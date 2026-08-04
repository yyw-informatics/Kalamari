#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Offline NCBI-taxonomy resolver for the Kalamari curation-triage tool (Tier 0).

This module is **pure and offline**: it parses NCBI ``taxdump`` files
(``nodes.dmp`` / ``names.dmp`` / ``merged.dmp``) with the same simple
``line.split`` idiom as ``bin/generate_sepia_reference.py`` and climbs the
lineage.  It never touches the network, so the taxonomy-resolution logic is
fully unit-testable against a tiny fixture dump.

Why a real taxdump and not the bundled ``src/taxonomy/build`` dump:
    The bundled ``build`` dump is only the *custom-modification overlay* (the 18
    synthetic ``9000000``-``9000017`` nodes).  The bundled ``deprecated`` dump is
    a stale, filtered snapshot that is missing 23 of the 268 real taxids in
    ``chromosomes.tsv`` (verified 2026-07).  So Tier 0 resolves against a full
    NCBI taxdump, with the ``build`` overlay layered on top so the synthetic
    Kalamari lineages (Listeria I-IV, Salmonella subspecies, C. botulinum
    groups) climb to their real parent species.

Resolution contract (see :meth:`TaxonomyResolver.resolve`):
    * apply ``merged.dmp`` (an old taxid that NCBI merged forwards is followed,
      e.g. ``46170`` Staphylococcus aureus subsp. aureus -> ``1280``);
    * climb parent pointers to the first ``species``-rank ancestor;
    * report the superkingdom (Bacteria / Archaea / Viruses / Eukaryota) by
      climbing to a known superkingdom taxid, not by trusting a rank string.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

# Superkingdom / domain anchor taxids.  We climb to one of these rather than
# matching a rank label, because NCBI has repeatedly reshuffled the rank of the
# "Viruses" node (superkingdom vs. acellular root vs. realm).
SUPERKINGDOM_TAXIDS: Dict[int, str] = {
    2: "Bacteria",
    2157: "Archaea",
    2759: "Eukaryota",
    10239: "Viruses",
    12884: "Viroids",
}

# Ranks we accept as "species-or-below" for the preflight guard.  Anything in
# this set (or literally rank == "species") satisfies the requirement that an
# eligible unit resolve to a real, species-or-below NCBI taxon.
SPECIES_OR_BELOW_RANKS = frozenset(
    {
        "species",
        "subspecies",
        "strain",
        "isolate",
        "serotype",
        "serogroup",
        "biotype",
        "forma specialis",
        "pathogroup",
        "varietas",
        "morph",
    }
)

# Ranks at genus level and above.  A declared taxon at one of these ranks is a
# legitimate higher-level Kalamari entry (e.g. a bare "Ruminococcus_sp." genus)
# that may soft-demote to review:genus-level rather than hard-failing.  Anything
# NOT in this set that fails to reach a species -- including the ambiguous
# "no rank" that NCBI uses for most strain leaf nodes, and None -- is treated as
# "should have resolved to a species" and HARD-FAILS the preflight (the safe
# default, so a strain/no-rank key can never silently pass).
GENUS_OR_ABOVE_RANKS = frozenset(
    {
        "genus", "subgenus", "section", "subsection", "series", "subseries",
        "species group", "species subgroup",
        "tribe", "subtribe",
        "family", "subfamily", "superfamily",
        "order", "suborder", "superorder", "infraorder", "parvorder",
        "class", "subclass", "infraclass", "superclass",
        "phylum", "subphylum", "superphylum",
        "kingdom", "subkingdom", "superkingdom",
        "domain", "realm", "clade", "cohort", "subcohort", "supercohort",
    }
)


def is_genus_or_above(rank) -> bool:
    return rank in GENUS_OR_ABOVE_RANKS

# Extra directories to search for a full NCBI taxdump, tried after ``--taxdump``
# and ``$NCBI_TAXDUMP``.
#
# **Deliberately empty in the repository, and it must stay empty.**  A hard-coded
# absolute path here would be wrong twice over: it publishes where one institution
# keeps its reference data, and it means nothing on anybody else's filesystem.  A
# site records its own locations OUTSIDE version control -- see
# :func:`site_taxdump_paths`.  This tuple is the programmatic hook a caller or a
# test can override.
DEFAULT_TAXDUMP_SEARCH: Tuple[str, ...] = ()

# A gitignored, machine-local list of taxdump directories: one per line, blank
# lines and ``#`` comments ignored.  It lives at the repo root so a developer sets
# it once per clone and never has to remember an environment variable.
SITE_TAXDUMP_FILE = ".taxdump-search"

# Same idea as an environment variable, for CI and one-off shells.  Several
# directories, separated the way ``$PATH`` separates them.
SITE_TAXDUMP_ENV = "NCBI_TAXDUMP_SEARCH"

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def site_taxdump_paths(repo_root: Optional[str] = None) -> List[str]:
    """Taxdump directories this machine has configured, in search order.

    Reads ``$NCBI_TAXDUMP_SEARCH`` first, then the gitignored ``.taxdump-search``
    file at the repo root.  Both are optional and a missing or unreadable file is
    simply no configuration -- with none of them set the resolver falls back to
    the bundled dumps, which is what a fresh clone does.
    """
    out: List[str] = []
    env = os.environ.get(SITE_TAXDUMP_ENV, "")
    out.extend(p for p in env.split(os.pathsep) if p.strip())
    path = os.path.join(repo_root or REPO_ROOT, SITE_TAXDUMP_FILE)
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.split("#", 1)[0].strip()
                if line:
                    out.append(line)
    except OSError:
        pass
    return out


@dataclass(frozen=True)
class Resolution:
    """Result of resolving one input taxid to an NCBI species."""

    input_taxid: int
    input_rank: Optional[str]  # rank of the (merge-followed) input node
    input_name: Optional[str]
    species_taxid: Optional[int]  # first species-rank ancestor, if any
    species_name: Optional[str]
    superkingdom: Optional[str]
    resolved_from: str  # provenance: "self" | "climb" | "merged" | "overlay" | "unknown"

    @property
    def is_species_or_below(self) -> bool:
        """True iff the input node itself is species-rank-or-below.

        Used by the preflight to tell a genuinely genus-level entry (soft
        review) apart from a species/strain key whose species ancestor is
        missing from the dump (hard failure)."""
        return self.input_rank in SPECIES_OR_BELOW_RANKS

    @property
    def resolves_to_species(self) -> bool:
        return self.species_taxid is not None


def _split_dmp(line: str) -> List[str]:
    """Parse one NCBI ``.dmp`` row into stripped fields.

    Rows look like ``9000000\t|\t1639\t|\tsubspecies\t|\t...`` — fields are
    separated by ``\\t|\\t`` and the row ends with ``\\t|``.  Splitting on the
    bare ``|`` and stripping whitespace is exactly what
    ``generate_sepia_reference.py`` does, so we match it.
    """
    return [field.strip() for field in line.rstrip("\n").split("|")]


class TaxonomyResolver:
    """Climb NCBI taxonomy offline, with a Kalamari custom-node overlay.

    Parameters
    ----------
    nodes, parents, ranks, names, merged:
        In-memory maps.  Most callers use :meth:`from_dumps` instead of building
        these by hand; the raw constructor exists so tests can inject a handful
        of nodes without touching disk.
    """

    def __init__(
        self,
        parents: Dict[int, int],
        ranks: Dict[int, str],
        names: Dict[int, str],
        merged: Optional[Dict[int, int]] = None,
    ) -> None:
        self.parents = parents
        self.ranks = ranks
        self.names = names
        self.merged = merged or {}

    # ------------------------------------------------------------------ load
    @classmethod
    def from_dumps(
        cls,
        base_taxdump: str,
        overlay_dirs: Optional[List[str]] = None,
    ) -> "TaxonomyResolver":
        """Build a resolver from a base taxdump plus optional overlay dumps.

        ``base_taxdump`` is a directory holding ``nodes.dmp``/``names.dmp`` (and
        optionally ``merged.dmp``).  ``overlay_dirs`` are layered on top in
        order, so the Kalamari ``build`` overlay (synthetic nodes) wins over the
        base dump.  Overlay dumps need not carry ``merged.dmp``.
        """
        parents: Dict[int, int] = {}
        ranks: Dict[int, str] = {}
        names: Dict[int, str] = {}
        merged: Dict[int, int] = {}

        def load_dir(path: str) -> None:
            nodes_path = os.path.join(path, "nodes.dmp")
            names_path = os.path.join(path, "names.dmp")
            merged_path = os.path.join(path, "merged.dmp")
            if os.path.exists(nodes_path):
                with open(nodes_path) as fh:
                    for line in fh:
                        f = _split_dmp(line)
                        if len(f) < 3 or not f[0]:
                            continue
                        taxid = int(f[0])
                        parents[taxid] = int(f[1])
                        ranks[taxid] = f[2]
            if os.path.exists(names_path):
                with open(names_path) as fh:
                    for line in fh:
                        if "scientific name" not in line:
                            continue
                        f = _split_dmp(line)
                        if len(f) < 2 or not f[0]:
                            continue
                        names[int(f[0])] = f[1]
            if os.path.exists(merged_path):
                with open(merged_path) as fh:
                    for line in fh:
                        f = _split_dmp(line)
                        if len(f) < 2 or not f[0]:
                            continue
                        merged[int(f[0])] = int(f[1])

        load_dir(base_taxdump)
        for overlay in overlay_dirs or []:
            load_dir(overlay)

        return cls(parents=parents, ranks=ranks, names=names, merged=merged)

    @staticmethod
    def find_base_taxdump(explicit: Optional[str] = None) -> Optional[str]:
        """Resolve which base taxdump directory to use.

        Priority: explicit path -> ``$NCBI_TAXDUMP`` -> this machine's configured
        search paths (``$NCBI_TAXDUMP_SEARCH``, then the gitignored
        ``.taxdump-search`` file) -> ``DEFAULT_TAXDUMP_SEARCH`` -> ``None`` (the
        caller should then fall back to the bundled dump).

        Nothing here is site-specific: a clone with no configuration finds no
        taxdump and says so, rather than silently depending on one institution's
        filesystem layout.
        """
        candidates: List[str] = []
        if explicit:
            candidates.append(explicit)
        env = os.environ.get("NCBI_TAXDUMP")
        if env:
            candidates.append(env)
        candidates.extend(site_taxdump_paths())
        candidates.extend(DEFAULT_TAXDUMP_SEARCH)
        for path in candidates:
            if path and os.path.exists(os.path.join(path, "nodes.dmp")):
                return path
        return None

    # -------------------------------------------------------------- resolve
    def _follow_merge(self, taxid: int) -> int:
        """Follow a chain of ``merged.dmp`` redirects (usually length 0 or 1)."""
        seen = set()
        while taxid in self.merged and taxid not in seen:
            seen.add(taxid)
            taxid = self.merged[taxid]
        return taxid

    def contains(self, taxid: int) -> bool:
        return self._follow_merge(taxid) in self.parents

    def lineage(self, taxid: int) -> List[int]:
        """Return the lineage from ``taxid`` up to the root (self first)."""
        taxid = self._follow_merge(taxid)
        out: List[int] = []
        seen = set()
        cur = taxid
        while cur not in seen:
            seen.add(cur)
            out.append(cur)
            parent = self.parents.get(cur)
            if parent is None or parent == cur or cur == 1:
                break
            cur = parent
        return out

    def superkingdom(self, taxid: int) -> Optional[str]:
        for node in self.lineage(taxid):
            if node in SUPERKINGDOM_TAXIDS:
                return SUPERKINGDOM_TAXIDS[node]
        return None

    def resolve(self, taxid: int) -> Resolution:
        """Resolve one input taxid to an NCBI species (see module docstring)."""
        merged_taxid = self._follow_merge(taxid)
        resolved_from = "merged" if merged_taxid != taxid else "self"

        if merged_taxid not in self.parents:
            return Resolution(
                input_taxid=taxid,
                input_rank=None,
                input_name=self.names.get(merged_taxid),
                species_taxid=None,
                species_name=None,
                superkingdom=None,
                resolved_from="unknown",
            )

        input_rank = self.ranks.get(merged_taxid)
        input_name = self.names.get(merged_taxid)
        superkingdom = self.superkingdom(merged_taxid)

        species_taxid: Optional[int] = None
        for node in self.lineage(merged_taxid):
            if self.ranks.get(node) == "species":
                species_taxid = node
                if node != merged_taxid and resolved_from == "self":
                    resolved_from = "climb"
                break

        return Resolution(
            input_taxid=taxid,
            input_rank=input_rank,
            input_name=input_name,
            species_taxid=species_taxid,
            species_name=self.names.get(species_taxid) if species_taxid else None,
            superkingdom=superkingdom,
            resolved_from=resolved_from,
        )

    def name(self, taxid: int) -> Optional[str]:
        return self.names.get(self._follow_merge(taxid))
