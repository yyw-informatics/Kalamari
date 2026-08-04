#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cached, mockable ncbi-datasets-cli layer for the Tier 1 J1b "better-rep" feed.

J1b asks, per *tracked* species: has a newer, NCBI-designated **reference /
representative** genome appeared since the last run?  The query is
``datasets summary genome taxon <taxid> --released-after <date> --as-json-lines``
(the plan's §4 signal).  ncbi-datasets-cli is not always on PATH and hitting NCBI is
non-deterministic, so -- exactly like ``ani_reports.py`` quarantines the ASSEMBLY_REPORTS
download -- this module quarantines the ``datasets`` call behind:

  * an **injectable fetcher** (default = the real subprocess; tests pass a lambda that
    returns fixture JSON-lines), and
  * a **query-scoped committed JSON cache** (facts distilled per species taxid), so a
    normal run and every unit test are fully offline.

The **pure gating** (which lanes may raise a better-rep flag, confounded suppression,
reference-vs-our-pick comparison) lives in ``leads.j1b_leads`` -- this module only
fetches + distils, keeping the network and the policy cleanly apart.

Epoch pinning: unlike the NCBI flat files (whose epoch is pinned only by header
asserts), ncbi-datasets-cli has a real version, so the cache ``meta`` records the
resolved ``datasets_version`` and the ``released_after`` window -- a version bump is a
tooling epoch (it re-opens dismissed J1b leads for re-review).

Query narrowing: the two feeds want DIFFERENT slices of NCBI, and asking for
the wide one everywhere is not viable -- an unfiltered ``taxon 28901`` (Salmonella
enterica) returns 135 000 records, and we track 151 species.  So:

  * **J1b** only ever fires on an NCBI-*designated* reference whose assembly base
    differs from our pick, so it asks with ``--reference`` (1-2 records per species) --
    the flag returns exactly the RefSeq-designated reference genome(s), which is the
    same set J1b would have filtered down to anyway.
  * **J2**'s wishlist uses ``_best_candidate`` (prefer a reference, else the newest
    genome), so a wishlist species needs the **unfiltered** query.  Those 19 taxids are
    small species; unfiltered they total ~3000 records, which is fine.

The per-taxid query mode is recorded in the cache and re-checked by ``covers()``, so a
reference-only entry can never silently satisfy a full-query need (that would make J2
miss every species whose only genome is not NCBI's reference).

The ``--released-after`` window is recorded the same way, and for the same reason.
``build()`` REPLACES the cache, so a windowed refresh is a *narrowing* rebuild, never a
top-up: every taxid's candidate list is re-fetched from inside the window alone.  A
narrowed entry is still *present*, and an empty candidate list is indistinguishable from
"NCBI has nothing" -- so without the recorded window a windowed cache would sail through
the ``--offline`` gate of a run that asked the unwindowed question and answer it with
silence (leads vanishing off the worklist rather than a coverage failure).  ``covers()``
therefore also demands a window at least as wide as the caller needs.  The window is a
J1b delta device only; J2's wishlist is always fetched unwindowed (see ``build``).
"""

from __future__ import annotations

import datetime
import json
import os
import subprocess
from dataclasses import asdict, dataclass, fields
from typing import Callable, Dict, Iterable, List, Optional, Set, Tuple

from ani_reports import accession_base   # reuse the version-stripped base helper

DEFAULT_DATASETS_BIN = "datasets"

# The two query shapes a cached taxid can have been built with (see module docstring).
# 'all' is the wide query and satisfies every need; 'reference' is the narrow
# ``--reference`` query and satisfies only a J1b-style need.
QUERY_ALL = "all"
QUERY_REFERENCE = "reference"

# The unwindowed query: no ``--released-after``, i.e. "everything NCBI has".
WINDOW_ALL = ""


# --------------------------------------------------------------------------- #
#  --released-after windows                                                    #
# --------------------------------------------------------------------------- #
def window_start(window: str) -> Optional[Tuple[int, int, int]]:
    """(year, month, day) of a ``--released-after`` value, or None if it is not a date.

    ncbi-datasets-cli documents MM/DD/YYYY and also accepts ISO YYYY-MM-DD, so both are
    understood rather than guessed at.
    """
    text = (window or "").strip()
    for fmt in ("%m/%d/%Y", "%Y-%m-%d"):
        try:
            when = datetime.datetime.strptime(text, fmt)
        except ValueError:
            continue
        return (when.year, when.month, when.day)
    return None


def window_covers(cached: str, needed: str) -> bool:
    """True when a cache fetched with window ``cached`` can answer window ``needed``.

    A window only ever REMOVES records, so the cached window has to start no later than
    the question being asked: the unwindowed cache answers everything, and an earlier
    start answers a later one -- never the reverse.  Two windows we cannot parse are
    trusted only when identical: refusing costs a rebuild, guessing costs an empty
    worklist that looks like good news.
    """
    cached, needed = (cached or "").strip(), (needed or "").strip()
    if cached == needed:
        return True
    if cached == WINDOW_ALL:
        return True                     # unwindowed cache = a superset of every window
    if needed == WINDOW_ALL:
        return False                    # a windowed cache cannot answer "everything"
    start_cached, start_needed = window_start(cached), window_start(needed)
    if start_cached is None or start_needed is None:
        return False
    return start_cached <= start_needed


# --------------------------------------------------------------------------- #
#  Distilled candidate record                                                 #
# --------------------------------------------------------------------------- #
@dataclass
class GenomeReport:
    """The datasets-summary facts J1b needs for one candidate assembly."""

    accession: str = ""              # versioned GCF_/GCA_
    refseq_category: str = ""        # 'reference genome' / 'representative genome' / ''
    release_date: str = ""           # YYYY-MM-DD
    species_taxid: Optional[int] = None
    organism_name: str = ""
    taxonomy_check_status: str = ""  # from the ANI block, when present
    source_db: str = ""              # 'refseq' / 'genbank' (from the accession prefix)

    @property
    def base(self) -> str:
        return accession_base(self.accession)

    @property
    def is_reference(self) -> bool:
        cat = (self.refseq_category or "").lower()
        return "reference" in cat or "representative" in cat


def _get(obj: dict, *path, default=""):
    """Safe nested lookup: _get(o, 'assembly_info', 'refseq_category')."""
    cur = obj
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur if cur is not None else default


def parse_datasets_report(obj: dict) -> GenomeReport:
    """One ``datasets summary genome`` JSON report -> GenomeReport (tolerant of missing keys)."""
    acc = _get(obj, "accession") or _get(obj, "current_accession")
    st = _get(obj, "organism", "tax_id", default=None)
    try:
        st = int(st) if st not in (None, "") else None
    except (TypeError, ValueError):
        st = None
    return GenomeReport(
        accession=acc,
        refseq_category=_get(obj, "assembly_info", "refseq_category"),
        release_date=_get(obj, "assembly_info", "release_date"),
        species_taxid=st,
        organism_name=_get(obj, "organism", "organism_name"),
        taxonomy_check_status=_get(obj, "average_nucleotide_identity",
                                   "taxonomy_check_status"),
        source_db="refseq" if str(acc).startswith("GCF_") else
                  ("genbank" if str(acc).startswith("GCA_") else ""),
    )


def distill_datasets(lines: Iterable[str]) -> List[GenomeReport]:
    """Pure: JSON-lines from ``datasets ... --as-json-lines`` -> candidate reports."""
    out: List[GenomeReport] = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        out.append(parse_datasets_report(json.loads(line)))
    return out


# --------------------------------------------------------------------------- #
#  Real fetcher (subprocess) -- injectable so tests never call it              #
# --------------------------------------------------------------------------- #
def datasets_command(taxid: int, released_after: str, datasets_bin: str,
                     reference_only: bool = False) -> List[str]:
    """The ncbi-datasets-cli invocation (centralised so the flags live in one place).

    ``reference_only`` adds ``--reference``, which restricts the answer to the
    RefSeq-designated reference genome(s) of the taxon.  That is the whole J1b
    question, and it is the difference between 1 record and 135 000 for a big species
    like Salmonella enterica -- see the module docstring.
    """
    cmd = [datasets_bin, "summary", "genome", "taxon", str(taxid), "--as-json-lines"]
    if reference_only:
        cmd += ["--reference"]
    if released_after:
        cmd += ["--released-after", released_after]   # MM/DD/YYYY per datasets-cli
    return cmd


def run_datasets(taxid: int, released_after: str,
                 datasets_bin: str = DEFAULT_DATASETS_BIN,
                 reference_only: bool = False) -> List[str]:
    """Call ``datasets`` and return its JSON-lines stdout (one report per line)."""
    cmd = datasets_command(taxid, released_after, datasets_bin, reference_only)
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)  # noqa: S603
    if proc.returncode != 0:
        raise RuntimeError("datasets failed for taxid %s (rc=%s): %s"
                           % (taxid, proc.returncode, proc.stderr.strip()[:400]))
    return proc.stdout.splitlines()


def datasets_version(datasets_bin: str = DEFAULT_DATASETS_BIN) -> str:
    """The pinned tool epoch string, e.g. 'ncbi-datasets-cli@16.15.0' (or @unknown)."""
    try:
        proc = subprocess.run([datasets_bin, "--version"], capture_output=True,
                              text=True, check=False)  # noqa: S603
        ver = (proc.stdout or proc.stderr).strip().split()[-1] if proc.returncode == 0 else "unknown"
    except (OSError, IndexError):
        ver = "unknown"
    return "ncbi-datasets-cli@%s" % ver


# --------------------------------------------------------------------------- #
#  Query-scoped committed cache                                                #
# --------------------------------------------------------------------------- #
def _report_from_dict(d: dict) -> GenomeReport:
    known = {f.name for f in fields(GenomeReport)}
    return GenomeReport(**{k: v for k, v in d.items() if k in known})


class DatasetsSummaryCache:
    """Distilled ``datasets summary`` candidates per species taxid, for a fixed
    ``released_after`` window.  Mirrors ``ani_reports.AssemblyReportsCache``:
    injectable fetcher, committed JSON cache, ``covers()`` gate for offline runs.
    """

    def __init__(
        self,
        cache_path: Optional[str] = None,
        fetcher: Callable[[int, str, bool], Iterable[str]] = None,
        datasets_bin: str = DEFAULT_DATASETS_BIN,
    ) -> None:
        self.cache_path = cache_path
        self.datasets_bin = datasets_bin
        # fetcher(taxid, released_after, reference_only) -> iterable of JSON-lines
        self._fetcher = fetcher or (
            lambda taxid, after, ref_only: run_datasets(taxid, after, datasets_bin, ref_only))
        self.by_species: Dict[int, List[GenomeReport]] = {}
        # taxid -> QUERY_ALL / QUERY_REFERENCE: which question we actually asked NCBI.
        self.query_mode: Dict[int, str] = {}
        # taxid -> the --released-after window that taxid was actually fetched with.
        self.window: Dict[int, str] = {}
        self.released_after: str = WINDOW_ALL
        self.meta: Dict[str, object] = {}
        if cache_path and os.path.exists(cache_path):
            self._load_cache(cache_path)

    def _load_cache(self, path: str) -> None:
        with open(path, encoding="utf-8") as fh:
            blob = json.load(fh)
        self.by_species = {int(k): [_report_from_dict(d) for d in v]
                           for k, v in blob.get("by_species", {}).items()}
        # A schema-1 cache carries no query_mode; every one of its entries came from
        # the unfiltered query, so 'all' is the truthful default.
        self.query_mode = {int(k): str(v)
                           for k, v in (blob.get("query_mode") or {}).items()}
        self.released_after = blob.get("released_after", WINDOW_ALL)
        # Schema <=2 recorded only the cache-level window, and every entry was fetched
        # with it, so ``window_of`` inherits it for any taxid with no window of its own.
        self.window = {int(k): str(v)
                       for k, v in (blob.get("window") or {}).items()}
        self.meta = blob.get("meta", {})

    def save(self) -> None:
        if not self.cache_path:
            return
        os.makedirs(os.path.dirname(os.path.abspath(self.cache_path)), exist_ok=True)
        payload = {
            "schema": 3,
            "meta": self.meta,
            "released_after": self.released_after,
            "query_mode": {str(k): self.mode_of(k) for k in sorted(self.by_species)},
            "window": {str(k): self.window_of(k) for k in sorted(self.by_species)},
            "by_species": {str(k): [asdict(r) for r in v]
                           for k, v in sorted(self.by_species.items())},
        }
        tmp = "%s.tmp.%d" % (self.cache_path, os.getpid())
        try:
            with open(tmp, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=1, sort_keys=True)
                fh.write("\n")
            os.replace(tmp, self.cache_path)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)

    def mode_of(self, taxid: int) -> str:
        """The query mode this taxid's cached candidate list was built with."""
        return self.query_mode.get(int(taxid), QUERY_ALL)

    def window_of(self, taxid: int) -> str:
        """The ``--released-after`` window this taxid's candidate list was fetched with.

        Falls back to the cache-level window, which is what a schema<=2 cache recorded
        (every entry there was fetched with it).  The committed cache has ``""``, so it
        keeps reading as unwindowed and still answers every question.
        """
        return self.window.get(int(taxid), self.released_after)

    def uncovered(self, species: Iterable[int], reference_only: Iterable[int] = (),
                  released_after: str = WINDOW_ALL) -> List[Tuple[int, str]]:
        """The requested taxids this cache cannot answer, as (taxid, why) pairs.

        Three distinct failures, deliberately reported apart:

          * a taxid the cache never queried ('absent');
          * one queried with ``--reference`` when the caller needs the wide answer,
            which would silently hide every non-reference genome from J2's wishlist
            feeder;
          * one fetched inside a ``--released-after`` window narrower than the question
            being asked.  ``build`` replaces rather than tops up, so a windowed refresh
            leaves entries that are present but pruned -- and a pruned entry answers
            "nothing found", which reads exactly like "NCBI has nothing".

        ``released_after`` is the window the caller needs for the *windowed* feed, i.e.
        J1b's ``reference_only`` species.  The wide-query taxids are J2's wishlist,
        whose question ("does this species have a genome at all?") is standing, not a
        delta, so they always need the unwindowed answer.
        """
        narrow: Set[int] = {int(t) for t in reference_only}
        out: List[Tuple[int, str]] = []
        for taxid in sorted({int(s) for s in species}):
            if taxid not in self.by_species:
                out.append((taxid, "absent"))
                continue
            if taxid not in narrow and self.mode_of(taxid) == QUERY_REFERENCE:
                out.append((taxid, "cached as reference-only, need the full query"))
                continue
            need = released_after if taxid in narrow else WINDOW_ALL
            have = self.window_of(taxid)
            if not window_covers(have, need):
                out.append((taxid, "cached only since %s, need %s"
                            % (have or "(all)", need or "(all)")))
        return out

    def covers(self, species: Iterable[int], reference_only: Iterable[int] = (),
               released_after: str = WINDOW_ALL) -> bool:
        """True when every requested taxid is cached AT LEAST as widely as it is needed.

        ``reference_only`` names the taxids the caller is happy to answer from the
        narrow ``--reference`` query (J1b's tracked species).  Everything else needs a
        ``QUERY_ALL`` entry -- a wide cache entry satisfies a narrow need, never the
        other way round.  ``released_after`` applies the same rule to the date window.
        """
        return not self.uncovered(species, reference_only, released_after)

    def build(self, species_taxids: Iterable[int], released_after: str,
              version_epoch: Optional[str] = None, log=lambda m: None,
              reference_only: Iterable[int] = ()) -> None:
        """Fetch + distil candidates for each species (via the injected fetcher).

        Taxids listed in ``reference_only`` are asked with ``--reference`` (J1b's
        narrow question); every other taxid gets the unfiltered query.  ``build``
        REPLACES the cache, so it must always be called over the whole union of taxids
        the enabled feeds need -- see ``tier1_metadata.datasets_cache``.

        ``released_after`` reaches only the ``reference_only`` taxids.  It is a J1b
        delta device ("has a newer reference appeared?"), while the wide-query taxids
        are J2's wishlist, which asks the standing question "does this species have a
        genome at all?".  A window would not narrow that answer, it would falsify it:
        a species whose only genome predates the window would read as "still no genome
        available" and re-raise a lead we already have a genome for.
        """
        narrow: Set[int] = {int(t) for t in reference_only}
        self.released_after = released_after
        self.by_species = {}
        self.query_mode = {}
        self.window = {}
        for taxid in sorted(set(int(t) for t in species_taxids)):
            ref_only = taxid in narrow
            window = released_after if ref_only else WINDOW_ALL
            reports = distill_datasets(self._fetcher(taxid, window, ref_only))
            self.by_species[taxid] = reports
            self.query_mode[taxid] = QUERY_REFERENCE if ref_only else QUERY_ALL
            self.window[taxid] = window
            log("datasets taxon %d [%s]: %d candidate genome(s) since %s"
                % (taxid, self.query_mode[taxid], len(reports), window or "(all)"))
        self.meta = {
            "datasets_version": version_epoch or "ncbi-datasets-cli@unset",
            "released_after": released_after,
            "n_species": len(self.by_species),
            "n_reference_only": sum(1 for m in self.query_mode.values()
                                    if m == QUERY_REFERENCE),
            "n_windowed": sum(1 for w in self.window.values() if w != WINDOW_ALL),
            "n_records": sum(len(v) for v in self.by_species.values()),
        }

    def candidates(self, taxid: Optional[int]) -> List[GenomeReport]:
        if taxid is None:
            return []
        return self.by_species.get(int(taxid), [])
