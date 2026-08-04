#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tier 1 lead model: turn the per-signal detections (J1a / J1b / J2) into uniform,
self-contained ``Lead`` records the ledger/worklist loop can diff, suppress, and
render (pure & offline).

A **lead** is one candidate SME work item.  The three signals feed it differently:

  * **J1a** (pick-health) -- reuses ``coverage.py`` unchanged: a lead per eligible
    pick whose ``CoverageRow.actionable_flags`` is non-empty (the frozen 6-flag
    worklist vocabulary).  The lead's identity is the ``unit_id`` and its evidence
    is the *set* of actionable flags, so the monthly delta can ask "did this pick's
    actionable-flag set GAIN a member since last run?" (see ``ledger.compute_delta``).
  * **J1b** (better-rep) -- a lead per (unit, newer candidate assembly) the
    source-gated feed surfaces; identity = ``unit_id#candidate``; delta = first time
    that candidate is seen.  [built in a later step]
  * **J2** (coverage feeder) -- a lead per newly-available genome for a species we
    want but don't have; identity = ``todo:<species>#<candidate>`` or
    ``lpsn:<species>``; delta = first appearance.  [built in a later step]

Two identity concepts drive the loop and are deliberately kept apart:

  * ``lead_key``            -- the *stable* identity of the lead across runs (what a
                              review row and the delta both key on).
  * ``evidence_fingerprint``-- the *content* hash of the lead's current evidence.  A
                              dismissal (``reviews.tsv``) is honoured only while the
                              fingerprint is unchanged, so a rejected lead re-surfaces
                              the moment its evidence genuinely changes (a new flag,
                              a newer candidate) -- "suppress until evidence changes".

``epoch`` is the *tooling/schema* epoch for the lead's signal (NOT the monthly run):
it moves only when NCBI changes a file layout or a pinned tool version bumps, so a
dismissal persists across ordinary monthly runs but a genuine recompute wave re-opens
every dismissed lead for re-review.  Keeping epoch OUT of the fingerprint and IN the
record lets the two concerns stay legible.

This module holds **no network and no I/O** -- ``coverage.py`` (+ its cache) produces
the rows, ``ledger.py`` persists/diffs the leads, and the CLI wires them together --
mirroring how Tier 0 and the coverage probe split pure logic from the cached download.
"""

from __future__ import annotations

import csv
import hashlib
from dataclasses import asdict, dataclass, field, fields
from typing import Dict, List, Optional, Sequence, Set

# --------------------------------------------------------------------------- #
#  Tooling / schema epochs.                                                    #
#                                                                              #
#  Bump ONLY on a real schema or pinned-tool change, NEVER for an ordinary     #
#  monthly data refresh -- a dismissal in reviews.tsv is scoped to its signal's #
#  epoch, so moving these re-opens every dismissed lead of that signal for      #
#  re-review (the "recompute wave = epoch bump" rule from the plan's §8).       #
# --------------------------------------------------------------------------- #
# J1a reads NCBI's ANI_report_prokaryotes + assembly_summary_{refseq,genbank};
# their column layout is pinned + header-verified in ani_reports.py ("2026-07").
# Move this when that verification date moves (an NCBI schema change).
ANI_EPOCH = "ncbi-assembly-reports@2026-07"
# Fallbacks for when the corresponding cache is absent, so a lead is never stamped with
# a blank epoch.  The REAL epochs come from the committed caches' meta blocks (the
# resolved ncbi-datasets-cli version, the LPSN release) -- see
# ``tier1_metadata.tool_versions`` and ``collect_j2``.
DATASETS_EPOCH = "ncbi-datasets-cli@unset"
LPSN_EPOCH = "lpsn@unset"
# Tier 2's ANI epoch (skani).  A marker lead's epoch is richer -- software + database
# + our interpretation rule -- and is composed per task by ``tier2.marker_epoch``.
SKANI_EPOCH = "skani@unset"

# Resolved run parameters stamped into every ledger record (plan §10 defaults).
# Tier-1 metadata itself uses none of the ANI thresholds, but stamping them keeps a
# record self-describing for the Tier-2 confirmation that consumes the survivors.
DEFAULT_PARAMS = {"species_ani": 95, "af_gate": 0.5}

SIGNAL_J1A = "J1a"
SIGNAL_J1B = "J1b"
SIGNAL_J2 = "J2"
SIGNAL_T2 = "T2"

# The evidence keys that define a Tier-2 lead's content (see ``Lead.fingerprint``).
# Deliberately the VERDICT and the call, never the raw ANI number: a skani version
# bump moves the last decimal place of every comparison, and busting every SME
# dismissal over 0.01% ANI would train reviewers to ignore the tool.  A genuine
# tooling change is handled the right way instead -- by the epoch, which is matched
# separately.
T2_FINGERPRINT_KEYS = ("t2_check", "verdict", "call", "nearest_label", "redundant_with")

# NCBI Datasets genome landing page (the SME's one-click "look at this" link).
GENOME_LINK_BASE = "https://www.ncbi.nlm.nih.gov/datasets/genome/"


def _fp(payload: str) -> str:
    """A short, stable content hash (first 16 hex of sha1) for an evidence payload."""
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:16]


def lpsn_epoch(release: str) -> str:
    """The LPSN tooling epoch string for a release tag, e.g. 'lpsn@2026-06-01'.

    The ``lpsn@`` prefix is never optional.  Stamping the bare release tag would let a
    run with no ``--lpsn-release`` write the naked literal ``unset`` -- an epoch that
    names no tool and sorts next to nothing.
    """
    return "lpsn@%s" % (str(release).strip() or "unset")


# --------------------------------------------------------------------------- #
#  The uniform lead record                                                     #
# --------------------------------------------------------------------------- #
@dataclass
class Lead:
    """One candidate SME work item, uniform across signals and self-contained.

    ``as_record()`` serialises it for the append-only ledger; ``fingerprint``
    derives the content hash that scopes a dismissal.
    """

    signal: str                 # J1a / J1b / J2
    lead_key: str               # stable cross-run identity (delta + review key)
    unit_id: str                # Kalamari unit ('' for a J2 species we have no unit for)
    scientific_name: str
    source_lane: str            # NCBI-REF / SME / FDA-ARGOS / NCTC3000 / ...
    epoch: str                  # tooling/schema epoch of this signal (see module docstring)
    actionable_flags: List[str] = field(default_factory=list)  # J1a: the frozen worklist flags
    flags: List[str] = field(default_factory=list)             # J1a: full flag context
    candidate_accession: str = ""     # J1b/J2: the newer/covering genome
    deposit_date: str = ""            # seq_rel_date of the pick (J1a) or candidate (J1b/J2)
    evidence: Dict[str, str] = field(default_factory=dict)     # human-facing verdict fields
    link: str = ""
    summary: str = ""                 # one-line human description for the worklist

    @property
    def fingerprint(self) -> str:
        """Content hash of the lead's *current evidence* (epoch-independent).

        A ``reviews.tsv`` dismissal is honoured only while this is unchanged, so the
        lead re-surfaces the instant its evidence genuinely changes -- for J1a when the
        actionable-flag SET changes, for J1b/J2 when the candidate assembly changes.
        Epoch is deliberately excluded (it is matched separately) so an ordinary
        monthly refresh never busts a dismissal.
        """
        if self.signal == SIGNAL_J1A:
            payload = "J1a|%s|%s" % (self.unit_id, ",".join(sorted(self.actionable_flags)))
        elif self.signal == SIGNAL_J1B:
            payload = "J1b|%s|%s" % (self.unit_id, self.candidate_accession)
        elif self.signal == SIGNAL_T2:
            # Tier 2: the verdict IS the evidence (the ANI number is context).
            payload = "T2|%s|%s" % (self.lead_key, ",".join(
                "%s=%s" % (k, self.evidence.get(k, "")) for k in T2_FINGERPRINT_KEYS))
        else:  # J2
            payload = "J2|%s|%s" % (self.lead_key, self.candidate_accession)
        return _fp(payload)

    def to_dict(self) -> dict:
        """Round-trippable plain dict (used for the sharded per-shard lead fragments)."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "Lead":
        known = {f.name for f in fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in known})

    def as_record(self, run_id: str, delta: str, params: Dict = None,
                  tool_versions: Dict = None) -> dict:
        """A self-contained NDJSON ledger record for this lead in a given run."""
        return {
            "schema": 1,
            "run_id": run_id,
            "signal": self.signal,
            "lead_key": self.lead_key,
            "epoch": self.epoch,
            "delta": delta,
            "unit_id": self.unit_id,
            "scientificName": self.scientific_name,
            "source_lane": self.source_lane,
            "candidate_accession": self.candidate_accession,
            "deposit_date": self.deposit_date,
            "actionable_flags": list(self.actionable_flags),
            "flags": list(self.flags),
            "evidence": dict(self.evidence),
            "evidence_fingerprint": self.fingerprint,
            "link": self.link,
            "summary": self.summary,
            "params": dict(params or DEFAULT_PARAMS),
            "tool_versions": dict(tool_versions or {}),
        }


# --------------------------------------------------------------------------- #
#  J1a -- pick-health leads (reuse coverage.py's CoverageRow unchanged)         #
# --------------------------------------------------------------------------- #
# Human-readable gloss for each actionable flag (worklist prose; keep terse).
_FLAG_GLOSS = {
    "tax_check_failed": "NCBI taxonomy-check FAILED",
    "tax_check_inconclusive": "NCBI taxonomy-check inconclusive",
    "ncbi_reclassified": "NCBI now files this assembly under a different species",
    "refseq_excluded_identity": "dropped from RefSeq for an identity/contamination reason",
    "assembly_suppressed": "assembly SUPPRESSED by NCBI",
    "assembly_replaced": "assembly REPLACED by NCBI",
}


def _j1a_summary(row) -> str:
    """One-line human description of a J1a pick-health lead."""
    p = row.pick
    reasons = "; ".join(_FLAG_GLOSS.get(f, f) for f in row.actionable_flags)
    context = []
    if "ncbi_reclassified" in row.actionable_flags and row.ncbi_species_name:
        context.append("NCBI species now %s" % row.ncbi_species_name)
    if row.best_match_species_name and "best_match_species_mismatch" in row.flags:
        context.append("best ANI match %s" % row.best_match_species_name)
    if row.excluded_from_refseq:
        context.append("excluded_from_refseq: %s" % row.excluded_from_refseq)
    tail = (" -- " + "; ".join(context)) if context else ""
    return "%s (%s, %s): %s%s" % (
        p.scientific_name, p.parent_assembly, p.source_lane, reasons, tail)


def j1a_leads(rows: Sequence["object"], epoch: str = ANI_EPOCH) -> List[Lead]:
    """Build J1a leads from the coverage rows: one per pick with actionable flags.

    ``rows`` are ``coverage.CoverageRow`` objects (imported lazily by the caller so
    this module stays dependency-light).  Precision is inherited wholesale from
    ``coverage.py``: we surface exactly the rows whose ``actionable_flags`` is
    non-empty -- the confounded-demotion and approved-mismatch rules already applied.
    """
    leads: List[Lead] = []
    for row in rows:
        actionable = list(row.actionable_flags)
        if not actionable:
            continue
        p = row.pick
        evidence = {
            "declared_species": p.species_name,
            "taxonomy_check": row.taxonomy_check,
            "best_match_status": row.best_match_status,
            "best_match_species": row.best_match_species_name,
            "ncbi_current_species": row.ncbi_species_name,
            "excluded_from_refseq": row.excluded_from_refseq,
            "relation_to_type_material": row.relation_to_type_material,
        }
        leads.append(Lead(
            signal=SIGNAL_J1A,
            lead_key=p.unit_id,
            unit_id=p.unit_id,
            scientific_name=p.scientific_name,
            source_lane=p.source_lane,
            epoch=epoch,
            actionable_flags=actionable,
            flags=list(row.flags),
            candidate_accession="",
            deposit_date=row.seq_rel_date,
            evidence=evidence,
            link="%s%s/" % (GENOME_LINK_BASE, p.parent_assembly),
            summary=_j1a_summary(row),
        ))
    return leads


# --------------------------------------------------------------------------- #
#  J1b -- source-gated better-rep feed                                          #
# --------------------------------------------------------------------------- #
# ONLY these lanes get a "a newer reference genome appeared" flag.  For SME /
# FDA-ARGOS / NCTC3000 picks, divergence from the type/reference strain is a
# DELIBERATE curatorial choice (e.g. MC58 for N. meningitidis), so J1b must never nag
# "a reference exists and you're not it" -- such a pick is surfaced only if it itself
# degrades, which is J1a's job.  (User-confirmed: FDA-ARGOS + NCTC3000 = like SME.)
J1B_LANES = frozenset({"NCBI-REF"})


def _j1b_summary(pick, cand) -> str:
    rel = (" released %s" % cand.release_date) if cand.release_date else ""
    return ("%s (%s, %s): NCBI's reference/representative genome is now %s%s, not our "
            "pick %s" % (pick.scientific_name, pick.parent_assembly, pick.source_lane,
                         cand.accession, rel, pick.parent_assembly))


def j1b_leads(picks: Sequence["object"], candidate_source, epoch: str = DATASETS_EPOCH):
    """Source-gated better-rep leads: a *new* NCBI reference genome that isn't our pick.

    ``candidate_source`` is any object with ``.candidates(taxid) -> [GenomeReport-like]``
    (e.g. ``datasets_summary.DatasetsSummaryCache``); its candidates are already scoped
    to genomes released since the last run, so any reference among them is genuinely new.

    Gating (precision-first): only NCBI-REF, non-confounded picks; a candidate must be
    NCBI-designated reference/representative AND a different assembly base than our pick.
    """
    out: List[Lead] = []
    for p in picks:
        if p.source_lane not in J1B_LANES:      # SME / FDA-ARGOS / NCTC3000: deliberate
            continue
        if getattr(p, "is_confounded", False):  # informational:confounded -> suppress
            continue
        if not p.species_taxid:
            continue
        for cand in candidate_source.candidates(p.species_taxid):
            if not cand.is_reference or cand.base == p.base:
                continue
            out.append(Lead(
                signal=SIGNAL_J1B,
                lead_key="%s#%s" % (p.unit_id, cand.base),
                unit_id=p.unit_id,
                scientific_name=p.scientific_name,
                source_lane=p.source_lane,
                epoch=epoch,
                candidate_accession=cand.accession,
                deposit_date=cand.release_date,
                evidence={
                    "declared_species": p.species_name,
                    "candidate_refseq_category": cand.refseq_category,
                    "candidate_taxonomy_check": cand.taxonomy_check_status,
                    "our_assembly": p.parent_assembly,
                },
                link="%s%s/" % (GENOME_LINK_BASE, cand.accession),
                summary=_j1b_summary(p, cand),
            ))
    return out


# --------------------------------------------------------------------------- #
#  J2 -- coverage feeders (chromosomes-todo.tsv wishlist + LPSN gaps)           #
# --------------------------------------------------------------------------- #
TODO_WISHLIST = "XXXXXX"          # nuccoreAcc sentinel: "a species we want but lack"


@dataclass
class TodoEntry:
    """One row of chromosomes-todo.tsv (the SME wishlist; 4 cols, no header)."""

    scientific_name: str          # as written (spaces), e.g. 'Arcobacter butzleri'
    nuccore_acc: str              # a real accession, or TODO_WISHLIST
    taxid: Optional[int]
    parent: str

    @property
    def is_wishlist(self) -> bool:
        return self.nuccore_acc.strip().upper() == TODO_WISHLIST


def normalize_species(name: str) -> str:
    """Canonical species key: lowercase, underscores->spaces, collapsed whitespace."""
    return " ".join((name or "").replace("_", " ").split()).lower()


def species_genus(name: str) -> str:
    n = normalize_species(name)
    return n.split(" ")[0] if n else ""


def load_todo(path: str) -> List[TodoEntry]:
    """Load chromosomes-todo.tsv: 4 columns, NO header, tab-separated."""
    out: List[TodoEntry] = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or not row[0].strip():
                continue
            name = row[0].strip()
            # A blank accession cell means the same thing as a missing column: "want it,
            # don't have it yet" -- route to the wishlist sentinel so it gets a genome
            # lookup, not a bogus zero-accession "ready to add" lead.
            acc = row[1].strip() if len(row) > 1 and row[1].strip() else TODO_WISHLIST
            tx = row[2].strip() if len(row) > 2 else ""
            parent = row[3].strip() if len(row) > 3 else ""
            out.append(TodoEntry(scientific_name=name, nuccore_acc=acc,
                                 taxid=int(tx) if tx.isdigit() else None, parent=parent))
    return out


def _best_candidate(cands):
    """Pick the most add-worthy candidate genome: reference/representative first, then
    the most recently released."""
    if not cands:
        return None
    refs = [c for c in cands if c.is_reference]
    pool = refs or list(cands)
    return sorted(pool, key=lambda c: (c.is_reference, c.release_date or ""), reverse=True)[0]


def valid_since(value: str, since: str) -> bool:
    """True when an LPSN valid-publication date is at or after the ``since`` cutoff.

    Compared as left-anchored ISO-ish text truncated to the shorter side, so the cutoff
    and the value need not have the same granularity.  In practice LPSN gives only a
    year (``lpsn.valid_publication_year``), so this is what lets a ``2024-06-01`` cutoff
    still be answerable against the bare year ``2024``.  A record with
    NO date is DROPPED once a cutoff is set: precision first -- an undated name cannot be
    shown to be new, and quietly treating it as new is how a "newly published species"
    feed turns into a dump of the whole nomenclature.
    """
    if not since:
        return True
    v = (value or "").strip()
    if not v:
        return False
    n = min(len(v), len(since))
    return v[:n] >= since[:n]


def j2_leads(todo_entries: Sequence[TodoEntry], candidate_source, lpsn_records,
             tracked_genera: Set[str], have_species: Set[str],
             todo_epoch: str = DATASETS_EPOCH, lpsn_epoch: str = LPSN_EPOCH,
             max_lpsn: Optional[int] = 25, since: Optional[str] = None,
             log=lambda m: None) -> List[Lead]:
    """Coverage-feeder leads: wishlist species that now have a genome + LPSN gaps.

    * **todo feeder** (highest precision -- explicit SME intent): one lead per wishlist
      species that now has an available genome (or a pre-identified accession).  Its
      candidate comes from ncbi-datasets-cli, so it carries the *datasets* tooling
      epoch (``todo_epoch``).
    * **LPSN feeder**: a validly-published species in a tracked genus that Kalamari
      lacks.  Carries the *LPSN release* epoch (``lpsn_epoch``).  Sorted NEWEST FIRST
      and capped at ``max_lpsn`` (with a logged count of any dropped) so a first run
      cannot flood the worklist; ``None`` disables the LPSN half entirely.

    This is where the LPSN gate lives, on purpose.  ``lpsn.py`` distils the download
    faithfully -- synonyms, undated rows and all -- and the two filters that decide what
    counts as a coverage gap are applied HERE, over the in-memory records:

      * **status**: only a *correct name* is a gap.  A ``synonym`` row is a name LPSN
        itself has retired; surfacing it tells an SME to add a species that does not
        exist, which is exactly the kind of false lead that trains reviewers to stop
        reading the worklist.
      * **date**: with ``since`` set, only names validly published on/after that cutoff.
        Without it, every name LPSN has ever published in a tracked genus is "new to
        Kalamari", so the first run is a backlog dump rather than a monthly signal.

    Keeping both in the pure layer means changing our mind costs a flag, not a
    re-download -- and the committed cache stays a complete, honest copy of LPSN.

    The two feeders carry SEPARATE epochs so a bump in one tool re-baselines only its
    own leads (a datasets-cli bump must not re-open LPSN gaps, and vice-versa).
    ``have_species`` / ``tracked_genera`` must be the WHOLE-catalog sets (never a shard
    slice) -- an LPSN gap is a global-coverage fact.
    """
    out: List[Lead] = []
    todo_species = {normalize_species(e.scientific_name) for e in todo_entries}

    # --- todo wishlist feeder -------------------------------------------------
    for e in todo_entries:
        sp_key = normalize_species(e.scientific_name)
        if e.is_wishlist:
            cand = _best_candidate(candidate_source.candidates(e.taxid)) if e.taxid else None
            if cand is None:
                continue                       # nothing to add yet
            acc, cat = cand.accession, cand.refseq_category or "genome"
            link = "%s%s/" % (GENOME_LINK_BASE, acc)
            summary = ("coverage feeder: wishlist species %s now has an available %s (%s)"
                       % (e.scientific_name, cat, acc))
        else:
            acc = e.nuccore_acc                # a pre-identified accession, ready to add
            link = ""
            summary = ("coverage feeder: wishlist species %s is ready to add "
                       "(targeted accession %s)" % (e.scientific_name, acc))
        out.append(Lead(
            signal=SIGNAL_J2,
            lead_key="todo:%s#%s" % (e.taxid if e.taxid else sp_key, acc),
            unit_id="",
            scientific_name=e.scientific_name,
            source_lane="todo",
            epoch=todo_epoch,
            candidate_accession=acc,
            evidence={"taxid": str(e.taxid or ""), "parent": e.parent,
                      "source": "chromosomes-todo.tsv"},
            link=link,
            summary=summary,
        ))

    # --- LPSN newly-valid-species gap feeder ---------------------------------
    if max_lpsn is not None and lpsn_records:
        in_genus = [r for r in lpsn_records if r.genus.strip().lower() in tracked_genera]
        correct = [r for r in in_genus if r.is_correct_name]
        dated = [r for r in correct if valid_since(getattr(r, "valid_date", ""), since or "")]
        gaps = [r for r in dated
                if normalize_species(r.species) not in have_species
                and normalize_species(r.species) not in todo_species]
        if len(correct) < len(in_genus) or len(dated) < len(correct):
            log("J2 LPSN: %d name(s) in tracked genera -> %d correct name(s) -> %d "
                "published since %s" % (len(in_genus), len(correct), len(dated),
                                        since or "(no cutoff)"))
        # NEWEST FIRST, and before the cap: this feeder exists to surface *newly
        # published* species, so the 25 the cap keeps must be the 25 most recent ones.
        # Sorted alphabetically instead, the real catalog gives 25 Acinetobacter out of
        # ~4000 gaps -- a cap that silently hides every genus after 'A'.  Undated names
        # sort last (an empty year), which is right: they cannot be shown to be new.
        # Two passes because Python's sort is stable: name ascending, then year
        # descending, so same-year gaps stay in a deterministic alphabetical order.
        gaps.sort(key=lambda r: normalize_species(r.species))
        gaps.sort(key=lambda r: getattr(r, "valid_date", "") or "", reverse=True)
        if len(gaps) > max_lpsn:
            log("J2 LPSN: %d gaps in tracked genera; showing the %d newest (dropped %d "
                "-- raise --max-lpsn to see more)"
                % (len(gaps), max_lpsn, len(gaps) - max_lpsn))
            gaps = gaps[:max_lpsn]
        for r in gaps:
            out.append(Lead(
                signal=SIGNAL_J2,
                lead_key="lpsn:%s" % normalize_species(r.species),
                unit_id="",
                scientific_name=r.species,
                source_lane="lpsn",
                epoch=lpsn_epoch,
                candidate_accession="",
                deposit_date=getattr(r, "valid_date", ""),
                evidence={"genus": r.genus, "status": r.status, "authority": r.authority,
                          "valid_date": getattr(r, "valid_date", ""), "source": "LPSN"},
                link="",
                summary=("coverage feeder: validly-published species %s (genus %s) is "
                         "in a tracked genus but not in Kalamari" % (r.species, r.genus)),
            ))
    return out
