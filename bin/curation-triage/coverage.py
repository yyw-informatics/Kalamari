#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tier 1 J1a coverage-probe core logic: join our picks to NCBI's ASSEMBLY_REPORTS,
type the match, and preview the pick-health flags (pure & offline).

Given (a) the eligible genome-units from Tier 0's ``taxon_policy.tsv``, (b) each
unit's source lane from ``chromosomes.tsv``, and (c) the distilled NCBI facts from
``ani_reports.AssemblyReportsCache``, this module answers the make-or-break question
the coverage probe exists to measure:

    For each pick, is there a *live* row for OUR assembly in
      (a) ANI_report_prokaryotes.txt  and  (b) assembly_summary_{refseq,genbank}.txt?
    And if so, what does NCBI's own verdict say about it?

It contains **no network and no I/O of its own** (the caller loads the inputs and
writes the outputs), so every rule here is unit-testable against tiny fixtures --
symmetric with how Tier 0 kept ``policy.py`` pure behind ``ncbi.py``.

Match typing (per file), most-informative first:
  * ``exact``          -- our assembly, our version, is a live row.
  * ``version-differs``-- our base assembly is live but at a *newer* version (our
                          recorded accession is stale; genome is essentially the same).
  * ``taxon-only``     -- our exact assembly is absent, but the species is present as
                          some *other* (newer) assembly (the thin-coverage case the
                          design memo warned about).
  * ``none``           -- neither the assembly nor the species is present.

J1a flags (honouring NCBI's ``approved-mismatch``, and NOT firing ANI-identity flags
on the confounded taxa, where an ANI/type-strain mismatch is *expected* and is
resolved by markers, never by ANI):
  * ``tax_check_failed`` / ``tax_check_inconclusive``  -- NCBI's taxonomy-check verdict.
  * ``best_match_species_mismatch``                    -- best type-strain ANI match is a
                                                          different species than declared.
  * ``ncbi_reclassified``                              -- NCBI now files this assembly under
                                                          a different species than declared.
  * ``refseq_excluded_identity`` / ``refseq_excluded_qc`` -- excluded_from_refseq populated,
                                                          split into identity vs annotation-QC.
  * ``assembly_superseded``                            -- our version is stale (see above).
  * ``assembly_suppressed`` / ``assembly_replaced``    -- version_status is not 'latest'.
  * ``genbank_only``                                   -- GenBank-only pick (never in RefSeq).
  * ``not_in_ani_report`` / ``not_in_assembly_summary``-- the signal is missing for this pick.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set

from ani_reports import (
    AssemblyReportsCache,
    AniReportRow,
    AssemblySummaryRow,
    accession_base,
    split_accession_version,
)

# ``best-match-status`` value NCBI uses to mark an assignment it has reviewed and
# approved despite an ANI mismatch -- we must NEVER flag these as mislabels.
APPROVED_MISMATCH = "approved-mismatch"

# Source-lane collapse: the many SME sub-sources become one "SME" lane; the design
# treats NCBI-REF / SME / FDA-ARGOS / NCTC differently (J1b), so keep them distinct.
SOURCE_LANE = {
    "NCBI-REF": "NCBI-REF",
    "NBCI-REF": "NCBI-REF",   # tolerate the known upstream typo
    "CDC SME": "SME",
    "UGA SME": "SME",
    "StaPH-B SME": "SME",
    "FDA-ARGOS": "FDA-ARGOS",
    "NCTC3000": "NCTC3000",
}

# Substrings in an ``excluded_from_refseq`` reason that mark an IDENTITY/taxonomy or
# provenance problem (high-value J1a: the organism may not be what it claims), as
# opposed to an annotation/assembly-quality one (``misassembled``, ``fragmented``,
# ``sequence duplications``, ``many frameshifted proteins`` -- kept as QC).  Terms are
# from NCBI's excluded_from_refseq controlled vocabulary.
IDENTITY_EXCLUSION_HINTS = (
    "unverified source organism",
    "contaminated",
    "chimeric",
    "mixed culture",
    "hybrid",
    "genus undefined",
    "untrustworthy as type",
    "metagenome",
    "derived from environmental",
    "derived from single cell",
)

MATCH_EXACT = "exact"
MATCH_VERSION = "version-differs"
MATCH_TAXON = "taxon-only"
MATCH_NONE = "none"


def is_real_assembly_base(base: str) -> bool:
    """True only for a real genome-assembly accession base (GCF_/GCA_).

    A Tier-0 unit whose assembly did not resolve carries a ``strain:<taxid>`` /
    ``acc:<x>`` proxy instead; such a base can never match an ASSEMBLY_REPORTS row and
    must be kept out of the join set (and surfaced), not silently mis-bucketed.
    """
    return base.startswith(("GCF_", "GCA_"))


def _names_differ(a: str, b: str) -> bool:
    """True unless both names are known and equal (case-insensitive)."""
    a, b = (a or "").strip().lower(), (b or "").strip().lower()
    return a != b if (a and b) else True


def _version_of(base: str, versioned_accs) -> Optional[int]:
    """Version int of the accession in ``versioned_accs`` whose base == ``base``."""
    for acc in versioned_accs:
        b, v = split_accession_version(acc)
        if b == base:
            return v
    return None


# --------------------------------------------------------------------------- #
#  Inputs                                                                      #
# --------------------------------------------------------------------------- #
@dataclass
class Pick:
    """One eligible Tier-0 genome-unit, as a join candidate."""

    unit_id: str
    scientific_name: str
    species_taxid: Optional[int]
    species_name: str
    parent_assembly: str          # versioned GCF/GCA
    accession_set: List[str]      # nuccore accessions (for source-lane lookup)
    source_lane: str = "unknown"
    regime: str = ""              # confounded regime A/B ('' if not confounded)
    resolver: str = ""            # marker / ANI

    @property
    def base(self) -> str:
        return accession_base(self.parent_assembly)

    @property
    def version(self) -> Optional[int]:
        return split_accession_version(self.parent_assembly)[1]

    @property
    def is_confounded(self) -> bool:
        return bool(self.regime)


def load_source_lanes(chromosomes_path: str) -> Dict[str, str]:
    """Map each nuccore accession base -> its source lane (collapsed)."""
    out: Dict[str, str] = {}
    with open(chromosomes_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            acc = accession_base(r["nuccoreAcc"])
            out[acc] = SOURCE_LANE.get(r["source"].strip(), r["source"].strip())
    return out


def load_picks(policy_path: str, source_lanes: Dict[str, str]) -> List[Pick]:
    """Load the eligible genome-units from the Tier-0 policy table as picks."""
    picks: List[Pick] = []
    with open(policy_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            if r["eligibility"] != "eligible":
                continue
            accs = [a for a in r["accession_set"].split(",") if a]
            lanes = {source_lanes.get(accession_base(a)) for a in accs}
            lanes.discard(None)
            lane = next(iter(lanes)) if len(lanes) == 1 else ("mixed" if lanes else "unknown")
            st = r.get("species_taxid", "").strip()
            picks.append(
                Pick(
                    unit_id=r["unit_id"],
                    scientific_name=r["scientificName"],
                    species_taxid=int(st) if st.isdigit() else None,
                    species_name=r.get("species_name", "").strip(),
                    parent_assembly=r["parent_assembly"].strip(),
                    accession_set=accs,
                    source_lane=lane,
                    regime=r.get("regime", "").strip(),
                    resolver=r.get("resolver", "").strip(),
                )
            )
    return picks


# --------------------------------------------------------------------------- #
#  Per-pick coverage record                                                    #
# --------------------------------------------------------------------------- #
@dataclass
class CoverageRow:
    pick: Pick
    ani_match: str = MATCH_NONE
    summary_match: str = MATCH_NONE
    in_refseq: bool = False
    newer_version_available: bool = False   # NCBI's live version > our recorded one
    # ANI-report verdict fields
    taxonomy_check: str = ""
    best_match_status: str = ""
    best_match_species_taxid: Optional[int] = None
    best_match_species_name: str = ""
    ncbi_species_taxid: Optional[int] = None
    ncbi_species_name: str = ""
    # assembly_summary fields
    version_status: str = ""
    excluded_from_refseq: str = ""
    relation_to_type_material: str = ""
    seq_rel_date: str = ""
    flags: List[str] = field(default_factory=list)

    @property
    def release_year(self) -> Optional[int]:
        y = self.seq_rel_date[:4]
        return int(y) if y.isdigit() else None

    @property
    def actionable_flags(self) -> List[str]:
        """Flags that by themselves put a pick on the SME worklist.

        Precision-first: a differing best-ANI-match (``best_match_species_mismatch``)
        is NOT actionable on its own -- NCBI's ``taxonomy-check-status`` already
        integrates the best match, and for a correctly-named-but-divergent reference
        (its type strain absent from the ANI panel) that status stays OK.  So we defer
        the yes/no to NCBI's verdict and keep the best-match only as the *explanation*.
        """
        return [f for f in self.flags if f in ACTIONABLE_FLAGS]

    def as_row(self) -> Dict[str, str]:
        p = self.pick
        return {
            "unit_id": p.unit_id,
            "scientificName": p.scientific_name,
            "source_lane": p.source_lane,
            "species_taxid": str(p.species_taxid or ""),
            "declared_species": p.species_name,
            "parent_assembly": p.parent_assembly,
            "ani_match": self.ani_match,
            "summary_match": self.summary_match,
            "in_refseq": "yes" if self.in_refseq else "no",
            "taxonomy_check": self.taxonomy_check,
            "best_match_status": self.best_match_status,
            "best_match_species": self.best_match_species_name,
            "ncbi_current_species": self.ncbi_species_name,
            "version_status": self.version_status,
            "excluded_from_refseq": self.excluded_from_refseq,
            "relation_to_type_material": self.relation_to_type_material,
            "seq_rel_date": self.seq_rel_date,
            "confounded": f"{self.pick.regime}/{self.pick.resolver}" if self.pick.is_confounded else "",
            "flags": ";".join(self.flags),
        }


REPORT_COLUMNS = [
    "unit_id", "scientificName", "source_lane", "species_taxid", "declared_species",
    "parent_assembly", "ani_match", "summary_match", "in_refseq", "taxonomy_check",
    "best_match_status", "best_match_species", "ncbi_current_species", "version_status",
    "excluded_from_refseq", "relation_to_type_material", "seq_rel_date", "confounded", "flags",
]

# Flags that by themselves put a pick on the SME worklist (drift/withdrawal signals
# where NCBI's own verdict -- or a hard status -- carries the call).  Everything
# else (best_match context, annotation-QC exclusions, stale versions, coverage
# notes, confounded/approved ambiguity) is informational.
ACTIONABLE_FLAGS = frozenset({
    "tax_check_failed",
    "tax_check_inconclusive",
    "ncbi_reclassified",
    "refseq_excluded_identity",
    "assembly_suppressed",
    "assembly_replaced",
})

# The demotions below are deliberately NARROW: each removes only the flags that are
# genuinely UNRELIABLE for its situation, and lets every orthogonal signal survive so a
# real problem is never hidden.
#
# Confounded taxa (ANI can't split them; markers do): the ANI/type-strain SIMILARITY
# verdicts are the expected artifact, so demote those.  But ``ncbi_reclassified`` (from
# NCBI's assigned taxid -- metadata, not ANI) and ``refseq_excluded_identity`` (an active
# curation withdrawal) are NOT ANI-similarity artifacts, so they stay actionable even
# for a confounded pick (else a reclassified or contaminated Salmonella/Listeria
# reference would be silently hidden as mere ANI ambiguity).
CONFOUNDED_DEMOTE_FLAGS = frozenset({
    "tax_check_failed", "tax_check_inconclusive", "best_match_species_mismatch",
})

# NCBI-approved mismatch: NCBI has reviewed and approved the assembly's SPECIES call
# despite a best-match discrepancy, so demote the species-identity flags.  The
# independent ``taxonomy-check-status`` verdict and any RefSeq withdrawal are separate
# axes and survive.
APPROVED_DEMOTE_FLAGS = frozenset({
    "best_match_species_mismatch", "ncbi_reclassified",
})


# --------------------------------------------------------------------------- #
#  Match typing                                                                #
# --------------------------------------------------------------------------- #
def _match_type(base: str, our_version: Optional[int], candidate_versioned: Sequence[str]) -> Optional[str]:
    """exact / version-differs for a base found among versioned candidate accessions."""
    for acc in candidate_versioned:
        b, v = split_accession_version(acc)
        if b == base:
            if our_version is not None and v is not None and our_version != v:
                return MATCH_VERSION
            return MATCH_EXACT
    return None


def type_ani_match(pick: Pick, cache: AssemblyReportsCache) -> tuple:
    # Fallback via the GCF<->GCA twin: an ANI row is computed per GenBank assembly and
    # its RefSeq accession column is occasionally blank, so a GCF pick can be present
    # only under its GCA twin base.
    twin = cache.twin_of.get(pick.base)
    row = cache.ani.get(pick.base) or (cache.ani.get(twin) if twin else None)
    if row is not None:
        # Version-type ONLY against our own accession family.  A twin-only match (the
        # blank-RefSeq-column case) matches on the GCA accession, whose version advances
        # INDEPENDENTLY of our GCF version -- comparing across the GCF<->GCA boundary
        # would fabricate a spurious 'version-differs', so a twin-only match is exact.
        mt = _match_type(pick.base, pick.version, row.versioned_accessions())
        return (mt or MATCH_EXACT), row
    if pick.species_taxid and cache.ani_species_present.get(pick.species_taxid, 0) > 0:
        return MATCH_TAXON, None
    return MATCH_NONE, None


def type_summary_match(pick: Pick, cache: AssemblyReportsCache) -> tuple:
    """Return (match_type, in_refseq, chosen_summary_row).

    Prefer the RefSeq row (authoritative version_status for a live pick); fall back
    to the GenBank row, which is where a pick that was dropped from RefSeq lives
    (with its ``excluded_from_refseq`` reason).
    """
    rs = cache.refseq.get(pick.base)
    gb = cache.genbank.get(pick.base)
    row = rs or gb
    if row is not None:
        cand = [row.assembly_accession, row.gbrs_paired_asm]
        mt = _match_type(pick.base, pick.version, cand) or MATCH_EXACT
        return mt, (rs is not None), row
    if pick.species_taxid and cache.summary_species_present.get(pick.species_taxid, 0) > 0:
        return MATCH_TAXON, False, None
    return MATCH_NONE, False, None


# --------------------------------------------------------------------------- #
#  Flag logic                                                                  #
# --------------------------------------------------------------------------- #
def classify_exclusion(reason: str) -> Optional[str]:
    """'' -> None; identity-type reason -> 'identity'; else -> 'qc'."""
    if not reason:
        return None
    low = reason.lower()
    if any(h in low for h in IDENTITY_EXCLUSION_HINTS):
        return "identity"
    return "qc"


def compute_flags(row: CoverageRow, ani: Optional[AniReportRow],
                  summary: Optional[AssemblySummaryRow]) -> List[str]:
    """The J1a pick-health flags for one coverage row (order = display order)."""
    flags: List[str] = []
    pick = row.pick

    # --- Tier-0 resolution gap: no real GCF/GCA to join on ----------------
    if not is_real_assembly_base(pick.base):
        flags.append("unresolved_assembly")

    # --- coverage gaps (informational for the probe) ----------------------
    if row.ani_match == MATCH_NONE:
        flags.append("not_in_ani_report")
    elif row.ani_match == MATCH_TAXON:
        flags.append("ani_taxon_only")
    if row.summary_match == MATCH_NONE:
        flags.append("not_in_assembly_summary")

    approved = bool(ani and ani.best_match_status == APPROVED_MISMATCH)

    # --- ANI-report identity verdict --------------------------------------
    if ani is not None:
        if ani.taxonomy_check_status == "Failed":
            flags.append("tax_check_failed")
        elif ani.taxonomy_check_status == "Inconclusive":
            flags.append("tax_check_inconclusive")
        if (ani.best_match_species_taxid and pick.species_taxid
                and ani.best_match_species_taxid != pick.species_taxid):
            flags.append("best_match_species_mismatch")
        if (ani.species_taxid and pick.species_taxid
                and ani.species_taxid != pick.species_taxid
                and _names_differ(ani.species_name, pick.species_name)):
            # taxid AND (when both known) name differ -> a real reclassification, not
            # a same-name taxid merge.
            flags.append("ncbi_reclassified")

    # --- assembly_summary withdrawal / provenance -------------------------
    excl_class = classify_exclusion(row.excluded_from_refseq)
    if excl_class == "identity":
        flags.append("refseq_excluded_identity")
    elif excl_class == "qc":
        flags.append("refseq_excluded_qc")

    if summary is not None:
        if summary.version_status == "suppressed":
            flags.append("assembly_suppressed")
        elif summary.version_status == "replaced":
            flags.append("assembly_replaced")
        # GenBank-only: matched only in GenBank and has no RefSeq twin.
        if (not row.in_refseq and summary.source_db == "genbank"
                and not summary.gbrs_paired_asm and excl_class is None):
            flags.append("genbank_only")

    if row.newer_version_available:
        flags.append("assembly_superseded")

    # --- demote EXPECTED identity ambiguity (narrowly) --------------------
    # Collapse only the flags that are genuinely unreliable for the situation; every
    # orthogonal signal (metadata reclassification, RefSeq withdrawal, and -- for an
    # approved mismatch -- NCBI's independent taxonomy-check verdict) survives.
    if pick.is_confounded and any(f in CONFOUNDED_DEMOTE_FLAGS for f in flags):
        flags = [f for f in flags if f not in CONFOUNDED_DEMOTE_FLAGS]
        flags.append("confounded_ani_ambiguity")
    elif approved and any(f in APPROVED_DEMOTE_FLAGS for f in flags):
        flags = [f for f in flags if f not in APPROVED_DEMOTE_FLAGS]
        flags.append("approved_mismatch")

    return flags


# --------------------------------------------------------------------------- #
#  Orchestration                                                               #
# --------------------------------------------------------------------------- #
def build_coverage(picks: Sequence[Pick], cache: AssemblyReportsCache) -> List[CoverageRow]:
    rows: List[CoverageRow] = []
    for pick in picks:
        ani_match, ani = type_ani_match(pick, cache)
        summary_match, in_refseq, summary = type_summary_match(pick, cache)
        row = CoverageRow(pick=pick, ani_match=ani_match,
                          summary_match=summary_match, in_refseq=in_refseq)
        if ani is not None:
            row.taxonomy_check = ani.taxonomy_check_status
            row.best_match_status = ani.best_match_status
            row.best_match_species_taxid = ani.best_match_species_taxid
            row.best_match_species_name = ani.best_match_species_name
            row.ncbi_species_taxid = ani.species_taxid
            row.ncbi_species_name = ani.species_name
        # authoritative excluded_from_refseq: whichever summary row carries it.
        excl = ""
        for r in (cache.refseq.get(pick.base), cache.genbank.get(pick.base)):
            if r and r.excluded_from_refseq:
                excl = r.excluded_from_refseq
                break
        row.excluded_from_refseq = excl
        if summary is not None:
            row.version_status = summary.version_status
            row.relation_to_type_material = summary.relation_to_type_material
            row.seq_rel_date = summary.seq_rel_date
        # 'superseded' only when NCBI's live version is strictly NEWER than ours
        # (a same-base version difference alone is direction-agnostic).
        if pick.version is not None:
            file_versions = []
            if ani is not None:
                file_versions.append(_version_of(pick.base, ani.versioned_accessions()))
            if summary is not None:
                file_versions.append(
                    _version_of(pick.base, [summary.assembly_accession, summary.gbrs_paired_asm]))
            row.newer_version_available = any(
                v is not None and v > pick.version for v in file_versions)
        row.flags = compute_flags(row, ani, summary)
        rows.append(row)
    rows.sort(key=lambda r: r.pick.scientific_name.lower())
    return rows


def write_report(rows: Sequence[CoverageRow], path: str) -> None:
    # LF terminators (not csv's default CRLF) so the trailing 'flags' field never
    # carries a stray '\r' that a naive `cut -f19` downstream would inherit.
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=REPORT_COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for r in rows:
            writer.writerow(r.as_row())


# --------------------------------------------------------------------------- #
#  Summary / breakdowns (pure -> a plain dict the CLI renders)                 #
# --------------------------------------------------------------------------- #
def _live(match: str) -> bool:
    """A live row for OUR assembly = exact or version-differs (same base assembly)."""
    return match in (MATCH_EXACT, MATCH_VERSION)


def _assembly_kind(row: "CoverageRow") -> str:
    b = row.pick.base
    if b.startswith("GCF_"):
        return "GCF"
    if b.startswith("GCA_"):
        return "GCA"
    return "other"


def _year_bucket(year: Optional[int]) -> str:
    if year is None:
        return "unknown"
    if year <= 2010:
        return "<=2010"
    if year <= 2015:
        return "2011-2015"
    if year <= 2020:
        return "2016-2020"
    return ">=2021"


def _rate_table(rows: Sequence[CoverageRow], key_fn) -> Dict[str, Dict[str, int]]:
    out: Dict[str, Dict[str, int]] = {}
    for r in rows:
        k = key_fn(r)
        d = out.setdefault(k, {"n": 0, "ani_live": 0, "summary_live": 0, "in_refseq": 0})
        d["n"] += 1
        d["ani_live"] += int(_live(r.ani_match))
        d["summary_live"] += int(_live(r.summary_match))
        d["in_refseq"] += int(r.in_refseq)
    return out


def summarize(rows: Sequence[CoverageRow]) -> dict:
    """Coverage headline + all mandated breakdowns + the J1a flag preview."""
    n = len(rows)
    ani_live = sum(_live(r.ani_match) for r in rows)
    summary_live = sum(_live(r.summary_match) for r in rows)

    def _match_counts(attr: str) -> Dict[str, int]:
        c: Dict[str, int] = {}
        for r in rows:
            m = getattr(r, attr)
            c[m] = c.get(m, 0) + 1
        return c

    # J1a flag preview: count each flag reason, split confounded vs not.
    flag_counts: Dict[str, int] = {}
    actionable_units = 0
    actionable_by_lane: Dict[str, int] = {}
    for r in rows:
        for f in r.flags:
            flag_counts[f] = flag_counts.get(f, 0) + 1
        if r.actionable_flags:
            actionable_units += 1
            actionable_by_lane[r.pick.source_lane] = actionable_by_lane.get(r.pick.source_lane, 0) + 1

    return {
        "n_picks": n,
        "headline": {
            "ani_live": ani_live,
            "ani_live_pct": round(100.0 * ani_live / n, 1) if n else 0.0,
            "summary_live": summary_live,
            "summary_live_pct": round(100.0 * summary_live / n, 1) if n else 0.0,
            "in_refseq": sum(r.in_refseq for r in rows),
        },
        "ani_match_types": _match_counts("ani_match"),
        "summary_match_types": _match_counts("summary_match"),
        "by_source_lane": _rate_table(rows, lambda r: r.pick.source_lane),
        "by_assembly_kind": _rate_table(rows, _assembly_kind),
        "by_year": _rate_table(rows, lambda r: _year_bucket(r.release_year)),
        "flag_counts": flag_counts,
        "actionable_units": actionable_units,
        "actionable_by_lane": actionable_by_lane,
        "confounded_units": sum(1 for r in rows if r.pick.is_confounded),
    }
