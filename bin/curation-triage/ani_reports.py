#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cached, mockable NCBI ASSEMBLY_REPORTS layer for the Tier 1 J1a coverage probe.

This is the *only* part of ``coverage_probe.py`` that needs the network, and it is
deliberately quarantined here so the join / match-typing / flag logic in
``coverage.py`` stays pure and offline -- exactly how Tier 0 quarantined the
network in ``ncbi.py`` behind ``policy.py``/``taxonomy.py``.

Three NCBI flat files carry NCBI's *own* precomputed verdict on every prokaryotic
assembly -- the signal the pick-health sentinel (J1a) reads:

    * ``ANI_report_prokaryotes.txt`` -- NCBI's nightly ANI-vs-type-strain check:
      taxonomy-check-status (OK/Failed/Inconclusive), best-match species, and the
      NCBI-blessed ``best-match-status`` values (incl. ``approved-mismatch``,
      which we must NEVER flag).
    * ``assembly_summary_refseq.txt`` / ``assembly_summary_genbank.txt`` --
      per-assembly ``version_status`` (latest/replaced/suppressed),
      ``excluded_from_refseq`` (the withdrawal reason), ``relation_to_type_material``,
      and ``seq_rel_date`` (assembly age).

Two traps these files already solve for us:

  * **GCF<->GCA equivalence.** The ANI report carries *both* the GenBank (GCA) and
    RefSeq (GCF) accession on every row, and each ``assembly_summary`` row carries
    its twin in ``gbrs_paired_asm`` -- so a pick stored as a GCF matches a row keyed
    on its GCA twin and vice-versa, for free.
  * **Accession version.** Accessions are matched on their *version-stripped* base
    (``GCF_000008685``), and the version we hold is compared against the version the
    live file reports so a superseded pick (our ``.1`` vs NCBI's ``.2``) is surfaced.

The files are large (~0.9 / 0.2 / 1.6 GB), so a fetch **streams** each file line by
line and keeps only the rows that touch our ~250 picks plus per-species presence
counts.  That *distilled* result is what gets committed to the JSON cache, so a
normal probe run and every unit test are fully offline and the committed artifact
stays small.  The stream openers are injectable (``opener=``) so tests never hit the
wire.
"""

from __future__ import annotations

import json
import os
import urllib.request
from dataclasses import dataclass, asdict, fields
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Set, Tuple

ASSEMBLY_REPORTS = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS"
ANI_REPORT_URL = f"{ASSEMBLY_REPORTS}/ANI_report_prokaryotes.txt"
SUMMARY_REFSEQ_URL = f"{ASSEMBLY_REPORTS}/assembly_summary_refseq.txt"
SUMMARY_GENBANK_URL = f"{ASSEMBLY_REPORTS}/assembly_summary_genbank.txt"

# "not applicable" sentinels NCBI uses in these files; treated as empty.
_NA = frozenset({"", "na", "NA", "n/a", "-"})


# --------------------------------------------------------------------------- #
#  Column layout of the real files (0-based).  Verified against the live       #
#  headers 2026-07.  ``expected`` header names are asserted at parse time so   #
#  an NCBI schema change fails loudly (an "epoch bump") instead of silently    #
#  reading the wrong column.                                                   #
# --------------------------------------------------------------------------- #
# ANI_report_prokaryotes.txt (25 cols, header line starts with '# genbank-accession')
ANI_COLS = {
    "genbank_accession": 0,
    "refseq_accession": 1,
    "species_taxid": 3,
    "organism_name": 4,
    "species_name": 5,
    "excluded_from_refseq": 8,   # anomaly note in this file (not authoritative exclusion)
    "best_match_species_taxid": 16,
    "best_match_species_name": 17,
    "best_match_status": 22,     # <- 'approved-mismatch' lives here
    "comment": 23,
    "taxonomy_check_status": 24,
}
ANI_HEADER_CHECK = {0: "genbank-accession", 1: "refseq-accession",
                    24: "taxonomy-check-status", 22: "best-match-status"}

# The SAME file also names the type-strain assembly NCBI compared each row against,
# with its ANI and both coverages.  Tier 2 needs that accession to know what to skani
# a pick against, so resolving a type strain costs no extra query -- and NCBI's own
# numbers become a free cross-check on ours.  Kept in a separate constant so the
# Tier-1 distillation and its committed cache are untouched by Tier 2's needs.
TYPE_COLS = {
    "declared_type_assembly": 9,       # type strain of the species the pick DECLARES
    "declared_type_organism_name": 10,
    "declared_type_category": 11,
    "declared_type_ani": 12,
    "declared_type_qcoverage": 13,
    "declared_type_scoverage": 14,
    "best_match_type_assembly": 15,    # type strain the pick actually matched best
}
TYPE_HEADER_CHECK = {9: "declared-type-assembly", 12: "declared-type-ANI",
                     15: "best-match-type-assembly"}

# assembly_summary_{refseq,genbank}.txt (>=38 cols; 2nd comment line is the header
# starting with '#assembly_accession')
SUMMARY_COLS = {
    "assembly_accession": 0,
    "refseq_category": 4,
    "taxid": 5,
    "species_taxid": 6,
    "organism_name": 7,
    "version_status": 10,        # latest / replaced / suppressed
    "seq_rel_date": 14,          # YYYY/MM/DD
    "gbrs_paired_asm": 17,       # the GCF<->GCA twin
    "excluded_from_refseq": 20,  # authoritative withdrawal reason
    "relation_to_type_material": 21,
}
SUMMARY_HEADER_CHECK = {0: "assembly_accession", 10: "version_status",
                        17: "gbrs_paired_asm", 20: "excluded_from_refseq"}


# --------------------------------------------------------------------------- #
#  Accession helpers                                                           #
# --------------------------------------------------------------------------- #
def split_accession_version(acc: str) -> Tuple[str, Optional[int]]:
    """``GCF_000008685.2`` -> (``GCF_000008685``, 2); ``GCF_1`` -> (``GCF_1``, None)."""
    acc = (acc or "").strip()
    if "." in acc:
        base, _, ver = acc.rpartition(".")
        if ver.isdigit():
            return base, int(ver)
    return acc, None


def accession_base(acc: str) -> str:
    return split_accession_version(acc)[0]


def _clean(value: str) -> str:
    """Normalise an NCBI 'na'-style empty field to ''."""
    v = (value or "").strip()
    return "" if v in _NA else v


# --------------------------------------------------------------------------- #
#  Distilled row records                                                       #
# --------------------------------------------------------------------------- #
@dataclass
class AniReportRow:
    """The ANI-report facts J1a needs for one assembly (version-stripped keys)."""

    genbank_accession: str = ""      # versioned, e.g. GCA_000008685.2
    refseq_accession: str = ""       # versioned, e.g. GCF_000008685.2 ('' if none)
    species_taxid: Optional[int] = None      # NCBI's *current* species for this assembly
    species_name: str = ""
    best_match_species_taxid: Optional[int] = None
    best_match_species_name: str = ""
    best_match_status: str = ""      # species-match / mismatch / approved-mismatch / ...
    taxonomy_check_status: str = ""  # OK / Failed / Inconclusive
    excluded_from_refseq: str = ""   # anomaly note (secondary to assembly_summary)
    comment: str = ""

    def versioned_accessions(self) -> List[str]:
        return [a for a in (self.genbank_accession, self.refseq_accession) if a]


@dataclass
class AssemblySummaryRow:
    """The assembly_summary facts J1a needs for one assembly."""

    assembly_accession: str = ""     # versioned
    refseq_category: str = ""        # 'reference genome' / 'representative genome' / ''
    taxid: Optional[int] = None
    species_taxid: Optional[int] = None
    version_status: str = ""         # latest / replaced / suppressed
    seq_rel_date: str = ""           # YYYY/MM/DD
    gbrs_paired_asm: str = ""        # twin accession (versioned)
    excluded_from_refseq: str = ""   # authoritative withdrawal reason ('' if none)
    relation_to_type_material: str = ""
    source_db: str = ""              # 'refseq' or 'genbank'

    @property
    def release_year(self) -> Optional[int]:
        y = self.seq_rel_date[:4]
        return int(y) if y.isdigit() else None


# --------------------------------------------------------------------------- #
#  Header validation + row parsing                                            #
# --------------------------------------------------------------------------- #
def _validate_header(fields_list: List[str], expected: Dict[int, str], what: str) -> None:
    for idx, name in expected.items():
        got = fields_list[idx].lstrip("#").strip() if idx < len(fields_list) else "<missing>"
        if got != name:
            raise ValueError(
                f"{what}: unexpected column layout -- expected {name!r} at index "
                f"{idx} but found {got!r}. NCBI may have changed the file schema "
                f"(treat as an epoch bump and update ANI_COLS/SUMMARY_COLS)."
            )


def _assembly_core(acc_base: str) -> str:
    """The numeric core of an assembly accession base: ``GCF_000283675`` -> ``000283675``.

    A GCF and its canonical GCA twin share this core; only the prefix differs.
    """
    for p in ("GCF_", "GCA_"):
        if acc_base.startswith(p):
            return acc_base[len(p):]
    return acc_base


def parse_ani_row(fields_list: List[str]) -> AniReportRow:
    c = ANI_COLS
    st = _clean(fields_list[c["species_taxid"]])
    bt = _clean(fields_list[c["best_match_species_taxid"]])
    return AniReportRow(
        genbank_accession=_clean(fields_list[c["genbank_accession"]]),
        refseq_accession=_clean(fields_list[c["refseq_accession"]]),
        species_taxid=int(st) if st.isdigit() else None,
        species_name=_clean(fields_list[c["species_name"]]),
        best_match_species_taxid=int(bt) if bt.isdigit() else None,
        best_match_species_name=_clean(fields_list[c["best_match_species_name"]]),
        best_match_status=_clean(fields_list[c["best_match_status"]]),
        taxonomy_check_status=_clean(fields_list[c["taxonomy_check_status"]]),
        excluded_from_refseq=_clean(fields_list[c["excluded_from_refseq"]]),
        comment=_clean(fields_list[c["comment"]]),
    )


def parse_summary_row(fields_list: List[str], source_db: str) -> AssemblySummaryRow:
    c = SUMMARY_COLS
    tx = _clean(fields_list[c["taxid"]])
    st = _clean(fields_list[c["species_taxid"]])
    return AssemblySummaryRow(
        assembly_accession=_clean(fields_list[c["assembly_accession"]]),
        refseq_category=_clean(fields_list[c["refseq_category"]]),
        taxid=int(tx) if tx.isdigit() else None,
        species_taxid=int(st) if st.isdigit() else None,
        version_status=_clean(fields_list[c["version_status"]]),
        seq_rel_date=_clean(fields_list[c["seq_rel_date"]]),
        gbrs_paired_asm=_clean(fields_list[c["gbrs_paired_asm"]]),
        excluded_from_refseq=_clean(fields_list[c["excluded_from_refseq"]]),
        relation_to_type_material=_clean(fields_list[c["relation_to_type_material"]]),
        source_db=source_db,
    )


# --------------------------------------------------------------------------- #
#  Streaming distillers (pure over an iterator of lines)                        #
# --------------------------------------------------------------------------- #
@dataclass
class AniDistill:
    rows: Dict[str, AniReportRow]        # wanted base accession -> row
    species_present: Dict[int, int]      # species_taxid -> count of assemblies seen
    n_scanned: int = 0


@dataclass
class SummaryDistill:
    rows: Dict[str, AssemblySummaryRow]  # wanted base accession -> row
    species_present: Dict[int, int]      # species_taxid -> count of 'latest' assemblies
    species_reference: Set[int]          # species with a reference/representative genome
    n_scanned: int = 0


def distill_type_strains(lines: Iterable[str],
                         wanted_bases: Set[str]) -> Dict[str, Dict[str, str]]:
    """Keep the type-strain columns of the ANI report for our picks (Tier 2).

    Returns ``{assembly_base: {declared_type_assembly, ..., species_taxid}}``.  A
    separate pass from ``distill_ani_report`` on purpose: Tier 1's committed cache is
    frozen, and widening it would force a rebuild of a file three tiers depend on just
    to add a column only Tier 2 reads.

    NCBI's own ``declared-type-ANI`` and coverages come along, which lets Tier 2's
    independently computed skani numbers be checked against the source they are meant
    to confirm.
    """
    out: Dict[str, Dict[str, str]] = {}
    header_seen = False
    min_cols = max(TYPE_HEADER_CHECK)
    for line in lines:
        if line.startswith("#"):
            if not header_seen:
                fields = line.rstrip("\n").split("\t")
                if len(fields) > min_cols:
                    _validate_header(fields, TYPE_HEADER_CHECK,
                                     "ANI_report_prokaryotes.txt (type columns)")
                    header_seen = True
            continue
        line = line.rstrip("\n")
        if not line or not header_seen:
            continue
        f = line.split("\t")
        if len(f) <= max(TYPE_COLS.values()):
            continue
        gbk = accession_base(_clean(f[ANI_COLS["genbank_accession"]]))
        rfq = accession_base(_clean(f[ANI_COLS["refseq_accession"]]))
        hits = [b for b in (gbk, rfq) if b and b in wanted_bases]
        if not hits:
            continue
        rec = {name: _clean(f[idx]) for name, idx in TYPE_COLS.items()}
        rec["species_taxid"] = f[ANI_COLS["species_taxid"]].strip()
        rec["species_name"] = _clean(f[ANI_COLS["species_name"]])
        for base in hits:
            out[base] = rec
    return out


def distill_ani_report(
    lines: Iterable[str],
    wanted_bases: Set[str],
    wanted_species: Set[int],
) -> AniDistill:
    """Keep ANI-report rows touching our picks; count assemblies per wanted species.

    A row matches a pick when *either* its GenBank or its RefSeq accession, stripped
    of version, is in ``wanted_bases`` (this is what makes GCF<->GCA equivalence
    free).  ``species_present`` counts, for each of *our* species taxids, how many
    assemblies the report holds -- used to classify a taxon-only match if a pick's
    exact assembly is ever absent.
    """
    rows: Dict[str, AniReportRow] = {}
    species_present: Dict[int, int] = {}
    n = 0
    header_seen = False
    min_cols = max(ANI_HEADER_CHECK)   # the real header has all columns; a '## note' does not
    for line in lines:
        if line.startswith("#"):
            # Validate the actual COLUMN header (the wide comment line), not a leading
            # '## description' comment -- and only that line.
            if not header_seen:
                fields = line.rstrip("\n").split("\t")
                if len(fields) > min_cols:
                    _validate_header(fields, ANI_HEADER_CHECK, "ANI_report_prokaryotes.txt")
                    header_seen = True
            continue
        line = line.rstrip("\n")
        if not line:
            continue
        if not header_seen:
            raise ValueError("ANI_report_prokaryotes.txt: reached data rows before a "
                             "validated column header (NCBI schema change?).")
        n += 1
        f = line.split("\t")
        if len(f) <= ANI_COLS["taxonomy_check_status"]:
            continue
        # match on the CLEANED base so an NCBI 'na'/'-' sentinel never matches.
        gbk = accession_base(_clean(f[ANI_COLS["genbank_accession"]]))
        rfq = accession_base(_clean(f[ANI_COLS["refseq_accession"]]))
        matched = (gbk and gbk in wanted_bases) or (rfq and rfq in wanted_bases)
        st = f[ANI_COLS["species_taxid"]].strip()
        if st.isdigit():
            sti = int(st)
            if sti in wanted_species:
                species_present[sti] = species_present.get(sti, 0) + 1
        if matched:
            row = parse_ani_row(f)
            # gbk and rfq are the two accessions of the SAME assembly, so no
            # cross-genome collision here (unlike gbrs_paired in assembly_summary).
            for base in (gbk, rfq):
                if base and base in wanted_bases:
                    rows[base] = row
    return AniDistill(rows=rows, species_present=species_present, n_scanned=n)


def distill_assembly_summary(
    lines: Iterable[str],
    wanted_bases: Set[str],
    wanted_species: Set[int],
    source_db: str,
) -> SummaryDistill:
    """Keep assembly_summary rows touching our picks; per-species presence counts.

    Matches on the assembly's own accession base or its ``gbrs_paired_asm`` twin
    base.  ``species_present`` counts only ``version_status == 'latest'`` assemblies
    for our species (a taxon-only match should point at a *live* alternative, not a
    replaced/suppressed one); ``species_reference`` records species that have a
    reference/representative genome available.
    """
    rows: Dict[str, AssemblySummaryRow] = {}
    own_keyed: Set[str] = set()
    species_present: Dict[int, int] = {}
    species_reference: Set[int] = set()
    n = 0
    header_seen = False
    min_cols = max(SUMMARY_HEADER_CHECK)  # the wide column header, not the '## See ...' note
    what = f"assembly_summary_{source_db}.txt"
    for line in lines:
        if line.startswith("#"):
            if not header_seen:
                fields = line.rstrip("\n").split("\t")
                if len(fields) > min_cols:
                    _validate_header(fields, SUMMARY_HEADER_CHECK, what)
                    header_seen = True
            continue
        line = line.rstrip("\n")
        if not line:
            continue
        if not header_seen:
            raise ValueError(f"{what}: reached data rows before a validated column "
                             f"header (NCBI schema change?).")
        n += 1
        f = line.split("\t")
        if len(f) <= SUMMARY_COLS["relation_to_type_material"]:
            continue
        acc = accession_base(_clean(f[SUMMARY_COLS["assembly_accession"]]))
        twin = accession_base(_clean(f[SUMMARY_COLS["gbrs_paired_asm"]]))
        matched = (acc and acc in wanted_bases) or (twin and twin in wanted_bases)
        st = f[SUMMARY_COLS["species_taxid"]].strip()
        vstatus = f[SUMMARY_COLS["version_status"]].strip()
        if st.isdigit() and int(st) in wanted_species:
            sti = int(st)
            if vstatus == "latest":
                species_present[sti] = species_present.get(sti, 0) + 1
            cat = f[SUMMARY_COLS["refseq_category"]].strip().lower()
            if "reference" in cat or "representative" in cat:
                species_reference.add(sti)
        if matched:
            row = parse_summary_row(f, source_db)
            # An own-accession match is authoritative and always wins; a gbrs_paired
            # twin match is only a fallback (NCBI re-pairing means a twin base can
            # belong to a DIFFERENT genome), so it must never overwrite an
            # own-accession row for the same base.
            if acc and acc in wanted_bases:
                rows[acc] = row
                own_keyed.add(acc)
            # A twin binding is only trusted when it is the CANONICAL pairing (same
            # numeric core, GCF_x <-> GCA_x); this rejects a foreign genome that names
            # our base as its paired accession via a stale/re-paired gbrs entry.
            if (twin and twin in wanted_bases and twin != acc and twin not in own_keyed
                    and _assembly_core(acc) == _assembly_core(twin)):
                rows.setdefault(twin, row)
    return SummaryDistill(rows=rows, species_present=species_present,
                          species_reference=species_reference, n_scanned=n)


# --------------------------------------------------------------------------- #
#  Line openers (network / file).  Injectable so tests never hit the wire.     #
# --------------------------------------------------------------------------- #
def open_url_lines(url: str, timeout: int = 300) -> Iterator[str]:
    """Stream a (plain-text) NCBI file line by line without buffering it all."""
    req = urllib.request.Request(url, headers={"User-Agent": "kalamari-curation-triage"})
    resp = urllib.request.urlopen(req, timeout=timeout)  # noqa: S310 - fixed NCBI host
    for raw in resp:
        yield raw.decode("utf-8", "replace")


def open_file_lines(path: str) -> Iterator[str]:
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            yield line


# --------------------------------------------------------------------------- #
#  Cache-backed resolver                                                       #
# --------------------------------------------------------------------------- #
class AssemblyReportsCache:
    """Distilled ANI-report + assembly_summary facts for a fixed set of picks.

    Mirrors ``ncbi.NuccoreResolver``: a JSON cache is loaded if present; a build
    streams the three live files (via injectable openers) and distils only the rows
    touching ``wanted_bases`` plus per-species presence counts, then persists.  In
    ``--offline`` mode nothing is streamed and a cache that does not cover the
    requested picks raises ``LookupError``.

    The cache is *query-scoped* (it is a function of the live files AND the pick
    set), so ``meta.wanted_bases`` records the set it was built for; a probe over a
    different pick set must rebuild it.
    """

    def __init__(
        self,
        cache_path: Optional[str] = None,
        ani_opener: Callable[[], Iterable[str]] = None,
        refseq_opener: Callable[[], Iterable[str]] = None,
        genbank_opener: Callable[[], Iterable[str]] = None,
    ) -> None:
        self.cache_path = cache_path
        self._ani_opener = ani_opener or (lambda: open_url_lines(ANI_REPORT_URL))
        self._refseq_opener = refseq_opener or (lambda: open_url_lines(SUMMARY_REFSEQ_URL))
        self._genbank_opener = genbank_opener or (lambda: open_url_lines(SUMMARY_GENBANK_URL))

        self.wanted_bases: Set[str] = set()
        self.wanted_species: Set[int] = set()
        self.twin_of: Dict[str, str] = {}   # pick base -> its GCF<->GCA twin base
        self.ani: Dict[str, AniReportRow] = {}
        self.refseq: Dict[str, AssemblySummaryRow] = {}
        self.genbank: Dict[str, AssemblySummaryRow] = {}
        self.ani_species_present: Dict[int, int] = {}
        self.summary_species_present: Dict[int, int] = {}
        self.summary_species_reference: Set[int] = set()
        self.meta: Dict[str, dict] = {}

        if cache_path and os.path.exists(cache_path):
            self._load_cache(cache_path)

    # ---- persistence -----------------------------------------------------
    def _load_cache(self, path: str) -> None:
        with open(path) as fh:
            blob = json.load(fh)
        self.wanted_bases = set(blob.get("wanted_bases", []))
        self.wanted_species = {int(x) for x in blob.get("wanted_species", [])}
        self.twin_of = dict(blob.get("twin_of", {}))
        self.ani = {b: _row_from_dict(AniReportRow, d) for b, d in blob.get("ani_rows", {}).items()}
        self.refseq = {b: _row_from_dict(AssemblySummaryRow, d)
                       for b, d in blob.get("summary_refseq", {}).items()}
        self.genbank = {b: _row_from_dict(AssemblySummaryRow, d)
                        for b, d in blob.get("summary_genbank", {}).items()}
        self.ani_species_present = {int(k): v for k, v in blob.get("ani_species_present", {}).items()}
        self.summary_species_present = {int(k): v
                                        for k, v in blob.get("summary_species_present", {}).items()}
        self.summary_species_reference = {int(x) for x in blob.get("summary_species_reference", [])}
        self.meta = blob.get("meta", {})

    def save(self) -> None:
        if not self.cache_path:
            return
        os.makedirs(os.path.dirname(os.path.abspath(self.cache_path)), exist_ok=True)
        payload = {
            "schema": 1,
            "meta": self.meta,
            "wanted_bases": sorted(self.wanted_bases),
            "wanted_species": sorted(self.wanted_species),
            "twin_of": dict(sorted(self.twin_of.items())),
            "ani_rows": {b: asdict(r) for b, r in sorted(self.ani.items())},
            "summary_refseq": {b: asdict(r) for b, r in sorted(self.refseq.items())},
            "summary_genbank": {b: asdict(r) for b, r in sorted(self.genbank.items())},
            "ani_species_present": {str(k): v for k, v in sorted(self.ani_species_present.items())},
            "summary_species_present": {str(k): v
                                        for k, v in sorted(self.summary_species_present.items())},
            "summary_species_reference": sorted(self.summary_species_reference),
        }
        tmp = f"{self.cache_path}.tmp.{os.getpid()}"
        try:
            with open(tmp, "w") as fh:
                json.dump(payload, fh, indent=1, sort_keys=True)
                fh.write("\n")
            os.replace(tmp, self.cache_path)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)

    # ---- build / resolve -------------------------------------------------
    def covers(self, bases: Iterable[str], species: Optional[Iterable[int]] = None) -> bool:
        """True iff the committed cache was distilled for (a superset of) these picks.

        Checks the species too: a Tier-0 retaxonomy can change a pick's
        ``species_taxid`` while its assembly base is unchanged, which would make the
        cached per-species presence counts stale -- so an offline run must rebuild.
        """
        if not set(bases).issubset(self.wanted_bases):
            return False
        if species is not None and not set(int(s) for s in species).issubset(self.wanted_species):
            return False
        return True

    def build(
        self,
        wanted_bases: Set[str],
        wanted_species: Set[int],
        include_genbank: bool = True,
        log=lambda m: None,
    ) -> None:
        """Stream the three live files and distil facts for the given pick set.

        Order matters: the assembly_summaries are streamed FIRST so we learn each
        pick's GCF<->GCA twin (via ``gbrs_paired_asm``), then the ANI report is
        streamed for the picks AND their twins -- so a GCF pick whose ANI row lists
        only the GenBank accession (RefSeq column blank) is still matched via its twin.
        """
        self.wanted_bases = set(wanted_bases)
        self.wanted_species = set(wanted_species)

        log("Streaming assembly_summary_refseq ...")
        rs = distill_assembly_summary(self._refseq_opener(), wanted_bases, wanted_species, "refseq")
        self.refseq = rs.rows
        log(f"  refseq summary: {rs.n_scanned} rows scanned, {len(rs.rows)} picks matched")

        gb = SummaryDistill(rows={}, species_present={}, species_reference=set())
        if include_genbank:
            log("Streaming assembly_summary_genbank (large) ...")
            gb = distill_assembly_summary(self._genbank_opener(), wanted_bases,
                                          wanted_species, "genbank")
            log(f"  genbank summary: {gb.n_scanned} rows scanned, {len(gb.rows)} picks matched")
        self.genbank = gb.rows

        # Learn each pick's GCF<->GCA twin from the summary pairing.
        self.twin_of = {}
        for base in wanted_bases:
            row = self.refseq.get(base) or self.genbank.get(base)
            if row is not None:
                twin = _twin_from_row(row, base)
                if twin and twin.startswith(("GCF_", "GCA_")):
                    self.twin_of[base] = twin
        ani_wanted = set(wanted_bases) | set(self.twin_of.values())

        log(f"Streaming ANI report ({len(wanted_bases)} picks, +{len(self.twin_of)} twins) ...")
        a = distill_ani_report(self._ani_opener(), ani_wanted, wanted_species)
        self.ani = a.rows
        self.ani_species_present = a.species_present
        log(f"  ANI report: {a.n_scanned} rows scanned, {len(a.rows)} rows matched")

        # Merge per-species presence across summaries; union the reference set.
        self.summary_species_present = dict(rs.species_present)
        for sp, c in gb.species_present.items():
            self.summary_species_present[sp] = max(self.summary_species_present.get(sp, 0), c)
        self.summary_species_reference = rs.species_reference | gb.species_reference

        self.meta = {
            "ani_report": {"n_scanned": a.n_scanned, "n_matched": len(a.rows)},
            "assembly_summary_refseq": {"n_scanned": rs.n_scanned, "n_matched": len(rs.rows)},
            "assembly_summary_genbank": {"n_scanned": gb.n_scanned, "n_matched": len(gb.rows)},
        }


def _twin_from_row(row: AssemblySummaryRow, base: str) -> Optional[str]:
    """The GCF<->GCA twin base of ``base`` given a summary row that matched it.

    The row matched ``base`` either on its own accession or on its ``gbrs_paired_asm``;
    the twin is whichever of the two is NOT ``base``.
    """
    own = accession_base(row.assembly_accession)
    paired = accession_base(row.gbrs_paired_asm)
    if own == base:
        return paired or None
    if paired == base:
        return own or None
    return None


def _row_from_dict(cls, d: dict):
    known = {f.name for f in fields(cls)}
    return cls(**{k: v for k, v in d.items() if k in known})
