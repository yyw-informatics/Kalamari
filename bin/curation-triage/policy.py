#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tier 0 core logic: build the genome-unit policy table (pure & offline).

Given (a) the parsed ``chromosomes.tsv`` rows, (b) a :class:`taxonomy.TaxonomyResolver`,
(c) a map of ``accession -> ncbi.NuccoreRecord``, and (d) a confounded-taxa
manifest, this module produces one :class:`GenomeUnit` per genome and runs the
preflight.  It contains **no network and no I/O of its own** (the caller loads
the inputs), so every rule here is unit-testable against tiny fixtures.

The five jobs, in order:

  1. **Group rows into genome-units by assembly, never by the Kalamari taxid.**
     The grouping key comes from the nuccore record (real GCF/GCA if known, else
     the record's own strain-level ``taxid``, else the accession).  This merges
     multi-replicon genomes (Aliivibrio's two chromosomes) into one unit and
     keeps distinct strains sharing a Kalamari taxid (the two Yersinia
     enterocolitica rows) as separate units.

  2. **Resolve a real NCBI species taxid** by climbing the taxonomy; if the
     declared taxid is genus-level (Ruminococcus_sp.) fall back to the genome's
     own organism taxid, which is usually a real species.

  3. **Classify eligibility with a programmatic guard** — organelle/oversized ->
     ``exclude:organelle``; virus -> ``exclude:virus``; bacteria/archaea ->
     ``eligible``; eukaryote -> ``eukaryote``; bacteria-with-no-species ->
     ``review:genus-level`` (surfaced, not silently passed).

  4. **Tag confounded taxa** from the committed manifest with ``regime`` (A/B)
     and ``resolver`` (marker/ANI).

  5. **Detect the shared-accession quirk** (one physical sequence labelled as
     two subspecies) and flag it instead of silently double-counting.

The preflight (:func:`run_preflight`) then HARD-FAILS if any ``eligible`` unit
lacks a species-rank taxid.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from taxonomy import TaxonomyResolver, is_genus_or_above
from ncbi import NuccoreRecord, strip_version, is_real_assembly

# Eligibility classes emitted in the policy table.
ELIG_ELIGIBLE = "eligible"
ELIG_ORGANELLE = "exclude:organelle"
ELIG_VIRUS = "exclude:virus"
ELIG_EUKARYOTE = "eukaryote"
ELIG_GENUS = "review:genus-level"
ELIG_UNRESOLVED = "review:unresolved"

# moltype substrings that indicate a segment / non-chromosomal element.
SEGMENT_MOLTYPE_HINTS = ("segment",)

# Default organelle/eukaryote size cap (bp).  A grouped unit summing above this
# is treated as an organelle/eukaryote assembly regardless of the genome field.
DEFAULT_SIZE_CAP_BP = 50_000_000


# --------------------------------------------------------------------------- #
#  Data structures                                                            #
# --------------------------------------------------------------------------- #
@dataclass
class ChromosomeRow:
    scientific_name: str
    accession: str
    taxid: int
    parent: int
    source: str


@dataclass
class ManifestEntry:
    trigger_taxid: int
    label: str
    regime: str
    resolver: str
    note: str = ""


@dataclass
class GenomeUnit:
    unit_id: str
    group_key: str
    scientific_names: List[str]
    kalamari_taxids: List[int]
    accession_set: List[str]
    parent_assembly: Optional[str]
    species_taxid: Optional[int] = None
    species_name: Optional[str] = None
    superkingdom: Optional[str] = None
    eligibility: str = ELIG_UNRESOLVED
    regime: str = ""
    resolver: str = ""
    total_length: Optional[int] = None
    n_replicons: int = 0
    flags: List[str] = field(default_factory=list)
    note: str = ""

    # ---- rendering helpers ------------------------------------------------
    def scientific_name_field(self) -> str:
        return ";".join(self.scientific_names)

    def kalamari_taxid_field(self) -> str:
        return ";".join(str(t) for t in self.kalamari_taxids)

    def accession_field(self) -> str:
        return ",".join(self.accession_set)

    def as_row(self) -> Dict[str, str]:
        """Render to the committed TSV schema (required 9 cols first)."""
        return {
            "unit_id": self.unit_id,
            "scientificName": self.scientific_name_field(),
            "kalamari_taxid": self.kalamari_taxid_field(),
            "species_taxid": str(self.species_taxid) if self.species_taxid else "",
            "accession_set": self.accession_field(),
            "parent_assembly": self.parent_assembly or "",
            "eligibility": self.eligibility,
            "regime": self.regime,
            "resolver": self.resolver,
            # --- extra columns to aid human review (not required by spec) ---
            "species_name": self.species_name or "",
            "superkingdom": self.superkingdom or "",
            "total_length_bp": str(self.total_length) if self.total_length is not None else "",
            "n_replicons": str(self.n_replicons),
            "flags": ";".join(self.flags),
            "note": self.note,
        }


TABLE_COLUMNS = [
    "unit_id",
    "scientificName",
    "kalamari_taxid",
    "species_taxid",
    "accession_set",
    "parent_assembly",
    "eligibility",
    "regime",
    "resolver",
    "species_name",
    "superkingdom",
    "total_length_bp",
    "n_replicons",
    "flags",
    "note",
]


# --------------------------------------------------------------------------- #
#  Loaders                                                                     #
# --------------------------------------------------------------------------- #
def load_chromosome_rows(path: str) -> List[ChromosomeRow]:
    rows: List[ChromosomeRow] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            rows.append(
                ChromosomeRow(
                    scientific_name=r["scientificName"].strip(),
                    accession=strip_version(r["nuccoreAcc"]),
                    taxid=int(r["taxid"]),
                    parent=int(r["parent"]),
                    source=r["source"].strip(),
                )
            )
    return rows


def load_confounded_manifest(path: str) -> List[ManifestEntry]:
    entries: List[ManifestEntry] = []
    with open(path, newline="") as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = [p.strip() for p in line.rstrip("\n").split("\t")]
            if parts[0] == "trigger_taxid":  # header
                continue
            if len(parts) < 4:
                continue
            entries.append(
                ManifestEntry(
                    trigger_taxid=int(parts[0]),
                    label=parts[1],
                    regime=parts[2],
                    resolver=parts[3],
                    note=parts[4] if len(parts) > 4 else "",
                )
            )
    return entries


# --------------------------------------------------------------------------- #
#  Grouping                                                                    #
# --------------------------------------------------------------------------- #
def _sanitize(text: str) -> str:
    return text.replace(":", "-").replace("/", "-").replace(" ", "_")


def _bucket_key(rec: Optional[NuccoreRecord], row: ChromosomeRow) -> str:
    """Key that decides which rows share a genome-unit.

    Priority: real GCF/GCA assembly (ground truth -- merges replicons and keeps
    distinct genomes apart) -> a strain proxy SCOPED TO THE KALAMARI IDENTITY ->
    the accession.  The identity-scoped proxy is deliberate: without a real
    assembly, only replicons agreeing on the Kalamari (scientificName, taxid)
    may merge, so distinct genomes that share a species-level organism taxid
    (the 12 Salmonella / 4 Listeria lineages all report their species taxid)
    can NEVER be merged into one chimera by the fallback.
    """
    if rec and is_real_assembly(rec.assembly_acc):
        return rec.assembly_acc  # type: ignore[return-value]
    if rec and rec.strain_taxid:
        return f"strain:{rec.strain_taxid}|{row.scientific_name}|{row.taxid}"
    return f"acc:{row.accession}"


def group_rows_into_units(
    rows: Sequence[ChromosomeRow],
    records: Dict[str, NuccoreRecord],
) -> List[GenomeUnit]:
    """Group chromosome rows into genome-units by assembly (never by taxid)."""
    # Detect the shared-accession quirk: an accession appearing on >1 row.
    acc_row_count: Dict[str, int] = {}
    for r in rows:
        acc_row_count[r.accession] = acc_row_count.get(r.accession, 0) + 1
    shared_accessions = {a for a, n in acc_row_count.items() if n > 1}

    # Bucket rows by their genome-group key.  A shared accession is one physical
    # sequence, so all its rows are forced into one unit (keyed by the accession)
    # regardless of their differing labels -- surfaced, never double-counted.
    buckets: Dict[str, List[ChromosomeRow]] = {}
    for r in rows:
        rec = records.get(r.accession)
        # A shared accession is one physical sequence, so its differently-labelled
        # rows must land in one unit. If the accession has a real assembly, the
        # assembly key already unifies them AND correctly keeps any sibling
        # replicons of that assembly together, so prefer it; only fall back to a
        # dedicated 'sharedacc:' key when there is no assembly to group on (which
        # otherwise the identity-scoped proxy would split by label).
        if r.accession in shared_accessions and not (rec and is_real_assembly(rec.assembly_acc)):
            key = f"sharedacc:{r.accession}"
        else:
            key = _bucket_key(rec, r)
        buckets.setdefault(key, []).append(r)

    units: List[GenomeUnit] = []
    used_ids: Dict[str, int] = {}
    for group_key, group_rows in sorted(buckets.items()):
        accessions = sorted({r.accession for r in group_rows})
        names: List[str] = []
        for r in group_rows:
            if r.scientific_name not in names:
                names.append(r.scientific_name)
        names.sort()
        taxids = sorted({r.taxid for r in group_rows})

        recs = [records[a] for a in accessions if a in records]
        lengths = [rec.length for rec in recs if rec.length is not None]
        total_length = sum(lengths) if lengths else None
        real_assembly = next(
            (rec.assembly_acc for rec in recs if is_real_assembly(rec.assembly_acc)), None
        )
        # Clean display proxy when no real assembly (strain taxid, else acc).
        strain = next((rec.strain_taxid for rec in recs if rec.strain_taxid), None)
        proxy = f"strain:{strain}" if strain else f"acc:{accessions[0]}"
        parent_assembly = real_assembly or proxy

        rep_name = names[0] if names else group_key
        base_id = f"{_sanitize(rep_name)}~{_sanitize(parent_assembly)}"
        n = used_ids.get(base_id, 0)
        used_ids[base_id] = n + 1
        unit_id = base_id if n == 0 else f"{base_id}#{n + 1}"

        flags: List[str] = []
        if any(a in shared_accessions for a in accessions):
            flags.append("shared_accession")
        if len(names) > 1:
            flags.append("label_conflict")
        if not is_real_assembly(real_assembly) and len(accessions) > 1:
            # multi-replicon unit grouped without a real assembly (strain proxy)
            # -- reviewer should confirm these replicons are one genome.
            flags.append("grouped_by_proxy")
        if any(a not in records for a in accessions):
            # a replicon NCBI never returned a record for -- its grouping and
            # organelle status could not be confirmed.
            flags.append("missing_record")

        units.append(
            GenomeUnit(
                unit_id=unit_id,
                group_key=group_key,
                scientific_names=names,
                kalamari_taxids=taxids,
                accession_set=accessions,
                parent_assembly=parent_assembly,
                total_length=total_length,
                n_replicons=len(accessions),
                flags=flags,
            )
        )
    return units


# --------------------------------------------------------------------------- #
#  Species resolution                                                          #
# --------------------------------------------------------------------------- #
def resolve_species(
    unit: GenomeUnit,
    resolver: TaxonomyResolver,
    records: Dict[str, NuccoreRecord],
) -> None:
    """Fill ``species_taxid`` / ``species_name`` / ``superkingdom`` on a unit.

    Primary source is the Kalamari declared taxid climbed to species.  If that
    yields no species (a genus-level 'sp.' entry) we fall back to the genome's
    own organism taxid from NCBI, which is usually a real species.
    """
    # Resolve every declared Kalamari taxid; prefer one that reaches a species.
    resolutions = [resolver.resolve(t) for t in unit.kalamari_taxids]
    species_res = next((r for r in resolutions if r.resolves_to_species), None)
    primary = species_res or resolutions[0]

    unit.superkingdom = next(
        (r.superkingdom for r in resolutions if r.superkingdom), primary.superkingdom
    )

    # Independently resolve the genome's OWN NCBI organism to a species -- this is
    # the ground truth NCBI records for the deposited sequence.
    organism_species = set()
    for acc in unit.accession_set:
        rec = records.get(acc)
        if rec and rec.strain_taxid:
            org = resolver.resolve(rec.strain_taxid)
            if org.resolves_to_species:
                organism_species.add((org.species_taxid, org.species_name))

    if species_res is not None:
        unit.species_taxid = species_res.species_taxid
        unit.species_name = species_res.species_name
        # Surface a mislabel: the Kalamari taxid climbs to a DIFFERENT species
        # than the deposited genome's own NCBI organism (e.g. the Kalamari row
        # for Yersinia_kristensenii carries taxid 631 = Y. intermedia).  This is
        # exactly the kind of drift the triage tool exists to flag, so keep the
        # declared species (faithful policy snapshot) but flag + note the NCBI
        # organism for SME review.
        # Only meaningful for the genomes that feed the ANI analysis (bacteria/
        # archaea); a host mitogenome's declared-vs-organism species is not a
        # curation signal.  The unit is clean only when the declared label(s) and
        # the genome organism(s) agree on ONE species; anything else -- the taxid
        # points at a different species, the labels disagree with each other, or
        # replicons report different organisms (a chimera hint) -- is flagged.
        org_taxids = {t for t, _ in organism_species}
        declared_species = {r.species_taxid for r in resolutions if r.resolves_to_species}
        clean = (
            len(declared_species) == 1
            and bool(org_taxids)
            and org_taxids == declared_species
        )
        mismatch = bool(org_taxids) and not clean
        if mismatch and unit.superkingdom in ("Bacteria", "Archaea"):
            unit.flags.append("taxid_organism_mismatch")
            org_str = ", ".join(sorted(f"{t} ({n})" for t, n in organism_species))
            declared_str = ", ".join(str(t) for t in sorted(declared_species))
            unit.note = (unit.note + "; " if unit.note else "") + (
                f"declared species {{{declared_str}}} "
                f"({species_res.species_name}) but genome organism is {org_str}"
            )
        return

    # Kalamari taxid did not reach a species — try the genome's organism taxid.
    for acc in unit.accession_set:
        rec = records.get(acc)
        if not rec or not rec.strain_taxid:
            continue
        org = resolver.resolve(rec.strain_taxid)
        if org.resolves_to_species:
            unit.species_taxid = org.species_taxid
            unit.species_name = org.species_name
            if org.superkingdom:
                unit.superkingdom = org.superkingdom
            unit.flags.append("species_from_organism")
            unit.note = (
                f"declared taxid {primary.input_taxid} is "
                f"{primary.input_rank or 'unresolvable'}; species taken from "
                f"genome organism {rec.strain_taxid}"
            )
            return

    # Nothing reached a species. Decide hard-fail vs. soft-review by rank:
    # demote to review:genus-level ONLY if EVERY declared label is genuinely
    # genus-rank-or-above (a legitimate higher-level entry).  If any label is
    # species/strain/subspecies -- or the ambiguous "no rank"/unknown -- it
    # "should have" resolved, so it stays eligible and the preflight hard-fails.
    unit._primary_rank = primary.input_rank  # type: ignore[attr-defined]
    unit._demote_to_genus = all(  # type: ignore[attr-defined]
        is_genus_or_above(r.input_rank) for r in resolutions
    )


# --------------------------------------------------------------------------- #
#  Eligibility                                                                 #
# --------------------------------------------------------------------------- #
def _is_non_nuclear(rec: NuccoreRecord) -> bool:
    """A replicon that is an organelle or a genome segment (not a nuclear/main
    chromosome).  A unit made ENTIRELY of these is an organelle/segment-only
    submission; a unit that also has real chromosomes is not."""
    if rec.is_organelle:
        return True
    if rec.moltype and any(h in rec.moltype.lower() for h in SEGMENT_MOLTYPE_HINTS):
        return True
    return False


def classify_eligibility(
    unit: GenomeUnit,
    records: Dict[str, NuccoreRecord],
    size_cap_bp: int = DEFAULT_SIZE_CAP_BP,
) -> None:
    """Programmatic eligibility guard (order matters; see module docstring)."""
    recs = [records[a] for a in unit.accession_set if a in records]
    sk = unit.superkingdom

    # 1. Viruses first -- a segmented virus (Andes hantavirus) is still a virus,
    #    not an "organelle/segment" exclusion.
    if sk in ("Viruses", "Viroids"):
        unit.eligibility = ELIG_VIRUS
        return

    # 2. Organelle-ONLY submission -> exclude:organelle (catches the Arripis
    #    fish mitogenome and every host mito/plastid).  Crucially this fires only
    #    when EVERY replicon is non-nuclear, so a eukaryote nuclear assembly that
    #    merely bundles its mitochondrion (Cryptococcus: 14 chromosomes + 1 mito)
    #    is NOT excluded here -- it falls through to the eukaryote track.  Require
    #    a record for EVERY replicon first, so a missing nuclear-chromosome record
    #    can't make a real assembly look organelle-only.
    all_present = len(recs) == len(unit.accession_set)
    non_nuclear = [rec for rec in recs if _is_non_nuclear(rec)]
    if recs and all_present and len(non_nuclear) == len(recs):
        unit.eligibility = ELIG_ORGANELLE
        kinds = sorted({(rec.genome or rec.moltype or "?") for rec in non_nuclear})
        unit.note = (unit.note + "; " if unit.note else "") + f"molecule={','.join(kinds)}"
        return

    # 3. Oversized grouped assembly -> exclude:organelle (huge eukaryote nuclear).
    #    Gated to NON-bacterial/archaeal domains: a bacterium/archaeon is never
    #    >50 Mbp, so an oversized prokaryotic unit is a data problem to keep
    #    eligible (and let the preflight/summary surface), not silently relabel
    #    as an organelle.
    if (
        unit.total_length is not None
        and unit.total_length > size_cap_bp
        and sk not in ("Bacteria", "Archaea")
    ):
        unit.eligibility = ELIG_ORGANELLE
        unit.note = (unit.note + "; " if unit.note else "") + f"length>{size_cap_bp}bp"
        return

    # 4. Classify by superkingdom of the resolved organism.
    if sk in ("Bacteria", "Archaea"):
        if unit.species_taxid is not None:
            unit.eligibility = ELIG_ELIGIBLE
        else:
            # Bacteria/archaea with no species: a genuine genus-or-above entry
            # soft-demotes to review:genus-level; anything more specific (or an
            # ambiguous "no rank") stays 'eligible' so the preflight HARD-FAILS.
            if getattr(unit, "_demote_to_genus", False):
                unit.eligibility = ELIG_GENUS
                unit.note = (unit.note + "; " if unit.note else "") + "no species-rank ancestor"
            else:
                unit.eligibility = ELIG_ELIGIBLE  # -> caught by preflight
    elif sk == "Eukaryota":
        unit.eligibility = ELIG_EUKARYOTE
    else:
        unit.eligibility = ELIG_UNRESOLVED


# --------------------------------------------------------------------------- #
#  Confounded tagging                                                          #
# --------------------------------------------------------------------------- #
def tag_confounded(
    unit: GenomeUnit,
    manifest: Sequence[ManifestEntry],
    resolver: TaxonomyResolver,
) -> None:
    """Set ``regime``/``resolver`` if the unit's lineage hits a manifest trigger."""
    lineage_taxids = set()
    if unit.species_taxid:
        lineage_taxids.update(resolver.lineage(unit.species_taxid))
    for t in unit.kalamari_taxids:
        lineage_taxids.update(resolver.lineage(t))

    for entry in manifest:
        if entry.trigger_taxid in lineage_taxids:
            unit.regime = entry.regime
            unit.resolver = entry.resolver
            if "confounded" not in unit.flags:
                unit.flags.append("confounded")
            tag = f"confounded:{entry.label}"
            unit.note = (unit.note + "; " if unit.note else "") + tag
            return


# --------------------------------------------------------------------------- #
#  Orchestration + preflight                                                   #
# --------------------------------------------------------------------------- #
def build_units(
    rows: Sequence[ChromosomeRow],
    resolver: TaxonomyResolver,
    records: Dict[str, NuccoreRecord],
    manifest: Sequence[ManifestEntry],
    size_cap_bp: int = DEFAULT_SIZE_CAP_BP,
) -> List[GenomeUnit]:
    """Run the full Tier 0 pipeline and return the genome-units, sorted."""
    units = group_rows_into_units(rows, records)
    for unit in units:
        resolve_species(unit, resolver, records)
        classify_eligibility(unit, records, size_cap_bp=size_cap_bp)
        tag_confounded(unit, manifest, resolver)
    units.sort(key=lambda u: (u.scientific_name_field().lower(), u.group_key))
    return units


class PreflightError(RuntimeError):
    """Raised when an eligible unit fails to resolve to a species taxon."""


def run_preflight(units: Sequence[GenomeUnit]) -> List[GenomeUnit]:
    """Return the offending units (eligible but no species). Empty = pass."""
    return [u for u in units if u.eligibility == ELIG_ELIGIBLE and u.species_taxid is None]


def assert_preflight(units: Sequence[GenomeUnit]) -> None:
    """Raise :class:`PreflightError` if any eligible unit lacks a species taxid."""
    offenders = run_preflight(units)
    if offenders:
        lines = [
            f"  {u.unit_id}: {u.scientific_name_field()} "
            f"(kalamari_taxid={u.kalamari_taxid_field()}) did not resolve to a "
            f"species-or-below NCBI taxon"
            for u in offenders
        ]
        raise PreflightError(
            f"Preflight FAILED: {len(offenders)} eligible unit(s) have no real "
            f"species taxid:\n" + "\n".join(lines)
        )


def write_table(units: Sequence[GenomeUnit], path: str) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TABLE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for u in units:
            writer.writerow(u.as_row())
