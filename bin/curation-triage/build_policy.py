#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_policy.py -- Tier 0 of the Kalamari curation-triage enhancement.

Reads ``src/chromosomes.tsv`` plus NCBI taxonomy and emits a committed,
human-reviewable ``src/curation-triage/taxon_policy.tsv`` -- one row per
**genome-unit** (replicons grouped by assembly, never by the Kalamari taxid
column).  Every downstream tier keys off this table.

Output columns
--------------
Required (docs/CURATION_TRIAGE_DESIGN.md section 4):

    unit_id          stable, readable id: "<repName>~<groupKey>"
    scientificName   Kalamari name(s) for the unit (';'-joined if >1)
    kalamari_taxid   Kalamari-declared taxid(s) (';'-joined if >1)
    species_taxid    resolved real NCBI species-rank taxid (blank if none)
    accession_set    ','-joined nuccore accessions in this genome-unit
    parent_assembly  real GCF/GCA if resolved, else the "strain:<taxid>" proxy
    eligibility      eligible | exclude:organelle | exclude:virus | eukaryote
                     | review:genus-level | review:unresolved
    regime           confounded regime A/B (blank if not confounded)
    resolver         marker | ANI (blank if not confounded)

Extra review columns follow: species_name, superkingdom, total_length_bp,
n_replicons, flags, note.

How it runs
-----------
Grouping and the organelle guard need three NCBI facts per accession
(parent genome, molecule type, length).  These come from one batched nuccore
``esummary`` call and are cached to ``src/curation-triage/nuccore_cache.json``.
Once that cache is committed, the build (and every unit test) is fully offline.

    # first build (populates + commits the cache):
    python3 bin/curation-triage/build_policy.py

    # offline rebuild from the committed cache (fails if anything is uncached):
    python3 bin/curation-triage/build_policy.py --offline

Taxonomy comes from a full NCBI taxdump, with the bundled ``src/taxonomy/build``
overlay layered on top so the synthetic Kalamari lineages climb to their real
parent species.  Point it at one with ``--taxdump DIR``, ``$NCBI_TAXDUMP``,
``$NCBI_TAXDUMP_SEARCH``, or a gitignored ``.taxdump-search`` file at the repo
root (one directory per line) so a machine records its own layout once and the
repository carries nobody's absolute paths.

Exit status is non-zero if the preflight finds an eligible unit with no real
species taxid (a synthetic/strain/genus key must never silently pass).
"""

from __future__ import annotations

import argparse
import collections
import os
import sys

# Allow "import taxonomy/ncbi/policy" whether run as a script or a module.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from taxonomy import TaxonomyResolver  # noqa: E402
from ncbi import NuccoreResolver, is_real_assembly  # noqa: E402
import policy  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEFAULT_CHROMOSOMES = os.path.join(REPO_ROOT, "src", "chromosomes.tsv")
DEFAULT_OUT = os.path.join(REPO_ROOT, "src", "curation-triage", "taxon_policy.tsv")
DEFAULT_CACHE = os.path.join(REPO_ROOT, "src", "curation-triage", "nuccore_cache.json")
DEFAULT_MANIFEST = os.path.join(REPO_ROOT, "src", "curation-triage", "confounded_manifest.tsv")
BUNDLED_OVERLAY = os.path.join(REPO_ROOT, "src", "taxonomy", "build")
BUNDLED_FALLBACK = os.path.join(REPO_ROOT, "src", "taxonomy", "deprecated")


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def resolve_taxonomy(args) -> TaxonomyResolver:
    base = TaxonomyResolver.find_base_taxdump(args.taxdump)
    overlays = [BUNDLED_OVERLAY]
    if base is None:
        log(
            "WARNING: no full NCBI taxdump found (--taxdump / $NCBI_TAXDUMP / "
            "$NCBI_TAXDUMP_SEARCH / .taxdump-search). Falling back to the bundled "
            "partial 'deprecated' dump; some taxids may fail to resolve."
        )
        base = BUNDLED_FALLBACK
    log(f"Taxonomy base dump: {base}")
    log(f"Taxonomy overlay(s): {overlays}")
    return TaxonomyResolver.from_dumps(base, overlay_dirs=overlays)


def print_summary(units, offenders) -> None:
    elig_counts = collections.Counter(u.eligibility for u in units)
    confounded = [u for u in units if u.regime]
    shared = [u for u in units if "shared_accession" in u.flags]
    organism_species = [u for u in units if "species_from_organism" in u.flags]
    real_asm = [u for u in units if is_real_assembly(u.parent_assembly)]
    proxy_multi = [u for u in units if "grouped_by_proxy" in u.flags]
    mislabels = [u for u in units if "taxid_organism_mismatch" in u.flags]

    print()
    print("=" * 72)
    print("Tier 0 policy table summary")
    print("=" * 72)
    print(f"Total genome-units: {len(units)}")
    print(f"Grouped by real GCF/GCA assembly: {len(real_asm)}/{len(units)}"
          f"  (multi-replicon units still on a strain-taxid proxy: {len(proxy_multi)})")
    print()
    print("Eligibility classes:")
    for cls, n in sorted(elig_counts.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"  {cls:<22} {n}")
    print()
    print(f"Confounded units (regime tagged): {len(confounded)}")
    by_regime = collections.Counter((u.regime, u.resolver) for u in confounded)
    for (regime, resolver), n in sorted(by_regime.items()):
        print(f"  regime {regime} / resolver {resolver:<7} {n}")
    for u in sorted(confounded, key=lambda u: u.scientific_name_field().lower()):
        print(f"    - {u.scientific_name_field()} "
              f"(species_taxid={u.species_taxid}, regime {u.regime}/{u.resolver})")
    print()
    if shared:
        print(f"Shared-accession quirks flagged: {len(shared)}")
        for u in shared:
            print(f"  - {u.accession_field()} shared by {u.scientific_name_field()} "
                  f"(kalamari_taxid={u.kalamari_taxid_field()})")
        print()
    if organism_species:
        print(f"Species taken from genome organism (genus-level declared taxid): "
              f"{len(organism_species)}")
        for u in organism_species:
            print(f"  - {u.scientific_name_field()} -> species_taxid {u.species_taxid} "
                  f"({u.species_name})")
        print()
    if mislabels:
        print(f"Taxid/organism species MISMATCH (declared taxid disagrees with the "
              f"genome's NCBI organism -- flag for SME review): {len(mislabels)}")
        for u in mislabels:
            print(f"  - {u.scientific_name_field()}: {u.note}")
        print()
    print(f"Preflight: {'FAILED' if offenders else 'passed'} "
          f"({len(offenders)} eligible unit(s) without a species taxid)")
    print("=" * 72)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Build the Tier 0 taxon policy table.")
    p.add_argument("--chromosomes", default=DEFAULT_CHROMOSOMES)
    p.add_argument("--out", default=DEFAULT_OUT)
    p.add_argument("--cache", default=DEFAULT_CACHE)
    p.add_argument("--manifest", default=DEFAULT_MANIFEST)
    p.add_argument("--taxdump", default=None,
                   help="Full NCBI taxdump dir (nodes/names/merged.dmp).")
    p.add_argument("--offline", action="store_true",
                   help="Never hit the network; fail if the cache is incomplete.")
    p.add_argument("--no-assemblies", action="store_true",
                   help="Skip elink nuccore->assembly resolution (grouping then "
                        "falls back to the strain-taxid proxy; less accurate).")
    p.add_argument("--size-cap-bp", type=int, default=policy.DEFAULT_SIZE_CAP_BP)
    p.add_argument("--no-write", action="store_true",
                   help="Compute + preflight + summarize but do not write the table.")
    args = p.parse_args(argv)

    rows = policy.load_chromosome_rows(args.chromosomes)
    log(f"Loaded {len(rows)} chromosome rows from {args.chromosomes}")

    manifest = policy.load_confounded_manifest(args.manifest)
    log(f"Loaded {len(manifest)} confounded-manifest trigger(s)")

    resolver = resolve_taxonomy(args)

    nuccore = NuccoreResolver(cache_path=args.cache)
    accessions = [r.accession for r in rows]
    records = nuccore.resolve_all(accessions, allow_fetch=not args.offline, log=log)
    if not args.offline:
        nuccore.save()  # persist esummary facts before the slower elink step
    if not args.offline and not args.no_assemblies:
        nuccore.resolve_assemblies(records, log=log)
        nuccore.save()
    if not args.offline:
        log(f"Wrote nuccore cache: {args.cache} ({len(nuccore.records)} records)")

    missing = [a for a in accessions if a not in records]
    if missing:
        log(f"WARNING: {len(missing)} accession(s) had no NCBI record: {missing}")

    units = policy.build_units(rows, resolver, records, manifest, size_cap_bp=args.size_cap_bp)
    offenders = policy.run_preflight(units)

    if not args.no_write:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        policy.write_table(units, args.out)
        log(f"Wrote policy table: {args.out} ({len(units)} units)")

    print_summary(units, offenders)

    if offenders:
        log("")
        policy_error = policy.PreflightError(
            f"{len(offenders)} eligible unit(s) without a species taxid"
        )
        log(f"PREFLIGHT FAILURE: {policy_error}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
