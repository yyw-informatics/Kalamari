#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_type_strains.py -- Tier 2: resolve the type-strain assembly for every pick.

Tier 2's ANI check asks "how far is our pick from the type strain of the species it
claims to be?", which first needs to know *which assembly* that type strain is.

The obvious way is a per-species ``datasets summary`` query.  The cheap way -- used
here -- is to notice that NCBI's ``ANI_report_prokaryotes.txt`` already names it: every
row carries ``declared-type-assembly`` (the type strain of the species the assembly
declares) alongside NCBI's own ANI and coverages against it.  Tier 1 already streams
that file, so resolving 249 type strains costs one more pass and no new dependency.

NCBI's numbers are kept next to ours deliberately.  Tier 2 computes its own ANI with
skani; having the source's value in the same record makes the two directly comparable,
so a disagreement shows up as a disagreement instead of hiding.

    python3 bin/curation-triage/build_type_strains.py            # stream from NCBI
    python3 bin/curation-triage/build_type_strains.py --ani-file ANI_report.txt
"""

from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ani_reports  # noqa: E402
import coverage  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")
DEFAULT_POLICY = os.path.join(SRC, "taxon_policy.tsv")
DEFAULT_CHROMOSOMES = os.path.join(REPO_ROOT, "src", "chromosomes.tsv")
DEFAULT_OUT = os.path.join(SRC, "type_strain_cache.json")


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def build(picks, lines) -> dict:
    """Map each pick's own assembly -> its declared type-strain assembly.

    Also keeps a species-taxid map as a fallback, but the per-assembly one is
    authoritative: a species with subspecies has a different declared type strain per
    subspecies, and eleven of our units share one species taxid.
    """
    wanted = {p.base for p in picks if p.base.startswith(("GCF_", "GCA_"))}
    by_base = ani_reports.distill_type_strains(lines, wanted)
    log("ANI report: type-strain rows for %d of %d picks" % (len(by_base), len(wanted)))

    by_assembly, by_species, no_type, self_type, subspecies_split = {}, {}, [], [], {}
    for pick in sorted(picks, key=lambda p: p.unit_id):
        rec = by_base.get(pick.base)
        if not rec:
            continue
        acc = rec.get("declared_type_assembly", "")
        # The column is not always an accession, and the two non-accession cases mean
        # opposite things:
        #   'same'    -> THIS assembly IS the declared type strain.  That is the
        #                strongest identity confirmation available, not a missing
        #                value, so it is recorded as such and skips the comparison.
        #   anything  -> NCBI has no type material for the declared species; there is
        #   else        nothing to compare against and Tier 2 must say so.
        # Handing either string to the downloader as if it were an accession would be
        # the obvious bug here.
        if acc.strip().lower() == "same":
            self_type.append(pick.unit_id)
            by_assembly[pick.base] = {
                "accession": pick.parent_assembly,
                "organism": rec.get("species_name", ""),
                "category": "self",
                "is_self_type": True,
                "source": "NCBI ANI_report_prokaryotes.txt declared-type-assembly='same'",
                "note": "this pick IS the declared type strain of its species",
            }
            continue
        if not acc.startswith(("GCA_", "GCF_")):
            no_type.append("%s (%s)" % (pick.unit_id, acc or "empty"))
            # Recorded rather than dropped.  "NCBI holds no type material for this
            # species" is an ANSWER, not a missing value -- and it is the reason
            # NCBI's own taxonomy check comes back Inconclusive for several of the
            # picks Tier 1 flagged.  Saying so turns an unexplained flag into an
            # explained one.
            by_assembly[pick.base] = {
                "accession": "",
                "organism": rec.get("species_name", ""),
                "category": acc or "none",
                "is_self_type": False,
                "has_no_type_material": True,
                "source": "NCBI ANI_report_prokaryotes.txt declared-type-assembly",
                "note": "NCBI holds no type material for this species, so no "
                        "type-strain comparison exists",
            }
            continue
        entry = {
            "accession": acc,
            "is_self_type": False,
            "organism": rec.get("declared_type_organism_name", ""),
            "category": rec.get("declared_type_category", ""),
            "source": "NCBI ANI_report_prokaryotes.txt declared-type-assembly",
            "ncbi_declared_type_ani": rec.get("declared_type_ani", ""),
            "ncbi_declared_type_qcoverage": rec.get("declared_type_qcoverage", ""),
            "ncbi_declared_type_scoverage": rec.get("declared_type_scoverage", ""),
            "best_match_type_assembly": rec.get("best_match_type_assembly", ""),
        }
        # Keyed by the pick's OWN assembly. NCBI names a declared type strain per
        # assembly, and for a species with subspecies those differ -- all eleven
        # Kalamari Salmonella units declare species 28901 but get eleven different,
        # subspecies-appropriate type strains. Keying only by species would compare
        # ten of them against the wrong genome.
        by_assembly[pick.base] = entry
        taxid = pick.species_taxid
        if taxid is None:
            continue
        prior = by_species.get(taxid)
        if prior and prior["accession"] != acc:
            subspecies_split.setdefault(taxid, set()).update([prior["accession"], acc])
        by_species.setdefault(taxid, entry)
    return {"by_assembly": by_assembly, "by_species": by_species, "no_type": no_type,
            "self_type": self_type, "subspecies_split": {k: sorted(v) for k, v in subspecies_split.items()}}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Resolve Tier-2 type-strain assemblies.")
    p.add_argument("--policy", default=DEFAULT_POLICY)
    p.add_argument("--chromosomes", default=DEFAULT_CHROMOSOMES)
    p.add_argument("--out", default=DEFAULT_OUT)
    p.add_argument("--ani-file", default=None,
                   help="A local copy of ANI_report_prokaryotes.txt (default: stream it).")
    args = p.parse_args(argv)

    lanes = coverage.load_source_lanes(args.chromosomes)
    picks = coverage.load_picks(args.policy, lanes)
    log("Loaded %d eligible picks" % len(picks))

    lines = (ani_reports.open_file_lines(args.ani_file) if args.ani_file
             else ani_reports.open_url_lines(ani_reports.ANI_REPORT_URL))
    result = build(picks, lines)
    by_species, by_assembly = result["by_species"], result["by_assembly"]

    payload = {
        "schema": 1,
        "meta": {
            "status": "populated",
            "source": "NCBI ASSEMBLY_REPORTS/ANI_report_prokaryotes.txt",
            "column": "declared-type-assembly (the type strain of the DECLARED species)",
            "n_assemblies": len(by_assembly),
            "n_species_fallback": len(by_species),
            "key": "by_assembly is authoritative; by_species is a fallback for a pick "
                   "with no row of its own",
            "picks_without_a_declared_type": len(result["no_type"]),
            "picks_that_are_themselves_the_type_strain": len(result["self_type"]),
            "note": "NCBI's own ANI and coverages against the same type strain are kept "
                    "alongside each entry so Tier 2's independent skani numbers can be "
                    "compared against the source they confirm.",
        },
        "by_assembly": dict(sorted(by_assembly.items())),
        "by_species": {str(k): v for k, v in sorted(by_species.items())},
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    tmp = "%s.tmp.%d" % (args.out, os.getpid())
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, args.out)

    log("Wrote %s: %d assemblies (%d species fallback) with a type strain"
        % (args.out, len(by_assembly), len(by_species)))
    if result["self_type"]:
        log("%d pick(s) ARE the declared type strain of their species (e.g. %s)"
            % (len(result["self_type"]), ", ".join(result["self_type"][:3])))
    if result["no_type"]:
        log("%d pick(s) have no declared type assembly (NCBI has none for that species): %s"
            % (len(result["no_type"]), ", ".join(result["no_type"][:5])))
    for taxid, accs in sorted(result["subspecies_split"].items()):
        log("species %s spans %d declared type assemblies (%s) -- expected where a "
            "species has subspecies; the per-assembly key handles it"
            % (taxid, len(accs), ", ".join(accs[:4])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
