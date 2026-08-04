#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
coverage_probe.py -- the one-shot coverage probe under Tier 1.

The **make-or-break measurement**.  Tier 1 rests on the "pick-health sentinel"
(J1a), and this is what shows the signal is really there: it joins the
eligible Tier-0 picks against NCBI's own precomputed verdicts and reports what
fraction have a *live row for OUR assembly* -- separately in the ANI report and in
assembly_summary -- plus the breakdowns and a J1a flag preview that show the real
precision, not just coverage.

    Inputs
      src/curation-triage/taxon_policy.tsv   (Tier-0 output; the join set = eligible units)
      src/chromosomes.tsv                   (source lane per accession)
      NCBI ASSEMBLY_REPORTS/{ANI_report_prokaryotes, assembly_summary_refseq,
                             assembly_summary_genbank}.txt  (streamed + distilled to cache)

    Outputs
      src/curation-triage/coverage_report.tsv       (one row per eligible unit)
      src/curation-triage/assembly_reports_cache.json (committed distilled cache)
      a printed summary: headline coverage + breakdowns + J1a flag preview

How it runs
-----------
    # first run: streams the three large NCBI files (~2.7 GB total), distils only
    # the rows touching our picks + per-species presence counts, commits the cache:
    python3 bin/curation-triage/coverage_probe.py

    # offline rerun from the committed cache (no network; fails if cache incomplete):
    python3 bin/curation-triage/coverage_probe.py --offline

The network is quarantined in ``ani_reports.py`` (streaming distillers + injectable
openers); ``coverage.py`` holds the pure join / match-typing / flag / summary logic,
so the whole thing is unit-testable offline -- exactly like Tier 0 split ``ncbi.py``
from ``policy.py``.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ani_reports  # noqa: E402
import coverage  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEFAULT_POLICY = os.path.join(REPO_ROOT, "src", "curation-triage", "taxon_policy.tsv")
DEFAULT_CHROMOSOMES = os.path.join(REPO_ROOT, "src", "chromosomes.tsv")
DEFAULT_OUT = os.path.join(REPO_ROOT, "src", "curation-triage", "coverage_report.tsv")
DEFAULT_CACHE = os.path.join(REPO_ROOT, "src", "curation-triage", "assembly_reports_cache.json")


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def _pct(num: int, den: int) -> str:
    return f"{(100.0 * num / den):.1f}%" if den else "n/a"


def _print_rate_table(title: str, table: dict, denom_key: str = "n") -> None:
    print(f"  {title}")
    for k in sorted(table, key=lambda x: (-table[x][denom_key], x)):
        d = table[k]
        print(f"    {k:<14} n={d['n']:<4} "
              f"ANI-live {d['ani_live']}/{d['n']} ({_pct(d['ani_live'], d['n'])})   "
              f"summary-live {d['summary_live']}/{d['n']} ({_pct(d['summary_live'], d['n'])})   "
              f"in-RefSeq {d['in_refseq']}/{d['n']}")


def print_summary(rows, s: dict) -> None:
    n = s["n_picks"]
    h = s["headline"]
    print()
    print("=" * 78)
    print("J1a coverage probe: does NCBI's verdict cover our picks?")
    print("=" * 78)
    print(f"Eligible genome-units probed: {n}")
    print()
    print("HEADLINE COVERAGE (live row for OUR assembly)")
    print(f"  ANI report        : {h['ani_live']}/{n}  ({h['ani_live_pct']}%)")
    print(f"  assembly_summary  : {h['summary_live']}/{n}  ({h['summary_live_pct']}%)")
    print(f"  live in RefSeq now: {h['in_refseq']}/{n}  ({_pct(h['in_refseq'], n)})")
    print()
    print("ANI-report match type:")
    for m, c in sorted(s["ani_match_types"].items(), key=lambda kv: -kv[1]):
        print(f"    {m:<16} {c}")
    print("assembly_summary match type:")
    for m, c in sorted(s["summary_match_types"].items(), key=lambda kv: -kv[1]):
        print(f"    {m:<16} {c}")
    print()
    print("BREAKDOWNS")
    _print_rate_table("by source lane:", s["by_source_lane"])
    _print_rate_table("by assembly kind (GCF/GCA):", s["by_assembly_kind"])
    _print_rate_table("by assembly release year:", s["by_year"])
    print()
    print("J1a FLAG PREVIEW (what the pick-health sentinel would surface)")
    print(f"  Units with >=1 actionable flag: {s['actionable_units']}/{n}  "
          f"(confounded units, ANI-ambiguity expected & marker-resolved: {s['confounded_units']})")
    if s["actionable_by_lane"]:
        lanes = ", ".join(f"{k}:{v}" for k, v in sorted(s["actionable_by_lane"].items()))
        print(f"  actionable units by source lane: {lanes}")
    print("  flag counts (all reasons; * = informational/context, not an SME action item):")
    for f, c in sorted(s["flag_counts"].items(), key=lambda kv: (-kv[1], kv[0])):
        mark = "" if f in coverage.ACTIONABLE_FLAGS else " *"
        print(f"    {f:<28} {c}{mark}")
    print()
    print("  Actionable picks (unit -> flags):")
    shown = 0
    for r in sorted(rows, key=lambda r: r.pick.scientific_name.lower()):
        if r.actionable_flags:
            shown += 1
            best = f" [best-match: {r.best_match_species_name}]" if r.best_match_species_name and \
                "best_match_species_mismatch" in r.flags else ""
            excl = f" [excluded: {r.excluded_from_refseq}]" if r.excluded_from_refseq else ""
            print(f"    - {r.pick.scientific_name} ({r.pick.parent_assembly}, {r.pick.source_lane}): "
                  f"{';'.join(r.actionable_flags)}{best}{excl}")
    if not shown:
        print("    (none)")
    print("=" * 78)


def build_or_load_cache(args, picks) -> ani_reports.AssemblyReportsCache:
    # Only real GCF/GCA bases can join an ASSEMBLY_REPORTS row; a 'strain:'/'acc:'
    # Tier-0 proxy (unresolved assembly) is kept in the report but out of the stream
    # set (it is surfaced via the 'unresolved_assembly' flag instead).
    wanted_bases = {p.base for p in picks if coverage.is_real_assembly_base(p.base)}
    wanted_species = {p.species_taxid for p in picks if p.species_taxid}
    unresolved = [p for p in picks if not coverage.is_real_assembly_base(p.base)]
    if unresolved:
        log(f"WARNING: {len(unresolved)} eligible pick(s) have no real GCF/GCA parent "
            f"assembly; joined as 'unresolved_assembly': {[p.unit_id for p in unresolved][:5]}")

    openers = {}
    if args.ani_file:
        openers["ani_opener"] = lambda: ani_reports.open_file_lines(args.ani_file)
    if args.refseq_file:
        openers["refseq_opener"] = lambda: ani_reports.open_file_lines(args.refseq_file)
    if args.genbank_file:
        openers["genbank_opener"] = lambda: ani_reports.open_file_lines(args.genbank_file)

    cache = ani_reports.AssemblyReportsCache(cache_path=args.cache, **openers)

    if args.offline:
        if not cache.covers(wanted_bases, wanted_species):
            missing = sorted(wanted_bases - cache.wanted_bases)[:5]
            missing_sp = sorted(wanted_species - cache.wanted_species)[:5]
            raise SystemExit(
                f"--offline: committed cache {args.cache} does not cover all "
                f"{len(wanted_bases)} picks / {len(wanted_species)} species "
                f"(e.g. missing bases {missing}, species {missing_sp}). "
                f"Rebuild without --offline."
            )
        log(f"Offline: loaded distilled cache for {len(cache.wanted_bases)} picks from {args.cache}")
        return cache

    log(f"Building assembly-reports cache for {len(wanted_bases)} picks "
        f"({len(wanted_species)} species) ...")
    cache.build(wanted_bases, wanted_species,
                include_genbank=not args.no_genbank, log=log)
    cache.save()
    log(f"Wrote distilled cache: {args.cache}")
    return cache


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="J1a coverage probe: does NCBI's verdict cover our picks?")
    p.add_argument("--policy", default=DEFAULT_POLICY)
    p.add_argument("--chromosomes", default=DEFAULT_CHROMOSOMES)
    p.add_argument("--out", default=DEFAULT_OUT)
    p.add_argument("--cache", default=DEFAULT_CACHE)
    p.add_argument("--offline", action="store_true",
                   help="Never hit the network; fail if the committed cache is incomplete.")
    p.add_argument("--no-genbank", action="store_true",
                   help="Skip the 1.6 GB genbank summary (loses GenBank-only picks + "
                        "the withdrawal reason for RefSeq-dropped picks).")
    p.add_argument("--no-write", action="store_true",
                   help="Compute + summarize but do not write the report table.")
    # Local pre-downloaded copies (dev / CI staging), bypassing the URL stream.
    p.add_argument("--ani-file", default=None, help="Local ANI_report_prokaryotes.txt.")
    p.add_argument("--refseq-file", default=None, help="Local assembly_summary_refseq.txt.")
    p.add_argument("--genbank-file", default=None, help="Local assembly_summary_genbank.txt.")
    args = p.parse_args(argv)

    source_lanes = coverage.load_source_lanes(args.chromosomes)
    picks = coverage.load_picks(args.policy, source_lanes)
    log(f"Loaded {len(picks)} eligible picks from {args.policy}")

    cache = build_or_load_cache(args, picks)

    rows = coverage.build_coverage(picks, cache)
    summary = coverage.summarize(rows)

    if not args.no_write:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        coverage.write_report(rows, args.out)
        log(f"Wrote coverage report: {args.out} ({len(rows)} rows)")

    print_summary(rows, summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
