#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tier1_metadata.py -- the monthly, delta-aware Tier 1 curation-triage net.

Turns the one-shot coverage probe into the recurring SME handoff loop:

    1. J1a  pick-health -- reuse coverage.py unchanged over the committed cache,
            then DIFF against the last ledger run so only *newly* actionable picks
            surface (not the standing set).
    2. J1b  source-gated better-rep feed   [collect_j1b]
    3. J2   coverage feeders               [collect_j2]
    4. Ledger + reviews + worklist: append one self-contained NDJSON record per lead
            per run to ledger.ndjson, honour reviews.tsv (accept/reject/defer) to
            suppress dismissed leads, and emit the surviving actionable rows as a TSV
            artifact + a $GITHUB_STEP_SUMMARY markdown block.

Informational, never gating.  Precision over recall.  The network is quarantined in
ani_reports.py (+ datasets_summary.py / lpsn.py for the feeders), so every run and
every test is fully offline from the committed caches -- exactly like Tier 0 and the
coverage probe.

    # offline monthly run from the committed cache (no network):
    python3 bin/curation-triage/tier1_metadata.py --offline --run-id 2026-08

    # refresh the cache from NCBI, then run (first run / cache rebuild):
    python3 bin/curation-triage/tier1_metadata.py --run-id 2026-08
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import coverage  # noqa: E402
import coverage_probe  # noqa: E402  (reuse its build_or_load_cache network helper)
import leads as leads_mod  # noqa: E402
import ledger as ledger_mod  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")
DEFAULT_POLICY = os.path.join(SRC, "taxon_policy.tsv")
DEFAULT_CHROMOSOMES = os.path.join(REPO_ROOT, "src", "chromosomes.tsv")
DEFAULT_CACHE = os.path.join(SRC, "assembly_reports_cache.json")
DEFAULT_LEDGER = os.path.join(SRC, "ledger.ndjson")
DEFAULT_REVIEWS = os.path.join(SRC, "reviews.tsv")
DEFAULT_WORKLIST = os.path.join(SRC, "worklist.tsv")
DEFAULT_DATASETS_CACHE = os.path.join(SRC, "datasets_summary_cache.json")
DEFAULT_TODO = os.path.join(REPO_ROOT, "src", "chromosomes-todo.tsv")
DEFAULT_LPSN_CACHE = os.path.join(SRC, "lpsn_cache.json")


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


# --------------------------------------------------------------------------- #
#  Lead collection (J1a now; J1b / J2 hook in here in later steps)              #
# --------------------------------------------------------------------------- #
def _j1b_tracked(picks) -> list:
    """The species taxids J1b queries: NCBI-REF, non-confounded, with a species taxid."""
    return sorted({p.species_taxid for p in picks
                   if p.source_lane in leads_mod.J1B_LANES
                   and not p.is_confounded and p.species_taxid})


def datasets_taxids(j1b_picks, args):
    """Which taxids each datasets-backed feed needs, and how widely to ask.

    Returns ``(needed, reference_only)``.  J1b's tracked species are asked with
    ``--reference`` -- it only ever fires on an NCBI-designated reference anyway, and
    the unfiltered query for a big species is unusable (Salmonella enterica alone
    returns 135 000 records).  J2's wishlist taxids need the FULL query, because
    ``_best_candidate`` falls back to the newest genome when a species has no
    designated reference at all -- which is the normal case for the rare species on the
    wishlist.  A taxid wanted by both feeds gets the full query: the wide answer
    satisfies J1b too, the narrow one would blind J2.
    """
    reference_only, full = set(), set()
    if args.with_j1b:
        reference_only |= set(_j1b_tracked(j1b_picks))
    if args.with_j2:
        full |= {e.taxid for e in leads_mod.load_todo(args.todo)
                 if e.is_wishlist and e.taxid}
    reference_only -= full
    return sorted(reference_only | full), reference_only


def datasets_cache(j1b_picks, args):
    """The shared ncbi-datasets cache, built ONCE over the UNION of every taxid the
    enabled datasets-backed feeds need (J1b tracked species + J2 wishlist taxids).

    J1b and J2 share one committed cache file; building each feed's slice separately
    would have each ``build()`` wipe the other's species (build replaces, not merges),
    so a subsequent offline run would fail its covers() gate.  Building the union once
    keeps a single coherent ``released_after`` window that covers both feeds.

    ``--released-after`` is a J1b device and reaches only J1b's taxids: the offline gate
    asks for exactly the window this run wants, so a cache rebuilt inside a NARROWER
    window fails loudly instead of answering with an empty candidate list (which reads
    as "NCBI has nothing" and quietly drops every standing lead off the worklist).
    """
    if not (args.with_j1b or args.with_j2 or args.with_lpsn):
        return None
    import datasets_summary as ds  # lazy: only when a datasets-backed feed is requested
    needed, reference_only = datasets_taxids(j1b_picks, args)
    window = args.released_after or ds.WINDOW_ALL
    cache = ds.DatasetsSummaryCache(cache_path=args.datasets_cache,
                                    datasets_bin=args.datasets_bin)
    if args.offline:
        gaps = cache.uncovered(needed, reference_only, window)
        if gaps:
            raise SystemExit(
                "--offline: datasets cache %s does not cover the %d datasets taxids "
                "for window %s: %d gap(s), e.g. %s. Rebuild without --offline (see "
                "--refresh-datasets-cache)."
                % (args.datasets_cache, len(needed), window or "(all)", len(gaps),
                   ["%d (%s)" % g for g in gaps[:5]]))
        log("datasets: offline cache covers %d taxid(s) (%d reference-only) for window %s"
            % (len(needed), len(reference_only), window or "(all)"))
    elif needed:
        cache.build(needed, args.released_after or "",
                    version_epoch=ds.datasets_version(args.datasets_bin), log=log,
                    reference_only=reference_only)
        cache.save()
        log("datasets: built shared cache over %d taxid(s) -- %d reference-only, "
            "%d full -- %d record(s) (%s)"
            % (len(needed), len(reference_only), len(needed) - len(reference_only),
               sum(len(v) for v in cache.by_species.values()), args.datasets_cache))
    return cache


def collect_j1b(picks, args, dcache) -> list:
    """J1b better-rep leads (gated behind --with-j1b; reads the shared datasets cache)."""
    if not args.with_j1b:
        return []
    tracked = _j1b_tracked(picks)
    epoch = dcache.meta.get("datasets_version") or leads_mod.DATASETS_EPOCH
    j1b = leads_mod.j1b_leads(picks, dcache, epoch=epoch)
    log("J1b: %d better-rep lead(s) over %d NCBI-REF tracked species" % (len(j1b), len(tracked)))
    return j1b


def collect_j2(all_picks, args, dcache) -> list:
    """J2 coverage-feeder leads: chromosomes-todo.tsv wishlist (--with-j2) + LPSN gaps
    (--with-lpsn).  Both gated off by default; reads the shared datasets cache + LPSN.

    ``all_picks`` MUST be the whole (unsharded) catalog: an LPSN "we don't carry this
    species" gap is a global-coverage fact, so computing it against a shard slice would
    falsely flag species held in other shards.  J2 leads are pick-slice-independent, so
    every shard emits the same set and collate dedups by lead_key (with the cap
    preserved because the capped sets are identical).
    """
    if not (args.with_j2 or args.with_lpsn):
        return []
    import lpsn as lpsn_mod  # lazy

    todo = leads_mod.load_todo(args.todo) if args.with_j2 else []
    # datasets-derived todo leads carry the datasets tooling epoch (like J1b).
    todo_epoch = (dcache.meta.get("datasets_version") if dcache else None) or leads_mod.DATASETS_EPOCH

    lpsn_records, max_lpsn, lpsn_epoch = [], None, leads_mod.LPSN_EPOCH
    if args.with_lpsn:
        lc = lpsn_mod.LpsnCache(
            cache_path=args.lpsn_cache,
            # The only working refresh path: an operator-downloaded complete-list CSV
            # (or the short-lived per-account link it came from).
            opener=((lambda: lpsn_mod.open_lpsn_source(args.lpsn_csv)) if args.lpsn_csv else None))
        # A LOCAL --lpsn-csv is offline-safe: distilling a file an operator already
        # downloaded touches no network.  That matters, because the alternative is to
        # refresh LPSN without --offline, which drags the 2.7 GB assembly-reports
        # rebuild along for a job that has nothing to do with it.  A download LINK is a
        # real fetch, so that one is still refused under --offline.
        if args.lpsn_csv and args.offline and lpsn_mod.is_remote_source(args.lpsn_csv):
            raise SystemExit("--offline: --lpsn-csv %s is a download link, not a file. "
                             "Save it first, or drop --offline." % args.lpsn_csv)
        if args.lpsn_csv:
            lc.build(args.lpsn_release or "unset", log=log)
            lc.save()
            log("LPSN: wrote cache %s" % args.lpsn_cache)
        elif not lc.records:
            # LPSN is the one feed with no unattended download, so a missing cache
            # cannot be fixed by "rebuild without --offline" the way the NCBI caches
            # can -- the message has to hand the operator the whole manual recipe.
            raise SystemExit("--with-lpsn: LPSN cache %s is missing or empty.\n%s"
                             % (args.lpsn_cache, lpsn_mod.LPSN_REFRESH_HELP))
        lpsn_records, max_lpsn = lc.records, args.max_lpsn
        lpsn_epoch = leads_mod.lpsn_epoch(lc.meta.get("lpsn_release") or lc.release)

    tracked_genera = {leads_mod.species_genus(p.species_name) for p in all_picks if p.species_name}
    have_species = {leads_mod.normalize_species(p.species_name) for p in all_picks if p.species_name}
    j2 = leads_mod.j2_leads(todo, dcache, lpsn_records, tracked_genera, have_species,
                            todo_epoch=todo_epoch, lpsn_epoch=lpsn_epoch,
                            max_lpsn=max_lpsn, since=args.lpsn_since, log=log)
    log("J2: %d coverage-feeder lead(s) (todo=%s, lpsn=%s)"
        % (len(j2), args.with_j2, args.with_lpsn))
    return j2


def collect_leads(rows, args, all_picks) -> list:
    """All Tier 1 leads for this run, across the enabled signals.

    ``rows`` are the coverage rows for this (possibly sharded) run; ``all_picks`` is the
    whole unsharded catalog, needed only by J2's global-coverage gap computation.
    """
    picks = [r.pick for r in rows]
    all_leads = list(leads_mod.j1a_leads(rows, epoch=leads_mod.ANI_EPOCH))
    log("J1a: %d picks actionable of %d probed" % (len(all_leads), len(rows)))
    dcache = datasets_cache(picks, args)             # shared, built ONCE over the union
    all_leads.extend(collect_j1b(picks, args, dcache))     # per-pick: sharded slice is correct
    all_leads.extend(collect_j2(all_picks, args, dcache))  # global: whole catalog
    return all_leads


def shard_picks(picks, shard: int, shards_total: int):
    """Deterministic 1-of-N slice of the picks (by sorted unit_id) for a CI shard.

    Sharding splits only the WORK; the committed cache is global, so an offline shard
    still resolves every pick in its slice.  Collate re-unions the shard leads before
    the delta, so standing/resolved are computed over the whole pick set, never
    per-shard (which would misfire on picks a shard never saw).
    """
    if shards_total <= 1:
        return picks
    ordered = sorted(picks, key=lambda p: p.unit_id)
    return [p for i, p in enumerate(ordered) if i % shards_total == shard]


def emit_leads_fragment(leads, path: str) -> None:
    """Write this shard's raw leads as an NDJSON fragment for collate.py to merge.

    Written atomically (tmp + os.replace) so a shard killed mid-write never leaves a
    truncated fragment that collate would have to treat as a failed shard.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = "%s.tmp.%d" % (path, os.getpid())
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            for lead in leads:
                fh.write(json.dumps(lead.to_dict(), sort_keys=True, ensure_ascii=False) + "\n")
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def cache_meta(path: str) -> dict:
    """The ``meta`` block of a committed cache file (empty when absent/unreadable).

    Deliberately a raw JSON read rather than constructing the cache object: this is
    called from ``collate.py`` too, which has no business importing the fetch layer of
    a tool it never runs.
    """
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh).get("meta") or {}
    except (OSError, ValueError):
        return {}


def tool_versions(datasets_cache_path: str = DEFAULT_DATASETS_CACHE,
                  lpsn_cache_path: str = DEFAULT_LPSN_CACHE) -> dict:
    """Pinned tool/data epochs stamped into every ledger record this run.

    Read from the COMMITTED caches' meta blocks, not from the module constants: the
    constants say ``@unset``, and a ledger record claiming an unset tool version is
    worse than no record at all -- it is the field an SME uses to decide whether a
    verdict is stale.  The constants remain the fallback for a cache that is absent or
    predates the meta block.

    ``collate.py`` recomputes these from the same committed files rather than reading
    them out of a shard fragment.  Recomputing was chosen because the caches are
    committed artifacts of the same tree the shards ran from, so both routes give the
    same answer -- and this one cannot be poisoned by a stale or hand-edited fragment,
    and needs no schema change to the fragment format.
    """
    ds_meta = cache_meta(datasets_cache_path)
    lpsn_meta = cache_meta(lpsn_cache_path)
    release = lpsn_meta.get("lpsn_release")
    return {
        "assembly_reports_epoch": leads_mod.ANI_EPOCH,
        "ncbi_datasets_cli": ds_meta.get("datasets_version") or leads_mod.DATASETS_EPOCH,
        "lpsn": leads_mod.lpsn_epoch(release) if release else leads_mod.LPSN_EPOCH,
    }


# --------------------------------------------------------------------------- #
#  Human summary                                                               #
# --------------------------------------------------------------------------- #
def print_summary(part) -> None:
    s = ledger_mod.summarize_partition(part)
    print()
    print("=" * 78)
    print("Tier 1 curation-triage worklist (run %s)" % part.run_id)
    print("=" * 78)
    print("NEW actionable leads : %d   (the SME worklist this run)" % s["new"])
    print("standing             : %d   (already surfaced, not yet actioned)" % s["standing"])
    print("accepted -> edits    : %d" % s["accepted"])
    print("dismissed            : %d   (reject/defer in reviews.tsv)" % s["dismissed"])
    print("resolved this run    : %d" % s["resolved"])
    print("ledger records added : %d" % s["records_appended"])
    if part.new:
        print()
        print("NEW leads:")
        for lead in part.new:
            print("  - [%s] %s" % (lead.signal, lead.summary))
    print("=" * 78)


def write_markdown(part, args) -> None:
    md = ledger_mod.render_markdown(part)
    if args.summary_md:
        os.makedirs(os.path.dirname(os.path.abspath(args.summary_md)), exist_ok=True)
        with open(args.summary_md, "w", encoding="utf-8") as fh:
            fh.write(md + "\n")
        log("Wrote worklist markdown: %s" % args.summary_md)
    # In CI, append straight to the job summary (the SME-facing surface).
    gh = os.environ.get("GITHUB_STEP_SUMMARY")
    if gh:
        with open(gh, "a", encoding="utf-8") as fh:
            fh.write(md + "\n")


# --------------------------------------------------------------------------- #
#  Driver                                                                      #
# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Tier 1 curation-triage net: the monthly, delta-aware run.")
    p.add_argument("--policy", default=DEFAULT_POLICY)
    p.add_argument("--chromosomes", default=DEFAULT_CHROMOSOMES)
    p.add_argument("--cache", default=DEFAULT_CACHE)
    p.add_argument("--ledger", default=DEFAULT_LEDGER)
    p.add_argument("--reviews", default=DEFAULT_REVIEWS)
    p.add_argument("--out", default=DEFAULT_WORKLIST, help="Per-run worklist TSV artifact.")
    p.add_argument("--summary-md", default=None,
                   help="Write the $GITHUB_STEP_SUMMARY markdown block to this file.")
    p.add_argument("--run-id", default=None,
                   help="Run tag (default: current YYYY-MM). Provenance only, not the "
                        "suppression epoch.")
    p.add_argument("--offline", action="store_true",
                   help="Never hit the network; fail if the committed cache is incomplete.")
    p.add_argument("--no-genbank", action="store_true",
                   help="Skip the large genbank summary when (re)building the cache.")
    p.add_argument("--no-write", action="store_true",
                   help="Compute + summarize but do not append the ledger or write files.")
    # J1b (source-gated better-rep feed) -- off by default; needs ncbi-datasets-cli
    # (or a committed datasets cache for --offline).
    p.add_argument("--with-j1b", action="store_true",
                   help="Enable the J1b better-rep feed (NCBI-REF lanes only).")
    p.add_argument("--datasets-cache", default=DEFAULT_DATASETS_CACHE)
    p.add_argument("--datasets-bin", default="datasets",
                   help="ncbi-datasets-cli executable (e.g. a pinned conda env's bin).")
    p.add_argument("--refresh-datasets-cache", action="store_true",
                   help="Rebuild ONLY the ncbi-datasets cache for the enabled feeds and "
                        "stop. Never touches the 2.7 GB assembly-reports cache, the "
                        "ledger or the worklist.")
    p.add_argument("--released-after", default=None,
                   help="J1b only: ask NCBI for genomes released after this date "
                        "(MM/DD/YYYY); default = all. J2's wishlist is always fetched "
                        "unwindowed -- 'has this species a genome at all?' is a standing "
                        "question. NOTE: this narrows the whole rebuilt cache, it does "
                        "NOT top it up, so a later run wanting a wider window must "
                        "rebuild.")
    # J2 (coverage feeders) -- also off by default.
    p.add_argument("--with-j2", action="store_true",
                   help="Enable the J2 chromosomes-todo.tsv wishlist feeder.")
    p.add_argument("--with-lpsn", action="store_true",
                   help="Also enable the J2 LPSN newly-valid-species gap feeder "
                        "(newest first, capped by --max-lpsn). Runs offline from the "
                        "committed src/curation-triage/lpsn_cache.json; refreshing that "
                        "cache is a manual step -- see --lpsn-csv.")
    p.add_argument("--todo", default=DEFAULT_TODO)
    p.add_argument("--lpsn-cache", default=DEFAULT_LPSN_CACHE)
    p.add_argument("--lpsn-csv", default=None,
                   help="Rebuild the LPSN cache from an operator-downloaded complete-list "
                        "CSV -- a local path or a direct http(s) link. This is the only "
                        "refresh path today: the public URL is a landing page and the "
                        "export needs a registered account. LPSN's REST API would make "
                        "it unattended (future work).")
    p.add_argument("--lpsn-release", default=None, help="LPSN release tag (epoch).")
    p.add_argument("--lpsn-since", default=None,
                   help="Only LPSN names validly published on/after this date "
                        "(YYYY or YYYY-MM-DD). Undated names are dropped.")
    p.add_argument("--max-lpsn", type=int, default=25,
                   help="Cap on LPSN gap leads per run (prevents a first-run flood).")
    # Sharded-CI seam: process a 1-of-N slice and emit its raw leads for collate.py
    # to merge (so the delta is computed over the whole pick set, not per shard).
    p.add_argument("--shard", type=int, default=0, help="This shard's index (0-based).")
    p.add_argument("--shards-total", type=int, default=1, help="Total number of shards.")
    p.add_argument("--emit-leads", default=None,
                   help="Write this run's raw leads as an NDJSON fragment and stop "
                        "(no delta/ledger/worklist). Used by CI shards; merged by collate.py.")
    # Local pre-downloaded copies (dev / CI staging), bypassing the URL stream.
    p.add_argument("--ani-file", default=None)
    p.add_argument("--refseq-file", default=None)
    p.add_argument("--genbank-file", default=None)
    args = p.parse_args(argv)

    run_id = args.run_id or datetime.date.today().strftime("%Y-%m")

    source_lanes = coverage.load_source_lanes(args.chromosomes)
    all_picks = coverage.load_picks(args.policy, source_lanes)   # whole catalog (J2 needs it)
    picks = shard_picks(all_picks, args.shard, args.shards_total)
    log("Loaded %d eligible picks from %s%s" % (
        len(picks), args.policy,
        "" if args.shards_total <= 1 else " (shard %d/%d)" % (args.shard, args.shards_total)))

    # Datasets-only refresh: the J1a assembly-reports cache is 2.7 GB of NCBI flat files
    # and has nothing to do with this feed, so rebuilding it as a side effect of "please
    # re-ask NCBI about our 170 taxids" would be absurd -- and the same reasoning says
    # this run must not append the ledger or write a worklist either.  Runs UNSHARDED
    # over the whole catalog on purpose: build() replaces the cache, it never merges.
    if args.refresh_datasets_cache:
        if args.offline:
            raise SystemExit("--refresh-datasets-cache needs the network; drop --offline.")
        if not (args.with_j1b or args.with_j2):
            raise SystemExit("--refresh-datasets-cache: nothing to fetch -- pass "
                             "--with-j1b and/or --with-j2.")
        datasets_cache(all_picks, args)
        return 0

    # Reuse the coverage probe's cached + mockable network layer verbatim.
    cache = coverage_probe.build_or_load_cache(args, picks)
    rows = coverage.build_coverage(picks, cache)

    leads = collect_leads(rows, args, all_picks)

    # Shard mode: emit raw leads for collate.py and stop (the delta/ledger/worklist
    # are done ONCE, globally, by collate over the merged fragments).
    if args.emit_leads:
        emit_leads_fragment(leads, args.emit_leads)
        log("Emitted %d leads: %s" % (len(leads), args.emit_leads))
        return 0

    prior = ledger_mod.load_ledger(args.ledger)
    reviews = ledger_mod.load_reviews(args.reviews)
    log("Prior ledger: %d records; reviews: %d rows" % (len(prior), len(reviews)))

    part = ledger_mod.partition(leads, prior, reviews, run_id=run_id,
                                params=leads_mod.DEFAULT_PARAMS,
                                tool_versions=tool_versions(args.datasets_cache,
                                                            args.lpsn_cache))

    if not args.no_write:
        ledger_mod.append_ledger(args.ledger, part.records)
        log("Appended %d ledger records: %s" % (len(part.records), args.ledger))
        ledger_mod.write_worklist(part, args.out)
        log("Wrote worklist: %s (%d rows)" % (args.out, len(ledger_mod.worklist_rows(part))))
        write_markdown(part, args)

    print_summary(part)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
