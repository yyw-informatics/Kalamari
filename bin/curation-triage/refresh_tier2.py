#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
refresh_tier2.py -- run a Tier-2 cache refresh in shards and merge the pieces.

``tier2_confirm.py`` can refresh its own caches, but only as one process that holds
every pinned tool at once.  That does not fit a CI runner: the eleven pinned
environments total about 13 GB.  So a refresh is **split by tool**, each shard writes
a cache *fragment*, and a merge step folds the fragments back into the committed
caches.

Why by tool and not by unit
---------------------------
The sub-species panels are cliques -- eleven Salmonella units compare against each
other, four Listeria lineages against each other -- so a per-unit split hands the
same clique to several shards and computes it several times over.  The tool axis has
no such overlap: every (tool, accession) job belongs to exactly one tool, and a
runner then installs exactly one environment.

Three modes, one plan
---------------------
All three build the **same plan** as a real Tier-2 run, by calling
``tier2_confirm.build_plan``.  That is the point: a CI cache key, a shard's job list
and the merge's coverage check must never disagree with what the monthly run needs.

    # 1. the accession list -- this is the CI genome-cache key
    python3 bin/curation-triage/refresh_tier2.py --print-accessions

    # 2. one tool shard (or the skani shard) -> a cache FRAGMENT
    python3 bin/curation-triage/refresh_tier2.py --tools sistr,seqsero2 \
        --typer-env-dir <env-dir> --out frag-sistr.json
    python3 bin/curation-triage/refresh_tier2.py --skani \
        --skani-bin <env-dir>/skani/bin/skani --out frag-skani.json

    # 3. fold the fragments into a full cache, over the cache they replace
    python3 bin/curation-triage/refresh_tier2.py --merge frag-*.json \
        --base src/curation-triage/tier2_typer_cache.json \
        --out src/curation-triage/tier2_typer_cache.json

The merge is a pure function (``merge_caches``) and it is **strict on purpose**: two
fragments that disagree about the same key are an error, not a last-writer-wins
overwrite.  A silent overwrite there would let a shard that ran a stale environment
quietly replace a good number, and the resulting cache is what every Tier-2 verdict
is then computed from.  For the same reason a merge that does not cover the whole
plan is refused unless ``--allow-partial`` says the gap is deliberate.

Why a merge also needs a BASE and a shrink guard
------------------------------------------------
Covering the plan is not the same thing as keeping the cache whole.  The plan routes
only the tools that decide a verdict, so the shard matrix has no ECTyper shard --
ECTyper is context-only and its 943 MiB download is opt-in.  The committed typer
cache nevertheless holds one ``ectyper|...`` result and its tool version, because a
past run with ``--with-ectyper`` put them there.  Merging the shard fragments alone
therefore produced a cache that covered every planned job and was still *smaller*
than the one it replaced, and nothing downstream could see it: context keys sit
outside the Tier-2 fingerprint, so the delta stayed ``recurring`` and both the
coverage replay and the four-point update gate passed on the shrunken cache.

Two things close that hole:

``--base FILE``
    A floor under the fragments.  It fills the keys no fragment refreshed and loses
    every key a fragment does answer, so a full refresh can still MOVE a call (the
    fragment wins) while carrying the results no shard produces (the base keeps
    them).  Bases are not peers of the fragments, so this does not soften the
    conflict-strictness above: two *fragments* that disagree are still an error.

``--out`` shrink guard
    Before writing, the merge compares itself against the cache it is about to
    overwrite and refuses to drop a result key, a ``meta.tool_versions`` entry or any
    other meta value it already holds.  ``--allow-drop`` is the deliberate escape and
    it names everything it drops.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import coverage  # noqa: E402
import genomes as genomes_mod  # noqa: E402
import ledger as ledger_mod  # noqa: E402
import markers as markers_mod  # noqa: E402
import policy as policy_mod  # noqa: E402
import skani_runner  # noqa: E402
import tier2_confirm  # noqa: E402
import typers as typers_mod  # noqa: E402

KIND_ANI = "ani"
KIND_TYPER = "typer"

# Recomputed from the merged results, never merged: a fragment's counts describe the
# fragment, and carrying them forward would describe the merged cache wrongly.
RECOMPUTED_META = ("n_results", "n_screened_out")


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


class MergeConflict(ValueError):
    """Two fragments disagree about the same key -- refuse rather than pick one."""


# --------------------------------------------------------------------------- #
#  Pure: cache identity, merging, coverage                                     #
# --------------------------------------------------------------------------- #
def cache_kind(blob: dict) -> str:
    """``ani`` or ``typer`` for a cache payload; '' when it cannot be told apart.

    Inferred rather than stamped into the file, so a fragment is byte-comparable with
    the committed caches and no new meta key appears in a refreshed commit.
    """
    meta = blob.get("meta") or {}
    if "skani_version" in meta or "screen" in meta:
        return KIND_ANI
    if "tool_versions" in meta or "amrfinder_database_version" in meta:
        return KIND_TYPER
    results = blob.get("results") or {}
    values = [v for v in results.values() if isinstance(v, dict)]
    if values and all(("ani" in v or "screened_out" in v) for v in values):
        return KIND_ANI
    if values:
        return KIND_TYPER
    return ""


def _fold_peers(blobs) -> dict:
    """Fold ``(label, blob)`` **peers** into raw merge parts; disagreement is an error.

    Peers, not layers: nothing here outranks anything else, so two blobs that answer
    the same key differently cannot be reconciled and must be refused.  Returns the
    parts (schema, results, meta, tool versions) rather than a finished cache, so the
    same strict fold can be run over the fragments and over the base files.
    """
    results, results_from = {}, {}
    meta, meta_from = {}, {}
    versions, versions_from = {}, {}
    schema, schema_from = None, ""
    want_screened = False

    for label, blob in blobs:
        blob_meta = dict(blob.get("meta") or {})
        blob_schema = blob.get("schema")
        if schema is None:
            schema, schema_from = blob_schema, label
        elif blob_schema != schema:
            raise MergeConflict("schema: %s says %r, %s says %r"
                                % (label, blob_schema, schema_from, schema))

        for key, val in (blob.get("results") or {}).items():
            if key in results and results[key] != val:
                raise MergeConflict(
                    "result %s: %s says %r, %s already says %r"
                    % (key, label, val, results_from[key], results[key]))
            results[key] = val
            results_from[key] = label

        for tool, ver in (blob_meta.pop("tool_versions", None) or {}).items():
            if tool in versions and versions[tool] != ver:
                raise MergeConflict(
                    "tool version %s: %s says %r, %s already says %r"
                    % (tool, label, ver, versions_from[tool], versions[tool]))
            versions[tool] = ver
            versions_from[tool] = label

        for key, val in blob_meta.items():
            if key in RECOMPUTED_META:
                want_screened = want_screened or key == "n_screened_out"
                continue
            if key in meta and meta[key] != val:
                raise MergeConflict("meta %s: %s says %r, %s already says %r"
                                    % (key, label, val, meta_from[key], meta[key]))
            meta[key] = val
            meta_from[key] = label

    return {"schema": schema, "schema_from": schema_from, "results": results,
            "meta": meta, "versions": versions, "want_screened": want_screened}


def merge_caches(fragments, base=()) -> dict:
    """Pure: fold ``(label, blob)`` fragments into one ``{schema, meta, results}``.

    Union the results and ``meta.tool_versions``; recompute the counts; carry every
    other meta value (the skani epoch, the screen thresholds, the single
    ``amrfinder_database_version``) forward.  Any disagreement between two fragments
    -- a differing ANI for the same pair, two skani versions, two AMRFinder database
    dates -- raises ``MergeConflict``.  A cache that silently kept one of two
    conflicting values would make every verdict derived from it unattributable.

    ``base`` is the cache the merge replaces, and it is a **floor, not a peer**: it
    supplies only the keys no fragment answers.  A fragment that answers a base key
    wins without a conflict, because that is exactly what a refresh is -- the shard
    just re-ran the tool.  What the base is for is the other direction: the keys no
    shard produces at all (the context-only ECTyper result and its tool version) would
    otherwise vanish from a cache that still covers the whole plan.
    """
    top = _fold_peers(fragments)
    under = _fold_peers(base)

    if top["schema"] is None:
        top["schema"], top["schema_from"] = under["schema"], under["schema_from"]
    elif under["schema"] is not None and under["schema"] != top["schema"]:
        # A base written under a different schema is not a stale number the refresh
        # supersedes; it is a file that means something else.
        raise MergeConflict("schema: base %s says %r, %s says %r"
                            % (under["schema_from"], under["schema"],
                               top["schema_from"], top["schema"]))
    for key, val in under["results"].items():
        top["results"].setdefault(key, val)
    for tool, ver in under["versions"].items():
        top["versions"].setdefault(tool, ver)
    for key, val in under["meta"].items():
        top["meta"].setdefault(key, val)

    meta, results = top["meta"], top["results"]
    if top["versions"]:
        meta["tool_versions"] = dict(sorted(top["versions"].items()))
    meta["n_results"] = len(results)
    screened = [v for v in results.values() if isinstance(v, dict) and "screened_out" in v]
    if screened or top["want_screened"] or under["want_screened"]:
        meta["n_screened_out"] = sum(1 for v in screened if v.get("screened_out"))
    return {"schema": top["schema"] if top["schema"] is not None else 1, "meta": meta,
            "results": dict(sorted(results.items()))}


def uncovered_pairs(blob: dict, pairs) -> list:
    """Planned ANI pairs a merged cache cannot answer (either direction counts)."""
    results = blob.get("results") or {}
    return [(q, r) for q, r in pairs
            if "%s|%s" % (q, r) not in results and "%s|%s" % (r, q) not in results]


def uncovered_jobs(blob: dict, jobs) -> list:
    """Planned ``(tool, accession)`` jobs a merged cache cannot answer."""
    results = blob.get("results") or {}
    return [(t, a) for t, a in jobs if "%s|%s" % (t, a) not in results]


def dropped_keys(previous: dict, merged: dict) -> list:
    """Everything ``merged`` would delete from ``previous``, named for an error message.

    Covering the plan is a weaker promise than keeping the cache whole: the plan lists
    only the tools that decide a verdict, so a merged cache can answer every planned
    job and still be missing entries the committed cache holds -- the context-only
    ECTyper result is exactly that.  Nothing re-derives a dropped key, and every check
    downstream reads only what IS in the cache, so the loss leaves no trace anywhere.
    A dropped ``meta`` value is the same kind of loss: a tool version and the AMRFinder
    database date are epochs, and losing an epoch is what stops a dismissed lead from
    re-opening when the tool moves.
    """
    prev_meta = (previous.get("meta") or {})
    new_meta = (merged.get("meta") or {})
    lost = sorted(set(previous.get("results") or {}) - set(merged.get("results") or {}))
    lost += ["meta.tool_versions.%s" % t
             for t in sorted(set(prev_meta.get("tool_versions") or {})
                             - set(new_meta.get("tool_versions") or {}))]
    # RECOMPUTED_META describes the file, not the evidence, so it may legitimately move.
    lost += ["meta.%s" % k for k in sorted(set(prev_meta) - set(new_meta))
             if k != "tool_versions" and k not in RECOMPUTED_META]
    return lost


def unversioned_tools(blob: dict, jobs) -> list:
    """Planned tools with no recorded version in the merged meta.

    A cache whose results are complete but whose meta is missing a tool has lost that
    tool's epoch, and the epoch is what re-opens a dismissed lead when the tool moves.
    """
    versions = (blob.get("meta") or {}).get("tool_versions") or {}
    return sorted({t for t, _ in jobs} - set(versions))


# --------------------------------------------------------------------------- #
#  Plan (identical to what a Tier-2 run computes)                              #
# --------------------------------------------------------------------------- #
def build_plan(args):
    """The same plan ``tier2_confirm.py`` builds, so the shards match the real run."""
    source_lanes = coverage.load_source_lanes(args.chromosomes)
    picks = coverage.load_picks(args.policy, source_lanes)
    unit_meta = tier2_confirm.load_unit_meta(args.policy)
    confounded_entries = policy_mod.load_confounded_manifest(args.confounded_manifest)
    manifest = markers_mod.load_manifest(args.marker_manifest)
    survivors = tier2_confirm.survivor_unit_ids(ledger_mod.load_ledger(args.ledger))
    triggers = tier2_confirm.trigger_units(picks, survivors, unit_meta)
    ts_cache = genomes_mod.TypeStrainCache(
        cache_path=os.path.join(args.cache_dir, tier2_confirm.TYPE_STRAIN_CACHE_NAME))
    plan = tier2_confirm.build_plan(triggers, manifest, confounded_entries, ts_cache,
                                    args.panel_dir)
    if args.with_ectyper:
        plan.typer_jobs = sorted(set(plan.typer_jobs) | set(context_jobs(plan)))
    return plan


def context_jobs(plan) -> list:
    """``(tool, accession)`` jobs for the context-only tools of the planned tasks.

    Off by default because ECTyper downloads a 943 MiB Zenodo mash sketch the first
    time it runs.  Nothing in the pinned pipeline needs it -- ECTyper decides no
    verdict -- so paying for that download has to be an explicit choice.
    """
    jobs = set()
    for task_id, context_tools in typers_mod.TASK_CONTEXT_TOOLS.items():
        # The accessions a task's decisive tools already run on are exactly the
        # accessions its context tool would add value for.
        decisive = set(typers_mod.TASK_TOOLS.get(task_id, []))
        accessions = {a for tool, a in plan.typer_jobs if tool in decisive}
        jobs |= {(tool, a) for tool in context_tools for a in accessions}
    return sorted(jobs)


# --------------------------------------------------------------------------- #
#  Fragment writing                                                            #
# --------------------------------------------------------------------------- #
def write_cache(path: str, payload: dict) -> None:
    """Write a cache/fragment in exactly the committed caches' format."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = "%s.tmp.%d" % (path, os.getpid())
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=1, sort_keys=True)
            fh.write("\n")
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def read_cache(path: str) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def genome_store(args, accessions):
    store = genomes_mod.GenomeStore(args.genome_cache, offline=args.offline,
                                    datasets_bin=args.datasets_bin)
    if not args.offline:
        store.prefetch(sorted(accessions), log=log)
    return store


def run_typer_shard(args, plan, tools) -> dict:
    """Run one tool shard and return its fragment payload."""
    jobs = [(t, a) for t, a in plan.typer_jobs if t in tools]
    if not jobs:
        raise SystemExit("no planned job for tool(s) %s; the plan covers %s"
                         % (",".join(sorted(tools)),
                            ",".join(sorted({t for t, _ in plan.typer_jobs}))))
    problems = typers_mod.preflight(tools, args.typer_env_dir)
    if problems:
        if not args.allow_unpinned:
            raise SystemExit("preflight failed -- refusing to refresh with unpinned "
                             "tools:\n  %s\nInstall the pinned environments (see "
                             "setup_envs.sh) or pass --allow-unpinned."
                             % "\n  ".join(problems))
        log("WARNING: preflight problems, continuing under --allow-unpinned:\n  %s"
            % "\n  ".join(problems))

    store = genome_store(args, [a for _, a in jobs])
    # A fresh, empty cache: a shard must RE-RUN its jobs, so it must not be seeded
    # from the committed cache (which would make every job a cache hit).
    cache = typers_mod.TyperCache(cache_path=None, env_dir=args.typer_env_dir)
    cache.build(jobs, store, args.workdir or os.path.dirname(os.path.abspath(args.out)),
                log=log, only_tools=tools)
    if cache.misses and not args.allow_misses:
        raise SystemExit("%d of %d job(s) produced no result: %s. A cache committed "
                         "with silently absent calls is worse than no refresh; fix the "
                         "environment or pass --allow-misses."
                         % (len(cache.misses), len(jobs), ", ".join(sorted(cache.misses))))
    log("typer shard %s: %d result(s), %d miss(es)"
        % (",".join(sorted(tools)), len(cache.results), len(cache.misses)))
    return {"schema": 1, "meta": cache.meta, "results": dict(sorted(cache.results.items()))}


def run_skani_shard(args, plan) -> dict:
    """Run the whole skani slice and return its fragment payload."""
    skani_bin = args.skani_bin
    resolved = skani_bin if os.sep in skani_bin else shutil.which(skani_bin)
    if not (resolved and os.path.exists(resolved) and os.access(resolved, os.X_OK)):
        if not args.allow_unpinned:
            raise SystemExit("preflight failed -- skani binary %r does not resolve to an "
                             "executable. Install the pinned environment (see "
                             "setup_envs.sh) or pass --allow-unpinned." % skani_bin)
        log("WARNING: skani binary %r does not resolve; continuing under --allow-unpinned"
            % skani_bin)

    store = genome_store(args, plan.accessions)
    cache = skani_runner.AniCache(cache_path=None, skani_bin=skani_bin)
    cache.build(plan.ani_pairs, store, labels=plan.labels, log=log)
    if cache.misses and not args.allow_misses:
        raise SystemExit("%d genome(s) could not be resolved for skani: %s. Fix the "
                         "genome cache or pass --allow-misses."
                         % (len(set(cache.misses)), ", ".join(sorted(set(cache.misses)))))
    log("skani shard: %d result(s) over %d planned pair(s)"
        % (len(cache.results), len(plan.ani_pairs)))
    return {"schema": 1, "meta": cache.meta, "results": dict(sorted(cache.results.items()))}


def read_previous(path: str):
    """The cache a merge is about to overwrite, or ``None`` when there is none.

    An unreadable file at ``--out`` is refused rather than ignored: "I could not tell
    what was there, so I overwrote it" is the same silent loss this guard exists for.
    """
    if not os.path.exists(path):
        return None
    try:
        return read_cache(path)
    except (OSError, ValueError) as exc:
        raise SystemExit("merge refused -- cannot read the cache it would overwrite "
                         "(%s): %s. Move that file aside if it is not a cache." % (path, exc))


def merge_mode(args, plan) -> dict:
    """Merge fragments, check the merged cache covers the plan, and write it."""
    fragments = [(os.path.basename(p), read_cache(p)) for p in args.merge]
    base = [(os.path.basename(p), read_cache(p)) for p in (args.base or [])]
    try:
        merged = merge_caches(fragments, base=base)
    except MergeConflict as exc:
        raise SystemExit("merge refused -- %s. Re-run the shard that produced the stale "
                         "value; do not pick one of the two by hand." % exc)
    kind = args.kind or cache_kind(merged)
    if not kind:
        raise SystemExit("cannot tell an ANI cache from a typer cache in %s; pass "
                         "--kind ani|typer." % ", ".join(args.merge))

    if kind == KIND_ANI:
        gaps = ["%s|%s" % p for p in uncovered_pairs(merged, plan.ani_pairs)]
        total, what = len(plan.ani_pairs), "ANI comparison(s)"
        extra = []
    else:
        gaps = ["%s|%s" % j for j in uncovered_jobs(merged, plan.typer_jobs)]
        total, what = len(plan.typer_jobs), "typer job(s)"
        # A tool whose results are all present but whose version is not recorded has
        # lost its epoch, which is what re-opens a dismissed lead when the tool moves.
        extra = ["%s (no recorded tool version)" % t
                 for t in unversioned_tools(merged, plan.typer_jobs)]
    if (gaps or extra) and not args.allow_partial:
        raise SystemExit("merge of %d fragment(s) covers %d of %d %s. Missing: %s. Run "
                         "the missing shard(s), or pass --allow-partial for a deliberate "
                         "partial cache."
                         % (len(fragments), total - len(gaps), total, what,
                            (gaps + extra)[:6]))
    if gaps or extra:
        log("WARNING: --allow-partial: %d of %d %s not covered%s"
            % (len(gaps), total, what,
               "; %d tool(s) without a recorded version" % len(extra) if extra else ""))

    if base:
        report_base_use(fragments, base)
    check_no_shrink(args, merged, kind, bool(base))

    log("merged %d fragment(s)%s into %d result(s) [%s]"
        % (len(fragments), " over %d base file(s)" % len(base) if base else "",
           len(merged["results"]), kind))
    return merged


def report_base_use(fragments, base) -> None:
    """Say out loud which keys the base carried and which ones the fragments moved.

    A carried key is evidence nobody re-ran this round, and a moved key is the event
    the whole refresh exists to produce.  Both are silent otherwise, because the base
    resolves without a conflict by design.
    """
    from_fragments = {}
    for _, blob in fragments:
        from_fragments.update(blob.get("results") or {})
    from_base = {}
    for _, blob in base:
        from_base.update(blob.get("results") or {})
    carried = sorted(set(from_base) - set(from_fragments))
    moved = sorted(k for k in set(from_base) & set(from_fragments)
                   if from_base[k] != from_fragments[k])
    if carried:
        log("base: carried %d result(s) forward that no fragment refreshed: %s"
            % (len(carried), name_some(carried)))
    if moved:
        log("base: %d result(s) MOVED -- the fragment supersedes the base: %s"
            % (len(moved), name_some(moved)))


def name_some(keys, limit: int = 10) -> str:
    """Name the keys in an operator-facing message, capped so it stays readable."""
    shown = ", ".join(keys[:limit])
    return shown if len(keys) <= limit else "%s (+%d more)" % (shown, len(keys) - limit)


def committed_cache_path(args, kind: str):
    """Where a merged cache of this kind is destined to land, or ``None`` if unknown."""
    name = {KIND_ANI: tier2_confirm.ANI_CACHE_NAME,
            KIND_TYPER: tier2_confirm.TYPER_CACHE_NAME}.get(kind)
    return os.path.join(args.cache_dir, name) if name else None


def check_no_shrink(args, merged: dict, kind: str, has_base: bool) -> None:
    """Refuse a merge that would delete something the cache at ``--out`` already holds."""
    previous = read_previous(args.out)
    if previous is None and not has_base:
        # A merge STAGED to a scratch path and moved into place afterwards (which is how
        # CI does it) has no file to compare against here, so this guard cannot see the
        # drop.  Say so rather than pass quietly: silence is what made the lost ECTyper
        # result invisible in the first place.
        destination = committed_cache_path(args, kind)
        if destination and os.path.exists(destination):
            log("WARNING: %s does not exist, so this merge cannot be checked against the "
                "cache it replaces. If it is staged onto %s, pass --base %s (or write "
                "--out straight there); otherwise a key only that cache holds is dropped "
                "silently." % (args.out, destination, destination))
    lost = dropped_keys(previous, merged) if previous else []
    if not lost:
        return
    if not args.allow_drop:
        raise SystemExit(
            "merge refused -- writing %s would DROP %d entr(y/ies) it already holds: "
            "%s. Nothing re-derives a dropped key and every check downstream reads only "
            "what IS in the cache, so this loss would be invisible. Pass --base %s to "
            "carry them forward (a fragment still wins over the base), or --allow-drop "
            "if the drop is deliberate."
            % (args.out, len(lost), name_some(lost), args.out))
    log("WARNING: --allow-drop: %s loses %d entr(y/ies): %s"
        % (args.out, len(lost), name_some(lost)))


# --------------------------------------------------------------------------- #
#  Driver                                                                      #
# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Sharded Tier-2 cache refresh (fragments + merge).")
    # Modes.
    p.add_argument("--print-accessions", action="store_true",
                   help="Print the plan's accessions, one per line, sorted "
                        "(the CI genome-cache key).")
    p.add_argument("--tools", default="",
                   help="Comma-separated tool shard to run (e.g. sistr,seqsero2).")
    p.add_argument("--skani", action="store_true", help="Run the skani (ANI) shard.")
    p.add_argument("--merge", nargs="+", default=None,
                   help="Cache fragments to merge into --out.")
    p.add_argument("--base", nargs="+", default=None,
                   help="Cache(s) the merge replaces, used as a floor: they supply only "
                        "the keys no fragment answers, so a shard still moves a call but "
                        "the results no shard produces (context-only tools) survive.")
    p.add_argument("--out", default=None, help="Fragment / merged cache to write.")
    p.add_argument("--kind", choices=(KIND_ANI, KIND_TYPER), default=None,
                   help="Force the cache kind when merging (normally inferred).")
    # Plan inputs (defaults match tier2_confirm.py).
    p.add_argument("--policy", default=tier2_confirm.DEFAULT_POLICY)
    p.add_argument("--chromosomes", default=tier2_confirm.DEFAULT_CHROMOSOMES)
    p.add_argument("--confounded-manifest", default=tier2_confirm.DEFAULT_CONFOUNDED)
    p.add_argument("--marker-manifest", default=tier2_confirm.DEFAULT_MANIFEST)
    p.add_argument("--panel-dir", default=tier2_confirm.DEFAULT_PANEL_DIR)
    p.add_argument("--ledger", default=tier2_confirm.DEFAULT_LEDGER)
    p.add_argument("--cache-dir", default=tier2_confirm.DEFAULT_CACHE_DIR)
    p.add_argument("--genome-cache", default=tier2_confirm.DEFAULT_GENOME_CACHE)
    p.add_argument("--workdir", default=None,
                   help="Scratch directory for typer runs (default: beside --out).")
    # Tools.
    p.add_argument("--datasets-bin", default=genomes_mod.DEFAULT_DATASETS_BIN)
    p.add_argument("--skani-bin", default=skani_runner.DEFAULT_SKANI_BIN)
    p.add_argument("--typer-env-dir", default=None,
                   help="Directory of per-tool conda envs (<dir>/<env>/bin/...).")
    p.add_argument("--with-ectyper", action="store_true",
                   help="Plan the context-only ECTyper jobs too. Off by default: "
                        "ECTyper downloads a 943 MiB Zenodo sketch and decides no verdict.")
    # Escapes, each one loud.
    p.add_argument("--offline", action="store_true",
                   help="Never download a genome; run only over the cached FASTAs.")
    p.add_argument("--allow-unpinned", action="store_true",
                   help="Downgrade the preflight failure to a warning.")
    p.add_argument("--allow-misses", action="store_true",
                   help="Write a fragment even when some jobs produced no result.")
    p.add_argument("--allow-partial", action="store_true",
                   help="Accept a merge that does not cover the whole plan.")
    p.add_argument("--allow-drop", action="store_true",
                   help="Write a merged cache even when it DELETES result keys or meta "
                        "values the file at --out already holds. Every dropped entry is "
                        "named in the log.")
    args = p.parse_args(argv)

    tools = sorted({t.strip() for t in args.tools.split(",") if t.strip()})
    modes = [bool(args.print_accessions), bool(tools or args.skani), bool(args.merge)]
    if sum(modes) != 1:
        p.error("pick exactly one mode: --print-accessions, --tools/--skani, or --merge")
    if tools and args.skani:
        p.error("--tools and --skani write different caches; run them separately")
    if (args.base or args.allow_drop) and not args.merge:
        p.error("--base and --allow-drop only mean something for --merge")
    unknown = [t for t in tools if t not in typers_mod.NORMALIZERS]
    if unknown:
        p.error("unknown tool(s): %s" % ", ".join(unknown))
    if not args.print_accessions and not args.out:
        p.error("--out is required for a shard or a merge")

    plan = build_plan(args)
    log("Plan: %d ANI pair(s), %d typer job(s), %d genome(s)"
        % (len(plan.ani_pairs), len(plan.typer_jobs), len(plan.accessions)))

    if args.print_accessions:
        for accession in plan.accessions:
            print(accession)
        return 0

    if args.merge:
        payload = merge_mode(args, plan)
    elif args.skani:
        payload = run_skani_shard(args, plan)
    else:
        payload = run_typer_shard(args, plan, set(tools))

    write_cache(args.out, payload)
    log("Wrote %s (%d result(s))" % (args.out, len(payload["results"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
