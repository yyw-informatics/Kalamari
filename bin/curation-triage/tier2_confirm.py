#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tier2_confirm.py -- Tier 2 targeted confirmation for the curation-triage tool.

Tier 1 reads NCBI's metadata verdict on all 249 eligible picks and hands
the SME a short worklist.  Tier 2 runs only on what that worklist and the policy
table single out, and answers at the sequence level:

    1. **Tier-1 survivors** -- our pick against the type strain of its declared
       species, ANI with the aligned fraction as a QC gate.  This turns "NCBI says
       inconclusive" into a number the SME can act on.
    2. **Confounded units (mandatory)** -- for these the metadata cannot answer the
       question at all, so they are confirmed every run whether or not Tier 1
       flagged them:
         * **regime A** (Shigella/E. coli, anthracis/cereus, pestis/pseudotb) --
           the pinned typer decides and the marker call OVERRIDES ANI; the ANI
           number is attached as context only.
         * **regime B** (Listeria lineages, Salmonella subspecies, botulinum
           groups) -- the sub-species panel checks placement and distinctness, and
           the typer supplies the label.

Tier 2 is a **new signal on the same ledger**: its leads carry ``signal='T2'`` and
flow through the same delta + reviews-suppression + worklist loop Tier 1 uses.  A
confirming Tier-2 verdict annotates the Tier-1 row and never removes it -- this is
an informational tool, so the SME still decides.

Everything external is cached and injectable (``genomes.py``, ``skani_runner.py``,
``typers.py``), so an offline run replays from the committed caches:

    # offline monthly run from the committed caches (no network, no tools):
    python3 bin/curation-triage/tier2_confirm.py --offline --run-id 2026-08

    # refresh: download the genomes, run skani + the typers, then evaluate.
    # --typer-env-dir is where the pinned per-tool conda envs live (setup_envs.sh);
    # a refresh that cannot find them stops before it runs anything:
    python3 bin/curation-triage/tier2_confirm.py --run-id 2026-08 \\
        --typer-env-dir .cache/curation-triage/envs
"""

from __future__ import annotations

import argparse
import csv
import datetime
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import coverage  # noqa: E402
import genomes as genomes_mod  # noqa: E402
import ledger as ledger_mod  # noqa: E402
import markers as markers_mod  # noqa: E402
import policy as policy_mod  # noqa: E402
import skani_runner  # noqa: E402
import tier2  # noqa: E402
import typers as typers_mod  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")
DEFAULT_POLICY = os.path.join(SRC, "taxon_policy.tsv")
DEFAULT_CHROMOSOMES = os.path.join(REPO_ROOT, "src", "chromosomes.tsv")
DEFAULT_CONFOUNDED = os.path.join(SRC, "confounded_manifest.tsv")
DEFAULT_MANIFEST = os.path.join(SRC, "marker_manifest.tsv")
DEFAULT_PANEL_DIR = os.path.join(SRC, "panels")
DEFAULT_LEDGER = os.path.join(SRC, "ledger.ndjson")
DEFAULT_REVIEWS = os.path.join(SRC, "reviews.tsv")
DEFAULT_WORKLIST = os.path.join(SRC, "worklist_tier2.tsv")
DEFAULT_CACHE_DIR = SRC
# Genome FASTA is megabytes per accession and is never committed (see genomes.py).
DEFAULT_GENOME_CACHE = os.path.join(REPO_ROOT, ".cache", "curation-triage", "genomes")

ANI_CACHE_NAME = "tier2_ani_cache.json"
TYPER_CACHE_NAME = "tier2_typer_cache.json"
TYPE_STRAIN_CACHE_NAME = "type_strain_cache.json"


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


# --------------------------------------------------------------------------- #
#  Trigger set                                                                 #
# --------------------------------------------------------------------------- #
def load_unit_meta(policy_path: str) -> dict:
    """unit_id -> the policy columns Tier 2 needs beyond ``coverage.Pick``."""
    out = {}
    with open(policy_path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            note = (r.get("note") or "").strip()
            out[r["unit_id"]] = {
                "kalamari_taxid": (r.get("kalamari_taxid") or "").strip(),
                "confounded_label": (note.split("confounded:", 1)[1].strip()
                                     if "confounded:" in note else ""),
            }
    return out


def survivor_unit_ids(ledger_records) -> set:
    """Units currently on the Tier-1 worklist (the leads Tier 2 confirms).

    Read from the latest record per lead key, so a lead that has since resolved is
    not re-confirmed and a standing lead still is.  Tier-2's own records are skipped
    -- Tier 2 must not feed itself.
    """
    out = set()
    for key, rec in ledger_mod.latest_by_key(ledger_records).items():
        if rec.get("signal") == "T2" or key.startswith("t2:"):
            continue
        if rec.get("delta") == "resolved":
            continue
        if not rec.get("actionable_flags"):
            continue
        if rec.get("unit_id"):
            out.add(rec["unit_id"])
    return out


def trigger_units(picks, survivors, meta, want_survivors=True, want_confounded=True):
    """The units Tier 2 processes, with the reason each one qualified."""
    triggers = []
    for pick in sorted(picks, key=lambda p: p.unit_id):
        reasons = []
        if want_confounded and pick.is_confounded:
            reasons.append("confounded regime %s" % pick.regime)
        if want_survivors and pick.unit_id in survivors:
            reasons.append("Tier-1 survivor")
        if reasons:
            triggers.append((pick, reasons, meta.get(pick.unit_id, {})))
    return triggers


def tasks_for_pick(pick, unit_meta, manifest, confounded_entries):
    """Marker task ids for a unit: its species taxid plus its group's trigger taxid.

    The group trigger matters for *Bacillus*: our units are *B. cereus* (1396) and
    *B. thuringiensis* (1428), but the manifest triggers on the *B. cereus group*
    (86661).  Tier 0 already resolved that lineage containment and recorded the group
    label, so it is looked up here rather than climbing the taxonomy again.
    """
    label = (unit_meta or {}).get("confounded_label", "")
    taxids = {pick.species_taxid} if pick.species_taxid else set()
    taxids |= {e.trigger_taxid for e in confounded_entries if e.label == label}
    return markers_mod.tasks_for_taxids(manifest, sorted(t for t in taxids if t))


# --------------------------------------------------------------------------- #
#  Planning: what has to be computed before anything can be decided             #
# --------------------------------------------------------------------------- #
class Plan:
    """The genomes, ANI pairs and typer jobs one run needs."""

    def __init__(self):
        self.ani_pairs = []        # (query_accession, reference_accession)
        self.labels = {}           # reference accession -> human label
        self.typer_jobs = []       # (tool, accession)
        self.panels = {}           # unit_id -> [PanelMember]
        self.type_strain = {}      # unit_id -> the type-strain cache entry ({} if none)
        self.tasks = {}            # unit_id -> [task_id]

    @property
    def accessions(self):
        out = {q for q, _ in self.ani_pairs} | {r for _, r in self.ani_pairs}
        out |= {acc for _, acc in self.typer_jobs}
        return sorted(a for a in out if a)


def build_plan(triggers, manifest, confounded_entries, ts_cache, panel_dir) -> Plan:
    plan = Plan()
    for pick, _reasons, meta in triggers:
        tasks = tasks_for_pick(pick, meta, manifest, confounded_entries)
        plan.tasks[pick.unit_id] = tasks

        # ANI vs the type strain: for regime A this is context only, but it is still
        # worth computing -- an SME reading a marker override wants to see the number
        # the override replaced.
        # By the pick's OWN assembly first: NCBI names a per-assembly declared type
        # strain, which differs between the subspecies of one species (see
        # genomes.TypeStrainCache).  The species taxid is only the fallback.
        #
        # The whole entry is kept, not just the accession: NCBI also tells us when the
        # pick IS the type strain, and when the species has no type material at all.
        # Both are answers that need no comparison, so neither adds an ANI pair.
        entry = ts_cache.get(pick.base, pick.species_taxid) or {}
        plan.type_strain[pick.unit_id] = entry
        ts_acc = entry.get("accession", "")
        if (ts_acc and not entry.get("is_self_type")
                and not entry.get("has_no_type_material")
                and ts_acc != pick.parent_assembly):
            plan.ani_pairs.append((pick.parent_assembly, ts_acc))
            plan.labels[ts_acc] = "type strain of %s" % (pick.species_name or "the species")

        # Sub-species panel (regime B only; regime A is resolved by markers).
        if pick.regime == markers_mod.REGIME_B:
            members = tier2.load_panel(tier2.panel_path(panel_dir, meta.get("kalamari_taxid")))
            plan.panels[pick.unit_id] = members
            for m in members:
                if m.accession and m.accession != pick.parent_assembly:
                    plan.ani_pairs.append((pick.parent_assembly, m.accession))
                    plan.labels[m.accession] = m.lineage_label or m.unit_id

        for task_id in tasks:
            plan.typer_jobs.extend(typers_mod.jobs_for_task(task_id, pick.parent_assembly))

    plan.ani_pairs = sorted(set(plan.ani_pairs))
    plan.typer_jobs = sorted(set(plan.typer_jobs))
    return plan


def uncached(plan: "Plan", ani_cache, typer_cache):
    """The planned ANI pairs and typer jobs the committed caches cannot answer."""
    pairs = [pair for pair in plan.ani_pairs if not ani_cache.has(*pair)]
    jobs = [job for job in plan.typer_jobs if typer_cache.get(*job) is None]
    return pairs, jobs


def check_typer_environments(plan: "Plan", typer_cache, ani_cache, args) -> None:
    """Refuse a refresh whose typers are not in their pinned environments.

    Only the tools with UNCACHED work are checked.  A tool whose jobs are all cached
    runs nothing this month, so its 3-GB environment does not have to be present just
    to replay the committed calls -- and demanding it would break the ordinary monthly
    refresh on any runner that holds one shard.

    Same shape as ``refresh_tier2.py``'s shard gate, deliberately: an operator who has
    seen one of these messages should recognise the other.
    """
    _pairs, pending = uncached(plan, ani_cache, typer_cache)
    tools = sorted({tool for tool, _ in pending})
    if not tools:
        log("Typers: all %d planned job(s) are already cached; no typer will run."
            % len(plan.typer_jobs))
        return
    problems = typers_mod.preflight(tools, args.typer_env_dir)
    if not problems:
        return
    if not args.allow_unpinned:
        raise SystemExit(
            "preflight failed -- refusing to run %d typer job(s) with unpinned tools:"
            "\n  %s\nInstall the pinned environments (see bin/curation-triage/"
            "setup_envs.sh) and pass --typer-env-dir, or --skip-typers to keep the "
            "cached calls, or --allow-unpinned to run anyway."
            % (len(pending), "\n  ".join(problems)))
    log("WARNING: preflight problems, continuing under --allow-unpinned:\n  %s"
        % "\n  ".join(problems))


def check_offline_coverage(plan: "Plan", ani_cache, typer_cache, cache_dir: str) -> None:
    """Fail loudly when an offline run's caches do not cover the plan.

    Without this an offline run against a missing or empty cache builds zero leads,
    appends zero ledger records and exits 0 -- the monthly job goes green while
    reporting nothing.  A curation tool that silently reports nothing is worse than
    one that fails, because the silence is indistinguishable from "no findings".

    The same shape as ``tier1_metadata.py``'s datasets gate, deliberately: an operator
    who has seen one of these messages should recognise the other.
    """
    pairs, jobs = uncached(plan, ani_cache, typer_cache)
    if not (pairs or jobs):
        return
    raise SystemExit(
        "--offline: the committed Tier-2 cache directory %s does not cover %d of %d ANI "
        "comparison(s) (e.g. %s) and %d of %d typer job(s) (e.g. %s). Rebuild without "
        "--offline, or pass --allow-uncached to evaluate only what is cached."
        % (cache_dir, len(pairs), len(plan.ani_pairs),
           ["%s|%s" % p for p in pairs[:5]] or "none",
           len(jobs), len(plan.typer_jobs),
           ["%s|%s" % j for j in jobs[:5]] or "none"))


# --------------------------------------------------------------------------- #
#  Evaluation                                                                  #
# --------------------------------------------------------------------------- #
def evaluate(triggers, plan, ani_cache, typer_cache, manifest, args, survivors=frozenset()):
    """Run the pure checks over the cached results and build the Tier-2 leads."""
    leads, resolutions, unevaluated = [], [], []
    for pick, _reasons, meta in triggers:
        unit_id = pick.unit_id
        entry = plan.type_strain.get(unit_id) or {}
        ts_acc = entry.get("accession", "")
        ani_result, conf = None, None
        has_evidence = False
        if entry.get("is_self_type"):
            # No comparison to make, and none needed: being the type strain is the
            # strongest identity confirmation available.
            conf = tier2.confirm_self_type(unit_id, pick.species_name, ts_acc)
            has_evidence = True
        elif entry.get("has_no_type_material"):
            conf = tier2.confirm_no_type_material(unit_id, pick.species_name)
            has_evidence = True
        elif ts_acc:
            ani_result = ani_cache.get(pick.parent_assembly, ts_acc,
                                       plan.labels.get(ts_acc, ""))
            conf = tier2.confirm_species(unit_id, pick.species_name, ani_result,
                                         species_ani=args.species_ani,
                                         af_gate=args.af_gate, epsilon=args.epsilon)
            has_evidence = ani_result is not None

        # --- panel (regime B) ------------------------------------------------
        placement = None
        members = plan.panels.get(unit_id) or []
        if members:
            results = []
            for m in members:
                res = ani_cache.get(pick.parent_assembly, m.accession,
                                    m.lineage_label or m.unit_id)
                if res is not None:
                    results.append(res)
            if results:
                has_evidence = True
                external = any(m.role != "sibling" for m in members)
                placement = tier2.place_against_panel(
                    unit_id, pick.scientific_name, results, derep=args.derep,
                    epsilon=args.epsilon, af_gate=args.af_gate, self_panel=not external)

        # --- markers ---------------------------------------------------------
        calls, pins_by_task, context_by_task = [], {}, {}
        for task_id in plan.tasks.get(unit_id, []):
            evidence = typers_mod.evidence_for_task(task_id, pick.parent_assembly, typer_cache)
            if not evidence:
                continue
            has_evidence = True
            call = markers_mod.resolve(task_id, evidence, pick.species_name,
                                       pick.scientific_name)
            call.overrides_ani = (pick.regime == markers_mod.REGIME_A)
            calls.append(call)
            pins_by_task[task_id] = markers_mod.task_pins(manifest, task_id)
            # Background context (ECTyper today): attached to the lead the SME already
            # reads, never routed as a task of its own.
            context_by_task[task_id] = typers_mod.context_for_task(
                task_id, pick.parent_assembly, typer_cache)

        if not has_evidence:
            unevaluated.append(unit_id)
            if not args.emit_unevaluated:
                continue

        res = tier2.resolve_unit(unit_id, pick.regime, ani=conf, panel=placement,
                                 marker_calls=calls)
        resolutions.append(res)

        # --- leads -----------------------------------------------------------
        # Cross-link to the Tier-1 lead this confirms.  A J1a lead's key IS its
        # unit_id, so the SME can line the two rows up; blank for a unit Tier 1
        # never flagged (the mandatory confounded ones).
        confirms = unit_id if unit_id in survivors else ""
        for call in calls:
            pins = pins_by_task[call.task_id]
            leads.append(tier2.marker_lead(
                pick, call, tier2.marker_epoch(pins), pins,
                ani=conf if res.overrode_ani else None, confirms_lead=confirms,
                context=context_by_task.get(call.task_id)))
        if placement is not None:
            leads.append(tier2.panel_lead(pick, placement, ani_cache.epoch(), confirms))
        # Regime A gets NO independent ANI verdict: the marker call is the verdict and
        # the ANI number already rides along inside the marker lead as context.
        if conf is not None and not res.overrode_ani and (
                ani_result is not None or conf.gate in (tier2.GATE_SELF_TYPE,
                                                        tier2.GATE_NO_TYPE)):
            leads.append(tier2.ani_lead(pick, conf, ani_cache.epoch(), confirms))

    return leads, resolutions, unevaluated


# --------------------------------------------------------------------------- #
#  Reporting                                                                   #
# --------------------------------------------------------------------------- #
def tool_versions(ani_cache, typer_cache, manifest) -> dict:
    """Every pinned tool/DB/rule version stamped into this run's ledger records."""
    rules = sorted({r.interpretation_rule_version for r in manifest
                    if r.interpretation_rule_version})
    return {
        "skani": ani_cache.epoch(),
        "typers": dict(typer_cache.meta.get("tool_versions") or {}),
        "marker_databases": sorted({"%s=%s" % (r.database, r.database_version)
                                    for r in manifest if r.database}),
        "interpretation_rules": rules,
    }


def print_summary(part, resolutions, unevaluated, plan, run_id) -> None:
    needs, confirms = tier2.split_by_verdict(
        list(part.new) + list(part.accepted) + list(part.standing))
    new_needs, _ = tier2.split_by_verdict(part.new)
    print()
    print("=" * 78)
    print("Tier 2 confirmation (run %s)" % run_id)
    print("=" * 78)
    print("units evaluated      : %d" % len(resolutions))
    print("units with no data   : %d   (no cached ANI or typer result yet)" % len(unevaluated))
    print("ANI comparisons      : %d planned" % len(plan.ani_pairs))
    print("typer jobs           : %d planned" % len(plan.typer_jobs))
    print("-" * 78)
    print("findings to review   : %d   (%d new this run)" % (len(needs), len(new_needs)))
    print("confirmations        : %d   (annotate the Tier-1 row, not work)" % len(confirms))
    print("dismissed            : %d   (reject/defer in reviews.tsv)" % len(part.dismissed))
    print("ledger records added : %d" % len(part.records))
    overrides = [r for r in resolutions if r.overrode_ani]
    if overrides:
        print("-" * 78)
        print("regime-A marker overrides (ANI is context only for these): %d" % len(overrides))
        for r in overrides:
            print("  - %s -> %s" % (r.unit_id, r.verdict))
    if new_needs:
        print("-" * 78)
        print("NEW Tier-2 findings:")
        for lead in new_needs:
            print("  - %s" % lead.summary)
    print("=" * 78)


# --------------------------------------------------------------------------- #
#  Driver                                                                      #
# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Tier 2 targeted confirmation.")
    p.add_argument("--policy", default=DEFAULT_POLICY)
    p.add_argument("--chromosomes", default=DEFAULT_CHROMOSOMES)
    p.add_argument("--confounded-manifest", default=DEFAULT_CONFOUNDED)
    p.add_argument("--marker-manifest", default=DEFAULT_MANIFEST)
    p.add_argument("--panel-dir", default=DEFAULT_PANEL_DIR)
    p.add_argument("--ledger", default=DEFAULT_LEDGER)
    p.add_argument("--reviews", default=DEFAULT_REVIEWS)
    p.add_argument("--out", default=DEFAULT_WORKLIST, help="Per-run Tier-2 worklist TSV.")
    p.add_argument("--summary-md", default=None,
                   help="Write the $GITHUB_STEP_SUMMARY markdown block to this file.")
    p.add_argument("--run-id", default=None, help="Run tag (default: current YYYY-MM).")
    # Caches: one directory holds all three committed result caches.
    p.add_argument("--cache-dir", default=DEFAULT_CACHE_DIR,
                   help="Directory holding the committed ANI / typer / type-strain caches.")
    p.add_argument("--genome-cache", default=DEFAULT_GENOME_CACHE,
                   help="Gitignored FASTA cache (sequence is never committed).")
    p.add_argument("--datasets-bin", default=genomes_mod.DEFAULT_DATASETS_BIN,
                   help="Path to the ncbi-datasets-cli binary.")
    p.add_argument("--skani-bin", default=skani_runner.DEFAULT_SKANI_BIN,
                   help="Path to the skani binary.")
    p.add_argument("--typer-env-dir", default=None,
                   help="Directory of per-tool conda envs (<dir>/<tool>/bin/...); "
                        "without it each typer must be on PATH.")
    p.add_argument("--allow-unpinned", action="store_true",
                   help="Downgrade the typer preflight failure to a warning. Without "
                        "it a refresh whose pinned environments are missing stops "
                        "before it runs, rather than recording every job as a miss.")
    p.add_argument("--offline", action="store_true",
                   help="Never fetch or run a tool; evaluate only what is cached.")
    p.add_argument("--allow-uncached", action="store_true",
                   help="With --offline, continue when the caches do not cover the whole "
                        "plan (for a deliberate partial run). Without it an uncovered "
                        "plan is an error, because 'no findings' and 'no data' must not "
                        "look the same.")
    p.add_argument("--no-write", action="store_true",
                   help="Compute and summarise, but do not append the ledger or write files.")
    # Scope.
    p.add_argument("--only-confounded", action="store_true",
                   help="Skip the Tier-1 survivors and confirm only the confounded units.")
    p.add_argument("--only-survivors", action="store_true",
                   help="Skip the mandatory confounded units and confirm only Tier-1 leads.")
    # Staged refresh: skani is cheap and the typers are not, so they can be rebuilt
    # independently instead of forcing an all-or-nothing refresh.
    p.add_argument("--skip-typers", action="store_true",
                   help="Do not run the typers when refreshing (keep the cached calls).")
    p.add_argument("--skip-ani", action="store_true",
                   help="Do not run skani when refreshing (keep the cached numbers).")
    p.add_argument("--emit-unevaluated", action="store_true",
                   help="Also emit a lead for units with no Tier-2 evidence yet "
                        "(off by default: 'not measured' is not a finding).")
    # Thresholds (design section 10 defaults).
    p.add_argument("--species-ani", type=float, default=tier2.DEFAULT_SPECIES_ANI)
    p.add_argument("--af-gate", type=float, default=tier2.DEFAULT_AF_GATE)
    p.add_argument("--derep", type=float, default=tier2.DEFAULT_DEREP)
    p.add_argument("--epsilon", type=float, default=tier2.DEFAULT_EPSILON)
    args = p.parse_args(argv)

    run_id = args.run_id or datetime.date.today().strftime("%Y-%m")
    ani_path = os.path.join(args.cache_dir, ANI_CACHE_NAME)
    typer_path = os.path.join(args.cache_dir, TYPER_CACHE_NAME)
    ts_path = os.path.join(args.cache_dir, TYPE_STRAIN_CACHE_NAME)

    # --- inputs -------------------------------------------------------------
    source_lanes = coverage.load_source_lanes(args.chromosomes)
    picks = coverage.load_picks(args.policy, source_lanes)
    unit_meta = load_unit_meta(args.policy)
    confounded_entries = policy_mod.load_confounded_manifest(args.confounded_manifest)
    manifest = markers_mod.load_manifest(args.marker_manifest)
    prior = ledger_mod.load_ledger(args.ledger)
    reviews = ledger_mod.load_reviews(args.reviews)

    survivors = survivor_unit_ids(prior)
    triggers = trigger_units(picks, survivors, unit_meta,
                             want_survivors=not args.only_confounded,
                             want_confounded=not args.only_survivors)
    log("Tier 2 triggers: %d unit(s) -- %d confounded, %d Tier-1 survivor(s)"
        % (len(triggers),
           sum(1 for p, _, _ in triggers if p.is_confounded),
           sum(1 for p, _, _ in triggers if p.unit_id in survivors)))

    ts_cache = genomes_mod.TypeStrainCache(cache_path=ts_path)
    plan = build_plan(triggers, manifest, confounded_entries, ts_cache, args.panel_dir)
    log("Plan: %d ANI pair(s), %d typer job(s), %d genome(s)"
        % (len(plan.ani_pairs), len(plan.typer_jobs), len(plan.accessions)))

    # --- external work (skipped entirely when offline) ----------------------
    ani_cache = skani_runner.AniCache(cache_path=ani_path, skani_bin=args.skani_bin)
    typer_cache = typers_mod.TyperCache(cache_path=typer_path, env_dir=args.typer_env_dir)
    if not args.offline:
        # Check the typer environments BEFORE anything downloads or runs, so a wrong
        # --typer-env-dir stops here instead of turning into a cache full of holes.
        # typers.tool_binary() falls back to a bare-name PATH lookup when the pinned
        # environment is absent, so a missing environment does not fail: every job
        # becomes a miss and the run still exits 0.  Checking first also means a
        # failure costs no skani work, since nothing has been computed yet.
        if not args.skip_typers:
            check_typer_environments(plan, typer_cache, ani_cache, args)
        store = genomes_mod.GenomeStore(args.genome_cache, offline=False,
                                        datasets_bin=args.datasets_bin)
        store.prefetch(plan.accessions, log=log)      # one bulk call, not N
        if not args.skip_ani:
            ani_cache.build(plan.ani_pairs, store, labels=plan.labels, log=log)
        if not args.skip_typers:
            typer_cache.build(plan.typer_jobs, store, args.genome_cache, log=log)
        if not args.no_write:
            if not args.skip_ani:
                ani_cache.save()
            if not args.skip_typers:
                typer_cache.save()
        if ani_cache.misses or typer_cache.misses:
            log("WARNING: %d genome(s) and %d typer job(s) could not be completed; "
                "their units stay unevaluated rather than being guessed."
                % (len(set(ani_cache.misses)), len(typer_cache.misses)))
    elif not args.allow_uncached:
        check_offline_coverage(plan, ani_cache, typer_cache, args.cache_dir)
        log("--offline: the committed caches cover all %d ANI comparison(s) and all %d "
            "typer job(s)." % (len(plan.ani_pairs), len(plan.typer_jobs)))
    else:
        missing_ani, missing_jobs = uncached(plan, ani_cache, typer_cache)
        if missing_ani or missing_jobs:
            log("--allow-uncached: %d of %d ANI comparison(s) and %d of %d typer job(s) "
                "are not cached; those checks are reported as unevaluated."
                % (len(missing_ani), len(plan.ani_pairs),
                   len(missing_jobs), len(plan.typer_jobs)))

    # --- pure evaluation ----------------------------------------------------
    leads, resolutions, unevaluated = evaluate(triggers, plan, ani_cache, typer_cache,
                                               manifest, args, survivors=survivors)
    log("Built %d Tier-2 lead(s) over %d evaluated unit(s); %d unit(s) had no evidence"
        % (len(leads), len(resolutions), len(unevaluated)))

    part = ledger_mod.partition(leads, prior, reviews, run_id=run_id,
                                params={"species_ani": args.species_ani,
                                        "af_gate": args.af_gate,
                                        "derep": args.derep,
                                        "epsilon": args.epsilon},
                                tool_versions=tool_versions(ani_cache, typer_cache, manifest),
                                # Tier 2 recomputes no Tier-1 signal, so it is entitled
                                # to mark nothing resolved: the absence of a J1a lead
                                # here means "not my job", not "fixed".
                                resolve_signals=())

    if not args.no_write:
        ledger_mod.append_ledger(args.ledger, part.records)
        log("Appended %d ledger records: %s" % (len(part.records), args.ledger))
        ledger_mod.write_worklist(part, args.out)
        log("Wrote Tier-2 worklist: %s" % args.out)
        md = tier2.render_markdown(part, run_id)
        if args.summary_md:
            os.makedirs(os.path.dirname(os.path.abspath(args.summary_md)), exist_ok=True)
            with open(args.summary_md, "w", encoding="utf-8") as fh:
                fh.write(md + "\n")
        gh = os.environ.get("GITHUB_STEP_SUMMARY")
        if gh:
            with open(gh, "a", encoding="utf-8") as fh:
                fh.write(md + "\n")

    print_summary(part, resolutions, unevaluated, plan, run_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
