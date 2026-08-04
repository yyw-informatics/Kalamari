# -*- coding: utf-8 -*-
"""Tier 2 end to end, offline, against the REAL committed caches.

Nothing here is synthetic.  The ANI numbers come from a real skani 0.3.2 run, the
marker calls from real runs of the pinned typers, and the type strains from NCBI's
own ANI report -- all committed under ``src/curation-triage/``.  The test therefore
exercises exactly what ships, and a cache refresh that changes a verdict will show up
here as a failing assertion, which for a curation tool is the point.
"""

import json
import os
import shutil

import pytest

import ledger as ledger_mod
import markers
import tier2
import tier2_confirm

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")

BACILLUS_CEREUS = "Bacillus_cereus~GCF_000283675.1"
BACILLUS_THUR = "Bacillus_thuringiensis~GCF_000008505.1"
ECOLI = "Escherichia_coli~GCF_003018575.1"
YERSINIA = "Yersinia_pseudotuberculosis~GCF_000834295.1"
LISTERIA_III = "Listeria_monocytogenes_III~GCF_013282645.1"
SALMONELLA_IIIA = "Salmonella_enterica_IIIa~GCF_000018625.1"
BUCHNERA = "Buchnera_aphidicola~GCF_000009605.1"
KORARCHAEUM = "Candidatus_Korarchaeum_cryptofilum~GCF_000019605.1"


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    """One offline run over the committed caches; returns its T2 ledger records."""
    tmp = tmp_path_factory.mktemp("tier2")
    ledger_path = tmp / "ledger.ndjson"
    shutil.copyfile(os.path.join(SRC, "ledger.ndjson"), ledger_path)
    # Keep the pre-run state so a test can prove the run only APPENDED to it. Snapshot
    # it here rather than re-reading the committed file later: what must not change is
    # what this run actually read.
    shutil.copyfile(os.path.join(SRC, "ledger.ndjson"), tmp / "ledger.before.ndjson")
    reviews_path = tmp / "reviews.tsv"
    reviews_path.write_text("\t".join(ledger_mod.REVIEW_COLUMNS) + "\n", encoding="utf-8")

    rc = tier2_confirm.main([
        "--offline", "--run-id", "2026-08",
        "--ledger", str(ledger_path), "--reviews", str(reviews_path),
        "--out", str(tmp / "worklist_tier2.tsv"),
        "--summary-md", str(tmp / "worklist_tier2.md"),
    ])
    assert rc == 0
    records = [json.loads(l) for l in open(ledger_path, encoding="utf-8")]
    return {r["lead_key"]: r for r in records if r["signal"] == "T2"}, tmp


def check(t2, unit, kind):
    return t2[next(k for k in t2 if k.startswith("t2:%s#%s" % (unit, kind)))]


# --------------------------------------------------------------------------- #
#  Coverage of the trigger set                                                 #
# --------------------------------------------------------------------------- #
def test_every_triggered_unit_is_actually_evaluated(run):
    """41 units, none left unmeasured -- the caches cover the whole trigger set."""
    t2, _ = run
    units = {k[len("t2:"):].split("#")[0] for k in t2}
    assert len(units) == 41
    assert len(t2) == 75          # 37 ani + 21 marker + 17 panel checks


def test_all_three_checks_are_represented(run):
    t2, _ = run
    kinds = {}
    for rec in t2.values():
        kinds[rec["evidence"]["t2_check"]] = kinds.get(rec["evidence"]["t2_check"], 0) + 1
    assert kinds == {"ani": 37, "marker": 21, "panel": 17}


# --------------------------------------------------------------------------- #
#  The regime-A override, on real marker calls                                 #
# --------------------------------------------------------------------------- #
def test_regime_a_units_get_no_competing_ani_verdict(run):
    """The override is structural: only the marker call is a verdict."""
    t2, _ = run
    for unit in (BACILLUS_CEREUS, BACILLUS_THUR, ECOLI, YERSINIA):
        assert "t2:%s#ani" % unit not in t2
        rec = check(t2, unit, "marker")
        assert rec["evidence"]["overrides_ani"] == "yes"


def test_the_thuringiensis_pick_is_nearest_to_anthracis_but_the_markers_say_otherwise(run):
    """The case regime A exists for: ANI points at anthracis, the markers rule it out."""
    t2, _ = run
    rec = check(t2, BACILLUS_THUR, "marker")
    ev = rec["evidence"]
    assert ev["marker_closest_type_strain"] == "anthracis"
    assert ev["marker_anthracis_chromosome"] == "no"
    assert ev["marker_pxo1_markers"] == "none" and ev["marker_pxo2_markers"] == "none"
    # The ANI number is carried as context only, never as a competing verdict.
    assert "context only" in ev["ani_context"]


def test_the_cereus_pick_keeps_its_label_but_records_the_paranthracis_proximity(run):
    t2, _ = run
    ev = check(t2, BACILLUS_CEREUS, "marker")["evidence"]
    assert ev["verdict"] == tier2.VERDICT_CONFIRMED
    assert ev["marker_closest_type_strain"] == "paranthracis"
    assert "paranthracis" in ev["marker_reasons"]


def test_the_e_coli_pick_is_not_shigella(run):
    t2, _ = run
    ev = check(t2, ECOLI, "marker")["evidence"]
    assert ev["verdict"] == tier2.VERDICT_CONFIRMED
    assert ev["marker_ipaH"] == "absent"


def test_the_yersinia_task_says_undecided_because_its_classifier_is_a_gap(run):
    """No cgMLST caller exists, so the honest answer is 'undecided', not a guess."""
    t2, _ = run
    ev = check(t2, YERSINIA, "marker")["evidence"]
    assert ev["verdict"] == tier2.VERDICT_INDETERMINATE
    assert "no usable cgMLST" in ev["marker_reasons"]


# --------------------------------------------------------------------------- #
#  Regime B: panel decides membership, the typer supplies the label            #
# --------------------------------------------------------------------------- #
def test_a_regime_b_unit_gets_all_three_checks(run):
    t2, _ = run
    kinds = sorted(k.split("#", 1)[1].split(":")[0]
                   for k in t2 if k.startswith("t2:%s#" % LISTERIA_III))
    assert kinds == ["ani", "marker", "panel"]
    assert check(t2, LISTERIA_III, "marker")["evidence"]["overrides_ani"] == "no"


def test_the_salmonella_subspecies_label_is_confirmed_by_two_typers(run):
    t2, _ = run
    ev = check(t2, SALMONELLA_IIIA, "marker")["evidence"]
    assert ev["verdict"] == tier2.VERDICT_CONFIRMED
    assert ev["marker_subspecies"] == "arizonae"
    assert ev["marker_serovar_agreement"] == "agree"


def test_the_listeria_lineage_stays_undecided_without_a_cgmlst_caller(run):
    t2, _ = run
    ev = check(t2, LISTERIA_III, "marker")["evidence"]
    assert ev["verdict"] == tier2.VERDICT_INDETERMINATE
    assert ev["marker_predicted_molecular_serogroup"]        # serogroup IS reported
    assert ev["marker_broad_lineage"] == "none"              # lineage is not invented


def test_the_lineage_panels_are_all_distinct(run):
    """No two Kalamari lineage references collapsed into the same genome."""
    t2, _ = run
    panels = [r for k, r in t2.items() if r["evidence"]["t2_check"] == "panel"]
    assert len(panels) == 17
    assert not [p for p in panels if p["evidence"]["redundant_with"]]


# --------------------------------------------------------------------------- #
#  The states that only real data produced                                     #
# --------------------------------------------------------------------------- #
def test_a_species_with_no_type_material_is_explained_not_left_blank(run):
    """This is also why NCBI's own taxonomy check on these picks is inconclusive."""
    t2, _ = run
    rec = check(t2, KORARCHAEUM, "ani")
    assert rec["evidence"]["gate"] == tier2.GATE_NO_TYPE
    assert rec["actionable_flags"] == ["t2_no_type_material"]
    assert "no type material" in rec["evidence"]["reasons"]


def test_a_pair_too_distant_for_skani_is_a_finding_not_a_gap(run):
    t2, _ = run
    rec = check(t2, BUCHNERA, "ani")
    assert rec["evidence"]["gate"] == tier2.GATE_SCREENED
    assert rec["evidence"]["verdict"] == tier2.VERDICT_REFUTED


def test_a_pick_that_is_its_own_type_strain_needs_no_comparison(run):
    t2, _ = run
    selfs = [r for r in t2.values() if r["evidence"].get("gate") == tier2.GATE_SELF_TYPE]
    assert selfs, "expected at least one pick that IS its species' type strain"
    assert selfs[0]["evidence"]["verdict"] == tier2.VERDICT_CONFIRMED


def test_the_aligned_fraction_gate_fires_on_real_comparisons(run):
    t2, _ = run
    unusable = [r for r in t2.values()
                if r["evidence"].get("gate") == tier2.GATE_UNUSABLE]
    assert unusable, "expected the af gate to veto at least one real comparison"
    for rec in unusable:
        assert rec["evidence"]["verdict"] == tier2.VERDICT_INDETERMINATE


# --------------------------------------------------------------------------- #
#  Provenance and the handoff loop                                             #
# --------------------------------------------------------------------------- #
def test_every_record_pins_the_tool_the_database_and_our_rule(run):
    t2, _ = run
    rec = check(t2, BACILLUS_CEREUS, "marker")
    assert rec["epoch"].startswith("btyper3@3.4.0")
    assert "rule:bacillus_anthracis_vs_b_cereus_group@1" in rec["epoch"]
    assert rec["tool_versions"]["skani"] == "skani@0.3.2"
    assert rec["tool_versions"]["typers"]["amrfinderplus"] == "amrfinderplus@4.2.7"
    assert any("2026-05-15.1" in d for d in rec["tool_versions"]["marker_databases"])


def test_the_tier1_baseline_is_untouched_by_a_tier2_run(run):
    """Tier 2 only APPENDS T2 records: every Tier-1 record it read comes out unchanged.

    Stated as before-vs-after rather than as a Tier-1 record count. The ledger is
    append-only and the monthly bot commits a longer copy of it every run, so a pinned
    count would fail on the bot's own push -- and it would not even test the right
    thing: what matters is that a Tier-2 run never rewrites or resolves Tier-1 work.
    """
    _, tmp = run
    before = [json.loads(l) for l in open(tmp / "ledger.before.ndjson", encoding="utf-8")]
    after = [json.loads(l) for l in open(tmp / "ledger.ndjson", encoding="utf-8")]
    tier1_before = [r for r in before if r["signal"] != "T2"]
    assert tier1_before, "no Tier-1 baseline to protect -- the check would be vacuous"
    assert [r for r in after if r["signal"] != "T2"] == tier1_before

    appended = after[len(before):]
    assert after[:len(before)] == before, "the run rewrote history instead of appending"
    assert appended and {r["signal"] for r in appended} == {"T2"}
    # A Tier-1 lead is never re-labelled -- in particular never marked resolved -- by a
    # tier that did not recompute the Tier-1 signals.
    assert not {r["lead_key"] for r in appended} & {r["lead_key"] for r in tier1_before}


def test_confirmations_are_recorded_but_do_not_become_work(run):
    t2, _ = run
    confirmed = [r for r in t2.values()
                 if r["evidence"]["verdict"] == tier2.VERDICT_CONFIRMED]
    assert len(confirmed) == 42
    md = open(os.path.join(str(run[1]), "worklist_tier2.md"), encoding="utf-8").read()
    assert "**33 Tier-2 finding(s) need review**" in md
    assert "Confirmations (42)" in md


def test_marker_calls_carry_their_interpretation_rule_version(run):
    """Every marker lead has a real rule; a context-only task never becomes a lead."""
    t2, _ = run
    for rec in t2.values():
        if rec["evidence"]["t2_check"] == tier2.CHECK_MARKER:
            assert rec["evidence"]["interpretation_rule_version"].endswith("@1")
    assert not any("escherichia_coli_serotype" in k for k in t2)


def test_context_only_tasks_are_not_routed_to_the_lead_path(run):
    """ECTyper is background for an E. coli pick, never a verdict about it."""
    manifest = markers.load_manifest(os.path.join(SRC, "marker_manifest.tsv"))
    assert "escherichia_coli_serotype_pathotype" not in markers.tasks_for_taxids(
        manifest, [562])
    assert "escherichia_coli_serotype_pathotype" in markers.tasks_for_taxids(
        manifest, [562], decisive_only=False)


def test_the_e_coli_lead_carries_the_real_ectyper_serotype_as_labelled_context(run):
    """The cached ECTyper result is shown to the SME without becoming a verdict."""
    t2, _ = run
    ev = check(t2, ECOLI, "marker")["evidence"]
    assert ev["ectyper_serotype"] == "O118/O151:H16"
    assert ev["ectyper_qc"] == "WARNING MIXED O-TYPE"
    assert ev["ectyper_db_version"].startswith("v1.0")
    assert "context" in ev["ectyper_role"]
    # ShigEiFinder's own (empty) serotype field is untouched: the two never merge.
    assert ev["marker_serotype"] == ""
    assert ev["verdict"] == tier2.VERDICT_CONFIRMED


def test_context_reaches_exactly_one_lead_and_adds_no_new_one(run):
    """Attaching context must not move any of the counts above."""
    t2, _ = run
    with_context = [k for k, r in t2.items() if "ectyper_serotype" in r["evidence"]]
    assert len(with_context) == 1
    assert with_context[0].startswith("t2:%s#marker" % ECOLI)


# --------------------------------------------------------------------------- #
#  Offline strictness                                                          #
# --------------------------------------------------------------------------- #
def offline_run(tmp_path, *extra):
    ledger_path = tmp_path / "ledger.ndjson"
    shutil.copyfile(os.path.join(SRC, "ledger.ndjson"), ledger_path)
    return tier2_confirm.main([
        "--offline", "--run-id", "2026-08", "--no-write",
        "--cache-dir", str(tmp_path / "empty-caches"),
        "--ledger", str(ledger_path), "--reviews", os.path.join(SRC, "reviews.tsv"),
    ] + list(extra))


def test_an_offline_run_the_caches_do_not_cover_fails_instead_of_reporting_nothing(tmp_path):
    """It used to build 0 leads, append 0 records and exit 0 -- silently green."""
    (tmp_path / "empty-caches").mkdir()
    with pytest.raises(SystemExit) as exc:
        offline_run(tmp_path)
    message = str(exc.value)
    assert "does not cover" in message
    assert "Rebuild without --offline" in message
    # It names what is missing, so the operator knows which shard to re-run.
    assert "sistr|" in message or "abricate|" in message


def test_a_deliberate_partial_offline_run_is_still_possible(tmp_path):
    (tmp_path / "empty-caches").mkdir()
    assert offline_run(tmp_path, "--allow-uncached") == 0


def test_the_committed_caches_pass_the_coverage_gate(run):
    """The gate must not fire on what ships: the module-scope run above is the proof."""
    t2, _ = run
    assert len(t2) == 75


# --------------------------------------------------------------------------- #
#  Trigger-set helpers                                                         #
# --------------------------------------------------------------------------- #
def test_a_resolved_tier1_lead_is_not_re_confirmed():
    resolved = [{"lead_key": "u", "unit_id": "u", "signal": "J1a", "delta": "resolved",
                 "actionable_flags": []}]
    assert tier2_confirm.survivor_unit_ids(resolved) == set()


def test_tier2_does_not_feed_itself():
    own = [{"lead_key": "t2:u#ani", "unit_id": "u", "signal": "T2", "delta": "new",
            "actionable_flags": ["t2_ani_refutes_species"]}]
    assert tier2_confirm.survivor_unit_ids(own) == set()


def test_the_type_strain_is_keyed_per_assembly_not_per_species():
    """All eleven Salmonella units share species 28901 but need different type strains."""
    import genomes
    cache = genomes.TypeStrainCache(os.path.join(SRC, "type_strain_cache.json"))
    iiia = cache.accession("GCF_000018625", 28901)
    iiib = cache.accession("GCA_013162305", 28901)
    assert iiia and iiib and iiia != iiib
