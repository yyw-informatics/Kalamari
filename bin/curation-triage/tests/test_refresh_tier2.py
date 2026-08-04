# -*- coding: utf-8 -*-
"""The sharded Tier-2 refresh -- accession list, fragments, merge.

The merge is the part worth guarding.  It is the last step before a refreshed cache
is committed, and every Tier-2 verdict is computed from that cache, so a merge that
quietly picked one of two conflicting numbers would produce findings nobody can
attribute.  Everything here runs offline over the committed caches; no tool and no
network is involved.
"""

import argparse
import json
import os

import pytest

import refresh_tier2
import tier2_confirm

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")
ANI_CACHE = os.path.join(SRC, "tier2_ani_cache.json")
TYPER_CACHE = os.path.join(SRC, "tier2_typer_cache.json")


def load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def frag(results, **meta):
    return {"schema": 1, "meta": dict(meta), "results": dict(results)}


def write(tmp_path, name, blob):
    path = tmp_path / name
    path.write_text(json.dumps(blob), encoding="utf-8")
    return str(path)


# --------------------------------------------------------------------------- #
#  The accession list (the CI genome-cache key)                                #
# --------------------------------------------------------------------------- #
def test_the_accession_list_is_sorted_deduplicated_and_covers_both_caches(capsys):
    """CI keys its genome cache on this list, so a missing accession is a cache miss."""
    assert refresh_tier2.main(["--print-accessions"]) == 0
    printed = capsys.readouterr().out.split()
    assert printed == sorted(set(printed))
    typer_accessions = {k.split("|", 1)[1] for k in load(TYPER_CACHE)["results"]}
    ani_accessions = set()
    for key in load(ANI_CACHE)["results"]:
        ani_accessions.update(key.split("|"))
    assert typer_accessions <= set(printed)
    assert ani_accessions <= set(printed)


def plan_args(**over):
    """The default plan inputs, as the CLI would resolve them."""
    args = dict(policy=tier2_confirm.DEFAULT_POLICY,
                chromosomes=tier2_confirm.DEFAULT_CHROMOSOMES,
                confounded_manifest=tier2_confirm.DEFAULT_CONFOUNDED,
                marker_manifest=tier2_confirm.DEFAULT_MANIFEST,
                panel_dir=tier2_confirm.DEFAULT_PANEL_DIR,
                ledger=tier2_confirm.DEFAULT_LEDGER,
                cache_dir=tier2_confirm.DEFAULT_CACHE_DIR, with_ectyper=False)
    args.update(over)
    return argparse.Namespace(**args)


def test_the_plan_matches_what_a_real_tier2_run_computes():
    """A shard that plans differently from the monthly run refreshes the wrong things."""
    plan = refresh_tier2.build_plan(plan_args())
    assert len(plan.ani_pairs) == 158
    assert len(plan.typer_jobs) == 37
    assert len(plan.accessions) == 64


def test_planning_ectyper_is_opt_in_because_of_its_943_mib_download():
    plan = refresh_tier2.build_plan(plan_args(with_ectyper=True))
    assert ("ectyper", "GCF_003018575.1") in plan.typer_jobs
    assert len(plan.typer_jobs) == 38
    # ... and it adds no genome: ECTyper types the E. coli pick we already download.
    assert len(plan.accessions) == 64


# --------------------------------------------------------------------------- #
#  The pure merge                                                              #
# --------------------------------------------------------------------------- #
def test_merging_fragments_unions_the_results_and_the_tool_versions():
    merged = refresh_tier2.merge_caches([
        ("a", frag({"sistr|GCF_1.1": {"subspecies": "arizonae"}},
                   tool_versions={"sistr": "sistr@1.1.3"}, n_results=1)),
        ("b", frag({"abricate|GCF_2.1": {"abricate_hit_count": "3"}},
                   tool_versions={"abricate": "abricate@1.4.0"}, n_results=1)),
    ])
    assert sorted(merged["results"]) == ["abricate|GCF_2.1", "sistr|GCF_1.1"]
    assert merged["meta"]["tool_versions"] == {"abricate": "abricate@1.4.0",
                                               "sistr": "sistr@1.1.3"}


def test_the_counts_are_recomputed_from_the_merged_results_not_carried_over():
    """A fragment's n_results describes the fragment; the merged cache is bigger."""
    merged = refresh_tier2.merge_caches([
        ("a", frag({"x|1": {"ani": 97.0, "screened_out": False}}, n_results=1,
                   n_screened_out=0)),
        ("b", frag({"y|2": {"ani": None, "screened_out": True}}, n_results=1,
                   n_screened_out=1)),
    ])
    assert merged["meta"]["n_results"] == 2
    assert merged["meta"]["n_screened_out"] == 1


def test_a_conflicting_result_between_two_fragments_is_an_error_not_an_overwrite():
    """Last-writer-wins here would let a stale shard replace a good number."""
    with pytest.raises(refresh_tier2.MergeConflict) as exc:
        refresh_tier2.merge_caches([
            ("good", frag({"x|1": {"ani": 97.0}})),
            ("stale", frag({"x|1": {"ani": 91.0}})),
        ])
    assert "x|1" in str(exc.value) and "stale" in str(exc.value)


def test_two_fragments_that_disagree_about_a_tool_version_are_an_error():
    """A tool version IS the epoch, so two of them in one cache is unattributable."""
    with pytest.raises(refresh_tier2.MergeConflict):
        refresh_tier2.merge_caches([
            ("a", frag({}, tool_versions={"sistr": "sistr@1.1.3"})),
            ("b", frag({}, tool_versions={"sistr": "sistr@1.1.2"})),
        ])


def test_two_amrfinder_database_versions_cannot_be_merged_into_one_cache():
    with pytest.raises(refresh_tier2.MergeConflict) as exc:
        refresh_tier2.merge_caches([
            ("a", frag({}, amrfinder_database_version="2026-05-15.1")),
            ("b", frag({}, amrfinder_database_version="2026-07-01.1")),
        ])
    assert "amrfinder_database_version" in str(exc.value)


def test_an_identical_value_in_two_fragments_is_not_a_conflict():
    """Overlapping shards are allowed; only DISagreement is refused."""
    merged = refresh_tier2.merge_caches([
        ("a", frag({"x|1": {"ani": 97.0}}, skani_version="skani@0.3.2")),
        ("b", frag({"x|1": {"ani": 97.0}}, skani_version="skani@0.3.2")),
    ])
    assert merged["meta"]["skani_version"] == "skani@0.3.2"
    assert merged["meta"]["n_results"] == 1


def test_a_mixed_schema_version_is_refused():
    with pytest.raises(refresh_tier2.MergeConflict):
        refresh_tier2.merge_caches([("a", {"schema": 1, "meta": {}, "results": {}}),
                                    ("b", {"schema": 2, "meta": {}, "results": {}})])


def test_the_cache_kind_is_told_apart_from_the_committed_caches():
    assert refresh_tier2.cache_kind(load(ANI_CACHE)) == refresh_tier2.KIND_ANI
    assert refresh_tier2.cache_kind(load(TYPER_CACHE)) == refresh_tier2.KIND_TYPER
    assert refresh_tier2.cache_kind({"meta": {}, "results": {}}) == ""


# --------------------------------------------------------------------------- #
#  Coverage of the plan                                                        #
# --------------------------------------------------------------------------- #
def test_the_committed_caches_merge_into_themselves_unchanged(tmp_path):
    """Round-trip on real data: merge is a no-op on a cache that is already whole."""
    for source in (ANI_CACHE, TYPER_CACHE):
        out = str(tmp_path / os.path.basename(source))
        assert refresh_tier2.main(["--merge", source, "--out", out]) == 0
        assert load(out) == load(source)


def test_a_merge_that_misses_a_planned_job_is_refused(tmp_path):
    blob = load(TYPER_CACHE)
    blob["results"].pop("sistr|GCF_000018625.1")
    path = write(tmp_path, "partial.json", blob)
    with pytest.raises(SystemExit) as exc:
        refresh_tier2.main(["--merge", path, "--out", str(tmp_path / "out.json")])
    assert "sistr|GCF_000018625.1" in str(exc.value)
    assert "--allow-partial" in str(exc.value)


def test_a_merge_that_misses_a_planned_ani_pair_is_refused(tmp_path):
    blob = load(ANI_CACHE)
    victim = sorted(blob["results"])[0]
    blob["results"].pop(victim)
    # The reverse comparison answers the same question, so drop it too if it is there.
    blob["results"].pop("|".join(reversed(victim.split("|"))), None)
    path = write(tmp_path, "partial.json", blob)
    with pytest.raises(SystemExit) as exc:
        refresh_tier2.main(["--merge", path, "--out", str(tmp_path / "out.json")])
    assert "ANI comparison(s)" in str(exc.value)


def test_a_tool_with_results_but_no_recorded_version_is_still_refused(tmp_path):
    """Complete results with a missing epoch is exactly the shard bug this guards."""
    blob = load(TYPER_CACHE)
    blob["meta"]["tool_versions"].pop("sistr")
    path = write(tmp_path, "no-version.json", blob)
    with pytest.raises(SystemExit) as exc:
        refresh_tier2.main(["--merge", path, "--out", str(tmp_path / "out.json")])
    assert "no recorded tool version" in str(exc.value)


def test_a_conflict_reaches_the_operator_as_a_message_not_a_traceback(tmp_path):
    good = load(TYPER_CACHE)
    stale = load(TYPER_CACHE)
    stale["results"]["abricate|GCF_000834295.1"]["abricate_hit_count"] = "999"
    with pytest.raises(SystemExit) as exc:
        refresh_tier2.main(["--merge", write(tmp_path, "good.json", good),
                            write(tmp_path, "stale.json", stale),
                            "--out", str(tmp_path / "out.json")])
    assert "merge refused" in str(exc.value)
    assert "abricate|GCF_000834295.1" in str(exc.value)
    assert not (tmp_path / "out.json").exists()


# --------------------------------------------------------------------------- #
#  Keeping the cache whole (a stronger promise than covering the plan)          #
# --------------------------------------------------------------------------- #
def nine_shard_fragments(tmp_path):
    """The fragments a FULL typer refresh produces -- one per tool the plan routes.

    Built by splitting the committed cache by tool, which is what the shards recompute.
    There is no ectyper fragment, and that is the point: ECTyper is context-only, the
    plan lists it only behind --with-ectyper, and the refresh workflow's matrix has no
    entry for it -- so the committed ectyper result is a key no shard ever produces.
    """
    committed = load(TYPER_CACHE)
    plan = refresh_tier2.build_plan(plan_args())
    paths = []
    for tool in sorted({t for t, _ in plan.typer_jobs}):
        results = {k: v for k, v in committed["results"].items()
                   if k.split("|", 1)[0] == tool}
        meta = {"tool_versions": {tool: committed["meta"]["tool_versions"][tool]},
                "n_results": len(results)}
        if tool == "amrfinderplus":
            meta["amrfinder_database_version"] = \
                committed["meta"]["amrfinder_database_version"]
        paths.append(write(tmp_path, "frag-typer-%s.json" % tool,
                           {"schema": 1, "meta": meta, "results": results}))
    return paths


def test_a_full_typer_refresh_that_would_delete_the_committed_ectyper_result_is_refused(tmp_path):
    """The regression: 37 of 37 planned jobs covered, and the cache still got smaller.

    Nothing downstream could see it -- context keys sit outside the Tier-2 fingerprint,
    so the delta stayed `recurring` and both the coverage replay and the update gate
    passed on a cache that had lost a result and a tool version.
    """
    frags = nine_shard_fragments(tmp_path)
    assert len(frags) == 9
    out = write(tmp_path, "tier2_typer_cache.json", load(TYPER_CACHE))
    with pytest.raises(SystemExit) as exc:
        refresh_tier2.main(["--merge"] + frags + ["--out", out])
    assert "ectyper|GCF_003018575.1" in str(exc.value)
    assert "meta.tool_versions.ectyper" in str(exc.value)
    assert "--allow-drop" in str(exc.value)
    # ... and the cache it refused to shrink is exactly as it was.
    assert load(out) == load(TYPER_CACHE)


def test_the_replaced_cache_as_a_base_carries_the_context_only_ectyper_result_forward(tmp_path):
    frags = nine_shard_fragments(tmp_path)
    out = write(tmp_path, "tier2_typer_cache.json", load(TYPER_CACHE))
    assert refresh_tier2.main(["--merge"] + frags
                              + ["--base", TYPER_CACHE, "--out", out]) == 0
    assert load(out) == load(TYPER_CACHE)


def test_allow_drop_is_the_only_way_to_shrink_a_cache_and_it_names_what_it_dropped(tmp_path, capsys):
    frags = nine_shard_fragments(tmp_path)
    out = write(tmp_path, "tier2_typer_cache.json", load(TYPER_CACHE))
    assert refresh_tier2.main(["--merge"] + frags + ["--out", out, "--allow-drop"]) == 0
    assert "ectyper|GCF_003018575.1" not in load(out)["results"]
    warning = capsys.readouterr().err
    assert "--allow-drop" in warning and "ectyper|GCF_003018575.1" in warning


def test_a_merge_that_would_shrink_the_tool_versions_of_the_cache_it_overwrites_is_refused(tmp_path):
    """A tool version IS an epoch: losing one stops a dismissed lead re-opening.

    ectyper is used here because it is not a planned tool, so the plan-coverage check
    has nothing to say about it -- only the shrink guard can catch this.
    """
    blob = load(TYPER_CACHE)
    blob["meta"]["tool_versions"].pop("ectyper")
    path = write(tmp_path, "no-ectyper-version.json", blob)
    out = write(tmp_path, "tier2_typer_cache.json", load(TYPER_CACHE))
    with pytest.raises(SystemExit) as exc:
        refresh_tier2.main(["--merge", path, "--out", out])
    assert "meta.tool_versions.ectyper" in str(exc.value)


def test_a_merge_that_would_drop_the_amrfinder_database_version_is_refused(tmp_path):
    """The AMRFinder database date is an epoch too, and it lives only in meta."""
    blob = load(TYPER_CACHE)
    blob["meta"].pop("amrfinder_database_version")
    path = write(tmp_path, "no-db-version.json", blob)
    out = write(tmp_path, "tier2_typer_cache.json", load(TYPER_CACHE))
    with pytest.raises(SystemExit) as exc:
        refresh_tier2.main(["--merge", path, "--out", out])
    assert "meta.amrfinder_database_version" in str(exc.value)


def test_a_merge_staged_to_a_scratch_path_says_it_cannot_check_for_dropped_keys(tmp_path, capsys):
    """CI merges into RUNNER_TEMP and moves the file afterwards, so --out is empty.

    The guard has nothing to compare against there. It has to SAY so: silence is what
    made the lost ECTyper result invisible in the first place.
    """
    frags = nine_shard_fragments(tmp_path)
    assert refresh_tier2.main(["--merge"] + frags
                              + ["--out", str(tmp_path / "staged.json")]) == 0
    err = capsys.readouterr().err
    assert "cannot be checked against the cache it replaces" in err
    assert "tier2_typer_cache.json" in err


def test_a_base_fills_only_the_keys_no_fragment_answers_so_a_refresh_can_move_a_call():
    """A base is a floor, not a peer: the fragment is the fresh run and it wins."""
    merged = refresh_tier2.merge_caches(
        [("fresh", frag({"sistr|GCF_1.1": {"serovar": "Typhi"}},
                        tool_versions={"sistr": "sistr@1.1.3"}))],
        base=[("committed", frag({"sistr|GCF_1.1": {"serovar": "Paratyphi"},
                                  "ectyper|GCF_2.1": {"ectyper_serotype": "O157:H7"}},
                                 tool_versions={"sistr": "sistr@1.1.2",
                                                "ectyper": "ectyper@2.0.0"}))])
    assert merged["results"]["sistr|GCF_1.1"] == {"serovar": "Typhi"}
    assert merged["results"]["ectyper|GCF_2.1"] == {"ectyper_serotype": "O157:H7"}
    assert merged["meta"]["tool_versions"] == {"ectyper": "ectyper@2.0.0",
                                               "sistr": "sistr@1.1.3"}
    assert merged["meta"]["n_results"] == 2


def test_two_fragments_still_conflict_even_when_a_base_answers_the_same_key():
    """A base must not soften the strictness between peers."""
    with pytest.raises(refresh_tier2.MergeConflict) as exc:
        refresh_tier2.merge_caches(
            [("good", frag({"x|1": {"ani": 97.0}})),
             ("stale", frag({"x|1": {"ani": 91.0}}))],
            base=[("committed", frag({"x|1": {"ani": 97.0}}))])
    assert "x|1" in str(exc.value)


def test_base_and_allow_drop_are_rejected_outside_a_merge(tmp_path):
    for argv in (["--print-accessions", "--base", TYPER_CACHE],
                 ["--tools", "sistr", "--out", str(tmp_path / "f.json"), "--allow-drop"]):
        with pytest.raises(SystemExit):
            refresh_tier2.main(argv)


def test_a_deliberate_partial_merge_is_allowed_when_it_says_so(tmp_path):
    blob = load(TYPER_CACHE)
    blob["results"] = {"abricate|GCF_000834295.1": blob["results"]["abricate|GCF_000834295.1"]}
    path = write(tmp_path, "one.json", blob)
    out = str(tmp_path / "out.json")
    assert refresh_tier2.main(["--merge", path, "--out", out, "--allow-partial"]) == 0
    assert load(out)["meta"]["n_results"] == 1


# --------------------------------------------------------------------------- #
#  Shard mode guard rails                                                      #
# --------------------------------------------------------------------------- #
def test_a_shard_refuses_to_run_a_tool_that_is_not_in_a_pinned_environment(tmp_path):
    """Falling back to a bare-name PATH lookup would refresh with an unpinned build."""
    with pytest.raises(SystemExit) as exc:
        refresh_tier2.main(["--tools", "abricate", "--out", str(tmp_path / "f.json"),
                            "--typer-env-dir", str(tmp_path / "no-envs")])
    assert "preflight failed" in str(exc.value)


def test_a_shard_for_a_tool_with_no_planned_job_says_so_rather_than_writing_nothing(tmp_path):
    with pytest.raises(SystemExit) as exc:
        refresh_tier2.main(["--tools", "ectyper", "--out", str(tmp_path / "f.json")])
    assert "no planned job" in str(exc.value)


def test_the_modes_are_mutually_exclusive(tmp_path):
    for argv in (["--print-accessions", "--skani", "--out", str(tmp_path / "f.json")],
                 ["--tools", "sistr", "--skani", "--out", str(tmp_path / "f.json")],
                 ["--tools", "sistr"]):
        with pytest.raises(SystemExit):
            refresh_tier2.main(argv)


def test_an_unknown_tool_name_is_rejected_before_anything_runs(tmp_path):
    with pytest.raises(SystemExit):
        refresh_tier2.main(["--tools", "sistr,not_a_typer", "--out", str(tmp_path / "f.json")])
