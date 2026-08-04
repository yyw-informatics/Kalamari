# -*- coding: utf-8 -*-
"""The four-point gate a pin change has to pass.

The gate's promise is that no pinned thing can move quietly and no Tier-2 call can
change without something explaining it.  These tests attack that promise from both
ends: they tamper with a pinned artifact and check the gate notices, and they move a
call with nothing pinned having moved and check the gate refuses to file it under a
category.  The last test is the one CI depends on -- ``--check`` has to pass on the
committed tree exactly as it stands.

Everything runs offline over committed files.  The tampering happens on COPIES in
``tmp_path``; nothing here writes into the repo.
"""

import json
import os
import shutil

import pytest

import tier2
import update_gate

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")
PINS = os.path.join(SRC, "tool_pins.tsv")
CHANGELOG = os.path.join(SRC, "pin_changelog.tsv")
ENVS = os.path.join(SRC, "envs")


@pytest.fixture(scope="module")
def state():
    return update_gate.pin_state()


@pytest.fixture(scope="module")
def panel():
    return update_gate.replay_panel()


@pytest.fixture(scope="module")
def changelog():
    return update_gate.load_changelog()


def fake_repo(tmp_path, edit_pins=None, edit_lockfile=None):
    """A throwaway repo root holding a copy of the pin table and the lockfiles.

    The digest check joins every lockfile path against a repo root, so tampering has
    to happen inside a copied tree; the committed files are never touched.
    """
    root = tmp_path / "repo"
    (root / "src" / "curation-triage").mkdir(parents=True)
    shutil.copytree(ENVS, str(root / "src" / "curation-triage" / "envs"))
    text = open(PINS, encoding="utf-8").read()
    if edit_pins:
        text = edit_pins(text)
    (root / "src" / "curation-triage" / "tool_pins.tsv").write_text(text, encoding="utf-8")
    shutil.copy(os.path.join(SRC, "marker_manifest.tsv"),
                str(root / "src" / "curation-triage" / "marker_manifest.tsv"))
    if edit_lockfile:
        name, transform = edit_lockfile
        path = root / "src" / "curation-triage" / "envs" / name
        path.write_text(transform(path.read_text(encoding="utf-8")), encoding="utf-8")
    return root


def pins_of(root):
    return update_gate.load_pins(str(root / "src" / "curation-triage" / "tool_pins.tsv"))


def working_src(tmp_path, name="src"):
    """A throwaway copy of the committed inputs, for tests that edit curation data."""
    path = str(tmp_path / name)
    shutil.copytree(SRC, path)
    return path


# --------------------------------------------------------------------------- #
#  Point 1 -- the immutable digest                                             #
# --------------------------------------------------------------------------- #
def test_the_digest_check_catches_a_tampered_lockfile(tmp_path):
    """One appended line is enough: the digest is the artifact's whole identity."""
    root = fake_repo(tmp_path, edit_lockfile=(
        "skani.linux-64.lock",
        lambda text: text + "https://conda.anaconda.org/conda-forge/linux-64/"
                            "sneaky-1.0-0.conda#%s\n" % ("d" * 32)))
    problems = update_gate.digest_problems(pins_of(root), str(root))
    assert len(problems) == 1
    assert "skani.linux-64.lock" in problems[0] and "sha256" in problems[0]
    # The failure has to carry the fix, not just the fact.
    assert "conda list --explicit --md5" in problems[0]
    assert "pin_changelog.tsv" in problems[0]


def test_the_digest_check_passes_on_the_untouched_lockfiles(tmp_path):
    """The same copied tree with nothing edited must produce no finding at all."""
    root = fake_repo(tmp_path)
    assert update_gate.digest_problems(pins_of(root), str(root)) == []


def test_a_build_recorded_only_in_the_pin_table_is_rejected(tmp_path):
    """The build is half the changelog identity, so the lockfile must install it."""
    root = fake_repo(tmp_path, edit_pins=lambda t: t.replace("0.3.2-h79ce301_0",
                                                             "0.3.3-h79ce301_0"))
    problems = update_gate.digest_problems(pins_of(root), str(root))
    assert any("installs no such package" in p for p in problems)


def test_the_dated_amrfinder_database_is_checked_as_a_string_not_as_a_lockfile(tmp_path):
    """No lockfile can capture it, so its shape and its second recording are the check."""
    root = fake_repo(tmp_path, edit_pins=lambda t: t.replace(
        "\tpost-install\t2026-05-15.1\t", "\tpost-install\tlatest\t"))
    problems = update_gate.digest_problems(pins_of(root), str(root))
    assert any("post-install" in p and "database_version" in p for p in problems)
    assert any("marker_manifest.tsv" in p for p in problems)


def test_a_bundled_database_with_no_recorded_identity_is_rejected(tmp_path):
    """A bundled database that names itself nothing pins nothing."""
    root = fake_repo(tmp_path, edit_pins=lambda t: t.replace(
        "\tbundled\tbundled 2020-06-04\t", "\tbundled\t\t"))
    problems = update_gate.digest_problems(pins_of(root), str(root))
    assert any("database_version is empty" in p for p in problems)


def test_no_panel_call_may_depend_on_a_tool_whose_data_is_a_runtime_download(panel):
    """ECTyper's 943 MiB Zenodo sketch is pinned by nothing, so it decides nothing."""
    runtime = {"tool:%s" % r["env"] for r in update_gate.load_pins()
               if r["database_kind"] == "runtime-download"}
    assert runtime == {"tool:ectyper"}, "the runtime-download set moved: %s" % runtime
    assert not [c for c in panel.calls if runtime & set(c.components)]
    assert update_gate.digest_problems(update_gate.load_pins(), REPO_ROOT, panel) == []


def test_a_panel_call_that_depended_on_ectyper_would_fail_the_gate(panel):
    """The guard has to bite, so make a call depend on it and watch point 1 fail."""
    hijacked = update_gate.Panel(
        calls=[update_gate.Call(key="t2:X#marker:y", components=["tool:ectyper"])],
        units=["X"])
    problems = update_gate.digest_problems(update_gate.load_pins(), REPO_ROOT, hijacked)
    assert any("fetched at run time" in p for p in problems)


# --------------------------------------------------------------------------- #
#  Point 2 -- the validation panel                                             #
# --------------------------------------------------------------------------- #
def test_the_panel_is_the_committed_tier2_evidence_replayed_offline(panel):
    """41 trigger units, 75 calls: the same numbers the real Tier-2 run produces."""
    assert len(panel.units) == 41
    assert len(panel.calls) == 75
    assert panel.unevaluated == []
    by_check = {}
    for call in panel.calls:
        by_check[call.check] = by_check.get(call.check, 0) + 1
    assert by_check == {tier2.CHECK_ANI: 37, tier2.CHECK_MARKER: 21, tier2.CHECK_PANEL: 17}


def test_the_call_table_records_no_ani_number(panel):
    """A table that reports a difference on every skani bump teaches everyone to
    ignore it, so the measurement stays out and only the call goes in."""
    for call in panel.calls:
        for field in (call.verdict, call.call, call.flags):
            assert not any(part.replace(".", "", 1).isdigit() and "." in part
                           for part in field.split()), \
                "%s carries a number: %r" % (call.key, field)


def test_a_written_call_table_reads_back_with_the_same_digest(tmp_path, panel):
    """--compare works on files, so the digest has to survive the round trip."""
    path = str(tmp_path / "panel.tsv")
    update_gate.write_panel(panel, path)
    assert update_gate.panel_digest(update_gate.read_panel(path)) == panel.digest


def test_every_panel_call_names_the_components_it_depends_on(panel, state):
    """A call with no components can never be attributed, so it can never be
    explained -- every one of them must reduce to pins that exist."""
    for call in panel.calls:
        assert call.components, "%s depends on nothing" % call.key
        for component in call.components:
            assert component in state, "%s depends on unpinned %s" % (call.key, component)
            assert update_gate.component_class(component) in update_gate.CHANGE_CLASSES


def test_a_marker_call_depends_on_its_task_tools_and_an_ani_call_on_skani(panel):
    by_key = {c.key: c for c in panel.calls}
    salmonella = by_key["t2:Salmonella_enterica_I~GCF_000006945.2"
                        "#marker:salmonella_subspecies_and_serovar"]
    assert set(salmonella.components) == {
        "tool:sistr", "db:sistr", "tool:seqsero2", "db:seqsero2",
        "rule:salmonella_subspecies_and_serovar",
        "threshold:salmonella_subspecies_and_serovar"}
    ani = next(c for c in panel.calls if c.check == tier2.CHECK_ANI)
    assert set(ani.components) == {"tool:skani", update_gate.ANI_THRESHOLD_COMPONENT}


def test_a_panel_that_lost_a_unit_fails_rather_than_validating_less(changelog):
    """A shrinking panel is the quiet failure: it still passes, over less evidence."""
    shrunk = update_gate.Panel(calls=[], units=["a", "b"], unevaluated=["a", "b"])
    problems = update_gate.panel_problems(shrunk, changelog)
    assert any("replayed 0 calls" in p for p in problems)
    assert any("no Tier-2 evidence" in p for p in problems)


def test_the_panel_scope_is_frozen_in_the_gate_and_covers_both_kinds_of_unit(panel):
    """The 41 units are named here, not read from a file curation rewrites."""
    assert len(update_gate.PANEL_SCOPE) == 41
    assert set(panel.units) == set(update_gate.PANEL_SCOPE)
    assert panel.missing == []
    # Both halves are represented: the mandatory confounded regimes and the Tier-1
    # leads. If one half ever emptied, the panel would stop validating a whole path.
    assert len(update_gate.PANEL_CONFOUNDED_UNITS) == 21
    assert len(update_gate.PANEL_LEAD_UNITS) == 20


def test_resolving_a_tier1_lead_does_not_move_the_panel(tmp_path, panel, capsys):
    """Ordinary curation must not fail the gate.

    The panel used to be built from the ledger, so the monthly collate appending a
    `resolved` record dropped that unit, shrank the panel from 75 calls to 74 and
    moved the digest. Point 2 then failed with no row anyone could legally add: no
    pin had moved, and a changelog row may not record from == to. A resolved lead is
    a correct outcome, so it has to leave the panel exactly where it was.
    """
    src = working_src(tmp_path)
    unit = update_gate.PANEL_LEAD_UNITS[0]
    ledger_path = os.path.join(src, "ledger.ndjson")
    with open(ledger_path, encoding="utf-8") as fh:
        records = [json.loads(line) for line in fh if line.strip()]
    standing = [r for r in records
                if r.get("unit_id") == unit and r.get("signal") != "T2"
                and r.get("actionable_flags")]
    assert standing, "%s is no longer a standing Tier-1 lead" % unit
    resolved = dict(standing[-1], delta="resolved", run_id="2026-09",
                    summary="an SME reviewed this lead and closed it")
    with open(ledger_path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(resolved, sort_keys=True) + "\n")

    replayed = update_gate.replay_panel(src=src)
    assert len(replayed.units) == len(panel.units)
    assert len(replayed.calls) == len(panel.calls)
    assert replayed.digest == panel.digest
    # ...and the gate CI runs still passes over that tree.
    assert update_gate.main(["--check", "--src", src]) == 0
    assert "all four points pass." in capsys.readouterr().out


def test_a_scope_unit_the_policy_table_no_longer_has_is_reported_not_dropped(tmp_path):
    """A silently smaller panel is the failure mode a frozen scope must not create."""
    panel = update_gate.replay_panel(
        scope=list(update_gate.PANEL_SCOPE) + ["Nothing_atall~GCF_000000000.1"])
    assert panel.missing == ["Nothing_atall~GCF_000000000.1"]
    problems = update_gate.panel_problems(panel, [])
    assert any("frozen panel scope" in p and "PANEL_SCOPE" in p for p in problems)


# -- point 2, second half: the caches have to come from the pinned tools ------ #
def test_the_committed_caches_say_they_were_produced_by_the_pinned_tools(panel):
    """The digest only certifies a pin state if that pin state produced the evidence."""
    assert update_gate.provenance_problems(panel, update_gate.load_pins()) == []
    assert panel.provenance["skani"] == "skani@0.3.2"
    assert panel.provenance[update_gate.AMRFINDER_DB_COMPONENT] == "2026-05-15.1"


def test_a_panel_replayed_from_the_previous_pins_caches_fails_point_two(panel):
    """The hole this closes: bump the lockfile, the pin table and the manifest, add a
    changelog row, and every one of the 75 validated calls still came from the OLD
    tool -- while the changelog permanently records the new pin as validated."""
    stale = update_gate.Panel(calls=panel.calls, units=panel.units,
                              provenance=dict(panel.provenance, skani="skani@0.3.1"))
    problems = update_gate.provenance_problems(stale, update_gate.load_pins())
    assert len(problems) == 1
    assert "skani@0.3.1" in problems[0] and "skani@0.3.2" in problems[0]


def test_a_cache_built_with_another_amrfinder_database_fails_point_two(panel):
    """The one database that moves without its package moving, so the one that can
    slip through every lockfile check."""
    stale = update_gate.Panel(
        calls=panel.calls, units=panel.units,
        provenance=dict(panel.provenance,
                        **{update_gate.AMRFINDER_DB_COMPONENT: "2025-01-01.1"}))
    problems = update_gate.provenance_problems(stale, update_gate.load_pins())
    assert len(problems) == 1 and "2025-01-01.1" in problems[0]
    assert "2026-05-15.1" in problems[0]


def test_a_tool_that_cannot_report_its_version_is_accepted_only_as_unknown(panel):
    """`shigeifinder --version` exits 2, so `shigeifinder@unknown` IS its identity --
    and a version string turning up there means the cache came from somewhere else."""
    assert panel.provenance["shigeifinder"] == "shigeifinder@unknown"
    assert update_gate.provenance_problems(panel, update_gate.load_pins()) == []
    claimed = update_gate.Panel(
        calls=panel.calls, units=panel.units,
        provenance=dict(panel.provenance, shigeifinder="shigeifinder@1.3.4"))
    problems = update_gate.provenance_problems(claimed, update_gate.load_pins())
    assert len(problems) == 1 and "cannot report its own version" in problems[0]


def test_caches_that_record_no_versions_at_all_certify_nothing(panel):
    blank = update_gate.Panel(calls=panel.calls, units=panel.units, provenance={})
    problems = update_gate.provenance_problems(blank, update_gate.load_pins())
    assert len(problems) == 1 and "record no tool versions" in problems[0]


def test_a_tool_a_call_depends_on_must_appear_in_the_caches_provenance(panel):
    """Deleting the entry must not be a way to pass the check the entry exists for."""
    thin = update_gate.Panel(
        calls=panel.calls, units=panel.units,
        provenance={k: v for k, v in panel.provenance.items() if k != "sistr"})
    problems = update_gate.provenance_problems(thin, update_gate.load_pins())
    assert any("tool:sistr" in p and "no version" in p for p in problems)


def test_a_moved_panel_fails_against_the_accepted_digest(panel, changelog):
    """The accepted call table lives in the changelog's validation column."""
    moved = update_gate.Panel(
        calls=[update_gate.Call(unit_id="u", check="ani", key="t2:u#ani",
                                verdict="confirmed", components=["tool:skani"])],
        units=["u"])
    problems = update_gate.panel_problems(moved, changelog)
    assert any("the validation panel moved" in p for p in problems)
    assert any("panel_calls_total" in p for p in problems)


# --------------------------------------------------------------------------- #
#  Point 3 -- classification                                                   #
# --------------------------------------------------------------------------- #
def call_on(*components):
    return update_gate.Call(unit_id="u", check="marker", key="t2:u#marker:t",
                            verdict="confirmed", components=list(components))


BASE_STATE = {
    "tool:sistr": "1.1.3+lock:aaaaaaaaaaaa",
    "db:sistr": "SISTR_V_1.1.3_db",
    "rule:salmonella_subspecies_and_serovar": "salmonella_subspecies_serovar@1",
    "threshold:salmonella_subspecies_and_serovar": "serovar_cgmlst[min_identity=90]",
}
DEPENDS = tuple(BASE_STATE)


def moved_state(component, value):
    after = dict(BASE_STATE)
    after[component] = value
    return after


@pytest.mark.parametrize("component,value,expected", [
    ("tool:sistr", "1.1.4+lock:bbbbbbbbbbbb", update_gate.CLASS_SOFTWARE),
    ("db:sistr", "SISTR_V_1.1.4_db", update_gate.CLASS_DATABASE),
    ("threshold:salmonella_subspecies_and_serovar", "serovar_cgmlst[min_identity=95]",
     update_gate.CLASS_THRESHOLD),
    ("rule:salmonella_subspecies_and_serovar", "salmonella_subspecies_serovar@2",
     update_gate.CLASS_RULE),
])
def test_classification_names_the_pinned_component_that_moved(component, value, expected):
    """One component moves at a time: each of the four classes, from a real shape."""
    before, after = call_on(*DEPENDS), call_on(*DEPENDS)
    after.verdict = "refuted"
    assert update_gate.classify_difference(
        before, after, BASE_STATE, moved_state(component, value)) == expected


def test_a_call_that_changed_with_nothing_moving_is_unattributable_and_fails(tmp_path):
    """The whole point: same pins, different answer, so the pipeline is not
    reproducible and there is no category that makes that acceptable."""
    before, after = call_on(*DEPENDS), call_on(*DEPENDS)
    after.verdict = "refuted"
    assert update_gate.classify_difference(before, after, BASE_STATE, dict(BASE_STATE)) \
        == update_gate.CLASS_UNATTRIBUTABLE

    before_path, after_path = str(tmp_path / "b.tsv"), str(tmp_path / "a.tsv")
    update_gate.write_panel(update_gate.Panel(calls=[before], units=["u"]), before_path)
    update_gate.write_panel(update_gate.Panel(calls=[after], units=["u"]), after_path)
    assert update_gate.main(["--compare", before_path, after_path,
                             "--pins-before", PINS]) == 1


def test_a_change_a_moved_pin_explains_is_reported_and_does_not_fail(tmp_path, panel):
    """A real, attributable pin change: the gate reports it and exits 0."""
    before_path, after_path = str(tmp_path / "b.tsv"), str(tmp_path / "a.tsv")
    update_gate.write_panel(panel, before_path)
    moved = update_gate.Panel(
        calls=[update_gate.Call(**dict(vars(c), verdict="refuted"))
               if c.key.endswith("#marker:bacillus_anthracis_vs_b_cereus_group")
               else c for c in panel.calls],
        units=panel.units)
    update_gate.write_panel(moved, after_path)
    older = str(tmp_path / "pins_before.tsv")
    with open(older, "w", encoding="utf-8") as fh:
        fh.write(open(PINS, encoding="utf-8").read().replace("3.4.0-pyhdfd78af_0",
                                                             "3.3.0-pyhdfd78af_0"))
    assert update_gate.main(["--compare", before_path, after_path,
                             "--pins-before", older]) == 0


def test_a_call_that_appeared_or_disappeared_is_classified_the_same_way():
    """A dropped call is a changed call; None on either side must not crash."""
    after = call_on(*DEPENDS)
    software = moved_state("tool:sistr", "1.1.4+lock:bbbbbbbbbbbb")
    assert update_gate.classify_difference(None, after, BASE_STATE, software) \
        == update_gate.CLASS_SOFTWARE
    assert update_gate.classify_difference(after, None, BASE_STATE, software) \
        == update_gate.CLASS_SOFTWARE
    assert update_gate.classify_difference(after, None, BASE_STATE, dict(BASE_STATE)) \
        == update_gate.CLASS_UNATTRIBUTABLE


def test_a_component_the_call_does_not_depend_on_cannot_explain_it():
    """Attribution is per call: a bumped lissero says nothing about a sistr call."""
    before, after = call_on(*DEPENDS), call_on(*DEPENDS)
    after.verdict = "refuted"
    unrelated = dict(BASE_STATE)
    unrelated["tool:lissero"] = "0.4.11+lock:cccccccccccc"
    assert update_gate.classify_difference(before, after, BASE_STATE, unrelated) \
        == update_gate.CLASS_UNATTRIBUTABLE


def test_our_own_edits_outrank_a_vendor_bump_when_both_moved():
    """Reported class is the root cause a reviewer can revert fastest; the row still
    lists every component that moved, so the ambiguity stays visible."""
    before, after = call_on(*DEPENDS), call_on(*DEPENDS)
    after.verdict = "refuted"
    both = moved_state("tool:sistr", "1.1.4+lock:bbbbbbbbbbbb")
    both["rule:salmonella_subspecies_and_serovar"] = "salmonella_subspecies_serovar@2"
    assert update_gate.classify_difference(before, after, BASE_STATE, both) \
        == update_gate.CLASS_RULE
    assert update_gate.moved_components(DEPENDS, BASE_STATE, both) == [
        "rule:salmonella_subspecies_and_serovar", "tool:sistr"]


def test_a_bundled_database_moving_with_its_package_is_reported_as_software():
    """A bundled database moves only because its package did, so the package is the
    root cause; only a database that moved on its own is a `database` change."""
    before, after = call_on(*DEPENDS), call_on(*DEPENDS)
    after.verdict = "refuted"
    both = moved_state("tool:sistr", "1.1.4+lock:bbbbbbbbbbbb")
    both["db:sistr"] = "SISTR_V_1.1.4_db"
    assert update_gate.classify_difference(before, after, BASE_STATE, both) \
        == update_gate.CLASS_SOFTWARE


def test_a_call_the_scope_added_or_removed_is_a_panel_change_not_a_failure():
    """A unit entering or leaving the frozen scope explains a call that appeared or
    disappeared. Calling that `unattributable` was a false alarm -- it told a curator
    the pipeline was not reproducible when nothing about it had changed."""
    call = call_on(*DEPENDS)
    before = dict(BASE_STATE, **{update_gate.PANEL_SCOPE_COMPONENT: "40 units+scope:aaa"})
    after = dict(BASE_STATE, **{update_gate.PANEL_SCOPE_COMPONENT: "41 units+scope:bbb"})
    assert update_gate.classify_difference(None, call, before, after) \
        == update_gate.CLASS_PANEL
    assert update_gate.classify_difference(call, None, before, after) \
        == update_gate.CLASS_PANEL
    # The scope is on the difference row, so `moved` never reads as an empty `-`.
    differences = update_gate.diff_panels([], [call], before, after)
    assert [d.moved for d in differences] == [[update_gate.PANEL_SCOPE_COMPONENT]]


def test_a_verdict_that_moved_during_a_scope_change_is_still_unattributable():
    """The scope decides WHICH units are replayed, never what a call says, so it must
    not become a blanket excuse for a call that changed on both sides."""
    before_call, after_call = call_on(*DEPENDS), call_on(*DEPENDS)
    after_call.verdict = "refuted"
    before = dict(BASE_STATE, **{update_gate.PANEL_SCOPE_COMPONENT: "40 units+scope:aaa"})
    after = dict(BASE_STATE, **{update_gate.PANEL_SCOPE_COMPONENT: "41 units+scope:bbb"})
    assert update_gate.classify_difference(before_call, after_call, before, after) \
        == update_gate.CLASS_UNATTRIBUTABLE


def test_a_real_scope_change_classifies_every_call_it_moved_as_panel(tmp_path, state,
                                                                     panel):
    """End to end over the committed evidence: drop one unit from the scope, replay,
    and every difference has to be attributable to the scope alone."""
    dropped = update_gate.PANEL_LEAD_UNITS[0]
    smaller = [u for u in update_gate.PANEL_SCOPE if u != dropped]
    before = update_gate.replay_panel(scope=smaller)
    assert len(before.calls) < len(panel.calls)
    before_state = dict(state)
    before_state[update_gate.PANEL_SCOPE_COMPONENT] = \
        update_gate.panel_scope_identity(smaller)
    differences = update_gate.diff_panels(before.calls, panel.calls, before_state, state)
    assert differences
    assert {d.classification for d in differences} == {update_gate.CLASS_PANEL}
    assert all(d.key.startswith("t2:%s#" % dropped) for d in differences)


def test_diff_panels_ignores_calls_that_did_not_move(panel, state):
    """75 identical calls must produce zero differences, or every run cries wolf."""
    assert update_gate.diff_panels(panel.calls, panel.calls, state, state) == []


# --------------------------------------------------------------------------- #
#  Point 4 -- the changelog                                                    #
# --------------------------------------------------------------------------- #
def write_changelog(tmp_path, rows, name="pin_changelog.tsv"):
    path = tmp_path / name
    lines = ["# a test changelog", "\t".join(update_gate.CHANGELOG_COLUMNS)]
    lines += ["\t".join(row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def row(component, to_value, from_value="-", class_name=None, date="2026-08-04",
        total="75", changed="0", validation="panel-replay:4c7b9e0f8066",
        reviewer="initial-pin", note="a note"):
    return (date, component, class_name or update_gate.component_class(component),
            from_value, to_value, total, changed, validation, reviewer, note)


def test_a_pin_change_with_no_changelog_row_fails_the_gate(tmp_path, state, panel):
    """The gate in one comparison: the tree says one thing, the record another."""
    recorded = [row(component, identity) for component, identity in state.items()]
    path = write_changelog(tmp_path, recorded)
    moved = dict(state)
    moved["tool:sistr"] = "1.1.4-pyhdc42f0e_0+lock:6a1a53a28066"
    problems = update_gate.pin_change_problems(moved, update_gate.load_changelog(path),
                                               panel, path)
    assert len(problems) == 1
    assert "tool:sistr changed" in problems[0]
    # The message has to be a recipe, not a complaint.
    assert "component=tool:sistr" in problems[0] and "class=software" in problems[0]
    assert "panel_calls_total=75" in problems[0]
    assert "--compare" in problems[0]


def test_a_newly_pinned_component_with_no_row_at_all_fails(tmp_path, state, panel):
    """A pin nobody recorded is as unexplainable as a pin that moved unrecorded."""
    recorded = [row(c, i) for c, i in state.items() if c != "tool:skani"]
    path = write_changelog(tmp_path, recorded)
    problems = update_gate.pin_change_problems(state, update_gate.load_changelog(path),
                                               panel, path)
    assert len(problems) == 1 and "has no row for it" in problems[0]


def test_a_component_that_left_the_tree_needs_a_removal_row(tmp_path, state, panel):
    recorded = [row(c, i) for c, i in state.items()]
    recorded.append(row("tool:gone", "1.0+lock:aaaaaaaaaaaa"))
    path = write_changelog(tmp_path, recorded)
    problems = update_gate.pin_change_problems(state, update_gate.load_changelog(path),
                                               panel, path)
    assert len(problems) == 1 and "nothing in the tree pins it any more" in problems[0]

    recorded.append(row("tool:gone", update_gate.NO_PIN,
                        from_value="1.0+lock:aaaaaaaaaaaa", date="2026-09-01"))
    path = write_changelog(tmp_path, recorded, name="removed.tsv")
    assert update_gate.pin_change_problems(state, update_gate.load_changelog(path),
                                           panel, path) == []


def test_unattributable_is_rejected_as_a_recorded_class(tmp_path):
    """It is a failure, not a fifth category, so it cannot be filed away as one."""
    path = write_changelog(tmp_path, [row("tool:sistr", "1.1.3+lock:aaaaaaaaaaaa",
                                          class_name="unattributable")])
    problems = update_gate.changelog_problems(update_gate.load_changelog(path), path)
    assert any("FAILURE, not a category" in p for p in problems)


def test_a_row_whose_class_contradicts_its_component_is_rejected(tmp_path):
    path = write_changelog(tmp_path, [row("db:sistr", "SISTR_V_1.1.3_db",
                                          class_name="software")])
    problems = update_gate.changelog_problems(update_gate.load_changelog(path), path)
    assert any("is a database pin but the row says software" in p for p in problems)


def test_a_broken_chain_of_identities_is_rejected(tmp_path):
    """A pin state nobody recorded means a call from that state can never be explained."""
    path = write_changelog(tmp_path, [
        row("tool:sistr", "1.1.3+lock:aaaaaaaaaaaa"),
        row("tool:sistr", "1.1.5+lock:cccccccccccc", from_value="1.1.4+lock:bbbbbbbbbbbb",
            date="2026-09-01", changed="2"),
    ])
    problems = update_gate.chain_problems(update_gate.load_changelog(path), path)
    assert len(problems) == 1 and "the chain cannot skip" in problems[0]


def test_the_first_row_of_a_component_must_start_from_no_pin(tmp_path):
    path = write_changelog(tmp_path, [row("tool:sistr", "1.1.3+lock:aaaaaaaaaaaa",
                                          from_value="1.1.2+lock:zzzzzzzzzzzz")])
    problems = update_gate.chain_problems(update_gate.load_changelog(path), path)
    assert any("must start from -" in p for p in problems)


def test_a_row_without_a_reviewer_or_a_panel_digest_is_rejected(tmp_path):
    path = write_changelog(tmp_path, [
        row("tool:sistr", "1.1.3+lock:aaaaaaaaaaaa", reviewer=""),
        row("tool:skani", "0.3.2+lock:bbbbbbbbbbbb", validation="looked fine to me"),
    ])
    problems = update_gate.changelog_problems(update_gate.load_changelog(path), path)
    assert any("reviewer is empty" in p for p in problems)
    assert any("panel-replay:" in p and "validation is" in p for p in problems)


def test_a_row_that_moved_calls_without_saying_why_is_rejected(tmp_path):
    """panel_calls_changed=3 with an empty note is a number nobody can act on."""
    path = write_changelog(tmp_path, [row("tool:sistr", "1.1.3+lock:aaaaaaaaaaaa",
                                          changed="3", note="")])
    assert any("the note is empty" in p
               for p in update_gate.classification_problems(
                   update_gate.load_changelog(path)))


def test_the_changelog_must_be_append_ordered_in_time(tmp_path):
    path = write_changelog(tmp_path, [
        row("tool:sistr", "1.1.3+lock:aaaaaaaaaaaa", date="2026-09-01"),
        row("tool:skani", "0.3.2+lock:bbbbbbbbbbbb", date="2026-08-04"),
    ])
    problems = update_gate.changelog_problems(update_gate.load_changelog(path), path)
    assert any("append-ordered" in p for p in problems)


def test_a_scope_change_is_recordable_as_a_panel_row(tmp_path, state, panel):
    """The gate must never be a dead end. A scope move is not a pin move, so it needs
    a class of its own -- and that row has to pass every point-4 check."""
    smaller = update_gate.panel_scope_identity(
        [u for u in update_gate.PANEL_SCOPE if u != update_gate.PANEL_LEAD_UNITS[0]])
    scope = update_gate.PANEL_SCOPE_COMPONENT
    recorded = [row(c, i) for c, i in state.items() if c != scope]
    recorded.append(row(scope, smaller))
    recorded.append(row(scope, state[scope], from_value=smaller, date="2026-09-01",
                        changed="1", note="put the unit back in the panel"))
    path = write_changelog(tmp_path, recorded)
    rows = update_gate.load_changelog(path)
    assert update_gate.changelog_problems(rows, path) == []
    assert update_gate.chain_problems(rows, path) == []
    assert update_gate.classification_problems(rows) == []
    assert update_gate.pin_change_problems(state, rows, panel, path) == []
    assert update_gate.latest_by_component(rows)[scope].class_name \
        == update_gate.CLASS_PANEL


def test_previous_state_is_the_last_recorded_identity_not_the_one_before_it(tmp_path,
                                                                            state):
    """--compare runs at step 3 of the documented loop: the pin has already moved and
    its new row does not exist yet. So the before-state is each component's last
    recorded `to`. Rolling back to `from` instead is one change too far -- against the
    baseline rows (`from` = `-`) it returns the CURRENT state, so `--compare` reports
    that nothing moved and files every real pin bump as unattributable."""
    recorded = [row(c, i) for c, i in state.items()]
    path = write_changelog(tmp_path, recorded)
    moved = dict(state)                       # the tree, after the bump
    moved["tool:sistr"] = "1.1.4+lock:dddddddddddd"
    earlier = update_gate.previous_state(moved, update_gate.load_changelog(path))
    assert earlier["tool:sistr"] == state["tool:sistr"]
    assert earlier["tool:skani"] == state["tool:skani"]
    assert update_gate.moved_components(["tool:sistr", "tool:skani"], earlier, moved) \
        == ["tool:sistr"]


def test_previous_state_rolls_back_to_the_last_change_of_a_component_that_moved_before(
        tmp_path, state):
    """The same rule when the component already has a change row: the before-state is
    where that row ENDED, not where it started."""
    recorded = [row(c, i) for c, i in state.items()]
    recorded.append(row("tool:sistr", "1.1.4+lock:dddddddddddd",
                        from_value=state["tool:sistr"], date="2026-09-01",
                        changed="1", note="bump"))
    path = write_changelog(tmp_path, recorded)
    moved = dict(state)
    moved["tool:sistr"] = "1.1.5+lock:eeeeeeeeeeee"
    earlier = update_gate.previous_state(moved, update_gate.load_changelog(path))
    assert earlier["tool:sistr"] == "1.1.4+lock:dddddddddddd"


def test_previous_state_drops_a_component_whose_last_row_removed_it(tmp_path, state):
    """A removal row says there was no such pin, and a component with no identity
    would read as one that existed and was called `-`."""
    recorded = [row(c, i) for c, i in state.items()]
    recorded.append(row("tool:gone", "1.0+lock:aaaaaaaaaaaa", date="2026-09-01"))
    recorded.append(row("tool:gone", update_gate.NO_PIN, date="2026-09-02",
                        from_value="1.0+lock:aaaaaaaaaaaa", changed="1", note="dropped"))
    path = write_changelog(tmp_path, recorded)
    earlier = update_gate.previous_state(dict(state), update_gate.load_changelog(path))
    assert "tool:gone" not in earlier


# --------------------------------------------------------------------------- #
#  The committed tree                                                          #
# --------------------------------------------------------------------------- #
def test_check_passes_on_the_committed_tree(capsys):
    """This is what CI runs. It has to be green on the tree as committed."""
    assert update_gate.main(["--check"]) == 0
    printed = capsys.readouterr().out
    assert "all four points pass." in printed
    for point in ("1. immutable digest", "2. validation panel",
                  "3. classification", "4. changelog"):
        assert "%-22s PASS" % point in printed


def test_a_missing_changelog_fails_with_a_recipe_rather_than_a_traceback():
    """The first thing a new checkout can get wrong, so it gets a sentence."""
    with pytest.raises(SystemExit) as excinfo:
        update_gate.main(["--check", "--changelog", os.path.join("no", "such.tsv")])
    assert "--print-state" in str(excinfo.value)


def test_the_panel_and_print_state_modes_write_what_a_reviewer_needs(tmp_path, capsys,
                                                                     state):
    """--panel prints the digest a changelog row needs; --print-state prints the
    identities those rows record."""
    out = str(tmp_path / "panel.tsv")
    assert update_gate.main(["--panel", "--out", out]) == 0
    printed = capsys.readouterr().out.strip().split("\t")
    assert printed[0].startswith(update_gate.VALIDATION_PREFIX) and printed[1] == "75"
    assert update_gate.panel_digest(update_gate.read_panel(out)) \
        == printed[0][len(update_gate.VALIDATION_PREFIX):]

    assert update_gate.main(["--print-state"]) == 0
    rows = [line.split("\t") for line in capsys.readouterr().out.splitlines()]
    assert {r[0]: r[2] for r in rows} == state


def test_the_committed_changelog_records_every_pinned_component(state, changelog):
    """38 components, 38 baseline rows: a pin with no row is a pin nobody can explain.

    38 rather than 37 since the panel's scope became a pinned component of its own."""
    assert len(state) == 38
    recorded = update_gate.latest_by_component(changelog)
    assert set(recorded) == set(state)
    for component, identity in state.items():
        assert recorded[component].to_value == identity


def test_the_committed_baseline_records_the_sistr_environment_finding(changelog):
    """A1 found the installed prefix disagreeing with its own lockfile; the outcome
    was that the environment was corrected and no pin moved, which is exactly the
    kind of decision a changelog exists to keep."""
    row_for_sistr = update_gate.latest_by_component(changelog)["tool:sistr"]
    assert row_for_sistr.from_value == update_gate.NO_PIN
    assert row_for_sistr.panel_calls_changed == "0"
    assert "setuptools" in row_for_sistr.note and "lockfile" in row_for_sistr.note


def test_the_amrfinder_database_identity_names_its_bytes_not_only_its_label(changelog, state):
    """The one database no lockfile can capture must be identified by content, not by name.

    `amrfinder --database_version` prints a dated LABEL, and any directory can be given that
    name. Since the freeze the gate's identity carries the roll-up digest of the 162 files as
    well, so a database swapped underneath an unchanged version string fails the gate instead
    of only failing setup_envs.sh.
    """
    row_for_db = update_gate.latest_by_component(changelog)["db:amrfinderplus"]
    assert state["db:amrfinderplus"].startswith("2026-05-15.1+db:")
    assert row_for_db.to_value == state["db:amrfinderplus"]
    assert row_for_db.class_name == update_gate.CLASS_DATABASE


def test_a_database_with_no_recorded_digest_keeps_its_bare_version():
    # The bundled databases need no digest of their own: the lockfile already pins their
    # bytes, and point 1 checks the lockfile. Only the post-install one carries a digest.
    bundled = {"database_version": "vfdb 4592 seqs", "database_digest": "in-lockfile"}
    frozen = {"database_version": "2026-05-15.1",
              "database_digest": "sha256:" + "a" * 64}
    assert update_gate.database_identity(bundled) == "vfdb 4592 seqs"
    assert update_gate.database_identity(frozen) == "2026-05-15.1+db:" + "a" * 12


def test_the_committed_changelog_is_well_formed_and_chained(changelog):
    assert update_gate.changelog_problems(changelog) == []
    assert update_gate.chain_problems(changelog) == []
    assert update_gate.classification_problems(changelog) == []


def test_the_gate_names_no_local_filesystem_path(state, panel, tmp_path):
    """Every message a reviewer reads has to work on any machine."""
    written = str(tmp_path / "panel.tsv")
    update_gate.write_panel(panel, written)
    texts = [open(os.path.join(REPO_ROOT, "bin", "curation-triage",
                               "update_gate.py"), encoding="utf-8").read(),
             open(CHANGELOG, encoding="utf-8").read(),
             open(written, encoding="utf-8").read()]
    for text in texts:
        for leak in ("/scicomp/", "/home/", "/Users/"):
            assert leak not in text
