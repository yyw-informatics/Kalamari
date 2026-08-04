# -*- coding: utf-8 -*-
"""worklist_from_ledger.py -- rendering the release worklist from the committed ledger.

The release asset must show what a curator actually saw. These tests pin the two ways
that can go wrong: recomputing the delta (which would relabel every lead) and losing the
per-row run id (which would misreport when the evidence was last computed).

One rule for everything below that touches ``src/curation-triage/ledger.ndjson``: assert
only what SURVIVES AN APPEND. That ledger is append-only and the monthly bot commits a
grown copy of it every run, while this test file is a BLOCKING check that fires on
``src/curation-triage/**``. Anything pinned to the seed history -- a record count, the
exact set of run ids or signals present -- would therefore go red on the bot's own first
push and block every later merge, turning an informational tool into a gate.
"""

import csv
import json
import os
import re

import pytest

import ledger as ledger_mod
import worklist_from_ledger as wfl
from leads import SIGNAL_J1A, SIGNAL_J1B, SIGNAL_J2, SIGNAL_T2

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
COMMITTED_LEDGER = os.path.join(REPO_ROOT, "src", "curation-triage", "ledger.ndjson")
COMMITTED_REVIEWS = os.path.join(REPO_ROOT, "src", "curation-triage", "reviews.tsv")

KNOWN_SIGNALS = {SIGNAL_J1A, SIGNAL_J1B, SIGNAL_J2, SIGNAL_T2}
KNOWN_DELTAS = {"new", "recurring", "resolved"}
RUN_ID_RE = re.compile(r"^\d{4}-\d{2}$")      # run ids are YYYY-MM, by convention


def _rec(lead_key, run_id, delta, signal="J1a", name="Escherichia_coli",
         flags=("tax_check_failed",), epoch="ncbi-assembly-reports@2026-07"):
    return {
        "schema": 1, "run_id": run_id, "signal": signal, "lead_key": lead_key,
        "epoch": epoch, "delta": delta, "unit_id": lead_key,
        "scientificName": name, "source_lane": "NCBI-REF",
        "candidate_accession": "", "deposit_date": "2020-01-01",
        "actionable_flags": list(flags), "flags": list(flags), "evidence": {},
        "evidence_fingerprint": "deadbeefdeadbeef", "link": "", "summary": "why it fired",
        "params": {}, "tool_versions": {},
    }


def _write(path, records):
    with open(path, "w", encoding="utf-8") as fh:
        for r in records:
            fh.write(json.dumps(r, sort_keys=True) + "\n")
    return str(path)


def test_the_species_name_survives_the_record_to_lead_round_trip():
    # Lead.from_dict drops keys it does not know, and the record spells the field
    # 'scientificName' while the dataclass spells it 'scientific_name'. Without the
    # mapping every rendered row would have a blank species.
    lead = wfl.lead_from_record(_rec("u1", "2026-07", "new", name="Listeria_monocytogenes"))
    assert lead.scientific_name == "Listeria_monocytogenes"
    assert lead.lead_key == "u1"
    assert lead.actionable_flags == ["tax_check_failed"]


def test_the_bucket_comes_from_the_recorded_delta_and_is_never_recomputed():
    # This is the whole point: re-running the pipeline at tag time would diff against a
    # ledger that already holds these records and call every one of them 'recurring'.
    records = [_rec("a", "2026-08", "new"), _rec("b", "2026-08", "recurring")]
    part = wfl.partition_from_ledger(records, [], run_id="2026-08")
    assert [l.lead_key for l in part.new] == ["a"]
    assert [l.lead_key for l in part.standing] == ["b"]


def test_a_resolved_record_is_reported_as_resolved_and_never_as_a_worklist_row():
    records = [_rec("a", "2026-08", "new"), _rec("gone", "2026-08", "resolved")]
    part = wfl.partition_from_ledger(records, [], run_id="2026-08")
    assert len(part.resolved) == 1
    assert part.resolved[0]["lead_key"] == "gone"
    assert "gone" not in [r["lead_key"] for r in ledger_mod.worklist_rows(part)]


def test_a_dismissal_recorded_after_the_run_still_suppresses_the_lead():
    # reviews.tsv is applied at render time, so an SME decision taken after the monthly
    # run does not reappear on the release asset.
    records = [_rec("a", "2026-08", "new")]
    reviews = [{"lead_key": "a", "decision": "reject", "evidence_fingerprint": "",
                "epoch": "", "reviewer": "sme", "date": "2026-08-02", "note": ""}]
    part = wfl.partition_from_ledger(records, reviews, run_id="2026-08")
    assert part.new == []
    assert [l.lead_key for l in part.dismissed] == ["a"]


def test_an_accept_stays_on_the_worklist_as_a_proposed_edit():
    records = [_rec("a", "2026-08", "recurring")]
    reviews = [{"lead_key": "a", "decision": "accept", "evidence_fingerprint": "",
                "epoch": "", "reviewer": "sme", "date": "2026-08-02", "note": ""}]
    part = wfl.partition_from_ledger(records, reviews, run_id="2026-08")
    assert [l.lead_key for l in part.accepted] == ["a"]
    assert part.standing == []


def test_the_snapshot_keeps_only_the_latest_record_per_lead():
    # An append-only log records the same lead every run; the snapshot is the current
    # state, not the whole history.
    records = [_rec("a", "2026-07", "new"), _rec("a", "2026-08", "recurring"),
               _rec("b", "2026-08", "new")]
    part = wfl.partition_from_ledger(records, [], run_id=None)
    assert [l.lead_key for l in part.new] == ["b"]
    assert [l.lead_key for l in part.standing] == ["a"]


def test_the_snapshot_spans_signals_that_carry_different_run_ids():
    # Tier 1 and Tier 2 are appended by different jobs, so "the last run" would show
    # only one tier. The snapshot must show both.
    records = [_rec("j1a-lead", "2026-07", "new"),
               _rec("t2:x#ani", "2026-08", "new", signal="T2", epoch="skani@0.3.2")]
    part = wfl.partition_from_ledger(records, [], run_id=None)
    assert {l.signal for l in part.new} == {"J1a", "T2"}


def test_each_row_carries_its_own_records_run_id_not_the_newest_one(tmp_path):
    records = [_rec("old", "2026-07", "new"), _rec("fresh", "2026-08", "new")]
    part = wfl.partition_from_ledger(records, [], run_id=None)
    out = tmp_path / "worklist.tsv"
    n = wfl.write_worklist(part, str(out), wfl._run_id_of(records))
    assert n == 2
    with open(out, encoding="utf-8") as fh:
        rows = {r["lead_key"]: r["run_id"] for r in csv.DictReader(fh, delimiter="\t")}
    assert rows == {"old": "2026-07", "fresh": "2026-08"}


def test_the_snapshot_heading_never_claims_to_be_a_tier_1_run(tmp_path):
    records = [_rec("a", "2026-08", "new")]
    part = wfl.partition_from_ledger(records, [], run_id=None)
    md = wfl.render(part, None, tag="v3.9.2")
    assert md.startswith("## Kalamari curation-triage — worklist snapshot at v3.9.2")
    assert "Tier 1 worklist" not in md


def test_render_markdown_still_defaults_to_the_tier_1_heading():
    # The monthly job passes no title; its heading must not have moved.
    part = ledger_mod.Partition(run_id="2026-08")
    assert ledger_mod.render_markdown(part).startswith(
        "## Kalamari curation-triage — Tier 1 worklist (run 2026-08)")


def test_an_unknown_run_id_fails_loudly_instead_of_writing_an_empty_worklist(tmp_path):
    path = _write(tmp_path / "ledger.ndjson", [_rec("a", "2026-08", "new")])
    assert wfl.main(["--ledger", path, "--run-id", "1999-01",
                     "--out", str(tmp_path / "w.tsv")]) == 1
    assert not (tmp_path / "w.tsv").exists()


def test_a_missing_ledger_fails_loudly(tmp_path):
    assert wfl.main(["--ledger", str(tmp_path / "nope.ndjson")]) == 1


# --------------------------------------------------------------------------- #
#  The committed ledger: invariants that survive the bot's monthly append       #
# --------------------------------------------------------------------------- #
def render_and_check_invariants(ledger_path, reviews_path, out_dir, tag="v3.9.2"):
    """Render ``ledger_path`` and assert what must be true of ANY state of the ledger.

    Stated as a mapping rather than as counts: the snapshot renders one row per live
    lead -- latest record per lead_key, minus the resolved ones and the ones an SME
    dismissed -- and every row carries its own record's run id. That holds on the seed
    file and on the same file after any number of monthly appends, which is the point
    (see the module docstring). Returns (rows, records) so a caller can add assertions
    about the particular ledger it passed in.
    """
    records = ledger_mod.load_ledger(ledger_path)
    assert records, "%s must parse into at least one record" % ledger_path
    out = os.path.join(str(out_dir), "worklist.tsv")
    md = os.path.join(str(out_dir), "worklist.md")
    assert wfl.main(["--ledger", str(ledger_path), "--reviews", str(reviews_path),
                     "--out", out, "--summary-md", md, "--tag", tag]) == 0
    with open(out, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    latest = ledger_mod.latest_by_key(records)
    live = {k for k, rec in latest.items() if rec.get("delta") != "resolved"}
    # A lead an SME touched may legitimately be suppressed; which of the three decisions
    # suppresses it is ledger.review_decision's business, tested above on synthetic data.
    reviewed = {(r.get("lead_key") or "").strip()
                for r in ledger_mod.load_reviews(str(reviews_path))
                if (r.get("decision") or "").strip().lower() in ledger_mod.REVIEW_DECISIONS}

    keys = [r["lead_key"] for r in rows]
    assert len(keys) == len(set(keys)), "a lead may appear on the worklist only once"
    assert set(keys) <= live, "a resolved or unknown lead reached the worklist"
    assert live - reviewed <= set(keys), "a live, unreviewed lead is missing from the worklist"

    for row in rows:
        key = row["lead_key"]
        # Per-row provenance: a snapshot mixes runs, so a row must carry the run id of
        # the record it was rendered from, never the newest one in the file.
        assert row["run_id"] == latest[key].get("run_id"), "row %s lost its run id" % key
        assert RUN_ID_RE.match(row["run_id"] or ""), "row %s has no run id" % key
        assert row["scientificName"], "row %s renders without a species" % key
        assert row["signal"] in KNOWN_SIGNALS
    with open(md, encoding="utf-8") as fh:
        assert fh.read().startswith(
            "## Kalamari curation-triage — worklist snapshot at %s" % tag)
    return rows, records


def test_the_committed_ledger_renders_every_live_lead_as_exactly_one_worklist_row(tmp_path):
    # The release asset is built from this file, so a rendering change must fail here.
    rows, _ = render_and_check_invariants(COMMITTED_LEDGER, COMMITTED_REVIEWS, tmp_path)
    assert rows, "the committed ledger must render a non-empty worklist"


def test_every_record_in_the_committed_ledger_has_the_shape_the_renderer_needs():
    """The ledger's shape is still pinned; only its length and history are not."""
    records = ledger_mod.load_ledger(COMMITTED_LEDGER)
    assert records
    for rec in records:
        assert rec["schema"] == 1
        assert rec["lead_key"]
        # J2 is the only signal about a species Kalamari does not carry yet ('todo:'
        # leads), so it is the only one allowed to have no unit id.
        assert rec["unit_id"] or rec["signal"] == SIGNAL_J2, rec["lead_key"]
        assert RUN_ID_RE.match(rec.get("run_id") or ""), rec["lead_key"]
        assert rec["signal"] in KNOWN_SIGNALS
        assert rec["delta"] in KNOWN_DELTAS
        assert rec["scientificName"], rec["lead_key"]
        assert rec["evidence_fingerprint"], rec["lead_key"]
    # Append-only, so both tiers stay in the file forever: the release snapshot has to
    # span them, and a seed that lost one of them would be a real regression.
    assert {SIGNAL_J1A, SIGNAL_T2} <= {rec["signal"] for rec in records}


def _next_run_id(records):
    """The month after the newest run in the file (run ids are YYYY-MM)."""
    year, month = (int(x) for x in max(r.get("run_id", "") for r in records).split("-"))
    return "%04d-%02d" % (year + month // 12, month % 12 + 1)


def _append_a_month(src_ledger, dst):
    """Grow a copy of ``src_ledger`` the way one monthly bot run grows the real one.

    That run re-records surviving leads (delta ``recurring``), records what newly fired
    (``new``), and records what stopped firing (``resolved``) -- then commits the file.
    Only a *copy* is ever grown here; the committed ledger is read-only for tests.

    Half the live leads are left alone on purpose, so the grown snapshot mixes run ids
    the way the real one does (Tier 1 and Tier 2 are appended by different jobs and do
    not always both fire). Returns (path, this month's run id, the keys left alone).
    """
    records = ledger_mod.load_ledger(src_ledger)
    run_id = _next_run_id(records)
    latest = ledger_mod.latest_by_key(records)
    live = [rec for rec in latest.values() if rec.get("delta") != "resolved"]
    grown = list(records)
    for rec in live[:len(live) // 2]:
        grown.append(dict(rec, run_id=run_id, delta="recurring"))
    grown.append(dict(live[-1], run_id=run_id, delta="resolved", actionable_flags=[]))
    grown.append(_rec("Newly_flagged~GCF_999999999.1", run_id, "new",
                      name="Newly_flagged"))
    grown.append(_rec("t2:Newly_flagged~GCF_999999999.1#ani", run_id, "new",
                      signal=SIGNAL_T2, name="Newly_flagged", epoch="skani@0.3.2"))
    untouched = {rec["lead_key"] for rec in live[len(live) // 2:-1]}
    return _write(dst, grown), run_id, untouched


def test_the_invariants_still_hold_after_the_monthly_append_the_bot_itself_commits(tmp_path):
    """Regression: the bot's own ledger commit must not turn this blocking check red.

    ``ledger.partition()`` appends a record for every surviving lead on every run, and
    curation-triage.yml commits the result, which fires curation-triage-tests.yml. So
    the assertions above run against a ledger that is longer than the seed, carries run
    ids the seed never had and can carry signals the seed never had. If any of them ever
    pins the seed history again, this test goes red here instead of in the required
    check on the bot's push.
    """
    grown, run_id, untouched = _append_a_month(COMMITTED_LEDGER, tmp_path / "grown.ndjson")
    rows, records = render_and_check_invariants(grown, COMMITTED_REVIEWS, tmp_path)
    seed = ledger_mod.load_ledger(COMMITTED_LEDGER)
    assert len(records) > len(seed), "the simulated month did not actually grow the ledger"
    by_key = {r["lead_key"]: r for r in rows}
    assert "Newly_flagged~GCF_999999999.1" in by_key         # this month's new lead renders
    assert "t2:Newly_flagged~GCF_999999999.1#ani" in by_key  # from either tier
    resolved = {r["lead_key"] for r in records if r["delta"] == "resolved"}
    assert resolved and not resolved & set(by_key)           # a resolved lead leaves the list
    # The leads this month did not re-record keep the run id they were last computed
    # under: a snapshot mixes runs and must not restamp them.
    assert untouched
    assert all(by_key[k]["run_id"] != run_id for k in untouched)


def test_a_ledger_that_lost_its_per_row_run_ids_still_fails_the_invariants(tmp_path):
    """The invariants replaced pinned counts, not real coverage: they must still bite."""
    records = ledger_mod.load_ledger(COMMITTED_LEDGER)
    no_run_id = _write(tmp_path / "no_run_id.ndjson",
                       [dict(r, run_id="") for r in records])
    with pytest.raises(AssertionError):
        render_and_check_invariants(no_run_id, COMMITTED_REVIEWS, tmp_path)
    no_species = _write(tmp_path / "no_species.ndjson",
                        [dict(r, scientificName="") for r in records])
    with pytest.raises(AssertionError):
        render_and_check_invariants(no_species, COMMITTED_REVIEWS, tmp_path)
