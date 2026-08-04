# -*- coding: utf-8 -*-
"""collate.py: merge sharded lead fragments into one worklist + ledger, and report
shard health (loud-but-non-blocking) -- all offline."""

import csv
import json

import collate
from leads import Lead, ANI_EPOCH


def _lead(unit, flags):
    return Lead(signal="J1a", lead_key=unit, unit_id=unit, scientific_name=unit,
                source_lane="NCBI-REF", epoch=ANI_EPOCH,
                actionable_flags=list(flags), flags=list(flags), summary=unit)


def _fragment(path, leads):
    with open(path, "w", encoding="utf-8") as fh:
        for l in leads:
            fh.write(json.dumps(l.to_dict(), sort_keys=True) + "\n")


def test_lead_roundtrip():
    l = _lead("U1", ["tax_check_failed"])
    assert Lead.from_dict(l.to_dict()) == l


def _statuses(path):
    with open(path, newline="") as fh:
        return {r["scientificName"]: r["status"] for r in csv.DictReader(fh, delimiter="\t")}


def test_collate_merges_shards(tmp_path):
    d = tmp_path / "frags"
    d.mkdir()
    _fragment(d / "leads.0.ndjson", [_lead("A", ["tax_check_failed"])])
    _fragment(d / "leads.1.ndjson", [_lead("B", ["assembly_suppressed"])])
    ledger = tmp_path / "ledger.ndjson"
    out = tmp_path / "worklist.tsv"
    rc = collate.main(["--leads-dir", str(d), "--shards-total", "2", "--run-id", "2026-08",
                       "--ledger", str(ledger), "--reviews", str(tmp_path / "none.tsv"),
                       "--out", str(out)])
    assert rc == 0                                   # both shards present
    assert _statuses(str(out)) == {"A": "new", "B": "new"}
    assert ledger.exists()


def test_collate_reports_missing_shard(tmp_path, capsys):
    d = tmp_path / "frags"
    d.mkdir()
    _fragment(d / "leads.0.ndjson", [_lead("A", ["tax_check_failed"])])
    # shard 1 crashed -> no leads.1.ndjson
    rc = collate.main(["--leads-dir", str(d), "--shards-total", "2", "--run-id", "2026-08",
                       "--ledger", str(tmp_path / "l.ndjson"),
                       "--reviews", str(tmp_path / "none.tsv"),
                       "--out", str(tmp_path / "w.tsv"), "--no-write"])
    assert rc == 1                                   # loud: a shard was missing
    assert "shards_failed=1/2" in capsys.readouterr().out


def test_collate_missing_shard_does_not_append_or_resolve(tmp_path):
    # A missing fragment (crashed shard) must NOT append the ledger -- otherwise the
    # crashed shard's units would be falsely marked 'resolved', corrupting the
    # append-only baseline. Collate must fail loud (exit 1) and touch nothing.
    import ledger as ledger_mod
    d = tmp_path / "frags"
    d.mkdir()
    _fragment(d / "leads.0.ndjson", [_lead("A", ["tax_check_failed"])])
    # shard 1 crashed -> no leads.1.ndjson
    ledger = tmp_path / "ledger.ndjson"
    ledger_mod.append_ledger(str(ledger),
                             [_lead("B", ["tax_check_inconclusive"]).as_record("2026-07", "new")])
    before = ledger.read_text()
    out = tmp_path / "w.tsv"
    rc = collate.main(["--leads-dir", str(d), "--shards-total", "2", "--run-id", "2026-08",
                       "--ledger", str(ledger), "--reviews", str(tmp_path / "none.tsv"),
                       "--out", str(out)])
    assert rc == 1
    assert ledger.read_text() == before      # append-only baseline untouched
    assert not out.exists()                  # no misleading partial worklist


def test_collate_corrupt_fragment_is_failed_shard(tmp_path, capsys):
    # A truncated/corrupt fragment must be counted as a failed shard, not crash collate.
    d = tmp_path / "frags"
    d.mkdir()
    _fragment(d / "leads.0.ndjson", [_lead("A", ["tax_check_failed"])])
    (d / "leads.1.ndjson").write_text('{"signal":"J1a","lead_key":"B" TRUNCATED\n')
    ledger = tmp_path / "ledger.ndjson"
    rc = collate.main(["--leads-dir", str(d), "--shards-total", "2", "--run-id", "2026-08",
                       "--ledger", str(ledger), "--reviews", str(tmp_path / "none.tsv"),
                       "--out", str(tmp_path / "w.tsv")])
    assert rc == 1
    assert "shards_failed=1/2" in capsys.readouterr().out
    assert not ledger.exists()               # nothing appended on an incomplete run


def test_collate_unreadable_fragment_is_failed_shard(tmp_path, capsys):
    # A fragment that EXISTS but cannot be read (bad permissions on the runner) passes
    # the missing-file check, so it has to be caught where it is opened -- otherwise
    # collate crashes, or worse, silently treats the run as complete.
    import os
    d = tmp_path / "frags"
    d.mkdir()
    _fragment(d / "leads.0.ndjson", [_lead("A", ["tax_check_failed"])])
    bad = d / "leads.1.ndjson"
    _fragment(bad, [_lead("B", ["assembly_suppressed"])])
    os.chmod(bad, 0o000)
    if os.access(str(bad), os.R_OK):          # running as root: chmod proves nothing
        os.chmod(bad, 0o644)
        return
    try:
        ledger = tmp_path / "ledger.ndjson"
        rc = collate.main(["--leads-dir", str(d), "--shards-total", "2",
                           "--run-id", "2026-08", "--ledger", str(ledger),
                           "--reviews", str(tmp_path / "none.tsv"),
                           "--out", str(tmp_path / "w.tsv")])
        assert rc == 1
        assert "shards_failed=1/2" in capsys.readouterr().out
        assert not ledger.exists()
    finally:
        os.chmod(bad, 0o644)


# --------------------------------------------------------------------------- #
#  --shards-total is a required input, not a guess                             #
# --------------------------------------------------------------------------- #
def test_collate_requires_an_explicit_shard_count(tmp_path):
    d = tmp_path / "frags"
    d.mkdir()
    _fragment(d / "leads.0.ndjson", [_lead("A", ["tax_check_failed"])])
    try:
        collate.main(["--leads-dir", str(d), "--run-id", "2026-08",
                      "--ledger", str(tmp_path / "l.ndjson"),
                      "--reviews", str(tmp_path / "none.tsv"),
                      "--out", str(tmp_path / "w.tsv")])
        assert False, "expected collate to refuse an unstated shard count"
    except SystemExit as e:
        assert "--shards-total is required" in str(e)


def test_an_unstated_shard_count_can_no_longer_falsely_resolve_the_other_shards(tmp_path):
    # The defect this guards: --shards-total used to default to 1, so pointing collate
    # at a 2-shard directory read only leads.0.ndjson, saw nothing missing, called the
    # run complete, and wrote 'resolved' records for every unit the other shard held --
    # permanent false entries in an append-only ledger.
    import ledger as ledger_mod
    d = tmp_path / "frags"
    d.mkdir()
    _fragment(d / "leads.0.ndjson", [_lead("A", ["tax_check_failed"])])
    _fragment(d / "leads.1.ndjson", [_lead("B", ["assembly_suppressed"])])
    ledger = tmp_path / "ledger.ndjson"
    ledger_mod.append_ledger(str(ledger),
                             [_lead("B", ["assembly_suppressed"]).as_record("2026-07", "new")])
    before = ledger.read_text()
    try:
        collate.main(["--leads-dir", str(d), "--run-id", "2026-08",
                      "--ledger", str(ledger), "--reviews", str(tmp_path / "none.tsv"),
                      "--out", str(tmp_path / "w.tsv")])
        assert False, "expected collate to refuse rather than read only shard 0"
    except SystemExit:
        pass
    assert ledger.read_text() == before        # B was never falsely 'resolved'


def test_collate_cross_checks_an_explicit_fragment_list_against_the_shard_count(tmp_path):
    d = tmp_path / "frags"
    d.mkdir()
    _fragment(d / "leads.0.ndjson", [_lead("A", ["tax_check_failed"])])
    try:
        collate.main(["--leads", str(d / "leads.0.ndjson"), "--shards-total", "4",
                      "--run-id", "2026-08", "--ledger", str(tmp_path / "l.ndjson"),
                      "--reviews", str(tmp_path / "none.tsv"),
                      "--out", str(tmp_path / "w.tsv")])
        assert False, "expected collate to refuse a mismatched --shards-total"
    except SystemExit as e:
        assert "1 fragment path" in str(e)


def test_collate_accepts_an_explicit_fragment_list_without_a_shard_count(tmp_path):
    # With --leads the path count IS the shard count, so stating it is optional.
    d = tmp_path / "frags"
    d.mkdir()
    _fragment(d / "leads.0.ndjson", [_lead("A", ["tax_check_failed"])])
    out = tmp_path / "worklist.tsv"
    rc = collate.main(["--leads", str(d / "leads.0.ndjson"), "--run-id", "2026-08",
                       "--ledger", str(tmp_path / "l.ndjson"),
                       "--reviews", str(tmp_path / "none.tsv"), "--out", str(out)])
    assert rc == 0
    assert _statuses(str(out)) == {"A": "new"}


def test_collate_stamps_the_same_resolved_epochs_as_a_single_process_run(tmp_path):
    # collate recomputes tool_versions from the committed caches; a shard fragment
    # carries no epochs, so this is what keeps the two paths' ledger records identical.
    import datasets_summary as ds
    import ledger as ledger_mod
    import tier1_metadata
    dcache = tmp_path / "ds.json"
    c = ds.DatasetsSummaryCache(cache_path=str(dcache), fetcher=lambda t, a, r: [])
    c.build([100], released_after="", version_epoch="ncbi-datasets-cli@18.34.0")
    c.save()
    d = tmp_path / "frags"
    d.mkdir()
    _fragment(d / "leads.0.ndjson", [_lead("A", ["tax_check_failed"])])
    ledger = tmp_path / "ledger.ndjson"
    rc = collate.main(["--leads-dir", str(d), "--shards-total", "1", "--run-id", "2026-08",
                       "--ledger", str(ledger), "--reviews", str(tmp_path / "none.tsv"),
                       "--out", str(tmp_path / "w.tsv"), "--datasets-cache", str(dcache),
                       "--lpsn-cache", str(tmp_path / "no-lpsn.json")])
    assert rc == 0
    (rec,) = ledger_mod.load_ledger(str(ledger))
    assert rec["tool_versions"] == tier1_metadata.tool_versions(
        str(dcache), str(tmp_path / "no-lpsn.json"))
    assert rec["tool_versions"]["ncbi_datasets_cli"] == "ncbi-datasets-cli@18.34.0"


def test_collate_global_delta_not_per_shard(tmp_path):
    # A unit absent from THIS run's fragments but actionable last run must be
    # 'resolved' by collate (the global view), which a per-shard delta could not see.
    d = tmp_path / "frags"
    d.mkdir()
    _fragment(d / "leads.0.ndjson", [_lead("A", ["tax_check_failed"])])
    _fragment(d / "leads.1.ndjson", [])              # shard ran, found nothing
    ledger = tmp_path / "ledger.ndjson"
    # seed a prior run where B was actionable
    import ledger as ledger_mod
    ledger_mod.append_ledger(str(ledger), [_lead("B", ["tax_check_inconclusive"]).as_record("2026-07", "new")])
    collate.main(["--leads-dir", str(d), "--shards-total", "2", "--run-id", "2026-08",
                  "--ledger", str(ledger), "--reviews", str(tmp_path / "none.tsv"),
                  "--out", str(tmp_path / "w.tsv")])
    recs = ledger_mod.load_ledger(str(ledger))
    resolved = [r for r in recs if r["run_id"] == "2026-08" and r["delta"] == "resolved"]
    assert [r["lead_key"] for r in resolved] == ["B"]
