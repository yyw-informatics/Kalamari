# -*- coding: utf-8 -*-
"""tier1_metadata.main() end-to-end, fully offline.

Drives the whole monthly loop through the real CLI: build the cache from local
NCBI-shaped files (no network), run J1a, diff against the growing ledger, apply
reviews, and render the worklist -- across successive runs so the delta, the review
suppression, and the resolved transition are exercised together.  Tiny synthetic
fixtures reused from the coverage-probe test builders."""

import csv
import os

import tier1_metadata
from test_ani_reports import ani_line, summary_line, ANI_HEADER, SUMMARY_HEADER_1, SUMMARY_HEADER_2
from test_coverage_cli import POLICY_HEADER, _policy_row


def _chrom(tmp_path):
    chrom = tmp_path / "chromosomes.tsv"
    chrom.write_text(
        "scientificName\tnuccoreAcc\ttaxid\tparent\tsource\n"
        "Foo\tCP1\t100\t1\tNCBI-REF\n"
        "Bar\tCP2\t200\t1\tCDC SME\n"
        "Baz\tCP3\t300\t1\tNCBI-REF\n"
    )
    return chrom


def _policy(tmp_path):
    policy = tmp_path / "taxon_policy.tsv"
    policy.write_text("\n".join([
        POLICY_HEADER,
        _policy_row("Foo~GCF_1.1", "Foo", 100, "CP1", "GCF_1.1"),
        _policy_row("Bar~GCF_2.1", "Bar", 200, "CP2", "GCF_2.1"),
        _policy_row("Baz~GCF_3.1", "Baz", 300, "CP3", "GCF_3.1", regime="B", resolver="ANI"),
    ]) + "\n")
    return policy


def _summaries(tmp_path):
    refseq = tmp_path / "refseq.txt"
    refseq.write_text("".join([
        SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
        summary_line(acc="GCF_1.1", species_taxid=100, paired="GCA_1.1"),
        summary_line(acc="GCF_2.1", species_taxid=200, paired="GCA_2.1"),
        summary_line(acc="GCF_3.1", species_taxid=300, paired="GCA_3.1"),
    ]))
    genbank = tmp_path / "genbank.txt"
    genbank.write_text("".join([
        SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
        summary_line(acc="GCA_1.1", species_taxid=100, paired="GCF_1.1", source_db="genbank"),
        summary_line(acc="GCA_2.1", species_taxid=200, paired="GCF_2.1", source_db="genbank"),
        summary_line(acc="GCA_3.1", species_taxid=300, paired="GCF_3.1", source_db="genbank"),
    ]))
    return refseq, genbank


def _ani(tmp_path, name, bar_check="Failed", foo_check="OK", baz_check="Failed"):
    ani = tmp_path / name
    ani.write_text("".join([
        ANI_HEADER + "\n",
        ani_line(genbank="GCA_1.1", refseq="GCF_1.1", species_taxid=100, best_taxid=100,
                 tax_check=foo_check),
        ani_line(genbank="GCA_2.1", refseq="GCF_2.1", species_taxid=200, best_taxid=200,
                 tax_check=bar_check),
        # Baz is confounded (regime B): a Failed tax-check demotes to non-actionable.
        ani_line(genbank="GCA_3.1", refseq="GCF_3.1", species_taxid=300, best_taxid=888,
                 best_name="Sister sp.", best_status="mismatch", tax_check=baz_check),
    ]))
    return ani


def _run(tmp_path, ani, run_id, ledger, reviews=None):
    cache = tmp_path / "cache.json"
    out = tmp_path / ("worklist_%s.tsv" % run_id)
    refseq, genbank = _summaries(tmp_path)
    argv = ["--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
            "--cache", str(cache), "--ledger", str(ledger), "--out", str(out),
            "--run-id", run_id,
            "--ani-file", str(ani), "--refseq-file", str(refseq), "--genbank-file", str(genbank)]
    if reviews is not None:
        argv += ["--reviews", str(reviews)]
    rc = tier1_metadata.main(argv)
    assert rc == 0
    return out


def _statuses(worklist_path):
    with open(worklist_path, newline="") as fh:
        return {r["scientificName"]: r["status"] for r in csv.DictReader(fh, delimiter="\t")}


def _worklist(worklist_path):
    with open(worklist_path, newline="") as fh:
        return {r["scientificName"]: r for r in csv.DictReader(fh, delimiter="\t")}


# --------------------------------------------------------------------------- #
def test_first_run_surfaces_actionable_only(tmp_path):
    ledger = tmp_path / "ledger.ndjson"
    out = _run(tmp_path, _ani(tmp_path, "ani1.txt"), "2026-08", ledger)
    st = _statuses(out)
    # Foo clean -> absent; Bar Failed -> new; Baz confounded+Failed -> demoted, absent.
    assert st == {"Bar": "new"}
    assert ledger.exists()


def test_second_run_is_standing_not_new(tmp_path):
    ledger = tmp_path / "ledger.ndjson"
    _run(tmp_path, _ani(tmp_path, "ani1.txt"), "2026-08", ledger)
    out2 = _run(tmp_path, _ani(tmp_path, "ani2.txt"), "2026-09", ledger)
    assert _statuses(out2) == {"Bar": "standing"}      # same evidence -> not new


def test_newly_actionable_pick_surfaces_as_new(tmp_path):
    ledger = tmp_path / "ledger.ndjson"
    _run(tmp_path, _ani(tmp_path, "ani1.txt", foo_check="OK"), "2026-08", ledger)
    # next run: Foo's tax-check flips OK -> Failed => Foo is newly actionable.
    out2 = _run(tmp_path, _ani(tmp_path, "ani2.txt", foo_check="Failed"), "2026-09", ledger)
    st = _statuses(out2)
    assert st.get("Foo") == "new"
    assert st.get("Bar") == "standing"


def test_resolved_transition(tmp_path):
    ledger = tmp_path / "ledger.ndjson"
    _run(tmp_path, _ani(tmp_path, "ani1.txt", bar_check="Failed"), "2026-08", ledger)
    # Bar's tax-check recovers to OK => no longer actionable => resolved, off the worklist.
    out2 = _run(tmp_path, _ani(tmp_path, "ani2.txt", bar_check="OK"), "2026-09", ledger)
    assert "Bar" not in _statuses(out2)
    import ledger as ledger_mod
    recs = ledger_mod.load_ledger(str(ledger))
    resolved = [r for r in recs if r["run_id"] == "2026-09" and r["delta"] == "resolved"]
    assert [r["lead_key"] for r in resolved] == ["Bar~GCF_2.1"]


def test_j1b_offline_surfaces_better_rep(tmp_path):
    # Prebuild a datasets cache offline (injected fetcher) for the tracked species, then
    # run the driver --with-j1b --offline: Foo (NCBI-REF, taxid 100) gets a better-rep
    # lead for a newly-designated reference genome.
    import datasets_summary as ds
    from test_datasets_summary import report_line
    dcache = tmp_path / "ds.json"
    fetcher = lambda taxid, after, ref_only: (
        [report_line("GCF_88.1", 100, category="reference genome")] if taxid == 100 else [])
    c = ds.DatasetsSummaryCache(cache_path=str(dcache), fetcher=fetcher)
    # tracked species = NCBI-REF, non-confounded picks: Foo(100). (Bar is SME; Baz confounded.)
    # J1b asks narrowly, so the committed cache records taxid 100 as reference-only.
    c.build([100], released_after="", version_epoch="ncbi-datasets-cli@16.0.0",
            reference_only={100})
    c.save()

    ledger = tmp_path / "ledger.ndjson"
    cache = tmp_path / "cache.json"
    out = tmp_path / "worklist.tsv"
    refseq, genbank = _summaries(tmp_path)
    base = ["--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
            "--cache", str(cache), "--ledger", str(ledger), "--out", str(out),
            "--run-id", "2026-08"]
    # Step 1: build the ANI cache from the local files (no ledger written).
    tier1_metadata.main(base + ["--no-write", "--ani-file", str(_ani(tmp_path, "a.txt")),
                                "--refseq-file", str(refseq), "--genbank-file", str(genbank)])
    # Step 2: fully offline run (ANI cache + datasets cache both committed) with J1b.
    rc = tier1_metadata.main(base + ["--offline", "--with-j1b", "--datasets-cache", str(dcache)])
    assert rc == 0
    wl = _worklist(out)
    assert "Foo" in wl and wl["Foo"]["signal"] == "J1b"
    assert wl["Foo"]["candidate_accession"] == "GCF_88.1"


def test_j2_offline_surfaces_wishlist(tmp_path):
    # Prebuild a datasets cache for a wishlist taxid, then run --with-j2 --offline.
    import datasets_summary as ds
    from test_datasets_summary import report_line
    dcache = tmp_path / "ds.json"
    c = ds.DatasetsSummaryCache(
        cache_path=str(dcache),
        fetcher=lambda taxid, after, ref_only:
            [report_line("GCF_777.1", taxid, category="reference genome")])
    c.build([28197], released_after="", version_epoch="ncbi-datasets-cli@16.0.0")
    c.save()
    todo = tmp_path / "todo.tsv"
    todo.write_text("Arcobacter butzleri\tXXXXXX\t28197\t28196\n")

    ledger = tmp_path / "ledger.ndjson"
    cache = tmp_path / "cache.json"
    out = tmp_path / "worklist.tsv"
    refseq, genbank = _summaries(tmp_path)
    base = ["--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
            "--cache", str(cache), "--ledger", str(ledger), "--out", str(out), "--run-id", "2026-08"]
    tier1_metadata.main(base + ["--no-write", "--ani-file", str(_ani(tmp_path, "a.txt")),
                                "--refseq-file", str(refseq), "--genbank-file", str(genbank)])
    rc = tier1_metadata.main(base + ["--offline", "--with-j2",
                                     "--todo", str(todo), "--datasets-cache", str(dcache)])
    assert rc == 0
    wl = _worklist(out)
    assert "Arcobacter butzleri" in wl
    assert wl["Arcobacter butzleri"]["signal"] == "J2"
    assert wl["Arcobacter butzleri"]["candidate_accession"] == "GCF_777.1"


def test_datasets_cache_union_covers_both_j1b_and_j2(tmp_path, monkeypatch):
    # A combined non-offline --with-j1b --with-j2 build must produce ONE datasets cache
    # covering BOTH J1b's tracked species AND J2's wishlist taxids (the two feeds share
    # the cache file; a per-feed build would clobber the other).
    import datasets_summary as ds
    from test_datasets_summary import report_line
    asked = {}

    def fake_run(taxid, after, datasets_bin="datasets", reference_only=False):
        asked[taxid] = reference_only
        return [report_line("GCF_%d.1" % taxid, taxid, category="reference genome")]

    monkeypatch.setattr(ds, "run_datasets", fake_run)
    monkeypatch.setattr(ds, "datasets_version",
                        lambda datasets_bin="datasets": "ncbi-datasets-cli@test")

    todo = tmp_path / "todo.tsv"
    todo.write_text("Wanted species\tXXXXXX\t999\t1\n")   # J2 wishlist taxid 999
    dcache = tmp_path / "ds.json"
    cache = tmp_path / "cache.json"
    refseq, genbank = _summaries(tmp_path)
    # Foo(100) is NCBI-REF -> a J1b tracked species; Bar(200) SME, Baz(300) confounded.
    tier1_metadata.main([
        "--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
        "--cache", str(cache), "--ledger", str(tmp_path / "ledger.ndjson"),
        "--out", str(tmp_path / "w.tsv"), "--run-id", "2026-08",
        "--with-j1b", "--with-j2", "--todo", str(todo), "--datasets-cache", str(dcache),
        "--ani-file", str(_ani(tmp_path, "a.txt")),
        "--refseq-file", str(refseq), "--genbank-file", str(genbank),
    ])
    c = ds.DatasetsSummaryCache(cache_path=str(dcache))
    assert c.covers([100, 999], reference_only={100})   # BOTH feeds' taxids in one cache
    # ...and each feed got the query IT needs: narrow for J1b, wide for J2's wishlist.
    assert asked == {100: True, 999: False}
    assert c.mode_of(100) == ds.QUERY_REFERENCE and c.mode_of(999) == ds.QUERY_ALL


def test_a_taxid_wanted_by_both_feeds_gets_the_wide_query(tmp_path, monkeypatch):
    # J1b would be happy with --reference, J2's wishlist would not: the wide answer
    # serves both, the narrow one blinds J2. So the full query must win the overlap.
    import datasets_summary as ds
    from test_datasets_summary import report_line
    asked = {}

    def fake_run(taxid, after, datasets_bin="datasets", reference_only=False):
        asked[taxid] = reference_only
        return [report_line("GCF_%d.1" % taxid, taxid, category="reference genome")]

    monkeypatch.setattr(ds, "run_datasets", fake_run)
    monkeypatch.setattr(ds, "datasets_version",
                        lambda datasets_bin="datasets": "ncbi-datasets-cli@test")
    todo = tmp_path / "todo.tsv"
    todo.write_text("Foo bar\tXXXXXX\t100\t1\n")     # same taxid as the NCBI-REF pick Foo
    dcache = tmp_path / "ds.json"
    refseq, genbank = _summaries(tmp_path)
    tier1_metadata.main([
        "--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
        "--cache", str(tmp_path / "cache.json"), "--ledger", str(tmp_path / "l.ndjson"),
        "--out", str(tmp_path / "w.tsv"), "--run-id", "2026-08", "--no-write",
        "--with-j1b", "--with-j2", "--todo", str(todo), "--datasets-cache", str(dcache),
        "--ani-file", str(_ani(tmp_path, "a.txt")),
        "--refseq-file", str(refseq), "--genbank-file", str(genbank),
    ])
    assert asked == {100: False}
    assert ds.DatasetsSummaryCache(cache_path=str(dcache)).mode_of(100) == ds.QUERY_ALL


def test_offline_refuses_a_datasets_cache_that_is_only_reference_wide(tmp_path):
    # The regression this guards: a cache built for J1b silently answering J2, which
    # would report "no genome available yet" for every wishlist species whose only
    # genome is not NCBI's designated reference.
    import datasets_summary as ds
    dcache = tmp_path / "ds.json"
    c = ds.DatasetsSummaryCache(cache_path=str(dcache), fetcher=lambda t, a, r: [])
    c.build([28197], released_after="", reference_only={28197})
    c.save()
    todo = tmp_path / "todo.tsv"
    todo.write_text("Arcobacter butzleri\tXXXXXX\t28197\t28196\n")
    refseq, genbank = _summaries(tmp_path)
    base = ["--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
            "--cache", str(tmp_path / "cache.json"), "--ledger", str(tmp_path / "l.ndjson"),
            "--out", str(tmp_path / "w.tsv"), "--run-id", "2026-08"]
    tier1_metadata.main(base + ["--no-write", "--ani-file", str(_ani(tmp_path, "a.txt")),
                                "--refseq-file", str(refseq), "--genbank-file", str(genbank)])
    try:
        tier1_metadata.main(base + ["--offline", "--with-j2", "--todo", str(todo),
                                    "--datasets-cache", str(dcache)])
        assert False, "expected --offline to refuse a reference-only cache entry"
    except SystemExit as e:
        assert "reference-only" in str(e)


def test_offline_refuses_a_datasets_cache_narrowed_by_a_released_after_window(tmp_path):
    # The regression this guards, seen for real on the committed cache: a
    # --released-after refresh REPLACES every candidate list with the window's contents,
    # so the same 170 taxids are still present but nearly empty. The old gate checked
    # only presence, so the run exited 0, printed "cache covers 170 taxid(s)", and
    # dropped 77 standing J1b/J2 leads off the worklist with no 'resolved' record --
    # zero leads instead of a coverage failure.
    import datasets_summary as ds
    from test_datasets_summary import report_line
    dcache = tmp_path / "ds.json"
    c = ds.DatasetsSummaryCache(
        cache_path=str(dcache),
        # Inside the window NCBI returns nothing; the real reference is older than it.
        fetcher=lambda taxid, after, ref_only:
            [] if after else [report_line("GCF_88.1", 100, category="reference genome")])
    c.build([100], released_after="07/01/2026", reference_only={100})
    c.save()
    assert c.candidates(100) == []          # the narrowing really did gut the entry

    refseq, genbank = _summaries(tmp_path)
    base = ["--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
            "--cache", str(tmp_path / "cache.json"), "--ledger", str(tmp_path / "l.ndjson"),
            "--out", str(tmp_path / "w.tsv"), "--run-id", "2026-08",
            "--datasets-cache", str(dcache)]
    tier1_metadata.main(base + ["--no-write", "--ani-file", str(_ani(tmp_path, "a.txt")),
                                "--refseq-file", str(refseq), "--genbank-file", str(genbank)])
    try:
        tier1_metadata.main(base + ["--offline", "--with-j1b"])
        assert False, "expected --offline to refuse a window-narrowed cache"
    except SystemExit as e:
        assert "cached only since 07/01/2026" in str(e)
        assert "Rebuild without --offline" in str(e)
    # ...and the same cache still answers the run that asked for exactly that window.
    assert tier1_metadata.main(base + ["--offline", "--with-j1b", "--no-write",
                                       "--released-after", "07/01/2026"]) == 0


def test_a_wishlist_taxid_is_fetched_unwindowed_even_when_released_after_is_given(
        tmp_path, monkeypatch):
    # "Does this species have a genome at all?" is a standing fact, not a delta, so a
    # window gives J2 a WRONG answer rather than a narrower one. The window belongs to
    # J1b's reference feed alone.
    import datasets_summary as ds
    from test_datasets_summary import report_line
    asked = {}

    def fake_run(taxid, after, datasets_bin="datasets", reference_only=False):
        asked[taxid] = after
        return [report_line("GCF_%d.1" % taxid, taxid, category="reference genome")]

    monkeypatch.setattr(ds, "run_datasets", fake_run)
    monkeypatch.setattr(ds, "datasets_version",
                        lambda datasets_bin="datasets": "ncbi-datasets-cli@test")
    todo = tmp_path / "todo.tsv"
    todo.write_text("Wanted species\tXXXXXX\t999\t1\n")      # J2 wishlist taxid 999
    dcache = tmp_path / "ds.json"
    refseq, genbank = _summaries(tmp_path)
    tier1_metadata.main([
        "--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
        "--cache", str(tmp_path / "cache.json"), "--ledger", str(tmp_path / "l.ndjson"),
        "--out", str(tmp_path / "w.tsv"), "--run-id", "2026-08", "--no-write",
        "--with-j1b", "--with-j2", "--todo", str(todo), "--datasets-cache", str(dcache),
        "--released-after", "07/01/2026",
        "--ani-file", str(_ani(tmp_path, "a.txt")),
        "--refseq-file", str(refseq), "--genbank-file", str(genbank),
    ])
    assert asked == {100: "07/01/2026", 999: ""}    # J1b windowed, J2 wishlist not
    c = ds.DatasetsSummaryCache(cache_path=str(dcache))
    assert c.window_of(999) == ds.WINDOW_ALL
    # ...so a later J2-only offline run (which never passes a window) still passes.
    assert c.covers([999])


def test_tool_versions_stamp_the_resolved_cache_epochs(tmp_path):
    # A ledger record claiming 'ncbi-datasets-cli@unset' is worse than none: it is the
    # field an SME reads to decide whether a verdict is stale.
    import datasets_summary as ds
    dcache = tmp_path / "ds.json"
    c = ds.DatasetsSummaryCache(cache_path=str(dcache), fetcher=lambda t, a, r: [])
    c.build([100], released_after="", version_epoch="ncbi-datasets-cli@18.34.0")
    c.save()
    tv = tier1_metadata.tool_versions(str(dcache), str(tmp_path / "no-lpsn.json"))
    assert tv["ncbi_datasets_cli"] == "ncbi-datasets-cli@18.34.0"
    assert tv["lpsn"] == "lpsn@unset"          # no cache -> the honest fallback
    # An absent datasets cache falls back too, rather than crashing the run.
    assert tier1_metadata.tool_versions(str(tmp_path / "nope.json"),
                                        str(tmp_path / "nope2.json"))[
        "ncbi_datasets_cli"] == "ncbi-datasets-cli@unset"


def test_the_ledger_record_carries_the_resolved_datasets_epoch(tmp_path):
    import datasets_summary as ds
    import ledger as ledger_mod
    from test_datasets_summary import report_line
    dcache = tmp_path / "ds.json"
    c = ds.DatasetsSummaryCache(
        cache_path=str(dcache),
        fetcher=lambda t, a, r: [report_line("GCF_88.1", 100, category="reference genome")])
    c.build([100], released_after="", version_epoch="ncbi-datasets-cli@18.34.0",
            reference_only={100})
    c.save()
    ledger = tmp_path / "ledger.ndjson"
    refseq, genbank = _summaries(tmp_path)
    base = ["--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
            "--cache", str(tmp_path / "cache.json"), "--ledger", str(ledger),
            "--out", str(tmp_path / "w.tsv"), "--run-id", "2026-08",
            "--datasets-cache", str(dcache)]
    tier1_metadata.main(base + ["--no-write", "--ani-file", str(_ani(tmp_path, "a.txt")),
                                "--refseq-file", str(refseq), "--genbank-file", str(genbank)])
    tier1_metadata.main(base + ["--offline", "--with-j1b"])
    recs = ledger_mod.load_ledger(str(ledger))
    assert all(r["tool_versions"]["ncbi_datasets_cli"] == "ncbi-datasets-cli@18.34.0"
               for r in recs)


def test_refresh_datasets_cache_writes_only_the_datasets_cache(tmp_path, monkeypatch):
    # The datasets refresh must not drag the 2.7 GB assembly-reports rebuild, the
    # ledger append or the worklist along with it.
    import datasets_summary as ds
    from test_datasets_summary import report_line
    monkeypatch.setattr(ds, "run_datasets",
                        lambda taxid, after, datasets_bin="datasets", reference_only=False:
                            [report_line("GCF_%d.1" % taxid, taxid, category="reference genome")])
    monkeypatch.setattr(ds, "datasets_version",
                        lambda datasets_bin="datasets": "ncbi-datasets-cli@test")
    dcache = tmp_path / "ds.json"
    acache = tmp_path / "assembly.json"
    ledger = tmp_path / "ledger.ndjson"
    out = tmp_path / "w.tsv"
    rc = tier1_metadata.main([
        "--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
        "--cache", str(acache), "--ledger", str(ledger), "--out", str(out),
        "--refresh-datasets-cache", "--with-j1b", "--datasets-cache", str(dcache)])
    assert rc == 0
    assert dcache.exists()
    assert not acache.exists() and not ledger.exists() and not out.exists()


def test_refresh_datasets_cache_refuses_offline(tmp_path):
    try:
        tier1_metadata.main([
            "--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
            "--cache", str(tmp_path / "a.json"), "--refresh-datasets-cache",
            "--with-j1b", "--offline"])
        assert False, "expected --refresh-datasets-cache --offline to be refused"
    except SystemExit as e:
        assert "needs the network" in str(e)


def test_with_lpsn_offline_explains_how_to_build_the_missing_cache(tmp_path):
    # LPSN has no committed cache and no unattended download, so the failure has to
    # hand the operator the whole recipe instead of "rebuild without --offline".
    refseq, genbank = _summaries(tmp_path)
    base = ["--policy", str(_policy(tmp_path)), "--chromosomes", str(_chrom(tmp_path)),
            "--cache", str(tmp_path / "cache.json"), "--ledger", str(tmp_path / "l.ndjson"),
            "--out", str(tmp_path / "w.tsv"), "--run-id", "2026-08",
            "--lpsn-cache", str(tmp_path / "absent-lpsn.json")]
    tier1_metadata.main(base + ["--no-write", "--ani-file", str(_ani(tmp_path, "a.txt")),
                                "--refseq-file", str(refseq), "--genbank-file", str(genbank)])
    try:
        tier1_metadata.main(base + ["--offline", "--with-lpsn"])
        assert False, "expected --with-lpsn without a cache to be refused"
    except SystemExit as e:
        msg = str(e)
        assert "--lpsn-csv" in msg and "sign in" in msg and "--lpsn-release" in msg


def test_j2_lpsn_sharding_uses_global_universe(tmp_path):
    # Two same-genus picks split across shards; a species Kalamari HOLDS in the other
    # shard must NOT be reported as an LPSN gap (the gap is a global-coverage fact).
    import json
    import lpsn as lpsn_mod
    from test_j2_feeders import LPSN_HEADER

    policy = tmp_path / "policy.tsv"
    policy.write_text("\n".join([
        POLICY_HEADER,
        _policy_row("Ab~GCF_1.1", "Arcobacter_butzleri", 100, "CP1", "GCF_1.1",
                    species_name="Arcobacter butzleri"),
        _policy_row("Ac~GCF_2.1", "Arcobacter_cryaerophilus", 200, "CP2", "GCF_2.1",
                    species_name="Arcobacter cryaerophilus"),
    ]) + "\n")
    chrom = tmp_path / "chromosomes.tsv"
    chrom.write_text("scientificName\tnuccoreAcc\ttaxid\tparent\tsource\n"
                     "Arcobacter_butzleri\tCP1\t100\t1\tNCBI-REF\n"
                     "Arcobacter_cryaerophilus\tCP2\t200\t1\tNCBI-REF\n")
    ani = tmp_path / "ani.txt"
    ani.write_text("".join([
        ANI_HEADER + "\n",
        ani_line(genbank="GCA_1.1", refseq="GCF_1.1", species_taxid=100, best_taxid=100),
        ani_line(genbank="GCA_2.1", refseq="GCF_2.1", species_taxid=200, best_taxid=200),
    ]))
    refseq = tmp_path / "refseq.txt"
    refseq.write_text("".join([SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
                               summary_line(acc="GCF_1.1", species_taxid=100, paired="GCA_1.1"),
                               summary_line(acc="GCF_2.1", species_taxid=200, paired="GCA_2.1")]))
    genbank = tmp_path / "genbank.txt"
    genbank.write_text("".join([SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
                                summary_line(acc="GCA_1.1", species_taxid=100, paired="GCF_1.1", source_db="genbank"),
                                summary_line(acc="GCA_2.1", species_taxid=200, paired="GCF_2.1", source_db="genbank")]))

    # LPSN cache: two HELD species + one genuine gap (venerupis).
    lcache = tmp_path / "lpsn.json"
    rows = [["Arcobacter", "butzleri", "", "correct name", "x", "1991-01-01"],
            ["Arcobacter", "cryaerophilus", "", "correct name", "x", "1992-01-01"],
            ["Arcobacter", "venerupis", "", "correct name", "x", "2012-01-01"]]
    lc = lpsn_mod.LpsnCache(cache_path=str(lcache),
                            opener=lambda: [LPSN_HEADER + "\n"] + [",".join(r) + "\n" for r in rows])
    lc.build(release="LPSN-2026-06")
    lc.save()

    cache = tmp_path / "cache.json"
    base = ["--policy", str(policy), "--chromosomes", str(chrom), "--cache", str(cache),
            "--ledger", str(tmp_path / "ledger.ndjson"),
            "--run-id", "2026-08", "--shard", "0", "--shards-total", "2"]
    # Step 1: build the ANI cache for shard 0's slice.
    tier1_metadata.main(base + ["--no-write", "--ani-file", str(ani),
                                "--refseq-file", str(refseq), "--genbank-file", str(genbank)])
    # Step 2: offline shard-0 run with the LPSN feeder, emitting raw leads.
    frag = tmp_path / "leads.0.ndjson"
    rc = tier1_metadata.main(base + ["--offline", "--with-lpsn", "--lpsn-cache", str(lcache),
                                     "--emit-leads", str(frag)])
    assert rc == 0
    lpsn_species = {json.loads(l)["scientific_name"]
                    for l in frag.read_text().splitlines() if l.strip()
                    and json.loads(l)["source_lane"] == "lpsn"}
    # cryaerophilus is held (in shard 1) -> must NOT be a gap; venerupis genuinely is.
    assert "Arcobacter cryaerophilus" not in lpsn_species
    assert "Arcobacter venerupis" in lpsn_species

    # ...and the LPSN leads stamp 'lpsn@<release>', never the bare release tag.
    epochs = {json.loads(l)["epoch"] for l in frag.read_text().splitlines()
              if l.strip() and json.loads(l)["source_lane"] == "lpsn"}
    assert epochs == {"lpsn@LPSN-2026-06"}


def test_lpsn_since_narrows_the_gap_feeder_from_the_cli(tmp_path):
    # Same setup as the sharding test, but with a cutoff: only venerupis (2012) is new
    # enough, and the older gap is left out of the worklist.
    import json
    import lpsn as lpsn_mod
    from test_j2_feeders import LPSN_HEADER

    policy = tmp_path / "policy.tsv"
    policy.write_text("\n".join([
        POLICY_HEADER,
        _policy_row("Ab~GCF_1.1", "Arcobacter_butzleri", 100, "CP1", "GCF_1.1",
                    species_name="Arcobacter butzleri"),
    ]) + "\n")
    chrom = tmp_path / "chromosomes.tsv"
    chrom.write_text("scientificName\tnuccoreAcc\ttaxid\tparent\tsource\n"
                     "Arcobacter_butzleri\tCP1\t100\t1\tNCBI-REF\n")
    ani = tmp_path / "ani.txt"
    ani.write_text("".join([ANI_HEADER + "\n",
                            ani_line(genbank="GCA_1.1", refseq="GCF_1.1",
                                     species_taxid=100, best_taxid=100)]))
    refseq = tmp_path / "refseq.txt"
    refseq.write_text("".join([SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
                               summary_line(acc="GCF_1.1", species_taxid=100, paired="GCA_1.1")]))
    genbank = tmp_path / "genbank.txt"
    genbank.write_text("".join([SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
                                summary_line(acc="GCA_1.1", species_taxid=100,
                                             paired="GCF_1.1", source_db="genbank")]))

    lcache = tmp_path / "lpsn.json"
    rows = [["Arcobacter", "cryaerophilus", "", "correct name", "x", "1992-01-01"],
            ["Arcobacter", "venerupis", "", "correct name", "x", "2012-01-01"],
            ["Arcobacter", "retired", "", "synonym", "x", "2020-01-01"]]
    lc = lpsn_mod.LpsnCache(cache_path=str(lcache),
                            opener=lambda: [LPSN_HEADER + "\n"] + [",".join(r) + "\n" for r in rows])
    lc.build(release="LPSN-2026-06")
    lc.save()

    base = ["--policy", str(policy), "--chromosomes", str(chrom),
            "--cache", str(tmp_path / "cache.json"),
            "--ledger", str(tmp_path / "ledger.ndjson"), "--run-id", "2026-08"]
    tier1_metadata.main(base + ["--no-write", "--ani-file", str(ani),
                                "--refseq-file", str(refseq), "--genbank-file", str(genbank)])
    frag = tmp_path / "leads.0.ndjson"
    rc = tier1_metadata.main(base + ["--offline", "--with-lpsn", "--lpsn-cache", str(lcache),
                                     "--lpsn-since", "2000", "--emit-leads", str(frag)])
    assert rc == 0
    names = {json.loads(l)["scientific_name"] for l in frag.read_text().splitlines()
             if l.strip() and json.loads(l)["source_lane"] == "lpsn"}
    assert names == {"Arcobacter venerupis"}    # 1992 too old, 'retired' is a synonym


def test_reject_suppresses_lead(tmp_path):
    ledger = tmp_path / "ledger.ndjson"
    out = _run(tmp_path, _ani(tmp_path, "ani1.txt"), "2026-08", ledger)
    bar = _worklist(out)["Bar"]
    reviews = tmp_path / "reviews.tsv"
    with open(reviews, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["lead_key", "decision", "evidence_fingerprint",
                                           "epoch", "reviewer", "date", "note", "proposed_edit"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerow({"lead_key": bar["lead_key"], "decision": "reject",
                    "evidence_fingerprint": bar["evidence_fingerprint"], "epoch": bar["epoch"],
                    "reviewer": "t", "date": "2026-08-10", "note": "known", "proposed_edit": ""})
    ledger2 = tmp_path / "ledger2.ndjson"
    out2 = _run(tmp_path, _ani(tmp_path, "ani1.txt"), "2026-08", ledger2, reviews=reviews)
    assert "Bar" not in _statuses(out2)                # dismissed, off the worklist
