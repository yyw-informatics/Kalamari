# -*- coding: utf-8 -*-
"""J1b better-rep feed: the datasets-summary network scaffold (distiller, injected
fetcher, cache roundtrip) and the pure source-gating in leads.j1b_leads -- all offline."""

import json
import os

import coverage
import datasets_summary as ds
import leads
import tier1_metadata

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")


# --------------------------------------------------------------------------- #
#  Fixture builders (real datasets --as-json-lines report shape)               #
# --------------------------------------------------------------------------- #
def report_line(acc, taxid, category="", release="2023-05-01", tax_check="OK", name="Sp."):
    obj = {
        "accession": acc,
        "organism": {"tax_id": taxid, "organism_name": name},
        "assembly_info": {"refseq_category": category, "release_date": release},
        "average_nucleotide_identity": {"taxonomy_check_status": tax_check},
    }
    return json.dumps(obj) + "\n"


def _pick(unit="Foo~GCF_1.1", parent="GCF_1.1", lane="NCBI-REF", species_taxid=100,
          regime="", resolver=""):
    return coverage.Pick(unit_id=unit, scientific_name="Foo", species_taxid=species_taxid,
                         species_name="Foo bar", parent_assembly=parent, accession_set=["CP1"],
                         source_lane=lane, regime=regime, resolver=resolver)


# --------------------------------------------------------------------------- #
#  Distiller                                                                   #
# --------------------------------------------------------------------------- #
def test_distill_parses_fields():
    lines = [report_line("GCF_000002.1", 100, category="reference genome"),
             report_line("GCA_000003.1", 100, category="")]
    reports = ds.distill_datasets(lines)
    assert reports[0].base == "GCF_000002"
    assert reports[0].is_reference is True
    assert reports[0].source_db == "refseq"
    assert reports[0].release_date == "2023-05-01"
    assert reports[1].is_reference is False
    assert reports[1].source_db == "genbank"


def test_distill_tolerates_missing_keys():
    (r,) = ds.distill_datasets(['{"accession":"GCF_9.1"}\n'])
    assert r.base == "GCF_9" and r.species_taxid is None and r.refseq_category == ""


# --------------------------------------------------------------------------- #
#  Cache: injected fetcher, roundtrip, offline covers()                        #
# --------------------------------------------------------------------------- #
def _fetcher(mapping):
    return lambda taxid, after, ref_only: list(mapping.get(taxid, []))


def test_cache_build_roundtrip_and_covers(tmp_path):
    mapping = {100: [report_line("GCF_000002.1", 100, category="reference genome")],
               200: []}
    path = tmp_path / "ds.json"
    c = ds.DatasetsSummaryCache(cache_path=str(path), fetcher=_fetcher(mapping))
    c.build([100, 200], released_after="01/01/2026", version_epoch="ncbi-datasets-cli@16.0.0")
    c.save()
    assert c.candidates(100)[0].is_reference

    c2 = ds.DatasetsSummaryCache(cache_path=str(path))          # fresh load
    assert c2.covers([100, 200])
    assert not c2.covers([100, 999])
    assert c2.meta["datasets_version"] == "ncbi-datasets-cli@16.0.0"
    assert c2.candidates(100)[0].accession == "GCF_000002.1"


def test_cache_build_never_calls_default_subprocess(tmp_path):
    # Injecting a fetcher must not shell out to the real `datasets`.
    called = {"real": False}

    def boom(taxid, after, ref_only):
        called["real"] = True
        raise AssertionError("subprocess used")

    c = ds.DatasetsSummaryCache(cache_path=None, fetcher=boom)
    # build() DOES call the fetcher -- but our injected one, never run_datasets.
    try:
        c.build([100], "")
    except AssertionError:
        pass
    assert called["real"] is True   # our injected fetcher ran; the real path was never wired


# --------------------------------------------------------------------------- #
#  Query narrowing: --reference for J1b, the full query for J2                 #
# --------------------------------------------------------------------------- #
def test_the_datasets_command_asks_for_the_reference_genome_only_when_told_to():
    wide = ds.datasets_command(28901, "", "datasets")
    narrow = ds.datasets_command(28901, "", "datasets", reference_only=True)
    assert "--reference" not in wide            # an unfiltered S. enterica = 135 000 records
    assert "--reference" in narrow
    assert narrow[:5] == ["datasets", "summary", "genome", "taxon", "28901"]


def test_the_released_after_window_survives_the_reference_flag():
    cmd = ds.datasets_command(562, "01/01/2026", "datasets", reference_only=True)
    assert cmd[cmd.index("--released-after") + 1] == "01/01/2026"
    assert "--reference" in cmd


def test_build_asks_narrowly_only_for_the_taxids_named_reference_only(tmp_path):
    asked = {}

    def fetcher(taxid, after, ref_only):
        asked[taxid] = ref_only
        return []

    c = ds.DatasetsSummaryCache(cache_path=str(tmp_path / "ds.json"), fetcher=fetcher)
    c.build([100, 200, 300], released_after="", reference_only={100, 300})
    assert asked == {100: True, 200: False, 300: True}
    assert c.query_mode == {100: ds.QUERY_REFERENCE, 200: ds.QUERY_ALL,
                            300: ds.QUERY_REFERENCE}
    assert c.meta["n_reference_only"] == 2


def test_a_reference_only_cache_entry_cannot_serve_a_full_query_need(tmp_path):
    # The whole point of recording the query mode: J2's wishlist feeder falls back to
    # the newest genome when a species has no designated reference, so answering it
    # from a --reference query would silently report "no genome available yet".
    path = tmp_path / "ds.json"
    c = ds.DatasetsSummaryCache(cache_path=str(path), fetcher=_fetcher({}))
    c.build([100, 200], released_after="", reference_only={100})
    c.save()

    c2 = ds.DatasetsSummaryCache(cache_path=str(path))
    assert c2.covers([100, 200], reference_only={100})     # asked as built -> fine
    assert not c2.covers([100, 200])                        # 100 now needs the wide answer
    assert c2.uncovered([100, 200]) == [
        (100, "cached as reference-only, need the full query")]
    # ...and the wide entry happily serves a narrow need (a superset always does).
    assert c2.covers([200], reference_only={200})


def test_uncovered_separates_an_absent_taxid_from_a_too_narrow_one(tmp_path):
    path = tmp_path / "ds.json"
    c = ds.DatasetsSummaryCache(cache_path=str(path), fetcher=_fetcher({}))
    c.build([100], released_after="", reference_only={100})
    c.save()
    reasons = dict(ds.DatasetsSummaryCache(cache_path=str(path)).uncovered([100, 999]))
    assert reasons[999] == "absent"
    assert "reference-only" in reasons[100]


# --------------------------------------------------------------------------- #
#  --released-after windows: build() replaces, so a window is a NARROWING        #
# --------------------------------------------------------------------------- #
def test_an_unwindowed_cache_answers_a_windowed_question_but_never_the_reverse():
    # A window only ever removes records, so width flows one way: the cache that asked
    # NCBI for everything can answer "and what about since July?", while a cache that
    # only ever saw July cannot answer "everything".
    assert ds.window_covers(ds.WINDOW_ALL, "07/01/2026") is True
    assert ds.window_covers("07/01/2026", ds.WINDOW_ALL) is False
    assert ds.window_covers("01/01/2026", "07/01/2026") is True    # earlier start = wider
    assert ds.window_covers("07/01/2026", "01/01/2026") is False
    assert ds.window_covers("07/01/2026", "07/01/2026") is True
    assert ds.window_covers("2026-01-01", "07/01/2026") is True    # ISO understood too
    # An unparseable window is trusted only when identical: refusing costs a rebuild,
    # guessing costs an empty worklist that looks like good news.
    assert ds.window_covers("last tuesday", "last tuesday") is True
    assert ds.window_covers("last tuesday", "07/01/2026") is False


def test_a_cache_narrowed_by_a_window_cannot_answer_the_unwindowed_question(tmp_path):
    # THE regression that matters: build() REPLACES the cache, so a --released-after refresh
    # re-fetches every taxid from inside the window -- it is not a top-up. The pruned
    # entry is still present and its empty candidate list reads exactly like "NCBI has
    # nothing", so without this check the offline gate waves it through and the leads
    # simply vanish off the worklist.
    path = tmp_path / "ds.json"
    c = ds.DatasetsSummaryCache(cache_path=str(path), fetcher=_fetcher({}))
    c.build([100], released_after="07/01/2026", reference_only={100})
    c.save()

    c2 = ds.DatasetsSummaryCache(cache_path=str(path))
    assert c2.covers([100], reference_only={100}, released_after="07/01/2026")
    assert not c2.covers([100], reference_only={100})            # the unwindowed question
    assert c2.uncovered([100], reference_only={100}) == [
        (100, "cached only since 07/01/2026, need (all)")]
    # A window at least as wide as the cache's is still answerable.
    assert c2.covers([100], reference_only={100}, released_after="09/01/2026")


def test_build_never_windows_a_wishlist_taxid_even_when_a_window_is_given(tmp_path):
    # J2's wishlist asks a STANDING question ("has this species a genome at all?"), so a
    # window does not narrow the answer, it falsifies it: a species whose only genome
    # predates the window would read as "still nothing" and re-raise a solved lead.
    asked = {}

    def fetcher(taxid, after, ref_only):
        asked[taxid] = after
        return []

    c = ds.DatasetsSummaryCache(cache_path=str(tmp_path / "ds.json"), fetcher=fetcher)
    c.build([100, 200], released_after="07/01/2026", reference_only={100})
    assert asked == {100: "07/01/2026", 200: ds.WINDOW_ALL}
    assert c.window_of(100) == "07/01/2026" and c.window_of(200) == ds.WINDOW_ALL
    assert c.meta["n_windowed"] == 1
    # ...so the wishlist half of that same cache still answers its own feed.
    assert c.covers([200])


def test_a_cache_entry_with_no_window_of_its_own_inherits_the_cache_level_one(tmp_path):
    # Schema<=2 wrote only the cache-level window and every entry was fetched with it,
    # so inheriting is the truthful reading -- and it is what catches an older windowed
    # cache that predates the per-taxid field.
    path = tmp_path / "old.json"
    path.write_text(json.dumps({
        "schema": 2, "meta": {}, "released_after": "07/01/2026",
        "query_mode": {"100": ds.QUERY_REFERENCE},
        "by_species": {"100": []}}))
    c = ds.DatasetsSummaryCache(cache_path=str(path))
    assert c.window_of(100) == "07/01/2026"
    assert not c.covers([100], reference_only={100})


def test_a_cache_written_before_query_modes_existed_reads_as_fully_queried(tmp_path):
    # Schema-1 caches were built by the unfiltered query only, so 'all' is the honest
    # default -- an over-cautious 'reference' would refuse a cache that is in fact fine.
    path = tmp_path / "old.json"
    path.write_text(json.dumps({
        "schema": 1, "meta": {}, "released_after": "",
        "by_species": {"100": [{"accession": "GCF_1.1"}]}}))
    c = ds.DatasetsSummaryCache(cache_path=str(path))
    assert c.mode_of(100) == ds.QUERY_ALL
    assert c.covers([100])


# --------------------------------------------------------------------------- #
#  Pure source-gating (leads.j1b_leads)                                        #
# --------------------------------------------------------------------------- #
class _Cands:
    """Minimal candidate_source: taxid -> [GenomeReport]."""
    def __init__(self, mapping):
        self.by_species = mapping

    def candidates(self, taxid):
        return self.by_species.get(taxid, [])


def _ref(acc, taxid=100):
    (r,) = ds.distill_datasets([report_line(acc, taxid, category="reference genome")])
    return r


def _nonref(acc, taxid=100):
    (r,) = ds.distill_datasets([report_line(acc, taxid, category="")])
    return r


def test_j1b_flags_new_reference_for_ncbi_ref_lane():
    picks = [_pick(parent="GCF_1.1", lane="NCBI-REF", species_taxid=100)]
    cands = _Cands({100: [_ref("GCF_2.1")]})            # a different reference genome
    out = leads.j1b_leads(picks, cands)
    assert len(out) == 1
    assert out[0].signal == "J1b"
    assert out[0].candidate_accession == "GCF_2.1"
    assert out[0].lead_key == "Foo~GCF_1.1#GCF_2"


def test_j1b_ignores_our_own_assembly_as_reference():
    picks = [_pick(parent="GCF_1.1", lane="NCBI-REF", species_taxid=100)]
    cands = _Cands({100: [_ref("GCF_1.2")]})            # same base GCF_1 = our pick
    assert leads.j1b_leads(picks, cands) == []


def test_j1b_ignores_non_reference_candidates():
    picks = [_pick(parent="GCF_1.1", lane="NCBI-REF", species_taxid=100)]
    cands = _Cands({100: [_nonref("GCF_2.1")]})
    assert leads.j1b_leads(picks, cands) == []


def test_j1b_never_nags_sme_fda_nctc_lanes():
    cands = _Cands({100: [_ref("GCF_2.1")]})
    for lane in ("SME", "FDA-ARGOS", "NCTC3000"):
        picks = [_pick(parent="GCF_1.1", lane=lane, species_taxid=100)]
        assert leads.j1b_leads(picks, cands) == [], lane


def test_j1b_suppresses_confounded_units():
    picks = [_pick(parent="GCF_1.1", lane="NCBI-REF", species_taxid=100,
                   regime="B", resolver="ANI")]                 # confounded
    cands = _Cands({100: [_ref("GCF_2.1")]})
    assert leads.j1b_leads(picks, cands) == []


def test_j1b_skips_pick_without_species_taxid():
    picks = [_pick(parent="GCF_1.1", lane="NCBI-REF", species_taxid=None)]
    cands = _Cands({100: [_ref("GCF_2.1")]})
    assert leads.j1b_leads(picks, cands) == []


# --------------------------------------------------------------------------- #
#  The committed cache (built for real against NCBI; replayed offline in CI)    #
# --------------------------------------------------------------------------- #
class _RefreshArgs:
    """Just enough of the CLI namespace for ``tier1_metadata.datasets_taxids``."""

    def __init__(self):
        self.with_j1b = True
        self.with_j2 = True
        self.todo = os.path.join(REPO_ROOT, "src", "chromosomes-todo.tsv")


def _committed_needs():
    lanes = coverage.load_source_lanes(tier1_metadata.DEFAULT_CHROMOSOMES)
    picks = coverage.load_picks(tier1_metadata.DEFAULT_POLICY, lanes)
    return tier1_metadata.datasets_taxids(picks, _RefreshArgs())


def test_the_committed_datasets_cache_answers_both_feeds_at_the_width_they_need():
    # This is what makes the monthly CI job an offline replay: if the cache stops
    # covering the tracked species (a policy change adds one, say), the run must fail
    # here in a test rather than in a cron job nobody is watching.
    needed, reference_only = _committed_needs()
    cache = ds.DatasetsSummaryCache(cache_path=tier1_metadata.DEFAULT_DATASETS_CACHE)
    assert cache.uncovered(needed, reference_only) == []
    assert len(needed) == len(cache.by_species)         # no stale extra taxids either


def test_the_committed_datasets_cache_is_unwindowed_so_it_answers_every_run():
    # It was built with no --released-after, and the monthly offline run asks the
    # unwindowed question, so it must pass the gate exactly as committed. Being
    # unwindowed it also answers a run that DOES pass a window (a superset always does).
    needed, reference_only = _committed_needs()
    cache = ds.DatasetsSummaryCache(cache_path=tier1_metadata.DEFAULT_DATASETS_CACHE)
    assert cache.released_after == ds.WINDOW_ALL
    assert all(cache.window_of(t) == ds.WINDOW_ALL for t in needed)
    assert cache.uncovered(needed, reference_only, ds.WINDOW_ALL) == []
    assert cache.covers(needed, reference_only, "07/01/2026")


def test_the_committed_datasets_cache_records_its_resolved_tool_epoch():
    cache = ds.DatasetsSummaryCache(cache_path=tier1_metadata.DEFAULT_DATASETS_CACHE)
    assert cache.meta["datasets_version"].startswith("ncbi-datasets-cli@")
    assert "unset" not in cache.meta["datasets_version"]


def test_the_committed_datasets_cache_kept_the_big_species_narrow():
    # The reason the narrowing exists: unfiltered, Salmonella enterica alone returns
    # 135 000 records. Every J1b-tracked species must be stored reference-only, and a
    # reference-only entry is 0-2 records -- so a fat one means the flag went missing.
    needed, reference_only = _committed_needs()
    cache = ds.DatasetsSummaryCache(cache_path=tier1_metadata.DEFAULT_DATASETS_CACHE)
    fat = {t: len(cache.candidates(t)) for t in reference_only
           if len(cache.candidates(t)) > 5}
    assert fat == {}
    assert all(cache.mode_of(t) == ds.QUERY_REFERENCE for t in reference_only)
    assert all(cache.mode_of(t) == ds.QUERY_ALL
               for t in set(needed) - set(reference_only))
