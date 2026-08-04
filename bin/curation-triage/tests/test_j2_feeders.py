# -*- coding: utf-8 -*-
"""J2 coverage feeders: chromosomes-todo.tsv loading, the LPSN network scaffold, and
the pure j2_leads matching (wishlist + LPSN gaps) -- all offline."""

import os

import lpsn
import leads
from leads import TodoEntry
from test_datasets_summary import _Cands, _ref, _nonref, report_line
import datasets_summary as ds

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))       # tests -> curation-triage -> bin -> repo
COMMITTED_LPSN_CACHE = os.path.join(REPO_ROOT, "src", "curation-triage", "lpsn_cache.json")


# --------------------------------------------------------------------------- #
#  load_todo: 4 cols, no header, spaces in names                               #
# --------------------------------------------------------------------------- #
def test_load_todo(tmp_path):
    p = tmp_path / "todo.tsv"
    p.write_text("Arcobacter butzleri\tXXXXXX\t28197\t28196\n"
                 "Caulobacter\tCP000927\t366602\t2648921\n")
    entries = leads.load_todo(str(p))
    assert entries[0].scientific_name == "Arcobacter butzleri"
    assert entries[0].is_wishlist is True and entries[0].taxid == 28197
    assert entries[1].is_wishlist is False and entries[1].nuccore_acc == "CP000927"


def test_normalize_species_and_genus():
    assert leads.normalize_species("Arcobacter_butzleri") == "arcobacter butzleri"
    assert leads.species_genus("Escherichia_coli") == "escherichia"


def test_load_todo_blank_accession_is_wishlist(tmp_path):
    # A blank accession cell must behave like the wishlist sentinel, not a bogus
    # "ready to add" lead with an empty accession.
    p = tmp_path / "todo.tsv"
    p.write_text("Arcobacter butzleri\t\t28197\t28196\n")
    (e,) = leads.load_todo(str(p))
    assert e.is_wishlist is True
    out = leads.j2_leads([e], _Cands({}), [], set(), set(), max_lpsn=None)
    assert out == []                          # no candidate yet -> no lead (not a blank one)


# --------------------------------------------------------------------------- #
#  LPSN distiller + cache                                                      #
# --------------------------------------------------------------------------- #
# The header the synthetic-row tests write: a SUBSET of the real LPSN header, using the
# real column names.  csv.DictReader keys by name, so a subset in a different order is a
# legitimate fixture -- ``test_the_real_lpsn_header_parses`` below covers the full thing.
# The sixth column is ``reference``, not a date: LPSN ships NO date column, and a year in
# the citation is the documented fallback when ``authors`` carries none.
LPSN_HEADER = "genus_name,sp_epithet,subsp_epithet,status,authors,reference"

# Five rows copied VERBATIM out of the real complete-list export (release 2026-08-04),
# quoting and all, so the header/parse tests are checked against what LPSN really ships
# rather than against our idea of it.  Long lines on purpose: reflowing them would stop
# them being verbatim.
LPSN_REAL_HEADER = ("genus_name,sp_epithet,subsp_epithet,reference,status,authors,"
                    "address,risk_grp,nomenclatural_type,record_no,record_lnk")
LPSN_REAL_ROWS = [
    # plain valid publication -> 2018
    'Abditibacterium,utsteinense,,"Oren A, Garrity GM. Validation list no. 184. List of new names and new combinations previously effectively, but not validly, published. Int J Syst Evol Microbiol 2018; 68:3379-3393.","VL; sp. nov.; validly published under the ICNP; correct name","Tahon et al. 2018",https://lpsn.dsmz.de/species/abditibacterium-utsteinense,,"DSM 105287; LMG 29911; R-68213",797965,',
    # Approved Lists -> 1980 (and a 'correct name, recommended for medical use' status)
    'Listeria,monocytogenes,,"Skerman VBD, McGowan V, Sneath PHA. Approved lists of bacterial names. Int J Syst Bacteriol 1980; 30:225-420.","VL; comb. nov.; validly published under the ICNP; correct name, recommended for medical use","(Murray et al. 1926) Pirie 1940 (Approved Lists 1980)",https://lpsn.dsmz.de/species/listeria-monocytogenes,2,"ATCC 15313; CCUG 15526; CIP 82.110; DSM 20600; NBIMCC 8669; NCTC 10357; SLCC 53",777590,',
    # new combination -> the COMBINATION year 1995, not the basionym's 1989
    'Abiotrophia,adiacens,,"Kawamura Y, Hou XG, Sultana F, Liu S, Yamamoto H, Ezaki T. Transfer of Streptococcus adjacens and Streptococcus defectivus to Abiotrophia gen. nov. as Abiotrophia adiacens comb. nov. and Abiotrophia defectiva comb. nov., respectively. Int J Syst Bacteriol 1995; 45:798-803.","VP; comb. nov.; validly published under the ICNP; synonym, not recommended for medical use","(Bouvet et al. 1989) Kawamura et al. 1995",https://lpsn.dsmz.de/species/abiotrophia-adiacens,2,"ATCC 49175; CCUG 27809; CIP 103.243; DSM 9848; GaD; LMG 14496; NCTC 13000",772466,776611',
    # emendation -> 2011, NOT the 2026 emendation
    'Acinetobacter,brisouii,,"Euzeby JP. Validation list no. 140. List of new names and new combinations previously effectively, but not validly, published. Int J Syst Evol Microbiol 2011; 61:1499-1501.","VL; sp. nov.; validly published under the ICNP; correct name, recommended for medical use","Anandham et al. 2011 emend. Radolfova-Krizova et al. 2026",https://lpsn.dsmz.de/species/acinetobacter-brisouii,2,"5YN5-8; ANC 4119; CCUG 61636; CIP 110357T; DSM 18516; KACC 11602",785346,',
    # 'inappropriate correction; misspelling' -- the status trap
    'Actinospongiicola,halichondriae,,"Oren A, Goker M. Validation list no. 227. List of new names and new combinations previously effectively, but not validly, published. Int J Syst Evol Microbiol 2026; 76:6989.","VL; sp. nov.; validly published under the ICNP, inappropriate correction; misspelling","Huang et al. 2026",https://lpsn.dsmz.de/species/actinospongiicola-halichondriae,,"DSM 114536; Hal317; LMG 32795",70412,70465',
]


def _lpsn_csv(*rows):
    return [LPSN_HEADER + "\n"] + [",".join(r) + "\n" for r in rows]


def _real_csv():
    return [LPSN_REAL_HEADER + "\n"] + [r + "\n" for r in LPSN_REAL_ROWS]


# --------------------------------------------------------------------------- #
#  The real header + the derived valid-publication year                        #
# --------------------------------------------------------------------------- #
def test_the_real_lpsn_header_parses_into_species_records():
    # The pinned LPSN_COLS were read off the real 2026-08-04 export; this is the
    # regression test that says so, using rows copied verbatim from it.
    recs = {r.species: r for r in lpsn.distill_lpsn(_real_csv())}
    assert set(recs) == {"Abditibacterium utsteinense", "Listeria monocytogenes",
                         "Abiotrophia adiacens", "Acinetobacter brisouii",
                         "Actinospongiicola halichondriae"}
    assert recs["Abditibacterium utsteinense"].authority == "Tahon et al. 2018"


def test_the_real_export_has_no_date_column_at_all():
    # Every 'valid_publication_date'-style name we once guessed at is absent, which is
    # exactly why the year is derived from the authority instead of read from a column.
    assert "date" not in lpsn.LPSN_COLS
    assert lpsn.column_candidates("authority") == ("authors",)
    for guess in ("valid_publication_date", "date_of_valid_publication",
                  "valid_publication_year", "valid_publication"):
        assert guess not in LPSN_REAL_HEADER


def test_a_plain_authority_year_is_the_valid_publication_year():
    assert lpsn.valid_publication_year("Tahon et al. 2018") == "2018"


def test_a_new_combination_dates_from_the_combination_not_the_basionym():
    # The basionym year comes FIRST and in parentheses; the current name dates from the
    # combination, so the LAST year wins.
    assert lpsn.valid_publication_year("(Bouvet et al. 1989) Kawamura et al. 1995") == "1995"


def test_an_emendation_does_not_move_the_valid_publication_year():
    # 2026 amended the description of a name published in 2011; it did not publish it.
    assert lpsn.valid_publication_year(
        "Anandham et al. 2011 emend. Radolfova-Krizova et al. 2026") == "2011"
    assert lpsn.valid_publication_year(
        "Jiang et al. 2016 emend. Madhaiyan et al. 2020 emend. Dobritsa 2020") == "2016"
    # 'corrig.' (a corrected spelling) is the same class of later annotation.
    assert lpsn.valid_publication_year("Smith 1999 corrig. Jones 2005") == "1999"


def test_an_approved_lists_name_dates_from_the_approved_lists():
    # Pre-1980 names gained standing under the code only with the Approved Lists, and
    # LPSN writes that year last -- which is why 'last year wins' is the right rule.
    assert lpsn.valid_publication_year(
        "(Murray et al. 1926) Pirie 1940 (Approved Lists 1980)") == "1980"


def test_the_year_falls_back_to_the_reference_when_the_authority_has_none():
    assert lpsn.valid_publication_year(
        "", "Oren A. Validation list no. 184. Int J Syst Evol Microbiol 2018; 68:3379-3393.") \
        == "2018"
    # ...and a journal volume or page range can never be mistaken for a year.
    assert lpsn.valid_publication_year("", "Int J Syst Bacteriol 30:225-420") == ""


def test_the_derived_year_lands_in_valid_date_for_every_real_row():
    recs = {r.species: r.valid_date for r in lpsn.distill_lpsn(_real_csv())}
    assert recs == {"Abditibacterium utsteinense": "2018",
                    "Listeria monocytogenes": "1980",
                    "Abiotrophia adiacens": "1995",
                    "Acinetobacter brisouii": "2011",
                    "Actinospongiicola halichondriae": "2026"}


# --------------------------------------------------------------------------- #
#  Taxonomic standing (the 'correct name' test)                                #
# --------------------------------------------------------------------------- #
def test_an_inappropriate_correction_is_a_misspelling_not_a_correct_name():
    # The trap: 12 rows in the 2026-08-04 export read '... validly published under the
    # ICNP, inappropriate correction; misspelling'. A substring test for "correct"
    # matches them, so it would hand an SME 10 misspelled species to add.
    (rec,) = [r for r in lpsn.distill_lpsn(_real_csv())
              if r.species == "Actinospongiicola halichondriae"]
    assert "correct" in rec.status.lower()          # the substring really is in there
    assert rec.is_correct_name is False             # ...and it must not count


def test_the_standing_is_the_last_semicolon_clause_minus_the_medical_use_rider():
    standing = lpsn.taxonomic_standing
    assert standing("VL; sp. nov.; validly published under the ICNP; correct name") \
        == "correct name"
    assert standing("VL; comb. nov.; validly published under the ICNP; "
                    "correct name, recommended for medical use") == "correct name"
    assert standing("VP; comb. nov.; validly published under the ICNP; "
                    "synonym, not recommended for medical use") == "synonym"
    assert standing("VL; sp. nov.; validly published under the ICNP, "
                    "inappropriate correction; misspelling") == "misspelling"
    assert standing("") == ""


def test_a_medical_use_rider_never_hides_a_correct_name():
    # 11 152 of the correct names carry ', recommended for medical use' -- dropping them
    # would blind the feeder to most of the pathogens it exists for.
    rec = lpsn.LpsnRecord(status="VL; comb. nov.; validly published under the ICNP; "
                                 "correct name, recommended for medical use")
    assert rec.is_correct_name is True


def test_distill_lpsn_correct_names_and_binomial():
    lines = _lpsn_csv(
        ["Arcobacter", "butzleri", "", "correct name", "Kiehlbauch 1991", "1991-01-01"],
        ["Arcobacter", "skirrowii", "", "synonym", "x", "1992-01-01"],
        ["Escherichia", "coli", "", "correct name", "Castellani 1919", "1919-01-01"],
    )
    recs = lpsn.distill_lpsn(lines)
    corrects = [r for r in recs if r.is_correct_name]
    assert {r.species for r in corrects} == {"Arcobacter butzleri", "Escherichia coli"}


def test_distill_lpsn_correct_name_wins_over_synonym_regardless_of_order():
    # A synonym/subspecies row for a binomial listed BEFORE the correct-name species row
    # must not shadow it (else a real coverage gap is silently dropped).
    lines = _lpsn_csv(
        ["Arcobacter", "butzleri", "thereius", "synonym", "x", "2011-01-01"],
        ["Arcobacter", "butzleri", "", "correct name", "y", "1991-01-01"],
    )
    corrects = [r for r in lpsn.distill_lpsn(lines) if r.is_correct_name]
    assert [r.species for r in corrects] == ["Arcobacter butzleri"]


def test_distill_lpsn_schema_change_raises():
    bad = ["genus_name,WRONG,status\n", "Arcobacter,butzleri,correct name\n"]
    try:
        lpsn.distill_lpsn(bad)
        assert False, "expected schema-change ValueError"
    except ValueError as e:
        assert "missing expected column" in str(e)


def test_the_distiller_keeps_synonyms_so_the_cache_stays_a_faithful_copy():
    # Gating is policy and lives in leads.j2_leads; the cache must not pre-filter, or
    # changing our mind about synonyms would need an account-gated re-download.
    lines = _lpsn_csv(
        ["Arcobacter", "butzleri", "", "correct name", "y", "1991-01-01"],
        ["Arcobacter", "skirrowii", "", "synonym", "x", "1992-01-01"],
    )
    recs = lpsn.distill_lpsn(lines)
    assert {r.species: r.status for r in recs} == {
        "Arcobacter butzleri": "correct name", "Arcobacter skirrowii": "synonym"}


def test_the_distiller_derives_the_valid_publication_year_into_its_own_field():
    # This test used to assert a full date read out of a guessed 'valid_publication_date'
    # column. The real export has no such column, so the expectation is now the YEAR
    # derived from the authority -- the only place LPSN records one.
    (rec,) = lpsn.distill_lpsn(_lpsn_csv(
        ["Arcobacter", "butzleri", "", "correct name", "Kiehlbauch et al. 1991", ""]))
    assert rec.valid_date == "1991"
    assert rec.authority == "Kiehlbauch et al. 1991"   # kept verbatim as SME context


def test_a_missing_authority_column_fails_loudly_instead_of_blanking_the_year():
    # The year drives the --lpsn-since gate and now comes from 'authors'; a silently
    # blank one would drop every gap the gate was meant to surface, so its absence must
    # be a schema error. (Before the real header was read, this guarded the guessed
    # date column -- same contract, now pointed at the column that really carries it.)
    no_authors = ["genus_name,sp_epithet,subsp_epithet,status,reference\n",
                  "Arcobacter,butzleri,,correct name,x\n"]
    try:
        lpsn.distill_lpsn(no_authors)
        assert False, "expected a missing-authority-column ValueError"
    except ValueError as e:
        assert "authority" in str(e) and "missing expected column" in str(e)


def test_a_renamed_column_is_pinnable_through_a_list_of_candidates(monkeypatch):
    # LPSN can still rename a column; re-pinning it must stay a one-line edit, so
    # LPSN_COLS keeps accepting a tuple of spellings.
    monkeypatch.setitem(lpsn.LPSN_COLS, "authority", ("authors", "authority_string"))
    lines = ["genus_name,sp_epithet,subsp_epithet,status,authority_string\n",
             "Arcobacter,butzleri,,correct name,Kiehlbauch et al. 1991\n"]
    (rec,) = lpsn.distill_lpsn(lines)
    assert rec.valid_date == "1991"


def test_pinning_a_column_to_a_bare_string_works_like_the_tuple_form(monkeypatch):
    # ...and "one-line change" must mean a bare string too, not only a 1-tuple.
    monkeypatch.setitem(lpsn.LPSN_COLS, "authority", "pinned_authors")
    lines = ["genus_name,sp_epithet,subsp_epithet,status,pinned_authors\n",
             "Arcobacter,butzleri,,correct name,Kiehlbauch et al. 1991\n"]
    (rec,) = lpsn.distill_lpsn(lines)
    assert rec.valid_date == "1991"


def test_lpsn_cache_roundtrip(tmp_path):
    path = tmp_path / "lpsn.json"
    c = lpsn.LpsnCache(
        cache_path=str(path),
        opener=lambda: _lpsn_csv(
            ["Arcobacter", "butzleri", "", "correct name", "Kiehlbauch et al. 1991", ""]))
    c.build(release="LPSN-2026-06-01")
    c.save()
    c2 = lpsn.LpsnCache(cache_path=str(path))
    assert c2.release == "LPSN-2026-06-01"
    assert [r.species for r in c2.in_genera({"Arcobacter"})] == ["Arcobacter butzleri"]
    assert c2.in_genera({"Escherichia"}) == []
    assert c2.records[0].valid_date == "1991"      # the derived year survives the roundtrip


def test_lpsn_cache_never_calls_default_opener():
    called = {"url": False}

    def boom():
        called["url"] = True
        raise AssertionError("network used")

    c = lpsn.LpsnCache(
        cache_path=None,
        opener=lambda: _lpsn_csv(["G", "s", "", "correct name", "", "2000-01-01"]))
    c.build(release="r")
    assert called["url"] is False


def test_the_default_opener_explains_that_lpsn_needs_an_account(tmp_path):
    # There is no unattended download: the public URL is a landing page and the
    # complete-list export is account-gated. Say that, rather than fetching HTML and
    # failing later with a confusing "missing expected column".
    c = lpsn.LpsnCache(cache_path=str(tmp_path / "nope.json"))
    try:
        c.build(release="r")
        assert False, "expected the default opener to refuse"
    except RuntimeError as e:
        assert "account" in str(e) and "--lpsn-csv" in str(e)


def test_an_operator_downloaded_csv_is_the_working_refresh_path(tmp_path):
    csv_path = tmp_path / "lpsn.csv"
    csv_path.write_text("".join(_lpsn_csv(
        ["Arcobacter", "butzleri", "", "correct name", "Kiehlbauch et al. 1991", ""])))
    c = lpsn.LpsnCache(cache_path=str(tmp_path / "lpsn.json"),
                       opener=lambda: lpsn.open_lpsn_source(str(csv_path)))
    c.build(release="LPSN-2026-06-01")
    assert [r.species for r in c.records] == ["Arcobacter butzleri"]
    assert c.meta["n_correct_names"] == 1 and c.meta["n_with_valid_date"] == 1


def test_a_local_csv_is_offline_safe_but_a_download_link_is_not():
    # tier1_metadata lets --offline rebuild the cache from a downloaded FILE (no network)
    # while still refusing a short-lived per-account download LINK.
    assert lpsn.is_remote_source("https://lpsn.dsmz.de/downloads/abc") is True
    assert lpsn.is_remote_source("HTTP://example.org/x.csv") is True
    assert lpsn.is_remote_source("lpsn_gss_2026-08-04.csv") is False
    assert lpsn.is_remote_source("") is False


def test_the_refresh_help_names_the_api_as_the_future_unattended_path():
    # The account now exists and the CSV path works; what is still missing is an
    # unattended refresh, and the next developer should not have to guess where it goes.
    assert lpsn.LPSN_API_ROOT == "https://api.lpsn.dsmz.de/"
    assert lpsn.LPSN_API_ROOT in lpsn.LPSN_REFRESH_HELP
    assert "--lpsn-csv" in lpsn.LPSN_REFRESH_HELP and "sign in" in lpsn.LPSN_REFRESH_HELP


# --------------------------------------------------------------------------- #
#  The COMMITTED cache (built from the real 2026-08-04 download)               #
# --------------------------------------------------------------------------- #
def test_the_committed_lpsn_cache_loads_with_the_shape_the_feeder_expects():
    # The whole point of committing it: --with-lpsn now runs offline, in CI, with no
    # account and no download. If this cache ever stops loading, that stops being true.
    c = lpsn.LpsnCache(cache_path=COMMITTED_LPSN_CACHE)
    assert c.release == "2026-08-04"
    assert c.meta["lpsn_release"] == "2026-08-04"
    assert len(c.records) == c.meta["n_records"] > 28000
    assert c.meta["n_correct_names"] > 22000
    # Every record must carry the fields j2_leads reads, and a derived year for each.
    assert all(r.genus and r.species and r.status and r.valid_date for r in c.records)


def test_the_committed_cache_answers_the_real_gap_question_for_a_tracked_genus():
    c = lpsn.LpsnCache(cache_path=COMMITTED_LPSN_CACHE)
    yersinia = {r.species: r for r in c.in_genera({"Yersinia"})}
    assert "Yersinia fenwicki" in yersinia                  # named 2026, Kalamari lacks it
    assert yersinia["Yersinia fenwicki"].valid_date == "2026"
    assert yersinia["Yersinia pestis"].valid_date == "1980"  # an Approved-Lists name
    # in_genera returns correct names only, so a retired name can never become a lead.
    assert all(r.is_correct_name for r in yersinia.values())


def test_the_committed_cache_kept_the_synonyms_it_does_not_report():
    # Faithful copy: the gating is policy and lives in leads.j2_leads, so the cache must
    # still contain the names it will never surface.
    c = lpsn.LpsnCache(cache_path=COMMITTED_LPSN_CACHE)
    standings = {lpsn.taxonomic_standing(r.status) for r in c.records}
    assert {"correct name", "synonym", "misspelling"} <= standings
    assert sum(1 for r in c.records if not r.is_correct_name) > 5000


# --------------------------------------------------------------------------- #
#  j2_leads: wishlist feeder                                                   #
# --------------------------------------------------------------------------- #
def _todo(name, acc="XXXXXX", taxid=100, parent="1"):
    return TodoEntry(scientific_name=name, nuccore_acc=acc, taxid=taxid, parent=parent)


def test_j2_wishlist_surfaces_when_genome_available():
    todo = [_todo("Arcobacter butzleri", acc="XXXXXX", taxid=28197)]
    cands = _Cands({28197: [_ref("GCF_500.1", taxid=28197)]})
    out = leads.j2_leads(todo, cands, lpsn_records=[], tracked_genera=set(),
                         have_species=set(), max_lpsn=None)
    assert len(out) == 1 and out[0].signal == "J2"
    assert out[0].candidate_accession == "GCF_500.1"
    assert out[0].lead_key == "todo:28197#GCF_500.1"


def test_j2_wishlist_silent_when_no_genome_yet():
    todo = [_todo("Arcobacter cloacae", acc="XXXXXX", taxid=1054034)]
    cands = _Cands({})                              # nothing available yet
    out = leads.j2_leads(todo, cands, [], set(), set(), max_lpsn=None)
    assert out == []


def test_j2_wishlist_ready_to_add_when_accession_known():
    todo = [_todo("Caulobacter", acc="CP000927", taxid=366602)]
    out = leads.j2_leads(todo, _Cands({}), [], set(), set(), max_lpsn=None)
    assert len(out) == 1
    assert out[0].candidate_accession == "CP000927"
    assert "ready to add" in out[0].summary


# --------------------------------------------------------------------------- #
#  j2_leads: LPSN gap feeder                                                   #
# --------------------------------------------------------------------------- #
def _rec(genus, species, status="correct name", valid_date="2000-01-01"):
    return lpsn.LpsnRecord(genus=genus, species=species, status=status,
                           valid_date=valid_date)


def test_j2_lpsn_gap_in_tracked_genus_not_held():
    recs = [_rec("Arcobacter", "Arcobacter butzleri"),   # tracked genus, not held -> gap
            _rec("Arcobacter", "Arcobacter cryaerophilus"),
            _rec("Streptomyces", "Streptomyces coelicolor")]  # untracked genus -> ignored
    out = leads.j2_leads(
        todo_entries=[], candidate_source=_Cands({}), lpsn_records=recs,
        tracked_genera={"arcobacter"},
        have_species={"arcobacter butzleri"},          # we already have this one
        max_lpsn=25)
    names = {l.scientific_name for l in out}
    assert names == {"Arcobacter cryaerophilus"}       # butzleri held, coelicolor untracked


def test_j2_lpsn_cap_and_dropped_logged():
    recs = [_rec("Arcobacter", "Arcobacter sp%02d" % i) for i in range(10)]
    logs = []
    out = leads.j2_leads([], _Cands({}), recs, {"arcobacter"}, set(),
                         max_lpsn=3, log=logs.append)
    lpsn_leads = [l for l in out if l.source_lane == "lpsn"]
    assert len(lpsn_leads) == 3
    assert any("dropped 7" in m for m in logs)


def test_the_cap_keeps_the_newest_gaps_not_the_alphabetically_first_ones():
    # The real catalog has ~4000 LPSN gaps. Sorted by name, a cap of 25 shows 25
    # Acinetobacter and hides every genus after 'A' -- useless for a feeder whose whole
    # point is "newly published". Newest first makes the cap mean what it says.
    recs = [_rec("Arcobacter", "Arcobacter aardvarki", valid_date="1985"),
            _rec("Arcobacter", "Arcobacter beta", valid_date="2001"),
            _rec("Arcobacter", "Arcobacter zulu", valid_date="2026")]
    out = leads.j2_leads([], _Cands({}), recs, {"arcobacter"}, set(), max_lpsn=2)
    assert [l.scientific_name for l in out] == ["Arcobacter zulu", "Arcobacter beta"]


def test_same_year_gaps_stay_in_a_deterministic_alphabetical_order():
    # Ties must not depend on the order LPSN happened to write the rows, or the monthly
    # worklist would churn for no reason.
    recs = [_rec("Arcobacter", "Arcobacter charlie", valid_date="2026"),
            _rec("Arcobacter", "Arcobacter alpha", valid_date="2026"),
            _rec("Arcobacter", "Arcobacter bravo", valid_date="2026")]
    out = leads.j2_leads([], _Cands({}), recs, {"arcobacter"}, set(), max_lpsn=25)
    assert [l.scientific_name for l in out] == [
        "Arcobacter alpha", "Arcobacter bravo", "Arcobacter charlie"]


def test_an_undated_gap_sorts_last_so_the_cap_never_spends_a_slot_on_it():
    recs = [_rec("Arcobacter", "Arcobacter aaaa", valid_date=""),
            _rec("Arcobacter", "Arcobacter zzzz", valid_date="1985")]
    out = leads.j2_leads([], _Cands({}), recs, {"arcobacter"}, set(), max_lpsn=1)
    assert [l.scientific_name for l in out] == ["Arcobacter zzzz"]


def test_j2_lpsn_disabled_when_max_none():
    recs = [_rec("Arcobacter", "Arcobacter butzleri")]
    out = leads.j2_leads([], _Cands({}), recs, {"arcobacter"}, set(), max_lpsn=None)
    assert out == []


def test_a_synonym_name_is_never_reported_as_a_coverage_gap():
    # A synonym is a name LPSN itself has retired; telling an SME to add it is exactly
    # the kind of false lead that trains reviewers to stop reading the worklist.
    recs = [_rec("Arcobacter", "Arcobacter butzleri", status="synonym"),
            _rec("Arcobacter", "Arcobacter venerupis", status="correct name")]
    out = leads.j2_leads([], _Cands({}), recs, {"arcobacter"}, set(), max_lpsn=25)
    assert [l.scientific_name for l in out] == ["Arcobacter venerupis"]


def test_the_since_cutoff_keeps_only_newly_published_names():
    recs = [_rec("Arcobacter", "Arcobacter old", valid_date="1991-06-15"),
            _rec("Arcobacter", "Arcobacter new", valid_date="2024-03-01")]
    out = leads.j2_leads([], _Cands({}), recs, {"arcobacter"}, set(),
                         max_lpsn=25, since="2024-01-01")
    assert [l.scientific_name for l in out] == ["Arcobacter new"]


def test_a_year_only_cutoff_matches_a_full_date():
    # LPSN's date granularity is unverified, so the cutoff compares left-anchored and
    # truncates to the shorter side: '2024' must accept '2024-03-01'.
    assert leads.valid_since("2024-03-01", "2024") is True
    assert leads.valid_since("2023-12-31", "2024") is False
    assert leads.valid_since("2024", "2024-03-01") is True     # can't disprove -> keep


def test_an_undated_name_is_dropped_once_a_cutoff_is_set():
    recs = [_rec("Arcobacter", "Arcobacter undated", valid_date="")]
    assert leads.j2_leads([], _Cands({}), recs, {"arcobacter"}, set(),
                          max_lpsn=25, since="2024") == []
    # ...but with no cutoff it is still a gap (we are not asking about dates at all).
    out = leads.j2_leads([], _Cands({}), recs, {"arcobacter"}, set(), max_lpsn=25)
    assert [l.scientific_name for l in out] == ["Arcobacter undated"]


def test_the_lpsn_lead_carries_its_valid_date_as_evidence():
    recs = [_rec("Arcobacter", "Arcobacter venerupis", valid_date="2012-05-01")]
    (lead,) = leads.j2_leads([], _Cands({}), recs, {"arcobacter"}, set(), max_lpsn=25)
    assert lead.evidence["valid_date"] == "2012-05-01"
    assert lead.deposit_date == "2012-05-01"


def test_the_lpsn_epoch_always_names_the_tool():
    # It used to stamp the bare release tag, so a run without --lpsn-release produced
    # the naked literal 'unset' as an epoch.
    assert leads.lpsn_epoch("2026-06-01") == "lpsn@2026-06-01"
    assert leads.lpsn_epoch("") == "lpsn@unset"
    assert leads.lpsn_epoch("unset") == leads.LPSN_EPOCH


def test_j2_todo_and_lpsn_carry_separate_epochs():
    # A datasets-cli bump must re-baseline todo leads, an LPSN bump the LPSN leads --
    # so the two feeders must stamp DISTINCT tooling epochs.
    todo = [_todo("Arcobacter butzleri", acc="XXXXXX", taxid=28197)]
    cands = _Cands({28197: [_ref("GCF_1.1", taxid=28197)]})
    recs = [_rec("Arcobacter", "Arcobacter cryaerophilus")]
    out = leads.j2_leads(todo, cands, recs, {"arcobacter"}, set(),
                         todo_epoch="ncbi-datasets-cli@9", lpsn_epoch="lpsn@9", max_lpsn=25)
    by_lane = {l.source_lane: l for l in out}
    assert by_lane["todo"].epoch == "ncbi-datasets-cli@9"
    assert by_lane["lpsn"].epoch == "lpsn@9"
