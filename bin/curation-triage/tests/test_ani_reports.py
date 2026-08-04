# -*- coding: utf-8 -*-
"""The ASSEMBLY_REPORTS network layer: streaming distillers, header validation,
accession/version helpers, and cache roundtrip -- all offline (no wire)."""

import pytest

import ani_reports
from ani_reports import (
    AssemblyReportsCache,
    accession_base,
    split_accession_version,
    distill_ani_report,
    distill_assembly_summary,
    ANI_COLS,
    SUMMARY_COLS,
)


# --------------------------------------------------------------------------- #
#  Raw-line builders (real column layout; 'na' everywhere we don't set)        #
# --------------------------------------------------------------------------- #
ANI_HEADER = ("# genbank-accession\trefseq-accession\ttaxid\tspecies-taxid\t"
              "organism-name\tspecies-name\tassembly-name\tassembly-type-category\t"
              "excluded-from-refseq\tdeclared-type-assembly\tdeclared-type-organism-name\t"
              "declared-type-category\tdeclared-type-ANI\tdeclared-type-qcoverage\t"
              "declared-type-scoverage\tbest-match-type-assembly\tbest-match-species-taxid\t"
              "best-match-species-name\tbest-match-type-category\tbest-match-type-ANI\t"
              "best-match-type-qcoverage\tbest-match-type-scoverage\tbest-match-status\t"
              "comment\ttaxonomy-check-status")

SUMMARY_HEADER_1 = "##  See ftp://.../README_assembly_summary.txt for a description."
SUMMARY_HEADER_2 = ("#assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\t"
                    "taxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\t"
                    "version_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\t"
                    "asm_name\tasm_submitter\tgbrs_paired_asm\tpaired_asm_comp\tftp_path\t"
                    "excluded_from_refseq\trelation_to_type_material\tasm_not_live_date\t"
                    "assembly_type\tgroup\tgenome_size")


def ani_line(**kw):
    f = ["na"] * 25
    f[ANI_COLS["genbank_accession"]] = kw.get("genbank", "na")
    f[ANI_COLS["refseq_accession"]] = kw.get("refseq", "na")
    f[ANI_COLS["species_taxid"]] = str(kw.get("species_taxid", "na"))
    f[ANI_COLS["species_name"]] = kw.get("species_name", "na")
    f[ANI_COLS["best_match_species_taxid"]] = str(kw.get("best_taxid", "na"))
    f[ANI_COLS["best_match_species_name"]] = kw.get("best_name", "na")
    f[ANI_COLS["best_match_status"]] = kw.get("best_status", "species-match")
    f[ANI_COLS["taxonomy_check_status"]] = kw.get("tax_check", "OK")
    f[ANI_COLS["excluded_from_refseq"]] = kw.get("excluded", "na")
    return "\t".join(f) + "\n"


def summary_line(**kw):
    f = ["na"] * 26
    f[SUMMARY_COLS["assembly_accession"]] = kw.get("acc", "na")
    f[SUMMARY_COLS["refseq_category"]] = kw.get("category", "na")
    f[SUMMARY_COLS["taxid"]] = str(kw.get("taxid", "na"))
    f[SUMMARY_COLS["species_taxid"]] = str(kw.get("species_taxid", "na"))
    f[SUMMARY_COLS["version_status"]] = kw.get("version_status", "latest")
    f[SUMMARY_COLS["seq_rel_date"]] = kw.get("seq_rel_date", "2018/01/01")
    f[SUMMARY_COLS["gbrs_paired_asm"]] = kw.get("paired", "na")
    f[SUMMARY_COLS["excluded_from_refseq"]] = kw.get("excluded", "na")
    f[SUMMARY_COLS["relation_to_type_material"]] = kw.get("relation", "na")
    return "\t".join(f) + "\n"


# --------------------------------------------------------------------------- #
#  Accession helpers                                                           #
# --------------------------------------------------------------------------- #
def test_split_accession_version():
    assert split_accession_version("GCF_000008685.2") == ("GCF_000008685", 2)
    assert split_accession_version("GCF_1") == ("GCF_1", None)
    assert accession_base("GCA_000003135.1") == "GCA_000003135"


# --------------------------------------------------------------------------- #
#  ANI-report distiller                                                        #
# --------------------------------------------------------------------------- #
def test_distill_ani_matches_via_either_accession():
    lines = [
        ANI_HEADER + "\n",
        ani_line(genbank="GCA_1.1", refseq="GCF_1.1", species_taxid=100, best_taxid=100),
        ani_line(genbank="GCA_2.1", refseq="GCF_2.1", species_taxid=200, best_taxid=200),  # our species, other assembly
        ani_line(genbank="GCA_3.1", refseq="GCF_3.1", species_taxid=300, best_taxid=300),  # unrelated
    ]
    # We want the GCF_1 pick (by its RefSeq accession) and count species 100 & 200.
    d = distill_ani_report(lines, wanted_bases={"GCF_1"}, wanted_species={100, 200})
    assert set(d.rows) == {"GCF_1"}
    assert d.rows["GCF_1"].genbank_accession == "GCA_1.1"
    assert d.species_present == {100: 1, 200: 1}   # 300 not counted
    assert d.n_scanned == 3


def test_distill_ani_matches_gca_pick_via_genbank_column():
    lines = [ANI_HEADER + "\n",
             ani_line(genbank="GCA_9.1", refseq="GCF_9.1", species_taxid=100, best_taxid=100)]
    d = distill_ani_report(lines, wanted_bases={"GCA_9"}, wanted_species={100})
    assert set(d.rows) == {"GCA_9"}   # matched on the GenBank column


def test_ani_header_schema_change_raises():
    bad = "# genbank-accession\tSOMETHING-ELSE\t" + "\t".join(["x"] * 23)
    with pytest.raises(ValueError, match="unexpected column layout"):
        distill_ani_report([bad + "\n"], wanted_bases=set(), wanted_species=set())


def test_ani_leading_description_comment_does_not_break_header():
    # A '## note' before the real column header must be skipped, not validated as the header.
    lines = ["## ANI report, generated nightly\n", ANI_HEADER + "\n",
             ani_line(genbank="GCA_1.1", refseq="GCF_1.1", species_taxid=100, best_taxid=100)]
    d = distill_ani_report(lines, wanted_bases={"GCF_1"}, wanted_species={100})
    assert set(d.rows) == {"GCF_1"}


def test_ani_data_before_validated_header_raises():
    # A malformed/absent column header (data arrives first) must fail LOUD, not silently
    # parse against the wrong layout.
    with pytest.raises(ValueError, match="before a validated column header"):
        distill_ani_report([ani_line(genbank="GCA_1.1", refseq="GCF_1.1", species_taxid=100)],
                           wanted_bases={"GCF_1"}, wanted_species={100})


def test_summary_data_before_validated_header_raises():
    with pytest.raises(ValueError, match="before a validated column header"):
        distill_assembly_summary([summary_line(acc="GCF_1.1", species_taxid=100)],
                                 wanted_bases={"GCF_1"}, wanted_species={100}, source_db="refseq")


def test_summary_twin_binding_rejects_foreign_repairing():
    # A foreign GenBank row (GCA_...Y) that names our GCF_...X as its paired accession
    # must NOT hijack our pick's row; only the canonical same-core twin binds.
    lines = [
        SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
        # foreign genome Y claims GCF_000000X as its RefSeq twin (stale/re-pairing)
        summary_line(acc="GCA_000000Y.1", species_taxid=100, paired="GCF_000000X.1",
                     excluded="derived from metagenome", source_db="genbank"),
        # the legitimate same-core row for our pick
        summary_line(acc="GCA_000000X.1", species_taxid=100, paired="GCF_000000X.1",
                     source_db="genbank"),
    ]
    d = distill_assembly_summary(lines, wanted_bases={"GCF_000000X"}, wanted_species={100},
                                 source_db="genbank")
    assert set(d.rows) == {"GCF_000000X"}
    # bound to the canonical GCA_...X row, NOT the foreign GCA_...Y row
    assert d.rows["GCF_000000X"].assembly_accession == "GCA_000000X.1"
    assert d.rows["GCF_000000X"].excluded_from_refseq == ""


# --------------------------------------------------------------------------- #
#  assembly_summary distiller                                                  #
# --------------------------------------------------------------------------- #
def test_distill_summary_matches_via_paired_twin_and_counts_latest():
    lines = [
        SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
        summary_line(acc="GCF_1.1", species_taxid=100, paired="GCA_1.1", category="reference genome"),
        summary_line(acc="GCF_2.1", species_taxid=100, paired="GCA_2.1", version_status="suppressed"),
        summary_line(acc="GCF_3.1", species_taxid=999, paired="GCA_3.1"),
    ]
    # Pick given as GCA_1 must match via gbrs_paired_asm; species-100 latest count = 1
    # (the suppressed GCF_2 is not counted as a live alternative).
    d = distill_assembly_summary(lines, wanted_bases={"GCA_1"}, wanted_species={100}, source_db="refseq")
    assert set(d.rows) == {"GCA_1"}
    assert d.rows["GCA_1"].refseq_category == "reference genome"
    assert d.species_present == {100: 1}
    assert 100 in d.species_reference


def test_summary_header_schema_change_raises():
    bad = "#assembly_accession\tbioproject\tNOPE\t" + "\t".join(["x"] * 23)
    with pytest.raises(ValueError, match="unexpected column layout"):
        distill_assembly_summary([SUMMARY_HEADER_1 + "\n", bad + "\n"],
                                 wanted_bases=set(), wanted_species=set(), source_db="refseq")


# --------------------------------------------------------------------------- #
#  Cache build (injected openers) + roundtrip + offline coverage               #
# --------------------------------------------------------------------------- #
def _openers():
    ani = [ANI_HEADER + "\n",
           ani_line(genbank="GCA_1.1", refseq="GCF_1.1", species_taxid=100, best_taxid=100)]
    refseq = [SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
              summary_line(acc="GCF_1.1", species_taxid=100, paired="GCA_1.1")]
    genbank = [SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
               summary_line(acc="GCA_1.1", species_taxid=100, paired="GCF_1.1", source_db="genbank")]
    return (lambda: iter(ani), lambda: iter(refseq), lambda: iter(genbank))


def test_cache_build_roundtrip_and_offline(tmp_path):
    a, r, g = _openers()
    path = tmp_path / "cache.json"
    c = AssemblyReportsCache(cache_path=str(path), ani_opener=a, refseq_opener=r, genbank_opener=g)
    c.build(wanted_bases={"GCF_1"}, wanted_species={100})
    # The ANI row is keyed under the pick AND its learned GCA twin (twin enrichment).
    assert set(c.ani) == {"GCF_1", "GCA_1"}
    assert set(c.refseq) == {"GCF_1"}
    assert c.twin_of == {"GCF_1": "GCA_1"}
    c.save()
    assert path.exists()

    # Fresh resolver loads the cache; covers() checks bases AND species.
    c2 = AssemblyReportsCache(cache_path=str(path))
    assert c2.covers({"GCF_1"}, {100})
    assert not c2.covers({"GCF_1", "GCF_zzz"})
    assert not c2.covers({"GCF_1"}, {999})       # species not covered -> stale cache
    assert c2.twin_of == {"GCF_1": "GCA_1"}
    assert c2.ani["GCF_1"].genbank_accession == "GCA_1.1"


def test_ani_row_matched_via_twin_when_refseq_column_blank(tmp_path):
    # NCBI's ANI row lists only the GenBank accession (refseq column blank). A GCF
    # pick must still be found via its GCA twin learned from the summary pairing.
    ani = [ANI_HEADER + "\n",
           ani_line(genbank="GCA_5.1", refseq="na", species_taxid=100, best_taxid=100)]
    refseq = [SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
              summary_line(acc="GCF_5.1", species_taxid=100, paired="GCA_5.1")]
    genbank = [SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
               summary_line(acc="GCA_5.1", species_taxid=100, paired="GCF_5.1", source_db="genbank")]
    c = AssemblyReportsCache(cache_path=None,
                             ani_opener=lambda: iter(ani),
                             refseq_opener=lambda: iter(refseq),
                             genbank_opener=lambda: iter(genbank))
    c.build(wanted_bases={"GCF_5"}, wanted_species={100})
    assert c.twin_of == {"GCF_5": "GCA_5"}
    assert "GCA_5" in c.ani            # ANI row stored under the twin base
    # and the coverage layer resolves the GCF pick through the twin (see test_coverage).


def test_build_never_touches_default_network_opener(tmp_path):
    # Sanity: constructing with injected openers must not call the URL default.
    called = {"url": False}

    def boom():
        called["url"] = True
        raise AssertionError("network used")

    a, r, g = _openers()
    c = AssemblyReportsCache(cache_path=None, ani_opener=a, refseq_opener=r, genbank_opener=g)
    c.build(wanted_bases={"GCF_1"}, wanted_species={100})
    assert called["url"] is False
