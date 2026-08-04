# -*- coding: utf-8 -*-
"""Pure J1a coverage logic: join, GCF<->GCA equivalence, version handling,
taxon-only fallback, and the precision-tuned flag rules (offline, no network)."""

from ani_reports import AniReportRow, AssemblySummaryRow, AssemblyReportsCache
import coverage
from coverage import Pick


# --------------------------------------------------------------------------- #
#  Builders                                                                    #
# --------------------------------------------------------------------------- #
def ani_row(genbank="", refseq="", species_taxid=None, species_name="",
            best_taxid=None, best_name="", best_status="species-match",
            tax_check="OK", excluded="", comment=""):
    return AniReportRow(
        genbank_accession=genbank, refseq_accession=refseq,
        species_taxid=species_taxid, species_name=species_name,
        best_match_species_taxid=best_taxid, best_match_species_name=best_name,
        best_match_status=best_status, taxonomy_check_status=tax_check,
        excluded_from_refseq=excluded, comment=comment,
    )


def summary_row(acc, category="", taxid=None, species_taxid=None,
                version_status="latest", seq_rel_date="2018/01/01", paired="",
                excluded="", relation="", source_db="refseq"):
    return AssemblySummaryRow(
        assembly_accession=acc, refseq_category=category, taxid=taxid,
        species_taxid=species_taxid, version_status=version_status,
        seq_rel_date=seq_rel_date, gbrs_paired_asm=paired,
        excluded_from_refseq=excluded, relation_to_type_material=relation,
        source_db=source_db,
    )


def pick(name, parent, species_taxid=100, accs=None, lane="NCBI-REF", regime="", resolver=""):
    return Pick(
        unit_id=f"{name}~{parent}", scientific_name=name, species_taxid=species_taxid,
        species_name="", parent_assembly=parent, accession_set=accs or ["CP0001"],
        source_lane=lane, regime=regime, resolver=resolver,
    )


def cache_with(ani=None, refseq=None, genbank=None, ani_species=None, summary_species=None):
    c = AssemblyReportsCache(cache_path=None)
    c.ani = dict(ani or {})
    c.refseq = dict(refseq or {})
    c.genbank = dict(genbank or {})
    c.ani_species_present = dict(ani_species or {})
    c.summary_species_present = dict(summary_species or {})
    c.wanted_bases = set(c.ani) | set(c.refseq) | set(c.genbank)
    return c


def one(cache, p):
    return coverage.build_coverage([p], cache)[0]


# --------------------------------------------------------------------------- #
#  Join + match typing                                                         #
# --------------------------------------------------------------------------- #
def test_exact_match_on_refseq_accession():
    c = cache_with(
        ani={"GCF_1": ani_row(genbank="GCA_1.1", refseq="GCF_1.1", species_taxid=100,
                              best_taxid=100)},
        refseq={"GCF_1": summary_row("GCF_1.1", species_taxid=100, paired="GCA_1.1")},
    )
    r = one(c, pick("Foo", "GCF_1.1", species_taxid=100))
    assert r.ani_match == coverage.MATCH_EXACT
    assert r.summary_match == coverage.MATCH_EXACT
    assert r.in_refseq is True
    assert r.actionable_flags == []


def test_gcf_gca_equivalence_pick_stored_as_gca():
    # Our pick is a GCA; the ANI row carries its GCF twin, the summary its paired GCF.
    c = cache_with(
        ani={"GCA_9": ani_row(genbank="GCA_9.1", refseq="GCF_9.1", species_taxid=100,
                              best_taxid=100)},
        refseq={"GCA_9": summary_row("GCF_9.1", species_taxid=100, paired="GCA_9.1")},
    )
    r = one(c, pick("Bar", "GCA_9.1", species_taxid=100))
    assert r.ani_match == coverage.MATCH_EXACT
    assert r.summary_match == coverage.MATCH_EXACT
    assert r.in_refseq is True


def test_version_differs_flags_superseded():
    # We hold .1 but NCBI's live row is .2 -> same base assembly, stale version.
    c = cache_with(
        ani={"GCF_2": ani_row(genbank="GCA_2.2", refseq="GCF_2.2", species_taxid=100,
                              best_taxid=100)},
        refseq={"GCF_2": summary_row("GCF_2.2", species_taxid=100, paired="GCA_2.2")},
    )
    r = one(c, pick("Baz", "GCF_2.1", species_taxid=100))
    assert r.ani_match == coverage.MATCH_VERSION
    assert r.summary_match == coverage.MATCH_VERSION
    assert "assembly_superseded" in r.flags
    assert "assembly_superseded" not in r.actionable_flags  # informational


def test_taxon_only_when_exact_absent_but_species_present():
    c = cache_with(ani={}, refseq={}, ani_species={100: 7}, summary_species={100: 7})
    r = one(c, pick("Qux", "GCF_404.1", species_taxid=100))
    assert r.ani_match == coverage.MATCH_TAXON
    assert r.summary_match == coverage.MATCH_TAXON
    assert "ani_taxon_only" in r.flags


def test_none_when_absent_and_species_absent():
    c = cache_with(ani={}, refseq={})
    r = one(c, pick("Nope", "GCF_500.1", species_taxid=999))
    assert r.ani_match == coverage.MATCH_NONE
    assert "not_in_ani_report" in r.flags
    assert "not_in_assembly_summary" in r.flags


# --------------------------------------------------------------------------- #
#  Flag logic (precision-tuned)                                                #
# --------------------------------------------------------------------------- #
def test_tax_check_failed_is_actionable():
    c = cache_with(
        ani={"GCF_3": ani_row(refseq="GCF_3.1", species_taxid=100, best_taxid=200,
                              best_name="Other sp.", best_status="mismatch", tax_check="Failed")},
        refseq={"GCF_3": summary_row("GCF_3.1", species_taxid=100)},
    )
    r = one(c, pick("Fail", "GCF_3.1", species_taxid=100))
    assert "tax_check_failed" in r.actionable_flags


def test_tax_check_inconclusive_is_actionable():
    c = cache_with(
        ani={"GCF_4": ani_row(refseq="GCF_4.1", species_taxid=100, best_taxid=100,
                              tax_check="Inconclusive")},
        refseq={"GCF_4": summary_row("GCF_4.1", species_taxid=100)},
    )
    r = one(c, pick("Incon", "GCF_4.1", species_taxid=100))
    assert "tax_check_inconclusive" in r.actionable_flags


def test_best_match_mismatch_with_ok_status_is_context_not_actionable():
    # The precision-first rule: a divergent best-match with taxonomy-check OK defers
    # to NCBI's verdict -- captured as context, never on the worklist.
    c = cache_with(
        ani={"GCF_5": ani_row(refseq="GCF_5.1", species_taxid=100, best_taxid=200,
                              best_name="Sister sp.", best_status="species-match", tax_check="OK")},
        refseq={"GCF_5": summary_row("GCF_5.1", species_taxid=100)},
    )
    r = one(c, pick("Diverge", "GCF_5.1", species_taxid=100))
    assert "best_match_species_mismatch" in r.flags
    assert r.actionable_flags == []


def test_ncbi_reclassified_is_actionable():
    # NCBI now files the assembly under a different species than we declare.
    c = cache_with(
        ani={"GCF_6": ani_row(refseq="GCF_6.1", species_taxid=222, species_name="New sp.",
                              best_taxid=222, tax_check="OK")},
        refseq={"GCF_6": summary_row("GCF_6.1", species_taxid=222)},
    )
    r = one(c, pick("Reclass", "GCF_6.1", species_taxid=111))
    assert "ncbi_reclassified" in r.actionable_flags


def test_excluded_from_refseq_identity_vs_qc():
    assert coverage.classify_exclusion("unverified source organism") == "identity"
    assert coverage.classify_exclusion("contaminated") == "identity"
    assert coverage.classify_exclusion("derived from metagenome") == "identity"
    assert coverage.classify_exclusion("many frameshifted proteins") == "qc"
    assert coverage.classify_exclusion("misassembled") == "qc"          # assembly QC, not identity
    assert coverage.classify_exclusion("sequence duplications") == "qc"
    assert coverage.classify_exclusion("") is None

    # identity exclusion is actionable; qc exclusion is informational.
    c = cache_with(
        ani={"GCF_7": ani_row(refseq="GCF_7.1", species_taxid=100, best_taxid=100)},
        genbank={"GCF_7": summary_row("GCA_7.1", species_taxid=100, source_db="genbank",
                                      excluded="unverified source organism", paired="GCF_7.1")},
    )
    r = one(c, pick("Excl", "GCF_7.1", species_taxid=100))
    assert r.in_refseq is False               # not in the refseq summary
    assert "refseq_excluded_identity" in r.actionable_flags

    c2 = cache_with(
        ani={"GCF_8": ani_row(refseq="GCF_8.1", species_taxid=100, best_taxid=100)},
        genbank={"GCF_8": summary_row("GCA_8.1", species_taxid=100, source_db="genbank",
                                      excluded="many frameshifted proteins", paired="GCF_8.1")},
    )
    r2 = one(c2, pick("ExclQC", "GCF_8.1", species_taxid=100))
    assert "refseq_excluded_qc" in r2.flags
    assert r2.actionable_flags == []


def test_approved_mismatch_demotes_species_flags_only():
    # NCBI approved the SPECIES mismatch, so the species-identity flags are demoted --
    # but NCBI's INDEPENDENT taxonomy-check verdict is a separate axis and survives.
    # Normal case: approved + taxonomy-check OK -> nothing actionable ("never flag").
    c_ok = cache_with(
        ani={"GCF_9": ani_row(refseq="GCF_9.1", species_taxid=100, best_taxid=200,
                              best_name="Other", best_status=coverage.APPROVED_MISMATCH,
                              tax_check="OK")},
        refseq={"GCF_9": summary_row("GCF_9.1", species_taxid=100)},
    )
    r = one(c_ok, pick("Approved", "GCF_9.1", species_taxid=100))
    assert "approved_mismatch" in r.flags
    assert "best_match_species_mismatch" not in r.flags
    assert r.actionable_flags == []

    # Pathological case: approved + taxonomy-check Failed -> the species flag is demoted
    # but the independent Failed verdict is NOT erased.
    c_fail = cache_with(
        ani={"GCF_9b": ani_row(refseq="GCF_9b.1", species_taxid=100, best_taxid=200,
                               best_name="Other", best_status=coverage.APPROVED_MISMATCH,
                               tax_check="Failed")},
        refseq={"GCF_9b": summary_row("GCF_9b.1", species_taxid=100)},
    )
    r2 = one(c_fail, pick("ApprovedFail", "GCF_9b.1", species_taxid=100))
    assert "best_match_species_mismatch" not in r2.flags   # species mismatch approved
    assert "tax_check_failed" in r2.actionable_flags        # independent verdict survives


def test_confounded_keeps_reclassification_actionable():
    # ncbi_reclassified is from NCBI's assigned taxid (metadata), not ANI similarity, so
    # it is NOT the "ANI can't split" artifact and must survive the confounded demotion.
    c = cache_with(
        ani={"GCF_rc": ani_row(refseq="GCF_rc.1", species_taxid=999, species_name="Other sp.",
                               best_taxid=999, tax_check="Inconclusive")},
        refseq={"GCF_rc": summary_row("GCF_rc.1", species_taxid=999)},
    )
    r = one(c, pick("Confounded_reclassed", "GCF_rc.1", species_taxid=1396,
                    regime="A", resolver="marker"))
    assert "ncbi_reclassified" in r.actionable_flags        # metadata signal survives
    assert "confounded_ani_ambiguity" in r.flags            # the ANI part still demoted
    assert "tax_check_inconclusive" not in r.flags


def test_twin_only_match_is_exact_not_version_differs():
    # Regression for the review finding: a GCF pick found only via its GCA twin must NOT
    # be version-typed across the GCF<->GCA boundary (independent version numbers).
    c = cache_with(
        ani={"GCA_x": ani_row(genbank="GCA_x.3", refseq="", species_taxid=100, best_taxid=100)},
        refseq={"GCF_x": summary_row("GCF_x.2", species_taxid=100, paired="GCA_x.3")},
    )
    c.twin_of = {"GCF_x": "GCA_x"}
    r = one(c, pick("TwinVer", "GCF_x.2", species_taxid=100))   # our v2 vs GCA v3
    assert r.ani_match == coverage.MATCH_EXACT                  # NOT version-differs
    assert "assembly_superseded" not in r.flags


def test_genbank_only_pick_is_flagged():
    # A GenBank-only pick (no RefSeq twin, no exclusion reason) -> informational note.
    c = cache_with(
        ani={"GCA_g": ani_row(genbank="GCA_g.1", refseq="", species_taxid=100, best_taxid=100)},
        genbank={"GCA_g": summary_row("GCA_g.1", species_taxid=100, source_db="genbank", paired="")},
    )
    r = one(c, pick("GbOnly", "GCA_g.1", species_taxid=100))
    assert r.in_refseq is False
    assert "genbank_only" in r.flags
    assert r.actionable_flags == []


def test_confounded_identity_ambiguity_demoted():
    # A confounded (marker-resolved) taxon: the expected ANI/type-strain mismatch
    # collapses to one informational note, never the worklist.
    c = cache_with(
        ani={"GCF_10": ani_row(refseq="GCF_10.1", species_taxid=1396, best_taxid=999,
                               best_name="Bacillus paranthracis", best_status="mismatch",
                               tax_check="Failed")},
        refseq={"GCF_10": summary_row("GCF_10.1", species_taxid=1396)},
    )
    r = one(c, pick("Bacillus_cereus", "GCF_10.1", species_taxid=1396,
                    regime="A", resolver="marker"))
    assert "confounded_ani_ambiguity" in r.flags
    assert r.actionable_flags == []
    assert "tax_check_failed" not in r.flags


def test_confounded_identity_exclusion_stays_actionable():
    # Regression for the review finding: a REAL withdrawal from RefSeq for a
    # contamination/provenance reason must survive the confounded demotion (it is
    # orthogonal to "ANI can't split this pair"), so a contaminated Salmonella
    # reference is not silently hidden as ANI ambiguity.
    c = cache_with(
        ani={"GCF_12": ani_row(refseq="GCF_12.1", species_taxid=28901, best_taxid=28901,
                               tax_check="Inconclusive")},
        genbank={"GCF_12": summary_row("GCA_12.1", species_taxid=28901, source_db="genbank",
                                       excluded="contaminated", paired="GCF_12.1")},
    )
    r = one(c, pick("Salmonella_lineage", "GCF_12.1", species_taxid=28901, regime="B", resolver="ANI"))
    assert "refseq_excluded_identity" in r.actionable_flags   # NOT demoted away
    assert "confounded_ani_ambiguity" in r.flags               # the ANI part still demoted
    assert "tax_check_inconclusive" not in r.flags


def test_assembly_suppressed_and_replaced_are_actionable():
    # version_status is 100% 'latest' in today's data, so these paths are only guarded
    # by this test -- and they feed the make-or-break actionable count.
    c_sup = cache_with(
        ani={"GCF_s": ani_row(refseq="GCF_s.1", species_taxid=100, best_taxid=100)},
        refseq={"GCF_s": summary_row("GCF_s.1", species_taxid=100, version_status="suppressed")},
    )
    assert "assembly_suppressed" in one(c_sup, pick("Sup", "GCF_s.1", species_taxid=100)).actionable_flags

    c_rep = cache_with(
        ani={"GCF_r": ani_row(refseq="GCF_r.1", species_taxid=100, best_taxid=100)},
        refseq={"GCF_r": summary_row("GCF_r.1", species_taxid=100, version_status="replaced")},
    )
    assert "assembly_replaced" in one(c_rep, pick("Rep", "GCF_r.1", species_taxid=100)).actionable_flags


def test_superseded_only_when_file_is_newer():
    # We hold .1, file has .2 -> superseded.
    c_new = cache_with(
        ani={"GCF_v": ani_row(genbank="GCA_v.2", refseq="GCF_v.2", species_taxid=100, best_taxid=100)},
        refseq={"GCF_v": summary_row("GCF_v.2", species_taxid=100, paired="GCA_v.2")},
    )
    assert "assembly_superseded" in one(c_new, pick("New", "GCF_v.1", species_taxid=100)).flags
    # We hold .2, file has .1 (a lagging row) -> version differs but NOT superseded.
    c_old = cache_with(
        ani={"GCF_w": ani_row(genbank="GCA_w.1", refseq="GCF_w.1", species_taxid=100, best_taxid=100)},
        refseq={"GCF_w": summary_row("GCF_w.1", species_taxid=100, paired="GCA_w.1")},
    )
    r = one(c_old, pick("Old", "GCF_w.2", species_taxid=100))
    assert r.ani_match == coverage.MATCH_VERSION
    assert "assembly_superseded" not in r.flags


def test_ncbi_reclassified_ignores_same_name_taxid_merge():
    # Different taxid but the SAME species name -> a taxid merge, not a reclassification.
    c = cache_with(
        ani={"GCF_m": ani_row(refseq="GCF_m.1", species_taxid=222, species_name="Foo bar",
                              best_taxid=222, tax_check="OK")},
        refseq={"GCF_m": summary_row("GCF_m.1", species_taxid=222)},
    )
    p = pick("Merge", "GCF_m.1", species_taxid=111)
    p.species_name = "Foo bar"
    r = one(c, p)
    assert "ncbi_reclassified" not in r.flags


def test_ani_matched_via_twin_fallback():
    # ANI row keyed only under the GCA twin (RefSeq column was blank); the GCF pick
    # resolves through twin_of instead of being mislabeled taxon-only.
    c = cache_with(
        ani={"GCA_t": ani_row(genbank="GCA_t.1", refseq="", species_taxid=100, best_taxid=100)},
        refseq={"GCF_t": summary_row("GCF_t.1", species_taxid=100, paired="GCA_t.1")},
    )
    c.twin_of = {"GCF_t": "GCA_t"}
    r = one(c, pick("Twin", "GCF_t.1", species_taxid=100))
    assert r.ani_match == coverage.MATCH_EXACT
    assert "not_in_ani_report" not in r.flags


def test_unresolved_assembly_proxy_is_flagged_and_bucketed_other():
    # A Tier-0 'strain:NNN' proxy parent cannot join; surface it, don't mis-bucket as GCA.
    c = cache_with()
    p = pick("Unresolved", "strain:12345", species_taxid=100)
    r = one(c, p)
    assert "unresolved_assembly" in r.flags
    assert r.ani_match == coverage.MATCH_NONE
    s = coverage.summarize([r])
    assert "other" in s["by_assembly_kind"] and "GCA" not in s["by_assembly_kind"]


def test_missing_species_taxid_is_graceful_no_taxon_fallback():
    # Blank species_taxid: an exact assembly match still resolves; identity flags stay
    # silent; and with no taxid there is deliberately no taxon-only fallback.
    c_hit = cache_with(
        ani={"GCF_n": ani_row(refseq="GCF_n.1", species_taxid=100, best_taxid=100)},
        refseq={"GCF_n": summary_row("GCF_n.1", species_taxid=100)},
    )
    r = one(c_hit, pick("NoTax", "GCF_n.1", species_taxid=None))
    assert r.ani_match == coverage.MATCH_EXACT
    assert r.actionable_flags == []

    c_miss = cache_with(ani={}, refseq={}, ani_species={100: 5}, summary_species={100: 5})
    r2 = one(c_miss, pick("NoTaxMiss", "GCF_z.1", species_taxid=None))
    assert r2.ani_match == coverage.MATCH_NONE       # no taxid -> no taxon-only fallback


def test_summarize_empty_input_is_safe():
    s = coverage.summarize([])
    assert s["n_picks"] == 0
    assert s["headline"]["ani_live_pct"] == 0.0
    assert s["by_source_lane"] == {}
    assert s["actionable_units"] == 0


def test_source_lane_mixed_and_unknown(tmp_path):
    chrom = tmp_path / "chromosomes.tsv"
    chrom.write_text(
        "scientificName\tnuccoreAcc\ttaxid\tparent\tsource\n"
        "M\tCP1\t100\t1\tCDC SME\n"
        "M\tCP2\t100\t1\tFDA-ARGOS\n"      # same unit, two lanes -> 'mixed'
    )
    policy = tmp_path / "policy.tsv"
    policy.write_text(
        "unit_id\tscientificName\tkalamari_taxid\tspecies_taxid\taccession_set\tparent_assembly\t"
        "eligibility\tregime\tresolver\tspecies_name\tsuperkingdom\ttotal_length_bp\tn_replicons\tflags\tnote\n"
        "M~GCF_1.1\tM\t100\t100\tCP1,CP2\tGCF_1.1\teligible\t\t\tSp.\tBacteria\t3000000\t2\t\t\n"
        "U~GCF_2.1\tU\t100\t100\tCP9\tGCF_2.1\teligible\t\t\tSp.\tBacteria\t3000000\t1\t\t\n"
    )
    lanes = coverage.load_source_lanes(str(chrom))
    picks = coverage.load_picks(str(policy), lanes)
    by_name = {p.scientific_name: p for p in picks}
    assert by_name["M"].source_lane == "mixed"
    assert by_name["U"].source_lane == "unknown"     # CP9 has no chromosomes.tsv row


def test_confounded_qc_exclusion_survives_demotion():
    # A confounded taxon still surfaces an orthogonal annotation-QC exclusion note
    # (it isn't an identity signal, so demotion leaves it alone) -- but still
    # non-actionable.
    c = cache_with(
        ani={"GCF_11": ani_row(refseq="GCF_11.1", species_taxid=28901, best_taxid=28901)},
        genbank={"GCF_11": summary_row("GCA_11.1", species_taxid=28901, source_db="genbank",
                                       excluded="many frameshifted proteins", paired="GCF_11.1")},
    )
    r = one(c, pick("Salmonella_lineage", "GCF_11.1", species_taxid=28901, regime="B", resolver="ANI"))
    assert "refseq_excluded_qc" in r.flags
    assert r.actionable_flags == []


# --------------------------------------------------------------------------- #
#  Source lanes + summary                                                      #
# --------------------------------------------------------------------------- #
def test_source_lane_collapse(tmp_path):
    chrom = tmp_path / "chromosomes.tsv"
    chrom.write_text(
        "scientificName\tnuccoreAcc\ttaxid\tparent\tsource\n"
        "A\tCP1\t100\t1\tCDC SME\n"
        "B\tCP2\t100\t1\tNCBI-REF\n"
        "C\tCP3\t100\t1\tFDA-ARGOS\n"
        "D\tCP4\t100\t1\tUGA SME\n"
    )
    lanes = coverage.load_source_lanes(str(chrom))
    assert lanes["CP1"] == "SME"
    assert lanes["CP2"] == "NCBI-REF"
    assert lanes["CP3"] == "FDA-ARGOS"
    assert lanes["CP4"] == "SME"


def test_summarize_headline_and_breakdowns():
    c = cache_with(
        ani={
            "GCF_a": ani_row(refseq="GCF_a.1", species_taxid=100, best_taxid=100),
            "GCF_b": ani_row(refseq="GCF_b.1", species_taxid=100, best_taxid=100, tax_check="Failed"),
        },
        refseq={
            "GCF_a": summary_row("GCF_a.1", species_taxid=100),
            "GCF_b": summary_row("GCF_b.1", species_taxid=100),
        },
        ani_species={200: 3}, summary_species={200: 3},
    )
    picks = [
        pick("A", "GCF_a.1", species_taxid=100, lane="NCBI-REF"),
        pick("B", "GCF_b.1", species_taxid=100, lane="SME"),
        pick("C", "GCF_c.1", species_taxid=200, lane="SME"),  # taxon-only
    ]
    rows = coverage.build_coverage(picks, c)
    s = coverage.summarize(rows)
    assert s["n_picks"] == 3
    assert s["headline"]["ani_live"] == 2          # A + B exact; C only taxon-only
    assert s["headline"]["ani_live_pct"] == 66.7
    assert s["ani_match_types"]["taxon-only"] == 1
    assert s["by_source_lane"]["SME"]["n"] == 2
    assert s["actionable_units"] == 1              # only B (tax_check_failed)
    assert s["flag_counts"]["tax_check_failed"] == 1
