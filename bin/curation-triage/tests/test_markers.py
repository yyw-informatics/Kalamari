# -*- coding: utf-8 -*-
"""Tier 2: the marker manifest and our versioned interpretation rules.

The rules are where a tool's output becomes a judgement, so each test names the
mistake it exists to prevent.
"""

import os

import pytest

import markers

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
MANIFEST_PATH = os.path.join(REPO_ROOT, "src", "curation-triage", "marker_manifest.tsv")


@pytest.fixture(scope="module")
def manifest():
    return markers.load_manifest(MANIFEST_PATH)


# --------------------------------------------------------------------------- #
#  The committed manifest                                                      #
# --------------------------------------------------------------------------- #
def test_every_pinned_row_names_a_tool_and_a_rule(manifest):
    assert manifest, "manifest did not load"
    for row in manifest:
        assert row.interpretation_rule_version, "%s/%s has no rule version" % (
            row.task_id, row.component)
        if not row.is_gap:
            assert row.tool, "%s/%s is pinned but names no tool" % (row.task_id, row.component)
            assert row.database_version, "%s/%s pins no database version" % (
                row.task_id, row.component)


def test_the_four_researched_tasks_are_pinned_to_the_2026_08_02_decisions(manifest):
    """Three of the four moved away from abricate as the classifier; keep them there."""
    def component(task_id, name):
        return next(r for r in markers.rows_for_task(manifest, task_id) if r.component == name)

    assert component("bacillus_anthracis_vs_b_cereus_group",
                     "group_taxonomy_and_marker_profile").conda_build.startswith("3.4.0")
    yersinia = component("yersinia_pestis_vs_y_pseudotuberculosis", "cgmlst_species_assignment")
    # The research pass named BIGSdb cgMLST the classifier, but nothing can run it:
    # no packaged CLI, no pinned allele caller.  It is recorded as a gap so the task
    # returns "undecided" instead of falling back to plasmid presence.
    assert "BIGSdb" in yersinia.tool and yersinia.is_gap
    assert component("listeria_monocytogenes_serogroup_lineage_cc",
                     "serogroup").conda_build.startswith("0.4.10")
    bont = component("clostridium_botulinum_group_and_bont_type", "bont_type")
    assert bont.conda_build.startswith("4.2.7") and bont.database_version == "2026-05-15.1"


def test_abricate_survives_only_as_a_pre_screen_for_yersinia(manifest):
    rows = markers.rows_for_task(manifest, "yersinia_pestis_vs_y_pseudotuberculosis")
    abricate = next(r for r in rows if r.tool == "abricate")
    assert abricate.role == markers.ROLE_PRE_SCREEN


def test_the_botulinum_group_classifier_is_recorded_as_a_gap(manifest):
    """An unpinned tool must be visible as a gap, not silently missing."""
    rows = markers.rows_for_task(manifest, "clostridium_botulinum_group_and_bont_type")
    gap = next(r for r in rows if r.component == "genome_group")
    assert gap.is_gap and not gap.tool


def test_listeria_keeps_three_separate_components(manifest):
    rows = markers.rows_for_task(manifest, "listeria_monocytogenes_serogroup_lineage_cc")
    assert {r.component for r in rows} == {"serogroup", "mlst_st_cc", "cgmlst_lin"}


def test_a_unit_routes_to_its_task_by_species_or_group_taxid(manifest):
    assert markers.tasks_for_taxids(manifest, [1639]) == \
        ["listeria_monocytogenes_serogroup_lineage_cc"]
    # Our Bacillus units are B. cereus (1396) / B. thuringiensis (1428); the manifest
    # triggers on the B. cereus GROUP (86661), which Tier 0 resolved for us.
    assert markers.tasks_for_taxids(manifest, [1396]) == []
    assert markers.tasks_for_taxids(manifest, [1396, 86661]) == \
        ["bacillus_anthracis_vs_b_cereus_group"]


def test_the_generic_fallback_never_triggers_on_its_own(manifest):
    assert "generic_marker_fallback" not in markers.tasks_for_taxids(manifest, [0, 562])


def test_pins_carry_software_database_and_rule_together(manifest):
    pins = markers.task_pins(manifest, "bacillus_anthracis_vs_b_cereus_group")
    comp = pins["components"][0]
    assert comp["method"]["conda_package"] == "btyper3"
    assert comp["database"]["version_or_scheme_id"]
    assert pins["interpretation_rule_version"] == "bacillus_anthracis_vs_b_cereus_group@1"
    assert pins["regime"] == markers.REGIME_A


# --------------------------------------------------------------------------- #
#  Rule: Shigella / EIEC vs E. coli                                            #
# --------------------------------------------------------------------------- #
def test_ipah_absent_confirms_an_e_coli_pick():
    call = markers.call_shigella_eiec({"ipaH": "absent"}, "Escherichia coli")
    assert call.verdict == markers.VERDICT_CONFIRMED and call.overrides_ani


def test_an_e_coli_pick_that_types_as_shigella_is_refuted():
    call = markers.call_shigella_eiec(
        {"ipaH": "present", "organism": "Shigella sonnei", "cluster": "C3"},
        "Escherichia coli")
    assert call.verdict == markers.VERDICT_REFUTED
    assert call.call == "Shigella sonnei"


def test_ipah_without_a_cluster_does_not_guess_between_shigella_and_eiec():
    call = markers.call_shigella_eiec({"ipaH": "present"}, "Escherichia coli")
    assert call.verdict == markers.VERDICT_INDETERMINATE


# --------------------------------------------------------------------------- #
#  Rule: B. anthracis vs the B. cereus group                                   #
# --------------------------------------------------------------------------- #
def test_anthracis_chromosome_with_both_plasmids_is_anthracis():
    call = markers.call_bacillus_anthracis(
        {"genomospecies": "Bacillus anthracis", "anthracis_chromosome": "yes",
         "pxo1_genes": ["pagA", "lef", "cya"], "pxo2_genes": ["capA", "capB", "capC"]},
        "Bacillus anthracis")
    assert call.call == markers.ANTHRACIS
    assert call.verdict == markers.VERDICT_CONFIRMED


def test_plasmid_cured_anthracis_is_still_anthracis():
    """Losing pXO1/pXO2 does not make a select agent into B. cereus."""
    call = markers.call_bacillus_anthracis(
        {"genomospecies": "Bacillus anthracis", "anthracis_chromosome": "yes",
         "pxo1_genes": [], "pxo2_genes": []}, "Bacillus anthracis")
    assert call.call == markers.ANTHRACIS_CURED
    assert call.verdict == markers.VERDICT_CONFIRMED


def test_a_cereus_carrying_anthrax_plasmids_is_never_labelled_plain_cereus():
    call = markers.call_bacillus_anthracis(
        {"genomospecies": "Bacillus cereus", "anthracis_chromosome": "no",
         "pxo1_genes": ["pagA", "lef"], "pxo2_genes": []}, "Bacillus cereus")
    assert call.call == markers.CEREUS_ANTHRAX_TOXIN
    assert call.verdict == markers.VERDICT_REFUTED


def test_a_different_subspecies_refutes_the_declared_label():
    """BTyper3 3.4 carries the familiar name in the SUBSPECIES field."""
    call = markers.call_bacillus_anthracis(
        {"genomospecies": "B. mosaicus subsp. paranthracis",
         "btyper3_subspecies": "paranthracis", "anthracis_chromosome": "no",
         "pxo1_genes": [], "pxo2_genes": []}, "Bacillus cereus")
    assert call.verdict == markers.VERDICT_REFUTED


def test_no_assigned_subspecies_cannot_confirm_the_label_but_still_rules_out_anthracis():
    """Kalamari's B. thuringiensis reference (Al Hakam) carries no Bt toxin genes, so
    BTyper3 assigns it no biovar. 'Cannot confirm' is the honest answer, not 'wrong'."""
    call = markers.call_bacillus_anthracis(
        {"genomospecies": "B. mosaicus", "btyper3_subspecies": "",
         "closest_type_strain": "anthracis", "anthracis_chromosome": "no",
         "pxo1_genes": [], "pxo2_genes": [], "bt_genes": []},
        "Bacillus thuringiensis")
    assert call.verdict == markers.VERDICT_INDETERMINATE
    reasons = "; ".join(call.reasons)
    assert "no subspecies or biovar" in reasons
    assert "rule out anthracis" in reasons


def test_no_chromosome_evidence_gives_no_verdict():
    call = markers.call_bacillus_anthracis({}, "Bacillus cereus")
    assert call.verdict == markers.VERDICT_INDETERMINATE


# --------------------------------------------------------------------------- #
#  Rule: Y. pestis vs Y. pseudotuberculosis                                    #
# --------------------------------------------------------------------------- #
def test_cgmlst_over_enough_loci_makes_the_yersinia_call():
    call = markers.call_yersinia_pestis(
        {"cgmlst_species_assignment": "Yersinia pseudotuberculosis", "loci_called": "487"},
        "Yersinia pseudotuberculosis")
    assert call.verdict == markers.VERDICT_CONFIRMED and call.confidence == "high"


def test_plasmid_markers_alone_never_make_a_pestis_call():
    """The mistake the 2026-08-02 pinning pass removed: plasmid presence is not species."""
    call = markers.call_yersinia_pestis(
        {"plasmid_markers": ["caf1", "pst", "pla"]}, "Yersinia pseudotuberculosis")
    assert call.verdict == markers.VERDICT_INDETERMINATE
    assert "never a species call" in "; ".join(call.reasons)


def test_too_few_called_loci_falls_back_to_chromosome_targets():
    call = markers.call_yersinia_pestis(
        {"cgmlst_species_assignment": "Yersinia pestis", "loci_called": "120",
         "chromosome_markers": ["ypo2088"]}, "Yersinia pseudotuberculosis")
    assert call.call == markers.PESTIS
    assert call.verdict == markers.VERDICT_REFUTED
    assert call.confidence == "medium"


def test_conflicting_chromosome_targets_give_no_verdict():
    call = markers.call_yersinia_pestis(
        {"chromosome_markers": ["ypo2088", "opgG"]}, "Yersinia pestis")
    assert call.verdict == markers.VERDICT_INDETERMINATE


# --------------------------------------------------------------------------- #
#  Rule: Listeria serogroup / lineage / CC                                     #
# --------------------------------------------------------------------------- #
def test_the_listeria_fields_stay_separate():
    call = markers.call_listeria_lineage(
        {"serogroup": "4b", "mlst_st": "1", "mlst_cc": "CC1", "cgmlst_lin": "L1.2",
         "broad_lineage": "I"}, "Listeria_monocytogenes_I")
    assert call.verdict == markers.VERDICT_CONFIRMED
    for key in ("predicted_molecular_serogroup", "mlst_st", "mlst_clonal_complex",
                "lin_code", "broad_lineage"):
        assert key in call.fields


def test_a_serogroup_alone_is_never_a_lineage():
    call = markers.call_listeria_lineage({"serogroup": "4b"}, "Listeria_monocytogenes_I")
    assert call.verdict == markers.VERDICT_INDETERMINATE
    assert "serogroup alone cannot decide" in "; ".join(call.reasons)


def test_a_wrong_lineage_is_refuted_and_never_overrides_ani():
    call = markers.call_listeria_lineage(
        {"broad_lineage": "II", "mlst_cc": "CC7"}, "Listeria_monocytogenes_III")
    assert call.verdict == markers.VERDICT_REFUTED
    assert call.overrides_ani is False


# --------------------------------------------------------------------------- #
#  Rule: Salmonella subspecies / serovar                                       #
# --------------------------------------------------------------------------- #
def test_salmonella_subspecies_is_checked_against_the_kalamari_suffix():
    call = markers.call_salmonella(
        {"subspecies": "arizonae", "sistr_serovar": "IIIa 41:z4,z23:-",
         "seqsero2_serovar": "IIIa 41:z4,z23:-"}, "Salmonella_enterica_IIIa")
    assert call.verdict == markers.VERDICT_CONFIRMED


def test_a_kalamari_subspecies_with_no_accepted_name_abstains():
    """IIa/IIb/VIII/IX/X have no name for SISTR to return; do not force a match."""
    call = markers.call_salmonella({"subspecies": "enterica"}, "Salmonella_enterica_IX")
    assert call.verdict == markers.VERDICT_INDETERMINATE


def test_two_typers_disagreeing_is_recorded_but_does_not_refute_the_subspecies():
    call = markers.call_salmonella(
        {"subspecies": "arizonae", "sistr_serovar": "IIIa 41:z4,z23:-",
         "seqsero2_serovar": "IIIb 61:k:1,5,(7)"}, "Salmonella_enterica_IIIa")
    assert call.verdict == markers.VERDICT_CONFIRMED
    assert call.confidence == "medium"
    assert "disagree" in "; ".join(call.reasons)


# --------------------------------------------------------------------------- #
#  Rule: C. botulinum group + bont type                                        #
# --------------------------------------------------------------------------- #
def test_the_toxin_type_is_recorded_but_the_group_is_not_inferred_from_it():
    call = markers.call_botulinum({"bont_types": ["B"], "bont_subtype": "B1"},
                                  "Clostridium_botulinum_groupI")
    assert call.fields["bont_gene_type"] == "B"
    assert call.verdict == markers.VERDICT_INDETERMINATE
    assert "no pinned genome-wide group classifier" in "; ".join(call.reasons)


def test_a_botulinum_reference_with_no_toxin_gene_is_worth_a_look():
    call = markers.call_botulinum({"bont_types": []}, "Clostridium_botulinum_groupII")
    assert "no bont gene detected" in "; ".join(call.reasons)


# --------------------------------------------------------------------------- #
#  Dispatch                                                                    #
# --------------------------------------------------------------------------- #
def test_dispatch_routes_name_keyed_rules_on_the_kalamari_name():
    call = markers.resolve("listeria_monocytogenes_serogroup_lineage_cc",
                           {"broad_lineage": "IV"}, "Listeria monocytogenes",
                           "Listeria_monocytogenes_IV")
    assert call.verdict == markers.VERDICT_CONFIRMED


def test_a_task_with_no_rule_yet_degrades_to_no_verdict():
    call = markers.resolve("some_future_task", {"anything": 1}, "X", "Y")
    assert call.verdict == markers.VERDICT_INDETERMINATE
    assert call.rule_version == "none"


def test_evidence_is_flattened_for_the_ledger():
    call = markers.call_shigella_eiec({"ipaH": "absent"}, "Escherichia coli")
    ev = call.as_evidence()
    assert ev["marker_verdict"] == markers.VERDICT_CONFIRMED
    assert ev["overrides_ani"] == "yes"
    assert ev["interpretation_rule_version"] == markers.RULE_SHIGELLA
