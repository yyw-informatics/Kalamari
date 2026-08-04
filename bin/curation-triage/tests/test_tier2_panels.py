# -*- coding: utf-8 -*-
"""Tier 2: the sub-species panels -- generation, loading, and placement."""

import os

import build_panels
import tier2

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")
POLICY = os.path.join(SRC, "taxon_policy.tsv")
PANEL_DIR = os.path.join(SRC, "panels")


# --------------------------------------------------------------------------- #
#  The committed panels                                                        #
# --------------------------------------------------------------------------- #
def test_the_committed_panels_match_what_build_panels_generates():
    """The panels are generated, so a policy-table edit must not silently drift."""
    generated = build_panels.build_panels(POLICY)
    on_disk = sorted(f for f in os.listdir(PANEL_DIR) if f.endswith(".acclist"))
    assert on_disk == sorted(generated)
    for name, text in generated.items():
        with open(os.path.join(PANEL_DIR, name), encoding="utf-8") as fh:
            assert fh.read() == text, "%s is stale: rerun build_panels.py --write" % name


def test_every_regime_b_unit_has_a_panel_and_no_regime_a_unit_does():
    """Regime A is exactly where similarity is the wrong ruler, so it gets markers."""
    groups = build_panels.load_panel_groups(POLICY)
    members = [m for ms in groups.values() for m in ms]
    assert len(members) == 17, "expected the 17 regime-B confounded units"
    assert {m["regime"] for m in members} == {"B"}
    for m in members:
        assert os.path.exists(tier2.panel_path(PANEL_DIR, m["kalamari_taxid"]))


def test_a_panel_lists_the_siblings_and_never_the_unit_itself():
    listeria_iii = tier2.load_panel(tier2.panel_path(PANEL_DIR, 9000002))
    labels = {m.lineage_label for m in listeria_iii}
    assert labels == {"Listeria_monocytogenes_I", "Listeria_monocytogenes_II",
                      "Listeria_monocytogenes_IV"}
    assert all(m.accession.startswith(("GCF_", "GCA_")) for m in listeria_iii)
    assert all(m.role == "sibling" for m in listeria_iii)


def test_the_synthetic_lineages_all_have_a_panel():
    """The 13 metadata-blind units are the reason panels exist."""
    synthetic = [f for f in os.listdir(PANEL_DIR) if f.startswith("90000")]
    assert len(synthetic) == 13


def test_comment_lines_and_the_header_are_not_read_as_members():
    members = tier2.load_panel(tier2.panel_path(PANEL_DIR, 9000005))
    assert len(members) == 1                       # C. botulinum has one sibling
    assert members[0].lineage_label == "Clostridium_botulinum_groupII"


def test_a_missing_panel_file_is_an_empty_panel_not_a_crash():
    assert tier2.load_panel(tier2.panel_path(PANEL_DIR, 4242424)) == []


# --------------------------------------------------------------------------- #
#  Placement                                                                   #
# --------------------------------------------------------------------------- #
def r(ref, ani, af=0.9, label=""):
    return tier2.AniResult(query="GCF_Q.1", reference=ref, ani=ani, af_query=af,
                           af_reference=af, reference_label=label or ref)


def test_a_unit_distinct_from_all_its_siblings_is_consistent():
    place = tier2.place_against_panel(
        "u", "Listeria_monocytogenes_III",
        [r("A", 98.9, label="lin_I"), r("B", 98.7, label="lin_II"), r("C", 98.6, label="lin_IV")])
    assert place.verdict == tier2.VERDICT_CONFIRMED
    assert place.nearest_label == "lin_I"
    assert place.margin == 0.2


def test_a_sibling_above_the_derep_threshold_is_a_curation_defect():
    """Two Kalamari lineage references that are the same genome is a real finding."""
    place = tier2.place_against_panel("u", "Salmonella_enterica_IIIa",
                                      [r("X", 99.7, label="Salmonella_enterica_IIIb")])
    assert place.verdict == tier2.VERDICT_REFUTED
    assert place.not_distinct == ["Salmonella_enterica_IIIb"]
    assert "not distinct" in "; ".join(place.reasons)


def test_the_derep_threshold_is_inclusive_and_overridable():
    assert tier2.place_against_panel("u", "d", [r("X", 99.5)]).verdict == \
        tier2.VERDICT_REFUTED
    assert tier2.place_against_panel("u", "d", [r("X", 99.4)]).verdict == \
        tier2.VERDICT_CONFIRMED
    assert tier2.place_against_panel("u", "d", [r("X", 99.4)], derep=99.0).verdict == \
        tier2.VERDICT_REFUTED


def test_panel_comparisons_below_the_aligned_fraction_gate_are_dropped():
    place = tier2.place_against_panel("u", "d", [r("X", 99.9, af=0.2), r("Y", 98.0)])
    assert place.nearest_accession == "Y"
    assert "dropped" in "; ".join(place.reasons)


def test_no_usable_comparison_gives_no_verdict():
    place = tier2.place_against_panel("u", "d", [r("X", 99.9, af=0.1)])
    assert place.verdict == tier2.VERDICT_INDETERMINATE
    assert place.nearest_label == ""


def test_a_self_panel_reports_the_neighbourhood_but_does_not_assign_a_lineage():
    """The unit's own lineage is not in a sibling panel, so 'nearest' is not 'is'."""
    place = tier2.place_against_panel("u", "lineage_III", [r("A", 98.9, label="lineage_I")])
    assert place.verdict == tier2.VERDICT_CONFIRMED     # distinct, which is the question asked
    assert place.nearest_label == "lineage_I"           # not treated as a mismatch


def test_an_external_reference_panel_does_assign_a_lineage():
    place = tier2.place_against_panel(
        "u", "lineage_III", [r("A", 99.0, label="lineage_I"), r("B", 97.0, label="lineage_III")],
        self_panel=False)
    assert place.verdict == tier2.VERDICT_REFUTED
    assert place.nearest_label == "lineage_I"


def test_two_equally_close_external_references_are_ambiguous_not_a_coin_flip():
    place = tier2.place_against_panel(
        "u", "lineage_I", [r("A", 99.00, label="lineage_I"), r("B", 98.95, label="lineage_II")],
        self_panel=False)
    assert place.verdict == tier2.VERDICT_INDETERMINATE
    assert place.ambiguous
