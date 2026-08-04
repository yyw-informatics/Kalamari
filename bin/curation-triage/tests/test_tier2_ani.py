# -*- coding: utf-8 -*-
"""Tier 2: the ANI aligned-fraction gate and the type-strain confirmation."""

import pytest

import skani_runner
import tier2
from skani_runner import parse_skani


def res(ani=98.0, af_q=0.9, af_r=0.9, ref="TYPESTRAIN_1"):
    return tier2.AniResult(query="GCF_1.1", reference=ref, ani=ani,
                           af_query=af_q, af_reference=af_r)


# --------------------------------------------------------------------------- #
#  The gate                                                                    #
# --------------------------------------------------------------------------- #
def test_high_ani_over_good_alignment_is_same_species():
    assert tier2.ani_gate(res(ani=98.0)) == tier2.GATE_SAME


def test_low_ani_over_good_alignment_is_a_different_species():
    assert tier2.ani_gate(res(ani=91.0)) == tier2.GATE_DIFFERENT


def test_aligned_fraction_below_the_gate_vetoes_a_high_ani():
    """A high ANI over a sliver of the genome is not weak evidence, it is none."""
    assert tier2.ani_gate(res(ani=99.9, af_q=0.2, af_r=0.2)) == tier2.GATE_UNUSABLE


def test_the_weaker_aligned_fraction_decides():
    """Containment: the query aligns fully, but covers little of the reference."""
    r = res(ani=99.0, af_q=0.98, af_r=0.20)
    assert r.af == pytest.approx(0.20)
    assert tier2.ani_gate(r) == tier2.GATE_UNUSABLE


def test_ani_inside_the_epsilon_band_is_too_close_to_call():
    assert tier2.ani_gate(res(ani=95.1)) == tier2.GATE_BORDERLINE
    assert tier2.ani_gate(res(ani=94.9)) == tier2.GATE_BORDERLINE
    # Just outside the band, the call is made.
    assert tier2.ani_gate(res(ani=95.3)) == tier2.GATE_SAME
    assert tier2.ani_gate(res(ani=94.7)) == tier2.GATE_DIFFERENT


def test_a_missing_number_is_unusable_not_zero():
    assert tier2.ani_gate(tier2.AniResult(ani=None, af_query=0.9)) == tier2.GATE_UNUSABLE
    assert tier2.ani_gate(tier2.AniResult(ani=98.0, af_query=None,
                                          af_reference=None)) == tier2.GATE_UNUSABLE


def test_thresholds_are_overridable():
    assert tier2.ani_gate(res(ani=96.0), species_ani=97.0) == tier2.GATE_DIFFERENT
    assert tier2.ani_gate(res(ani=99.0, af_q=0.4, af_r=0.4), af_gate=0.3) == tier2.GATE_SAME


# --------------------------------------------------------------------------- #
#  confirm_species                                                             #
# --------------------------------------------------------------------------- #
def test_confirming_a_pick_against_its_type_strain():
    conf = tier2.confirm_species("u", "Mycobacterium leprae", res(ani=98.9, af_q=0.87))
    assert conf.verdict == tier2.VERDICT_CONFIRMED
    ev = conf.as_evidence()
    assert ev["t2_check"] == tier2.CHECK_ANI
    assert ev["ani"] == "98.90" and ev["aligned_fraction"] == "0.870"


def test_a_pick_below_the_species_boundary_is_refuted():
    conf = tier2.confirm_species("u", "Citrobacter freundii", res(ani=88.0))
    assert conf.verdict == tier2.VERDICT_REFUTED
    assert "below the 95.0% species boundary" in "; ".join(conf.reasons)


def test_an_unusable_alignment_gives_no_verdict_rather_than_a_wrong_one():
    conf = tier2.confirm_species("u", "Buchnera aphidicola", res(ani=96.5, af_q=0.31, af_r=0.29))
    assert conf.verdict == tier2.VERDICT_INDETERMINATE
    assert conf.gate == tier2.GATE_UNUSABLE


def test_no_type_strain_means_unresolved_never_a_guess():
    conf = tier2.confirm_species("u", "Whatever sp.", None)
    assert conf.verdict == tier2.VERDICT_INDETERMINATE
    assert "type strain" in "; ".join(conf.reasons)
    assert conf.as_evidence()["type_strain_assembly"] == "unresolved"


# --------------------------------------------------------------------------- #
#  skani output parsing                                                        #
# --------------------------------------------------------------------------- #
SKANI_OUT = (
    "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
    "/c/TYPESTRAIN_1396.fna\t/c/GCF_000283675.1.fna\t97.80\t84.10\t85.30\tref\tqry\n"
    "/c/GCF_000019305.1.fna\t/c/GCF_000283675.1.fna\t91.02\t70.00\t71.50\tref\tqry\n"
)


def test_skani_output_is_parsed_by_column_name_and_converted_to_fractions():
    rows = parse_skani(SKANI_OUT)
    assert [r.reference for r in rows] == ["TYPESTRAIN_1396", "GCF_000019305.1"]
    assert rows[0].query == "GCF_000283675.1"
    assert rows[0].ani == pytest.approx(97.80)
    # skani reports the aligned fraction as a percentage; the gate works in 0-1.
    assert rows[0].af_query == pytest.approx(0.853)
    assert rows[0].af == pytest.approx(0.841)


def test_a_renamed_skani_column_fails_loudly_instead_of_parsing_garbage():
    with pytest.raises(ValueError):
        parse_skani("Ref_file\tQuery_file\tIDENTITY\n/a.fna\t/b.fna\t97.8\n")


def test_empty_skani_output_is_empty_not_an_error():
    assert parse_skani("") == []


# --------------------------------------------------------------------------- #
#  The cache epoch                                                             #
# --------------------------------------------------------------------------- #
def test_a_failed_version_probe_never_overwrites_the_recorded_epoch():
    """A missing skani is a broken runner, not a new skani.

    Recording it as one would bump the epoch, which re-baselines every Tier-2 ANI
    lead and voids every dismissal in reviews.tsv.
    """
    cache = skani_runner.AniCache(runner=lambda q, refs: "")
    cache.meta = {"skani_version": "skani@0.3.2"}
    cache.build([], store=None, version_epoch="skani@unknown")
    assert cache.meta["skani_version"] == "skani@0.3.2"


def test_a_real_version_change_is_recorded():
    cache = skani_runner.AniCache(runner=lambda q, refs: "")
    cache.meta = {"skani_version": "skani@0.3.2"}
    cache.build([], store=None, version_epoch="skani@0.4.0")
    assert cache.meta["skani_version"] == "skani@0.4.0"
