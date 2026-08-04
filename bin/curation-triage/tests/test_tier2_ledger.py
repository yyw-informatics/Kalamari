# -*- coding: utf-8 -*-
"""Tier 2 attaching to the shared ledger: the override, the delta, the suppression.

Tier 2 is a new signal on the SAME ledger, so these tests guard the seam: it must
diff and suppress like every other signal, and it must not disturb Tier 1's state.
"""

import ledger as ledger_mod
import markers
import tier2
from coverage import Pick
from leads import Lead, SIGNAL_J1A, SIGNAL_T2


def pick(unit_id="Bacillus_cereus~GCF_1.1", regime="A", name="Bacillus_cereus",
         species="Bacillus cereus"):
    return Pick(unit_id=unit_id, scientific_name=name, species_taxid=1396,
                species_name=species, parent_assembly="GCF_1.1", accession_set=["CP1"],
                source_lane="SME", regime=regime, resolver="marker" if regime == "A" else "ANI")


def ani_conf(ani=97.8, verdict=tier2.VERDICT_CONFIRMED):
    res = tier2.AniResult(query="GCF_1.1", reference="TYPESTRAIN_1396", ani=ani,
                          af_query=0.85, af_reference=0.84)
    return tier2.confirm_species("u", "Bacillus cereus", res)


def marker_call(verdict=markers.VERDICT_REFUTED, call="Bacillus paranthracis"):
    return markers.MarkerCall(task_id="bacillus_anthracis_vs_b_cereus_group",
                              rule_version="bacillus_anthracis_vs_b_cereus_group@1",
                              call=call, verdict=verdict, confidence="high",
                              overrides_ani=True)


# --------------------------------------------------------------------------- #
#  The regime-A override                                                       #
# --------------------------------------------------------------------------- #
def test_a_regime_a_marker_call_decides_and_records_that_it_overrode_ani():
    res = tier2.resolve_unit("u", markers.REGIME_A, ani=ani_conf(),
                             marker_calls=[marker_call()])
    assert res.verdict == tier2.VERDICT_REFUTED       # the ANI said confirmed
    assert res.deciding_check == tier2.CHECK_MARKER
    assert res.overrode_ani is True


def test_an_indeterminate_regime_a_marker_does_not_hand_the_decision_back_to_ani():
    """For these taxa the ANI answer was never admissible, so 'undecided' is honest."""
    res = tier2.resolve_unit("u", markers.REGIME_A, ani=ani_conf(),
                             marker_calls=[marker_call(verdict=markers.VERDICT_INDETERMINATE)])
    assert res.verdict == tier2.VERDICT_INDETERMINATE


def test_regime_b_lets_the_panel_decide_membership():
    place = tier2.PanelPlacement(unit_id="u", verdict=tier2.VERDICT_REFUTED,
                                 not_distinct=["sibling"])
    res = tier2.resolve_unit("u", markers.REGIME_B, ani=ani_conf(), panel=place,
                             marker_calls=[marker_call(verdict=markers.VERDICT_CONFIRMED,
                                                       call="ok")])
    assert res.verdict == tier2.VERDICT_REFUTED
    assert res.overrode_ani is False


def test_a_wrong_regime_b_label_is_a_defect_even_on_a_well_placed_genome():
    place = tier2.PanelPlacement(unit_id="u", verdict=tier2.VERDICT_CONFIRMED)
    res = tier2.resolve_unit("u", markers.REGIME_B, panel=place,
                             marker_calls=[marker_call(verdict=markers.VERDICT_REFUTED)])
    assert res.verdict == tier2.VERDICT_REFUTED
    assert res.deciding_check == tier2.CHECK_MARKER


def test_an_ordinary_unit_is_decided_by_ani():
    res = tier2.resolve_unit("u", "", ani=ani_conf())
    assert res.deciding_check == tier2.CHECK_ANI
    assert res.verdict == tier2.VERDICT_CONFIRMED


# --------------------------------------------------------------------------- #
#  Leads                                                                       #
# --------------------------------------------------------------------------- #
def pins():
    return {"task_id": "bacillus_anthracis_vs_b_cereus_group", "regime": "A",
            "interpretation_rule_version": "bacillus_anthracis_vs_b_cereus_group@1",
            "components": [{"component": "group_taxonomy_and_marker_profile",
                            "method": {"name": "BTyper3", "conda_package": "btyper3",
                                       "software_version": "3.4.0-0"},
                            "database": {"version_or_scheme_id": "bundled-with-3.4.0-0"},
                            "status": "pinned"}]}


def test_a_marker_epoch_pins_software_database_and_rule_together():
    epoch = tier2.marker_epoch(pins())
    assert epoch == ("btyper3@3.4.0-0+db:bundled-with-3.4.0-0"
                     "+rule:bacillus_anthracis_vs_b_cereus_group@1")


def test_for_regime_a_the_ani_number_rides_along_as_context_only():
    lead = tier2.marker_lead(pick(), marker_call(), "e", pins(), ani=ani_conf(),
                             confirms_lead="Bacillus_cereus~GCF_1.1")
    assert "context only" in lead.evidence["ani_context"]
    assert lead.evidence["verdict"] == tier2.VERDICT_REFUTED
    assert lead.evidence["confirms_lead"] == "Bacillus_cereus~GCF_1.1"
    assert lead.lead_key == "t2:Bacillus_cereus~GCF_1.1#marker:bacillus_anthracis_vs_b_cereus_group"


def test_tier2_lead_keys_live_in_their_own_namespace():
    p = pick()
    keys = {tier2.marker_lead(p, marker_call(), "e", pins()).lead_key,
            tier2.ani_lead(p, ani_conf(), "e").lead_key}
    assert all(k.startswith("t2:") for k in keys)
    assert len(keys) == 2          # one unit can carry several checks without colliding


def test_background_context_rides_along_without_changing_the_fingerprint():
    """ECTyper's serotype is context: a new serotype must not re-open a dismissed lead."""
    p = pick()
    plain = tier2.marker_lead(p, marker_call(), "e", pins())
    with_context = tier2.marker_lead(p, marker_call(), "e", pins(),
                                     context={"ectyper_serotype": "O118/O151:H16",
                                              "ectyper_role": "background context only"})
    assert with_context.evidence["ectyper_serotype"] == "O118/O151:H16"
    assert with_context.fingerprint == plain.fingerprint


def test_a_context_field_can_never_overwrite_a_verdict_field():
    """Context is written first, so the marker call always wins the shared keys."""
    lead = tier2.marker_lead(pick(), marker_call(), "e", pins(),
                             context={"verdict": tier2.VERDICT_CONFIRMED,
                                      "marker_call": "something else"})
    assert lead.evidence["verdict"] == tier2.VERDICT_REFUTED
    assert lead.evidence["marker_call"] == "Bacillus paranthracis"


def test_confirmations_are_recorded_but_are_not_work():
    p = pick()
    needs, confirms = tier2.split_by_verdict([
        tier2.marker_lead(p, marker_call(), "e", pins()),
        tier2.ani_lead(p, ani_conf(), "e"),
    ])
    assert [l.evidence["verdict"] for l in needs] == [tier2.VERDICT_REFUTED]
    assert [l.evidence["verdict"] for l in confirms] == [tier2.VERDICT_CONFIRMED]


# --------------------------------------------------------------------------- #
#  Delta + suppression on the shared ledger                                    #
# --------------------------------------------------------------------------- #
def t2_lead(verdict=tier2.VERDICT_CONFIRMED, call="Bacillus cereus", epoch="e1"):
    return Lead(signal=SIGNAL_T2, lead_key="t2:u#marker:t", unit_id="u",
                scientific_name="U", source_lane="SME", epoch=epoch,
                actionable_flags=["t2_marker_confirms_label"],
                evidence={"t2_check": "marker", "verdict": verdict, "call": call})


def test_a_repeated_identical_verdict_is_standing_not_new():
    lead = t2_lead()
    prior = [lead.as_record("2026-07", "new")]
    assert ledger_mod.compute_delta([lead], prior) == {lead.lead_key: "recurring"}


def test_a_changed_verdict_surfaces_as_new_work():
    prior = [t2_lead().as_record("2026-07", "new")]
    flipped = t2_lead(verdict=tier2.VERDICT_REFUTED, call="Bacillus paranthracis")
    assert ledger_mod.compute_delta([flipped], prior) == {flipped.lead_key: "new"}


def test_the_fingerprint_ignores_ani_jitter_but_not_the_verdict():
    """A skani bump moves the last decimal of every number; that must not bust reviews."""
    a = t2_lead()
    b = t2_lead()
    a.evidence["ani"] = "97.80"
    b.evidence["ani"] = "97.81"
    assert a.fingerprint == b.fingerprint
    b.evidence["verdict"] = tier2.VERDICT_REFUTED
    assert a.fingerprint != b.fingerprint


def test_an_epoch_bump_reopens_a_dismissed_tier2_lead():
    """A database refresh is a tooling epoch: the call must be re-reviewed."""
    lead = t2_lead()
    reviews = [{"lead_key": lead.lead_key, "decision": "reject",
                "evidence_fingerprint": lead.fingerprint, "epoch": "e1"}]
    assert ledger_mod.review_decision(lead, reviews) == "reject"
    bumped = t2_lead(epoch="e2")
    assert ledger_mod.review_decision(bumped, reviews) is None


def test_tier2_never_marks_a_tier1_lead_resolved():
    """A Tier-2 run says nothing about whether a Tier-1 lead is still actionable."""
    j1a = Lead(signal=SIGNAL_J1A, lead_key="unit_a", unit_id="unit_a", scientific_name="A",
               source_lane="NCBI-REF", epoch="ncbi@1", actionable_flags=["tax_check_failed"])
    prior = [j1a.as_record("2026-07", "new")]
    part = ledger_mod.partition([t2_lead()], prior, [], run_id="2026-08",
                                resolve_signals=())
    assert part.resolved == []
    assert [r["lead_key"] for r in part.records] == ["t2:u#marker:t"]


def test_tier1_still_resolves_its_own_leads_by_default():
    j1a = Lead(signal=SIGNAL_J1A, lead_key="unit_a", unit_id="unit_a", scientific_name="A",
               source_lane="NCBI-REF", epoch="ncbi@1", actionable_flags=["tax_check_failed"])
    prior = [j1a.as_record("2026-07", "new")]
    part = ledger_mod.partition([], prior, [], run_id="2026-08")
    assert [r["delta"] for r in part.resolved] == ["resolved"]
