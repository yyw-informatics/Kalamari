# -*- coding: utf-8 -*-
"""Lead model: J1a lead building from CoverageRow, fingerprints, and record shape.

Uses real ``coverage.Pick`` / ``coverage.CoverageRow`` objects (no network) so the
reuse contract with the coverage layer is exercised, not mocked away."""

import coverage
import leads
from leads import Lead


def _pick(unit="Foo~GCF_1.1", name="Foo", parent="GCF_1.1", lane="NCBI-REF",
          species_taxid=100, species_name="Foo bar", regime="", resolver=""):
    return coverage.Pick(
        unit_id=unit, scientific_name=name, species_taxid=species_taxid,
        species_name=species_name, parent_assembly=parent, accession_set=["CP1"],
        source_lane=lane, regime=regime, resolver=resolver)


def _row(flags, **kw):
    p = _pick(**kw)
    r = coverage.CoverageRow(pick=p, flags=list(flags))
    return r


# --------------------------------------------------------------------------- #
#  j1a_leads: only actionable rows, evidence + link + summary populated         #
# --------------------------------------------------------------------------- #
def test_j1a_leads_only_actionable():
    rows = [
        _row(["tax_check_failed"], unit="A", name="Aa"),                 # actionable
        _row(["best_match_species_mismatch"], unit="B", name="Bb"),      # NOT actionable
        _row([], unit="C", name="Cc"),                                    # clean
        _row(["refseq_excluded_identity"], unit="D", name="Dd"),         # actionable
    ]
    ls = leads.j1a_leads(rows)
    assert {l.unit_id for l in ls} == {"A", "D"}
    assert all(l.signal == "J1a" for l in ls)


def test_j1a_lead_fields_and_link():
    r = _row(["tax_check_failed"], unit="Cit~GCF_9.1", name="Citrobacter_freundii",
             parent="GCF_9.1", lane="NCBI-REF")
    r.taxonomy_check = "Failed"
    r.best_match_species_name = "Citrobacter cronae"
    (lead,) = leads.j1a_leads([r])
    assert lead.lead_key == "Cit~GCF_9.1"
    assert lead.actionable_flags == ["tax_check_failed"]
    assert lead.link == "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_9.1/"
    assert lead.evidence["taxonomy_check"] == "Failed"
    assert "FAILED" in lead.summary


# --------------------------------------------------------------------------- #
#  Fingerprint: content hash of the actionable-flag SET (order-independent)      #
# --------------------------------------------------------------------------- #
def test_j1a_fingerprint_set_based():
    a = Lead(signal="J1a", lead_key="U", unit_id="U", scientific_name="U",
             source_lane="NCBI-REF", epoch="e", actionable_flags=["a", "b"])
    b = Lead(signal="J1a", lead_key="U", unit_id="U", scientific_name="U",
             source_lane="NCBI-REF", epoch="e", actionable_flags=["b", "a"])   # reordered
    c = Lead(signal="J1a", lead_key="U", unit_id="U", scientific_name="U",
             source_lane="NCBI-REF", epoch="e", actionable_flags=["a"])        # smaller set
    assert a.fingerprint == b.fingerprint          # order-independent
    assert a.fingerprint != c.fingerprint          # set changed -> fingerprint changed


def test_j1a_fingerprint_epoch_independent():
    a = Lead(signal="J1a", lead_key="U", unit_id="U", scientific_name="U",
             source_lane="NCBI-REF", epoch="e1", actionable_flags=["a"])
    b = Lead(signal="J1a", lead_key="U", unit_id="U", scientific_name="U",
             source_lane="NCBI-REF", epoch="e2", actionable_flags=["a"])
    assert a.fingerprint == b.fingerprint          # epoch is matched separately, not hashed


def test_j1b_fingerprint_tracks_candidate():
    a = Lead(signal="J1b", lead_key="U#GCA_9", unit_id="U", scientific_name="U",
             source_lane="NCBI-REF", epoch="e", candidate_accession="GCA_9")
    b = Lead(signal="J1b", lead_key="U#GCA_10", unit_id="U", scientific_name="U",
             source_lane="NCBI-REF", epoch="e", candidate_accession="GCA_10")
    assert a.fingerprint != b.fingerprint


# --------------------------------------------------------------------------- #
#  Record serialisation                                                        #
# --------------------------------------------------------------------------- #
def test_as_record_self_contained():
    lead = Lead(signal="J1a", lead_key="U", unit_id="U", scientific_name="Sp",
                source_lane="SME", epoch="e", actionable_flags=["tax_check_failed"],
                flags=["tax_check_failed", "best_match_species_mismatch"],
                deposit_date="2019/10/16", link="http://x", summary="why")
    rec = lead.as_record(run_id="2026-08", delta="new",
                         tool_versions={"ncbi_datasets_cli": "16.0.0"})
    assert rec["run_id"] == "2026-08"
    assert rec["delta"] == "new"
    assert rec["actionable_flags"] == ["tax_check_failed"]
    assert rec["evidence_fingerprint"] == lead.fingerprint
    assert rec["tool_versions"]["ncbi_datasets_cli"] == "16.0.0"
    assert rec["params"]["af_gate"] == 0.5      # stamped default
