# -*- coding: utf-8 -*-
"""The Tier 1 handoff loop: delta semantics, review suppression, partition,
ledger round-trip, and worklist render -- all pure and offline."""

import ledger
from leads import Lead, ANI_EPOCH, DATASETS_EPOCH


def j1a(unit, flags, epoch=ANI_EPOCH, name=None):
    return Lead(signal="J1a", lead_key=unit, unit_id=unit,
                scientific_name=name or unit, source_lane="NCBI-REF", epoch=epoch,
                actionable_flags=list(flags), flags=list(flags),
                summary="%s: %s" % (name or unit, ";".join(flags)))


def j1b(unit, candidate, epoch=DATASETS_EPOCH):
    return Lead(signal="J1b", lead_key="%s#%s" % (unit, candidate), unit_id=unit,
                scientific_name=unit, source_lane="NCBI-REF", epoch=epoch,
                candidate_accession=candidate)


def prior_from(leads, run_id="2026-07", delta="new"):
    """A prior ledger (list of records) as if these leads were recorded last run."""
    return [l.as_record(run_id, delta) for l in leads]


# --------------------------------------------------------------------------- #
#  Delta: J1a set-diff                                                         #
# --------------------------------------------------------------------------- #
def test_delta_new_when_no_prior():
    leads = [j1a("U1", ["tax_check_failed"])]
    d = ledger.compute_delta(leads, prior_records=[])
    assert d["U1"] == "new"


def test_delta_recurring_when_flag_set_unchanged():
    leads = [j1a("U1", ["tax_check_inconclusive"])]
    prior = prior_from(leads)
    assert ledger.compute_delta(leads, prior)["U1"] == "recurring"


def test_delta_new_when_flag_set_gains_member():
    prior = prior_from([j1a("U1", ["tax_check_inconclusive"])])
    now = [j1a("U1", ["tax_check_inconclusive", "refseq_excluded_identity"])]
    assert ledger.compute_delta(now, prior)["U1"] == "new"      # gained a member


def test_delta_recurring_when_flag_set_shrinks_but_still_actionable():
    prior = prior_from([j1a("U1", ["tax_check_inconclusive", "refseq_excluded_identity"])])
    now = [j1a("U1", ["tax_check_inconclusive"])]               # lost one, still actionable
    assert ledger.compute_delta(now, prior)["U1"] == "recurring"   # did not GAIN -> not new


def test_delta_epoch_bump_rebaselines_to_new():
    prior = prior_from([j1a("U1", ["tax_check_failed"], epoch="ncbi@2026-07")])
    now = [j1a("U1", ["tax_check_failed"], epoch="ncbi@2027-01")]   # same flags, new epoch
    assert ledger.compute_delta(now, prior)["U1"] == "new"


def test_delta_j1b_first_appearance():
    now = [j1b("U1", "GCA_999")]
    assert ledger.compute_delta(now, prior_records=[])["U1#GCA_999"] == "new"
    prior = prior_from(now)
    assert ledger.compute_delta(now, prior)["U1#GCA_999"] == "recurring"
    # a different candidate is a different key -> new again
    other = [j1b("U1", "GCA_1000")]
    assert ledger.compute_delta(other, prior)["U1#GCA_1000"] == "new"


# --------------------------------------------------------------------------- #
#  Suppression: "until evidence changes" + tolerant blank cells                 #
# --------------------------------------------------------------------------- #
def _review(lead, decision, fingerprint=None, epoch=None, **extra):
    r = {"lead_key": lead.lead_key, "decision": decision,
         "evidence_fingerprint": lead.fingerprint if fingerprint is None else fingerprint,
         "epoch": lead.epoch if epoch is None else epoch}
    r.update(extra)
    return r


def test_suppressed_while_fingerprint_and_epoch_match():
    lead = j1a("U1", ["tax_check_inconclusive"])
    reviews = [_review(lead, "reject")]
    assert ledger.review_decision(lead, reviews) == "reject"


def test_resurfaces_when_flag_set_changes():
    reviewed = j1a("U1", ["tax_check_inconclusive"])
    reviews = [_review(reviewed, "reject")]                     # fingerprint of the OLD set
    worse = j1a("U1", ["tax_check_inconclusive", "refseq_excluded_identity"])
    assert ledger.review_decision(worse, reviews) is None      # evidence changed -> re-surfaces


def test_blank_fingerprint_suppresses_until_epoch_bump():
    lead = j1a("U1", ["tax_check_inconclusive"])
    reviews = [_review(lead, "reject", fingerprint="")]        # blank fp = any evidence
    worse = j1a("U1", ["tax_check_inconclusive", "refseq_excluded_identity"])
    assert ledger.review_decision(worse, reviews) == "reject"  # still suppressed (same epoch)
    bumped = j1a("U1", ["tax_check_inconclusive"], epoch="ncbi@2099-01")
    assert ledger.review_decision(bumped, reviews) is None      # epoch bump busts it


def test_blank_fingerprint_and_epoch_suppresses_forever():
    lead = j1a("U1", ["tax_check_inconclusive"])
    reviews = [_review(lead, "reject", fingerprint="", epoch="")]
    forever = j1a("U1", ["x", "y"], epoch="whatever")
    assert ledger.review_decision(forever, reviews) == "reject"


def test_reject_overrides_accept():
    lead = j1a("U1", ["tax_check_failed"])
    reviews = [_review(lead, "accept"), _review(lead, "reject")]
    assert ledger.review_decision(lead, reviews) == "reject"


def test_unrelated_review_does_not_match():
    lead = j1a("U1", ["tax_check_failed"])
    reviews = [_review(j1a("OTHER", ["tax_check_failed"]), "reject")]
    assert ledger.review_decision(lead, reviews) is None


# --------------------------------------------------------------------------- #
#  Partition: buckets + resolved + records                                     #
# --------------------------------------------------------------------------- #
def test_partition_buckets():
    prior = prior_from([j1a("U_std", ["tax_check_inconclusive"])])   # standing
    leads = [
        j1a("U_new", ["tax_check_failed"]),                          # new
        j1a("U_std", ["tax_check_inconclusive"]),                    # recurring/standing
        j1a("U_acc", ["refseq_excluded_identity"]),                  # accepted
        j1a("U_rej", ["assembly_suppressed"]),                       # rejected
    ]
    reviews = [_review(leads[2], "accept", proposed_edit="add row X"),
               _review(leads[3], "reject")]
    part = ledger.partition(leads, prior, reviews, run_id="2026-08")
    assert [l.unit_id for l in part.new] == ["U_new"]
    assert [l.unit_id for l in part.standing] == ["U_std"]
    assert [l.unit_id for l in part.accepted] == ["U_acc"]
    assert [l.unit_id for l in part.dismissed] == ["U_rej"]
    # every current lead recorded (state log), regardless of review state
    assert len([r for r in part.records if r["delta"] != "resolved"]) == 4


def test_partition_resolved_transition():
    prior = prior_from([j1a("U1", ["tax_check_failed"]),
                        j1a("U2", ["tax_check_inconclusive"])])
    now = [j1a("U1", ["tax_check_failed"])]        # U2 no longer actionable
    part = ledger.partition(now, prior, reviews=[], run_id="2026-08")
    assert [r["lead_key"] for r in part.resolved] == ["U2"]
    assert part.resolved[0]["delta"] == "resolved"
    assert part.resolved[0]["actionable_flags"] == []


def test_resolved_not_re_emitted():
    # once a unit has a 'resolved' record, a later run must not re-resolve it.
    prior = [j1a("U2", ["tax_check_inconclusive"]).as_record("2026-07", "new"),
             {"schema": 1, "run_id": "2026-08", "signal": "J1a", "lead_key": "U2",
              "epoch": ANI_EPOCH, "delta": "resolved", "actionable_flags": []}]
    part = ledger.partition([], prior, reviews=[], run_id="2026-09")
    assert part.resolved == []


# --------------------------------------------------------------------------- #
#  Ledger round-trip + worklist render                                         #
# --------------------------------------------------------------------------- #
def test_ledger_append_and_load_roundtrip(tmp_path):
    path = str(tmp_path / "ledger.ndjson")
    recs = [j1a("U1", ["tax_check_failed"]).as_record("2026-08", "new")]
    ledger.append_ledger(path, recs)
    ledger.append_ledger(path, [j1a("U2", ["assembly_suppressed"]).as_record("2026-09", "new")])
    loaded = ledger.load_ledger(path)
    assert [r["lead_key"] for r in loaded] == ["U1", "U2"]
    assert ledger.latest_by_key(loaded)["U2"]["run_id"] == "2026-09"


def test_load_reviews_missing_file_is_empty(tmp_path):
    assert ledger.load_reviews(str(tmp_path / "nope.tsv")) == []


def test_write_worklist_and_markdown(tmp_path):
    prior = prior_from([j1a("U_std", ["tax_check_inconclusive"])])
    leads = [j1a("U_new", ["tax_check_failed"], name="New sp"),
             j1a("U_std", ["tax_check_inconclusive"], name="Std sp")]
    part = ledger.partition(leads, prior, reviews=[], run_id="2026-08")
    out = str(tmp_path / "worklist.tsv")
    ledger.write_worklist(part, out)
    with open(out) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        body = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]
    assert header == ledger.WORKLIST_COLUMNS
    statuses = {r["scientificName"]: r["status"] for r in body}
    assert statuses == {"New sp": "new", "Std sp": "standing"}

    md = ledger.render_markdown(part)
    assert "**1 new**" in md
    assert "New actionable leads" in md
    assert "New sp" in md
