#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Tier 1 handoff loop: ledger (append-only NDJSON) + reviews (SME decisions) +
the monthly delta and the worklist render (pure & offline).

This is the piece the plan specifies only in prose, so the schemas and semantics
live here, documented, and everything is a pure function over in-memory records so
the delta/suppression logic is unit-testable with zero network and zero real files.

Delta semantics ("what is *newly* actionable", locked design decision)
----------------------------------------------------------------------
Against the most recent prior ledger record for a lead:
  * **J1a** surfaces as ``new`` when its actionable-flag SET *gains* a member
    (``current - prior`` is non-empty); an unchanged set is ``recurring`` (standing),
    and a set that shrank to empty is ``resolved`` (an informational note, not work).
  * **J1b / J2** surface as ``new`` on the first appearance of the candidate (the
    ``lead_key`` includes the candidate accession, so a fresh candidate is a fresh key).
  * On a signal's **epoch bump** the prior is treated as empty (re-baselined), so a
    genuine NCBI recompute / tool-version change re-opens the standing set for review.

Suppression semantics ("until evidence changes", locked design decision)
-------------------------------------------------------------------------
A ``reviews.tsv`` row suppresses a lead from the *new* worklist while BOTH hold:
  * its ``evidence_fingerprint`` matches the lead's (a changed flag set / candidate
    busts it), AND
  * its ``epoch`` matches the lead's (an epoch bump busts it).
The two columns are *tolerant*: leaving ``evidence_fingerprint`` blank suppresses the
lead until the next epoch bump regardless of evidence; leaving ``epoch`` blank too
suppresses it forever.  So one mechanism gives the SME per-lead choice of
until-evidence-changes (fill both -- the recommended default), until-epoch (blank
fingerprint), or forever (blank both).  ``accept`` routes the lead to the
proposed-edit list instead of dismissing it; ``reject`` / ``defer`` dismiss it.

Ledger record = one self-contained NDJSON line per *surfaced or resolved* lead per
run (NOT one per taxon per run -- a per-taxon log would bury the signal under ~3000
no-op rows a year).  Every actionable lead of the run is recorded (regardless of
review state) so the next run always has a diff baseline; suppression only affects
what the worklist *shows*.
"""

from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from leads import Lead, SIGNAL_J1A, SIGNAL_T2, _fp

REVIEW_DECISIONS = {"accept", "reject", "defer"}

# reviews.tsv schema.  ``lead_key`` + ``decision`` are the only required cells; the
# rest scope the dismissal and record who/why.  See module docstring for the tolerant
# fingerprint/epoch matching that gives per-lead until-evidence / until-epoch / forever.
REVIEW_COLUMNS = [
    "lead_key", "decision", "evidence_fingerprint", "epoch",
    "reviewer", "date", "note", "proposed_edit",
]

# The per-run worklist artifact (the surviving actionable rows).
WORKLIST_COLUMNS = [
    "run_id", "status", "signal", "lead_key", "unit_id", "scientificName",
    "source_lane", "actionable_flags", "candidate_accession", "deposit_date",
    "summary", "link", "evidence_fingerprint", "epoch",
]


# --------------------------------------------------------------------------- #
#  Ledger I/O (append-only NDJSON)                                             #
# --------------------------------------------------------------------------- #
def load_ledger(path: str) -> List[dict]:
    """Read the append-only ledger; a missing file is an empty history."""
    if not path or not os.path.exists(path):
        return []
    out: List[dict] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            out.append(json.loads(line))
    return out


def append_ledger(path: str, records: Sequence[dict]) -> None:
    """Append records as NDJSON (one compact JSON object per line)."""
    if not records:
        return
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "a", encoding="utf-8") as fh:
        for rec in records:
            fh.write(json.dumps(rec, sort_keys=True, ensure_ascii=False) + "\n")


def latest_by_key(records: Sequence[dict]) -> Dict[str, dict]:
    """The most recent ledger record per ``lead_key`` (file order = chronological)."""
    out: Dict[str, dict] = {}
    for rec in records:
        key = rec.get("lead_key")
        if key:
            out[key] = rec       # later append wins
    return out


# --------------------------------------------------------------------------- #
#  Reviews I/O                                                                 #
# --------------------------------------------------------------------------- #
def load_reviews(path: str) -> List[dict]:
    """Read reviews.tsv; a missing file means no decisions yet."""
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def review_decision(lead: Lead, reviews: Sequence[dict]) -> Optional[str]:
    """The SME decision that currently suppresses ``lead``, or None.

    Tolerant match (see module docstring): a review row applies when its ``lead_key``
    matches and -- if present -- its ``evidence_fingerprint`` and ``epoch`` match the
    lead's.  A blank cell means "don't constrain on this axis".  If several rows match,
    the most decisive wins: ``reject`` > ``defer`` > ``accept`` (a hard no is never
    overridden by a softer earlier decision).
    """
    found: List[str] = []
    for r in reviews:
        if (r.get("lead_key") or "").strip() != lead.lead_key:
            continue
        decision = (r.get("decision") or "").strip().lower()
        if decision not in REVIEW_DECISIONS:
            continue
        fp = (r.get("evidence_fingerprint") or "").strip()
        if fp and fp != lead.fingerprint:
            continue
        ep = (r.get("epoch") or "").strip()
        if ep and ep != lead.epoch:
            continue
        found.append(decision)
    for pref in ("reject", "defer", "accept"):
        if pref in found:
            return pref
    return None


# --------------------------------------------------------------------------- #
#  Delta                                                                       #
# --------------------------------------------------------------------------- #
def compute_delta(leads: Sequence[Lead], prior_records: Sequence[dict]) -> Dict[str, str]:
    """Map each lead's ``lead_key`` -> 'new' | 'recurring' (see module docstring).

    A lead whose signal epoch differs from its prior record's epoch is re-baselined
    (prior treated as empty) so an epoch bump re-opens it as ``new``.
    """
    prior = latest_by_key(prior_records)
    out: Dict[str, str] = {}
    for lead in leads:
        pr = prior.get(lead.lead_key)
        same_epoch = pr is not None and pr.get("epoch") == lead.epoch
        if lead.signal == SIGNAL_J1A:
            prior_set = set(pr.get("actionable_flags") or []) if same_epoch else set()
            gained = set(lead.actionable_flags) - prior_set
            out[lead.lead_key] = "new" if gained else "recurring"
        elif lead.signal == SIGNAL_T2:
            # Tier 2 keys are stable per (unit, check), so the delta is "did the
            # VERDICT change?" -- a flip from confirmed to refuted is the whole point
            # of re-running, and must surface as new work rather than as a standing row.
            prior_fp = pr.get("evidence_fingerprint") if same_epoch else None
            out[lead.lead_key] = "recurring" if prior_fp == lead.fingerprint else "new"
        else:
            # J1b / J2: first appearance of this candidate (lead_key includes it).
            out[lead.lead_key] = "recurring" if same_epoch else "new"
    return out


def _resolved_records(leads: Sequence[Lead], prior_records: Sequence[dict],
                      run_id: str) -> List[dict]:
    """One 'resolved' record for each J1a unit that WAS actionable and no longer is.

    Emitted once on the transition (the resolved record itself has empty flags, so it
    is not re-emitted on subsequent runs).

    Only ever call this when the run actually recomputed J1a over the whole pick set.
    "Absent from this run's leads" means "resolved" only if the run looked; a run that
    never computed J1a (Tier 2) or saw only part of it (a crashed shard) would
    otherwise mark every standing lead resolved and destroy the delta baseline.
    ``partition(resolve_signals=...)`` is the switch that enforces this.
    """
    prior = latest_by_key(prior_records)
    current_j1a = {lead.lead_key for lead in leads if lead.signal == SIGNAL_J1A}
    resolved: List[dict] = []
    for key, pr in prior.items():
        if pr.get("signal") != SIGNAL_J1A:
            continue
        if pr.get("delta") == "resolved" or not pr.get("actionable_flags"):
            continue                     # already resolved / never actionable
        if key in current_j1a:
            continue                     # still actionable
        resolved.append({
            "schema": 1,
            "run_id": run_id,
            "signal": SIGNAL_J1A,
            "lead_key": key,
            "epoch": pr.get("epoch", ""),
            "delta": "resolved",
            "unit_id": pr.get("unit_id", ""),
            "scientificName": pr.get("scientificName", ""),
            "source_lane": pr.get("source_lane", ""),
            "candidate_accession": "",
            "deposit_date": pr.get("deposit_date", ""),
            "actionable_flags": [],
            "flags": [],
            "evidence": {},
            "evidence_fingerprint": _fp("J1a|%s|" % key),
            "link": pr.get("link", ""),
            "summary": "%s: no longer actionable (was %s)" % (
                pr.get("scientificName", key), ";".join(pr.get("actionable_flags") or [])),
            "params": pr.get("params", {}),
            "tool_versions": pr.get("tool_versions", {}),
        })
    return resolved


# --------------------------------------------------------------------------- #
#  Partition: turn this run's leads + history + reviews into the worklist       #
# --------------------------------------------------------------------------- #
@dataclass
class Partition:
    run_id: str
    new: List[Lead] = field(default_factory=list)         # newly actionable, undismissed
    standing: List[Lead] = field(default_factory=list)    # still actionable, not new
    accepted: List[Lead] = field(default_factory=list)    # SME accepted -> proposed edit
    dismissed: List[Lead] = field(default_factory=list)   # reject/defer -> suppressed
    resolved: List[dict] = field(default_factory=list)    # resolved-this-run records
    records: List[dict] = field(default_factory=list)     # everything to append this run
    delta: Dict[str, str] = field(default_factory=dict)


def partition(leads: Sequence[Lead], prior_records: Sequence[dict],
              reviews: Sequence[dict], run_id: str,
              params: Dict = None, tool_versions: Dict = None,
              resolve_signals: Sequence[str] = (SIGNAL_J1A,)) -> Partition:
    """Diff ``leads`` against history + reviews into worklist buckets + ledger records.

    Every current lead is recorded (state log); the *new* bucket is the SME's work
    (newly actionable and not dismissed).  ``accept`` routes to ``accepted``;
    ``reject``/``defer`` to ``dismissed``.

    ``resolve_signals`` names the signals this run recomputed in full, and so is
    entitled to mark resolved.  It defaults to J1a (what the Tier-1 driver and
    collate always compute over the whole pick set).  A run that computes a
    different signal -- Tier 2 -- must pass ``()``: its leads say nothing about
    whether a Tier-1 lead is still actionable, and treating their absence as
    "resolved" would silently wipe the Tier-1 baseline.
    """
    delta = compute_delta(leads, prior_records)
    part = Partition(run_id=run_id, delta=delta)
    for lead in leads:
        d = delta[lead.lead_key]
        part.records.append(lead.as_record(run_id, d, params, tool_versions))
        decision = review_decision(lead, reviews)
        if decision == "accept":
            part.accepted.append(lead)
        elif decision in ("reject", "defer"):
            part.dismissed.append(lead)
        elif d == "new":
            part.new.append(lead)
        else:
            part.standing.append(lead)
    if SIGNAL_J1A in (resolve_signals or ()):
        part.resolved = _resolved_records(leads, prior_records, run_id)
        part.records.extend(part.resolved)
    return part


def summarize_partition(part: Partition) -> dict:
    return {
        "run_id": part.run_id,
        "new": len(part.new),
        "standing": len(part.standing),
        "accepted": len(part.accepted),
        "dismissed": len(part.dismissed),
        "resolved": len(part.resolved),
        "records_appended": len(part.records),
    }


# --------------------------------------------------------------------------- #
#  Render: worklist artifact (TSV) + $GITHUB_STEP_SUMMARY markdown              #
# --------------------------------------------------------------------------- #
def _worklist_row(lead: Lead, run_id: str, status: str) -> Dict[str, str]:
    return {
        "run_id": run_id,
        "status": status,
        "signal": lead.signal,
        "lead_key": lead.lead_key,
        "unit_id": lead.unit_id,
        "scientificName": lead.scientific_name,
        "source_lane": lead.source_lane,
        "actionable_flags": ";".join(lead.actionable_flags),
        "candidate_accession": lead.candidate_accession,
        "deposit_date": lead.deposit_date,
        "summary": lead.summary,
        "link": lead.link,
        "evidence_fingerprint": lead.fingerprint,
        "epoch": lead.epoch,
    }


def worklist_rows(part: Partition) -> List[Dict[str, str]]:
    """The surviving actionable rows (new + standing + accepted), most-urgent first."""
    rows: List[Dict[str, str]] = []
    for lead in part.new:
        rows.append(_worklist_row(lead, part.run_id, "new"))
    for lead in part.accepted:
        rows.append(_worklist_row(lead, part.run_id, "accepted"))
    for lead in part.standing:
        rows.append(_worklist_row(lead, part.run_id, "standing"))
    return rows


def write_worklist(part: Partition, path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=WORKLIST_COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in worklist_rows(part):
            writer.writerow(row)


def _md_cell(s: str) -> str:
    return (s or "").replace("|", "\\|").replace("\n", " ")


def _lead_table(leads: Sequence[Lead]) -> List[str]:
    out = ["| species | lane | signal | flags / candidate | why | link |",
           "| --- | --- | --- | --- | --- | --- |"]
    for lead in leads:
        detail = ";".join(lead.actionable_flags) or lead.candidate_accession
        out.append("| %s | %s | %s | %s | %s | %s |" % (
            _md_cell(lead.scientific_name), _md_cell(lead.source_lane),
            lead.signal, _md_cell(detail), _md_cell(lead.summary),
            ("[genome](%s)" % lead.link) if lead.link else ""))
    return out


def render_markdown(part: Partition, title: str = None) -> str:
    """A $GITHUB_STEP_SUMMARY-style markdown block; each row self-contained.

    ``title`` overrides the heading.  The default names the Tier-1 run, which is right
    for the monthly job; a caller rendering a different view -- the release snapshot
    built from the committed ledger, say -- passes its own so the heading does not
    claim to be something it is not.
    """
    s = summarize_partition(part)
    lines: List[str] = []
    lines.append(title or
                 "## Kalamari curation-triage — Tier 1 worklist (run %s)" % part.run_id)
    lines.append("")
    lines.append(
        "**%d new** actionable lead(s) · %d standing · %d accepted (pending edit) · "
        "%d resolved · %d dismissed" % (
            s["new"], s["standing"], s["accepted"], s["resolved"], s["dismissed"]))
    lines.append("")
    lines.append("_Informational triage, never a pass/fail gate. "
                 "Record decisions in `reviews.tsv` (accept / reject / defer)._")
    lines.append("")

    lines.append("### New actionable leads")
    if part.new:
        lines.extend(_lead_table(part.new))
    else:
        lines.append("_None this run._")
    lines.append("")

    if part.accepted:
        lines.append("### Accepted → proposed edits")
        lines.extend(_lead_table(part.accepted))
        lines.append("")

    if part.resolved:
        lines.append("### Resolved since last run")
        for rec in part.resolved:
            lines.append("- %s" % _md_cell(rec.get("summary", rec.get("lead_key", ""))))
        lines.append("")

    if part.standing:
        lines.append("<details><summary>Standing (%d) — already surfaced, not yet "
                     "actioned</summary>" % len(part.standing))
        lines.append("")
        lines.extend(_lead_table(part.standing))
        lines.append("")
        lines.append("</details>")
        lines.append("")

    return "\n".join(lines)
