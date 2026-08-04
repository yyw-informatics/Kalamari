#!/usr/bin/env python3
"""Render a worklist from the committed ledger, for the release asset.

Why this exists.  The monthly job writes ``worklist.tsv``/``worklist.md`` as *per-run
output* and never commits them (.gitignore), so at a ``v*.*.*`` tag there is no worklist
file in the checkout to attach to the GitHub Release.  Two ways to fix that are wrong:

* Committing the worklist every month duplicates a *view* of the ledger next to the
  ledger itself, and the two can drift.
* Re-running the pipeline at tag time recomputes the delta against a ledger that
  ALREADY contains this run's records, so every lead comes back ``recurring`` and the
  release document would announce "no new leads" for a month that had plenty.

So this reads the ledger and renders it.  The ledger already records each lead's delta
(``new`` / ``recurring`` / ``resolved``) as the run computed it, so the labels are the
ones a curator actually saw -- nothing is recomputed and nothing can drift.

Two views:

* ``--run-id 2026-08`` -- exactly what that run recorded.
* default (no ``--run-id``) -- the **snapshot**: the latest record per lead, across every
  signal.  This is the right view for a release, because Tier 1 and Tier 2 are appended
  by different jobs and can carry different run ids; a snapshot shows the state of the
  worklist as of the tag rather than the output of whichever tier ran last.

``reviews.tsv`` is applied the same way the live pipeline applies it, so a lead an SME
dismissed does not reappear on the release asset.

Offline, stdlib only.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import Dict, List, Optional, Sequence

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ledger as ledger_mod                                        # noqa: E402
from leads import Lead                                             # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")
DEFAULT_LEDGER = os.path.join(SRC, "ledger.ndjson")
DEFAULT_REVIEWS = os.path.join(SRC, "reviews.tsv")

# The ledger record spells the species field the JSON way; the Lead dataclass spells it
# the Python way.  Lead.from_dict() silently drops keys it does not know, so without this
# mapping every rendered row would have a blank species name.
RECORD_TO_LEAD = {"scientificName": "scientific_name"}


def lead_from_record(rec: dict) -> Lead:
    """Rebuild the Lead a ledger record was serialised from."""
    d = {RECORD_TO_LEAD.get(k, k): v for k, v in rec.items()}
    return Lead.from_dict(d)


def select_records(records: Sequence[dict], run_id: Optional[str] = None) -> List[dict]:
    """The records to render: one run, or the latest record per lead (the snapshot)."""
    if run_id:
        return [r for r in records if r.get("run_id") == run_id]
    # latest_by_key keeps the LAST record per lead_key in file order, which is exactly
    # "the current state" for an append-only log.
    latest = ledger_mod.latest_by_key(records)
    return [latest[k] for k in sorted(latest)]


def partition_from_ledger(records: Sequence[dict], reviews: Sequence[dict],
                          run_id: Optional[str] = None) -> ledger_mod.Partition:
    """Bucket committed records into the same Partition the live pipeline produces.

    The bucket comes from the record's OWN recorded delta -- it is never recomputed --
    except that a review decision still wins, so a dismissal recorded after the run is
    honoured here too.  ``part.records`` stays empty: rendering appends nothing.
    """
    selected = select_records(records, run_id)
    label = run_id or _latest_run_id(selected)
    part = ledger_mod.Partition(run_id=label)
    for rec in selected:
        delta = rec.get("delta", "")
        if delta == "resolved":
            part.resolved.append(rec)
            continue
        lead = lead_from_record(rec)
        part.delta[lead.lead_key] = delta
        decision = ledger_mod.review_decision(lead, reviews)
        if decision == "accept":
            part.accepted.append(lead)
        elif decision in ("reject", "defer"):
            part.dismissed.append(lead)
        elif delta == "new":
            part.new.append(lead)
        else:
            part.standing.append(lead)
    return part


def _latest_run_id(records: Sequence[dict]) -> str:
    """The newest run id present.  Run ids are YYYY-MM, so lexical order is date order."""
    return max((r.get("run_id", "") for r in records), default="")


def _run_id_of(records: Sequence[dict]) -> Dict[str, str]:
    """lead_key -> the run id of that lead's rendered record (for per-row provenance)."""
    return {r.get("lead_key", ""): r.get("run_id", "") for r in records}


def write_worklist(part: ledger_mod.Partition, path: str,
                   run_id_of: Dict[str, str] = None) -> int:
    """Write the worklist TSV, stamping each row with its OWN record's run id.

    A snapshot mixes runs -- a Tier-1 lead last recorded in 2026-07 next to a Tier-2 lead
    from 2026-08 -- so stamping every row with the newest run id would misreport when the
    evidence was last computed.
    """
    rows = ledger_mod.worklist_rows(part)
    for row in rows:
        if run_id_of:
            row["run_id"] = run_id_of.get(row["lead_key"], row["run_id"])
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=ledger_mod.WORKLIST_COLUMNS,
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return len(rows)


def render(part: ledger_mod.Partition, run_id: Optional[str], tag: str = "") -> str:
    if run_id:
        title = ("## Kalamari curation-triage — worklist for run %s" % run_id)
    else:
        title = ("## Kalamari curation-triage — worklist snapshot%s"
                 % ((" at %s" % tag) if tag else ""))
    return ledger_mod.render_markdown(part, title=title)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ledger", default=DEFAULT_LEDGER)
    ap.add_argument("--reviews", default=DEFAULT_REVIEWS)
    ap.add_argument("--run-id", default=None,
                    help="render exactly this run; default is the latest record per "
                         "lead across every signal (the snapshot a release wants)")
    ap.add_argument("--tag", default="",
                    help="label for the snapshot heading, e.g. the git tag")
    ap.add_argument("--out", default=None, help="worklist TSV to write")
    ap.add_argument("--summary-md", default=None, help="markdown to write")
    args = ap.parse_args(argv)

    records = ledger_mod.load_ledger(args.ledger)
    if not records:
        sys.stderr.write("worklist_from_ledger: %s is empty or missing -- nothing to "
                         "render. A release asset needs a committed ledger.\n"
                         % args.ledger)
        return 1
    reviews = ledger_mod.load_reviews(args.reviews)
    part = partition_from_ledger(records, reviews, args.run_id)
    if args.run_id and not (part.new or part.standing or part.accepted or
                            part.dismissed or part.resolved):
        sys.stderr.write("worklist_from_ledger: no ledger record carries run_id %s "
                         "(present: %s)\n"
                         % (args.run_id,
                            ", ".join(sorted({r.get("run_id", "") for r in records}))))
        return 1

    n_rows = 0
    if args.out:
        n_rows = write_worklist(part, args.out,
                                _run_id_of(select_records(records, args.run_id)))
        print("worklist: %d row(s) -> %s" % (n_rows, args.out))
    md = render(part, args.run_id, args.tag)
    if args.summary_md:
        os.makedirs(os.path.dirname(os.path.abspath(args.summary_md)), exist_ok=True)
        with open(args.summary_md, "w", encoding="utf-8") as fh:
            fh.write(md + "\n")
        print("markdown: %s" % args.summary_md)
    if not args.out and not args.summary_md:
        print(md)
    s = ledger_mod.summarize_partition(part)
    print("rendered %s: %d new, %d standing, %d accepted, %d dismissed, %d resolved"
          % (part.run_id or "(snapshot)", s["new"], s["standing"], s["accepted"],
             s["dismissed"], s["resolved"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
