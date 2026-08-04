#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collate.py -- merge the sharded Tier 1 lead fragments into ONE monthly worklist.

In CI the pick set is processed by a sharded matrix (mirroring unit-testing.yml):
each shard runs ``tier1_metadata.py --shard i --shards-total N --emit-leads
leads.i.ndjson`` and writes only its raw leads.  collate then does the delta, the
review suppression, the ledger append, and the render ONCE over the merged union --
so ``standing`` / ``resolved`` are computed against the whole pick set, never
per-shard (a shard that never saw a unit must not conclude it "resolved").

Failure handling honours the plan's §8 "loud-but-non-blocking" rule: a shard whose
fragment is MISSING (the job crashed) is counted, ``shards_failed=K/N`` is printed,
and collate exits non-zero so the (non-required) CI job goes red -- but an empty
fragment (a shard that ran and simply found nothing) is a success, not a failure.

``--shards-total`` is REQUIRED with ``--leads-dir``, and never defaults.  It used to
default to 1, which turned the safety net above inside out: point collate at a
directory of 8 fragments and it would read only ``leads.0.ndjson``, see no missing
file, conclude the run was complete, and mark every unit the other 7 shards held as
``resolved`` -- appending false records to an append-only ledger (reproduced: 17 of
them).  The count is the one thing collate cannot infer, so it must be stated.

    python3 bin/curation-triage/collate.py --leads-dir run/ --shards-total 8 \
        --run-id 2026-08 --ledger src/curation-triage/ledger.ndjson
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import leads as leads_mod  # noqa: E402
import ledger as ledger_mod  # noqa: E402
import tier1_metadata  # noqa: E402  (reuse print_summary / write_markdown / tool_versions)

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")
DEFAULT_LEDGER = os.path.join(SRC, "ledger.ndjson")
DEFAULT_REVIEWS = os.path.join(SRC, "reviews.tsv")
DEFAULT_WORKLIST = os.path.join(SRC, "worklist.tsv")


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def fragment_path(leads_dir: str, prefix: str, shard: int) -> str:
    return os.path.join(leads_dir, "%s.%d.ndjson" % (prefix, shard))


def read_fragment(path: str):
    """Rebuild Lead objects from a shard fragment (or None if the file is missing)."""
    if not os.path.exists(path):
        return None
    out = []
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    out.append(leads_mod.Lead.from_dict(json.loads(line)))
    except (json.JSONDecodeError, ValueError, TypeError) as exc:
        # A truncated/corrupt fragment (e.g. a shard OOM-killed mid-write) must be
        # treated as a FAILED shard, not crash the whole collate -- so the surviving
        # shards' worklist still renders. Return None so gather_leads counts it.
        log("CORRUPT shard fragment (failed shard): %s (%s)" % (path, exc))
        return None
    except (PermissionError, OSError) as exc:
        # An UNREADABLE fragment is the same class of accident as a missing one: the
        # file exists, so the missing-file check above passed, but we cannot see its
        # leads. Letting this escape would crash collate; ignoring it would be worse
        # still -- the run would look complete and falsely 'resolve' that shard's units.
        log("UNREADABLE shard fragment (failed shard): %s (%s)" % (path, exc))
        return None
    return out


def gather_leads(paths):
    """Merge shard fragments; return (leads, n_failed). A missing file = failed shard."""
    merged = []
    seen = set()
    n_failed = 0
    for path in paths:
        frag = read_fragment(path)
        if frag is None:
            n_failed += 1
            log("MISSING shard fragment (failed shard): %s" % path)
            continue
        for lead in frag:
            if lead.lead_key in seen:      # disjoint shards shouldn't overlap; be defensive
                continue
            seen.add(lead.lead_key)
            merged.append(lead)
    return merged, n_failed


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Collate sharded Tier 1 lead fragments.")
    p.add_argument("--leads-dir", default=None,
                   help="Directory of shard fragments named <prefix>.<i>.ndjson.")
    p.add_argument("--frag-prefix", default="leads")
    p.add_argument("--shards-total", type=int, default=None,
                   help="How many shards the matrix ran. REQUIRED with --leads-dir; "
                        "with --leads it is optional but cross-checked against the "
                        "number of paths given.")
    p.add_argument("--leads", nargs="*", default=None,
                   help="Explicit fragment paths (alternative to --leads-dir).")
    p.add_argument("--ledger", default=DEFAULT_LEDGER)
    p.add_argument("--reviews", default=DEFAULT_REVIEWS)
    p.add_argument("--out", default=DEFAULT_WORKLIST)
    p.add_argument("--summary-md", default=None)
    p.add_argument("--run-id", default=None)
    p.add_argument("--datasets-cache", default=tier1_metadata.DEFAULT_DATASETS_CACHE)
    p.add_argument("--lpsn-cache", default=tier1_metadata.DEFAULT_LPSN_CACHE)
    p.add_argument("--no-write", action="store_true")
    args = p.parse_args(argv)

    run_id = args.run_id or datetime.date.today().strftime("%Y-%m")

    # Resolve the fragment list.  The shard count is a REQUIRED input, never a guess:
    # see the module docstring for what the old default of 1 did to the ledger.
    if args.leads:
        paths = list(args.leads)
        if args.shards_total is not None and args.shards_total != len(paths):
            raise SystemExit(
                "collate: --shards-total %d but %d fragment path(s) given with --leads. "
                "One of them is wrong, and guessing would risk marking the unlisted "
                "shards' units 'resolved'." % (args.shards_total, len(paths)))
    elif args.leads_dir:
        if args.shards_total is None:
            raise SystemExit(
                "collate: --shards-total is required with --leads-dir. Pass the SAME "
                "number the tier1 matrix ran with -- collate cannot tell a 1-shard run "
                "from 7 crashed shards, and would mark the missing shards' units "
                "'resolved' in the append-only ledger.")
        if args.shards_total < 1:
            raise SystemExit("collate: --shards-total must be >= 1 (got %d)."
                             % args.shards_total)
        paths = [fragment_path(args.leads_dir, args.frag_prefix, i)
                 for i in range(args.shards_total)]
    else:
        raise SystemExit("collate: give --leads-dir (+ --shards-total) or --leads.")

    leads, n_failed = gather_leads(paths)
    n_total = len(paths)
    log("Merged %d leads from %d/%d shards" % (len(leads), n_total - n_failed, n_total))

    # A partial union is NOT safe to persist: the delta would see the crashed shard's
    # picks as absent and falsely mark them 'resolved', permanently corrupting the
    # append-only ledger.  So on ANY failed/missing/corrupt shard, touch nothing and
    # fail loud (the job is non-required -- see curation-triage.yml).  Rerun when whole.
    if n_failed:
        log("INCOMPLETE run: %d/%d shard fragment(s) missing or corrupt. NOT appending "
            "the ledger or writing the worklist (an incomplete union would falsely mark "
            "unseen picks resolved). Rerun after all shards succeed." % (n_failed, n_total))
        print("shards_failed=%d/%d" % (n_failed, n_total))
        return 1

    prior = ledger_mod.load_ledger(args.ledger)
    reviews = ledger_mod.load_reviews(args.reviews)
    # Recomputed from the same committed caches the shards read, so a collate record
    # stamps the identical resolved epochs a single-process run would (see
    # tier1_metadata.tool_versions for why this is recomputed, not threaded through).
    part = ledger_mod.partition(leads, prior, reviews, run_id=run_id,
                                params=leads_mod.DEFAULT_PARAMS,
                                tool_versions=tier1_metadata.tool_versions(
                                    args.datasets_cache, args.lpsn_cache))

    if not args.no_write:
        ledger_mod.append_ledger(args.ledger, part.records)
        ledger_mod.write_worklist(part, args.out)
        tier1_metadata.write_markdown(part, args)
        log("Appended %d ledger records; wrote worklist %s" % (len(part.records), args.out))

    tier1_metadata.print_summary(part)
    print("shards_failed=0/%d" % n_total)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
