# Tier 1 — the monthly triage net that produces the worklist

Tier 1 is a **monthly, delta-aware net**. It hands subject-matter experts (SMEs) a
short, precise worklist of newly deposited public genomes worth reviewing, plus a
decision loop so dismissed leads do not re-surface. It is **informational, never a
pass/fail gate**, and it optimises for **precision** (few, high-value flags) because
SME time is the scarce resource.

It shares [`coverage.py`](../bin/curation-triage/coverage.py) and
[`ani_reports.py`](../bin/curation-triage/ani_reports.py) with the coverage probe (see
[Tier 1 coverage](CURATION_TRIAGE_TIER1_COVERAGE.md)), and keeps every network step
cached + mockable behind an injectable fetcher, so the whole loop and its tests run
offline.

## The three signals

| Signal | Question | Source | Gating |
| --- | --- | --- | --- |
| **J1a** pick-health | Has one of *our* picks gone bad (mislabel / withdrawal / reclass)? | NCBI's own verdict via `coverage.py` (`ANI_report_prokaryotes` + `assembly_summary`) | The frozen 6-flag `ACTIONABLE_FLAGS` set; confounded/approved demotions already applied |
| **J1b** better-rep | Has a newer, NCBI-designated *reference* genome appeared? | `ncbi-datasets-cli` (`datasets summary … --reference`, cached) | **NCBI-REF lanes only**; never nags SME / FDA-ARGOS / NCTC3000 (deliberate divergence); confounded suppressed |
| **J2** coverage feeder | Is there a species we should add but don't? | `chromosomes-todo.tsv` wishlist + LPSN validly-published species in tracked genera | Wishlist = explicit SME intent (highest precision); LPSN gaps capped per run |

J1a is the anchor. Every eligible pick has a live NCBI verdict — **249 / 249**, measured
in [Tier 1 coverage](CURATION_TRIAGE_TIER1_COVERAGE.md) — so Tier 1 metadata carries the
drift/mislabel half of the curation job on its own, and no self-run ANI is needed for
coverage. J1b and J2 stay **off by default on the command line** (`--with-j1b`,
`--with-j2`, `--with-lpsn`); the repository commits a real `datasets_summary_cache.json`
and a real `lpsn_cache.json`, so the monthly CI run enables **all three** and replays them
offline like J1a.

## Delta semantics — what is *newly* actionable

Each run diffs the current leads against the last ledger record for the same lead:

- **J1a** surfaces as `new` when its **actionable-flag set gains a member** vs last run
  (a set that only shrinks is `recurring`; a set that empties is `resolved`). So a pick
  that flips `OK → Failed`, or is newly excluded/withdrawn/reclassified, surfaces — the
  standing 21 do not re-surface every month.
- **J1b / J2** surface on the **first appearance** of a candidate (the `lead_key`
  embeds the candidate accession, so a fresh candidate is a fresh key).
- On a signal's **epoch bump** (a real NCBI schema change or a pinned tool-version
  change — *not* an ordinary monthly refresh) the prior is re-baselined, so the standing
  set is re-opened for review.

## Suppression — "until evidence changes"

A row in `reviews.tsv` suppresses a lead from the *new* worklist while **both** its
`evidence_fingerprint` and its `epoch` match the lead's. The two columns are tolerant:

| `evidence_fingerprint` | `epoch` | Effect |
| --- | --- | --- |
| filled | filled | suppress **until the evidence changes** (recommended default) |
| blank | filled | suppress **until the next epoch bump**, regardless of evidence |
| blank | blank | suppress **forever** |

`accept` routes the lead to the *proposed-edits* section (an accepted lead becomes a
proposed edit to [`src/chromosomes.tsv`](../src/chromosomes.tsv) — **no auto-PR in
v1**); `reject` / `defer` dismiss it. A later `reject` outranks an earlier `accept`.

## Files

| File | What |
| --- | --- |
| [`ledger.ndjson`](../src/curation-triage/ledger.ndjson) | Append-only, one self-contained JSON record per surfaced/resolved lead per run. The delta baseline. |
| [`reviews.tsv`](../src/curation-triage/reviews.tsv) | SME decisions (accept / reject / defer). Human-edited. |
| `worklist.tsv` | The per-run artifact: surviving actionable rows (new + standing + accepted). Regenerated each run, uploaded as a CI artifact — not version-controlled. |
| `worklist.md` | The `$GITHUB_STEP_SUMMARY` markdown block. Regenerated each run. |

### `ledger.ndjson` record (fields)

```
schema, run_id, signal (J1a/J1b/J2), lead_key, epoch, delta (new/recurring/resolved),
unit_id, scientificName, source_lane, candidate_accession, deposit_date,
actionable_flags[], flags[], evidence{...}, evidence_fingerprint, link, summary,
params{species_ani, af_gate}, tool_versions{assembly_reports_epoch, ncbi_datasets_cli, lpsn}
```

`lead_key` — J1a: `unit_id`; J1b: `unit_id#<candidate_base>`; J2: `todo:<taxid>#<acc>`
or `lpsn:<species>`. `epoch` is the **tooling/schema** epoch of the signal (not the run
month); it moves only on a schema/tool change.

### `reviews.tsv` columns

```
lead_key   decision   evidence_fingerprint   epoch   reviewer   date   note   proposed_edit
```

Worked example — copy the `lead_key` / `evidence_fingerprint` / `epoch` straight from
the worklist artifact:

```
Citrobacter_freundii~GCF_000648515.1   reject   9f3c…   ncbi-assembly-reports@2026-07   jdoe   2026-08-15   known contaminated deposit, leave as-is
Clavibacter_michiganensis_sepedonicus~GCF_000069225.1   accept   1a7e…   ncbi-assembly-reports@2026-07   jdoe   2026-08-15   reclassified to C. sepedonicus   Clavibacter_sepedonicus<TAB>GCF_000069225.1<TAB>31964<TAB>…<TAB>NCBI-REF
```

## Running it

```bash
# Offline monthly run (J1a only, from the committed cache — no network):
python3 bin/curation-triage/tier1_metadata.py --offline --run-id 2026-08

# Refresh the ANI cache from NCBI, then run (first run / cache rebuild):
python3 bin/curation-triage/tier1_metadata.py --run-id 2026-08

# Enable the network feeders (needs ncbi-datasets-cli / an LPSN cache):
python3 bin/curation-triage/tier1_metadata.py --run-id 2026-08 \
    --with-j1b --with-j2

# Sharded (CI): each shard emits raw leads; collate merges + diffs + renders ONCE.
python3 bin/curation-triage/tier1_metadata.py --offline --shard 0 --shards-total 8 \
    --emit-leads leads.0.ndjson
python3 bin/curation-triage/collate.py --leads-dir run --shards-total 8 --run-id 2026-08

# Rebuild the committed datasets cache (needs network + the pinned ncbi-datasets-cli).
# Deliberately no --released-after: see below.
python3 bin/curation-triage/tier1_metadata.py --refresh-datasets-cache \
    --with-j1b --with-j2 --datasets-bin <env-dir>/datasets/bin/datasets
```

**`--released-after` narrows the rebuild; it does not top the cache up.** `--refresh-datasets-cache`
REPLACES `datasets_summary_cache.json`, so a window leaves a file holding only the J1b genomes
released inside it. The monthly run asks with no window, and the cache records the window per
taxid, so an `--offline` run whose question is wider than the cached answer **fails** and names each
gap (`cached only since 07/01/2026, need (all)`) instead of quietly reporting an empty worklist.
So: either every later offline run passes the same window, or the cache is rebuilt unwindowed.
Nothing in CI passes a window — the refresh workflow has no such input — and the committed cache is
unwindowed. The window applies to **J1b only**; J2's wishlist taxids are always fetched unwindowed,
because "does this species have a genome at all?" is a standing fact, and windowing it would report
a species whose only genome predates the window as still uncovered.

`--shards-total` is **required** when collate reads a directory of fragments, and it is
cross-checked when the fragments are named individually. It never defaults: it is the
only way collate can tell "this pick produced no lead" from "the shard that owned this
pick never reported". A count of 1 against an 8-shard matrix reads only
`leads.0.ndjson`, finds no missing file, concludes the run was complete, and marks every
pick the other seven shards held `resolved` — false records appended to an append-only
ledger.

The delta/gating/render logic is pure ([`leads.py`](../bin/curation-triage/leads.py),
[`ledger.py`](../bin/curation-triage/ledger.py)); the network is quarantined in
[`ani_reports.py`](../bin/curation-triage/ani_reports.py) (J1a),
[`datasets_summary.py`](../bin/curation-triage/datasets_summary.py) (J1b) and
[`lpsn.py`](../bin/curation-triage/lpsn.py) (J2), each with an injectable fetcher + a
committed JSON cache + a `covers()`/`--offline` gate. CI is the non-gating
[`curation-triage.yml`](../.github/workflows/curation-triage.yml) (kept **out of required
status checks**; a broken shard goes red but never blocks a merge).

## What happens to a Tier-1 lead next

[Tier 2](CURATION_TRIAGE_TIER2.md) confirms every lead on this worklist at the sequence
level: skani against the type strain of its declared species, plus sub-species panels and
marker resolvers for the 21 confounded units, where the marker call overrides ANI. A
confirming Tier-2 verdict **annotates** the Tier-1 row rather than removing it,
cross-linked by `confirms_lead`, so the curator still decides. Tier 2 marks nothing
resolved — it says nothing about whether a Tier-1 lead is still actionable.

## CI and the release worklist

`ncbi-datasets-cli` is pinned at 18.34.0 by a committed lockfile, and
[`datasets_summary_cache.json`](../src/curation-triage/datasets_summary_cache.json) holds
170 real NCBI queries (3208 records), so the monthly job replays J1b and J2 offline. The
monthly [`curation-triage.yml`](../.github/workflows/curation-triage.yml) runs
`generate-matrix → tier1-shard ×8 → collate → tier2-replay` and makes **one** ledger
commit carrying both tiers, via fetch-rebase-push. A weekly job defeats GitHub's 60-day
auto-disable of scheduled workflows. See
[reproducibility](CURATION_TRIAGE_REPRODUCIBILITY.md).

A `v*.*.*` tag attaches the worklist to the GitHub Release. The worklist files themselves
are gitignored — they are a view of the ledger, and committing a view next to its source
lets the two drift — so the tag job **renders** the worklist from the committed ledger with
[`worklist_from_ledger.py`](../bin/curation-triage/worklist_from_ledger.py). Because the
ledger stores the delta each run computed, the release shows the labels the curator saw;
re-running the pipeline instead would diff against a ledger that already holds those
records and report every lead as standing. See
[reproducibility](CURATION_TRIAGE_REPRODUCIBILITY.md) §4.

The release step **tolerates** a race with upstream's `build-sketch` job rather than
serialising against it. A `concurrency` block serialises a workflow against **itself**,
and upstream's workflow declares none, so no `concurrency` group here can de-conflict the
two; doing it properly would mean editing an upstream file. The step waits for the release
object to appear and updates it instead of creating a second one.

## LPSN valid dates — live, on a derived date

The J2 LPSN feeder filters to correct names, so a synonym or a misspelling cannot fire as
a coverage gap. It takes an optional `--lpsn-since` cutoff and sorts gaps **newest first**
before `--max-lpsn` truncates — so the feed means "the most recently named species in a
genus we track that we do not have", not an alphabetical slice.
[`lpsn_cache.json`](../src/curation-triage/lpsn_cache.json) holds 28794 names from LPSN
release **2026-08-04**, so the feeder replays offline like the others.

The limits are worth knowing. LPSN publishes **no date column**, so the valid-publication
date is derived from the authority string (cut at `emend.`, then the last year in what
remains — which gives *E. coli* its 1980 Approved-Lists year, not the year of a later
emendation). Kalamari lacks about 4000 validly-published species in its tracked genera, so
the cap, not the cutoff, keeps the worklist short — and a dismissed lead still occupies a
cap slot. Refreshing the cache needs an LPSN account: `--lpsn-csv` + `--lpsn-release`.
Detail in [reproducibility](CURATION_TRIAGE_REPRODUCIBILITY.md) §5.

## Status and open items

Built and run offline on the committed tree: the three signals, the delta and suppression
logic, the append-only ledger, the 8-shard run plus `collate.py`, and the worklist
renderer. [`ledger.ndjson`](../src/curation-triage/ledger.ndjson) holds 96 records from
real runs — the 21 J1a leads and the 75 Tier-2 verdicts that confirm them.

Open:

- **No auto-PR.** An accepted lead becomes a *proposed edit* on the worklist; a human
  applies it to [`src/chromosomes.tsv`](../src/chromosomes.tsv).
- **The workflow has never run on a GitHub runner** — no runner is available in the
  development environment. Its step scripts are extracted from the YAML and run locally
  instead, and `actionlint` + `shellcheck` report zero findings; see
  [reproducibility](CURATION_TRIAGE_REPRODUCIBILITY.md) §6.
- **The unattended LPSN refresh.** Rebuilding `lpsn_cache.json` needs an account-gated CSV
  download. The REST API path is not implemented.
