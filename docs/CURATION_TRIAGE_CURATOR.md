# Curation triage — for curators: reading the worklist and recording a decision

This is the only document you need to action a worklist.

- Why it works this way: [design](CURATION_TRIAGE_DESIGN.md)
- How to run it: [running it](CURATION_TRIAGE_RUNNING.md)

---

## 1. What you get, and what it will never do

Once a month the tool hands you a short list of **leads**. A lead is one curated Kalamari
genome, or one missing species, worth a look. It answers two questions: has one of our
current picks gone bad, and is there a species we should track but do not?

**It is informational. It is never a pass/fail gate.** It flags and ranks; it does not
grade. Kalamari's genomes were never chosen for being statistically average, so a score
would fail genomes that are correctly curated. It also opens no pull requests — an edit
to [`src/chromosomes.tsv`](../src/chromosomes.tsv) is your decision.

It is built for **precision**, not recall. Your time is the scarce resource, so it would
rather show you 21 leads that are all worth reading than 200 that are not.

## 2. What gets flagged

| Signal | The question it asks | Source |
| --- | --- | --- |
| **J1a** | Has one of *our* picks gone bad? | NCBI's own verdict on our exact assembly |
| **J1b** | Has NCBI named a different reference genome? | NCBI's reference designation |
| **J2** | Is there a species we should add but do not have? | The [wishlist](../src/chromosomes-todo.tsv), and newly named species from LPSN |
| **T2** | What does the sequence itself say? | Our own ANI and marker genes — see §4 |

J1a is the anchor, and it is complete: **all 249 checkable picks have a live NCBI
verdict**. Six flags make a J1a lead:

| Flag | Meaning |
| --- | --- |
| `tax_check_failed` | NCBI's check contradicts the name we give the genome |
| `tax_check_inconclusive` | NCBI could not confirm the name |
| `ncbi_reclassified` | NCBI now files this assembly under a different species |
| `refseq_excluded_identity` | RefSeq dropped it for an identity or contamination reason |
| `assembly_suppressed` / `assembly_replaced` | NCBI suppressed or replaced the assembly |

On the committed run, 21 of 249 picks carry at least one: 17 inconclusive, 2 failed, 3
identity exclusions, 1 reclassified. A pick can carry more than one.

**J1b only fires for picks in the NCBI-REF lane.** For a pick an expert chose, being
different from the type strain is *deliberate*, so the tool never nags about it. Those
picks are still watched by J1a.

### What is deliberately not flagged

| State | Count | Why it is not a lead |
| --- | --- | --- |
| Best ANI match differs, but NCBI's check says OK | 41 | A correctly-named but divergent reference. We follow NCBI's verdict, which already weighs the best match. |
| Dropped from RefSeq for **quality** | 18 | Frameshifted proteins, fragmented assembly. Not about identity. |
| Our accession version is stale | 3 | A newer version of the same genome. |
| Expected ANI ambiguity in a confounded taxon | 4 | Similarity is the wrong ruler here; markers answer these (§4). |

The report still records all four, so you can filter for them.

**One exception is deliberate: a RefSeq withdrawal for an identity reason stays
actionable even for a confounded taxon.** A contaminated *Bacillus cereus* is a real
problem, unrelated to the ambiguity markers resolve, so it must not be hidden.

## 3. How to read a lead

| Column | What to do with it |
| --- | --- |
| `status` | `new`, `standing`, `accepted`, `resolved` or `dismissed` |
| `signal` | J1a / J1b / J2 / T2 |
| `lead_key` | The lead's stable name. **Copy into `reviews.tsv`.** |
| `actionable_flags`, `summary` | Why it is on the list, short and long |
| `link` | Straight to the genome at NCBI |
| `evidence_fingerprint`, `epoch` | **Copy both into `reviews.tsv`** — they decide how long your decision lasts (§6) |

**`new`** is the work. For J1a it means the flag set **gained** a member since last run: a
pick that went OK → Failed, or was newly withdrawn or reclassified. For the others it
means a candidate or verdict that was not there before.

**`standing`** already surfaced and is not yet actioned. A flag set that only shrinks
stays standing, so the same 21 picks do not shout at you every month. **`resolved`** means
the flag set is now empty.

## 4. What a Tier-2 line adds

Rows with signal `T2` come from the sequence itself: our own ANI against the type strain,
a panel comparison for sub-species lineages, and marker genes where similarity cannot
split a taxon. Four rules govern how you read them.

- **A finding is work; a confirmation is not.** A Tier-2 verdict that *agrees* with the
  label is recorded as an annotation, cross-linked to the Tier-1 row it confirms. It never
  removes that row. The committed run has 33 findings and 42 confirmations.
- **For three pairs of taxa the marker call is the verdict and ANI is only context** —
  *Shigella* vs *E. coli*, *B. anthracis* vs the *B. cereus* group, *Y. pestis* vs
  *Y. pseudotuberculosis*. No ANI number can contradict a marker call there, because only
  one of the two is ever a verdict. The run has 4 such overrides.
- **"Undecided" is an answer.** Three marker components have no runnable tool yet, so
  their tasks say so rather than falling back to weaker evidence. Guessing from a plasmid
  or a toxin gene gives a confident wrong answer.
- **"Not measured" is not a finding.** A unit with no cached evidence produces no row.

An example of the override on real data: *Bacillus thuringiensis* (`GCF_000008505.1`) has
*B. anthracis* as its nearest type strain at 97.84% ANI, but the markers rule anthrax out
— no anthracis chromosome signature, no pXO1, no pXO2.

## 5. Recording a decision

Write one row per decision into [`reviews.tsv`](../src/curation-triage/reviews.tsv):

```
lead_key   decision   evidence_fingerprint   epoch   reviewer   date   note   proposed_edit
```

Copy the first, third and fourth straight from the worklist row. `decision` is one of:

| Decision | Effect |
| --- | --- |
| `accept` | Real. Moves to the *proposed edits* section carrying your `proposed_edit` text; a human then edits `chromosomes.tsv`. |
| `reject` | Not a problem. Suppressed. |
| `defer` | Not now. Suppressed the same way. |

A later `reject` outranks an earlier `accept`, so you change your mind by appending a row.

```
Citrobacter_freundii~GCF_000648515.1   reject   9f3c…   ncbi-assembly-reports@2026-07   jdoe   2026-08-15   known contaminated deposit
Clavibacter_michiganensis_sepedonicus~GCF_000069225.1   accept   1a7e…   ncbi-assembly-reports@2026-07   jdoe   2026-08-15   reclassified to C. sepedonicus   Clavibacter_sepedonicus<TAB>GCF_000069225.1<TAB>31964<TAB>…
```

## 6. How long a decision lasts

How you fill the two evidence columns decides this. The first row is the recommended
default.

| `evidence_fingerprint` | `epoch` | Suppressed… |
| --- | --- | --- |
| filled | filled | **until the evidence changes** (recommended) |
| blank | filled | until the next epoch bump, whatever the evidence says |
| blank | blank | forever |

So two things bring a dismissed lead back:

1. **The evidence changed** — NCBI's verdict moved, or a new candidate appeared. That is a
   genuinely different question, so you should see it again.
2. **The epoch bumped.** The epoch is the version of the tools and data behind a signal,
   not the run month. It moves only when NCBI changes a report's shape, or a pinned tool
   or database version changes. An ordinary monthly rerun does **not** bump it.

For Tier-2 leads the fingerprint holds the **verdict**, not the ANI number. A new skani
release moves the last decimal of almost every comparison; if that re-opened every
dismissal, the tool would train you to ignore it.

One more thing about J2: the LPSN feed shows the newest missing species first, capped at
25 a run. Kalamari is missing about 4000 validly-published species in the genera it
tracks, so the cap is what keeps the list short. **A dismissed lead still occupies a cap
slot**, so dismissing the top 25 does not reveal the next 25 — raise `--max-lpsn`.

## 7. Where the worklist comes from

| Where | When |
| --- | --- |
| The monthly workflow's job summary | Every month, on the 1st |
| `worklist.tsv` + `worklist.md` as a run artifact | Every month; artifacts expire after 90 days |
| Attached to the GitHub Release for a `v*.*.*` tag | Each release — the durable copy |

The release copy is rendered from the committed ledger, so it shows the labels you
actually saw in that run, even a year later.

To reproduce a worklist locally, with no network and no tools installed:

```bash
python3 bin/curation-triage/tier1_metadata.py --offline --run-id "$(date -u +%Y-%m)"
python3 bin/curation-triage/tier2_confirm.py  --offline --run-id "$(date -u +%Y-%m)"
```

Both read the committed caches, so they cannot produce a new number.
