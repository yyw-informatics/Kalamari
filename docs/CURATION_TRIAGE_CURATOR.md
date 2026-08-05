# Curation triage — for curators: reading the worklist and recording a decision

This is the only document you need to action a worklist. It explains what the tool
flags, how to read one line of it, and how to write your decision down so the lead
does not come back next month.

- Why it works this way: [design](CURATION_TRIAGE_DESIGN.md)
- How to run it: [running it](CURATION_TRIAGE_RUNNING.md)

---

## 1. What you get, and what it will never do

Once a month the tool hands you a short list of **leads**. A lead is one curated
Kalamari genome, or one missing species, that is worth a look. It answers two
questions:

- Has one of our current picks gone bad — mislabelled, withdrawn, or reclassified?
- Is there a species we should track but do not?

**It is informational. It is never a pass/fail gate.** It flags and ranks; it does not
grade. Kalamari's genomes were never chosen for being statistically average, so a score
would fail genomes that are correctly curated. The tool also opens no pull requests: an
edit to [`src/chromosomes.tsv`](../src/chromosomes.tsv) is your decision, not its.

It is built for **precision**, not recall. Your time is the scarce resource, so the tool
would rather show you 21 leads that are all worth reading than 200 that are not.

## 2. What gets flagged

Three signals feed the worklist. The `signal` column on every row says which one found
the lead.

| Signal | The question it asks | Where the answer comes from |
| --- | --- | --- |
| **J1a** | Has one of *our* picks gone bad? | NCBI's own verdict on our exact assembly |
| **J1b** | Has NCBI named a different reference genome for a species we track? | NCBI's reference designation |
| **J2** | Is there a species we should add but do not have? | The [`chromosomes-todo.tsv`](../src/chromosomes-todo.tsv) wishlist, and newly named species from LPSN in genera we already track |
| **T2** | What does the sequence itself say? | Our own ANI measurement and marker genes — see §4 |

J1a is the anchor, and it is complete: **all 249 checkable picks have a live NCBI
verdict**. Six flags make a J1a lead:

| Flag | What it means |
| --- | --- |
| `tax_check_failed` | NCBI's own check contradicts the name we give the genome |
| `tax_check_inconclusive` | NCBI could not confirm the name |
| `ncbi_reclassified` | NCBI now files this assembly under a different species |
| `refseq_excluded_identity` | RefSeq dropped it for an identity or contamination reason |
| `assembly_suppressed` | NCBI suppressed the assembly |
| `assembly_replaced` | NCBI replaced the assembly with a different one |

On the committed run, 21 of 249 picks carry at least one: 17 inconclusive, 2 failed,
3 identity exclusions and 1 reclassified. A pick can carry more than one flag.

J1b only fires for picks in the **NCBI-REF** lane. For a pick an expert chose, being
different from the type strain is *deliberate*, so the tool never nags "a type strain
exists" about an SME, FDA-ARGOS or NCTC3000 pick. Those picks are still watched by J1a.

### What is deliberately *not* flagged

These four states are real, and the report records them, but they are kept off the
worklist on purpose. Each one would cost more attention than it returns.

| State | Count | Why it is not a lead |
| --- | --- | --- |
| Best ANI match is a different species, but NCBI's check says OK | 41 | A correctly-named but divergent reference. We follow NCBI's verdict, which already weighs the best match. |
| Dropped from RefSeq for **quality** (`refseq_excluded_qc`) | 18 | Frameshifted proteins or a fragmented assembly. Not a question about identity. |
| Our recorded accession version is stale (`assembly_superseded`) | 3 | A newer version of the same genome. Not a different organism. |
| The expected ANI ambiguity in a confounded taxon | 4 | Similarity is the wrong ruler here. Markers answer these instead (§4). |

The first row is the largest group, and it is worth understanding: a divergent best match
with an OK verdict usually means the species' type strain is simply not in NCBI's
comparison panel. The report still records all four states, so you can filter for them if
you want to look.

One exception is deliberate and worth knowing. **A RefSeq withdrawal for an identity
reason stays actionable even for a confounded taxon.** A contaminated *Bacillus cereus*
is a real problem. It has nothing to do with the ambiguity that markers resolve, so it
must not be hidden.

## 3. How to read a lead

The worklist has one row per lead, in `worklist.tsv` and in the markdown summary. The
columns you use are:

| Column | What to do with it |
| --- | --- |
| `status` | `new`, `standing`, `accepted`, `resolved` or `dismissed` — see below |
| `signal` | J1a / J1b / J2 / T2 |
| `lead_key` | The lead's stable name. **Copy this into `reviews.tsv`.** |
| `scientificName`, `source_lane` | Which pick, and which lane it came from |
| `actionable_flags` | Why it is on the list, in short form |
| `summary` | The same thing in a sentence, with the numbers |
| `link` | Straight to the genome at NCBI |
| `evidence_fingerprint` | A short hash of the evidence behind the lead. **Copy this too.** |
| `epoch` | Which version of the tools and data produced the lead. **Copy this too.** |

`status` tells you whether this is new work:

- **`new`** — worth reading now. For J1a it means the flag set **gained** a member since
  last month: a pick that went from OK to Failed, or was newly withdrawn or
  reclassified. For J1b, J2 and T2 it means a candidate or a verdict that was not there
  before.
- **`standing`** — already surfaced in an earlier run and not yet actioned. A flag set
  that only shrinks stays standing, so the same 21 picks do not shout at you every month.
- **`resolved`** — the flag set is now empty. NCBI changed its mind, or the pick was
  replaced.

## 4. What a Tier-2 line adds

Rows with signal `T2` come from the sequence itself. Three kinds of evidence feed them:
our own ANI measurement against the type strain of the declared species, a small panel
comparison for sub-species lineages, and marker genes for the taxa that similarity cannot
split.

Four rules govern how you read them:

- **A finding is work; a confirmation is not.** A Tier-2 verdict that *agrees* with the
  current label is recorded as an annotation under the finding list, cross-linked to the
  Tier-1 row it confirms. It never removes that row — you still decide. On the committed
  run there are 33 findings and 42 confirmations.
- **For three pairs of taxa the marker call is the verdict, and the ANI number is only
  context.** The pairs are *Shigella* vs *E. coli*, *B. anthracis* vs the *B. cereus*
  group, and *Y. pestis* vs *Y. pseudotuberculosis*. No ANI number can contradict a
  marker call for them, because only one of the two is ever a verdict. The committed run
  has 4 such overrides.
- **"Undecided" is an answer.** Three marker components have no tool that can run them
  yet, so their tasks say `t2_marker_indeterminate` rather than falling back to weaker
  evidence. Guessing from a plasmid or a toxin gene would give a confident wrong answer.
- **"Not measured" is not a finding.** A unit with no cached evidence produces no row at
  all. It is reported as a count, so an unmeasured unit never looks like a clean result.

A worked example from the committed run shows why the override matters.
*Bacillus thuringiensis* (`GCF_000008505.1`) has *B. anthracis* as its nearest type
strain by ANI, at 97.84%. The markers rule anthrax out: no anthracis chromosome
signature, no pXO1, no pXO2. The label is still reported as refuted, because BTyper3
places the genome in a group that does not contain the declared species — but it is not
refuted *as anthrax*.

## 5. Recording a decision

Write one row per decision into [`src/curation-triage/reviews.tsv`](../src/curation-triage/reviews.tsv).
It is a plain tab-separated file, edited by hand:

```
lead_key   decision   evidence_fingerprint   epoch   reviewer   date   note   proposed_edit
```

Copy `lead_key`, `evidence_fingerprint` and `epoch` straight from the worklist row.
`decision` is one of:

| Decision | Effect |
| --- | --- |
| `accept` | The lead is real. It moves to the *proposed edits* section of the worklist, carrying your `proposed_edit` text. A human then edits `src/chromosomes.tsv`. |
| `reject` | Not a problem. The lead is suppressed. |
| `defer` | Not now. The lead is suppressed the same way. |

A later `reject` outranks an earlier `accept`, so you can change your mind by appending
a row rather than editing history.

Two example rows:

```
Citrobacter_freundii~GCF_000648515.1   reject   9f3c…   ncbi-assembly-reports@2026-07   jdoe   2026-08-15   known contaminated deposit, leave as-is
Clavibacter_michiganensis_sepedonicus~GCF_000069225.1   accept   1a7e…   ncbi-assembly-reports@2026-07   jdoe   2026-08-15   reclassified to C. sepedonicus   Clavibacter_sepedonicus<TAB>GCF_000069225.1<TAB>31964<TAB>…<TAB>NCBI-REF
```

## 6. How long a decision lasts

How you fill the two evidence columns decides how long the lead stays quiet. The
recommended default is the first row.

| `evidence_fingerprint` | `epoch` | The lead stays suppressed… |
| --- | --- | --- |
| filled | filled | **until the evidence changes** (recommended) |
| blank | filled | until the next epoch bump, whatever the evidence says |
| blank | blank | forever |

So two things bring a dismissed lead back:

1. **The evidence changed.** The fingerprint no longer matches — NCBI's verdict moved,
   or a new candidate appeared. That is a genuinely different question, and you should
   see it again.
2. **The epoch bumped.** The epoch is the version of the tools and data behind a signal,
   not the run month. It moves only when NCBI changes the shape of a report, or when a
   pinned tool or database version changes. A bump re-opens the standing set, because
   every verdict was recomputed by something different. An ordinary monthly rerun does
   **not** bump it.

For Tier-2 leads the fingerprint deliberately holds the **verdict**, not the ANI number.
A new skani release moves the last decimal of almost every comparison; if that re-opened
every dismissal, the tool would train you to ignore it. A real tooling change is handled
by the epoch instead.

One more thing to know about J2. The LPSN feed shows the **newest** missing species
first, capped at 25 a run (`--max-lpsn`). Kalamari is missing about 4000
validly-published species in the genera it tracks, so the cap, not any date cutoff, is
what keeps the list short. **A dismissed lead still occupies a cap slot**, so dismissing
the top 25 does not reveal the next 25 — raise the cap to see further.

## 7. Where the worklist comes from

| Where | When |
| --- | --- |
| The job summary of the monthly workflow run | Every month, on the 1st |
| `worklist.tsv` + `worklist.md`, uploaded as a run artifact | Every month; artifacts expire after 90 days |
| Attached to the GitHub Release for a `v*.*.*` tag | On each release — this is the durable copy |

The release copy is rendered from the committed ledger, so it shows the labels you
actually saw in that run, even a year later.

To reproduce a worklist locally, with no network and no tools installed:

```bash
python3 bin/curation-triage/tier1_metadata.py --offline --run-id "$(date -u +%Y-%m)"
python3 bin/curation-triage/tier2_confirm.py --offline --run-id "$(date -u +%Y-%m)"
```

Both read the committed caches, so they cannot produce a new number — they show you the
same evidence the monthly job saw. Everything else about running the tool is in
[running it](CURATION_TRIAGE_RUNNING.md).
