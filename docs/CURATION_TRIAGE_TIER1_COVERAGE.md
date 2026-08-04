# Tier 1 coverage check — does NCBI already have a verdict on every pick?

Tier 1's **pick-health sentinel** (signal J1a) reads NCBI's own precomputed verdict on
each curated Kalamari pick. That works only if the verdict exists for our exact assembly,
so the whole design rests on the number below. This probe measures it: it joins every
eligible pick to NCBI's precomputed reports and states what fraction have a **live row
for our assembly**. The result is evidence about the data, so re-run the probe whenever
the pick set changes or the NCBI reports change shape.

- **Code:** [`bin/curation-triage/`](../bin/curation-triage/) —
  [`coverage_probe.py`](../bin/curation-triage/coverage_probe.py) (CLI),
  [`coverage.py`](../bin/curation-triage/coverage.py) (pure join/flag/summary),
  [`ani_reports.py`](../bin/curation-triage/ani_reports.py) (cached, mockable NCBI layer)
- **Input:** [`src/curation-triage/taxon_policy.tsv`](../src/curation-triage/taxon_policy.tsv)
  ([Tier 0](CURATION_TRIAGE_TIER0.md); the `eligible` rows are the join set, keyed on
  `parent_assembly`)
- **Outputs:** [`src/curation-triage/coverage_report.tsv`](../src/curation-triage/coverage_report.tsv)
  (one row per eligible pick) and
  [`assembly_reports_cache.json`](../src/curation-triage/assembly_reports_cache.json) (distilled cache)
- **Design:** [`CURATION_TRIAGE_DESIGN.md`](CURATION_TRIAGE_DESIGN.md) §4 (Tier 1 J1a)

## Headline result — the signal is there

**249 / 249 (100%) of eligible picks have a live row in NCBI's ANI report**, and the
same 249/249 have a live `assembly_summary` row.

Thin coverage is the obvious risk, and the reason it does not bite is worth knowing.
The ANI report keeps only the *latest* assembly per taxon, and many Kalamari picks are
decades old. The join key decides the outcome: joining on the INSDC replicon accession
a pick records (e.g. `AE002098`, year 2000) finds little, while joining on the **real
GCF/GCA parent assembly** that [Tier 0](CURATION_TRIAGE_TIER0.md) resolves finds
everything, because NCBI keeps an ANI-report row for those assemblies regardless of
age. Never join this report on replicon accessions.

| Coverage (live row for *our* assembly) | count | % |
|---|---|---|
| ANI report (`ANI_report_prokaryotes.txt`) | 249/249 | 100% |
| &nbsp;&nbsp;· exact version | 246 | |
| &nbsp;&nbsp;· superseded to a newer version | 3 | |
| `assembly_summary` (refseq ∪ genbank) | 249/249 | 100% |
| currently live in RefSeq | 228/249 | 91.6% |

Coverage is 100% across **every** breakdown — source lane (NCBI-REF / SME /
FDA-ARGOS / NCTC3000), assembly release year (even the ≤2010 cohort), and GCF vs GCA.
The 21 picks not "live in RefSeq" are still present in `assembly_summary_genbank`
(all `version_status=latest`) with an `excluded_from_refseq` reason — a signal in its
own right (see below), not a coverage gap.

## What the tool does

1. **Reads the eligible picks** from `taxon_policy.tsv` and their source lane from
   `chromosomes.tsv`.
2. **Streams + distils** the three NCBI ASSEMBLY_REPORTS files (~2.7 GB total),
   keeping only the rows that touch our picks plus per-species presence counts. The
   distilled result is the committed `assembly_reports_cache.json`, so re-runs and
   all tests are offline. Two NCBI traps are handled by the files themselves: the ANI
   report carries **both** the GCA and GCF accession per row, and each
   `assembly_summary` row carries its twin in `gbrs_paired_asm`, so **GCF↔GCA
   equivalence is free**; accessions are matched on their version-stripped base and
   the held version is compared so a **superseded** pick is surfaced.
3. **Types each match** — `exact` / `version-differs` / `taxon-only` / `none`.
4. **Previews the J1a flags** honouring NCBI's `approved-mismatch` and *not* firing
   the identity flags on the confounded taxa (where an ANI/type-strain mismatch is
   expected and is resolved by markers, never by ANI).

```bash
# first run: streams the 3 NCBI files, writes the committed cache + report
python3 bin/curation-triage/coverage_probe.py
# offline rerun from the committed cache (no network):
python3 bin/curation-triage/coverage_probe.py --offline
```

## J1a flag preview — real precision, not just coverage

Coverage being 100% is necessary but not sufficient; the probe also previews how many
picks the sentinel would actually *flag*. Deferring to NCBI's own
`taxonomy-check-status` (rather than raw best-match differences) keeps precision high:

**21 / 249 picks carry an actionable flag** on the committed full-database sweep:

| flag | count | meaning |
|---|---|---|
| `tax_check_inconclusive` | 17 | NCBI could not confirm the label |
| `tax_check_failed` | 2 | NCBI's ANI check contradicts the label |
| `refseq_excluded_identity` | 3 | dropped from RefSeq for an identity reason (*unverified source organism*) |
| `ncbi_reclassified` | 1 | NCBI now files the assembly under a different species |

The two hardest catches are real curation leads: *Citrobacter freundii* (best ANI
match *C. cronae*, "unverified source organism", dropped from RefSeq) and
*Lysinibacillus sphaericus* (best match *L. pinottii*, taxonomy-check Failed); plus
*Clavibacter michiganensis* subsp. *sepedonicus*, which NCBI has elevated to the
species *C. sepedonicus*. A RefSeq withdrawal for an identity/contamination reason
stays actionable **even for a confounded taxon** (e.g. *Bacillus cereus*, withdrawn
as "unverified source organism") — it is orthogonal to the ANI-can't-split ambiguity
that markers resolve, so it must not be hidden.

**Deliberately *not* actionable** (captured as context in the report, never on the
worklist):
- **`best_match_species_mismatch` when taxonomy-check is OK** (41 picks). A divergent
  best-ANI-match with an OK verdict is a correctly-named-but-divergent reference
  (often its type strain simply isn't in the ANI panel). NCBI's verdict already
  integrates the best match, so we defer to it. An SME can still filter the report
  for these.
- **`refseq_excluded_qc`** (18) — dropped from RefSeq for annotation/assembly quality
  (*many frameshifted proteins*, *fragmented assembly*), not identity.
- **`assembly_superseded`** (3) — our recorded accession version is stale; the genome
  is essentially the same.
- **`confounded_ani_ambiguity`** (4 of the 21 confounded units) — the expected
  ANI/type-strain ambiguity for the marker-resolved taxa.

## Output columns (`coverage_report.tsv`)

`unit_id`, `scientificName`, `source_lane`, `species_taxid`, `declared_species`,
`parent_assembly`, `ani_match`, `summary_match`, `in_refseq`, `taxonomy_check`,
`best_match_status`, `best_match_species`, `ncbi_current_species`, `version_status`,
`excluded_from_refseq`, `relation_to_type_material`, `seq_rel_date`, `confounded`,
`flags`.

## Design: pure logic vs. the network

Exactly like Tier 0 split `ncbi.py` from `policy.py`, the network is quarantined:
- [`ani_reports.py`](../bin/curation-triage/ani_reports.py) — the **only** network
  code: streaming distillers + injectable line openers + the JSON cache. Column
  layout is asserted at parse time, so an NCBI schema change fails loudly (an "epoch
  bump") instead of silently reading the wrong column.
- [`coverage.py`](../bin/curation-triage/coverage.py) — pure join, match-typing, flag
  logic, and summary; no I/O, unit-testable against tiny fixtures.

## Tests

```bash
cd bin/curation-triage && python3 -m pytest tests/ -q
```

`tests/test_coverage.py` (pure flag/join logic), `tests/test_ani_reports.py`
(distillers, header-schema guard, cache roundtrip, GCF↔GCA + version), and
`tests/test_coverage_cli.py` (the whole path end-to-end via injected local files, no
network) — all offline.

## What this measurement settles

With **100% ANI-report coverage**, **Tier 1's J1a sentinel carries the drift/mislabel
half of the curation job on its own** — no self-run type-strain ANI fallback is needed
*for coverage*. Tier 2's own-ANI compute is therefore reserved for what metadata is
structurally blind to: the sub-species synthetic lineages, and confirmation of the
handful of Tier-1 survivors. See [Tier 1](CURATION_TRIAGE_TIER1.md) and
[Tier 2](CURATION_TRIAGE_TIER2.md).
