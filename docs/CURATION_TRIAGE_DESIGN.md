# Curation triage — how it works, and why to trust it

This is the reference for anyone reviewing the add-on or taking it over. It covers the
architecture, the reasoning behind each rule, the one measurement the whole approach
rests on, and an honest account of what is not proven.

Every rule below is written with its reason. A rule without its reason gets deleted by
the next person who reads it.

- To action a worklist: [for curators](CURATION_TRIAGE_CURATOR.md)
- To run, refresh or re-pin anything: [running it](CURATION_TRIAGE_RUNNING.md)

---

## 1. What the tool is, and what it refuses to be

It gives Kalamari's subject-matter-expert (SME) curators a monthly **triage worklist**
of publicly deposited genomes worth examining. It turns curatorial choices that were
made once into choices that are watched continuously. There are two jobs:

- **J1 — representative drift**, for a species we already track: has our pick gone bad
  (mislabelled, withdrawn, reclassified), or has a better *trusted* genome appeared?
- **J2 — coverage gap**: is there a species we should add but do not have?

**It is informational, never pass/fail.** Being the statistical centre of a species was
never how these genomes were chosen, so the tool flags and ranks for human review — it
does not grade. SME time is the scarce resource, so the whole design optimises for
**precision**: few flags, each worth reading.

## 2. The three ideas that shape everything

**1. Kalamari picks on trusted provenance, not on centrality.** The `source` column of
[`src/chromosomes.tsv`](../src/chromosomes.tsv) holds 179 NCBI-REF rows, 130 SME rows
across three SME lanes, 10 FDA-ARGOS and 4 NCTC3000. These genomes were chosen for
completeness and trust. So the triage signal has to track provenance and correctness,
not statistical centrality.

**2. Do not recompute what NCBI already computes.** NCBI publishes a roughly nightly
verdict for every prokaryotic assembly, comparing it against the type strain of its
species using ANI (average nucleotide identity, a genome-wide percent-similarity
measure). Reading *its verdict on our own picks* is nearly free, and it catches the case
a similarity computation is structurally blind to: our pick was reclassified or
withdrawn.

**3. For the confounded taxa, genome-wide similarity is the wrong ruler.** For the
bacteria the Kalamari paper highlights — *Shigella*, anthrax, plague, and the
*Listeria* / *Salmonella* / botulinum splits — the meaningful difference lives in
specific genes: a virulence plasmid, a toxin gene. So the resolver is **targeted
marker-gene detection**, which is interpretable and already validated in public-health
laboratories. Not ANI, and not a learned embedding.

Together they give the shape of the whole thing: **a cheap metadata net finds the few
candidates → expensive confirmation runs only on those → the confounded taxa are
resolved by markers, not by similarity.**

## 3. The three tiers

### Tier 0 — the policy table

[`build_policy.py`](../bin/curation-triage/build_policy.py) reads `src/chromosomes.tsv`
plus NCBI taxonomy and writes the committed
[`taxon_policy.tsv`](../src/curation-triage/taxon_policy.tsv): one row per
**genome-unit**, 270 units built from 323 replicon rows. It decides what each Kalamari
entry really is and how it should be handled. Every later tier keys off it.

**Group replicons into units by their parent assembly — never by the Kalamari `taxid`
column.** The taxid column mixes ranks and is not trustworthy for grouping. Taxid `630`
covers **two different *Yersinia enterocolitica* strains**, which grouping by taxid would
merge into one chimera. In the other direction, the two chromosomes of one *Aliivibrio*
genome would look like two organisms. The real GCF/GCA assembly, resolved from NCBI, is
the ground truth: it merges the replicons of one genome and keeps distinct genomes apart
even when they share a species. The eleven *Salmonella enterica* lineage entries all
resolve to species taxid `28901`, and they are eleven different assemblies.

Where no real assembly resolves, the fallback grouping key is scoped to the Kalamari
identity (`scientificName` plus `taxid`), so the fallback can never merge two distinct
genomes that happen to share a species-level organism taxid.

**Resolve a real NCBI species-rank taxid.** The 13 synthetic Kalamari taxids
(`9000000`–`9000017`, such as the four *Listeria monocytogenes* lineages) do not exist
at NCBI; they climb through the bundled custom-taxonomy overlay to their real parent
species. Strain and subspecies taxids climb to species, and merged taxids are followed.
A genus-level entry (`Ruminococcus_sp.`, taxid `1263`) falls back to the deposited
genome's own organism taxid, which is a real species. **A preflight check hard-fails**
if any eligible unit lacks a real species-or-below taxon, so a synthetic, strain-level
or genus-level key can never pass silently.

**Decide eligibility with a programmatic guard, not a hand-written list.** A hand list
misses what nobody thought of; a guard catches it.

| Class | Units | Meaning |
| --- | --- | --- |
| `eligible` | 249 | Bacteria and Archaea with a real species taxid. These are checked. |
| `exclude:organelle` | 14 | An organelle-only submission, or an assembly over ~50 Mbp. The guard is what catches the *Arripis trutta* fish mitogenome. |
| `exclude:virus` | 4 | Viruses and viroids, including multi-segment genomes. |
| `eukaryote` | 3 | Eukaryote pathogens (*Cryptosporidium*, *Cryptococcus*). Recorded, outside the pipeline. |
| `review:genus-level` | 0 | Bacteria that genuinely resolve no further than genus. No unit lands here today; the class exists so such a unit is surfaced rather than passing silently. |

The organelle guard only fires when **every** replicon of a unit is non-nuclear. A
eukaryote nuclear assembly that merely bundles its mitochondrion (*Cryptococcus*: 14
chromosomes plus 1 mitochondrion) is therefore `eukaryote`, not `exclude:organelle`.

**Tag the confounded taxa — 21 units** — from the committed
[`confounded_manifest.tsv`](../src/curation-triage/confounded_manifest.tsv), which sets
a `regime` (A or B) and a `resolver` (marker or ANI). A confounded unit **stays
eligible**. It is held out of the better-representative feed, not dropped: these units
still need pick-health monitoring, they just must not be nagged about diverging from a
type strain, because for them that divergence is expected (§5).

### Tier 1 — the monthly metadata net

Tier 1 needs no genome downloads. It produces the worklist. Its three signals, its delta
rules and its suppression rules are described from the curator's side in
[for curators](CURATION_TRIAGE_CURATOR.md); what follows is why each one is shaped the
way it is.

- **J1a, the pick-health sentinel, is the primary signal.** It bulk-diffs our ~270
  parent assemblies against NCBI's `ANI_report_prokaryotes.txt` and
  `assembly_summary_{refseq,genbank}.txt`, and flags a frozen set of six conditions. It
  defers to NCBI's own `taxonomy-check-status` rather than to raw best-match
  differences, because NCBI's verdict already weighs the best match and a raw difference
  does not. It respects NCBI's `approved-mismatch` marker and never flags those.
- **J1b, the better-representative feed, is source-gated.** For NCBI-REF units it fires
  when NCBI's reference designation changes. For SME, FDA-ARGOS and NCTC3000 units it
  never fires on type-strain divergence, because that divergence is a deliberate
  curatorial choice; those picks are covered by J1a instead. Confounded units are
  suppressed here.
- **J2, the coverage feeder, is fed by intent rather than by scanning.** The
  [`chromosomes-todo.tsv`](../src/chromosomes-todo.tsv) wishlist is explicit curator
  intent, which is the highest-precision source there is. LPSN supplies newly
  validly-published species, limited to genera we already track.
- **The ledger is the delta baseline.** Each run appends to the append-only
  [`ledger.ndjson`](../src/curation-triage/ledger.ndjson) and reads
  [`reviews.tsv`](../src/curation-triage/reviews.tsv) to suppress dismissed leads.

### Tier 2 — targeted confirmation

Tier 2 answers what metadata cannot, at the sequence level, and only for the units that
need it. It is **confirmation and sub-species residue, never coverage**: §4 shows Tier 1
already sees every pick, so nothing here exists to find leads Tier 1 missed.

Two trigger groups, 41 units in total: the 21 Tier-1 survivors on the worklist, plus the
21 confounded units whether flagged or not, overlapping in one unit
(`Bacillus_cereus~GCF_000283675.1`). Three checks run:

| Check | Runs on | The question |
| --- | --- | --- |
| `ani` | Tier-1 survivors and regime B | How far is our pick from the type strain of its declared species? |
| `panel` | Regime B (17 units) | Is this lineage reference still distinct from its siblings, and where does it sit? |
| `marker` | All confounded units | What do the defining genes say? |

**The aligned fraction is a veto, not a tiebreaker.** An ANI computed over a sliver of a
genome is not weak evidence, it is no evidence. So the gate checks aligned fraction
first, and takes the **smaller** of the two sides: a small genome fully contained in a
larger one aligns well on its own side while covering little of the other, which is
exactly the case the gate exists to catch. Below the gate the verdict is "cannot
decide", never "different species". Within `epsilon` of the species boundary the answer
is "too close to call" — a precision-first tool has to be allowed to say that.

**The sub-species panels compare Kalamari against itself.** The 13 units with synthetic
taxids are metadata-blind: NCBI knows nothing about a lineage only Kalamari defines, so
Tier 1's whole signal is unavailable for them. Their panel members are Kalamari's own
sibling units — real accessions already curated in `chromosomes.tsv`, generated into
[`panels/*.acclist`](../src/curation-triage/panels/) (17 files) by
[`build_panels.py`](../bin/curation-triage/build_panels.py), with a test asserting the
committed files match so a policy edit cannot silently change them. That answers two
questions offline, without inventing reference accessions nobody has verified:

- **Is the panel still a panel?** A sibling at or above the dereplication line (99.5%)
  means two of Kalamari's own lineage references are effectively the same genome. That
  is a curation defect nothing else in the tool can see.
- **Where does a genome sit?** The nearest member, and the margin to the runner-up.

Because the unit's own lineage is not in its panel, "nearest" is a neighbourhood, not a
lineage assignment, and the code says so rather than over-claiming.

## 4. The measurement the design rests on

J1a works only if NCBI publishes a live verdict row for *our exact assembly*. That is an
assumption about somebody else's data, so it is measured rather than assumed.
[`coverage_probe.py`](../bin/curation-triage/coverage_probe.py) joins every eligible pick
to NCBI's precomputed reports and writes
[`coverage_report.tsv`](../src/curation-triage/coverage_report.tsv).

**249 of 249 eligible picks (100%) have a live row in NCBI's ANI report**, and the same
249 have a live `assembly_summary` row.

| Coverage — a live row for *our* assembly | Count | % |
| --- | --- | --- |
| ANI report (`ANI_report_prokaryotes.txt`) | 249/249 | 100% |
| · at the exact version we record | 246 | |
| · superseded to a newer version | 3 | |
| `assembly_summary` (refseq ∪ genbank) | 249/249 | 100% |
| currently live in RefSeq | 228/249 | 91.6% |

Coverage is 100% across every breakdown: source lane (152 NCBI-REF, 83 SME, 10
FDA-ARGOS, 4 NCTC3000), release year including the pre-2010 cohort, and GCF versus GCA.
The 21 picks not live in RefSeq are all present in `assembly_summary_genbank` with an
`excluded_from_refseq` reason — that is a signal in its own right, not a coverage gap.

**Thin coverage is the obvious risk, and the join key is why it does not bite. Never
join this report on replicon accessions; join on the parent assembly.** NCBI's ANI
report keeps only the latest assembly per taxon, and many Kalamari picks are decades
old. Joining on the INSDC replicon accession a pick records — `AE002098`, from the year
2000 — finds almost nothing. Joining on the real GCF/GCA parent assembly that Tier 0
resolves finds everything, because NCBI keeps a row for those assemblies regardless of
age.

Two NCBI conveniences make the join cheap, and both are properties of the files rather
than of our code: the ANI report carries **both** the GCA and the GCF accession on every
row, and each `assembly_summary` row carries its twin in `gbrs_paired_asm`, so GCF↔GCA
equivalence is free. Accessions are matched on their version-stripped base and the held
version is then compared, so a superseded pick is surfaced rather than missed.

**What the measurement settles.** Tier 1 metadata carries the drift and mislabel half of
the curation job on its own. No self-run type-strain ANI is needed *for coverage*. That
is what frees Tier 2's own ANI for the two things metadata is structurally blind to: the
sub-species synthetic lineages, and confirmation on the handful of Tier-1 survivors.

Re-run the probe whenever the pick set changes or NCBI changes the shape of its reports.
Its network layer asserts the column layout at parse time, so a schema change fails
loudly — an epoch bump — instead of silently reading the wrong column.

## 5. The confounded taxa: markers, not similarity

These 21 units split into two regimes, and only regime A is a case of "ANI fails".

| Regime | Taxa | Why ANI is the wrong ruler | Resolver |
| --- | --- | --- | --- |
| **A — ANI cannot split them** | *Shigella* ↔ *E. coli*; *B. anthracis* ↔ the *B. cereus* group; *Y. pestis* ↔ *Y. pseudotuberculosis* | Near-identical backbone; the signal is a small mobile or accessory gene set | **Marker genes**, which override ANI |
| **B — ANI resolves or over-splits** | *L. monocytogenes* lineages; *Salmonella* subspecies and serovars; *C. botulinum* groups | Core-genome divergence is already captured; the residual problem is *labelling* | An ANI panel places the genome; a typer supplies the label |

**The regime-A override is structural, not a preference.** A regime-A unit gets no
independent ANI verdict at all: the ANI number is folded into the marker lead as context
and the marker call is the verdict. There is no code path in which an ANI number can
contradict a regime-A marker call, because only one of them is ever a verdict. If the
markers cannot decide, the answer is "undecided" — it does **not** fall back to ANI,
because for these taxa the ANI answer was never admissible in the first place.

For regime B the panel decides membership and the typer supplies the label. A wrong
label on a correctly-placed genome is still a defect, so a refuting label is reported.

### How marker information is obtained

Wrap validated tools; do not hand-curate gene lists. Three layers feed the calls:
curated gene and virulence databases (VFDB, the NCBI reference gene catalog behind
AMRFinderPlus, CGE VirulenceFinder and PlasmidFinder); validated in-silico typers that
bundle their own database and interpretation (SISTR, SeqSero2, ShigEiFinder, ECTyper,
BTyper3, LisSero); and allele or scheme databases for the labelling splits (PubMLST and
Institut Pasteur BIGSdb).

The committed [`marker_manifest.tsv`](../src/curation-triage/marker_manifest.tsv) maps
each confounded taxon to its resolver tool, database, database version and defining
marker genes. It holds **15 components over 8 identification tasks**. A *task* is the
question ("is this anthracis?"); a *component* is one tool run that contributes evidence
to it. Twelve components are pinned and runnable; three are gaps.

| Task | Resolver (installed build) | Defining markers |
| --- | --- | --- |
| `shigella_eiec_vs_escherichia_coli` | ShigEiFinder `1.3.4-pyhdfd78af_0` | ipaH, cluster-specific genes, O/H (wzx/wzy, fliC) |
| `bacillus_anthracis_vs_b_cereus_group` | BTyper3 `3.4.0-pyhdfd78af_0` **plus our rule** | pXO1: pagA, lef, cya · pXO2: capA, capB, capC |
| `yersinia_pestis_vs_y_pseudotuberculosis` | **gap** — no runnable cgMLST caller; `mlst 2.23.0` is context only | chromosome: ypo2088, yihN (*pestis*), opgG (*pseudotuberculosis*) · pMT1: caf1 · pPCP1: pst |
| `listeria_monocytogenes_serogroup_lineage_cc` | LisSero `0.4.10-pyhdfd78af_0` + `mlst 2.23.0`; cgMLST/LIN is a **gap** | Doumith pattern: prs, lmo0737, ORF2819, ORF2110, lmo1118 |
| `salmonella_subspecies_and_serovar` | SISTR `1.1.3-pyhdc42f0e_2` + SeqSero2 `1.3.2-pyhdfd78af_0` | O (wzx/wzy), H (fliC/fljB), cgMLST330 |
| `clostridium_botulinum_group_and_bont_type` | AMRFinderPlus `4.2.7-hf69ffd2_0`, database `2026-05-15.1` | bont/A–G, with subtype where the evidence supports it |
| `escherichia_coli_serotype_pathotype` | ECTyper `2.0.0-pyhdfd78af_4` — context only, never a verdict | O/H alleles, stx subtype |
| `generic_marker_fallback` | abricate `1.4.0` over VFDB and PlasmidFinder; AMRFinderPlus `4.2.7` | any curated virulence or plasmid gene |

### Why these tools, and not the obvious ones

Four rows reject the tool an experienced reader would reach for first. Three of the four
reject `abricate` as a **classifier**. Each rejection is a rule, not a taste.

| Task | The obvious choice | Why it is wrong | What is used instead |
| --- | --- | --- | --- |
| *Y. pestis* vs *pseudotuberculosis* | abricate over VFDB / PlasmidFinder | **Plasmid presence alone misclassifies plasmid-cured *pestis* and plasmid-bearing *pseudotuberculosis*.** A plasmid hit must never make the species call. | BIGSdb-Pasteur cgMLST scheme 1 — a gap, so the task returns "undecided"; abricate stays as a cheap pre-screen |
| *C. botulinum* toxin call | abricate with a bont gene set | abricate's databases carry no version of their own, so a changed call cannot be attributed to anything. | AMRFinderPlus `--plus`, the one tool with a first-class dated database version, so a changed call reads as a database update rather than as biology |
| *L. monocytogenes* lineage | merge serogroup, MLST and cgMLST into one "lineage" field | Serogroup, classical MLST ST/CC and cgMLST/LIN answer **different questions**; merging them destroys the distinction. | The same tools, with the outputs kept as separate fields |
| *B. anthracis* vs *B. cereus* | trust BTyper3's species call | BTyper3 gives **group taxonomy and a marker profile**, not an anthracis verdict. Turning that into "is this anthracis" is our judgement. | BTyper3 plus our own interpretation rule, versioned as ours |

Two more rules the marker logic follows, each because breaking it produces a confident
wrong answer. Plasmid-cured *B. anthracis* is still *B. anthracis*, and a *B. cereus*
carrying anthrax plasmids is never labelled plain *B. cereus*. And "anthracis" is matched
as a species epithet, never as a substring: *Bacillus paranthracis* contains the string,
is a different organism, and is what NCBI's best ANI match reports for one of Kalamari's
own picks.

### The three gap components

Three components have **no runnable, pinnable caller**. They are recorded as
`status = gap` in the manifest rather than left as pins nothing can execute. A gap makes
its task return **"undecided"** — deliberately, because the alternative is falling back
to weaker evidence.

| Component | Why it cannot run | Consequence |
| --- | --- | --- |
| Yersinia → `cgmlst_species_assignment` | BIGSdb-Pasteur cgMLST scheme 1 (500 loci) has no packaged CLI and no pinned allele caller | The Yersinia task returns "undecided". Classical 7-locus MLST does run, but its ST **cannot** separate the two: *pestis* is a clone nested inside *pseudotuberculosis*, so it is context only. VFDB and PlasmidFinder contain none of `pst`, `ypo2088`, `yihN` or `opgG`, so the abricate pre-screen adds no discriminating power either. |
| Listeria → `cgmlst_lin` | Same: no packaged CLI and no pinned allele caller for scheme 15 (1748 loci) | LisSero's serogroup and classical MLST both run, but the broad lineage does not, so the four *Listeria* lineage units are checked by the ANI panel instead. Classical MLST reports no clonal complex and no lineage, and no ST→CC→lineage table is pinned, so the field stays empty rather than being inferred. |
| C. botulinum → `genome_group` | No genome-wide group I–IV classifier is pinned | The group verdict for the two *C. botulinum* group units comes from the ANI panel, **never** from the toxin type: group I and group II strains can carry overlapping bont types. |

### Version the interpretation rule separately

A marker call is the product of three things that change independently: the
**software**, the **database**, and the **rule we apply to the output**. Recording only
the first two makes a changed call unattributable — you cannot tell a new BTyper3 build
from a new reading of the same output. So every ledger record carries all three, with
`interpretation_rule_version` as its own field. The rule is ours; it lives in
[`markers.py`](../bin/curation-triage/markers.py) as code.

A change to any of the three has to pass the **update gate** before it lands. The gate
is described in [running it](CURATION_TRIAGE_RUNNING.md); the design requirement it
implements is: a candidate artifact must have an immutable digest, a fixed validation
panel must be re-run, every call difference must be classified as a software, database,
threshold or interpretation-rule effect, and the decision record and changelog must be
updated. A call that changed while nothing it depends on moved is **unattributable**,
which is a failure and not a fifth category.

## 6. Facts the data forces on the design

These are properties of NCBI's data and of the tools, not of any single run. Each one
changed a design decision.

- **The type strain is keyed per assembly, not per species.** NCBI names a declared type
  strain per assembly, and it differs between the subspecies of one species. The eleven
  Kalamari *Salmonella* units all share species taxid 28901 and resolve to **seven
  different** type strains. Keying by species would compare most of them against the
  wrong genome.
- **"No type material" and "this pick is the type strain" are answers, not gaps.** 30
  picks are their own species' type strain, and 7 species have no type material at all —
  which is also why NCBI's taxonomy check on several flagged picks reads Inconclusive.
- **skani returns no row for a distant pair, rather than a low score.** Its default
  screen drops pairs below roughly 80% identity or 15% aligned fraction, and silence
  reads as "not measured" when it means the opposite. So the screen is widened and a
  screened-out pair is recorded as a **result in its own right**. On the committed cache
  4 of 157 rows are screened-out pairs.
- **A gap component returns "undecided"** rather than falling back to weaker evidence.
- **LPSN publishes no date column at all.** The valid-publication date is therefore
  **derived from the authority string**: cut at `emend.` or `corrig.`, then take the last
  four-digit year in what remains. That gives *E. coli* 1980, its Approved-Lists year,
  rather than the year of a later emendation. The rule and its measured effect are in
  [running it](CURATION_TRIAGE_RUNNING.md).
- **A digest proves sameness, not availability.** The AMRFinderPlus database
  `2026-05-15.1` is frozen by a per-file digest manifest, so anyone can prove offline
  that a copy is the same database. But `amrfinder -u` only ever fetches the latest
  database, so nothing can make NCBI serve that dated version again.

## 7. What the tool deliberately does not do

Naming the exclusions is part of the design, because each one protects the precision
target in §1.

| Not done | Why |
| --- | --- |
| Pass/fail grading of a genome | Kalamari's picks were never chosen by centrality, so a score would fail genomes that are correctly curated. |
| Automatic pull requests against `src/chromosomes.tsv` | A curation change is an SME decision. The tool writes a worklist and a ledger; a human makes the edit. |
| Scanning all of NCBI for species we do not track | Precision collapses. J2 is fed by curator intent and by newly named species in genera we already track. |
| Eukaryote pathogens | Tier 0 marks them and stops. NCBI's type-strain ANI report is prokaryote-only, so the primary signal does not exist for them. |
| Viruses and organelles | Excluded permanently by the Tier-0 guard — the comparison unit does not apply. |
| A learned embedding as a resolver | §2, idea 3. An embedding can only ever be a novelty radar. |
| Recomputing ANI for every pick every month | §2, idea 2. Tier 2 runs only on Tier-1 survivors and on the metadata-blind sub-species lineages. |

### The research track

A novelty radar is a research idea, and it never blocks anything. Sketch-space novelty
comes first: embed references and new deposits as Mash sketches, which CI already
builds, and let the nearest-neighbour distance flag genomes far from everything. UMAP is
for looking at, not for deciding with — keep the high-dimensional distance as the stable
metric, because UMAP is unstable on refit. **Bacformer** (gene-content, Apache-2.0) is
the one learned pilot worth trying, because it is the only candidate on an axis
orthogonal to the sketch baseline. DNABERT-S and Genos-m are redundant with sketches,
Evo2 has no native genome vector and is heavy, vPST is viral-only, and NT-v2 and ProkBERT
carry non-commercial licences. None of them is ever a resolver.

## 8. Reference: the committed tables

### `taxon_policy.tsv` — one row per genome-unit

Required columns, in order: `unit_id` (a readable id, `<repName>~<groupKey>`),
`scientificName`, `kalamari_taxid`, `species_taxid` (the resolved real NCBI one),
`accession_set`, `parent_assembly` (a real GCF/GCA, or a `strain:<taxid>` proxy),
`eligibility`, `regime`, `resolver`. Review columns follow: `species_name`,
`superkingdom`, `total_length_bp`, `n_replicons`, `flags`, `note`.

Useful `flags` values:

| Flag | Meaning |
| --- | --- |
| `shared_accession` | One physical sequence labelled two ways, for example `CP015575` |
| `label_conflict` | The unit's rows disagree about the name |
| `grouped_by_proxy` | A multi-replicon unit grouped without a real assembly — confirm by hand |
| `missing_record` | NCBI returned no record for a replicon |
| `species_from_organism` | The species came from the deposited genome, because the declared taxid was genus-level |
| `taxid_organism_mismatch` | The Kalamari taxid climbs to a different species than the genome's own NCBI organism — a mislabel surfaced for review. The declared species is kept and the organism is recorded in `note`. |
| `confounded` | Similarity cannot split this taxon; see §5 |

### `coverage_report.tsv` — one row per eligible pick

`unit_id`, `scientificName`, `source_lane`, `species_taxid`, `declared_species`,
`parent_assembly`, `ani_match`, `summary_match`, `in_refseq`, `taxonomy_check`,
`best_match_status`, `best_match_species`, `ncbi_current_species`, `version_status`,
`excluded_from_refseq`, `relation_to_type_material`, `seq_rel_date`, `confounded`,
`flags`.

`ani_match` and `summary_match` are typed as `exact`, `version-differs`, `taxon-only` or
`none`.

### `ledger.ndjson` — one self-contained JSON record per lead per run

```
schema, run_id, signal (J1a/J1b/J2/T2), lead_key, epoch, delta (new/recurring/resolved),
unit_id, scientificName, source_lane, candidate_accession, deposit_date,
actionable_flags[], flags[], evidence{...}, evidence_fingerprint, link, summary,
params{...}, tool_versions{...}
```

`lead_key` is `unit_id` for J1a; `unit_id#<candidate_base>` for J1b; `todo:<taxid>#<acc>`
or `lpsn:<species>` for J2; and `t2:<unit_id>#<check>` for Tier 2. `epoch` is the
tooling or schema epoch of the signal, not the run month: it moves only on a schema or
tool change.

The committed ledger holds **96 records**: the 21 standing J1a records from run
`2026-07`, plus 75 Tier-2 records from the `2026-08` replay (37 `ani`, 17 `panel`, 21
`marker`; 21 of them cross-linked to a Tier-1 lead by `confirms_lead`).

Tier 2 is a **new signal on the same ledger**, not a new ledger, so its leads flow
through the same delta, suppression and rendering as Tier-1 leads. Two rules keep the
two tiers from interfering. **Tier 2 marks nothing resolved** — its verdicts say nothing
about whether a Tier-1 lead is still actionable. And **context rides on a lead but can
never move it**: ECTyper's cached serotype sits on the *E. coli* marker lead under
`ectyper_*` keys, written *before* the verdict fields so context cannot overwrite a
call, and absent from the fingerprint so a changed serotype cannot re-open a dismissed
lead. The distinct key prefix carries weight: ShigEiFinder and ECTyper both emit a plain
`serotype`, and a merge that kept the bare key would drop ECTyper's value into
ShigEiFinder's slot.

## 9. What the tests protect

531 checks run offline in about 25 seconds, from 31 test files. They are not a coverage
exercise; each one pins a decision that would otherwise rot. The named regressions are:

- Synthetic taxids resolve to their parent species.
- The two *Yersinia enterocolitica* strains stay separate.
- The *Aliivibrio* multi-replicon genome concatenates into one unit.
- The *Arripis* fish mitogenome is caught by the eligibility guard, not by a hand list.
- The `CP015575` shared accession is flagged, not double-counted.
- Preflight fails on a deliberately unresolvable taxon.
- The committed panels match what `build_panels.py` regenerates.
- `tool_pins.tsv` and `marker_manifest.tsv` agree about every typer, so a bumped build
  cannot be recorded in one file and forgotten in the other.
- The documented mlst scheme count is the one the pinned binary reports.
- No committed file carries a developer's absolute filesystem path.

Ledger tests deserve their own note, because the ledger is a **seed file that the bot's
own commits grow**: a monthly run appends about 201 records (21 J1a, 59 J1b, 46 J2 and
75 Tier-2). A test that pinned the record count, the exact run ids, or "21 Tier-1
records" would turn the *blocking* test workflow red on the first bot commit. What is
asserted instead are rules that survive an append: every live lead renders as exactly one
row carrying its own record's run id, and a Tier-2 run appends only Tier-2 records and
leaves every Tier-1 record byte-identical. A companion test grows a copy of the ledger by
one bot month and re-checks those rules; another feeds a deliberately broken ledger, so
the checks cannot rot into vacuous ones.

## 10. Honest limits

This is the canonical list. Everything here is a real gap, not a caveat for form's sake.

1. **No workflow has ever run on a real GitHub runner.** Every workflow is lint-clean —
   `actionlint` 1.7.12 with `shellcheck` 0.11.0 report zero findings on the three
   fork-owned files — and step scripts have been extracted from the parsed YAML and run
   locally. But GitHub Actions itself has never executed one. The job to watch first is
   the AMRFinderPlus environment: its `amrfinder -u` branch is unexercised, because
   running it locally would destroy the pinned database.
2. **The AMRFinderPlus archive has no durable home.** The database is frozen and
   verifiable, but the only copy of the bytes is a single tarball on a single machine,
   and NCBI cannot serve that dated version again. Two copies in institutional storage,
   with their digests checked on a schedule, is the whole task. **This is the only open
   item whose cost, if skipped, is unrecoverable.**
3. **Three marker components have no runnable caller** (§5): an allele caller for the
   BIGSdb *Yersinia* and *Listeria* cgMLST schemes, and a genome-wide *C. botulinum*
   group I–IV classifier. Closing any of them means pinning a **new** tool, which changes
   the marker manifest, moves calls, and has to pass the update gate on its way in. That
   is a different kind of change from "pin what exists and prove it repeats", which is
   why the two are kept apart.
4. **The LPSN refresh is not unattended.** The committed cache was downloaded by a
   registered curator. The REST API at `https://api.lpsn.dsmz.de/` is the intended
   replacement and is not implemented: it needs the account credentials as a CI secret, a
   token exchange, and paging over roughly 34000 records.
5. **There is no automatic pull request.** An accepted lead becomes a *proposed edit* on
   the worklist, and a human applies it.
6. **The embedding novelty radar is research, not code** (§7).

## 11. The code

24 Python modules plus `setup_envs.sh` live in
[`bin/curation-triage/`](../bin/curation-triage/). Eleven are entry points; the rest are
libraries they import.

```
build_policy.py  policy.py  taxonomy.py  ncbi.py        # Tier 0
coverage_probe.py  coverage.py  ani_reports.py          # Tier 1 J1a + the coverage probe
tier1_metadata.py  leads.py  ledger.py  collate.py      # the Tier 1 net
datasets_summary.py  lpsn.py                            # the J1b and J2 feeders
tier2_confirm.py  tier2.py  markers.py                  # Tier 2 driver and pure logic
genomes.py  skani_runner.py  typers.py                  # Tier 2 cached, injectable layers
build_panels.py  build_type_strains.py                  # regenerate the committed Tier-2 inputs
refresh_tier2.py  setup_envs.sh  update_gate.py         # refresh, pin, gate
worklist_from_ledger.py                                 # the release asset's worklist
```

There is no single `run.sh`: each tier is a self-contained command-line tool.

**One structural rule runs through all of it: the network is quarantined in one module
per tier, and everything else is pure.** [`ncbi.py`](../bin/curation-triage/ncbi.py) is
the only network code in Tier 0, [`ani_reports.py`](../bin/curation-triage/ani_reports.py),
[`datasets_summary.py`](../bin/curation-triage/datasets_summary.py) and
[`lpsn.py`](../bin/curation-triage/lpsn.py) in Tier 1, and
[`genomes.py`](../bin/curation-triage/genomes.py),
[`skani_runner.py`](../bin/curation-triage/skani_runner.py) and
[`typers.py`](../bin/curation-triage/typers.py) in Tier 2. Each has an injectable
fetcher and a committed cache. That is why Tier 2 is five modules rather than one, and
it is why the whole test suite runs with no network and no conda environment.

Committed inputs live in [`src/curation-triage/`](../src/curation-triage/): the policy
table, the caches, the ledger, `reviews.tsv`, the marker and confounded manifests, 17
panel files, 11 environment lockfiles, `tool_pins.tsv` and `pin_changelog.tsv`. Genome
FASTA is **never** committed; it lives in the gitignored `.cache/curation-triage/genomes/`.
Only the small derived results are kept.

Default thresholds, all overridable: `species_ani=95`, `af_gate=0.5`, `derep=99.5`,
`epsilon≈0.2` ANI, `size_cap=50 Mbp`, and each typer's own percent-identity and coverage
defaults.

## 12. The upstream data fixes

`src/chromosomes.tsv` carried one `source` value misspelled `NBCI-REF` upstream, which
is why the corrected file has 179 `NCBI-REF` rows rather than 178. The Tier-0 preflight
also caught three taxid mislabels. All four fixes are applied here and captured in
[`chromosomes-data-fixes.patch`](../chromosomes-data-fixes.patch) for submission
upstream. Tier 0 never reads the `source` column, and the lane table in
[`coverage.py`](../bin/curation-triage/coverage.py) still maps `NBCI-REF` to `NCBI-REF`,
so an unpatched checkout keeps working.
