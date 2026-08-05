# Curation triage — how it works, and why to trust it

The reference for anyone reviewing the add-on or taking it over. Every rule here is
written with its reason: a rule without its reason gets deleted by the next person.

- To action a worklist: [for curators](CURATION_TRIAGE_CURATOR.md)
- To run, refresh or re-pin anything: [running it](CURATION_TRIAGE_RUNNING.md)

## The rules, in one screen

Six rules about the data. Break any of them and the tool gives a confident wrong answer.
The operational rules are in [running it](CURATION_TRIAGE_RUNNING.md).

| Rule | Why |
| --- | --- |
| Group replicons by their **parent assembly**, never by the Kalamari `taxid` column | Taxid `630` covers two different *Y. enterocolitica* strains (§3) |
| Join NCBI's ANI report on the **parent assembly**, never on replicon accessions | The report keeps only the latest assembly per taxon; replicon accessions find almost nothing (§4) |
| A type strain is keyed **per assembly**, not per species | It differs between subspecies: 11 *Salmonella* units → 7 type strains (§6) |
| A missing skani row means **"very distant"**, not "not measured" | skani drops a distant pair silently, so a screened-out pair is recorded as a result (§6) |
| The LPSN valid-publication date is **derived**, not published | LPSN has no date column at all; it comes from the authority string (§6) |
| A digest proves **sameness**, not availability | It cannot make NCBI serve a dated database again (§6) |

**Contents**

1. [What the tool is, and what it refuses to be](#1-what-the-tool-is-and-what-it-refuses-to-be) ·
2. [The three ideas](#2-the-three-ideas-that-shape-everything) ·
3. [The three tiers](#3-the-three-tiers) ·
4. [The measurement it rests on](#4-the-measurement-the-design-rests-on) ·
5. [The confounded taxa](#5-the-confounded-taxa-markers-not-similarity) ·
6. [Facts the data forces](#6-facts-the-data-forces-on-the-design) ·
7. [Deliberately out of scope](#7-what-the-tool-deliberately-does-not-do) ·
8. [Honest limits](#8-honest-limits) ·
9. [The code, the tables, the tests](#9-the-code-the-tables-and-the-tests)

---

## 1. What the tool is, and what it refuses to be

It gives Kalamari's subject-matter-expert (SME) curators a monthly **triage worklist** of
publicly deposited genomes worth examining. It turns curatorial choices made once into
choices watched continuously. Two jobs:

- **J1 — representative drift**: has our pick gone bad (mislabelled, withdrawn,
  reclassified), or has a better *trusted* genome appeared?
- **J2 — coverage gap**: is there a species we should add but do not have?

**It is informational, never pass/fail.** Being the statistical centre of a species was
never how these genomes were chosen, so the tool flags and ranks — it does not grade. SME
time is the scarce resource, so the whole design optimises for **precision**.

## 2. The three ideas that shape everything

**1. Kalamari picks on trusted provenance, not centrality.** The `source` column holds
179 NCBI-REF rows, 130 SME rows across three lanes, 10 FDA-ARGOS and 4 NCTC3000. These
genomes were chosen for completeness and trust, so the signal must track provenance and
correctness, not statistical centrality.

**2. Do not recompute what NCBI already computes.** NCBI publishes a roughly nightly
verdict for every prokaryotic assembly against its species' type strain, using ANI
(average nucleotide identity, a genome-wide percent-similarity measure). Reading *its
verdict on our own picks* is nearly free, and it catches what a similarity computation is
blind to: our pick was reclassified or withdrawn.

**3. For the confounded taxa, genome-wide similarity is the wrong ruler.** For *Shigella*,
anthrax, plague and the *Listeria*, *Salmonella* and botulinum splits, the meaningful
difference lives in specific genes — a virulence plasmid, a toxin gene. So the resolver is
**targeted marker-gene detection**, already validated in public-health laboratories. Not
ANI, and not a learned embedding.

The shape follows: **a cheap metadata net finds the few candidates → expensive
confirmation runs only on those → the confounded taxa are resolved by markers.**

## 3. The three tiers

### Tier 0 — the policy table

[`build_policy.py`](../bin/curation-triage/build_policy.py) turns `src/chromosomes.tsv`
plus NCBI taxonomy into the committed
[`taxon_policy.tsv`](../src/curation-triage/taxon_policy.tsv): one row per
**genome-unit**, 270 units from 323 replicon rows. Every later tier keys off it.

**Group replicons by parent assembly, never by the `taxid` column.** That column mixes
ranks and fails both ways: taxid `630` covers **two different *Yersinia enterocolitica*
strains**, which grouping by taxid would merge into one chimera, while the two chromosomes
of one *Aliivibrio* genome would look like two organisms. The real GCF/GCA assembly is the
ground truth — it merges one genome's replicons and keeps distinct genomes apart even when
they share a species (the eleven *Salmonella enterica* entries all resolve to taxid
`28901` and are eleven different assemblies). Where no assembly resolves, the fallback key
is scoped to the Kalamari identity, so it can never build a chimera either.

**Resolve a real NCBI species taxid.** The 13 synthetic Kalamari taxids
(`9000000`–`9000017`, such as the four *Listeria monocytogenes* lineages) climb through
the bundled taxonomy overlay to their real parent species. A **preflight hard-fails** if
any eligible unit lacks a real species-or-below taxon, so a synthetic, strain-level or
genus-level key can never pass silently.

**Decide eligibility with a programmatic guard, not a hand list.** A hand list misses what
nobody thought of.

| Class | Units | Meaning |
| --- | --- | --- |
| `eligible` | 249 | Bacteria and Archaea with a real species taxid — these are checked |
| `exclude:organelle` | 14 | Organelle-only submission, or over ~50 Mbp. The guard is what catches the *Arripis trutta* fish mitogenome. |
| `exclude:virus` | 4 | Viruses and viroids |
| `eukaryote` | 3 | Eukaryote pathogens — recorded, outside the pipeline |
| `review:genus-level` | 0 | Resolves no further than genus. Empty today; the class exists so such a unit surfaces rather than passing silently. |

The organelle guard fires only when **every** replicon is non-nuclear, so a eukaryote
nuclear assembly that bundles its mitochondrion (*Cryptococcus*: 14 chromosomes + 1
mitochondrion) is `eukaryote`, not `exclude:organelle`.

**Tag the confounded taxa — 21 units** — from
[`confounded_manifest.tsv`](../src/curation-triage/confounded_manifest.tsv), which sets a
`regime` (A/B) and a `resolver` (marker/ANI). A confounded unit **stays eligible**. It is
held out of the better-representative feed, not dropped: it still needs pick-health
monitoring, it just must not be nagged about type-strain divergence, which for it is
expected (§5).

### Tier 1 — the monthly metadata net

No genome downloads. It produces the worklist. The signals and the suppression rules are
described from the curator's side in [for curators](CURATION_TRIAGE_CURATOR.md); the
design points are:

- **J1a defers to NCBI's `taxonomy-check-status`**, not to raw best-match differences,
  because NCBI's verdict already weighs the best match and a raw difference does not. It
  respects NCBI's `approved-mismatch` marker.
- **J1b is source-gated.** For SME, FDA-ARGOS and NCTC3000 units, divergence from a type
  strain is a deliberate curatorial choice, so J1b never fires on it. Those picks are
  covered by J1a.
- **J2 is fed by intent, not by scanning.** The wishlist is explicit curator intent, the
  highest-precision source there is; LPSN is limited to genera already tracked.

### Tier 2 — targeted confirmation

Only for units that need it. It is **confirmation and sub-species residue, never
coverage**: §4 shows Tier 1 already sees every pick. Two trigger groups make 41 units —
the 21 Tier-1 survivors plus the 21 confounded units, overlapping in one. Three checks:
`ani` (survivors and regime B), `panel` (regime B, 17 units), `marker` (all confounded).

**The aligned fraction is a veto, not a tiebreaker.** An ANI computed over a sliver of a
genome is not weak evidence, it is no evidence. So the gate checks aligned fraction first
and takes the **smaller** of the two sides — a small genome fully contained in a larger one
aligns well on its own side while covering little of the other, which is exactly the case
the gate exists to catch. Below the gate the verdict is "cannot decide", never "different
species"; within `epsilon` of the boundary it is "too close to call".

**The sub-species panels compare Kalamari against itself.** The 13 synthetic-taxid units
are metadata-blind — NCBI knows nothing about a lineage only Kalamari defines — so their
panel members are Kalamari's own sibling units, generated into 17
[`panels/*.acclist`](../src/curation-triage/panels/) files with a test asserting the
committed files match. That answers two questions offline without inventing reference
accessions nobody has verified: **is the panel still a panel?** (a sibling at or above the
99.5% dereplication line means two lineage references are effectively the same genome — a
curation defect nothing else can see), and **where does a genome sit?** Because the unit's
own lineage is absent from its panel, "nearest" is a neighbourhood, not a lineage
assignment, and the code says so rather than over-claiming.

## 4. The measurement the design rests on

J1a works only if NCBI publishes a live verdict for *our exact assembly*. That is an
assumption about someone else's data, so it is measured.
[`coverage_probe.py`](../bin/curation-triage/coverage_probe.py) joins every eligible pick
to NCBI's reports and writes
[`coverage_report.tsv`](../src/curation-triage/coverage_report.tsv).

**249 of 249 eligible picks (100%) have a live row in NCBI's ANI report.**

| Coverage — a live row for *our* assembly | Count |
| --- | --- |
| ANI report — at the exact version we record | 246 |
| ANI report — superseded to a newer version | 3 |
| `assembly_summary` (refseq ∪ genbank) | 249/249 |
| currently live in RefSeq | 228/249 |

It holds across every breakdown: source lane, release year including the pre-2010 cohort,
and GCF versus GCA. The 21 picks not live in RefSeq are all in
`assembly_summary_genbank` with an `excluded_from_refseq` reason — a signal in its own
right, not a gap.

**Thin coverage is the obvious risk, and the join key is why it does not bite. Never join
this report on replicon accessions; join on the parent assembly.** The report keeps only
the latest assembly per taxon, and many picks are decades old. Joining on the INSDC
replicon accession a pick records — `AE002098`, from the year 2000 — finds almost nothing.
Joining on the real GCF/GCA parent assembly finds everything, because NCBI keeps a row for
those regardless of age.

Two NCBI conveniences make it cheap, and both are properties of the files: the ANI report
carries both the GCA and GCF accession per row, and each `assembly_summary` row carries
its twin in `gbrs_paired_asm`, so GCF↔GCA equivalence is free. Accessions match on their
version-stripped base and the held version is then compared, which is what surfaces a
superseded pick.

**What it settles.** Tier 1 carries the drift and mislabel half of the job on its own. No
self-run type-strain ANI is needed *for coverage*, which frees Tier 2's own ANI for what
metadata is blind to: the sub-species lineages, and confirmation on the few survivors.

## 5. The confounded taxa: markers, not similarity

The 21 units split into two regimes, and only regime A is a case of "ANI fails".

| Regime | Taxa | Resolver |
| --- | --- | --- |
| **A — ANI cannot split them** | *Shigella* ↔ *E. coli* · *B. anthracis* ↔ the *B. cereus* group · *Y. pestis* ↔ *Y. pseudotuberculosis* | **Marker genes**, which override ANI |
| **B — ANI resolves or over-splits** | *L. monocytogenes* lineages · *Salmonella* subspecies and serovars · *C. botulinum* groups | An ANI panel places the genome; a typer supplies the label |

In regime A the two organisms share a near-identical backbone, and what separates them is
a small mobile or accessory gene set that a genome-wide average cannot see. In regime B
core divergence is real and a sketch captures it; the residual problem is *labelling*.

**The regime-A override is structural.** A regime-A unit gets no independent ANI verdict
at all: the number is folded into the marker lead as context, and the marker call is the
verdict. There is no code path where an ANI number can contradict a regime-A marker call,
because only one of them is ever a verdict. If markers cannot decide, the answer is
"undecided" — it does **not** fall back to ANI, because for these taxa the ANI answer was
never admissible.

For regime B the panel decides membership and the typer supplies the label. A wrong label
on a correctly-placed genome is still a defect, so a refuting label is reported.

### The marker manifest

Wrap validated tools; do not hand-curate gene lists. Three layers feed the calls: curated
gene and virulence databases (VFDB, the NCBI catalog behind AMRFinderPlus, CGE
VirulenceFinder and PlasmidFinder); validated in-silico typers that bundle their own
database and interpretation; and allele or scheme databases for the labelling splits
(PubMLST, Institut Pasteur BIGSdb).

[`marker_manifest.tsv`](../src/curation-triage/marker_manifest.tsv) holds **15 components
over 8 identification tasks**. A *task* is the question ("is this anthracis?"); a
*component* is one tool run contributing evidence to it. Twelve are pinned and runnable;
three are gaps.

| Task | Resolver | Defining markers |
| --- | --- | --- |
| `shigella_eiec_vs_escherichia_coli` | ShigEiFinder 1.3.4 | ipaH, cluster-specific genes, O/H |
| `bacillus_anthracis_vs_b_cereus_group` | BTyper3 3.4.0 **plus our rule** | pXO1: pagA, lef, cya · pXO2: capA, capB, capC |
| `yersinia_pestis_vs_y_pseudotuberculosis` | **gap**; `mlst 2.23.0` is context only | ypo2088, yihN, opgG, caf1, pst |
| `listeria_monocytogenes_serogroup_lineage_cc` | LisSero 0.4.10 + `mlst`; cgMLST/LIN is a **gap** | Doumith pattern: prs, lmo0737, ORF2819, ORF2110, lmo1118 |
| `salmonella_subspecies_and_serovar` | SISTR 1.1.3 + SeqSero2 1.3.2 | O (wzx/wzy), H (fliC/fljB), cgMLST330 |
| `clostridium_botulinum_group_and_bont_type` | AMRFinderPlus 4.2.7, DB `2026-05-15.1` | bont/A–G |
| `escherichia_coli_serotype_pathotype` | ECTyper 2.0.0 — context only, never a verdict | O/H alleles, stx subtype |
| `generic_marker_fallback` | abricate 1.4.0; AMRFinderPlus 4.2.7 | any curated virulence or plasmid gene |

A marker call is the product of three things that move independently: the **software**, the
**database**, and **the rule we apply to the output**. Recording only the first two makes a
changed call unattributable — you cannot tell a new BTyper3 build from a new reading of the
same output. So every ledger record carries all three, with `interpretation_rule_version`
as its own field. The rule is ours, and lives in
[`markers.py`](../bin/curation-triage/markers.py) as code. Any change to the three passes
the [update gate](CURATION_TRIAGE_RUNNING.md#7-the-update-gate).

### Why these tools, and not the obvious ones

Four tasks reject the tool an experienced reader would reach for first. Three reject
`abricate` as a **classifier**. Each is a rule, not a taste.

- **Yersinia: not abricate over VFDB/PlasmidFinder.** Plasmid presence alone misclassifies
  both directions — plasmid-cured *Y. pestis* exists, and so does plasmid-bearing
  *Y. pseudotuberculosis*. **A plasmid hit must never make the species call.** abricate
  stays as a cheap pre-screen.
- **C. botulinum toxin: not abricate with a bont gene set.** abricate's databases carry no
  version of their own, so a changed call could not be attributed to anything.
  AMRFinderPlus has a first-class dated database version, so a changed call reads as a
  database update rather than as biology.
- **Listeria: do not merge serogroup, MLST and cgMLST into one "lineage" field.** They
  answer **different questions**; merging destroys the distinction. Outputs stay separate.
- **Bacillus: do not trust BTyper3's species call.** It gives group taxonomy and a marker
  profile, not an anthracis verdict. Turning that into "is this anthracis" is our
  judgement, so it is versioned as ours.

Two more rules, each because breaking it gives a confident wrong answer.
**Plasmid-cured *B. anthracis* is still *B. anthracis***, and a *B. cereus* carrying
anthrax plasmids is never labelled plain *B. cereus*. And **"anthracis" is matched as a
species epithet, never as a substring** — *Bacillus paranthracis* contains the string, is a
different organism, and is what NCBI's best ANI match reports for one of Kalamari's own
picks.

### The three gaps

Three components have **no runnable, pinnable caller**, so the manifest records them as
`gap` rather than as pins nothing can execute. A gap makes its task return
**"undecided"** — deliberately, because the alternative is falling back to weaker evidence.

- **Yersinia cgMLST** (scheme 1, 500 loci) has no packaged CLI and no pinned allele caller.
  Classical 7-locus MLST runs but **cannot** separate the two species — *pestis* is a clone
  nested inside *pseudotuberculosis* — so it is context only. The abricate pre-screen adds
  nothing either: VFDB and PlasmidFinder contain none of `pst`, `ypo2088`, `yihN`, `opgG`.
- **Listeria cgMLST/LIN** (scheme 15, 1748 loci): same reason. Serogroup and classical
  MLST still run, so the four lineage units fall to the ANI panel. No ST→CC→lineage table
  is pinned, so the field stays empty rather than being inferred.
- **C. botulinum genome group**: no genome-wide group I–IV classifier is pinned, so the
  group comes from the ANI panel and **never** from the toxin type — group I and group II
  strains can carry overlapping bont types.

## 6. Facts the data forces on the design

Properties of NCBI's data and of the tools, not of any one run. Each changed a decision.

- **The type strain is keyed per assembly, not per species.** NCBI names one per assembly
  and it differs between subspecies: the eleven *Salmonella* units share species taxid
  28901 and resolve to **seven different** type strains. Keying by species would compare
  most of them against the wrong genome.
- **"No type material" and "this pick is the type strain" are answers, not gaps.** 30
  picks are their own species' type strain; 7 species have no type material at all. The
  second group is also why NCBI's check on several flagged picks reads Inconclusive.
- **skani returns no row for a distant pair, rather than a low score.** Its default screen
  drops pairs below roughly 80% identity or 15% aligned fraction, and silence reads as
  "not measured" when it means the opposite. So the screen is widened and a screened-out
  pair is recorded as a **result in its own right** — 4 of the 157 committed rows.
- **A gap component returns "undecided"** rather than falling back to weaker evidence.
- **LPSN publishes no date column at all**, so the valid-publication date the J2 feeder
  uses is *derived* from the authority string. It is an authority year, not a full date.
  The rule is in [running it](CURATION_TRIAGE_RUNNING.md#3-the-committed-caches).
- **A digest proves sameness, not availability.** The AMRFinderPlus database
  `2026-05-15.1` is frozen by a per-file digest manifest, so anyone can prove offline that
  a copy is the same database. But `amrfinder -u` only ever fetches the latest, so nothing
  can make NCBI serve that dated version again.

## 7. What the tool deliberately does not do

Each exclusion protects the precision target in §1.

| Not done | Why |
| --- | --- |
| Pass/fail grading of a genome | The picks were never chosen by centrality, so a score would fail genomes that are correctly curated |
| Automatic pull requests | A curation change is an SME decision |
| Scanning all of NCBI for untracked species | Precision collapses; J2 is fed by intent instead |
| Eukaryote pathogens | NCBI's type-strain ANI report is prokaryote-only, so the primary signal does not exist |
| Viruses and organelles | The comparison unit does not apply |
| A learned embedding as a resolver | §2, idea 3 — it can only ever be a novelty radar |
| Recomputing ANI for every pick monthly | §2, idea 2 |

**The research track.** A novelty radar never blocks anything. Sketch-space novelty comes
first: embed references and new deposits as Mash sketches, which CI already builds, and
let nearest-neighbour distance flag genomes far from everything. Keep the
high-dimensional distance as the stable metric — UMAP is for looking at, not deciding
with, because it is unstable on refit. **Bacformer** (gene-content, Apache-2.0) is the one
learned pilot worth trying, being the only candidate on an axis orthogonal to the sketch
baseline; DNABERT-S and Genos-m are redundant with sketches, Evo2 has no native genome
vector, vPST is viral-only, and NT-v2 and ProkBERT carry non-commercial licences. None of
them is ever a resolver.

## 8. Honest limits

The canonical list. All real gaps, not caveats for form's sake.

1. **No workflow has ever run on a real GitHub runner.** They are lint-clean and their
   step scripts have been run locally, but GitHub Actions itself has never executed one.
   Watch the AMRFinderPlus job first: its `amrfinder -u` branch is unexercised, because
   running it locally would destroy the pinned database.
2. **The AMRFinderPlus archive has no durable home.** The database is frozen and
   verifiable, but the only copy of the bytes is a single tarball on a single machine, and
   NCBI cannot serve that dated version again. Two copies in institutional storage with
   scheduled digest checks is the whole task. **This is the only open item whose cost, if
   skipped, is unrecoverable.**
3. **Three marker components have no runnable caller** (§5). Closing any means pinning a
   **new** tool, which changes the manifest, moves calls, and must pass the update gate —
   a different kind of change from "pin what exists and prove it repeats".
4. **The LPSN refresh is not unattended.** The committed cache was downloaded by a
   registered curator; the REST API path needs credentials as a CI secret, a token
   exchange and paging over roughly 34000 records.
5. **There is no automatic pull request.** An accepted lead becomes a proposed edit, and a
   human applies it.
6. **The embedding novelty radar is research, not code** (§7).

## 9. The code, the tables, and the tests

24 Python modules plus `setup_envs.sh` live in
[`bin/curation-triage/`](../bin/curation-triage/); eleven are entry points, the rest are
libraries. There is no single `run.sh` — each tier is a self-contained command-line tool.

**One structural rule runs through all of it: the network is quarantined, and everything
else is pure.**

| Tier | The only modules that reach the network |
| --- | --- |
| 0 | [`ncbi.py`](../bin/curation-triage/ncbi.py) |
| 1 | [`ani_reports.py`](../bin/curation-triage/ani_reports.py), [`datasets_summary.py`](../bin/curation-triage/datasets_summary.py), [`lpsn.py`](../bin/curation-triage/lpsn.py) |
| 2 | [`genomes.py`](../bin/curation-triage/genomes.py), [`skani_runner.py`](../bin/curation-triage/skani_runner.py), [`typers.py`](../bin/curation-triage/typers.py) |

Each has an injectable fetcher and a committed cache. That is why Tier 2 is five modules
rather than one, and why the whole test suite runs with no network and no conda
environment. Genome FASTA is **never** committed; only the small derived results are.
Default thresholds, all overridable: `species_ani=95`, `af_gate=0.5`, `derep=99.5`,
`epsilon≈0.2`, `size_cap=50 Mbp`.

**The committed tables.** [`taxon_policy.tsv`](../src/curation-triage/taxon_policy.tsv)
carries one row per genome-unit; its `flags` column is where Tier 0 surfaces trouble —
`shared_accession` (one physical sequence labelled twice, e.g. `CP015575`),
`grouped_by_proxy` (grouped without a real assembly, confirm by hand),
`taxid_organism_mismatch` (the Kalamari taxid climbs to a different species than the
deposited organism — the declared species is **kept** and the organism recorded in `note`),
and `confounded`. [`coverage_report.tsv`](../src/curation-triage/coverage_report.tsv) has
one row per eligible pick; its `ani_match` and `summary_match` columns (`exact`,
`version-differs`, `taxon-only`, `none`) are what the 100% figure counts.

[`ledger.ndjson`](../src/curation-triage/ledger.ndjson) is append-only, one self-contained
JSON record per lead per run, holding the signal, `lead_key`, `epoch`, delta,
`evidence_fingerprint`, the evidence, and the tool versions and parameters that produced
it. `lead_key` is `unit_id` for J1a, `unit_id#<candidate>` for J1b, `todo:`/`lpsn:` for
J2, `t2:<unit_id>#<check>` for Tier 2. The committed file holds **96 records**: 21 J1a
from run `2026-07`, plus 75 Tier-2 from `2026-08` (37 `ani`, 17 `panel`, 21 `marker`; 21
cross-linked to a Tier-1 lead).

Tier 2 is a **new signal on the same ledger**, not a new ledger. Two rules keep the tiers
from interfering. **Tier 2 marks nothing resolved** — its verdicts say nothing about
whether a Tier-1 lead is still actionable. And **context rides on a lead but can never
move it**: ECTyper's serotype sits under `ectyper_*` keys, written *before* the verdict
fields so it cannot overwrite a call, and absent from the fingerprint so it cannot re-open
a dismissed lead. The prefix carries weight — ShigEiFinder and ECTyper both emit a plain
`serotype`, so a merge keeping the bare key would drop one into the other's slot.

**The tests.** 531 checks run offline in about 25 seconds. Each pins a decision that would
otherwise rot: synthetic taxids resolve; the two *Y. enterocolitica* strains stay
separate; *Aliivibrio* concatenates into one unit; the *Arripis* mitogenome is caught by
the guard; `CP015575` is flagged, not double-counted; preflight fails on an unresolvable
taxon; the committed panels match what `build_panels.py` regenerates; `tool_pins.tsv` and
`marker_manifest.tsv` agree about every typer; the documented mlst scheme count matches
the pinned binary; and no committed file carries an absolute path.

Ledger tests need care, because the ledger is a **seed file the bot's own commits grow** —
a monthly run appends about 201 records (21 J1a, 59 J1b, 46 J2, 75 Tier-2). A test pinning
the record count or the run ids would turn the *blocking* workflow red on the first bot
commit. What is asserted instead survives an append: every live lead renders as exactly
one row carrying its own record's run id, and a Tier-2 run appends only Tier-2 records and
leaves every Tier-1 record byte-identical.

**Upstream data fixes.** `src/chromosomes.tsv` carried one `source` value misspelled
`NBCI-REF`, which is why the corrected file has 179 `NCBI-REF` rows and not 178; the
preflight also caught three taxid mislabels. All four are captured in
[`chromosomes-data-fixes.patch`](../chromosomes-data-fixes.patch) for submission upstream.
Tier 0 never reads the `source` column, and
[`coverage.py`](../bin/curation-triage/coverage.py) still maps `NBCI-REF` to `NCBI-REF`,
so an unpatched checkout keeps working.
