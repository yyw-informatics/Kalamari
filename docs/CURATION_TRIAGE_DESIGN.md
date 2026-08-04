# Curation triage — design

*The authoritative design document for the curation-triage add-on: why it works the way it
does. Start with §1–§3 for the reasoning, §4 for the architecture, §11 for what is built, and §12 for
what is not.*

---

## 1. Purpose

Give Kalamari's subject-matter-expert (SME) curators an automated, monthly **triage worklist of newly
deposited public genomes worth examining** — turning the database's static, one-time curatorial choices
into continuously monitored ones. Two curation jobs:

- **J1 — representative drift** (for a species we already track): has our pick gone bad (mislabeled,
  withdrawn, reclassified), or has a better *trusted* genome appeared?
- **J2 — coverage gap**: is there a species/lineage we should add but don't?

**It is informational, never pass/fail.** Being the "centroid" was never how these genomes were chosen,
so the tool flags and ranks for human review — it does not grade. SME time is the scarce resource, so the
entire design optimizes for **precision**: few, high-value flags.

## 2. Design principles (the three ideas that shape everything)

1. **Kalamari picks on trusted provenance, not centrality.** The `source` column (179 NCBI-REF, 130
   SME across three SME lanes, 10 FDA-ARGOS, 4 NCTC3000) shows genomes were chosen for completeness +
   trust. So the triage signal must track provenance/correctness, not statistical centrality.
2. **Don't recompute what NCBI already computes.** NCBI publishes a ~nightly ANI-vs-type-strain verdict
   for every prokaryotic assembly (ANI = average nucleotide identity, a genome-wide percent-similarity
   measure). Reading *its verdict on our own picks* is nearly free and catches the "our pick was
   reclassified/withdrawn" case that a centroid computation is structurally blind to.
3. **For the confounded taxa, genome-wide similarity is the wrong ruler.** For the bacteria the paper
   highlights (Shigella, anthrax, plague, the Listeria/Salmonella/botulinum splits), the meaningful
   distinction lives in specific genes (a virulence plasmid, a toxin gene), so the resolver is
   **targeted marker-gene detection** — interpretable and already validated in public-health labs —
   not ANI and not a learned embedding.

Net architecture: **a cheap metadata net finds the few candidates → expensive confirmation runs only on
those → the confounded taxa are resolved by markers, not similarity.**

## 3. Locked decisions

- **Scope:** v1 = J1 (rep drift) + J2 metadata feeders. Outward missing-species scanning of all of NCBI
  and the learned-embedding novelty radar are a research track, not v1 (§6, §7).
- **Comparison unit:** compare Kalamari's curated entry *as shipped* (its nuccore replicon set), using
  ANI (size-robust) with **aligned fraction as a QC gate**; resolve the parent assembly only to dedup
  the pick's own genome out of any candidate panel.
- **Primary J1 signal:** NCBI's metadata type-strain verdict (the "pick-health sentinel"). This rests on
  one assumption — that NCBI publishes a live verdict row for our picks — and that assumption is
  measured, not assumed: **249/249 eligible picks have one** (§11).
- **Confounded taxa:** resolved by **marker genes**, which **override** the ANI verdict. Embeddings are
  never the resolver.
- **Embedding layer:** secondary, out of v1. Sketch-space (Mash) novelty first; **Bacformer**
  (gene-content, Apache-2.0) is the one learned pilot if pursued — not DNABERT-S — because it's the only
  candidate on an axis orthogonal to the sketch baseline.
- **Handoff:** per-run artifact + `$GITHUB_STEP_SUMMARY` worklist + committed append-only `ledger.ndjson`
  with a `reviews.tsv` SME decision loop. No auto-PRs in v1.
- **Delivery:** PR-ready add-on inside the Kalamari repo layout, rather than a standalone repo, so it
  stays in step with the curated tables (`src/chromosomes.tsv`, `src/taxonomy/`) it reads.

## 4. Architecture — the tiers

### Tier 0 — RESOLVE (committed policy table)
`build_policy.py` reads `src/chromosomes.tsv` + bundled taxonomy → committed `taxon_policy.tsv`, one row
per genome-unit: `unit_id, scientificName, kalamari_taxid, species_taxid (resolved), accession_set,
parent_assembly, eligibility, regime, resolver, species_name, superkingdom, total_length_bp,
n_replicons, flags, note`.
- **Group replicons by parent assembly**, never the taxid column (avoids the Y. enterocolitica 2-strain
  chimera and the CP015575 shared-accession trap; correctly merges Vibrio/Cryptococcus multi-replicon).
- **Resolve a real NCBI species taxid** (climb taxonomy; map the 13 synthetic `9000000`–`9000017` taxids
  and strain-rank taxids via the `parent` column). **Preflight hard-fails** on any unresolvable eligible unit.
- **Eligibility via programmatic guard** (not a hand list): Bacteria/Archaea → `eligible`; molecule type
  = mito/plastid/segment OR assembly > ~50 Mbp → `exclude:organelle` (auto-catches Arripis); viruses →
  `exclude:virus`; eukaryote pathogens → `eukaryote` (not in v1).
- **The confounded set stays `eligible`** and is marked by a `confounded` value in the `flags` column plus
  a `regime` (A/B) and a `resolver` (marker/ANI). It is *held out of J1b*, not dropped: these units still
  need pick-health monitoring, they just must not be nagged about type-strain divergence.

### Tier 1 — cheap metadata net (monthly, no genome downloads)
- **J1a · pick-health sentinel (primary signal).** Bulk-diff our ~270 parent assemblies against NCBI's
  `ASSEMBLY_REPORTS/ANI_report_prokaryotes.txt` + `assembly_summary_{refseq,genbank}.txt`. **Flag** when a
  pick's row shows: taxonomy-check → Failed/Inconclusive, best-match-species ≠ declared species,
  `excluded_from_refseq` populated, assembly superseded/suppressed, or taxid reassigned. Respect NCBI's
  `approved-mismatch` (never flag those). Catches the drift/mislabel/withdrawal half of J1.
- **J1b · source-gated better-rep feed.** Delta query per species (`datasets summary … --released-after`).
  For **NCBI-REF** units: flag when NCBI's reference designation changes. For **SME** units: divergence
  from type strain is *deliberate* — never nag "a type strain exists"; only surface if the pick itself
  degrades (via J1a). Confounded units are suppressed.
  **The delta is the ledger's, not the query's.** `--released-after` is a fetch-time window only: it makes
  the rebuild ask NCBI for less, it does not filter leads, and because the rebuild REPLACES the cache, a
  window leaves behind a cache the unwindowed monthly run cannot use. The committed cache is unwindowed
  and CI passes no window; see `docs/CURATION_TRIAGE_TIER1.md`.
- **J2 · coverage feeders.** New type-strain genomes matching `chromosomes-todo.tsv` (SME intent =
  highest precision) + newly validly-published species (LPSN) within tracked pathogen genera.
- **Ledger + reviews + worklist.** Append `ledger.ndjson` per taxon per run; read `reviews.tsv`
  (accept/reject/defer) to suppress dismissed leads; write the actionable rows to `$GITHUB_STEP_SUMMARY`
  + artifact.

### Tier 2 — targeted confirmation (only flagged survivors)
- **skani ANI** on candidate-vs-pick-vs-type-strain (aligned-fraction gate) to attach a precise number.
- **Sub-species mini-panels:** for the synthetic-taxid lineages, metadata is blind → a small
  ANI-to-lineage-reference panel is the only signal.
- **Marker-gene resolvers for confounded taxa (§5)** — the marker call **overrides** ANI.

### Side-track — J2 embedding novelty radar (research only, §7)
Sketch-space novelty (reuses the CI Mash sketch) + optional Bacformer probe. Never a resolver; a
"is this something we have no marker for" radar only.

## 5. Confounded-taxa handling — markers, not similarity

These split into two regimes; only Regime A is "ANI fails":

| Regime | Taxa | Why ANI is wrong | Resolver |
|---|---|---|---|
| **A — ANI can't split** | Shigella↔E. coli; B. anthracis↔B. cereus grp; Y. pestis↔Y. pseudotuberculosis | Near-identical backbone; signal is a small mobile/accessory gene set | **Marker genes** (below) — override ANI |
| **B — ANI already resolves / over-splits** | Listeria monocytogenes lineages; Salmonella subspecies/serovars; C. botulinum groups | Core-genome divergence sketches already capture; residual problem is *labeling* | ANI/sketch cluster + **custom-taxonomy reference genome supplies the name**; serovar/lineage typer for the label |

**How marker info is obtained — wrap validated tools, don't hand-curate.** Three layers:
1. **Curated gene/virulence databases** (the reference sequences): VFDB, NCBI reference gene catalog
   (AMRFinderPlus), CGE VirulenceFinder / PlasmidFinder.
2. **Validated in-silico typers** (bundle their own DB + interpretation): SISTR & SeqSero2 (Salmonella
   serovar), ShigEiFinder / ShigaTyper / ShigaPass (Shigella vs EIEC), ECTyper (E. coli), BTyper3
   (Bacillus cereus group incl. anthrax-marker profile).
3. **Allele/scheme databases** for the labeling splits: PubMLST / Institut Pasteur BIGSdb (Listeria
   cgMLST), EnteroBase.

**Detection mechanism:** BLAST/k-mer presence-absence + allele matching against the curated DB at fixed
%identity/%coverage; `abricate` is the generic wrapper for custom markers (e.g. a specific plasmid gene)
over VFDB/PlasmidFinder/VirulenceFinder/NCBI.

**In our tool:** a committed **marker manifest** (`src/curation-triage/marker_manifest.tsv`) maps each
confounded taxon → `{resolver_tool, database, database_version, defining_marker_genes}`; Tier 2 invokes
the typer, records the interpreted call, and **pins both the tool version AND the database version** in
the ledger (a DB update = a tooling epoch, so a changed call reads as a DB change, not biology).

The manifest holds **15 components over 8 identification tasks**. A *task* is the question ("is this
anthracis?"); a *component* is one tool run that contributes evidence to it. Twelve components are
pinned and runnable; three are `gap` (below).

| Taxon | Regime | Resolver tool (install) | Defining markers | DB pinning |
|---|---|---|---|---|
| Shigella / EIEC | A | **ShigEiFinder** (`bioconda::shigeifinder` 1.3.4) | ipaH + cluster-specific + O/H (wzx/wzy, fliC) | DB bundled → pin tool version (class a) |
| Salmonella serovar/subsp | B | **SISTR** (`bioconda::sistr_cmd` 1.1.3) + **SeqSero2** (`bioconda::seqsero2` 1.3.2) | O (wzx/wzy), H (fliC/fljB), cgMLST330 | SISTR DB on **Zenodo → pin DOI**; SeqSero2 bundled → pin tool (b / a) |
| E. coli serotype/pathotype | — | **ECTyper** (`bioconda::ectyper` 2.0.0) | O/H alleles, stx subtype | alleles bundled (DB ver in report); MASH species sketch on **Zenodo → pin** (a / b) |
| B. anthracis vs B. cereus grp | A | **BTyper3** (`bioconda::btyper3` **3.4.0-0**) — it supplies group taxonomy + a marker profile, *not* a standalone anthracis call | pXO1: pagA, lef, cya · pXO2: capA, capB, capC | bundled `seq_*_db` + a local PubMLST copy (downloaded 2023-07-07) → exact build or container digest + bundled-DB SHA-256 (class a) |
| Y. pestis vs Y. pseudotuberculosis | A | **BIGSdb-Pasteur *Yersinia* cgMLST scheme 1** (500 loci, `pubmlst_yersinia_seqdef`) — **`gap`**, see below. `mlst` 2.23.0 (classical 7-locus) runs as context; abricate is a pre-screen only, never the classifier | chromosome: ypo2088, yihN (*pestis*), opgG (*pseudotuberculosis*) · pMT1: caf1 · pPCP1: pst (pla only as corroboration, never alone) · lcrV = pCD1/pYV context, not a species discriminator | complete frozen scheme snapshot + SHA-256; pin the allele caller and its thresholds (class e→d) |
| Listeria monocytogenes lineage | B | **LisSero** (`bioconda::lissero` **0.4.10-0**) for serogroup + `mlst` 2.23.0 classical MLST for ST; BIGSdb-Lm cgMLST1748_v2/LIN (scheme 15, 1748 loci) for lineage — **`gap`**, see below. **Report these as separate, un-merged fields** | Doumith pattern: prs, lmo0737, ORF2819, ORF2110, lmo1118 | LisSero bundled DB (2020-06-04) → exact build + bundled-DB SHA-256 + thresholds (class a); schemes → independent full snapshot hashes, pinned allele caller, local regression panel (class e→d) |
| C. botulinum group/toxin | A/B | **NCBI AMRFinderPlus `--plus`** (`bioconda::ncbi-amrfinderplus` **4.2.7-0**) for the toxin call — *not* abricate. Genomic group I–IV needs a separate genome-wide classifier — **`gap`**, see below | bont/A–G, with subtype where the evidence supports it | dated DB **`2026-05-15.1`**: `amrfinder -V` + dated-DB-directory SHA-256 (class c+d) |
| *generic fallback* | — | **abricate** (`bioconda::abricate` 1.4.0) over VFDB/PlasmidFinder; **AMRFinderPlus** (`bioconda::ncbi-amrfinderplus` 4.2.7) | any curated VF/plasmid gene | AMRFinder **dated DB (`YYYY-MM-DD.N`) = best anchor**; abricate → freeze `$ABRICATE_DB` (VFDB 4592 seqs, PlasmidFinder 488 seqs, both dated 2026-Apr-3) (c+d / d) |

### Why these tools and not the obvious ones

Four rows reject the tool an experienced reader would reach for first. Three of the four reject
`abricate` as a **classifier**. Each rejection is a rule, not a preference:

| Task | The obvious choice | Why it is wrong | What is used instead |
|---|---|---|---|
| *Y. pestis* vs *pseudotuberculosis* | abricate over VFDB/PlasmidFinder | **Plasmid presence alone misclassifies plasmid-cured *pestis* and plasmid-bearing *pseudotuberculosis*.** A plasmid hit must never make the species call. | BIGSdb-Pasteur cgMLST scheme 1 — a `gap`, so the task returns "undecided"; abricate stays as a cheap pre-screen |
| *C. botulinum* toxin call | abricate + a bont gene set | abricate's databases carry no version of their own, so a changed call cannot be attributed. | AMRFinderPlus `--plus`, the one tool with a first-class dated DB version, so a changed call reads as a DB update rather than as biology |
| *L. monocytogenes* lineage | merge serogroup + MLST + cgMLST into one "lineage" field | Serogroup, classical MLST ST/CC and cgMLST/LIN answer **different questions**; merging them destroys the distinction. | the same tools, outputs kept as separate fields |
| *B. anthracis* vs *B. cereus* | trust BTyper3's species call | BTyper3 gives **group taxonomy and a marker profile**, not an anthracis verdict. Turning that into "is this anthracis" is our judgement. | BTyper3 plus our own interpretation rule, versioned as ours |

Two properties of BTyper3 3.4 that its output depends on: it reports the **revised** *B. cereus* group
nomenclature (`species(ANI)` is the genomospecies, e.g. *mosaicus*, and the familiar name appears in
`subspecies(ANI)`), and it **does not report `rpoB`** — so `rpoB` is not in its marker list even though
it is a published SNP discriminator. It writes to a file, not stdout.

### The three `gap` components

Three components have **no runnable, pinnable caller**, so they are recorded as `status = gap` in the
manifest rather than left as pins that nothing can execute. A `gap` component makes its task return
**"undecided"** — deliberately, because the alternative is falling back to weaker evidence:

| Component | Why it cannot run | Consequence |
|---|---|---|
| `yersinia_pestis_vs_y_pseudotuberculosis` → `cgmlst_species_assignment` | BIGSdb-Pasteur cgMLST has no packaged CLI and no pinned allele caller | The Yersinia task returns "undecided". Classical 7-locus MLST does run, but its ST **cannot** separate the two — *pestis* is a clone nested inside *pseudotuberculosis* — so it is context only. VFDB + PlasmidFinder contain none of `pst`, `ypo2088`, `yihN` or `opgG`, so the abricate pre-screen adds no discriminating power either. |
| `listeria_monocytogenes_serogroup_lineage_cc` → `cgmlst_lin` | Same: no packaged CLI, no pinned allele caller for scheme 15 | LisSero's serogroup and classical MLST both run, but the broad lineage does not, so the four Listeria lineage units are checked by the ANI panel instead. Classical MLST reports **no** clonal complex and no lineage, and no ST→CC→lineage table is pinned, so the lineage field stays empty rather than being inferred. |
| `clostridium_botulinum_group_and_bont_type` → `genome_group` | No genome-wide group I–IV classifier is pinned | The group verdict for the two C. botulinum group units comes from the ANI panel, **never** from the toxin type: group I and group II strains can carry overlapping bont types. |

`docs/CURATION_TRIAGE_TIER2.md` records what each installed tool actually reports, tool by tool.

**DB-pinning by tool class:** (a) *bundled-in-package* DBs — SeqSero2, ShigEiFinder, ShigaTyper, ECTyper
alleles, LisSero, `mlst` — pinning the conda version pins the DB; (b) *Zenodo-fetched* — SISTR DB,
ECTyper's MASH species sketch — pin the Zenodo DOI/version or pre-stage in CI; (c) *AMRFinderPlus* — the
only tool with a first-class dated DB version (`amrfinder --database_version`), the best reproducibility
anchor; (d) *abricate* — bundles static snapshots with no semver, so freeze the whole `$ABRICATE_DB` dir
+ checksum and record `abricate --list` (per-DB sequence count + date); (e) *ShigaPass* — not packaged,
pin a git commit if used. `e→d` means a scheme database with no version number of its own, so it is
pinned the way abricate is: freeze the whole snapshot and checksum it.

### Provenance: version the interpretation rule separately

A marker call is the product of three things that change independently: the **software**, the
**database**, and the **rule we apply to the output**. Recording only the first two makes a changed
call unattributable. So Tier 2 records all three, with `interpretation_rule_version` as its own field:

```yaml
method:   { name, software_version, conda_package, conda_build, container_digest }
database: { name, version_or_scheme_id, retrieved_utc, manifest_sha256 }
analysis: { command_line, thresholds, interpretation_rule_version, raw_output_paths }
```

Stable task identifiers, used as keys in the marker manifest and the ledger:
`bacillus_anthracis_vs_b_cereus_group`, `yersinia_pestis_vs_y_pseudotuberculosis`,
`listeria_monocytogenes_serogroup_lineage_cc`, `clostridium_botulinum_group_and_bont_type`.

**Update gate.** A new tool or database pin does not enter production until: (1) the candidate artifact
has an immutable digest; (2) a fixed validation panel has been re-run; (3) every call difference is
classified as a software, database, threshold, or interpretation-rule effect; and (4) the versioned
decision record and changelog are updated.

`bin/curation-triage/update_gate.py` implements all four as real checks over 38 pinned components
(11 software, 9 database, 8 rule, 9 threshold, and `panel:scope`) recorded in
`src/curation-triage/tool_pins.tsv` and `src/curation-triage/pin_changelog.tsv`. Three properties make
the gate mean something:

- **The validation panel is not a fixture.** It is the committed Tier-2 evidence replayed offline
  through the same functions the monthly run uses.
- **Point 2 also checks that the caches' own recorded tool versions are the pinned ones.** A panel
  replayed from the previous pin's evidence signs off on nothing.
- **The panel's scope is frozen in the gate (41 units) and pinned as `panel:scope`.** Read from the
  ledger instead, ordinary curation would silently move the panel's scope underneath the gate.

A call that changed while nothing it depends on moved is `unattributable`, which is a **failure**, not a
category — it means the same pins produced two different answers. See
`docs/CURATION_TRIAGE_REPRODUCIBILITY.md` §3.

## 6. Deliberately out of scope

Naming what the tool does *not* do is part of the design, because each exclusion protects the precision
target in §1.

| Not done | Why |
|---|---|
| **Pass/fail grading of a genome** | Kalamari's picks were never chosen by centrality, so a score would fail genomes that are correctly curated. The tool ranks for human review. |
| **Auto-PRs against `src/chromosomes.tsv`** | A curation change is an SME decision. The tool writes a worklist and a ledger; a human makes the edit. |
| **Outward scanning of all of NCBI for species we do not track** | Precision collapses. J2 is fed by SME intent (`chromosomes-todo.tsv`) and by newly validly-published names in genera we already track. |
| **Eukaryote pathogens** | Tier 0 marks them `eukaryote` and stops. NCBI's type-strain ANI report is prokaryote-only, so the primary signal does not exist for them. |
| **Viruses and organelles** | Excluded permanently by the Tier-0 eligibility guard — the comparison unit does not apply. |
| **A learned embedding as a resolver** | §3, principle 3. Embeddings can only ever be a novelty radar (§7). |
| **Recomputing ANI for every pick every month** | §3, principle 2. Tier 2 runs only on Tier-1 survivors and on the metadata-blind sub-species lineages. |

## 7. Research track — embedding novelty radar

Not in v1, and it never blocks v1.

- **Sketch-space novelty first:** embed refs + new deposits as Mash sketches (already built in CI),
  nearest-neighbor distance flags genomes far from everything; UMAP for SME visualization (keep the
  high-dim NN distance as the stable metric — UMAP is unstable on refit).
- **Optional Bacformer probe:** does a gene-content embedding separate the Regime-A confounded pairs
  better than ANI? Hypothesis test on a labeled 6-group panel (silhouette/ARI). Research only; deploy
  markers regardless.
- **Rejected for this use:** DNABERT-S/Genos-m (nucleotide axis = redundant with sketches), Evo2 (no
  native genome vector, heavy), vPST (viral-only), NT-v2 & ProkBERT (non-commercial licenses).

## 8. CI, reproducibility, and the ledger

`docs/CURATION_TRIAGE_REPRODUCIBILITY.md` is the operational account.

- Non-gating monthly workflow (`curation-triage.yml`): `schedule` + `workflow_dispatch` + `v*.*.*` tags,
  **out of required status checks** (that is what makes it non-gating — *not* blanket
  `continue-on-error`; real infra failures fail the shard loudly-but-non-blocking; collate prints
  `shards_failed=N/total`). It reuses `unit-testing.yml`'s sharded-matrix pattern:
  `generate-matrix → tier1-shard ×8 → collate → tier2-replay`, then **one** ledger commit carrying
  both tiers via fetch-rebase-push.
- **The monthly run installs no tools and touches no network.** It replays committed caches, so it
  cannot produce a new number. Refreshing the caches is a separate `workflow_dispatch` workflow
  (`curation-triage-refresh.yml`) that installs one pinned environment per runner and **shards by
  tool**, because the eleven environments total ~13 GB and the sub-species panels are cliques that a
  per-unit split would recompute several times. Genomes are cached by accession, keyed on the hash of
  the plan's accession list, and warmed by a single downloader job.
- **Pin everything**: 11 committed conda lockfiles + `tool_pins.tsv` (build, binary, version command,
  lockfile digest, database kind) + `setup_envs.sh`, which recreates each environment and then proves
  it by running every binary and comparing what it reports. Versions and resolved params are stamped
  into every ledger record. An NCBI file-schema change or a recompute wave is an **epoch bump**, so
  correlated re-flags read as tooling, not biology.
- A `v*.*.*` tag attaches the worklist as a **Release asset** (durable; run artifacts expire at 90
  days), with `files` only so it never rewrites the notes upstream's `build-sketch` job authors.
  The worklist is **rendered from the committed ledger** (`worklist_from_ledger.py`), not read off
  disk (the file is gitignored: it is a view of the ledger, and a committed view drifts) and not
  recomputed (re-running would diff against a ledger that already holds this run's records, so
  every lead would read as standing). The ledger stores each run's delta, so the release shows the
  labels the curator saw — see `docs/CURATION_TRIAGE_REPRODUCIBILITY.md` §4.
- **A shared `concurrency` group cannot de-conflict the tag path with upstream's `build-sketch` job.**
  `concurrency` serialises a workflow against itself, upstream's workflow declares none, and making
  them serialise would mean editing an upstream file. The release step tolerates the race instead — it
  waits for the release object and updates it.
- **Keepalive**: a weekly job re-enables the workflow through the API with least-privilege
  `permissions: actions: write`, defeating GitHub's 60-day scheduled-workflow auto-disable. It
  prevents that state; it cannot escape it, since a disabled workflow's crons stop firing too.
- The **update gate** (§5, four points) runs as its own CI job, offline and toolless:
  `update_gate.py --check`.

## 9. Testing

Unit tests against the real `src/chromosomes.tsv`: synthetic taxids resolve, the two Y. enterocolitica
strains stay separate, Arripis is excluded by the guard, multi-replicon units concatenate, preflight
fails on an unresolvable taxon. Plus a golden test for the metric/status logic and a PR smoke test on a
tiny pilot subset (enteric core: Salmonella, Listeria, E. coli, Campylobacter, Vibrio, Neisseria),
analogous to the existing `KALAMARI_DEBUG` path.

## 10. File layout

```
bin/curation-triage/             # 24 modules + setup_envs.sh + tests/
    build_policy.py policy.py taxonomy.py ncbi.py          # Tier 0
    coverage_probe.py coverage.py ani_reports.py           # Tier 1 J1a
    tier1_metadata.py leads.py ledger.py collate.py        # Tier 1 net
    datasets_summary.py lpsn.py                            # J1b / J2 feeders
    tier2_confirm.py tier2.py markers.py                   # Tier 2 driver + pure logic
    genomes.py skani_runner.py typers.py                   # Tier 2 cached/injectable layers
    build_panels.py build_type_strains.py                  # regenerate the committed Tier-2 inputs
    refresh_tier2.py setup_envs.sh update_gate.py          # refresh, pin, gate
    worklist_from_ledger.py                                # the release asset's worklist
src/curation-triage/             # committed policy table, caches, ledger, reviews,
                                 #   marker_manifest.tsv, panels/*.acclist (17),
                                 #   envs/*.lock (11), tool_pins.tsv, pin_changelog.tsv,
                                 #   datasets_summary_cache.json, lpsn_cache.json,
                                 #   amrfinderplus_db_manifest.tsv
.github/workflows/curation-triage.yml          # monthly triage + release asset + keepalive
.github/workflows/curation-triage-refresh.yml  # workflow_dispatch cache refresh, sharded by tool
.github/workflows/curation-triage-tests.yml    # pytest + update gate on push/PR (gating)
pyproject.toml                                 # pytest rootdir + testpaths
docs/CURATION_TRIAGE_TIER0.md
docs/CURATION_TRIAGE_TIER1.md
docs/CURATION_TRIAGE_TIER1_COVERAGE.md
docs/CURATION_TRIAGE_TIER2.md
docs/CURATION_TRIAGE_REPRODUCIBILITY.md
```

There is no single `run.sh` entry point: each tier is a self-contained CLI. Tier 2 is five modules, not
one, so the pure logic (`tier2.py`, `markers.py`) stays apart from the three external layers it needs
(`genomes.py`, `skani_runner.py`, `typers.py`) — the same separation Tier 0/1 keep between `policy.py`
and `ncbi.py`. Only the pure modules are exercised by the offline test suite; the external layers are
injectable so the tests never need a network or a conda environment.

Config defaults (all overridable): `species_ani=95`, `af_gate=0.5`, `derep=99.5`, `N_floor=5`,
`N_target=12`, `epsilon≈0.2 ANI`, `size_cap=50 Mbp`, marker `%id/%cov` per typer defaults.

## 11. Status

Everything below is committed and reproduces offline from the committed caches. Every number comes from
a committed artifact.

| Component | State | What is committed |
|---|---|---|
| **Tier 0 policy table** | built | `taxon_policy.tsv`: 270 units — 249 eligible (21 of them confounded: 4 regime A, 17 regime B), 14 organelle, 4 virus, 3 eukaryote. Every eligible unit resolves to a real GCF/GCA parent assembly; preflight passes. |
| **Tier 1 J1a coverage** | built, measured | **249/249** eligible picks have a live row in NCBI's `ANI_report_prokaryotes.txt` (246 the exact version, 3 a different version of the same assembly), and 249/249 in `assembly_summary`. `coverage_report.tsv`. |
| **Tier 1 net + SME handoff** | built | `ledger.ndjson` with 21 J1a records (run `2026-07`), the `reviews.tsv` decision loop with "suppress until the evidence changes" semantics, and the sharded collate path — all offline. |
| **J1b / J2 feeders** | built | `datasets_summary_cache.json` (170 taxid queries, 3208 records, built with the pinned `ncbi-datasets-cli 18.34.0`) and `lpsn_cache.json` (28794 names, LPSN release `2026-08-04`). Both replay offline; the monthly job runs `--with-lpsn`. |
| **Tier 2 confirmation** | built and run | Trigger set **41 units** (21 Tier-1 survivors + 21 confounded, overlapping in one). A full offline run evaluates **41/41** → 33 findings, 42 confirmations, 4 regime-A marker overrides, 75 ledger records (run `2026-08`). Caches: 157 skani rows over 158 planned comparisons (one pair is served by its reverse key), 38 typer results over 37 planned jobs (ECTyper adds a context-only serotype), and a type strain for all 249 picks. |
| **Marker manifest** | 12 of 15 components pinned | `marker_manifest.tsv`: 15 components over 8 tasks, 17 `panels/*.acclist`. Three components are `gap` (§5). |
| **Reproducibility + update gate** | built, verified by running | 11 conda environments pinned by explicit lockfiles (1148 package URLs), `tool_pins.tsv`, and a `setup_envs.sh` that recreates *and verifies* each one. `update_gate.py --check` passes over 38 pinned components and a 75-call panel across 41 units, with no `CONDA_PREFIX`, a bare `PATH`, and both proxies pointed at a dead port. |
| **Test suite** | 531 checks pass | `python3 -m pytest` from the repo root, ~25 s, fully offline. |

**Two consequences of the coverage number.** Tier 1 carries the drift/mislabel half of J1 on its own, so
no self-run type-strain ANI is needed for coverage; Tier 2's own ANI is reserved for the metadata-blind
sub-species lineages and for confirmation on Tier-1 survivors. And the coverage holds because the join
keys on the real GCF/GCA parent assembly rather than on INSDC replicon accessions — NCBI's verdict then
covers every pick regardless of the pick's age.

**The 21 J1a leads** over 249 picks break down as 17 taxonomy-check Inconclusive, 2 Failed, 3 RefSeq
identity exclusions and 1 reclassified (a pick can carry more than one flag), with the 21 confounded
units held out of J1b.

### Facts the data forces on the design

These are properties of NCBI's data and of the tools, not of any one run:

- **The type strain is keyed per assembly, not per species.** NCBI names one per assembly and it differs
  between subspecies: the eleven Salmonella units all share species taxid 28901 but resolve to seven
  different type strains.
- **"No type material" and "this pick IS the type strain" are answers, not gaps.** 30 picks are their own
  species' type strain and 7 species have no type material at all — which is also why NCBI's taxonomy
  check on several flagged picks reads Inconclusive.
- **skani returns no row for a distant pair**, which reads as "not measured" when it means the opposite.
  The screen is widened and a screened-out pair is recorded as a result.
- **A `gap` component returns "undecided"** rather than falling back to weaker evidence (§5).
- **LPSN publishes no date column.** The valid-publication date is derived from the authority string (cut
  at `emend.`, then the last year in what remains), and the `--max-lpsn` cap — not a date cutoff — is what
  keeps the worklist short.
- **A digest proves sameness, not availability.** The AMRFinderPlus database `2026-05-15.1` is frozen:
  all 162 files are digested into `amrfinderplus_db_manifest.tsv` with roll-up `2fb16b7f3c57…`, checked
  file by file by `setup_envs.sh` and folded into the update gate's identity. But `amrfinder -u` only
  ever fetches the latest database, so an archived copy of that dated directory is the only way back to
  it. `docs/CURATION_TRIAGE_REPRODUCIBILITY.md` §5 covers both of these residues.

### Two CI jobs, two questions

`.github/workflows/curation-triage-tests.yml` runs pytest and the update gate as separate jobs on push
and pull request. They ask different questions: **pytest asks whether the code behaves; the gate asks
whether the pinned world moved and whether the move is recorded.**

```bash
python3 -m pytest                            # 531 checks; rootdir + testpaths from pyproject.toml
python3 bin/curation-triage/update_gate.py --check
```

## 12. Open items

1. **No workflow has run on a real GitHub runner.** Every workflow is lint-clean and every step script has
   been executed locally by extracting it from the parsed YAML, but GitHub Actions itself has never
   executed one. The job to watch is the AMRFinderPlus environment: its `amrfinder -u` branch is
   unexercised, because running it locally would destroy the pinned `2026-05-15.1` database.
2. **The AMRFinderPlus archive has no durable home.** The database is frozen and verifiable, but the only
   copy of the bytes is a single tarball on a single machine, and NCBI cannot serve that dated version
   again. Two copies in institutional storage with their digests checked on a schedule is the whole task.
   It is the only open item whose cost, if skipped, is unrecoverable.
3. **Three `gap` marker components have no pinned caller** (§5): an allele caller for the BIGSdb
   *Yersinia* and *Listeria* cgMLST schemes, and a genome-wide *C. botulinum* group I–IV classifier. Each
   would be a new pin, so it enters through the update gate and its moved calls get changelog rows.
4. **The embedding novelty radar is research, not code** (§7). Sketch-space (Mash) novelty is the primary
   baseline; Bacformer is the single learned pilot if it is pursued. Never a resolver.
5. **LPSN needs an account.** The committed cache was downloaded by a registered curator. The LPSN REST
   API is the unattended replacement, and it is not wired up.

Four upstream data fixes are applied to `src/chromosomes.tsv` and captured for upstream submission in
`chromosomes-data-fixes.patch`: one `source` value misspelled `NBCI-REF` (hence 179 `NCBI-REF` rows, not
178), and three taxid mislabels caught by the Tier-0 preflight.
