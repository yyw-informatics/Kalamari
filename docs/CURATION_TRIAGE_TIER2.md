# Curation triage — Tier 2: targeted confirmation

Tier 1 reads NCBI's own verdict on all 249 eligible picks and hands the curator a short
worklist. **Tier 2 answers the questions metadata cannot**, at the sequence level, and only
for the units that need it.

It is still informational. It never gates anything, it opens no pull requests, and a
Tier‑2 verdict that agrees with the current label *annotates* the Tier‑1 row rather than
removing it — the curator still decides.

---

## 1. What Tier 2 runs on

Two trigger groups, 41 units in total:

| Trigger | Units | Why |
|---|---|---|
| Tier‑1 survivors (the leads on the worklist) | 21 | attach a number to NCBI's metadata verdict |
| Confounded units (mandatory, flagged or not) | 21 | metadata cannot answer these at all |
| overlap (`Bacillus_cereus~GCF_000283675.1` is both) | −1 | |

Tier 2 is **confirmation and sub‑species residue, never coverage**. NCBI's ANI report covers
100% of the eligible picks (249/249), so Tier 1 already sees every pick; nothing here exists
to find leads Tier 1 missed.

## 2. The three checks

| Check | Runs on | Question |
|---|---|---|
| `ani` | Tier‑1 survivors, and regime B | how far is our pick from the type strain of its declared species? |
| `panel` | regime B (17 units) | is this lineage reference still distinct from its siblings, and where does it sit? |
| `marker` | all confounded units | what do the defining genes say? |

### The aligned fraction is a veto, not a tiebreaker

An ANI computed over a sliver of a genome is not weak evidence, it is no evidence. So the
gate checks the aligned fraction **first**, and takes the *smaller* of the two sides: a small
genome fully contained in a larger one aligns well on its own side while covering little of
the other, which is exactly the case the gate exists to catch. Below `af_gate` the verdict is
"cannot decide", never "different species".

Within `epsilon` of the species boundary the answer is "too close to call". A precision‑first
tool must be allowed to say that.

### The regime‑A override is structural

For *Shigella* vs *E. coli*, *B. anthracis* vs the *B. cereus* group, and *Y. pestis* vs
*Y. pseudotuberculosis*, genome‑wide similarity is the wrong ruler. So a regime‑A unit gets
**no independent ANI verdict at all**: the ANI number is folded into the marker lead as
`ani_context` and the marker call is the verdict. There is no code path in which an ANI number
can contradict a regime‑A marker call, because only one of them is ever a verdict.

If the markers cannot decide, the answer is "undecided" — the decision does *not* fall back to
ANI, because for these taxa the ANI answer was never admissible.

For regime B the panel decides membership and the typer supplies the label. A wrong label on a
correctly‑placed genome is still a defect, so a refuting label is reported.

## 3. The marker manifest

`src/curation-triage/marker_manifest.tsv`, one row per (task, component). A marker call is the
product of three things that move independently, so all three are pinned and stamped into every
ledger record:

```
software  x  database  x  the rule we apply to the output
```

The third is ours. It lives in `markers.py` as code and carries its own
`interpretation_rule_version`. Recording only the tool and the database makes a changed call
unattributable — you cannot tell a new BTyper3 build from a new reading of the same output.

Every build below is **installed and has been run**. The manifest records what ran, not what a
package's documentation advertises; where a documented default and the installed build differ,
the manifest names the installed one.

| Task | Classifier (installed build) | Database pin |
|---|---|---|
| `shigella_eiec_vs_escherichia_coli` | ShigEiFinder `1.3.4-pyhdfd78af_0` | bundled → conda build |
| `bacillus_anthracis_vs_b_cereus_group` | BTyper3 `3.4.0-pyhdfd78af_0` **+ our rule** | bundled + PubMLST copy 2023‑07‑07 |
| `yersinia_pestis_vs_y_pseudotuberculosis` | **gap** — no runnable cgMLST caller; `mlst 2.23.0` is context only | frozen scheme snapshot (not installed) |
| `listeria_monocytogenes_serogroup_lineage_cc` | LisSero `0.4.10-pyhdfd78af_0` + `mlst 2.23.0`; cgMLST/LIN is a **gap** | bundled 2020‑06‑04 + scheme snapshot |
| `salmonella_subspecies_and_serovar` | SISTR `1.1.3-pyhdc42f0e_2` + SeqSero2 `1.3.2-pyhdfd78af_0` | SISTR_V_1.1.3_db (Zenodo, 2025‑04‑24) |
| `clostridium_botulinum_group_and_bont_type` | AMRFinderPlus `4.2.7-hf69ffd2_0`, DB `2026-05-15.1` | dated DB (the best anchor) |

**What the rules refuse to do**, each because it would produce a confident wrong answer:

- A plasmid marker never makes a *Yersinia* species call on its own — plasmid‑cured *Y. pestis*
  and plasmid‑bearing *Y. pseudotuberculosis* both exist. (This is why abricate was dropped as
  the classifier; it survives only as a pre‑screen.)
- A *Listeria* serogroup is never reported as a lineage. Serogroup, MLST ST/CC and cgMLST/LIN
  answer different questions and stay separate fields.
- The *C. botulinum* genomic group is never inferred from the toxin type. No genome‑wide
  group I–IV classifier is pinned yet, and the manifest records that row as a `gap` rather
  than hiding it.
- Plasmid‑cured *B. anthracis* is still *B. anthracis*, and a *B. cereus* carrying anthrax
  plasmids is never labelled plain *B. cereus*.
- "anthracis" is matched as a species epithet, never as a substring: *Bacillus paranthracis*
  contains the string and is a different organism — and it is what NCBI's best ANI match
  reports for one of Kalamari's own picks.

## 4. The sub-species panels

`src/curation-triage/panels/<kalamari_taxid>.acclist`, generated by `build_panels.py` from the
policy table (a test asserts the committed files match, so a policy edit cannot silently drift).

The 13 units with synthetic taxids are **metadata‑blind**: NCBI knows nothing about a lineage
only Kalamari defines, so Tier 1's whole signal is unavailable for them. Their panel members are
**Kalamari's own sibling units** — real accessions already curated in `chromosomes.tsv`. That
answers two questions offline, without inventing reference accessions we cannot verify:

- **Is the panel still a panel?** A sibling at or above `derep` (99.5%) means two of Kalamari's
  own lineage references are effectively the same genome. That is a curation defect nothing else
  in the tool can see.
- **Where does a genome sit?** The nearest member and the margin to the runner‑up.

With a sibling panel the unit's own lineage is *not* in the panel, so "nearest" is a
neighbourhood, not a lineage assignment — the code says so rather than over‑claiming. Each file
keeps a commented `external:` block for curated outside references, to be filled only after a
network pass has verified each accession.

## 5. How a Tier-2 result reaches the curator

Tier 2 is a **new signal on the same ledger**, not a new ledger. Leads carry `signal="T2"` and a
`t2:<unit_id>#<check>` key, and flow through the same delta, `reviews.tsv` suppression and
worklist render as Tier‑1 leads.

- **Delta** = *did the verdict change?* A flip from confirmed to refuted is the whole point of
  re‑running, so it surfaces as new work. An unchanged verdict is standing.
- **The fingerprint holds the verdict, not the ANI number.** A skani bump moves the last decimal
  place of every comparison; busting every dismissal over 0.01% ANI would train curators to
  ignore the tool. A genuine tooling change is handled by the **epoch** instead — and a database
  refresh *is* an epoch, so a changed call is re‑reviewed as a database effect rather than read
  as new biology.
- **Confirmations are recorded but are not work.** They render as annotations under the finding
  list, cross‑linked to the Tier‑1 row via `confirms_lead`.
- **Tier 2 marks nothing resolved.** Its leads say nothing about whether a Tier‑1 lead is still
  actionable, so `partition(resolve_signals=())` stops it from touching the Tier‑1 baseline.
- **Context rides on a lead but can never move it.** ECTyper's cached fields sit on the
  *E. coli* marker lead's evidence under `ectyper_*` keys — serotype `O118/O151:H16`, QC
  `WARNING MIXED O-TYPE`, database `v1.0`. Three properties make that safe: they are written
  *before* the verdict fields, so context cannot overwrite a call; they are absent from
  `T2_FINGERPRINT_KEYS`, so a changed serotype cannot re‑open a lead a curator dismissed; and
  ECTyper is not routed as its own task, so it is never a verdict. The distinct key prefix is
  load‑bearing — ShigEiFinder and ECTyper both emit a plain `serotype`, and a merge that keeps
  the bare key drops ECTyper's value into ShigEiFinder's slot.

## 6. Running it

```bash
# offline, replaying the committed caches (no network, no tools needed)
python3 bin/curation-triage/tier2_confirm.py --offline --run-id 2026-08

# refresh: download genomes, run skani + the typers, then evaluate
python3 bin/curation-triage/tier2_confirm.py --run-id 2026-08 \
  --typer-env-dir <dir of per-tool conda envs>

# skani is cheap and the typers are not, so each can be refreshed on its own
python3 bin/curation-triage/tier2_confirm.py --run-id 2026-08 --skip-typers

# narrow the scope
python3 bin/curation-triage/tier2_confirm.py --offline --only-confounded
```

Thresholds are all overridable: `--species-ani 95 --af-gate 0.5 --derep 99.5 --epsilon 0.2`.

Units with no cached evidence produce **nothing**. "Not measured" is not a finding, so it is
reported as a count, not as 41 indeterminate worklist rows.

An `--offline` run whose caches do not cover the plan **fails and names what is missing**
(with `--allow-uncached` for a deliberate partial run). Exiting 0 there would be a silent
no-op that reads as a clean run: zero leads built, zero records appended, and nothing said
about why.

A refresh (no `--offline`) preflights the typers that have uncached work **before** it downloads
a genome or runs skani, so a wrong `--typer-env-dir` stops the run instead of turning every job
into a miss. `--allow-unpinned` downgrades that to a warning. `--offline` never reaches it.

To refresh the caches in shards instead of one process holding every tool at once, see
[`refresh_tier2.py`](../bin/curation-triage/refresh_tier2.py) and
[Reproducibility](CURATION_TRIAGE_REPRODUCIBILITY.md) §2. When those shards are merged back, the cache being
replaced is passed as `--base`: it supplies only the keys no fragment answers, so a shard can
still move a call, but a result no shard produces — ECTyper's, since it is context-only and has
no shard — is not deleted by a full refresh. A merge that would drop a key is refused by name.

## 7. Files

| Path | What it is |
|---|---|
| `bin/curation-triage/tier2.py` | pure: the gate, panel placement, the override, the leads |
| `bin/curation-triage/markers.py` | pure: the manifest reader and our versioned rules |
| `bin/curation-triage/genomes.py` | cached FASTA + the type-strain map (injectable) |
| `bin/curation-triage/skani_runner.py` | cached skani (injectable) |
| `bin/curation-triage/typers.py` | cached typer adapters (injectable) |
| `bin/curation-triage/tier2_confirm.py` | the driver |
| `bin/curation-triage/refresh_tier2.py` | sharded cache refresh: accession list, per-tool shards, strict merge |
| `bin/curation-triage/build_panels.py` | regenerates the panels from the policy table |
| `src/curation-triage/marker_manifest.tsv` | the committed pins |
| `src/curation-triage/tool_pins.tsv` | the pinned environment behind each of those pins |
| `src/curation-triage/panels/*.acclist` | 17 committed panels |
| `bin/curation-triage/build_type_strains.py` | resolves every pick's type strain from NCBI's ANI report |
| `src/curation-triage/tier2_ani_cache.json` | 157 real skani rows, answering all 158 planned comparisons (one pair is served by its reverse key); 4 of them record a pair skani screened out |
| `src/curation-triage/tier2_typer_cache.json` | 38 real normalised typer results: the 37 planned jobs plus ECTyper's context-only serotype, which no shard refreshes |
| `src/curation-triage/type_strain_cache.json` | 249/249 picks → their declared type strain |
| `bin/curation-triage/tests/fixtures/tier2/raw/` | verbatim output from each of the 9 tools |

Genome FASTA is **never committed** — it lives in the gitignored `.cache/curation-triage/genomes/`.
Only the small derived results are.

## 8. What has actually been run

Every tool named above is **installed and has been run for real** (2026-08-04), and the
three committed caches hold real results, not placeholders:

| Cache | Contents |
|---|---|
| `type_strain_cache.json` | 249/249 picks resolved from NCBI's ANI report |
| `tier2_ani_cache.json` | 157 skani 0.3.2 rows covering all 158 planned comparisons |
| `tier2_typer_cache.json` | 38 normalised typer results across 10 tools (37 planned + ECTyper as context) |

A full offline run evaluates **41 of 41 triggered units, with nothing left unmeasured**:
33 findings to review, 42 confirmations, 4 regime‑A marker overrides.

### Tool behaviour a naive integration gets wrong

Each of these produces a confident wrong answer if you trust a tool's documentation
instead of its output:

- **AMRFinderPlus** names the column `Element symbol`, not `Gene symbol`, and writes
  `bont_E3` with an underscore.
- **BTyper3 3.4** uses the revised nomenclature (`species(ANI)` = `mosaicus`, the
  familiar name in `subspecies(ANI)`) and writes to a file, not stdout. Four of the
  nine tools write files, each with its own path convention.
- **skani** drops pairs below ~80% identity or 15% aligned fraction and returns **no
  row**, which reads as "not measured" when it means "very distant". So the screen is
  widened and a screened‑out pair is recorded as a result in its own right.
- **`mlst`** prints no header at all; **LisSero** reports `FULL`/`PARTIAL`/`NONE`;
  **ShigEiFinder** has no organism column and its `ipaH` is a bare `+`/`-`.
- Two conda environments need `setuptools<81` added (modern setuptools removed
  `pkg_resources`), LisSero needs its own `bin` on PATH for `makeblastdb`, `mlst` needs
  it for its Perl libs, and ShigEiFinder needs its `--tmpdir` to exist already —
  otherwise it dies *after* writing a header, leaving a file that looks like a clean
  "no results".

Two properties of the real data drive the design:

- **The type strain is keyed per assembly, not per species.** NCBI names a declared
  type strain per assembly, and it differs between the subspecies of one species. All
  eleven Kalamari *Salmonella* units share species taxid 28901 but get eleven different,
  subspecies-appropriate type strains; keying by species compares ten of them against
  the wrong genome.
- **Two states are answers, not missing values.** 30 picks *are* their species' type
  strain (NCBI writes `same`), and 7 species have **no type material at all** — which is
  also why NCBI's own taxonomy check on several Tier‑1 flagged picks is Inconclusive.

### Two rows the manifest records as gaps

These two schemes look pinnable on paper, but nothing here can run them, so the
manifest records them as `gap` and the affected tasks return "undecided" rather than
falling back to weaker evidence:

- **Yersinia cgMLST scheme 1** — no packaged CLI, no pinned allele caller. Classical
  7‑locus MLST runs, but on real genomes its ST cannot separate *Y. pestis* from
  *Y. pseudotuberculosis*, so it is context only. Worse, VFDB and PlasmidFinder contain
  **none** of `pst`, `ypo2088`, `yihN` or `opgG`, so the abricate pre-screen adds no
  discriminating power either.
- **Listeria cgMLST1748_v2 / LIN** — same reason, so the four lineage units get no
  lineage verdict from markers; their check comes from the ANI panel instead.

### Findings worth an SME's eye

- ***B. thuringiensis*** (`GCF_000008505.1`): its nearest type strain by ANI is
  ***B. anthracis*** (97.84%), but the markers rule anthrax out — no anthracis
  chromosome signature, no pXO1, no pXO2. This is exactly the case regime A exists for,
  and the override holds on real data. BTyper3 assigns it no biovar because the strain
  (Al Hakam) carries no Bt toxin genes, so the label is reported as unconfirmed rather
  than wrong.
- ***B. cereus*** (`GCF_000283675.1`): 91.77% ANI to the *B. cereus* type strain, below
  the species boundary, with *B. paranthracis* the nearest type strain — which is what
  NCBI's taxonomy check has been unhappy about.
- Five more picks sit below the species boundary against their own type strain, and
  three are too distant for skani to report at all.

## 9. Open items

The environments behind every call above are pinned, and the monthly workflow replays both tiers
offline — see [Reproducibility](CURATION_TRIAGE_REPRODUCIBILITY.md).

What is still open is the **two `gap` marker rows**: the BIGSdb *Yersinia* and *Listeria* cgMLST
allele callers, a genome-wide *C. botulinum* group I–IV classifier, and the *Listeria* CC →
lineage table. This is a missing tool, not a missing pin. Closing it means pinning a **new**
classifier, which changes the marker manifest, moves calls, and has to pass the update gate on
its way in — so it is kept apart from any change whose job is "change nothing, prove it repeats",
where a moved call would be ambiguous. Until then the four *Listeria* lineage units and the two
*C. botulinum* group units return "undecided", and the ANI panel is their check.
