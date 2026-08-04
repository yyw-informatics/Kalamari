# Tier 0 — the policy table: which genomes get checked, and which are excluded

The Kalamari curation-triage add-on is an automated, monthly, *informational*
worklist that flags newly deposited public genomes for subject-matter-expert
(SME) review. Tier 0 is the foundation every later tier keys off — a committed,
human-reviewable **policy table** that says, for each Kalamari genome, what it
really is and how it should be handled.

- **Code:** [`bin/curation-triage/`](../bin/curation-triage/)
- **Output:** [`src/curation-triage/taxon_policy.tsv`](../src/curation-triage/taxon_policy.tsv)
- **Design:** [`CURATION_TRIAGE_DESIGN.md`](CURATION_TRIAGE_DESIGN.md)

## What it does

`build_policy.py` reads [`src/chromosomes.tsv`](../src/chromosomes.tsv) plus NCBI
taxonomy and emits one row per **genome-unit** — 270 units from the 323 replicon
rows, of which 249 are `eligible`:

1. **Groups replicons into genome-units by assembly — never by the Kalamari
   `taxid` column.** The `taxid` column is rank-heterogeneous and untrustworthy
   for grouping: taxid `630` covers two *different* Yersinia enterocolitica
   strains, while the two chromosomes of one Aliivibrio genome would otherwise
   look like two organisms. The real GCF/GCA assembly (resolved from NCBI) is the
   ground truth — it merges the replicons of one genome and keeps distinct
   genomes apart even when they resolve to one species (the eleven *Salmonella
   enterica* lineage entries all resolve to species taxid `28901` yet are eleven
   different assemblies).

2. **Resolves a real NCBI species-rank `species_taxid`.** The 13 synthetic
   Kalamari taxids (`9000000`–`9000017`, e.g. the four *Listeria monocytogenes*
   lineages) do not exist at NCBI; they climb through the bundled custom-taxonomy
   overlay to their real parent species. Strain/subspecies taxids climb to
   species; merged taxids are followed. A genus-level entry (`Ruminococcus_sp.`,
   taxid `1263`) with no species falls back to the genome's own organism taxid
   (which is a real species, *R. bovis*).

3. **Classifies `eligibility` with a programmatic guard, not a hand list:**
   | class | units | meaning |
   |---|---|---|
   | `eligible` | 249 | Bacteria/Archaea with a real species taxid — participates in the ANI pipeline |
   | `exclude:organelle` | 14 | an organelle-only submission (host mitogenome / plastid) or an oversized assembly (> ~50 Mbp). A guard catches what a hand list misses, such as the *Arripis trutta* fish mitogenome |
   | `exclude:virus` | 4 | Viruses/Viroids (including multi-segment genomes) |
   | `eukaryote` | 3 | eukaryote pathogens (*Cryptosporidium*, *Cryptococcus*) — recorded, but outside the ANI pipeline |
   | `review:genus-level` | 0 | Bacteria/Archaea that genuinely resolve no further than genus. No unit lands here on the committed table; the class exists so such a unit is surfaced for review instead of passing silently |

   The organelle guard fires only when **every** replicon of a unit is
   non-nuclear, so a eukaryote nuclear assembly that merely bundles its
   mitochondrion (*Cryptococcus*: 14 chromosomes + 1 mito) is correctly
   `eukaryote`, not `exclude:organelle`.

4. **Tags confounded taxa** — 21 units — from the committed
   [`confounded_manifest.tsv`](../src/curation-triage/confounded_manifest.tsv)
   with `regime` (A/B) and `resolver` (marker/ANI). These are the taxa where
   whole-genome similarity cannot split the species and a marker must; see the
   confounded-taxa section of
   [`CURATION_TRIAGE_DESIGN.md`](CURATION_TRIAGE_DESIGN.md).

5. **Preflight HARD-FAILS** (non-zero exit) if any `eligible` unit lacks a real
   species-or-below NCBI taxon — so a synthetic/strain/genus key can never
   silently pass.

## Output columns (`taxon_policy.tsv`)

Required columns, in this order:

| column | meaning |
|---|---|
| `unit_id` | stable, readable id: `<repName>~<groupKey>` |
| `scientificName` | Kalamari name(s) for the unit (`;`-joined if the unit spans rows) |
| `kalamari_taxid` | Kalamari-declared taxid(s) |
| `species_taxid` | resolved real NCBI species-rank taxid (blank if none) |
| `accession_set` | `,`-joined nuccore accessions in this genome-unit |
| `parent_assembly` | real GCF/GCA if resolved, else a `strain:<taxid>` proxy |
| `eligibility` | see the table above |
| `regime` | confounded regime `A`/`B` (blank if not confounded) |
| `resolver` | `marker` / `ANI` (blank if not confounded) |

Extra review columns follow: `species_name`, `superkingdom`, `total_length_bp`,
`n_replicons`, `flags`, `note`. Useful `flags` values: `shared_accession` (one
physical sequence labelled two ways — e.g. `CP015575`), `label_conflict`,
`grouped_by_proxy` (multi-replicon unit grouped without a real assembly — confirm
manually), `missing_record` (NCBI returned no record for a replicon),
`species_from_organism` (species taken from the genome's organism because the
declared taxid was genus-level), `taxid_organism_mismatch` (the Kalamari `taxid`
column climbs to a different species than the deposited genome's own NCBI
organism — a mislabel surfaced for SME review; the declared species is kept and
the organism is recorded in `note`), `confounded`.

When no real assembly is resolved, the fallback grouping key is scoped to the
Kalamari identity (`scientificName` + `taxid`), so distinct genomes that share a
species-level organism taxid can never be merged into a chimera by the fallback.

## How to run

Taxonomy comes from a full NCBI taxdump (`nodes.dmp`/`names.dmp`/`merged.dmp`),
with the bundled `src/taxonomy/build` overlay layered on top for the synthetic
Kalamari lineages. The bundled `deprecated` dump is a stale partial snapshot and
is only a last-resort fallback.

```bash
# First build: fetches NCBI facts (nuccore esummary + elink->assembly) and
# writes the committed cache src/curation-triage/nuccore_cache.json, then the
# policy table. Set EMAIL (and NCBI_API_KEY if you have one) to raise NCBI rate
# limits, as the Kalamari README describes.
EMAIL=you@example.org python3 bin/curation-triage/build_policy.py

# Offline rebuild from the committed cache (no network; fails if incomplete):
python3 bin/curation-triage/build_policy.py --offline

# Point at a specific taxdump:
python3 bin/curation-triage/build_policy.py --taxdump /path/to/taxdump
# or export NCBI_TAXDUMP=/path/to/taxdump
```

**Recording a taxdump location once per machine.** Instead of passing `--taxdump`
every time, write the directories into `.taxdump-search` at the repo root — one per
line, `#` comments allowed — and the first one holding `nodes.dmp` is used. That file
is **gitignored on purpose**: an absolute path in a committed file would publish one
institution's filesystem layout and would be wrong on every other machine, so the
repository ships no default search path at all. `$NCBI_TAXDUMP_SEARCH` does the same
job for CI and one-off shells, separated like `$PATH`. With none of them set the
resolver falls back to the bundled partial dump and says so.

Useful flags: `--no-assemblies` (skip elink; group on the strain-taxid proxy),
`--no-write` (compute + preflight + summary only), `--size-cap-bp N`.

## Design: offline logic vs. the network

Everything except one module is pure and offline, so the whole classification is
unit-testable without the network:

- [`taxonomy.py`](../bin/curation-triage/taxonomy.py) — taxdump parsing + lineage
  climb (offline).
- [`policy.py`](../bin/curation-triage/policy.py) — grouping, eligibility guard,
  confounded tagging, preflight, table writer (offline, pure).
- [`ncbi.py`](../bin/curation-triage/ncbi.py) — the **only** network code: a
  cached, mockable layer that fetches the three NCBI facts per accession
  (parent genome, molecule type, length) via `esummary` and the real assembly
  via `elink`. Results persist to the committed `nuccore_cache.json`; the fetch
  function is injectable so tests never touch the wire.

## Tests

```bash
cd bin/curation-triage && python3 -m pytest tests/ -q
```

The suite covers six cases that must not regress: synthetic taxids resolve to
their parent species; the two *Yersinia enterocolitica* strains stay separate;
*Aliivibrio* multi-replicon concatenates to one unit; the *Arripis* fish
mitogenome is caught by the guard; the `CP015575` shared-accession quirk is
flagged, not double-counted; and preflight fails on a deliberately unresolvable
taxon.

## Upstream data note

`src/chromosomes.tsv` carries one `source` value misspelled `NBCI-REF` upstream
(vs. 179 `NCBI-REF`). This fork corrects it, and
[`chromosomes-data-fixes.patch`](../chromosomes-data-fixes.patch) holds that fix
and three taxid corrections for submission upstream. Tier 0 never reads the
`source` column, and the lane table in
[`coverage.py`](../bin/curation-triage/coverage.py) still maps `NBCI-REF` to
`NCBI-REF`, so an unpatched checkout keeps working.
