# Curation triage — running it: stages, caches, pins, and the gate

This is the operator's document. It covers how to run each stage, how to refresh the
cached evidence, how the tool environments are pinned, and what a change to a pin has to
pass before it lands.

Running the tools once, on one machine, proves the calls are right *here, today*.
Reproducing them **somewhere else, later** is a different problem, and it is the one
that decides whether a verdict is evidence or an anecdote. Three properties make it
evidence:

1. Anyone can rebuild the exact environment that produced a call.
2. CI can re-run the whole loop with no tools and no network, and can refresh the caches
   on demand without holding 13 GB of conda environments at once.
3. When a pin moves, every call that moves with it is explained, in writing, before it
   lands.

- What the output means: [for curators](CURATION_TRIAGE_CURATOR.md)
- Why any of this is shaped the way it is: [design](CURATION_TRIAGE_DESIGN.md)

---

## 1. Two ways to run

The monthly job and the refresh do different jobs, and confusing them is the main way to
get a wrong number into the ledger.

| | Monthly replay | Manual refresh |
| --- | --- | --- |
| Workflow | `curation-triage.yml` | `curation-triage-refresh.yml` |
| Trigger | cron, 06:17 UTC on the 1st | `workflow_dispatch` only |
| Network | none | NCBI genome downloads, NCBI datasets API |
| Tools installed | none | one pinned environment per runner |
| Reads | the committed caches | the tools |
| Writes | `ledger.ndjson`, in one commit carrying both tiers | the caches |
| Can move a call | no | **yes** |

**The monthly run is a replay.** Because it runs no tools, it cannot produce a new
number. That is the point: a verdict changes only when the evidence changes, and the
evidence changes only in a refresh.

## 2. Running each stage

### Tier 0 — the policy table

Taxonomy comes from a full NCBI taxdump (`nodes.dmp`, `names.dmp`, `merged.dmp`), with
the bundled `src/taxonomy/build` overlay layered on top for the synthetic Kalamari
lineages. The bundled `deprecated` dump is a stale partial snapshot and is only a
last-resort fallback.

```bash
# First build: fetches NCBI facts (nuccore esummary + elink->assembly), writes the
# committed cache src/curation-triage/nuccore_cache.json, then the policy table.
# Set EMAIL (and NCBI_API_KEY if you have one) to raise NCBI's rate limits.
EMAIL=you@example.org python3 bin/curation-triage/build_policy.py

# Offline rebuild from the committed cache (no network; fails if it is incomplete):
python3 bin/curation-triage/build_policy.py --offline

# Point at a specific taxdump:
python3 bin/curation-triage/build_policy.py --taxdump /path/to/taxdump
```

Other useful flags: `--no-assemblies` (skip elink and group on the strain-taxid proxy),
`--no-write` (compute, preflight and summarise only), `--size-cap-bp N`.

**Recording a taxdump location once per machine.** Instead of passing `--taxdump` every
time, write the directories into `.taxdump-search` at the repository root — one per
line, `#` comments allowed — and the first one holding `nodes.dmp` is used. That file is
**gitignored on purpose**: an absolute path in a committed file would publish one
institution's filesystem layout and would be wrong on every other machine, so the
repository ships no default search path at all. `$NCBI_TAXDUMP_SEARCH` does the same job
for CI and one-off shells, separated like `$PATH`. With none of them set, the resolver
falls back to the bundled partial dump and says so.

### The coverage probe

This measures whether NCBI has a live verdict for each pick. Re-run it whenever the pick
set changes or NCBI changes the shape of its reports. The first run streams and distils
the three `ASSEMBLY_REPORTS` files, about 2.7 GB in total, keeping only the rows that
touch our picks.

```bash
python3 bin/curation-triage/coverage_probe.py            # streams NCBI, writes cache + report
python3 bin/curation-triage/coverage_probe.py --offline  # replay from the committed cache
```

### Tier 1 — the monthly net

```bash
# Offline monthly run, J1a only, from the committed cache:
python3 bin/curation-triage/tier1_metadata.py --offline --run-id 2026-08

# Refresh the ANI cache from NCBI, then run (first run, or a cache rebuild):
python3 bin/curation-triage/tier1_metadata.py --run-id 2026-08

# With the other two feeders enabled:
python3 bin/curation-triage/tier1_metadata.py --offline --run-id 2026-08 \
    --with-j1b --with-j2 --with-lpsn

# Sharded, which is what CI does: each shard emits raw leads, collate merges,
# diffs and renders ONCE.
python3 bin/curation-triage/tier1_metadata.py --offline --shard 0 --shards-total 8 \
    --emit-leads leads.0.ndjson
python3 bin/curation-triage/collate.py --leads-dir run --shards-total 8 --run-id 2026-08
```

J1b and J2 are off by default on the command line. CI enables all three, because the
repository commits a real datasets cache and a real LPSN cache, so all three replay
offline.

**`--shards-total` is required and it never defaults.** `collate.py` refuses to guess
it. It is the only way collate can tell "this pick produced no lead" from "the shard
that owned this pick never reported". A count of 1 against an 8-shard matrix reads only
`leads.0.ndjson`, finds no file missing, concludes the run was complete, and marks every
pick the other seven shards held as `resolved` — false records appended to an
append-only ledger. So it is passed explicitly everywhere, and it is cross-checked when
the fragments are named individually.

### Tier 2 — targeted confirmation

```bash
# Offline, replaying the committed caches (no network, no tools needed):
python3 bin/curation-triage/tier2_confirm.py --offline --run-id 2026-08

# Refresh: download genomes, run skani and the typers, then evaluate:
python3 bin/curation-triage/tier2_confirm.py --run-id 2026-08 \
    --typer-env-dir <dir of per-tool conda envs>

# skani is cheap and the typers are not, so each can be refreshed on its own:
python3 bin/curation-triage/tier2_confirm.py --run-id 2026-08 --skip-typers

# Narrow the scope:
python3 bin/curation-triage/tier2_confirm.py --offline --only-confounded
```

Thresholds are all overridable: `--species-ani 95 --af-gate 0.5 --derep 99.5
--epsilon 0.2`.

Three refusals are deliberate:

- **Units with no cached evidence produce nothing.** "Not measured" is not a finding, so
  it is reported as a count rather than as indeterminate worklist rows.
- **An `--offline` run whose caches do not cover the plan fails and names what is
  missing** (`--allow-uncached` forces a deliberate partial run). Exiting 0 there would
  be a silent no-op that reads as a clean run: zero leads, zero records, and nothing said
  about why.
- **A refresh preflights the typers that have uncached work before it downloads a genome
  or runs skani**, so a wrong `--typer-env-dir` costs nothing instead of turning every job
  into a miss. Only tools with pending work are checked, because a tool whose jobs are
  all cached runs nothing and its ~3 GB environment need not be present.
  `--allow-unpinned` downgrades the stop to a warning. `--offline` never reaches the
  preflight.

### Rendering the release worklist

```bash
# What the tag job runs; --run-id <YYYY-MM> renders one run instead of the snapshot
python3 bin/curation-triage/worklist_from_ledger.py --tag v3.9.2 \
    --out worklist.tsv --summary-md worklist.md
```

**The worklist is rendered from the committed ledger.** The two obvious alternatives are
both wrong:

- *Read the worklist off disk.* There is no file to read. `worklist.tsv` and
  `worklist.md` are gitignored, because they are a **view** of the ledger, and committing
  a view next to its source lets the two drift.
- *Re-run the pipeline at tag time.* It would diff this month's leads against a ledger
  that already contains them, so every lead would come back `recurring` and the release
  document would announce "no new leads" for a month that had plenty.

The ledger already stores the delta each run computed, so rendering it reproduces the
labels the curator actually saw. The default view is a **snapshot**: the latest record
per lead across both tiers, which matters because Tier 1 and Tier 2 are appended by
different jobs and can carry different run ids. Each row keeps its own record's `run_id`,
so a Tier-1 lead last computed in July is not relabelled as August work. `reviews.tsv` is
applied at render time, so a lead dismissed after the run does not reappear on the
release asset.

## 3. The committed caches

Several files under `src/curation-triage/` look like scratch but are deliberately
committed: the monthly run and the tests read them, which is what makes the whole loop
work with no network. The comment block in `.gitignore` lists them and says not to
ignore them.

| Cache | What it holds |
| --- | --- |
| `nuccore_cache.json` | The Tier-0 NCBI facts per accession |
| `assembly_reports_cache.json` | The distilled J1a rows for our 249 picks |
| `datasets_summary_cache.json` | 170 real NCBI queries, 3208 records, 785 KB |
| `lpsn_cache.json` | 28794 names from LPSN release 2026-08-04 |
| `type_strain_cache.json` | 249/249 picks resolved to their declared type strain |
| `tier2_ani_cache.json` | 157 real skani rows answering all 158 planned comparisons; 4 of them record a pair skani screened out |
| `tier2_typer_cache.json` | 38 real normalised typer results: the 37 planned jobs, plus ECTyper's context-only serotype |

### The datasets cache

Built and committed for real with the pinned `ncbi-datasets-cli 18.34.0`. The query is
narrow on purpose. Asking NCBI for every genome of *Salmonella enterica* returns about
135 000 records; asking with `--reference` returns one. So the 151 tracked-species taxids
that feed J1b are queried **reference-only** — J1b fires only for NCBI-designated
reference genomes, so that is exactly the question it asks — while the 19 wishlist taxids
that feed J2 need the full list and get it, 3059 records between them. The cache records
the mode per taxid, and a reference-only entry cannot silently serve a full-query need.

```bash
# Rebuild it (needs the network and the pinned binary). Deliberately no --released-after.
python3 bin/curation-triage/tier1_metadata.py --refresh-datasets-cache \
    --with-j1b --with-j2 --datasets-bin <env-dir>/datasets/bin/datasets
```

**`--released-after` narrows the rebuild; it does not top the cache up.** The rebuild
**replaces** the cache, so a window leaves a file holding only the J1b genomes released
inside it — the rest are gone. The monthly job asks with no window at all, so a windowed
cache cannot answer it. The cache therefore records the window **per taxid** as well as
the mode, and an `--offline` run whose question is wider than the cached answer **fails**
and names each gap as `cached only since 07/01/2026, need (all)`.

The cost of not having that check is measured: a windowed cache answering an unwindowed
question loses 78 leads and invents one, while the coverage line still reads "covers 170
taxid(s)". Two consequences follow:

- The refresh workflow has **no `released_after` input**. Nothing in CI passes a matching
  window at replay time, so the input would have no correct use. The command-line flag
  stays for a local operator who knows a plain unwindowed rebuild has to follow.
- **J2's wishlist taxids are always fetched unwindowed**, whatever window is given.
  "Does this species have a genome at all?" is a standing fact, not a delta. A window
  would not narrow that answer, it would falsify it: a species whose only genome predates
  the window would read as "still no genome available" and re-raise a lead we already
  have a genome for.

The committed cache is unwindowed, so it answers both a windowed and an unwindowed run,
and a test pins that.

### The LPSN cache

`lpsn_cache.json` holds 28794 names from release `2026-08-04`, so the J2 LPSN feed
replays offline like the others and the monthly job runs it.

Three properties of the real file govern how it is read:

- **LPSN publishes no date column at all.** The header is `genus_name, sp_epithet,
  subsp_epithet, reference, status, authors, address, risk_grp, nomenclatural_type,
  record_no, record_lnk`. The valid-publication date is therefore **derived from the
  authority string**, and the rule matters: cut at `emend.` or `corrig.`, then take the
  last four-digit year in what remains. That gives *E. coli* 1980 (its Approved-Lists
  year, when the name gained standing), *Abiotrophia adiacens* 1995 (the combination, not
  its 1989 basionym), and *Acinetobacter brisouii* 2011 (not its 2026 emendation). A
  naive "last year in the string" disagrees on 2883 of 29745 rows. Every species row
  yields a year.
- **`authority` is really the `authors` column**, and `reference` is a usable fallback.
- **A correct name is the status field's last clause, never the substring `correct`.**
  The substring also matches `...validly published under the ICNP, inappropriate
  correction; misspelling` — 12 rows, 10 of them species — so a substring test would let
  the feeder propose a misspelling as a species worth adding.

**What limits it.** The date is an *authority year*, not a full date, and it is derived
here rather than published by LPSN, so a name with an unusual authority string gets no
year and sorts last. And the cap, not any cutoff, is what keeps the worklist short:
Kalamari lacks about **4000** validly-published species in the genera it tracks, so
`--max-lpsn` (25 by default) is the real gate. Gaps are sorted newest first, so the cap
keeps the most recently named ones.

Refreshing it is a local command, not CI, because the download is account-gated:

```bash
python3 bin/curation-triage/tier1_metadata.py --offline --no-write --with-j2 --with-lpsn \
    --lpsn-csv <downloaded.csv> --lpsn-release <YYYY-MM-DD> \
    --lpsn-cache src/curation-triage/lpsn_cache.json
```

## 4. The sharded refresh

The refresh is where tools run, and it is sharded **by tool**. Not by unit: the
sub-species panels are cliques — eleven *Salmonella* units compared against each other,
four *Listeria* lineages against each other — so a per-unit split hands the same clique
to several runners. The tool axis has no such overlap, and it is the axis that matters
for cost, because a runner installs exactly one environment.

```
plan  →  warm-genomes  →  { 9 typer shards | skani shard | datasets refresh }  →  merge-and-commit
```

- **`plan`** computes the genome-cache key as a hash of
  `refresh_tier2.py --print-accessions` (64 accessions today), so the cache key tracks
  the plan.
- **`warm-genomes` is the only downloader and the only genome-cache writer.** One writer
  means one key and no save races. It fails loudly if any planned genome is still
  missing, so every shard afterwards runs `--offline` and a missing genome is an error
  rather than nine quiet re-downloads.
- **Each shard writes a cache fragment**, and `merge-and-commit` folds the fragments back
  together. The merge is a pure function and it is strict: two fragments that disagree
  about the same key are an error, not last-writer-wins. A silent overwrite there would
  let a shard running a stale environment replace a good number, and that number is what
  every Tier-2 verdict is computed from.
- **The datasets refresh runs once, unsharded, and unwindowed**, because the rebuild
  replaces the file rather than merging into it: shards would clobber each other, and a
  window would leave a cache the monthly job cannot use (§3).

**The merge also refuses to shrink the cache it replaces, and this needs `--base`.** "The
fragments cover every planned job" is a weaker promise than "the merged cache keeps
everything the old one held". The plan lists only the tools that decide a verdict — 37
jobs — while the committed typer cache holds 38 results. The 38th is ECTyper's
context-only serotype, and **ECTyper has no shard**. Without a floor, a full refresh
deletes that result and its recorded tool version, and nothing downstream notices,
because every later check reads only what *is* in the cache.

So both merge steps pass `--base <the cache being replaced>` and write straight over it.
A base is a **floor, not a peer**: it supplies only the keys no fragment answers, so a
shard can still move a call without a conflict, and every carried-forward key is named in
the log. Writing over the real file is what arms the guard — staging through a scratch
path leaves it with nothing to compare against, so it can only warn.

The refresh never passes `--allow-misses`, `--allow-unpinned`, `--allow-partial` or
`--allow-drop`. A failed tool stops the run instead of committing a cache with silently
absent calls, and a lost result stops it instead of committing a cache with silently
absent keys.

**A failed version probe never replaces a version already on record.** The probe returns
`@unknown` when the binary was not found, which is a *missing tool*, not a new one. A
tool version **is** the Tier-2 epoch, so writing `@unknown` over the pins re-baselines
every Tier-2 lead, voids every dismissal in `reviews.tsv`, and stamps false provenance
into append-only records. The probe therefore runs only for tools that actually did work,
and `shigeifinder@unknown` is the one genuine `@unknown` there is (§5). Neither CI
workflow reaches the preflight: both call `tier2_confirm.py --offline`, which runs
nothing.

## 5. The pinned environments

Eleven conda environments, one committed **explicit lockfile** each, under
[`src/curation-triage/envs/`](../src/curation-triage/envs/). An explicit lockfile is a
list of exact package URLs with md5 checksums: installing from it repeats a fixed set of
downloads instead of re-running the solver, which can land on a different build on a
different day.

[`tool_pins.tsv`](../src/curation-triage/tool_pins.tsv) is the index. One row per
**environment**, not per tool key — the two `mlst` tool keys (`mlst_listeria`,
`mlst_yersinia`) run the same binary out of the same prefix and differ only by
`--scheme`, so eleven rows cover twelve tool keys.

| Environment | Package + build | Database | How the database is pinned |
| --- | --- | --- | --- |
| `shigeifinder` | shigeifinder `1.3.4-pyhdfd78af_0` | bundled | the lockfile |
| `btyper3` | btyper3 `3.4.0-pyhdfd78af_0` | bundled + PubMLST copy 2023-07-07 | the lockfile |
| `lissero` | lissero `0.4.10-pyhdfd78af_0` | bundled 2020-06-04 | the lockfile |
| `sistr` | sistr_cmd `1.1.3-pyhdc42f0e_2` | SISTR_V_1.1.3_db, 2025-04-24 | the lockfile |
| `seqsero2` | seqsero2 `1.3.2-pyhdfd78af_0` | bundled | the lockfile |
| `ectyper` | ectyper `2.0.0-pyhdfd78af_4` | alleles bundled; MASH sketch fetched at run time | **not pinned** — see §6 |
| `amrfinderplus` | ncbi-amrfinderplus `4.2.7-hf69ffd2_0` | dated `2026-05-15.1` | not in the package — **frozen by digest**, see §6 |
| `abricate` | abricate `1.4.0-h05cac1d_0` | vfdb 4592 seqs, plasmidfinder 488 seqs, both dated 2026-Apr-3 | the lockfile |
| `mlst` | mlst `2.23.0-hdfd78af_1` | bundled PubMLST, 144 schemes | the lockfile |
| `skani` | skani `0.3.2-h79ce301_0` | none | — |
| `datasets` | ncbi-datasets-cli `18.34.0-ha770c72_0` (conda-forge) | none | — |

The eleven lockfiles together name **1148 package URLs**, all from conda-forge or
bioconda. They total 145 KB, so the pins are committed; the ~13 GB they install are not.

`tool_pins.tsv` and `marker_manifest.tsv` overlap on the nine typers, and a test asserts
they agree, so a bumped build cannot be recorded in one file and forgotten in the other.

### Recreating and verifying them

```bash
# create all eleven, provision the AMRFinderPlus database, then verify every binary
bin/curation-triage/setup_envs.sh --env-dir <dir>

# just check what is already installed — creates nothing, downloads nothing
bin/curation-triage/setup_envs.sh --env-dir <dir> --verify-only

# one environment at a time (this is what CI does, one per runner)
bin/curation-triage/setup_envs.sh --env-dir <dir> --only sistr
```

Verification is not a formality. The script checks each lockfile's sha256 against
`tool_pins.tsv`, then runs every binary **through the prefix's `bin` on `PATH`** and
compares what it reports against the pin.

**Running through `PATH` matters: `abricate` and `mlst` need their environment's `bin`
first on `PATH`.** Both are Perl scripts. Called by absolute path, `abricate` dies with
`Can't locate Path/Tiny.pm` and `mlst` cannot find its Perl libraries — which reads as a
broken install and is really a missing `PATH` entry. LisSero needs the same for
`makeblastdb`.

Two tools cannot be checked the ordinary way, and the script says so rather than
inventing an answer. `shigeifinder --version` exits 2 and prints usage, so its
`version_command` is `none` and its only checkable identity is the conda build.
AMRFinderPlus is checked by `--database_version` as well as `--version`, and a database
that is not `2026-05-15.1` stops the run unless `--allow-db-drift` is given.

A binary that runs but prints **no version-looking token at all** — the ordinary shape of
a broken conda install, for example `error while loading shared libraries: libgomp.so.1`
— is reported by name as `reports '<nothing>' but the pin is '0.3.2'`, and the loop
carries on to the next environment. The hazard behind that, worth respecting when editing
the script: under `set -o pipefail` an empty `grep` exits 1, that status escapes through
the assignment, and `errexit` then kills the whole loop mid-run with no message and no
summary. Three places have this shape — a missing version token, an AMRFinderPlus build
that prints no database-version line, and an unreadable lockfile — and all three must
land in the failure summary instead of ending the run silently.

## 6. Data no lockfile can pin

Three data sources sit outside the conda packages. Each is handled differently, and each
leaves a residue worth reading before trusting it. LPSN is covered in §3; the other two
follow.

### AMRFinderPlus — frozen, and provably the same database

AMRFinderPlus is the only tool here with a first-class dated database version, which is
exactly why it is the *C. botulinum* toxin caller: a changed call is attributable to a
database update rather than to biology. But the database lives **outside** the conda
package. The lockfile installs the software and nothing else, and `amrfinder -u`
downloads the **latest** database, which is not necessarily the pinned `2026-05-15.1`.

So the version string is a **label**, and any directory can be given that name. The
pinned database is therefore **frozen** — an archived copy plus a committed digest
manifest, which turns "is this the same database?" into a question anyone can answer
offline.

- **A committed digest manifest.**
  [`amrfinderplus_db_manifest.tsv`](../src/curation-triage/amrfinderplus_db_manifest.tsv)
  holds one SHA-256 per file for all **162 files (239.4 MB)**, plus a single roll-up
  `2fb16b7f3c57…` over the sorted (name, hash) pairs. 20 KB of evidence in the
  repository; the 239.4 MB stays out of it. The roll-up is defined so that plain
  coreutils reproduces it.
- **`setup_envs.sh` checks the bytes, not the label.** It compares every installed file
  against the manifest and, on a difference, names the files that differ rather than
  saying "mismatch". Measured cost: 239.4 MB verified in 1.6 s cold and 0.8 s warm —
  cheap enough that there is one check and no weaker fast path to get wrong.
  `--allow-db-drift` is the deliberate escape.
- **The gate checks it too.** `tool_pins.tsv` carries a `database_digest` column, and the
  update gate's identity for this database is `2026-05-15.1+db:2fb16b7f3c57`, so a
  database swapped under an unchanged version string is a gate failure and not only a
  setup-script failure. `marker_manifest.tsv` and the Tier-2 caches keep the plain dated
  string on purpose: they record what a tool printed.
- **`typers.tool_env()` sets `CONDA_PREFIX` and `AMRFINDER_DB`** to the pinned prefix,
  and nothing may rely on the ambient values. Without them AMRFinderPlus reads whatever
  database an interactive shell's exported `CONDA_PREFIX` points at — the base conda
  install, say — so the recorded database version would describe a database the pin never
  named. On a runner, where nothing exports it, the same code fails with "No valid
  AMRFinder database is found".

**What still limits it — the sentence a future maintainer needs. A digest proves
*sameness*, not *availability*.** It cannot make NCBI serve `2026-05-15.1` again;
`amrfinder -u` only ever fetches the latest. If the archived copy is lost after NCBI has
moved on, the pinned database is **unrecoverable**, and the manifest can then only prove
that whatever turns up is not it.

So the archive is the thing to look after. It is a 41.9 MiB gzipped tar
(`amrfinderplus-db-2026-05-15.1.tar.gz`, sha256 `83c6c5af9b50…`) kept outside the
repository, and a restored copy reproduces the same roll-up digest. Keep at least two
copies on separate machines or in institutional storage, check their digests on a
schedule, and record where they are somewhere that is not one person's laptop. To
restore:

```bash
tar -xzf amrfinderplus-db-2026-05-15.1.tar.gz \
    -C <env-dir>/amrfinderplus/share/amrfinderplus/data
ln -sfn 2026-05-15.1 <env-dir>/amrfinderplus/share/amrfinderplus/data/latest
bin/curation-triage/setup_envs.sh --env-dir <env-dir> --only amrfinderplus --verify-only
```

Do **not** run `amrfinder -u` on a restored prefix: it replaces the frozen database with
the latest one. Getting a newer database is an update-gate event, not a silent pass.

### ECTyper's sketch — unpinned, and allowed to be

ECTyper fetches a 943 MiB MASH sketch from Zenodo at run time, so nothing pins it — and
nothing needs to, because ECTyper is **context-only**. Task planning drops it from every
decisive call, so it is never a verdict and never a lead. Its cached serotype is attached
to the *E. coli* marker lead as labelled background under `ectyper_*` keys. The sketch is
deliberately not provisioned, and point 1 of the update gate asserts that no call depends
on a run-time-download database.

## 7. The update gate

```bash
python3 bin/curation-triage/update_gate.py --check
```

This is what CI runs. It is offline, runs no tools, and must pass on the committed tree.
It implements the design's four requirements as four real checks.

| Point | What it checks | What failure means |
| --- | --- | --- |
| 1. Immutable digest | every lockfile's sha256, recomputed from disk, matches `tool_pins.tsv`; the dated AMRFinderPlus version matches `marker_manifest.tsv`; no call depends on a run-time-download database | a pinned artifact is not the artifact recorded |
| 2. Validation panel | the committed Tier-2 evidence replayed through the same functions the monthly run uses — 41 units, 75 calls — digested to `panel-replay:<12 hex>`; **plus** the caches' own recorded `skani_version`, `tool_versions` and `amrfinder_database_version` must equal what `tool_pins.tsv` pins | a call moved, or the evidence was not produced by the pinned tools |
| 3. Classification | every moved call is attributed to `software` / `database` / `threshold` / `rule` / `panel`, by asking which of *that call's* components moved | `unattributable`: the same pins produced two different answers |
| 4. Changelog | every component's current identity equals the last `to` in `pin_changelog.tsv`, and each component's rows form an unbroken chain | a pin moved and nobody wrote down why |

**Thirty-eight components are pinned**: 11 software, 9 database, 8 rule, 9 threshold, and
`panel:scope`. The threshold set includes `threshold:tier2_ani` — the species boundary,
the aligned-fraction gate, the dereplication line, epsilon, and skani's own screen
floors. Without it, a moved species boundary would be an unattributable difference, which
is the exact failure the gate exists to catch.

Three properties make the gate mean something.

**The validation panel is not a fixture.** It is the committed Tier-2 evidence replayed
offline through the same functions the monthly run uses.

**Point 2 also checks that the caches' own recorded tool versions are the pinned ones.**
Replaying the caches proves the calls repeat, but not that the pinned tools produced
them; a panel replayed from the previous pin's evidence signs off on nothing. The one
recorded exception is handled explicitly: a tool whose `version_command` is `none`
(`shigeifinder`) may only appear as `@unknown`, and a real version string there is
rejected as evidence from somewhere else.

**The panel's scope is frozen in the gate, not read from the ledger.** If the trigger set
came from the ledger and the policy table, a curator resolving a Tier-1 lead would drop a
unit, shrink the panel from 75 calls to 74, move the digest and fail point 2 — ordinary
curation breaking the gate, with no legal changelog row to add. So the scope is frozen as
`PANEL_SCOPE` (21 confounded units + 20 Tier-1 lead units = 41) and is a pinned component
in its own right, which means **changing what the gate replays is a recorded decision**,
written as a `panel` row whose `from` and `to` are the two scope identities. A call that
appeared or disappeared while the scope moved classifies as `panel`; a call present on
*both* sides during a scope change is still `unattributable`, so the new class cannot
become a blanket excuse. The trade is deliberate and worth knowing: new Tier-1 leads and
new confounded units do **not** enter the panel on their own. A gate whose evidence can
move underneath it is not a gate.

The panel call table carries **no ANI number**, deliberately. A skani bump moves the last
decimal of nearly every comparison; a table containing `96.83` would report a difference
on every refresh and teach everyone to ignore the gate. What has to stay put is the call,
not the measurement behind it.

### What to do when the gate fires

The failure message is a recipe, not a complaint. It prints the row to add and the
commands that produce its numbers.

```bash
# 1. on the current tree, before touching anything
python3 bin/curation-triage/update_gate.py --panel --out before.tsv

# 2. change the pin (new lockfile, new database, new threshold, new rule), then
python3 bin/curation-triage/update_gate.py --panel --out after.tsv

# 3. classify every call that moved
python3 bin/curation-triage/update_gate.py --compare before.tsv after.tsv

# 4. add the printed row to src/curation-triage/pin_changelog.tsv, then
python3 bin/curation-triage/update_gate.py --check
```

`before.tsv` is not a committed file, because it is a property of the tree you are
leaving.

The order matters, and the tool assumes it. At step 3 the pin has already moved but its
changelog row does not exist yet, so `--compare` reconstructs the earlier state as **the
last identity each component was recorded at** — the last `to` in `pin_changelog.tsv`,
not the `from` of that row. Reading `from` there makes a legitimate skani bump report
"pinned components moved: 0" and then fail the very call it has just explained, as
unattributable. A component whose last row records `to = '-'` is treated as removed and
drops out.

**If step 3 reports `unattributable`, stop.** It means a call changed while nothing it
depends on moved, and the usual cause is a cache refreshed without its tool version
changing. Adding a changelog row does not fix that; finding out what actually moved does.
`unattributable` is a **failure, not a fifth category**.

If instead it reports `panel`, the panel's *scope* moved — you added or removed a unit
the gate replays. That is a legal change and it gets a `panel:scope` row, but it is a
decision about what the gate validates, not a tool difference, and it should be reviewed
as one.

One consequence is worth knowing in advance: **a legitimate refresh that moves a call
will fail the gate**, because the panel digest changes before any changelog row accepts
it. That is intended. So the refresh workflow uploads the merged caches as an artifact
*before* it runs the gate, and the operator lands them after adding the row. Refresh →
record → re-run.

## 8. CI

Three fork-owned workflows.

**`curation-triage.yml` — the monthly run.** It fires on `schedule`, on
`workflow_dispatch`, and on `v*.*.*` tags. It is **non-gating**, and what makes it
non-gating is that it is kept **out of the branch's required status checks** — *not* a
blanket `continue-on-error`. A real infrastructure failure fails its shard loudly but
without blocking a merge, and collate prints `shards_failed=N/total`. It reuses
`unit-testing.yml`'s sharded-matrix pattern:
`generate-matrix → tier1-shard ×8 → collate → tier2-replay`, then **one** ledger commit
carrying both tiers, via fetch-rebase-push.

**`curation-triage-refresh.yml` — the cache refresh.** `workflow_dispatch` only, sharded
by tool (§4). Genomes are cached by accession, keyed on the hash of the plan's accession
list, and warmed by a single downloader job. Environment prefixes are large — mlst 2.1
GB, abricate 1.9 GB, seqsero2 1.6 GB — so each is restored from the Actions cache, keyed
on its lockfile.

**`curation-triage-tests.yml` — the blocking gate.** It runs pytest and the update gate as
two separate jobs on push and pull request, because they ask different questions:
**pytest asks whether the code behaves; the gate asks whether the pinned world moved and
whether the move is recorded.**

```bash
python3 -m pytest                                  # 531 checks, ~25 s, fully offline
python3 bin/curation-triage/update_gate.py --check
```

**The release asset.** A `v*.*.*` tag attaches `worklist.tsv` and `worklist.md` to the
GitHub Release, because run artifacts expire after 90 days and a curator looking for a
worklist a year later needs it to still be there. The release step passes `files` only —
never `generate_release_notes` — so it cannot rewrite the notes upstream's own release
job authors. A tag runs *only* that job. The worklist is rendered from the ledger (§2).

**A shared `concurrency` group cannot de-conflict the tag path with upstream's
`build-sketch` job.** `concurrency` serialises a workflow against *itself*, upstream's
workflow declares none, and making them serialise would mean editing an upstream file.
The release step tolerates the race instead: it waits for the release object to appear
and updates it, rather than creating a second one.

**Keepalive.** GitHub disables a scheduled workflow after 60 days without repository
activity. A weekly job re-enables it through the API with least-privilege
`permissions: actions: write`, and commits no heartbeat file. The honest limit is in the
job's own comment: it prevents the disabled state, it cannot escape it, because once
schedules stop firing the weekly cron stops with them.

## 9. Tool behaviour a naive integration gets wrong

Each of these produces a confident wrong answer if you trust a tool's documentation
instead of its output. The code that handles each one is named, because that is where the
fix lives if a tool version changes.

| Tool | The behaviour | Handled in |
| --- | --- | --- |
| **AMRFinderPlus** | Names the column **`Element symbol`**, not `Gene symbol`, and writes symbols with an underscore (`bont_E3`). | [typers.py:388](../bin/curation-triage/typers.py#L388) |
| **skani** | Drops a distant pair and returns **no row**, which reads as "not measured" when it means "very distant". The screen is widened and a screened-out pair is recorded as a result. | [skani_runner.py:43](../bin/curation-triage/skani_runner.py#L43) |
| **ShigEiFinder** | Needs its `--tmpdir` to **exist already**. Otherwise it dies *after* writing its header, leaving a file that looks like a clean "no results". The default `--tmpdir` is also a relative path. | [typers.py:595](../bin/curation-triage/typers.py#L595) |
| **abricate**, **mlst** | Perl scripts: they need their environment's `bin` first on `PATH`, or they die in a way that reads as a broken install (§5). LisSero needs the same for `makeblastdb`. | [typers.py:516](../bin/curation-triage/typers.py#L516) |
| **BTyper3 3.4** | Uses the **revised** *B. cereus* group nomenclature — `species(ANI)` is the genomospecies (e.g. *mosaicus*) and the familiar name is in `subspecies(ANI)` — and it does **not** report `rpoB`, so `rpoB` is not in its marker list even though it is a published SNP discriminator. It writes to a file, not stdout. | [typers.py](../bin/curation-triage/typers.py), [markers.py](../bin/curation-triage/markers.py) |
| **mlst**, **LisSero**, **ShigEiFinder** | `mlst` prints no header at all; LisSero reports `FULL`/`PARTIAL`/`NONE`; ShigEiFinder has no organism column and its `ipaH` is a bare `+`/`-`. Four of the nine tools write files rather than stdout, each with its own path convention. | [typers.py](../bin/curation-triage/typers.py) |
| Two environments | Need `setuptools<81` added, because modern setuptools removed `pkg_resources`. | the lockfiles |

## 10. What is run for real, and what is only wired

The rule here is to run a tool before believing its documentation, so this separates the
two honestly.

**Run for real (2026-08-04):**

| What | Evidence |
| --- | --- |
| All 11 lockfiles produced from the live prefixes | regenerating them into a throwaway tree reproduces all 11 byte-for-byte |
| Two environments recreated **from the committed lockfiles** | `datasets` and `skani` installed into a scratch prefix and reported their pinned versions |
| The AMRFinderPlus database inside the pinned prefix | `amrfinder --database_version` returns `2026-05-15.1` with `CONDA_PREFIX` and `AMRFINDER_DB` scrubbed |
| The datasets cache | 170 real queries against NCBI, 2m32s, committed |
| Sharded refresh, end to end | 3 typer shards (`abricate`, `mlst_yersinia`, `amrfinderplus`) plus the skani shard run, merged, and the result is **byte-identical** to the committed caches |
| The merge path the workflow takes | the committed typer cache split into the 9 fragments the matrix produces, merged with `--base` straight over the cache: 38 results in, 38 out, byte-identical, and the update gate still green |
| The merge's refusals | a partial fragment set exits 1 naming the uncovered jobs; a conflicting value is an error; the same 9-fragment merge **without** `--base` exits 1 naming `ectyper|GCF_003018575.1` rather than dropping it |
| The windowed-cache gate | two narrowed datasets caches (the current shape, and the legacy shape an older refresh wrote) both make an unwindowed `--offline` run exit 1 naming the gap; the committed cache still passes |
| `setup_envs.sh` failure paths | unknown `--only`, missing environment, tampered lockfile digest, database drift with and without `--allow-db-drift`, and a binary that prints no version token |
| The update gate | `--check` green, including with no `CONDA_PREFIX`, a bare `PATH`, and both proxies pointed at a dead port; failure paths forced for all four points |
| Both tiers offline | 8 Tier-1 shards plus collate, then the Tier-2 replay that produced the committed ledger |
| **Two more bot months, simulated** | the whole monthly pipeline re-run twice on a scratch copy (`2026-09`, `2026-10`): every lead comes back `recurring`, nothing is resolved, no epoch and no evidence fingerprint moves, and the update gate still passes on the grown tree |
| A failed shard | deleting one of the 8 lead fragments makes `collate.py` exit 1 with `shards_failed=1/8`, append nothing and write no worklist |
| The workflows' own step scripts | five extracted from the parsed YAML and executed locally against the pinned environments; every `run:` block in the three fork-owned files parses and passes `bash -n` |
| YAML and shell lint | `actionlint` 1.7.12 with `shellcheck` 0.11.0: zero findings on the three fork-owned workflow files |

**Wired but never executed:**

| What | Why not |
| --- | --- |
| Any real GitHub Actions run | no runner is available here; the step scripts are run locally instead |
| The `amrfinder -u` branch of `setup_envs.sh` | running it downloads the latest database and destroys the pinned `2026-05-15.1` — the state the pin exists to protect |
| 6 of the 9 typer shards (`btyper3`, `lissero`, `sistr`, `seqsero2`, `shigeifinder`, `mlst_listeria`) | minutes of compute each, and they change no committed cache; the three cheapest plus skani prove the fragment and merge path |
| The ECTyper shard | needs the 943 MiB Zenodo sketch, which is deliberately not provisioned; planning it stays behind `--with-ectyper` |
| The LPSN REST API refresh | not implemented (§3); the cache is refreshed by the account-gated local download |
| The keepalive job | fires only on GitHub, on a weekly cron |
| The release-asset job | fires only on a `v*.*.*` tag; the renderer it calls was run for real against the committed ledger (96 rows, both tiers, per-row run ids) — only the upload itself is untested |

The first of these is the honest headline, and it is repeated in the
[design's limits](CURATION_TRIAGE_DESIGN.md#10-honest-limits): **no workflow has ever run
on a real GitHub runner.**

## 11. Files

| Path | What it is |
| --- | --- |
| [`envs/*.linux-64.lock`](../src/curation-triage/envs/) | 11 explicit conda lockfiles, 1148 package URLs |
| [`tool_pins.tsv`](../src/curation-triage/tool_pins.tsv) | the pin index: build, binary, version command, lockfile digest, database kind and digest |
| [`pin_changelog.tsv`](../src/curation-triage/pin_changelog.tsv) | every identity each pinned component has ever had — 38 baseline rows, one per pinned component, plus the AMRFinderPlus database-freeze row |
| [`amrfinderplus_db_manifest.tsv`](../src/curation-triage/amrfinderplus_db_manifest.tsv) | one SHA-256 per file for the frozen database |
| [`setup_envs.sh`](../bin/curation-triage/setup_envs.sh) | creates each environment from its lockfile, provisions databases, verifies versions |
| [`refresh_tier2.py`](../bin/curation-triage/refresh_tier2.py) | the sharded refresh driver: accession list, per-tool shards, strict merge |
| [`update_gate.py`](../bin/curation-triage/update_gate.py) | the four-point gate |
| [`worklist_from_ledger.py`](../bin/curation-triage/worklist_from_ledger.py) | renders the release worklist from the committed ledger |
| [`tests/fixtures/tier2/raw/`](../bin/curation-triage/tests/fixtures/tier2/) | verbatim output from each of the 9 tools |
| [`.github/workflows/`](../.github/workflows/) | the three fork-owned workflows |
