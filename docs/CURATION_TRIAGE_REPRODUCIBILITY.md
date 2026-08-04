# Curation triage — reproducibility: pinned environments, caches, and the update gate

Every Tier-2 call comes out of eleven conda environments and their databases. Running those
tools once, on one machine, proves the calls are right *here, today*. Reproducing them
**somewhere else, later** is a different problem, and it is the one that decides whether a
Tier-2 verdict is evidence or an anecdote.

Three properties make it evidence:

1. Anyone can rebuild the exact environment that produced a call.
2. CI can re-run the whole loop with no tools and no network, and can refresh the caches on
   demand without holding 13 GB of conda environments at once.
3. When a pin moves, every call that moves with it is explained, in writing, before it lands.

---

## 1. What is pinned

Eleven conda environments, one committed **explicit lockfile** each, under
`src/curation-triage/envs/`. An explicit lockfile is a list of exact package URLs with md5
checksums: installing from it repeats a fixed set of downloads instead of re-running the
solver, which can land on a different build on a different day.

`src/curation-triage/tool_pins.tsv` is the index. One row per **environment**, not per tool
key — the two `mlst` tool keys (`mlst_listeria`, `mlst_yersinia`) run the same binary out of
the same prefix, so eleven rows cover twelve tool keys.

| Environment | Package + build | Database | How the database is pinned |
|---|---|---|---|
| `shigeifinder` | shigeifinder `1.3.4-pyhdfd78af_0` | bundled | the lockfile |
| `btyper3` | btyper3 `3.4.0-pyhdfd78af_0` | bundled + PubMLST copy 2023-07-07 | the lockfile |
| `lissero` | lissero `0.4.10-pyhdfd78af_0` | bundled 2020-06-04 | the lockfile |
| `sistr` | sistr_cmd `1.1.3-pyhdc42f0e_2` | SISTR_V_1.1.3_db, 2025-04-24 | the lockfile |
| `seqsero2` | seqsero2 `1.3.2-pyhdfd78af_0` | bundled | the lockfile |
| `ectyper` | ectyper `2.0.0-pyhdfd78af_4` | alleles bundled; MASH sketch fetched at run time | **not pinned** — see §5 |
| `amrfinderplus` | ncbi-amrfinderplus `4.2.7-hf69ffd2_0` | dated `2026-05-15.1` | not in the package — **frozen by digest**, see §5 |
| `abricate` | abricate `1.4.0-h05cac1d_0` | vfdb 4592 seqs, plasmidfinder 488 seqs | the lockfile |
| `mlst` | mlst `2.23.0-hdfd78af_1` | bundled PubMLST, 144 schemes | the lockfile |
| `skani` | skani `0.3.2-h79ce301_0` | none | — |
| `datasets` | ncbi-datasets-cli `18.34.0-ha770c72_0` (conda-forge) | none | — |

The eleven lockfiles together name **1148 package URLs**, all from conda-forge or bioconda.
They total 145 KB, so the pins are committed; the ~13 GB they install are not.

`tool_pins.tsv` and `marker_manifest.tsv` overlap on the nine typers, and a test asserts they
agree. A bumped build cannot be recorded in one file and forgotten in the other.

### Recreating the environments

```bash
# create all eleven, provision the AMRFinderPlus database, then verify every binary
bin/curation-triage/setup_envs.sh --env-dir <dir>

# just check what is already installed — creates nothing, downloads nothing
bin/curation-triage/setup_envs.sh --env-dir <dir> --verify-only

# one environment at a time (this is what CI does, one per runner)
bin/curation-triage/setup_envs.sh --env-dir <dir> --only sistr
```

Verification is not a formality. The script checks each lockfile's sha256 against
`tool_pins.tsv`, then runs every binary **through the prefix `bin` on `PATH`** and compares
what it reports against the pin. Running through `PATH` matters: `abricate` and `mlst` are
Perl scripts that die with `Can't locate Path/Tiny.pm` when only the binary is called by
absolute path, which reads as a broken install and is really a missing `PATH` entry.

Two tools cannot be checked the ordinary way, and the script says so rather than inventing an
answer. `shigeifinder --version` exits 2 and prints usage, so its `version_command` is `none`
and its only checkable identity is the conda build. AMRFinderPlus is checked by
`--database_version` as well as `--version`, and a database that is not `2026-05-15.1` stops
the run unless `--allow-db-drift` is given.

A binary that runs but prints **no version-looking token at all** — the ordinary shape of a broken
conda install, e.g. `error while loading shared libraries: libgomp.so.1` — is reported by name
as `reports '<nothing>' but the pin is '0.3.2'`, and the loop carries on to the next environment.
The hazard behind that, worth respecting when editing the script: under `set -o pipefail` an empty
`grep` exits 1, that status escapes through the assignment, and `errexit` then kills the whole loop
mid-run with no message and no summary. Three places have this shape — a missing version token, an
AMRFinderPlus build that prints no "database version" line, and an unreadable lockfile — and all
three must land in the failure summary instead of ending the run silently.

## 2. Two ways to run

The monthly job and the refresh do different jobs and must not be confused.

| | Monthly replay | Manual refresh |
|---|---|---|
| Workflow | `curation-triage.yml` | `curation-triage-refresh.yml` |
| Trigger | cron, 06:17 UTC on the 1st | `workflow_dispatch` only |
| Network | none | NCBI genome downloads, NCBI datasets API |
| Tools installed | none | one pinned environment per runner |
| Reads | the committed caches | the tools |
| Writes | `ledger.ndjson` (one commit, both tiers) | the caches |
| Can move a call | no | **yes** |

**The monthly run is a replay.** It reads committed caches and computes verdicts from them:
`generate-matrix` → eight `tier1-shard` jobs (`--offline --with-j1b --with-j2 --with-lpsn`) → `collate` →
`tier2-replay` → one ledger commit carrying both tiers. `--shards-total` is passed explicitly
everywhere, and `collate.py` refuses to guess it: it is the only way collate can tell "this pick
produced no lead" from "the shard that owned this pick never reported", and a shard count that
disagrees with the matrix silently marks every unseen pick `resolved` and corrupts the
append-only ledger.

Because it runs no tools, the monthly job cannot produce a new number. That is the point. A
verdict changes only when the evidence changes, and the evidence changes only in a refresh.

**The refresh is where tools run**, and it is sharded **by tool**. Not by unit: the sub-species
panels are cliques — eleven *Salmonella* units compared against each other, four *Listeria*
lineages against each other — so a per-unit split hands the same clique to several runners.
The tool axis has no such overlap, and it is the axis that matters for cost, since a runner
installs exactly one environment.

```
plan  →  warm-genomes  →  { 9 typer shards | skani shard | datasets refresh }  →  merge-and-commit
```

- `plan` computes the genome-cache key as a hash of `refresh_tier2.py --print-accessions` (64
  accessions today), so the cache key tracks the plan.
- `warm-genomes` is the **only** downloader and the only genome-cache writer. One writer means
  one key and no save races. It fails loudly if any planned genome is still missing, so every
  shard afterwards runs `--offline` and a missing genome is an error rather than nine quiet
  re-downloads.
- Each shard writes a **cache fragment**, and `merge-and-commit` folds the fragments back
  together. The merge is a pure function and it is strict: two fragments that disagree about
  the same key are an error, not last-writer-wins. A silent overwrite there would let a shard
  running a stale environment replace a good number, and that number is what every Tier-2
  verdict is computed from.
- **The merge also refuses to shrink the cache it replaces.** "The fragments cover every
  planned job" is a weaker promise than "the merged cache keeps everything the old one held":
  the plan lists only the tools that decide a verdict (37 jobs), while the committed typer cache
  holds 38 results — the 38th is ECTyper's context-only serotype, and ECTyper has no shard.
  Without a floor, a full refresh deletes that result and its recorded tool version, and nothing
  downstream notices, because every later check reads only what *is* in the cache. Both merge
  steps therefore pass `--base <the cache being replaced>` and write straight over it. A base is
  a **floor, not a peer**: it supplies only the keys no fragment answers, so a shard can still
  move a call without a conflict, and every carried-forward key is named in the log. Writing over
  the real file is what arms the guard — staging through a scratch path leaves it with nothing to
  compare against, so it can only warn.
- The datasets refresh runs **once, unsharded, and unwindowed**, because
  `DatasetsSummaryCache.build` replaces the file rather than merging into it — shards would
  clobber each other, and a window would leave a cache the monthly job cannot use (below).

The refresh never passes `--allow-misses`, `--allow-unpinned`, `--allow-partial` or
`--allow-drop`. A failed tool stops the run instead of committing a cache with silently absent
calls, and a lost result stops it instead of committing a cache with silently absent keys.

**A refresh that cannot find its typer environments stops before it runs anything.**
`tier2_confirm.py` preflights the tools that have uncached work — before any genome downloads and
before skani runs, so a wrong `--typer-env-dir` costs nothing. Only tools with pending work are
checked, because a tool whose jobs are all cached runs nothing and its ~3 GB environment need not
be present. `--allow-unpinned` downgrades the stop to a warning, matching `refresh_tier2.py`'s
shard gate. Behind that gate sits a subtler rule: **a failed version probe never replaces a
version already on record.** A probe returns `@unknown` when the binary was not found, which is a
missing tool, not a new one — and a tool version *is* the Tier-2 epoch, so writing `@unknown` over
the pins re-baselines every Tier-2 lead, voids every `reviews.tsv` dismissal, and stamps false
provenance into append-only records. The probe therefore runs only for tools that actually did
work, and `shigeifinder@unknown` is the one genuine `@unknown` there is. Neither CI workflow
reaches the preflight: both call `tier2_confirm.py --offline`, which runs nothing.

### The datasets cache

`src/curation-triage/datasets_summary_cache.json` is built and committed for real: 170 taxid
queries, 3208 records, 785 KB, built with the pinned `ncbi-datasets-cli 18.34.0`. It is what
lets the monthly job run J1b and J2 offline.

The query is narrow on purpose. Asking NCBI for every genome of *Salmonella enterica* returns
135 000 records; asking with `--reference` returns one. So the 151 tracked-species taxids that
feed J1b are queried `--reference`-only — J1b fires only for NCBI-designated reference genomes,
so that is exactly the question it asks — while the 19 wishlist taxids that feed J2 need the
full list and get it (3059 records between them). The cache records the mode per taxid, and a
reference-only entry cannot silently serve a full-query need.

**`--released-after` is a fetch-time window, and it does not top the cache up.** It narrows what
the rebuild *asks NCBI for*, so the file that lands holds only the J1b genomes released inside
the window — the rest are gone, because `build` replaces the cache. The monthly job replays with
no window at all, so a windowed cache cannot answer it. The cache therefore records the window
**per taxid** as well as the mode, and `--offline` refuses a cache whose window is narrower than
the question being asked, naming each gap as `cached only since 07/01/2026, need (all)`. The cost
of missing that check is measured: a windowed cache answering an unwindowed question loses 78
leads and invents one, while the coverage line still reads "covers 170 taxid(s)". Two
consequences:

- The refresh workflow has **no `released_after` input**. Nothing in CI passes a matching window
  at replay time, so the input has no correct use. The CLI flag stays for a local operator who
  knows a plain unwindowed rebuild has to follow.
- J2's wishlist taxids are **always fetched unwindowed**, whatever window is given. "Does this
  species have a genome at all?" is a standing fact, not a delta; a window would not narrow that
  answer, it would falsify it — a species whose only genome predates the window would read as
  "still no genome available" and re-raise a lead we already have a genome for.

The committed cache is unwindowed, so it answers both a windowed and an unwindowed run, and a
test pins that.

## 3. The update gate

`bin/curation-triage/update_gate.py --check` is what CI runs. It is offline, runs no tools, and
must pass on the committed tree. It implements the four points of design §5 as four real checks.

| Point | What it checks | What failure means |
|---|---|---|
| 1. Immutable digest | every lockfile's sha256, recomputed from disk, matches `tool_pins.tsv`; the dated AMRFinderPlus version matches `marker_manifest.tsv`; no call depends on a `runtime-download` database | a pinned artifact is not the artifact recorded |
| 2. Validation panel | the committed Tier-2 evidence replayed through the same functions the monthly run uses — 41 units, 75 calls — digested to `panel-replay:<12 hex>`; **plus** the caches' own recorded `skani_version` / `tool_versions` / `amrfinder_database_version` must equal what `tool_pins.tsv` pins | a call moved, or the evidence was not produced by the pinned tools |
| 3. Classification | every moved call is attributed to `software` / `database` / `threshold` / `rule` / `panel` by asking which of *that call's* components moved | `unattributable`: the same pins produced two different answers |
| 4. Changelog | every component's current identity equals the last `to` in `pin_changelog.tsv`, and each component's rows form an unbroken chain | a pin moved and nobody wrote down why |

Thirty-eight components are pinned: 11 software, 9 database, 8 rule, 9 threshold, and
`panel:scope`. The threshold set includes `threshold:tier2_ani` — the species boundary, the
aligned-fraction gate, the derep line, epsilon and skani's own screen floors. Without it a moved
species boundary would be an unattributable difference, which is the exact failure the gate
exists to catch.

**`panel:scope` exists because the panel must not be built from curation data.** If its trigger
set came from the ledger and the policy table, an SME resolving a Tier-1 lead would drop a unit,
shrink the panel from 75 calls to 74, move the digest and fail point 2 — ordinary curation
breaking the gate, with no legal changelog row to add. The scope is frozen in the gate as
`PANEL_SCOPE` (21 confounded units + 20 Tier-1 lead units = 41) and is a pinned component in its
own right, so **changing what the gate replays is a recorded decision**, written as a `panel` row
whose `from`/`to` are the two scope identities. A call that appeared or disappeared while the
scope moved classifies as `panel`; a call present on *both* sides during a scope change is still
`unattributable`, so the new class cannot become a blanket excuse. The trade is deliberate and
worth knowing: new Tier-1 leads and new confounded units do **not** enter the panel on their own.
A gate whose evidence can move underneath it is not a gate.

Point 2's provenance half closes the other side of the same door: replaying the caches proves the
calls repeat, but not that the pinned tools produced them. The one recorded exception is handled
explicitly — a tool whose `version_command` is `none` (`shigeifinder`, which exits 2 on
`--version`) may only appear as `@unknown`, and a real version string there is rejected as
evidence from somewhere else.

The panel call table carries **no ANI number**, deliberately. A skani bump moves the last
decimal of nearly every comparison; a table containing `96.83` would report a difference on
every refresh and teach everyone to ignore the gate. What has to stay put is the call, not the
measurement behind it.

### What a curator does when the gate fires

The failure message is a recipe, not a complaint. It prints the row to add and the commands
that produce its numbers.

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

`before.tsv` is not a committed file, because it is a property of the tree you are leaving.

The order matters, and the tool assumes it. At step 3 the pin has already moved but its
changelog row does not exist yet, so `--compare` reconstructs the earlier state as **the last
identity each component was recorded at** — the last `to` in `pin_changelog.tsv`, not the `from`
of that row. Reading `from` there makes a legitimate skani bump report "pinned components moved:
0" and then fail the very call it has just explained, as `unattributable`. A component whose last
row records `to = '-'` is treated as removed and drops out.

If step 3 reports `unattributable`, stop. It means a call changed while nothing it depends on
moved, and the usual cause is a cache refreshed without its tool version changing. Adding a
changelog row does not fix that; finding out what actually moved does.

If instead it reports `panel`, the panel's *scope* moved — you added or removed a unit the gate
replays. That is a legal change and it gets a `panel:scope` row, but it is a decision about what
the gate validates, not a tool difference, and it should be reviewed as one.

One consequence is worth knowing in advance: **a legitimate refresh that moves a call will fail
the gate**, because the panel digest changes before any changelog row accepts it. That is
intended. So the refresh workflow uploads the merged caches as an artifact *before* it runs the
gate, and the operator lands them after adding the row. Refresh → record → re-run.

## 4. Release asset and keepalive

A `v*.*.*` tag attaches `worklist.tsv` and `worklist.md` to the GitHub Release, because run
artifacts expire after 90 days and a curator looking for a worklist a year later needs it to
still be there. The release step passes `files` only — never `generate_release_notes` — so it
cannot rewrite the notes upstream's own release job authors. A tag runs *only* that job.

**The worklist is rendered from the committed ledger**, by
[`worklist_from_ledger.py`](../bin/curation-triage/worklist_from_ledger.py). The two obvious
alternatives are both wrong:

- *Read the worklist off disk.* There is no file to read. `worklist.tsv` and `worklist.md` are
  gitignored, because they are a **view** of the ledger, and committing a view next to its
  source lets the two drift.
- *Re-run the pipeline at tag time.* It would diff this month's leads against a ledger that
  already contains them, so every lead would come back `recurring` and the release document
  would announce "no new leads" for a month that had plenty.

The ledger already stores the delta each run computed, so rendering it reproduces the labels the
curator actually saw — nothing is recomputed and nothing can drift. The default view is a
**snapshot**: the latest record per lead across both tiers, which matters because Tier 1 and
Tier 2 are appended by different jobs and can carry different run ids. Each row keeps its own
record's `run_id`, so a Tier-1 lead last computed in July is not relabelled as August work.
`reviews.tsv` is applied at render time, so a lead an SME dismissed after the run does not
reappear on the release asset.

```bash
# what the tag job runs; --run-id <YYYY-MM> renders one run instead of the snapshot
python3 bin/curation-triage/worklist_from_ledger.py --tag v3.9.2 \
  --out worklist.tsv --summary-md worklist.md
```

GitHub disables a scheduled workflow after 60 days without repository activity. A weekly job
re-enables it through the API with least-privilege `permissions: actions: write`, and commits no
heartbeat file. The honest limit is in the job's own comment: it prevents the disabled state, it
cannot escape it, because once schedules stop firing the weekly cron stops with them.

## 5. Data that no lockfile can pin

Three data sources sit outside the conda packages: LPSN's list of validly published names, the
AMRFinderPlus database, and ECTyper's MASH sketch. Each is handled differently, and each leaves
a residue worth reading before trusting it.

### LPSN — an offline cache on a derived date

The J2 coverage feeder has two halves: the `chromosomes-todo.tsv` wishlist (explicit curator
intent, highest precision) and newly validly-published species from LPSN.
`src/curation-triage/lpsn_cache.json` holds **28794 names from release 2026-08-04**, so
`--with-lpsn` replays offline like the other feeds, and the monthly job runs it.

Three properties of the real file govern how it is read:

- **LPSN publishes no date column at all.** The header is `genus_name, sp_epithet,
  subsp_epithet, reference, status, authors, address, risk_grp, nomenclatural_type, record_no,
  record_lnk`. The valid-publication date is therefore **derived from the authority string**,
  and the rule matters: cut at `emend.`/`corrig.`, then take the last four-digit year in what
  remains. That gives *E. coli* 1980 (its Approved-Lists year, when the name gained standing),
  *Abiotrophia adiacens* 1995 (the combination, not its 1989 basionym), and
  *Acinetobacter brisouii* 2011 (not its 2026 emendation). A naive "last year in the string"
  disagrees on 2883 of 29745 rows. Every species row yields a year.
- **`authority` is really the `authors` column**, and `reference` is a usable fallback.
- **A correct name is the status field's last clause, never the substring `correct`.** The
  substring also matches `...validly published under the ICNP, inappropriate correction;
  misspelling` — 12 rows, 10 of them species — so a substring test lets the feeder propose a
  misspelling as a species worth adding.

**What limits it.** The date is an *authority year*, not a full date, and it is derived here
rather than a field LPSN publishes — a name whose authority string is unusual gets no year and
sorts last. And the cap, not the cutoff, is what keeps the worklist short: Kalamari lacks about
**4000** validly-published species in the genera it tracks, so `--max-lpsn` (25 by default) is
the real gate, and the gaps are sorted newest-first so the cap keeps the most recently named
ones. One consequence to know: **a dismissed lead still occupies a cap slot**, so dismissing the
top 25 does not reveal the next 25 — raise `--max-lpsn` to see further.

**Refreshing it** is a local command, not CI, because the download is account-gated:

```bash
python3 bin/curation-triage/tier1_metadata.py --offline --no-write --with-j2 --with-lpsn \
    --lpsn-csv <downloaded.csv> --lpsn-release <YYYY-MM-DD> \
    --lpsn-cache src/curation-triage/lpsn_cache.json
```

The REST API at `https://api.lpsn.dsmz.de/` is the intended **unattended** path and is not
implemented. It needs the account credentials as a CI secret, a token exchange, and paging over
~34000 records. It slots in as a third opener beside the file and URL ones, with no change to
the distiller, the cache format or the epoch.

### AMRFinderPlus — frozen, and provably the same database

AMRFinderPlus is the only tool here with a first-class dated database version, which is exactly
why it is the *C. botulinum* toxin caller: a changed call is attributable to a database update
rather than to biology. But the database lives **outside** the conda package. The lockfile
installs the software and nothing else, and `amrfinder -u` downloads the **latest** database,
which is not necessarily the pinned `2026-05-15.1`.

So the version string is a **label**, and any directory can be given that name. The pinned
database is therefore **frozen** — an archived copy plus a committed digest manifest, which turns
"is this the same database?" into a question anyone can answer offline:

- **A committed digest manifest.** `src/curation-triage/amrfinderplus_db_manifest.tsv` holds one
  SHA-256 per file for all **162 files (239.4 MB)**, plus a single roll-up
  `2fb16b7f3c57…` over the sorted (name, hash) pairs. 20 KB of evidence in the repository; the
  239.4 MB stays out of it. The roll-up is defined so plain coreutils reproduces it.
- **`setup_envs.sh` checks the bytes, not the label.** It compares every installed file against
  the manifest and, on a difference, names the files that differ rather than saying "mismatch".
  Measured cost: 239.4 MB verified in 1.6 s cold, 0.8 s warm — cheap enough that there is one
  check and no weaker fast path to get wrong. `--allow-db-drift` is the deliberate escape.
- **The gate checks it too.** `tool_pins.tsv` carries a `database_digest` column, and the update
  gate's identity for this database is `2026-05-15.1+db:2fb16b7f3c57`, so a database swapped
  under an unchanged version string is a gate failure, not only a setup-script failure.
  `marker_manifest.tsv` and the Tier-2 caches keep the plain dated string on purpose: they
  record what a tool printed.
- **`typers.tool_env()` sets `CONDA_PREFIX` and `AMRFINDER_DB`** to the pinned prefix, and
  nothing may rely on the ambient values. Without them AMRFinderPlus reads whatever database an
  interactive shell's exported `CONDA_PREFIX` points at — the base conda install, say — so the
  recorded database version describes a database the pin never named. On a runner, where nothing
  exports it, the same code fails with "No valid AMRFinder database is found".

**What still limits it — the sentence a future maintainer needs.** A digest proves *sameness*,
not *availability*. It cannot make NCBI serve `2026-05-15.1` again; `amrfinder -u` only ever
fetches the latest. If the archived copy is lost after NCBI has moved on, the pinned database is
**unrecoverable**, and the manifest can then only prove that whatever turns up is not it.

So the archive is the thing to look after. It is a 41.9 MiB gzipped tar
(`amrfinderplus-db-2026-05-15.1.tar.gz`, sha256 `83c6c5af9b50…`) kept outside the repository, and
a restored copy reproduces the same roll-up digest. Keep at least two copies on separate machines
or in institutional storage, check their digests on a schedule, and record where they are
somewhere that is not one person's laptop. To restore:

```bash
tar -xzf amrfinderplus-db-2026-05-15.1.tar.gz \
    -C <env-dir>/amrfinderplus/share/amrfinderplus/data
ln -sfn 2026-05-15.1 <env-dir>/amrfinderplus/share/amrfinderplus/data/latest
bin/curation-triage/setup_envs.sh --env-dir <env-dir> --only amrfinderplus --verify-only
```

Do **not** run `amrfinder -u` on a restored prefix: it replaces the frozen database with the
latest one. Getting a newer database is an update-gate event, not a silent pass.

### ECTyper's sketch — unpinned, and allowed to be

ECTyper fetches a 943 MiB MASH sketch from Zenodo at run time, so nothing pins it — and nothing
needs to, because ECTyper is **context-only**. Decisive task planning drops it, so it is never a
verdict and never a lead; its cached serotype is attached to the *E. coli* marker lead as
labelled background under `ectyper_*` keys, written before the verdict fields and absent from the
lead fingerprint, so a changed serotype can neither overwrite a call nor re-open a dismissed
lead. The sketch is deliberately not provisioned, and point 1 of the update gate asserts that no
call depends on a `runtime-download` database.

## 6. What is run for real, and what is only wired

The rule here is to run a tool before believing its documentation, so this section separates the
two honestly.

**Run for real (2026-08-04):**

| What | Evidence |
|---|---|
| All 11 lockfiles produced from the live prefixes | regenerating them into a throwaway tree reproduces all 11 byte-for-byte |
| Two environments recreated **from the committed lockfiles** | `datasets` and `skani` installed into a scratch prefix and reported their pinned versions |
| The AMRFinderPlus database inside the pinned prefix | `amrfinder --database_version` returns `2026-05-15.1` with `CONDA_PREFIX` and `AMRFINDER_DB` scrubbed |
| The datasets cache | 170 real queries against NCBI, 2m32s, committed |
| Sharded refresh, end to end | 3 typer shards (`abricate`, `mlst_yersinia`, `amrfinderplus`) + the skani shard run, merged, and the result is **byte-identical** to the committed caches |
| The merge path the workflow takes | the committed typer cache split into the 9 fragments the matrix produces, merged with `--base` straight over the cache: 38 results in, 38 out, byte-identical, and the update gate still green |
| The merge's refusals | a partial fragment set exits 1 naming the uncovered jobs; a conflicting value is an error; the same 9-fragment merge **without** `--base` exits 1 naming `ectyper|GCF_003018575.1` and `meta.tool_versions.ectyper` rather than dropping them |
| The windowed-cache gate | two narrowed datasets caches (the current shape, and the legacy shape an older refresh wrote) both make an unwindowed `--offline` run exit 1 naming the gap; the committed cache still passes |
| `setup_envs.sh` failure paths | unknown `--only`, missing environment, tampered lockfile digest, database drift with and without `--allow-db-drift`, and a binary that prints no version token |
| The update gate | `--check` green, including with no `CONDA_PREFIX`, a bare `PATH`, and both proxies pointed at a dead port; failure paths forced for all four points |
| Both tiers offline | 8 Tier-1 shards + collate, then the Tier-2 replay that produced the committed ledger |
| **Two more bot months, simulated** | the whole monthly pipeline re-run twice on a scratch copy (`2026-09`, `2026-10`), ledger 96 → 272 → 448 records: every lead comes back `recurring`, nothing is resolved, no epoch and no evidence fingerprint moves, and the update gate still passes on the grown tree |
| A failed shard | deleting one of the 8 lead fragments makes `collate.py` exit 1 with `shards_failed=1/8`, append nothing and write no worklist |
| The workflows' own step scripts | five extracted from the parsed YAML and executed locally against the pinned environments; every `run:` block in the three fork-owned files parses and passes `bash -n` |
| YAML and shell lint | `actionlint` 1.7.12 with `shellcheck` 0.11.0: zero findings on the three fork-owned workflow files |

**Wired but never executed:**

| What | Why not |
|---|---|
| Any real GitHub Actions run | no runner is available here; the step scripts are run locally instead |
| The `amrfinder -u` branch of `setup_envs.sh` | running it downloads the latest database and destroys the pinned `2026-05-15.1` — the state the pin exists to protect |
| 6 of the 9 typer shards (`btyper3`, `lissero`, `sistr`, `seqsero2`, `shigeifinder`, `mlst_listeria`) | minutes of compute each and they change no committed cache; the three cheapest plus skani prove the fragment/merge path |
| The ECTyper shard | needs the 943 MiB Zenodo sketch, which is deliberately not provisioned; planning it stays behind `--with-ectyper` |
| The LPSN REST API refresh | not implemented — see §5; the cache is refreshed by the account-gated local download |
| The keepalive job | fires only on GitHub, on a weekly cron |
| The release-asset job | fires only on a `v*.*.*` tag; the renderer it calls was run for real against the committed ledger (96 rows, both tiers, per-row run ids) — only the upload itself is untested |

## 7. Files, and the ledger seed

| Path | What it is |
|---|---|
| `src/curation-triage/envs/*.linux-64.lock` | 11 explicit conda lockfiles, 1148 package URLs |
| `src/curation-triage/tool_pins.tsv` | the pin index: build, binary, version command, lockfile digest, database kind |
| `src/curation-triage/pin_changelog.tsv` | every identity each pinned component has ever had (38 baseline rows, one per pinned component, plus the AMRFinderPlus database-freeze row) |
| `src/curation-triage/datasets_summary_cache.json` | 170 real NCBI queries, so J1b and J2 replay offline |
| `src/curation-triage/lpsn_cache.json` | 28794 names from LPSN release 2026-08-04, so `--with-lpsn` replays offline |
| `bin/curation-triage/setup_envs.sh` | creates each environment from its lockfile, provisions databases, verifies versions |
| `bin/curation-triage/refresh_tier2.py` | the sharded refresh driver: accession list, per-tool shards, strict merge |
| `bin/curation-triage/update_gate.py` | the four-point gate |
| `bin/curation-triage/worklist_from_ledger.py` | renders the release worklist from the committed ledger |
| `.github/workflows/curation-triage-refresh.yml` | `workflow_dispatch` cache refresh |

The committed ledger holds **96 records**: the 21 standing Tier-1 J1a records from `2026-07`, plus
the 75 Tier-2 records from the `2026-08` replay (37 `ani`, 21 `marker`, 17 `panel`; 21 of them
cross-linked to a Tier-1 lead by `confirms_lead`). The Tier-2 append is not idempotent — running
`tier2_confirm.py` twice writes the same 75 records twice — so regenerating them means starting
from a ledger that does not already hold that run's Tier-2 block.

It is a **seed file, and the bot's own commit grows it**: `ledger.partition()` re-records every
surviving lead each run, so one monthly run adds ~176 records (21 J1a + 59 J1b + 21 J2 + 75 T2).
A test must therefore assert rules that survive an append, never the seed's own history: a test
that pins the record count, the exact run ids, or "21 Tier-1 records" turns the *blocking* test
workflow red on the first bot commit. What is asserted instead is that every live lead renders as
exactly one row carrying its own record's run id, and that a Tier-2 run appends only `T2` records
and leaves every Tier-1 record byte-identical. A companion test grows a copy of the ledger by one
bot month and re-checks those rules, and another feeds a deliberately broken ledger so they cannot
rot into vacuous checks. Compaction, if the file ever gets unwieldy, is a separate decision.

## 8. Status and open items

Built and green on the committed tree: the eleven pinned environments, the offline caches, the
sharded refresh, the four-point update gate, the monthly replay of both tiers, and the release
worklist renderer.

Open:

- **The two `gap` marker rows** — the BIGSdb *Yersinia* and *Listeria* cgMLST allele callers, and
  a genome-wide *C. botulinum* group I–IV classifier. Each needs a **new pinned tool**, which is
  a marker-manifest change that moves calls and must pass the update gate on its way in. That is
  a different kind of change from "pin what exists and prove it repeats", which is why the two
  are kept apart.
- **The unattended LPSN refresh** (§5): the REST API path needs credentials as a CI secret, a
  token exchange, and paging.
- **The embedding novelty radar** (design §7) is a research track and never blocks this.

```bash
# the whole suite: 531 checks, no network, no external tools
python3 -m pytest
```
