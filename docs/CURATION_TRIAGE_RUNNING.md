# Curation triage — running it: stages, caches, pins, and the gate

The operator's document: how to run each stage, refresh the cached evidence, and move a
pin.

Running the tools once, on one machine, proves the calls are right *here, today*.
Reproducing them **somewhere else, later** is what decides whether a verdict is evidence
or an anecdote.

- What the output means: [for curators](CURATION_TRIAGE_CURATOR.md)
- Why any of this is shaped the way it is: [design](CURATION_TRIAGE_DESIGN.md)

## The rules, in one screen

Ten operational rules. Each one exists because breaking it corrupts data silently rather
than failing loudly. The rules about the data itself are in the
[design](CURATION_TRIAGE_DESIGN.md#the-rules-in-one-screen).

| Rule | Why |
| --- | --- |
| `--shards-total` is **required** and never guessed | A wrong count marks every unseen pick `resolved` in an append-only ledger (§2) |
| A failed version probe **never** writes `@unknown` over a recorded version | A tool version *is* the Tier-2 epoch, so it voids every dismissal (§4) |
| `--released-after` narrows a rebuild; it does **not** top the cache up | The rebuild replaces the cache, so a window deletes the rest (§3) |
| The Tier-2 merge must pass `--base` | Without it a full refresh silently deletes the ECTyper result (§4) |
| The release worklist is **rendered from the committed ledger** | Not read off disk, not recomputed — either one lies about the run (§2) |
| The gate's panel scope is **frozen**, not read from the ledger | Otherwise ordinary curation breaks the gate (§7) |
| `unattributable` is a **failure**, not a fifth category | It means the same pins produced two different answers (§7) |
| `abricate` and `mlst` need their env's `bin` first on `PATH` | Off `PATH` they die in a way that reads as a broken install (§5) |
| ShigEiFinder's `--tmpdir` must **exist already** | Otherwise it dies *after* its header, leaving a clean-looking empty file (§9) |
| AMRFinderPlus names its column `Element symbol` | Not `Gene symbol` — reading the wrong column finds nothing (§9) |

**Contents**

1. [Two ways to run](#1-two-ways-to-run) ·
2. [Running each stage](#2-running-each-stage) ·
3. [The committed caches](#3-the-committed-caches) ·
4. [The sharded refresh](#4-the-sharded-refresh) ·
5. [The pinned environments](#5-the-pinned-environments) ·
6. [Data no lockfile can pin](#6-data-no-lockfile-can-pin) ·
7. [The update gate](#7-the-update-gate) ·
8. [CI](#8-ci) ·
9. [Tool quirks](#9-tool-quirks)

---

## 1. Two ways to run

Confusing these is the main way to get a wrong number into the ledger.

| | Monthly replay | Manual refresh |
| --- | --- | --- |
| Workflow | `curation-triage.yml` | `curation-triage-refresh.yml` |
| Trigger | cron, 06:17 UTC on the 1st | `workflow_dispatch` only |
| Network | none | NCBI downloads and API |
| Tools installed | none | one pinned environment per runner |
| Reads | the committed caches | the tools |
| Writes | `ledger.ndjson`, one commit, both tiers | the caches |
| Can move a call | no | **yes** |

**The monthly run is a replay.** Because it runs no tools, it cannot produce a new number.
That is the point: a verdict changes only when the evidence changes, and the evidence
changes only in a refresh.

## 2. Running each stage

Every command takes `--help`; only the flags that carry a trap are explained here.

```bash
# Tier 0 — the policy table. EMAIL raises NCBI's rate limits.
EMAIL=you@example.org python3 bin/curation-triage/build_policy.py
python3 bin/curation-triage/build_policy.py --offline    # rebuild from the committed cache

# The coverage probe. First run streams ~2.7 GB of NCBI reports and distils them.
python3 bin/curation-triage/coverage_probe.py [--offline]

# Tier 1 — the monthly net. J1b and J2 are off by default; CI enables all three.
python3 bin/curation-triage/tier1_metadata.py --offline --run-id 2026-08 \
    [--with-j1b --with-j2 --with-lpsn]

# Tier 1, sharded (what CI does): shards emit raw leads, collate merges and renders ONCE.
python3 bin/curation-triage/tier1_metadata.py --offline --shard 0 --shards-total 8 \
    --emit-leads leads.0.ndjson
python3 bin/curation-triage/collate.py --leads-dir run --shards-total 8 --run-id 2026-08

# Tier 2 — targeted confirmation.
python3 bin/curation-triage/tier2_confirm.py --offline --run-id 2026-08
python3 bin/curation-triage/tier2_confirm.py --run-id 2026-08 --typer-env-dir <dir>

# The release worklist.
python3 bin/curation-triage/worklist_from_ledger.py --tag v3.9.2 \
    --out worklist.tsv --summary-md worklist.md
```

**`--shards-total` is required and never defaults.** `collate.py` refuses to guess it,
because it is the only way to tell "this pick produced no lead" from "the shard that owned
this pick never reported". Work through a wrong count: 1 against an 8-shard matrix reads
only `leads.0.ndjson`, finds no file missing, concludes the run was complete, and marks
every pick the other seven shards held as `resolved` — false records appended to an
append-only ledger.

**Tier 0's taxdump location.** Taxonomy comes from a full NCBI taxdump with the bundled
`src/taxonomy/build` overlay layered on for the synthetic lineages. Rather than passing
`--taxdump` every time, list directories in `.taxdump-search` at the repository root; the
first one holding `nodes.dmp` wins. That file is **gitignored on purpose** — an absolute
path in a committed file would publish one institution's filesystem layout and be wrong
everywhere else — so the repository ships no default search path at all.
`$NCBI_TAXDUMP_SEARCH` does the same for CI.

**Tier 2's three refusals.** A unit with no cached evidence produces **nothing**, because
"not measured" is not a finding. An `--offline` run whose caches do not cover the plan
**fails and names what is missing** (`--allow-uncached` forces a partial run); exiting 0
there would be a silent no-op that reads as a clean run. And a refresh preflights the
typers with uncached work **before** downloading a genome, so a wrong `--typer-env-dir`
costs nothing.

**The release worklist is rendered from the committed ledger.** Both alternatives are
wrong. *Reading it off disk* is impossible — `worklist.tsv` and `worklist.md` are
gitignored, because they are a **view** of the ledger and committing a view next to its
source lets the two drift. *Re-running the pipeline at tag time* would diff this month's
leads against a ledger that already contains them, so every lead would come back
`recurring` and the release would announce "no new leads" for a busy month.

Rendering reproduces the labels the curator actually saw. The default view is a
**snapshot** — the latest record per lead across both tiers — because the two tiers are
appended by different jobs and can carry different run ids. Each row keeps its own
record's `run_id`, so a Tier-1 lead last computed in July is not relabelled as August
work. `reviews.tsv` is applied at render time.

## 3. The committed caches

Several files under `src/curation-triage/` look like scratch but are deliberately
committed: the monthly run and the tests read them, which is what makes the loop work
offline. The `.gitignore` comment block lists them and says not to ignore them.

| Cache | What it holds |
| --- | --- |
| `nuccore_cache.json` | Tier-0 NCBI facts per accession |
| `assembly_reports_cache.json` | The distilled J1a rows for our 249 picks |
| `datasets_summary_cache.json` | 170 real NCBI queries, 3208 records |
| `lpsn_cache.json` | 28794 names from LPSN release 2026-08-04 |
| `type_strain_cache.json` | 249/249 picks resolved to their declared type strain |
| `tier2_ani_cache.json` | 157 skani rows covering all 158 planned comparisons |
| `tier2_typer_cache.json` | 38 typer results: 37 planned jobs plus ECTyper's context serotype |

### The datasets cache

The query is narrow on purpose. Asking NCBI for every genome of *Salmonella enterica*
returns about 135 000 records; asking with `--reference` returns one. So the 151
tracked-species taxids that feed J1b are queried **reference-only** — J1b fires only for
NCBI-designated reference genomes, so that is exactly the question it asks. The 19
wishlist taxids that feed J2 need the full list and get it, 3059 records between them. The
cache records the mode per taxid, so a reference-only entry cannot silently serve a
full-query need.

```bash
# Rebuild it. Deliberately no --released-after.
python3 bin/curation-triage/tier1_metadata.py --refresh-datasets-cache \
    --with-j1b --with-j2 --datasets-bin <env-dir>/datasets/bin/datasets
```

**`--released-after` narrows the rebuild; it does not top the cache up.** The rebuild
**replaces** the cache, so a window leaves only the J1b genomes released inside it — the
rest are gone. The monthly job asks with no window, so a windowed cache cannot answer it.
The cache therefore records the window **per taxid**, and an `--offline` run whose question
is wider than the cached answer **fails**, naming each gap. The cost of not having that
check is measured: a windowed cache answering an unwindowed question loses 78 leads and
invents one, while the coverage line still reads "covers 170 taxid(s)".

Two consequences. The refresh workflow has **no `released_after` input**, because nothing
in CI passes a matching window at replay time, so the input would have no correct use. And
**J2's wishlist taxids are always fetched unwindowed**: "does this species have a genome at
all?" is a standing fact, and a window would not narrow that answer, it would falsify it.

### The LPSN cache

**LPSN publishes no date column at all.** The valid-publication date is therefore
**derived from the authority string**: cut at `emend.` or `corrig.`, then take the last
four-digit year in what remains. That gives *E. coli* 1980 (its Approved-Lists year, when
the name gained standing), *Abiotrophia adiacens* 1995 (the combination, not its 1989
basionym), and *Acinetobacter brisouii* 2011 (not its 2026 emendation). A naive "last year
in the string" disagrees on 2883 of 29745 rows. Two more properties: `authority` is really
the `authors` column, and **a correct name is the status field's last clause, never the
substring `correct`** — the substring also matches `inappropriate correction; misspelling`,
so a substring test would propose a misspelling as a species worth adding.

**What limits it.** The date is an authority year, so a name with an unusual authority
string gets no year and sorts last. And the cap, not any cutoff, keeps the worklist short:
Kalamari lacks about **4000** validly-published species in tracked genera, so `--max-lpsn`
(25 by default) is the real gate, with gaps sorted newest first.

Refreshing is a local command, not CI, because the download is account-gated:

```bash
python3 bin/curation-triage/tier1_metadata.py --offline --no-write --with-j2 --with-lpsn \
    --lpsn-csv <downloaded.csv> --lpsn-release <YYYY-MM-DD> \
    --lpsn-cache src/curation-triage/lpsn_cache.json
```

## 4. The sharded refresh

The refresh is where tools run, and it is sharded **by tool**. Not by unit: the sub-species
panels are cliques — eleven *Salmonella* units against each other, four *Listeria* lineages
against each other — so a per-unit split hands the same clique to several runners. The tool
axis has no such overlap, and it is the axis that matters for cost, because a runner
installs exactly one environment.

```
plan → warm-genomes → { 9 typer shards | skani shard | datasets refresh } → merge-and-commit
```

- **`plan`** hashes `refresh_tier2.py --print-accessions` (64 accessions today) into the
  genome-cache key, so the key tracks the plan.
- **`warm-genomes` is the only downloader and the only genome-cache writer.** One writer
  means one key and no save races. It fails loudly if a planned genome is missing, so every
  shard afterwards runs `--offline`.
- **Each shard writes a fragment**, and the merge folds them back. The merge is pure and
  strict: two fragments disagreeing about one key is an error, not last-writer-wins. A
  silent overwrite would let a shard on a stale environment replace a good number, and that
  number is what every Tier-2 verdict is computed from.
- **The datasets refresh runs once, unsharded and unwindowed**, because the rebuild
  replaces the file rather than merging into it.

**The merge must pass `--base`, or a full refresh silently deletes a result.** "The
fragments cover every planned job" is weaker than "the merged cache keeps everything the
old one held". The plan lists only the tools that decide a verdict — 37 jobs — while the
committed typer cache holds 38. The 38th is ECTyper's context-only serotype, and **ECTyper
has no shard**. Without a floor, a full refresh deletes that result and its recorded tool
version, and nothing downstream notices, because every later check reads only what *is* in
the cache.

A base is a **floor, not a peer**: it supplies only the keys no fragment answers, so a
shard can still move a call, and every carried-forward key is named in the log. Writing
over the real file is what arms the guard — staging through a scratch path leaves it
nothing to compare against.

The refresh never passes `--allow-misses`, `--allow-unpinned`, `--allow-partial` or
`--allow-drop`. A failed tool stops the run rather than committing a cache with silently
absent calls.

**A failed version probe never replaces a version already on record.** The probe returns
`@unknown` when the binary was not found, which is a *missing tool*, not a new one. A tool
version **is** the Tier-2 epoch, so writing `@unknown` over the pins re-baselines every
Tier-2 lead, voids every dismissal in `reviews.tsv`, and stamps false provenance into
append-only records. The probe therefore runs only for tools that did work, and
`shigeifinder@unknown` is the one genuine `@unknown` there is.

## 5. The pinned environments

Eleven conda environments, one committed **explicit lockfile** each, under
[`envs/`](../src/curation-triage/envs/). An explicit lockfile is a list of exact package
URLs with md5 checksums, so installing from it repeats a fixed set of downloads instead of
re-running the solver, which can land on a different build on a different day.

[`tool_pins.tsv`](../src/curation-triage/tool_pins.tsv) is the index, with one row per
**environment**, not per tool key: `mlst_listeria` and `mlst_yersinia` run the same binary
from the same prefix and differ only by `--scheme`, so eleven rows cover twelve tool keys.

| Environment | Package + build | Database |
| --- | --- | --- |
| `shigeifinder` | shigeifinder `1.3.4-pyhdfd78af_0` | bundled |
| `btyper3` | btyper3 `3.4.0-pyhdfd78af_0` | bundled + PubMLST copy 2023-07-07 |
| `lissero` | lissero `0.4.10-pyhdfd78af_0` | bundled 2020-06-04 |
| `sistr` | sistr_cmd `1.1.3-pyhdc42f0e_2` | SISTR_V_1.1.3_db, 2025-04-24 |
| `seqsero2` | seqsero2 `1.3.2-pyhdfd78af_0` | bundled |
| `ectyper` | ectyper `2.0.0-pyhdfd78af_4` | alleles bundled; MASH sketch fetched at run time — **not pinned** (§6) |
| `amrfinderplus` | ncbi-amrfinderplus `4.2.7-hf69ffd2_0` | dated `2026-05-15.1` — **frozen by digest** (§6) |
| `abricate` | abricate `1.4.0-h05cac1d_0` | vfdb 4592 seqs, plasmidfinder 488 seqs, both dated 2026-Apr-3 |
| `mlst` | mlst `2.23.0-hdfd78af_1` | bundled PubMLST, 144 schemes |
| `skani` | skani `0.3.2-h79ce301_0` | none |
| `datasets` | ncbi-datasets-cli `18.34.0-ha770c72_0` | none |

The eleven lockfiles name **1148 package URLs** and total 145 KB, so the pins are
committed; the ~13 GB they install are not. `tool_pins.tsv` and `marker_manifest.tsv`
overlap on the nine typers, and a test asserts they agree.

```bash
bin/curation-triage/setup_envs.sh --env-dir <dir>                 # create all eleven, then verify
bin/curation-triage/setup_envs.sh --env-dir <dir> --verify-only   # check only; creates nothing
bin/curation-triage/setup_envs.sh --env-dir <dir> --only sistr    # one at a time — what CI does
```

Verification is not a formality: the script checks each lockfile's sha256, then runs every
binary **through the prefix's `bin` on `PATH`** and compares what it reports against the
pin.

**Running through `PATH` matters: `abricate` and `mlst` need their environment's `bin`
first on `PATH`.** Both are Perl scripts. Called by absolute path, `abricate` dies with
`Can't locate Path/Tiny.pm` and `mlst` cannot find its libraries — which reads as a broken
install and is really a missing `PATH` entry. LisSero needs the same for `makeblastdb`.

Two tools cannot be checked the ordinary way, and the script says so rather than inventing
an answer. `shigeifinder --version` exits 2 and prints usage, so its only checkable
identity is the conda build. AMRFinderPlus is also checked by `--database_version`, and a
database that is not `2026-05-15.1` stops the run unless `--allow-db-drift` is given.

A binary that runs but prints **no version-looking token at all** — the shape of a broken
conda install — is reported by name and the loop carries on. The hazard behind that is
worth respecting when editing the script: under `set -o pipefail` an empty `grep` exits 1,
that status escapes through the assignment, and `errexit` kills the whole loop with no
message and no summary. Three places have this shape — a missing version token, an
AMRFinderPlus build printing no database-version line, and an unreadable lockfile — and all
three must land in the failure summary instead.

## 6. Data no lockfile can pin

### AMRFinderPlus — frozen, and provably the same database

AMRFinderPlus is the only tool here with a first-class dated database version, which is
exactly why it is the *C. botulinum* toxin caller. But the database lives **outside** the
conda package: the lockfile installs the software and nothing else, and `amrfinder -u`
downloads the **latest** database, not necessarily the pinned one. So the version string is
a **label**, and any directory can be given that name.

The pinned database is therefore **frozen**.
[`amrfinderplus_db_manifest.tsv`](../src/curation-triage/amrfinderplus_db_manifest.tsv)
holds one SHA-256 per file for all **162 files (239.4 MB)**, plus a roll-up
`2fb16b7f3c57…` over the sorted (name, hash) pairs — 20 KB of evidence in the repository,
with the 239.4 MB outside it. `setup_envs.sh` compares every installed file against the
manifest and names any that differ (239.4 MB verified in 1.6 s cold). The gate checks it
too: the update gate's identity is `2026-05-15.1+db:2fb16b7f3c57`, so a database swapped
under an unchanged version string is a gate failure, not only a setup failure. And
`typers.tool_env()` sets `CONDA_PREFIX` and `AMRFINDER_DB` to the pinned prefix, because
without them AMRFinderPlus reads whatever database an interactive shell points at.

**What still limits it: a digest proves *sameness*, not *availability*.** It cannot make
NCBI serve `2026-05-15.1` again. If the archived copy is lost after NCBI has moved on, the
pinned database is **unrecoverable**, and the manifest can then only prove that whatever
turns up is not it.

So the archive is the thing to look after: a 41.9 MiB gzipped tar
(`amrfinderplus-db-2026-05-15.1.tar.gz`, sha256 `83c6c5af9b50…`) kept outside the
repository. Keep at least two copies on separate machines, check their digests on a
schedule, and record where they are somewhere that is not one person's laptop.

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
decisive call, so it is never a verdict and never a lead. The sketch is deliberately not
provisioned, and point 1 of the gate asserts that no call depends on a run-time-download
database.

## 7. The update gate

```bash
python3 bin/curation-triage/update_gate.py --check
```

Offline, runs no tools, and must pass on the committed tree. Four checks:

| Point | What failure means |
| --- | --- |
| **1. Immutable digest** | a pinned artifact is not the artifact recorded |
| **2. Validation panel** | a call moved, or the evidence was not produced by the pinned tools |
| **3. Classification** | `unattributable` — the same pins produced two different answers |
| **4. Changelog** | a pin moved and nobody wrote down why |

Point 1 recomputes every lockfile's sha256 from disk, checks the dated AMRFinderPlus
version against `marker_manifest.tsv`, and asserts no call depends on a run-time-download
database. Point 2 replays the committed Tier-2 evidence — 41 units, 75 calls — through the
same functions the monthly run uses, digests it to `panel-replay:<12 hex>`, *and* checks
that the caches' own recorded tool versions equal what `tool_pins.tsv` pins. Point 3
attributes every moved call by asking which of *that call's* components moved. Point 4
checks each component's identity against the last `to` in `pin_changelog.tsv`.

**Thirty-eight components are pinned**: 11 software, 9 database, 8 rule, 9 threshold, and
`panel:scope`. The thresholds include the species boundary, the aligned-fraction gate, the
dereplication line, epsilon and skani's screen floors — without them a moved species
boundary would be an unattributable difference, the exact failure the gate exists to catch.

Three properties make the gate mean something.

**The validation panel is not a fixture.** It is the committed evidence replayed offline
through the production functions.

**Point 2 also checks provenance.** Replaying caches proves the calls repeat, but not that
the pinned tools produced them; a panel replayed from the previous pin's evidence signs off
on nothing. The one exception is explicit: `shigeifinder`, whose `version_command` is
`none`, may only appear as `@unknown`, and a real version string there is rejected.

**The panel's scope is frozen, not read from the ledger.** Suppose it came from the ledger
and the policy table. A curator resolving a Tier-1 lead would drop a unit, shrink the panel
from 75 calls to 74, move the digest and fail point 2 — ordinary curation breaking the gate,
with no legal changelog row to add. So the scope is frozen as `PANEL_SCOPE` (21 confounded
units plus 20 Tier-1 lead units, 41 in all) and is a pinned component in its own right, so
**changing what the gate replays is a recorded decision**. A call present on *both* sides
during a scope change is still `unattributable`, so the new class cannot become a blanket
excuse. The trade is deliberate: new leads do **not** enter the panel on their own. A gate
whose evidence can move underneath it is not a gate.

The panel table carries **no ANI number**, deliberately. A skani bump moves the last
decimal of nearly every comparison, and a table containing `96.83` would report a
difference on every refresh and teach everyone to ignore the gate.

### When the gate fires

The failure message is a recipe. It prints the row to add and the commands producing its
numbers.

```bash
python3 bin/curation-triage/update_gate.py --panel --out before.tsv   # 1. before touching anything
# 2. change the pin, then
python3 bin/curation-triage/update_gate.py --panel --out after.tsv
python3 bin/curation-triage/update_gate.py --compare before.tsv after.tsv   # 3. classify what moved
# 4. add the printed row to src/curation-triage/pin_changelog.tsv, then --check
```

The order matters. At step 3 the pin has moved but its changelog row does not exist yet, so
`--compare` reconstructs the earlier state as **the last identity each component was
recorded at** — the last `to` in the changelog, not the `from` of that row. Reading `from`
would make a legitimate skani bump report "pinned components moved: 0" and then fail the
very call it just explained.

**If step 3 reports `unattributable`, stop.** A call changed while nothing it depends on
moved, usually a cache refreshed without its tool version changing. Adding a changelog row
does not fix that; finding out what actually moved does. If it reports `panel`, the scope
moved — legal, but a decision about what the gate validates, and it should be reviewed as
one.

One consequence worth knowing in advance: **a legitimate refresh that moves a call will
fail the gate**, because the panel digest changes before any changelog row accepts it. That
is intended. The refresh workflow uploads the merged caches as an artifact *before* running
the gate, and the operator lands them after adding the row. Refresh → record → re-run.

## 8. CI

Three fork-owned workflows.

**`curation-triage.yml` — the monthly run**, on `schedule`, `workflow_dispatch` and
`v*.*.*` tags. It is **non-gating**, and what makes it non-gating is being kept **out of
the branch's required status checks** — *not* a blanket `continue-on-error`. A real
failure fails its shard loudly without blocking a merge, and collate prints
`shards_failed=N/total`. The flow is `generate-matrix → tier1-shard ×8 → collate →
tier2-replay`, ending in **one** ledger commit carrying both tiers.

**`curation-triage-refresh.yml` — the cache refresh**, `workflow_dispatch` only, sharded
by tool (§4). Environment prefixes are large (mlst 2.1 GB, abricate 1.9 GB), so each is
restored from the Actions cache keyed on its lockfile.

**`curation-triage-tests.yml` — the blocking gate.** Two jobs, because they ask different
questions: **pytest asks whether the code behaves; the gate asks whether the pinned world
moved, and whether the move is recorded.**

```bash
python3 -m pytest                                   # 531 checks, ~25 s, fully offline
python3 bin/curation-triage/update_gate.py --check
```

**The release asset.** A `v*.*.*` tag attaches the worklist to the GitHub Release, because
run artifacts expire after 90 days. The step passes `files` only, never
`generate_release_notes`, so it cannot rewrite the notes upstream's own release job
authors. **A shared `concurrency` group cannot de-conflict this with upstream's
`build-sketch` job** — `concurrency` serialises a workflow against *itself*, upstream
declares none, and doing it properly would mean editing an upstream file. The step
tolerates the race instead: it waits for the release object and updates it.

**Keepalive.** GitHub disables a scheduled workflow after 60 days without activity, so a
weekly job re-enables it through the API with least-privilege `actions: write`. The honest
limit is in the job's own comment: it prevents the disabled state but cannot escape it,
because once schedules stop firing the weekly cron stops too.

**What has actually run.** Both tiers replay offline; the sharded refresh has been run end
to end with 3 typer shards plus skani, merged, and the result is byte-identical to the
committed caches; the merge's refusals, the windowed-cache gate, `setup_envs.sh`'s failure
paths and all four gate points have had their failures forced; and two further bot months
were simulated on a scratch copy with nothing resolving and the gate still passing. **No
workflow has ever run on a real GitHub runner** — see the
[design's limits](CURATION_TRIAGE_DESIGN.md#8-honest-limits).

## 9. Tool quirks

Each produces a confident wrong answer if you trust the documentation instead of the
output. The code that handles it is named, because that is where the fix goes if a version
changes.

- **AMRFinderPlus names its column `Element symbol`**, not `Gene symbol`, and writes
  symbols with an underscore (`bont_E3`).
  → [typers.py:388](../bin/curation-triage/typers.py#L388)
- **skani drops a distant pair and returns no row**, which the design widens the screen to
  work around. → [skani_runner.py:43](../bin/curation-triage/skani_runner.py#L43)
- **ShigEiFinder needs its `--tmpdir` to exist already.** Otherwise it dies *after* writing
  its header, leaving a file that looks like a clean "no results". Its default is also a
  relative path. → [typers.py:595](../bin/curation-triage/typers.py#L595)
- **abricate and mlst need their environment's `bin` first on `PATH`** (§5).
  → [typers.py:516](../bin/curation-triage/typers.py#L516)
- **BTyper3 3.4 uses the revised *B. cereus* group nomenclature** — `species(ANI)` is the
  genomospecies (e.g. *mosaicus*), the familiar name is in `subspecies(ANI)` — and it does
  **not** report `rpoB`. It writes to a file, not stdout.
  → [markers.py](../bin/curation-triage/markers.py)
- **Three tools report in shapes nothing else uses**: `mlst` prints no header at all,
  LisSero reports `FULL`/`PARTIAL`/`NONE`, and ShigEiFinder has no organism column with a
  bare `+`/`-` for `ipaH`. Four of the nine write files rather than stdout.
  → [typers.py](../bin/curation-triage/typers.py)
- **Two environments need `setuptools<81`**, because modern setuptools removed
  `pkg_resources`. → the lockfiles
