# Kalamari

## Synopsis

[![DOI](https://zenodo.org/badge/181688179.svg)](https://doi.org/10.5281/zenodo.13900883)

Kalamari is a database of completed and public assemblies, backed by trusted institutions.
These assemblies can be further used in formatted databases such as Kraken or Blast.

Completed assemblies means that you do not have to worry about the database itself being contaminated with "rogue" contigs.
Additionally, most assemblies were obtained by subject matter experts (SMEs) at
Centers for Disease Control and Prevention (CDC).
Those not from CDC come from other trusted institutions or projects such as
FDA-ARGOS.
Most genomes are from species that are either studied or are common contaminants
in the Enteric Diseases Laboratory Branch (EDLB) at CDC.

Kalamari also comes with a custom taxonomy database such as defining
_Shigella_ as a subspecies of _Escherichia coli_
or defining the four lineages of _Listeria monocytogenes_.
These changes have been backed by trusted SMEs in EDLB.

---

## Installation Overview
To start using Kalamari, you'll need to complete the following steps:
1. [Export variables for NCBI API](#export-variables-for-ncbi-api)
2. Install Kalamari dependencies
    - Choose either [Conda](#installation-with-conda) (preferred) or [manual installation](#manual-installation).
3. Download the databases
4. Build and filter the taxonomy directory

---

### Export Variables for NCBI API
NCBI edirect requests run considerably more smoothly when the following environment variables are set:
- `NCBI_API_KEY` 
- `EMAIL` 

Follow [these instructions](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities) to 
obtain an NCBI API key.
This key associates your edirect requests with your username.
Without it, edirect requests might be buggy.
After obtaining an NCBI API key, add it to your environment with

    export NCBI_API_KEY=unique_api_key

where `unique_api_key` is a unique hexadecimal number with characters from 0-9 and a-f.

You should also set your email address in the
`EMAIL` environment variable as edirect tries to guess it, which is an error prone process.
Add this variable to your environment with

    export EMAIL=my@email.address

using your own email address instead of `my@email.address`.

---

### Installation with `conda`

1. Create the Kalamari conda environment, then activate it.
```bash
conda create -n kalamari -c conda-forge -c bioconda kalamari
conda activate kalamari
```

When Kalamari is installed via conda, all scripts are placed on your `$PATH`, and the package data directory is installed inside the conda environment.
With the environment activated, run:
```bash
echo "$CONDA_PREFIX"
```
to see the location of the install and the directories containing the scripts, source files, etc.

2. Download the databases.

This step downloads the reference genome FASTA files for the Kalamari database. Note that this step takes a while to complete.
The databases are downloaded using the information contained in `src/chromosomes.tsv` and `src/plasmids.tsv`.
These files represent the chromosome and plasmid databases, respectively.

To download both the chromosome and plasmid databases with default settings, run:
```bash
downloadKalamari.sh
```

To include incomplete assemblies, please see the download section under
[manual installation](#manual-installation).

Files will output to: `${CONDA_PREFIX}/share/kalamari-<version>/kalamari/`

For more control over database downloads when using a conda installation, such as selecting databases, specifying an output directory, or setting download parameters, see [DOWNLOAD_PL.md](docs/DOWNLOAD_PL.md).

3. Build and filter the taxonomy directory

The taxonomy directory contains a locally generated NCBI taxonomy dump that incorporates Kalamari-specific modifications. It includes filtered `nodes.dmp` and `names.dmp` files representing only the TaxIDs present in the Kalamari database (and their ancestors). This taxonomy is used by downstream tools such as Kraken when building formatted databases.

```bash
buildTaxonomy.sh
filterTaxonomy.sh
```
The taxonomy directory will be located at: `${CONDA_PREFIX}/share/kalamari-<version>/taxonomy/`

4. Congrats! You are done! For instructions using Kalamari with Kraken, Sepia, BLAST, ANI, etc., see [database formatting instructions](docs/DATABASES.md).

---

### Manual installation

Manual installation is viable but less preferred.

1. Clone this repo locally:
```bash
git clone https://github.com/lskatz/Kalamari.git
```

2. Install dependencies:
   - Perl (5.x)
   - `wget` (or `curl`)
     - Debian/Ubuntu: `apt-get install wget`
   - [NCBI Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) (`edirect`, `esearch`, etc.)
     - Install via your package manager
     - Debian/Ubuntu: `apt install ncbi-entrez-direct`
   - [taxonkit](https://github.com/shenwei356/taxonkit/releases)

3. Add the `bin/` directory to your `$PATH`:
```bash
cd Kalamari
export PATH="$PWD/bin:$PATH"
```
Confirm with:
```bash
which downloadKalamari.sh
```
To make this change persistent across sessions, add the `export` line to your shell profile (e.g., `~/.bashrc` or `~/.zshrc`).

4. Download the databases.

This step downloads the reference genome FASTA files for the Kalamari database. Note that this step takes a while to complete.
The databases are downloaded using the information contained in `src/chromosomes.tsv` and `src/plasmids.tsv`.
These files represent the chromosome and plasmid databases, respectively.

Optionally, you can include assemblies that are not complete (i.e., more than one contig per chromosome)
by including `src/chromosomes-incomplete.tsv` by using `KALAMARI_EXPERIMENTAL` as shown below.

To download both the chromosome and plasmid databases with default settings, run:

```bash
downloadKalamari.sh
```

To include incomplete assemblies, include `KALAMARI_EXPERIMENTAL` in the environment.
You can either export it or include it when you execute it like so:

```bash
# either
export KALAMARI_EXPERIMENTAL=1
# or
KALAMARI_EXPERIMENTAL=1 downloadKalamari.sh
```

Files will output to: `<Kalamari_cloned_repo>/share/kalamari-<version>/kalamari/`

For more control over the database download such as selecting databases, specifying an output directory, or setting download parameters, see [DOWNLOAD_PL.md](docs/DOWNLOAD_PL.md).

5. Build and filter the taxonomy directory.

The taxonomy directory contains a locally generated NCBI taxonomy dump that incorporates Kalamari-specific modifications. It includes filtered `nodes.dmp` and `names.dmp` files representing only the TaxIDs present in the Kalamari database (and their ancestors). This taxonomy is used by downstream tools such as Kraken when building formatted databases.

```bash
buildTaxonomy.sh
filterTaxonomy.sh
```
The taxonomy directory will be located at: `<Kalamari_cloned_repo>/share/kalamari-<version>/taxonomy/`

6. Congrats! You are done! For instructions using Kalamari with Kraken, Sepia, BLAST, ANI, etc., see [database formatting instructions](docs/DATABASES.md).


## Database formatting instructions

[How to format and query databases](docs/DATABASES.md)

## Curation triage

An add-on that gives curators a monthly, ranked worklist of newly deposited public
genomes worth reviewing. It answers two questions: has one of our current picks gone bad
(mislabelled, withdrawn, reclassified), and is there a species we should track but don't.

It is **informational, never a pass/fail gate** — it flags and ranks for human review
rather than grading. It optimises for precision, because curator time is the scarce
resource. The monthly run needs no genome downloads and works offline from a committed
cache.

### How it works, in one paragraph

Kalamari's genomes were chosen for trusted provenance, not for being statistically average,
so the triage tracks provenance and correctness. NCBI already publishes a near-nightly
verdict comparing every prokaryotic assembly against its species' type strain; reading
*that verdict on Kalamari's own picks* is nearly free, and it catches the case a
similarity computation cannot see — that a pick was reclassified or withdrawn. So the
design is a cheap metadata check that finds a few candidates, with expensive confirmation
reserved for those few.

### Vocabulary

| Term | Meaning |
| --- | --- |
| **genome-unit** | One curated Kalamari entry. Replicons are grouped by their parent assembly, never by the taxid column. 270 units, 249 of them checkable. |
| **Tier 0** | The committed policy table. Decides which units are checked and which are deliberately excluded (viruses, organelles, eukaryotes). |
| **Tier 1** | The monthly metadata check. Needs no genome downloads. |
| **Tier 2** | Confirmation on the few survivors, plus marker genes for the taxa that genome-wide similarity cannot split. 41 units are in scope. |
| **J1a** | Signal: has one of our picks gone bad? Reads NCBI's verdict. This is the anchor and it runs offline. |
| **J1b** | Signal: has a better trusted genome appeared? Off by default. |
| **J2** | Signal: is there a species we should track but don't? Off by default. |
| **lead** | One flagged item on the worklist. |
| **epoch** | The tooling/schema version of a signal. Bumping it re-opens previously dismissed leads. |

### The documentation

Three documents, one per reader.

| Doc | Who it is for |
| --- | --- |
| [For curators](docs/CURATION_TRIAGE_CURATOR.md) | You have a worklist and want to action it: what gets flagged, how to read a lead, how to record a decision |
| [Design](docs/CURATION_TRIAGE_DESIGN.md) | You are reviewing or taking over the add-on: the architecture, the reasoning, the measurement it rests on, and the honest limits |
| [Running it](docs/CURATION_TRIAGE_RUNNING.md) | You are running a stage, refreshing a cache, or moving a pin: commands, caches, pinned environments, the update gate, CI |

### The code

Eleven entry points live in `bin/curation-triage/`; the rest of the files there are libraries
they import.

| Script | What it does |
| --- | --- |
| `build_policy.py` | Builds the Tier 0 policy table from `src/chromosomes.tsv` |
| `coverage_probe.py` | One-shot: measures how many picks NCBI has a verdict for |
| `tier1_metadata.py` | **The monthly run.** Produces the worklist |
| `collate.py` | Merges the shards when the monthly run is split across CI runners |
| `tier2_confirm.py` | Sequence-level confirmation of the flagged and confounded units |
| `refresh_tier2.py` | Refreshes the Tier 2 caches in shards, one tool at a time, then merges them |
| `build_panels.py` | Regenerates the sub-species ANI panels from the policy table |
| `build_type_strains.py` | Resolves each pick's type strain from NCBI's ANI report |
| `setup_envs.sh` | Recreates each pinned tool environment from its lockfile, then verifies it |
| `update_gate.py` | The check a tool, database, threshold or rule change must pass |
| `worklist_from_ledger.py` | Renders the worklist from the committed ledger (the release asset) |

```bash
# The monthly run, fully offline from the committed cache:
python3 bin/curation-triage/tier1_metadata.py --offline --run-id "$(date -u +%Y-%m)"

# Tier 2, offline from the committed caches:
python3 bin/curation-triage/tier2_confirm.py --offline --run-id "$(date -u +%Y-%m)"

# Did any pinned tool, database, threshold or rule move without being recorded?
python3 bin/curation-triage/update_gate.py --check

# Tests (531 checks, no network needed):
python3 -m pytest
```

### Status

Tiers 0, 1 and 2 are built and tested. Tier 2 runs for real with skani and all nine typers
installed: 41 of 41 triggered units evaluated, with the results committed, so an offline run
reproduces them with no network and no tools. Every tool environment is pinned by a committed
conda lockfile. Both tiers run in the monthly workflow and end in one ledger commit.
Refreshing the caches is a separate on-demand workflow, sharded by tool. A change to any
pinned tool, database, threshold or rule has to pass the update gate.

What is **not** proven, and what is still open, is listed once, in
[the design's honest limits](docs/CURATION_TRIAGE_DESIGN.md#10-honest-limits). Read it
before relying on any of this.

Two operational notes that live here rather than in a doc:

- The monthly workflow (`.github/workflows/curation-triage.yml`) is informational and **must
  be kept out of the branch's required status checks** — a broken shard should go red loudly
  without blocking a merge. The separate test workflow is an ordinary blocking gate.
- Some files under `src/curation-triage/` look like caches but are deliberately committed:
  the monthly run and the tests read them, so the whole loop works with no network. See the
  comment block in `.gitignore`.

## Contributing

Please see [CONTRIBUTING.md](CONTRIBUTING.md)

## Citation

Katz LS, Griswold T, Lindsey RL, Lauer AC, Im MS, Williams G, Halpin JL, Gómez GA, Kucerova Z, Morrison S, Page A, Den Bakker HC, Carleton HA. 2025. "Kalamari: a representative set of genomes of public health concern." _Microbiol Resour Announc_ 14:e00963-24. <https://doi.org/10.1128/mra.00963-24>
