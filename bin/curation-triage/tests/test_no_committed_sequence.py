# -*- coding: utf-8 -*-
"""No genome sequence and no database archive may ever be committed.

This is one of the tool's load-bearing rules, not housekeeping. Tier 2 is the first
tier that needs sequence, and the design's answer is that only the DERIVED results
are committed -- the ANI numbers and the normalised typer output, a few tens of KB --
while the megabytes of FASTA stay local scratch. A single genome slipping into the
history is not removable by a later commit.

The rule used to rest on one path rule, `/.cache/`. That covers the default cache
only: `--genome-cache` takes any path, and a run pointed inside the repository would
leave committable sequence behind. The same goes for the frozen AMRFinderPlus archive
-- the docs tell maintainers to keep copies of a 44 MB tarball, which makes a copy
landing in the repo a plausible accident rather than a hypothetical one.

So these tests check the rule two ways: nothing of that kind is tracked *now*, and
`.gitignore` would catch it *next time*.
"""

import os
import subprocess

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

# Extensions that mean "sequence" or "bulk archive". `.gz` is deliberately absent:
# src/provenance/ ships two committed .gz tables, and a blanket rule would either
# untrack them or force a carve-out that invites more carve-outs.
SEQUENCE_SUFFIXES = (".fna", ".fa", ".fasta", ".fsa", ".fna.gz", ".fasta.gz")
ARCHIVE_SUFFIXES = (".tar.gz", ".tgz", ".tar", ".zip")


def git(*args):
    return subprocess.run(["git", "-C", REPO_ROOT, *args],
                          capture_output=True, text=True, check=False)


def tracked_files():
    out = git("ls-files").stdout
    return [line for line in out.splitlines() if line.strip()]


def is_git_repo():
    return git("rev-parse", "--git-dir").returncode == 0


pytestmark = pytest.mark.skipif(not is_git_repo(),
                                reason="not a git checkout; nothing to inspect")


# --------------------------------------------------------------------------- #
#  Nothing of the kind is tracked now                                          #
# --------------------------------------------------------------------------- #
def test_no_sequence_file_is_tracked():
    bad = [f for f in tracked_files() if f.lower().endswith(SEQUENCE_SUFFIXES)]
    assert not bad, ("sequence must never be committed; found:\n  " + "\n  ".join(bad))


def test_no_bulk_archive_is_tracked():
    bad = [f for f in tracked_files() if f.lower().endswith(ARCHIVE_SUFFIXES)]
    assert not bad, ("bulk archives must never be committed; found:\n  "
                     + "\n  ".join(bad))


# --------------------------------------------------------------------------- #
#  .gitignore would catch the next one                                         #
# --------------------------------------------------------------------------- #
# Paths a real mistake would take: a run writing sequence into the repo, a copy of
# the frozen database archive parked next to the manifest that describes it.
WOULD_BE_MISTAKES = [
    "src/curation-triage/GCF_000283675.1.fna",
    "src/curation-triage/panels/reference.fasta",
    "bin/curation-triage/genomes/GCF_000009045.1.fna",
    "genomes/scratch.fa",
    "amrfinderplus-db-2026-05-15.1.tar.gz",
    "src/curation-triage/amrfinderplus-db-2026-05-15.1.tar.gz",
    "workspace-copy.zip",
]


@pytest.mark.parametrize("path", WOULD_BE_MISTAKES)
def test_gitignore_would_catch_it(path):
    assert git("check-ignore", "-q", path).returncode == 0, (
        "%s would be committable; .gitignore has a hole" % path)


# --------------------------------------------------------------------------- #
#  ...without knocking out anything the tool needs                             #
# --------------------------------------------------------------------------- #
# The two committed .gz tables are the reason *.gz is not in the ignore list. If a
# future tidy-up adds it, this fails before the tables silently leave the build.
MUST_STAY_TRACKED = [
    "src/provenance/assembly-complete.tsv.gz",
    "src/provenance/ncbi_ref.acc.gz",
    "src/curation-triage/amrfinderplus_db_manifest.tsv",
    "src/curation-triage/tier2_ani_cache.json",
    "src/curation-triage/tier2_typer_cache.json",
    "src/curation-triage/type_strain_cache.json",
]


@pytest.mark.parametrize("path", MUST_STAY_TRACKED)
def test_a_needed_file_is_not_ignored(path):
    if not os.path.exists(os.path.join(REPO_ROOT, path)):
        pytest.skip("%s is not present in this checkout" % path)
    assert git("check-ignore", "-q", path).returncode != 0, (
        "%s is ignored but the tool needs it committed" % path)


def test_no_tracked_file_matches_an_ignore_rule():
    """Tracked-and-ignored is the inconsistency that makes a .gitignore misleading."""
    out = git("ls-files", "-ci", "--exclude-standard").stdout.strip()
    assert not out, "these files are tracked but match an ignore rule:\n" + out
