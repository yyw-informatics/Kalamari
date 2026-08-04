# -*- coding: utf-8 -*-
"""The pinned tool environments: lockfiles, digests, and the pin table.

A Tier-2 marker call is only reproducible if the environment that produced it is.
These tests guard the four ways that promise quietly breaks:

  * a lockfile is regenerated but ``tool_pins.tsv`` still records the old digest,
  * a pin is bumped in one table and forgotten in the other,
  * an absolute local path leaks into a committed artifact, so the pins only work
    on the machine that made them, and
  * the AMRFinderPlus database is replaced under an unchanged version string.

That last one is the reason ``amrfinderplus_db_manifest.tsv`` exists.  No lockfile
can capture that database, and ``amrfinder -u`` only ever fetches the LATEST one, so
``2026-05-15.1`` is frozen rather than re-downloadable and its 162 per-file digests
are the only evidence that a machine holds the same database the ledger was built
with.

Everything here reads committed files, except the tests that run setup_envs.sh itself
against a stub environment built in ``tmp_path``, and the one test that checks the
committed manifest against a real installed database -- which skips cleanly when
there is none, because a CI runner has no 240 MB database and must not go red for it.
Nothing runs conda, nothing downloads, and nothing needs a real tool installed, so
the suite stays offline.
"""

import csv
import hashlib
import os
import re
import shutil
import stat
import subprocess

import pytest

import typers

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")
PINS_PATH = os.path.join(SRC, "tool_pins.tsv")
ENVS_DIR = os.path.join(SRC, "envs")
SETUP_SCRIPT = os.path.join(REPO_ROOT, "bin", "curation-triage", "setup_envs.sh")

MANIFEST_PATH = os.path.join(SRC, "marker_manifest.tsv")
DB_MANIFEST_PATH = os.path.join(SRC, "amrfinderplus_db_manifest.tsv")

DATABASE_KINDS = {"bundled", "post-install", "runtime-download", "none"}

# What each database_kind is allowed to offer as evidence for its identity.  The
# vocabulary is small on purpose: every value is a different physical situation, and
# collapsing two of them is how an unpinned database starts looking pinned.
DATABASE_EVIDENCE = {
    "bundled": "in-lockfile",          # travels in the package; lockfile_sha256 pins it
    "runtime-download": "unpinned",    # fetched when the tool runs; nothing pins it
    "none": "n/a",                     # no database at all
}

# Where the pinned environments live is machine-local, the same way the NCBI taxdump
# directories are: an absolute path in a committed test would publish one site's
# filesystem layout and be wrong everywhere else.  So it is read from the environment.
ENV_DIR_VAR = "CURATION_TRIAGE_ENV_DIR"     # the --env-dir setup_envs.sh was given
DB_DIR_VAR = "AMRFINDER_DB"                 # AMRFinderPlus's own variable
DB_UNDER_ENV_DIR = os.path.join("amrfinderplus", "share", "amrfinderplus",
                                "data", "latest")

# How to rebuild a lockfile, quoted verbatim in the staleness message so the fix is
# in the failure rather than in a document somebody has to go and find.
REGEN = ("conda list --explicit --md5 -p <env-dir>/%s "
         "> src/curation-triage/envs/%s.linux-64.lock (keep the comment header), "
         "then update lockfile_sha256 in src/curation-triage/tool_pins.tsv")


def load_pins(path=PINS_PATH):
    """Read ``tool_pins.tsv``; leading ``#`` lines are documentation."""
    with open(path, encoding="utf-8", newline="") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return [r for r in csv.DictReader(lines, delimiter="\t") if (r.get("tool") or "").strip()]


def load_manifest_rows(path=MANIFEST_PATH):
    with open(path, encoding="utf-8", newline="") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return [r for r in csv.DictReader(lines, delimiter="\t") if (r.get("task_id") or "").strip()]


def sha256_of(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def load_db_manifest(path=DB_MANIFEST_PATH):
    """The frozen database's header fields and its per-file rows.

    The header carries ``# name: value`` fields among the prose.  A field name is a
    single token, which is what keeps a sentence containing a colon from being read
    as one.
    """
    fields, body = {}, []
    with open(path, encoding="utf-8", newline="") as fh:
        for line in fh:
            if line.startswith("#"):
                name, sep, value = line[1:].strip().partition(": ")
                if sep and name and " " not in name:
                    fields.setdefault(name, value.strip())
                continue
            body.append(line)
    rows = [r for r in csv.DictReader(body, delimiter="\t")
            if (r.get("path") or "").strip()]
    return fields, rows


def roll_up(rows):
    """The manifest's roll-up digest, computed the way its header defines it.

    Ascending by raw path bytes, and for each file the exact line ``sha256sum``
    prints -- ``<hash><two spaces><path>\\n``.  Sizes are deliberately not part of
    it: they are recorded for humans, and folding them in later would move the
    digest without any file having changed.
    """
    digest = hashlib.sha256()
    for row in sorted(rows, key=lambda r: r["path"].encode("utf-8")):
        digest.update(("%s  %s\n" % (row["sha256"], row["path"])).encode("utf-8"))
    return digest.hexdigest()


def frozen_database_dir():
    """The installed AMRFinderPlus database on THIS machine, or '' if none is configured."""
    direct = os.environ.get(DB_DIR_VAR, "").strip()
    if direct:
        return direct
    env_dir = os.environ.get(ENV_DIR_VAR, "").strip()
    return os.path.join(env_dir, DB_UNDER_ENV_DIR) if env_dir else ""


@pytest.fixture(scope="module")
def pins():
    return load_pins()


@pytest.fixture(scope="module")
def db_manifest():
    return load_db_manifest()


# --------------------------------------------------------------------------- #
#  The pin table itself                                                        #
# --------------------------------------------------------------------------- #
def test_the_table_holds_one_row_per_pinned_environment(pins):
    """Eleven environments: the nine typers plus skani and the datasets CLI."""
    assert len(pins) == 11
    envs = [r["env"] for r in pins]
    assert len(set(envs)) == len(envs), "an environment is pinned twice"
    assert set(envs) == {"shigeifinder", "btyper3", "lissero", "sistr", "seqsero2",
                         "ectyper", "amrfinderplus", "abricate", "mlst",
                         "skani", "datasets"}


def test_every_typers_environment_has_a_pin(pins):
    """typers.py is the consumer, so a tool it can run with no pin is a hole."""
    pinned = {r["env"] for r in pins}
    used = {typers.TOOL_ENVS.get(tool, tool)
            for tools in typers.TASK_TOOLS.values() for tool in tools}
    assert used <= pinned, "unpinned environments: %s" % sorted(used - pinned)


def test_the_two_mlst_tool_keys_share_the_single_mlst_row(pins):
    """Row granularity is per environment: one prefix, one pin, two tool keys."""
    assert typers.TOOL_ENVS[typers.MLST_LISTERIA] == "mlst"
    assert typers.TOOL_ENVS[typers.MLST_YERSINIA] == "mlst"
    assert len([r for r in pins if r["env"] == "mlst"]) == 1


def test_every_row_names_the_binary_typers_would_actually_call(pins):
    """A pin that verifies the wrong executable verifies nothing."""
    by_env = {r["env"]: r for r in pins}
    for tool in (typers.SHIGEIFINDER, typers.BTYPER3, typers.LISSERO, typers.SISTR,
                 typers.SEQSERO2, typers.ECTYPER, typers.AMRFINDERPLUS,
                 typers.ABRICATE, typers.MLST_LISTERIA, typers.MLST_YERSINIA):
        env = typers.TOOL_ENVS.get(tool, tool)
        assert by_env[env]["binary"] == typers.TOOL_BINARIES.get(tool, tool)


def test_database_kind_is_one_of_the_four_documented_values(pins):
    for row in pins:
        assert row["database_kind"] in DATABASE_KINDS, \
            "%s has database_kind %r" % (row["tool"], row["database_kind"])


def test_the_two_environments_that_a_lockfile_cannot_pin_are_recorded_as_such(pins):
    """Naming the gaps is the point; a lockfile that quietly misses a database lies."""
    by_env = {r["env"]: r for r in pins}
    # AMRFinderPlus keeps its dated database outside the conda package, so the version
    # string is a label and the file digests are the evidence behind it.
    assert by_env["amrfinderplus"]["database_kind"] == "post-install"
    assert by_env["amrfinderplus"]["database_version"] == "2026-05-15.1"
    assert re.match(r"^sha256:[0-9a-f]{64}$", by_env["amrfinderplus"]["database_digest"])
    # ECTyper fetches a 943 MiB MASH sketch from Zenodo at run time.  It is
    # context-only, never planned decisively, so the sketch stays unprovisioned --
    # and `unpinned` says so rather than pretending some digest covers it.
    assert by_env["ectyper"]["database_kind"] == "runtime-download"
    assert by_env["ectyper"]["database_digest"] == "unpinned"


def test_every_row_offers_the_database_evidence_that_its_kind_allows(pins):
    """A version string is a label; database_digest is what can be checked.

    Each kind gets one legal value, because the four kinds are four different
    physical situations. Letting a `runtime-download` row claim `in-lockfile` would
    make an unpinnable database look pinned, which is the exact confusion the column
    exists to prevent.
    """
    for row in pins:
        kind, digest = row["database_kind"], row["database_digest"]
        if kind == "post-install":
            assert re.match(r"^sha256:[0-9a-f]{64}$", digest), \
                "%s is post-install, so its digest must be a sha256 roll-up, not %r" % (
                    row["tool"], digest)
        else:
            assert digest == DATABASE_EVIDENCE[kind], \
                "%s has database_kind %s, so database_digest must be %r, not %r" % (
                    row["tool"], kind, DATABASE_EVIDENCE[kind], digest)


# --------------------------------------------------------------------------- #
#  Digests: the table and the lockfiles must not drift apart                   #
# --------------------------------------------------------------------------- #
def test_the_recorded_lockfile_digests_match_the_lockfiles_on_disk(pins):
    """The digest is the update gate's immutable artifact identity."""
    for row in pins:
        path = os.path.join(REPO_ROOT, row["lockfile"])
        assert os.path.exists(path), "%s names a missing lockfile" % row["tool"]
        assert sha256_of(path) == row["lockfile_sha256"], \
            "%s is stale: rerun %s" % (row["lockfile"], REGEN % (row["env"], row["env"]))


def test_every_lockfile_on_disk_belongs_to_a_row_and_the_other_way_round(pins):
    on_disk = sorted(f for f in os.listdir(ENVS_DIR) if f.endswith(".linux-64.lock"))
    named = sorted(os.path.basename(r["lockfile"]) for r in pins)
    assert on_disk == named


def test_the_lockfile_path_is_repo_relative_and_named_after_its_environment(pins):
    for row in pins:
        assert row["lockfile"] == "src/curation-triage/envs/%s.linux-64.lock" % row["env"]


# --------------------------------------------------------------------------- #
#  Lockfile contents                                                           #
# --------------------------------------------------------------------------- #
def test_every_lockfile_is_a_conda_explicit_file(pins):
    """`conda create --file` reads exactly this shape; anything else fails at install."""
    for row in pins:
        with open(os.path.join(REPO_ROOT, row["lockfile"]), encoding="utf-8") as fh:
            lines = [ln.rstrip("\n") for ln in fh]
        assert "@EXPLICIT" in lines, "%s has no @EXPLICIT marker" % row["lockfile"]
        assert "# platform: linux-64" in lines, "%s does not declare its platform" % row["lockfile"]
        body = lines[lines.index("@EXPLICIT") + 1:]
        urls = [ln for ln in body if ln.strip()]
        assert urls, "%s lists no packages" % row["lockfile"]
        for url in urls:
            assert url.startswith("https://"), \
                "%s has a non-URL line after @EXPLICIT: %r" % (row["lockfile"], url)
            # `--md5` appends '#<checksum>', which is what makes the install
            # verifiable rather than merely repeatable.
            assert re.search(r"#[0-9a-f]{32}$", url), \
                "%s lost its md5 on %r" % (row["lockfile"], url)


def test_every_expected_version_matches_the_primary_package_build_in_its_lockfile(pins):
    """The pin must be the build the lockfile installs, not a version read elsewhere."""
    for row in pins:
        with open(os.path.join(REPO_ROOT, row["lockfile"]), encoding="utf-8") as fh:
            text = fh.read()
        stamp = "/%s-%s." % (row["conda_package"], row["conda_build"])
        assert stamp in text, \
            "%s installs no %s: rerun %s" % (row["lockfile"], stamp.strip("/."),
                                             REGEN % (row["env"], row["env"]))
        assert row["conda_build"].split("-")[0] == row["expected_version"], \
            "%s expects version %s but its build is %s" % (
                row["tool"], row["expected_version"], row["conda_build"])


def test_the_datasets_cli_comes_from_conda_forge_not_bioconda(pins):
    """Every other pin is bioconda; getting this one's channel wrong finds nothing."""
    row = next(r for r in pins if r["env"] == "datasets")
    with open(os.path.join(REPO_ROOT, row["lockfile"]), encoding="utf-8") as fh:
        text = fh.read()
    assert "conda-forge/linux-64/ncbi-datasets-cli-18.34.0" in text


# --------------------------------------------------------------------------- #
#  The frozen AMRFinderPlus database                                           #
# --------------------------------------------------------------------------- #
def test_the_database_manifest_has_the_header_fields_that_identify_the_freeze(db_manifest):
    """Those fields are read by setup_envs.sh, so a missing one breaks the check."""
    fields, rows = db_manifest
    for name in ("database_version", "files", "total_bytes", "rollup_sha256", "frozen_on"):
        assert name in fields, "%s declares no `# %s:` field" % (
            os.path.basename(DB_MANIFEST_PATH), name)
    assert re.match(r"^[0-9a-f]{64}$", fields["rollup_sha256"])
    assert re.match(r"^\d{4}-\d{2}-\d{2}$", fields["frozen_on"])
    assert int(fields["files"]) == len(rows)
    assert int(fields["total_bytes"]) == sum(int(r["size_bytes"]) for r in rows)


def test_the_roll_up_digest_is_reproducible_from_the_manifests_own_rows(db_manifest):
    """One short string has to identify the whole database, the same way every time.

    The roll-up is the only part of the freeze that a person will ever quote, so its
    definition has to be executable, not descriptive. Recomputing it here from the
    committed rows proves the header's recipe and the recorded value agree, without
    the 240 MB database being present.
    """
    fields, rows = db_manifest
    assert roll_up(rows) == fields["rollup_sha256"]
    # And it is stable: the same rows in a different order give the same digest,
    # because the recipe sorts before it hashes.
    assert roll_up(list(reversed(rows))) == fields["rollup_sha256"]


def test_the_manifest_rows_are_stored_in_the_byte_order_the_roll_up_assumes(db_manifest):
    """A locale-sorted manifest reads the same and digests differently.

    `LC_ALL=C` is in the recipe for a reason: a locale collation reorders case and
    punctuation, so a manifest regenerated without it would still list every file
    correctly while producing a roll-up nobody can reproduce.
    """
    _fields, rows = db_manifest
    paths = [r["path"] for r in rows]
    assert paths == sorted(paths, key=lambda p: p.encode("utf-8"))
    assert len(set(paths)) == len(paths), "the manifest lists a file twice"


def test_the_manifest_names_only_plain_files_inside_the_database_directory(db_manifest):
    """Every path is resolved against a directory the script is pointed at.

    A row saying `../../bin/amrfinder` or an absolute path would make the check read
    somewhere else entirely, so it would pass while the database was wrong.
    """
    _fields, rows = db_manifest
    for row in rows:
        path = row["path"]
        assert not os.path.isabs(path), "%r is absolute" % path
        assert "/" not in path and "\\" not in path, \
            "%r is not a plain file name in the dated directory" % path
        assert path not in (".", ".."), "%r is not a file" % path
        assert re.match(r"^[0-9a-f]{64}$", row["sha256"]), \
            "%s has no sha256" % path
        assert int(row["size_bytes"]) >= 0


def test_the_pin_table_and_the_database_manifest_record_the_same_frozen_database(pins, db_manifest):
    """Two recordings of one fact; if they disagree, neither can be trusted."""
    fields, rows = db_manifest
    row = next(r for r in pins if r["env"] == "amrfinderplus")
    assert row["database_version"] == fields["database_version"], \
        "tool_pins.tsv pins database %s but the manifest was taken from %s" % (
            row["database_version"], fields["database_version"])
    assert row["database_digest"] == "sha256:" + fields["rollup_sha256"], \
        "tool_pins.tsv records %s but the manifest's roll-up is %s" % (
            row["database_digest"], fields["rollup_sha256"])
    assert rows, "the manifest lists no files, so it proves nothing"


def test_the_committed_manifest_matches_the_frozen_database_on_this_machine(db_manifest):
    """The freeze is only worth anything if it still describes a real database.

    Skipped rather than failed when no database is configured: the 240 MB directory is
    machine-local working state that is never committed, so a CI runner has none and
    must not go red for it. Point it at one with $AMRFINDER_DB (the database directory)
    or $CURATION_TRIAGE_ENV_DIR (the --env-dir given to setup_envs.sh).
    """
    db_dir = frozen_database_dir()
    if not db_dir or not os.path.isdir(db_dir):
        pytest.skip(
            "no AMRFinderPlus database to check: set $%s to the dated database "
            "directory, or $%s to the environment directory setup_envs.sh built "
            "(the database is then <env-dir>/%s). CI has no 240 MB database, which "
            "is why this skips instead of failing." % (
                DB_DIR_VAR, ENV_DIR_VAR, DB_UNDER_ENV_DIR))

    fields, rows = db_manifest
    on_disk = sorted(name for name in os.listdir(db_dir)
                     if os.path.isfile(os.path.join(db_dir, name)))
    recorded = {r["path"]: r["sha256"] for r in rows}
    assert on_disk == sorted(recorded), \
        "the database at %s does not hold the recorded file list: %s" % (
            DB_DIR_VAR, sorted(set(on_disk) ^ set(recorded))[:5])

    differing = []
    for name in on_disk:
        if sha256_of(os.path.join(db_dir, name)) != recorded[name]:
            differing.append(name)
    assert not differing, (
        "%d file(s) of the installed database differ from %s: %s. `amrfinder -u` only "
        "fetches the LATEST database, so this one cannot be re-downloaded -- restore it "
        "from the archive, or move the pin through the update gate."
        % (len(differing), os.path.basename(DB_MANIFEST_PATH), differing[:5]))
    assert roll_up([{"path": n, "sha256": recorded[n]} for n in on_disk]) \
        == fields["rollup_sha256"]


# --------------------------------------------------------------------------- #
#  Agreement with the marker manifest                                          #
# --------------------------------------------------------------------------- #
def test_the_nine_typer_pins_agree_with_the_marker_manifest(pins):
    """One bumped build recorded in one table and not the other is unattributable."""
    by_package = {}
    for row in load_manifest_rows():
        package = (row.get("conda_package") or "").strip()
        if not package:
            continue                      # a gap row names no packaged tool
        by_package.setdefault(package, set()).add(
            ((row.get("conda_build") or "").strip(),
             (row.get("database_version") or "").strip()))

    assert len(by_package) == 9, "expected nine packaged typers, got %s" % sorted(by_package)
    for package, seen in by_package.items():
        assert len(seen) == 1, "%s is pinned inconsistently inside the manifest: %s" % (
            package, sorted(seen))

    pinned = {r["conda_package"]: r for r in pins}
    for package, seen in by_package.items():
        build, database_version = seen.pop()
        assert package in pinned, "%s is in marker_manifest.tsv but not in tool_pins.tsv" % package
        assert pinned[package]["conda_build"] == build, \
            "%s: tool_pins.tsv says %s, marker_manifest.tsv says %s" % (
                package, pinned[package]["conda_build"], build)
        assert pinned[package]["database_version"] == database_version, \
            "%s: tool_pins.tsv says database %s, marker_manifest.tsv says %s" % (
                package, pinned[package]["database_version"], database_version)


def test_skani_and_the_datasets_cli_are_the_two_rows_the_manifest_lacks(pins):
    """They decide nothing about markers, but their versions still move calls."""
    manifest_packages = {(r.get("conda_package") or "").strip()
                         for r in load_manifest_rows()}
    extra = {r["conda_package"] for r in pins} - manifest_packages
    assert extra == {"skani", "ncbi-datasets-cli"}


# --------------------------------------------------------------------------- #
#  No committed artifact may leak a local filesystem path                      #
# --------------------------------------------------------------------------- #
LOCAL_PATH = re.compile(r"(/scicomp/|/home/|/Users/|file:///)")


def test_no_pinned_artifact_leaks_a_local_filesystem_path(pins):
    """Package URLs are remote by design; a prefix path would pin one machine."""
    targets = [PINS_PATH, SETUP_SCRIPT, DB_MANIFEST_PATH] + [
        os.path.join(REPO_ROOT, r["lockfile"]) for r in pins]
    for path in targets:
        with open(path, encoding="utf-8") as fh:
            for number, line in enumerate(fh, start=1):
                assert not LOCAL_PATH.search(line), \
                    "%s:%d leaks a local path: %s" % (os.path.basename(path),
                                                      number, line.strip())


def test_the_setup_script_is_executable_and_says_why_the_zenodo_sketch_is_skipped():
    """The skipped 943 MiB download is a decision, so the script must state it."""
    assert os.access(SETUP_SCRIPT, os.X_OK), "setup_envs.sh is not executable"
    with open(SETUP_SCRIPT, encoding="utf-8") as fh:
        text = fh.read()
    assert "set -euo pipefail" in text
    assert "943 MiB" in text and "context-only" in text
    # `amrfinder -u` fetching the LATEST database is the one silent-drift risk in
    # the whole setup, so the script has to spell it out.
    assert "amrfinder -u" in text and "LATEST" in text


def test_the_setup_script_says_what_freezing_the_database_can_and_cannot_prove():
    """The one thing a future maintainer must not misread.

    A digest proves the database is the SAME one. It cannot make NCBI serve
    2026-05-15.1 again. Someone who reads "frozen" as "recoverable" will not keep the
    archive, and the pinned database will then be gone for good, so the script says it.
    """
    with open(SETUP_SCRIPT, encoding="utf-8") as fh:
        text = fh.read()
    assert "amrfinderplus_db_manifest.tsv" in text
    assert "SAMENESS, not availability" in text
    assert "archive" in text


# --------------------------------------------------------------------------- #
#  setup_envs.sh must survive a broken binary and still report it              #
# --------------------------------------------------------------------------- #
def stub_env(root, env_name, binary, body):
    """A directory that looks enough like a conda prefix for setup_envs.sh.

    Only two things make the script treat a prefix as installed: a ``conda-meta``
    directory and an executable ``bin/<binary>``.  Stubbing both keeps this test
    offline -- no conda, no download, no real tool.
    """
    prefix = root / env_name
    (prefix / "conda-meta").mkdir(parents=True)
    (prefix / "bin").mkdir(parents=True)
    path = prefix / "bin" / binary
    path.write_text("#!/bin/sh\n" + body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return prefix


def run_verify(env_dir, only, *extra):
    bash = shutil.which("bash")
    if bash is None:
        pytest.skip("setup_envs.sh needs bash")
    return subprocess.run(
        [bash, SETUP_SCRIPT, "--env-dir", str(env_dir), "--only", only,
         "--verify-only"] + list(extra),
        cwd=REPO_ROOT, capture_output=True, text=True)


# A stub that reports exactly the pinned software version and the pinned database
# version. The point of the two tests below is that this is not enough: the strings
# are right and the database is still the wrong one.
AMRFINDER_STUB = ('case "$1" in\n'
                  '  --database_version) echo "Database version: 2026-05-15.1" ;;\n'
                  '  *) echo "amrfinder 4.2.7" ;;\n'
                  'esac\n')


def stub_amrfinder_env(root, files):
    """A stub amrfinderplus prefix whose dated database holds ``files``."""
    prefix = stub_env(root, "amrfinderplus", "amrfinder", AMRFINDER_STUB)
    data = prefix / "share" / "amrfinderplus" / "data"
    dated = data / "2026-05-15.1"
    dated.mkdir(parents=True)
    for name, text in files.items():
        (dated / name).write_text(text, encoding="utf-8")
    # A RELATIVE link, the way setup_envs.sh makes one, so the prefix stays movable.
    os.symlink("2026-05-15.1", str(data / "latest"))
    return prefix


# A binary that is present and executable but broken the ordinary conda way: the
# loader message carries no X.Y token at all.
BROKEN_LOADER = (
    'echo "skani: error while loading shared libraries: libgomp.so.1: '
    'cannot open shared object file" >&2\nexit 127\n')


def test_a_binary_whose_output_holds_no_version_token_is_reported_and_does_not_abort_the_run(tmp_path):
    """Regression: `grep | head` with no match used to kill the script silently.

    ``first_version_token`` ends in a pipeline, so under ``set -euo pipefail`` a
    no-match made grep exit 1, pipefail carried that out of the pipeline, and the
    plain assignment tripped errexit.  The script died mid-loop: no FAIL line, no
    summary, and every environment after the broken one left unchecked -- a red
    refresh job whose log does not say which tool failed.
    """
    stub_env(tmp_path, "skani", "skani", BROKEN_LOADER)
    stub_env(tmp_path, "datasets", "datasets", 'echo "datasets version: 18.34.0"\n')

    done = run_verify(tmp_path, "skani,datasets")

    # The broken environment is named, with the empty reading spelled out.
    assert "FAIL  skani: skani reports '<nothing>' but the pin is '0.3.2'" in done.stdout
    # The loop carried on: the environment AFTER the broken one was still checked.
    assert "datasets (18.34.0-ha770c72_0" in done.stdout
    assert "ver   18.34.0" in done.stdout
    # And the run ended the way a failed verification is supposed to end.
    assert "2 environment(s) checked, 0 created, 1 problem(s)" in done.stdout
    assert "The pinned environments are NOT usable as they stand:" in done.stdout
    assert done.returncode == 1


def test_a_database_that_is_not_the_frozen_one_is_reported_by_name(tmp_path, db_manifest):
    """The version string passes and the digest still catches it.

    This is the whole reason the manifest exists. The stub prints exactly
    `2026-05-15.1`, so a check that only asked the tool would call this database
    pinned. The per-file comparison names the files instead of saying "mismatch",
    because "changed AMR.LIB" and "162 of 162 differ" are different problems.
    """
    _fields, rows = db_manifest
    first, second = rows[0]["path"], rows[1]["path"]
    stub_amrfinder_env(tmp_path, {first: "not the frozen file\n"})
    stub_env(tmp_path, "skani", "skani", 'echo "skani 0.3.2"\n')

    done = run_verify(tmp_path, "amrfinderplus,skani")

    assert "db    2026-05-15.1 (pinned)" in done.stdout      # the label agreed
    assert "amrfinderplus_db_manifest.tsv" in done.stdout    # the evidence did not
    assert "changed %s" % first in done.stdout
    assert "missing %s" % second in done.stdout
    assert "%d of %d recorded file(s) differ" % (len(rows), len(rows)) in done.stdout
    # Restoring it is the action, and it has to be in the failure: this database
    # cannot be downloaded again.
    assert "cannot be re-downloaded" in done.stdout
    # The loop carried on and ended the way a failed verification is supposed to.
    assert "ver   0.3.2" in done.stdout
    assert "2 environment(s) checked, 0 created, 1 problem(s)" in done.stdout
    assert done.returncode == 1


def test_allow_db_drift_turns_the_digest_failure_into_a_reported_line(tmp_path, db_manifest):
    """The deliberate escape has to cover the digest, not only the version string."""
    _fields, rows = db_manifest
    stub_amrfinder_env(tmp_path, {rows[0]["path"]: "not the frozen file\n"})

    done = run_verify(tmp_path, "amrfinderplus", "--allow-db-drift")

    assert "allowed by --allow-db-drift" in done.stdout
    assert "1 environment(s) checked, 0 created, 0 problem(s)" in done.stdout
    # A run that waved the database through must not then claim it was verified.
    assert "file for file" not in done.stdout
    assert done.returncode == 0


def test_setup_envs_still_passes_when_every_stubbed_binary_reports_its_pinned_version(tmp_path):
    """The guard against a no-match must not blind the ordinary matching path."""
    stub_env(tmp_path, "skani", "skani", 'echo "skani 0.3.2"\n')
    stub_env(tmp_path, "datasets", "datasets", 'echo "datasets version: 18.34.0"\n')

    done = run_verify(tmp_path, "skani,datasets")

    assert "2 environment(s) checked, 0 created, 0 problem(s)" in done.stdout
    assert "All pinned tools resolve and report the version recorded" in done.stdout
    assert done.returncode == 0


# --------------------------------------------------------------------------- #
#  Numbers quoted about a pinned database must be what the binary reports      #
# --------------------------------------------------------------------------- #
DOC_PATH = os.path.join(REPO_ROOT, "docs", "CURATION_TRIAGE_REPRODUCIBILITY.md")
CHANGELOG_PATH = os.path.join(SRC, "pin_changelog.tsv")

# `mlst --list` from the pinned 2.23.0-hdfd78af_1 build prints 144 scheme names.
# The prefix's db/pubmlst directory holds 145 entries, but one of them (dbases.sh)
# is a helper script, not a scheme -- counting the directory is the off-by-one both
# files used to carry.
MLST_SCHEMES = "144"
SCHEME_COUNT = re.compile(r"(\d+)\s+(?:PubMLST\s+)?schemes")


def test_the_documented_mlst_scheme_count_is_the_one_the_pinned_binary_reports():
    """tool_pins.tsv promises every value was read from the installed environment.

    A count taken from ``ls db/pubmlst`` instead of from ``mlst --list`` breaks that
    promise, and it also hides a real future change of exactly one scheme inside the
    off-by-one.
    """
    for path in (DOC_PATH, CHANGELOG_PATH):
        with open(path, encoding="utf-8") as fh:
            counts = SCHEME_COUNT.findall(fh.read())
        assert counts, "%s no longer states an mlst scheme count" % os.path.basename(path)
        assert set(counts) == {MLST_SCHEMES}, \
            "%s says %s schemes; `mlst --list` reports %s" % (
                os.path.basename(path), sorted(set(counts)), MLST_SCHEMES)
