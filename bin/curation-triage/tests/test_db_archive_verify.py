# -*- coding: utf-8 -*-
"""verify_db_archive.sh -- the check that makes a backup of the frozen database real.

The AMRFinderPlus database ``2026-05-15.1`` is the one pinned artifact that cannot be
fetched again: ``amrfinder -u`` only ever downloads the latest, and NCBI does not serve
dated versions.  So a copy of the archive is the only route back to the calls already in
the ledger -- and a copy nobody can check is a backup only in name.

These tests build a **tiny synthetic archive and manifest** rather than touching the real
44 MB one, so they run offline in milliseconds and do not depend on a file that lives
outside the repository.  What they pin is the behaviour that matters: it must say no.
A verifier that only ever passes is worse than none, because it converts an unchecked
backup into a believed one.
"""

import hashlib
import os
import subprocess
import tarfile

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
SCRIPT = os.path.join(REPO_ROOT, "bin", "curation-triage", "verify_db_archive.sh")
REAL_MANIFEST = os.path.join(REPO_ROOT, "src", "curation-triage",
                             "amrfinderplus_db_manifest.tsv")

VERIFIED, FAILED, CANNOT_RUN = 0, 1, 2


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def build_db(tmp_path, files):
    """A stand-in dated directory: flat, like the real one."""
    db = tmp_path / "9999-01-01.1"
    db.mkdir()
    for name, content in files.items():
        (db / name).write_bytes(content)
    return db


def rollup_of(files):
    """The manifest's documented roll-up: sorted `sha256sum` lines, then hashed."""
    lines = "".join("%s  %s\n" % (sha256_bytes(content), name)
                    for name, content in sorted(files.items()))
    return sha256_bytes(lines.encode())


def write_manifest(tmp_path, files, archive_sha, rollup=None, name="m.tsv"):
    path = tmp_path / name
    rows = ["# a synthetic manifest, same header convention as the real one",
            "# database_version: 9999-01-01.1",
            "# files: %d" % len(files),
            "# rollup_sha256: %s" % (rollup or rollup_of(files)),
            "# archive_name: synthetic.tar.gz",
            "# archive_sha256: %s" % archive_sha,
            "path\tsize_bytes\tsha256"]
    for n, content in sorted(files.items()):
        rows.append("%s\t%d\t%s" % (n, len(content), sha256_bytes(content)))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def make_archive(tmp_path, db_dir, name="synthetic.tar.gz"):
    archive = tmp_path / name
    with tarfile.open(archive, "w:gz") as tf:
        tf.add(str(db_dir), arcname=db_dir.name)
    return archive


def run(*args):
    return subprocess.run(["bash", SCRIPT, *args], capture_output=True, text=True)


@pytest.fixture
def good(tmp_path):
    files = {"AMR.LIB": b"reference sequences", "version.txt": b"9999-01-01.1\n",
             "taxgroup.tsv": b"taxgroup\tdata\n"}
    db = build_db(tmp_path, files)
    archive = make_archive(tmp_path, db)
    manifest = write_manifest(tmp_path, files, sha256_file(str(archive)))
    return files, archive, manifest


# --------------------------------------------------------------------------- #
#  It passes on a good copy                                                    #
# --------------------------------------------------------------------------- #
def test_a_good_archive_verifies(good):
    _, archive, manifest = good
    r = run("--manifest", str(manifest), str(archive))
    assert r.returncode == VERIFIED, r.stderr
    assert "VERIFIED" in r.stdout
    assert "162" not in r.stdout          # reads the given manifest, not the real one


def test_the_quick_mode_checks_the_digest_without_unpacking(good):
    _, archive, manifest = good
    r = run("--quick", "--manifest", str(manifest), str(archive))
    assert r.returncode == VERIFIED
    assert "quick" in r.stdout
    assert "unpacked" not in r.stdout


# --------------------------------------------------------------------------- #
#  It says no -- the whole reason it exists                                    #
# --------------------------------------------------------------------------- #
def test_bit_rot_in_the_archive_is_caught(good, tmp_path):
    _, archive, manifest = good
    rotted = tmp_path / "rotted.tar.gz"
    data = bytearray(archive.read_bytes())
    data[len(data) // 2] ^= 0xFF
    rotted.write_bytes(bytes(data))
    r = run("--quick", "--manifest", str(manifest), str(rotted))
    assert r.returncode == FAILED
    assert "not the pinned archive" in r.stderr


def test_an_intact_archive_of_the_wrong_database_is_caught(tmp_path):
    """The case that matters most: it unpacks fine, but it is not our database.

    Reached by recording the tampered archive's own digest, so the cheap first check
    passes and the roll-up has to do the work.
    """
    files = {"AMR.LIB": b"reference sequences", "version.txt": b"9999-01-01.1\n"}
    tampered = dict(files, **{"version.txt": b"tampered\n"})
    db = build_db(tmp_path, tampered)
    archive = make_archive(tmp_path, db)
    # digest of the tampered archive, but the per-file digests of the REAL one
    manifest = write_manifest(tmp_path, files, sha256_file(str(archive)))
    r = run("--manifest", str(manifest), str(archive))
    assert r.returncode == FAILED
    assert "not the frozen database" in r.stderr


def test_a_changed_file_is_named_not_just_counted(tmp_path):
    """'digest mismatch' is not actionable; a file name is."""
    files = {"AMR.LIB": b"reference sequences", "version.txt": b"9999-01-01.1\n"}
    tampered = dict(files, **{"version.txt": b"tampered\n"})
    db = build_db(tmp_path, tampered)
    archive = make_archive(tmp_path, db)
    manifest = write_manifest(tmp_path, files, sha256_file(str(archive)))
    r = run("--manifest", str(manifest), str(archive))
    assert "differs: version.txt" in r.stderr


def test_a_missing_file_in_the_archive_is_caught(tmp_path):
    files = {"AMR.LIB": b"reference sequences", "version.txt": b"9999-01-01.1\n"}
    db = build_db(tmp_path, {"AMR.LIB": files["AMR.LIB"]})     # version.txt dropped
    archive = make_archive(tmp_path, db)
    manifest = write_manifest(tmp_path, files, sha256_file(str(archive)))
    r = run("--manifest", str(manifest), str(archive))
    assert r.returncode == FAILED
    assert "expected 2 files, found 1" in r.stderr


# --------------------------------------------------------------------------- #
#  It distinguishes "failed" from "could not run"                              #
# --------------------------------------------------------------------------- #
def test_a_missing_archive_is_cannot_run_not_failed(good, tmp_path):
    _, _, manifest = good
    r = run("--manifest", str(manifest), str(tmp_path / "absent.tar.gz"))
    assert r.returncode == CANNOT_RUN       # a cron job must tell these apart
    assert "no such file" in r.stderr


def test_no_archive_argument_is_cannot_run(good):
    _, _, manifest = good
    assert run("--manifest", str(manifest)).returncode == CANNOT_RUN


def test_help_works_without_an_archive():
    r = run("--help")
    assert r.returncode == VERIFIED and "verify_db_archive.sh" in r.stdout


# --------------------------------------------------------------------------- #
#  The committed manifest carries what the script needs                        #
# --------------------------------------------------------------------------- #
def read_field(name):
    prefix = "# %s: " % name
    with open(REAL_MANIFEST, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(prefix):
                return line[len(prefix):].strip()
    return ""


def test_the_real_manifest_records_the_archive_digest():
    sha = read_field("archive_sha256")
    assert len(sha) == 64 and all(c in "0123456789abcdef" for c in sha)


def test_the_real_manifest_records_the_archive_name_but_never_a_path():
    """Copies live in storage this repo cannot know about, so only the name is pinned."""
    name = read_field("archive_name")
    assert name.endswith(".tar.gz")
    assert "/" not in name


def test_the_script_is_executable():
    assert os.access(SCRIPT, os.X_OK), "verify_db_archive.sh must be executable"
