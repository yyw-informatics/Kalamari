# -*- coding: utf-8 -*-
"""No committed file may carry a developer's absolute filesystem path.

This repository is a fork meant to go upstream, so an absolute path in a committed
file is a leak, not just a portability wart: it publishes where one institution
keeps its data and which account did the work, and it is wrong on every other
machine. `taxonomy.py` once hard-coded two site paths as taxdump search locations;
they are now read from a gitignored `.taxdump-search` file instead, and this test
is what stops them coming back.

Untracked scratch is deliberately out of scope: the genome cache under `.cache/`
holds real tool output whose first column IS an absolute path, and it is never
committed.
"""

import os
import re
import subprocess

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

# A path that names a real user's home or a site's mount point. `/home/` and
# `/Users/` need the trailing segment so a bare mention of the directory itself
# ("under /home") does not trip the guard -- it is a real path we care about.
SITE_PATH = re.compile(rb"(/scicomp/|/home/[a-z]|/Users/[A-Za-z])")

# Files that must be allowed to contain the pattern, because detecting it IS their
# job. Keep this list tiny; every entry is a hole in the guard.
ALLOWED = {
    os.path.join("bin", "curation-triage", "tests", "test_no_absolute_paths.py"),
    os.path.join("bin", "curation-triage", "tests", "test_tool_pins.py"),
    os.path.join("bin", "curation-triage", "tests", "test_update_gate.py"),
}

# Compressed and binary payloads: scanning them finds nothing useful and a false
# match inside compressed bytes would be unactionable.
SKIP_SUFFIXES = (".gz", ".zip", ".pdf", ".pptx", ".png", ".jpg", ".msh", ".pyc")

SKIP_DIRS = {".git", ".cache", "__pycache__", ".pytest_cache"}


def _tracked_files():
    """Every file git would consider part of the repo (tracked + untracked-not-ignored).

    Falls back to a filesystem walk when git is unavailable, so the guard still runs
    in an exported tarball. The walk skips the same directories .gitignore does for
    scratch, which is why SKIP_DIRS exists.
    """
    try:
        out = subprocess.run(  # noqa: S603,S607
            ["git", "-C", REPO_ROOT, "ls-files", "-co", "--exclude-standard"],
            capture_output=True, text=True, check=False, timeout=60)
        if out.returncode == 0 and out.stdout.strip():
            return [p for p in out.stdout.splitlines() if p.strip()]
    except (OSError, subprocess.SubprocessError):
        pass
    found = []
    for dirpath, dirnames, filenames in os.walk(REPO_ROOT):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for name in filenames:
            full = os.path.join(dirpath, name)
            found.append(os.path.relpath(full, REPO_ROOT))
    return found


def test_no_committed_file_contains_an_absolute_developer_path():
    offenders = []
    for rel in _tracked_files():
        if rel in ALLOWED or rel.endswith(SKIP_SUFFIXES):
            continue
        full = os.path.join(REPO_ROOT, rel)
        if not os.path.isfile(full):
            continue
        try:
            with open(full, "rb") as fh:
                blob = fh.read()
        except OSError:
            continue
        for n, line in enumerate(blob.splitlines(), 1):
            m = SITE_PATH.search(line)
            if m:
                offenders.append("%s:%d: %s" % (
                    rel, n, line.decode("utf-8", "replace").strip()[:120]))
                break
    assert not offenders, (
        "committed file(s) contain an absolute developer path -- move the value to "
        "configuration outside version control (see .taxdump-search):\n  "
        + "\n  ".join(offenders))


def test_the_guard_actually_matches_the_paths_it_claims_to(tmp_path):
    # A guard that matches nothing passes forever. Pin the patterns it must catch,
    # and the near-misses it must not.
    assert SITE_PATH.search(b"/scicomp/example-mount/ncbi-taxonomy")
    assert SITE_PATH.search(b'    "/home/jdoe/work/dump",')
    assert SITE_PATH.search(b"/Users/Someone/Desktop/x")
    assert not SITE_PATH.search(b"a relative path: src/curation-triage/envs")
    assert not SITE_PATH.search(b"the repo root, not /home or $HOME")


def test_a_planted_absolute_path_would_be_caught(tmp_path, monkeypatch):
    # Prove the scan is wired to the pattern, without editing a real file.
    planted = tmp_path / "leaky.py"
    planted.write_text('TAXDUMP = "/scicomp/example-mount/taxonomy"\n', encoding="utf-8")
    monkeypatch.setattr(
        "test_no_absolute_paths.REPO_ROOT", str(tmp_path))
    monkeypatch.setattr(
        "test_no_absolute_paths._tracked_files", lambda: ["leaky.py"])
    try:
        test_no_committed_file_contains_an_absolute_developer_path()
    except AssertionError as exc:
        assert "leaky.py:1" in str(exc)
    else:
        raise AssertionError("the guard did not catch a planted absolute path")


def test_the_taxdump_search_file_is_gitignored():
    # The whole fix depends on this file never being committed.
    from taxonomy import SITE_TAXDUMP_FILE
    out = subprocess.run(  # noqa: S603,S607
        ["git", "-C", REPO_ROOT, "check-ignore", SITE_TAXDUMP_FILE],
        capture_output=True, text=True, check=False, timeout=60)
    if out.returncode == 128:          # not a git repo (exported tarball)
        with open(os.path.join(REPO_ROOT, ".gitignore"), encoding="utf-8") as fh:
            assert SITE_TAXDUMP_FILE in fh.read()
    else:
        assert out.returncode == 0, (
            "%s is not gitignored; a machine's taxdump paths would be committed"
            % SITE_TAXDUMP_FILE)
