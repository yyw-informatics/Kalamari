#!/usr/bin/env bash
#
# verify_db_archive.sh -- prove that a copy of the frozen AMRFinderPlus database
# archive is still the right bytes, and still restores to the pinned database.
#
# WHY THIS EXISTS
# The pinned database 2026-05-15.1 is the one part of the tool chain that cannot be
# fetched again: `amrfinder -u` only ever downloads the LATEST database, and NCBI does
# not serve dated versions.  So the archive is not a convenience copy, it is the only
# route back to the calls already recorded in the ledger -- and an archive nobody
# checks is a backup only in name.
#
# The committed manifest already proves that an *installed* database is the frozen one
# (setup_envs.sh does that on every run).  This script answers the other question, the
# one a backup raises: is the copy sitting in some other storage still good?  It needs
# nothing but coreutils and tar, so it can run wherever the copy lives -- no repo
# checkout of the environments, no conda, no AMRFinderPlus.
#
# WHAT IT CHECKS, cheapest first
#   1. the archive's own sha256 matches archive_sha256 in the manifest;
#   2. it unpacks;
#   3. the unpacked directory has exactly the recorded file count;
#   4. every file's sha256 matches the manifest, and the roll-up over all of them
#      matches rollup_sha256.
# Step 1 alone catches bit-rot and truncation. Steps 2-4 catch the case that matters
# more: an archive that is intact but is not this database.
#
#   ./verify_db_archive.sh <path-to-archive.tar.gz>        # full check
#   ./verify_db_archive.sh --quick <path-to-archive.tar.gz>  # digest only, no unpack
#
# Exit status is 0 only when every check asked for passed, so it can be a cron job:
#   0  verified          1  a check failed          2  could not run the check
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MANIFEST="${REPO_ROOT}/src/curation-triage/amrfinderplus_db_manifest.tsv"

QUICK=0
ARCHIVE=""
while [ $# -gt 0 ]; do
    case "$1" in
        --quick) QUICK=1; shift ;;
        --manifest) MANIFEST="${2:-}"; shift 2 ;;
        -h|--help) sed -n '2,32p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*) echo "verify_db_archive.sh: unknown option $1" >&2; exit 2 ;;
        *) ARCHIVE="$1"; shift ;;
    esac
done

die()  { echo "verify_db_archive.sh: $*" >&2; exit 2; }
fail() { echo "FAILED: $*" >&2; exit 1; }

[ -n "${ARCHIVE}" ] || die "give the path to a copy of the archive (--help for usage)"
[ -f "${ARCHIVE}" ] || die "no such file: ${ARCHIVE}"
[ -f "${MANIFEST}" ] || die "no manifest at ${MANIFEST}"

# `# name: value` header fields, the same convention setup_envs.sh reads.
field() { sed -n "s/^# $1: \(.*\)$/\1/p" "${MANIFEST}" | head -1; }

WANT_ARCHIVE_SHA="$(field archive_sha256)"
WANT_ROLLUP="$(field rollup_sha256)"
WANT_FILES="$(field files)"
DB_VERSION="$(field database_version)"
[ -n "${WANT_ARCHIVE_SHA}" ] || die "the manifest records no archive_sha256"
[ -n "${WANT_ROLLUP}" ] || die "the manifest records no rollup_sha256"

echo "database:  ${DB_VERSION}"
echo "archive:   ${ARCHIVE}"

# ---- 1. the archive's own digest ------------------------------------------------
GOT_ARCHIVE_SHA="$(sha256sum "${ARCHIVE}" | cut -d' ' -f1)"
if [ "${GOT_ARCHIVE_SHA}" != "${WANT_ARCHIVE_SHA}" ]; then
    fail "archive digest differs -- this is not the pinned archive
  expected ${WANT_ARCHIVE_SHA}
  got      ${GOT_ARCHIVE_SHA}"
fi
echo "  [ok] archive sha256 matches"

if [ "${QUICK}" -eq 1 ]; then
    echo "VERIFIED (quick): the archive is byte-identical to the pinned copy."
    exit 0
fi

# ---- 2. it unpacks ---------------------------------------------------------------
# mktemp -d, not a fixed path: this may run on storage shared with other jobs, and a
# 240 MB unpack must never land on top of someone else's directory.
WORK="$(mktemp -d "${TMPDIR:-/tmp}/amrfinder-db-verify.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT
tar -xzf "${ARCHIVE}" -C "${WORK}" || fail "the archive did not unpack"

DB_DIR="${WORK}/${DB_VERSION}"
[ -d "${DB_DIR}" ] || DB_DIR="$(find "${WORK}" -mindepth 1 -maxdepth 2 -type d | head -1)"
[ -d "${DB_DIR}" ] || fail "no database directory inside the archive"
echo "  [ok] unpacked"

# ---- 3. the file count -----------------------------------------------------------
GOT_FILES="$(find "${DB_DIR}" -type f | wc -l | tr -d ' ')"
if [ -n "${WANT_FILES}" ] && [ "${GOT_FILES}" != "${WANT_FILES}" ]; then
    fail "expected ${WANT_FILES} files, found ${GOT_FILES}"
fi
echo "  [ok] ${GOT_FILES} files"

# ---- 4. every file, and the roll-up ---------------------------------------------
# Recomputed exactly as the manifest header documents: paths relative to the dated
# directory, sorted by raw bytes under LC_ALL=C, `sha256sum` output concatenated and
# hashed.  Any other ordering gives a different roll-up for identical files.
GOT_ROLLUP="$(cd "${DB_DIR}" && find . -type f -printf '%P\n' \
    | LC_ALL=C sort | xargs -d '\n' sha256sum | sha256sum | cut -d' ' -f1)"
if [ "${GOT_ROLLUP}" != "${WANT_ROLLUP}" ]; then
    # Name the files that moved: "digest mismatch" is not actionable, a few names are.
    echo "roll-up differs; comparing file by file..." >&2
    grep -v '^#' "${MANIFEST}" | tail -n +2 | while IFS=$'\t' read -r path _size want; do
        [ -n "${path:-}" ] || continue
        if [ ! -f "${DB_DIR}/${path}" ]; then
            echo "  missing: ${path}" >&2
        else
            got="$(sha256sum "${DB_DIR}/${path}" | cut -d' ' -f1)"
            [ "${got}" = "${want}" ] || echo "  differs: ${path}" >&2
        fi
    done | head -5 >&2
    fail "the archive unpacks, but it is not the frozen database
  expected roll-up ${WANT_ROLLUP}
  got              ${GOT_ROLLUP}"
fi

echo "  [ok] roll-up sha256 matches"
echo "VERIFIED: this archive restores to AMRFinderPlus database ${DB_VERSION}."
