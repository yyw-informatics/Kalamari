#!/usr/bin/env bash
#
# Create and verify the pinned tool environments for the curation-triage add-on.
#
# Why a script and not a README paragraph.  Tier 2 records a marker call together
# with the software build, the database version and our interpretation rule, so a
# call is only reproducible if the environment is.  Recreating it by hand from a
# package list re-runs the solver and can land on a different build; recreating it
# from the committed explicit lockfiles repeats an exact set of URLs.  This script
# does that, then PROVES the result by running every binary and comparing what it
# reports against src/curation-triage/tool_pins.tsv.
#
# Two things the lockfiles cannot pin, both handled here on purpose:
#
#   * AMRFinderPlus keeps its dated database OUTSIDE the conda package.  The pinned
#     version is 2026-05-15.1 and it must live inside the pinned prefix, because a
#     database found through an ambient CONDA_PREFIX belongs to some other install
#     and would silently stamp the wrong version into the ledger.  If the database
#     is missing this script runs `amrfinder -u`, which downloads the LATEST
#     database -- that may not be the pinned one.  Getting a newer database is an
#     update-gate event, not a silent pass, so the script then compares the
#     installed version against the pin and stops unless --allow-db-drift is given.
#
#     The version string is only a LABEL, though: any directory can be named
#     2026-05-15.1, and no command anywhere asks NCBI for a dated version, so this
#     database is FROZEN rather than re-downloadable.  The evidence is
#     src/curation-triage/amrfinderplus_db_manifest.tsv -- a sha256 for each of the
#     162 files plus one roll-up digest over all of them -- and this script compares
#     the installed database against it FILE BY FILE, naming the files that differ.
#     Measured on the machine the database was frozen on: 239.4 MB verified in 1.6 s
#     with the page cache dropped and 0.8 s warm.  That is cheap enough that there is
#     no separate "quick" mode to get wrong, and it runs once per setup or verify,
#     never once per genome.
#
#     Freezing proves SAMENESS, not availability.  If the local copy is lost and NCBI
#     has moved on, this exact database cannot be fetched again and the manifest can
#     only say that whatever turns up is not it.  Keeping an archive of the dated
#     directory somewhere durable is part of the pin, not an optional extra.
#
#   * ECTyper downloads a 943 MiB MASH species sketch from Zenodo on first use.
#     It is deliberately NOT provisioned.  ECTyper is context-only: decisive task
#     planning drops it, so the pinned pipeline never reaches the code path that
#     needs the sketch.  Downloading a gigabyte for output we never treat as a
#     verdict would be pure cost.
#
# Usage:
#   setup_envs.sh --env-dir <dir> [--only <env>[,<env>...]] [--conda <exe>]
#                 [--verify-only] [--allow-db-drift]
#
#   --env-dir DIR       where the environments live: DIR/<env> per environment
#   --only A,B          act on these environments only (default: all in tool_pins.tsv)
#   --conda EXE         conda/mamba executable (default: $CONDA_EXE, else `conda`)
#   --verify-only       check what is installed; create nothing, download nothing
#   --allow-db-drift    do not fail when the AMRFinderPlus database is not the pin,
#                       by version string or by digest
#
# Exit status is 0 only when every selected environment exists, every binary runs,
# every reported version matches the pin, and the AMRFinderPlus database matches its
# committed digest manifest.

set -euo pipefail

# Resolve the repo root from this script's own location, so the script works from
# any working directory and no absolute path is ever baked into the repo.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
PINS="${REPO_ROOT}/src/curation-triage/tool_pins.tsv"

# The per-file digests of the frozen AMRFinderPlus database.  240 MB cannot be
# committed; 162 digest lines can, and they are what turns "the directory says
# 2026-05-15.1" into "this is byte for byte the database the ledger was built with".
DB_MANIFEST_REL="src/curation-triage/amrfinderplus_db_manifest.tsv"
DB_MANIFEST="${REPO_ROOT}/${DB_MANIFEST_REL}"

# How many differing file names a failure prints before it summarises the rest.
# "digest mismatch" is not actionable; a few names plus a count is -- five names say
# whether this is one corrupted file or a wholesale swap of the database.
DB_DIFF_SHOWN=5

ENV_DIR=""
ONLY=""
CONDA_EXE_ARG="${CONDA_EXE:-conda}"
VERIFY_ONLY="no"
ALLOW_DB_DRIFT="no"

die() {
    printf 'setup_envs.sh: %s\n' "$*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
setup_envs.sh -- create and verify the pinned curation-triage tool environments.

  setup_envs.sh --env-dir <dir> [--only <env>[,<env>...]] [--conda <exe>]
                [--verify-only] [--allow-db-drift]

  --env-dir DIR       where the environments live: DIR/<env> per environment
  --only A,B          act on these environments only (default: all of them)
  --conda EXE         conda/mamba executable (default: $CONDA_EXE, else `conda`)
  --verify-only       check what is installed; create nothing, download nothing
  --allow-db-drift    do not fail when the AMRFinderPlus database is not the pin
                      (neither its version string nor its digest)

Environments and expected versions come from src/curation-triage/tool_pins.tsv;
each one is built from its explicit lockfile in src/curation-triage/envs/.
The AMRFinderPlus database is frozen, not re-downloadable, so it is checked file
by file against src/curation-triage/amrfinderplus_db_manifest.tsv.
Exit status is 0 only when every selected environment exists, every binary runs,
every reported version matches the pin, and that database matches its manifest.
EOF
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        --env-dir)        [ "$#" -ge 2 ] || die "--env-dir needs a directory"
                          ENV_DIR="$2"; shift 2 ;;
        --only)           [ "$#" -ge 2 ] || die "--only needs an environment list"
                          ONLY="$2"; shift 2 ;;
        --conda)          [ "$#" -ge 2 ] || die "--conda needs an executable"
                          CONDA_EXE_ARG="$2"; shift 2 ;;
        --verify-only)    VERIFY_ONLY="yes"; shift ;;
        --allow-db-drift) ALLOW_DB_DRIFT="yes"; shift ;;
        -h|--help)        usage; exit 0 ;;
        *)                die "unknown option '$1' (try --help)" ;;
    esac
done

[ -n "${ENV_DIR}" ] || die "--env-dir is required (try --help)"
[ -f "${PINS}" ] || die "pin table not found at ${PINS}"

# The environments are big -- ~13 GB for all eleven -- and are working state, not
# repository content, so --env-dir normally points OUTSIDE the repo.
mkdir -p -- "${ENV_DIR}"
ENV_DIR="$(cd -- "${ENV_DIR}" && pwd)"

# Is this environment one the caller asked for?
selected() {
    local name="$1"
    [ -z "${ONLY}" ] && return 0
    local wanted
    # A plain comma split; the names are simple identifiers, so this is enough.
    IFS=',' read -r -a wanted <<< "${ONLY}"
    local w
    for w in "${wanted[@]}"; do
        [ "${w}" = "${name}" ] && return 0
    done
    return 1
}

# The first version-looking token of a tool's output.  Taking the LAST word instead
# is a real trap: ectyper prints "ectyper 2.0.0 running database version 1.0", where
# the last token is the database version, so recording it as the software version
# would pin the wrong thing.
#
# The trailing `|| true` is load-bearing, not decoration.  A binary that is present
# but broken the ordinary conda way -- "error while loading shared libraries:
# libgomp.so.1" -- prints no version-looking token, so grep exits 1; under
# `set -euo pipefail` that status travels out of the pipeline and kills the whole
# run at the caller's plain assignment.  The operator then gets a red job that names
# no tool: no FAIL line, no summary, and every later environment left unchecked.
# An unreadable version is a finding to REPORT, so return empty and let the caller's
# note_failure say which tool it was and keep going.
first_version_token() {
    grep -o -E '[0-9]+\.[0-9]+(\.[0-9]+)*' <<< "$1" | head -n 1 || true
}

FAILURES=()
CHECKED=0
CREATED=0
# Set only when the frozen database was actually verified this run.  Without it the
# closing summary would claim a database check that `--only skani` never performed.
DB_DIGEST_VERIFIED="no"

note_failure() {
    FAILURES+=("$1")
    printf '  FAIL  %s\n' "$1"
}

# --------------------------------------------------------------------------- #
#  AMRFinderPlus database                                                      #
# --------------------------------------------------------------------------- #

# One `# name: value` field from the manifest header.
db_manifest_field() {
    sed -n "s/^# $1: *//p" -- "${DB_MANIFEST}" | head -n 1
}

# Every file of the installed database, as `<sha256>  <path>` lines.
#
# The three details here ARE the roll-up's definition, spelled out in the manifest
# header: paths relative to the database directory, sorted by raw bytes (LC_ALL=C,
# not by a locale's collation, which reorders case and punctuation), and the exact
# byte layout `sha256sum` prints.  Change any one of them and every digest below
# still matches while the roll-up does not.
db_file_digests() {
    local db_dir="$1"
    ( cd -- "${db_dir}" && find . -type f -printf '%P\n' | LC_ALL=C sort \
        | xargs -d '\n' -r sha256sum )
}

# Is the installed database the frozen one?
#
# This is the check the version string cannot make.  `amrfinder --database_version`
# reports whatever the directory calls itself, and `amrfinder -u` only ever fetches
# the LATEST database, so nothing can re-download 2026-05-15.1 to compare against.
# The committed manifest is the only evidence there is.
#
# Cost, measured where the database was frozen: 162 files, 239.4 MB, 1.6 s with the
# page cache dropped and 0.8 s warm.  Too cheap to be worth a second, weaker mode --
# a "quick" path that skips reads would pass a corrupted file, and this runs once per
# setup, not once per genome.
verify_amrfinder_db_digest() {
    local db_dir="$1" recorded="$2"
    local expected="${recorded#sha256:}"

    if [ "${expected}" = "${recorded}" ] || [ -z "${expected}" ]; then
        note_failure "amrfinderplus: tool_pins.tsv records database_digest '${recorded:-<nothing>}', which is not a sha256:<hex> roll-up; without it nothing proves the installed database is the frozen one"
        return 0
    fi
    if [ ! -f "${DB_MANIFEST}" ]; then
        note_failure "amrfinderplus: ${DB_MANIFEST_REL} is missing, so nothing can say whether the installed database is the frozen 2026-05-15.1"
        return 0
    fi
    # The pin table and the manifest are two recordings of one fact; if they disagree
    # there is no way to tell which one the database was checked against.
    local manifest_rollup
    manifest_rollup="$(db_manifest_field rollup_sha256)"
    if [ "${manifest_rollup}" != "${expected}" ]; then
        note_failure "amrfinderplus: tool_pins.tsv records database_digest sha256:${expected} but ${DB_MANIFEST_REL} declares rollup_sha256 ${manifest_rollup:-<nothing>}; regenerate the manifest or fix the pin"
        return 0
    fi
    if ! command -v sha256sum > /dev/null 2>&1; then
        printf '  db    digest NOT checked: sha256sum is not on PATH\n'
        return 0
    fi

    local listing actual count
    listing="$(mktemp)"
    # `|| true` for the same reason as everywhere else here: an unreadable database
    # must become a named failure, not an errexit death with no output.  An empty
    # listing digests to something that cannot equal the pin, so it still fails.
    db_file_digests "${db_dir}" > "${listing}" 2>/dev/null || true
    actual="$(sha256sum < "${listing}" | cut -d' ' -f1 || true)"
    count="$(wc -l < "${listing}" | tr -d ' ')"

    if [ "${actual}" = "${expected}" ]; then
        printf '  db    digest %s -- %s file(s) match %s\n' \
            "${expected:0:12}" "${count}" "${DB_MANIFEST_REL}"
        DB_DIGEST_VERIFIED="yes"
        rm -f -- "${listing}"
        return 0
    fi

    # The roll-up only says "not the same database".  Which files differ is what tells
    # a maintainer whether they are looking at one corrupted file, a partial download,
    # or a different release entirely -- so name them.
    local summary
    summary="$(awk -v shown="${DB_DIFF_SHOWN}" '
        NR == FNR {
            split($0, field, "\t")
            want[field[1]] = field[3]
            wanted[++nwant] = field[1]
            next
        }
        {
            # `sha256sum` writes "<hash><two spaces><path>", and a path may hold
            # spaces, so take it as the remainder of the line rather than as $2.
            path = substr($0, length($1) + 3)
            got[path] = $1
            found[++ngot] = path
        }
        END {
            n = 0
            for (i = 1; i <= nwant; i++) {
                p = wanted[i]
                if (!(p in got))            { if (++n <= shown) named[n] = "missing " p }
                else if (got[p] != want[p]) { if (++n <= shown) named[n] = "changed " p }
            }
            for (i = 1; i <= ngot; i++) {
                p = found[i]
                if (!(p in want))           { if (++n <= shown) named[n] = "extra " p }
            }
            list = ""
            for (i = 1; i <= n && i <= shown; i++)
                list = list (list == "" ? "" : ", ") named[i]
            if (n > shown) list = list ", and " (n - shown) " more"
            printf "%d of %d recorded file(s) differ (%d file(s) on disk): %s",
                   n, nwant, ngot, list
        }' <(grep -v '^#' -- "${DB_MANIFEST}" | tail -n +2) "${listing}")"
    rm -f -- "${listing}"

    if [ "${ALLOW_DB_DRIFT}" = "yes" ]; then
        printf '  db    digest %s -- NOT the frozen %s; %s; allowed by --allow-db-drift\n' \
            "${actual:0:12}" "${expected:0:12}" "${summary}"
    else
        note_failure "amrfinderplus database does not match ${DB_MANIFEST_REL}: ${summary}. Its roll-up is ${actual:-<unreadable>}, the frozen one is ${expected}. 2026-05-15.1 cannot be re-downloaded (\`amrfinder -u\` only fetches the LATEST), so restore it from the archive; a deliberate move to another database goes through the update gate, or rerun with --allow-db-drift"
    fi
}

provision_amrfinder_db() {
    local prefix="$1" pinned="$2" binary="$3" digest="$4"
    local data_dir="${prefix}/share/amrfinderplus/data"

    if [ ! -e "${data_dir}/latest" ]; then
        if [ -d "${data_dir}/${pinned}" ]; then
            # The pinned database is already unpacked; just make `latest` point at
            # it, with a RELATIVE link so the prefix stays movable.
            ln -sfn "${pinned}" "${data_dir}/latest"
            printf '  db    linked latest -> %s\n' "${pinned}"
        elif [ "${VERIFY_ONLY}" = "yes" ]; then
            note_failure "amrfinderplus: no database at ${data_dir}/latest; rerun without --verify-only"
            return 0
        else
            printf '  db    downloading with `amrfinder -u` -- this fetches the LATEST database,\n'
            printf '        which may not be the pinned %s; the check below decides.\n' "${pinned}"
            mkdir -p -- "${data_dir}"
            CONDA_PREFIX="${prefix}" AMRFINDER_DB="${data_dir}/latest" \
                PATH="${prefix}/bin:${PATH}" "${prefix}/bin/${binary}" -u \
                || { note_failure "amrfinderplus: 'amrfinder -u' failed"; return 0; }
        fi
    fi

    # Ask the tool, not the filesystem: the whole point is that the version the
    # ledger records is the version the tool actually reads.
    local out installed
    out="$(CONDA_PREFIX="${prefix}" AMRFINDER_DB="${data_dir}/latest" \
           PATH="${prefix}/bin:${PATH}" "${prefix}/bin/${binary}" --database_version 2>&1)" \
        || { note_failure "amrfinderplus: --database_version failed"; return 0; }
    # Same pipeline hazard as first_version_token: if the tool succeeds but prints no
    # "database version" line, grep exits 1 and errexit would abort the run silently.
    # An unreadable database version is a mismatch to report, not a reason to stop.
    installed="$(grep -i 'database version' <<< "${out}" | tail -n 1 | sed 's/.*: *//' || true)"

    if [ "${installed}" = "${pinned}" ]; then
        printf '  db    %s (pinned)\n' "${installed}"
    elif [ "${ALLOW_DB_DRIFT}" = "yes" ]; then
        printf '  db    %s -- NOT the pinned %s, allowed by --allow-db-drift\n' \
            "${installed}" "${pinned}"
    else
        note_failure "amrfinderplus database is '${installed:-<nothing>}' but the pin is '${pinned}'; a database change moves marker calls, so update the pin through the update gate or rerun with --allow-db-drift"
    fi

    # The version above is what the database calls itself; the digest below is what it
    # IS.  Both are reported, and the digest runs even when the version already
    # disagreed -- "162 of 162 files differ" is the sentence that tells the operator
    # they are looking at a different database rather than a mislabelled one.
    verify_amrfinder_db_digest "${data_dir}/latest" "${digest}"
}

# --------------------------------------------------------------------------- #
#  One environment                                                             #
# --------------------------------------------------------------------------- #
handle_env() {
    local tool="$1" env_name="$2" build="$3" binary="$4" version_command="$5"
    local expected="$6" lockfile="$7" sha="$8" db_kind="$9" db_version="${10}"
    local db_digest="${11}"
    local prefix="${ENV_DIR}/${env_name}"
    local lock_path="${REPO_ROOT}/${lockfile}"

    printf '%s (%s, db %s)\n' "${tool}" "${build}" "${db_kind}"
    CHECKED=$((CHECKED + 1))

    if [ ! -f "${lock_path}" ]; then
        note_failure "${tool}: lockfile ${lockfile} is missing"
        return 0
    fi

    # A lockfile that no longer matches its recorded digest means the pin table and
    # the pinned artifact disagree; installing from it would build an environment
    # nobody signed off on.
    if command -v sha256sum > /dev/null 2>&1; then
        local actual
        # `|| true` for the same reason as first_version_token: an unreadable
        # lockfile must land in note_failure with the tool's name, not abort the run
        # before anything is printed.  An empty digest can never equal the pin.
        actual="$(sha256sum -- "${lock_path}" | cut -d' ' -f1 || true)"
        if [ "${actual}" != "${sha}" ]; then
            note_failure "${tool}: ${lockfile} sha256 is ${actual:-<unreadable>} but tool_pins.tsv records ${sha}"
            return 0
        fi
    fi

    if [ ! -d "${prefix}/conda-meta" ]; then
        if [ "${VERIFY_ONLY}" = "yes" ]; then
            note_failure "${tool}: no environment at ${env_name} under the given --env-dir"
            return 0
        fi
        printf '  create from %s\n' "${lockfile}"
        "${CONDA_EXE_ARG}" create -y -p "${prefix}" --file "${lock_path}" > /dev/null \
            || { note_failure "${tool}: conda create failed from ${lockfile}"; return 0; }
        CREATED=$((CREATED + 1))
    else
        printf '  present, not recreated\n'
    fi

    if [ ! -x "${prefix}/bin/${binary}" ]; then
        note_failure "${tool}: ${env_name}/bin/${binary} is missing or not executable"
        return 0
    fi

    if [ "${db_kind}" = "post-install" ]; then
        provision_amrfinder_db "${prefix}" "${db_version}" "${binary}" "${db_digest}"
    fi

    if [ "${version_command}" = "none" ]; then
        # ShigEiFinder exits 2 on --version and prints usage, so it can never confirm
        # its own version.  Saying so is honest; guessing a version would not be.
        printf '  ver   not checkable: %s reports no version, pin is the conda build %s\n' \
            "${binary}" "${expected}"
        return 0
    fi

    # Call through the prefix bin on PATH, not just by absolute path.  abricate and
    # mlst are Perl scripts that need the prefix perl, and LisSero shells out to
    # makeblastdb; called by absolute path alone they die with errors that look like
    # a broken install and are really a missing PATH entry.
    local stdout_txt stderr_txt reported
    local -a version_argv
    # version_command is written binary-first for readability; drop that first word
    # and pass the rest as arguments.
    read -r -a version_argv <<< "${version_command}"
    stderr_txt="$(mktemp)"
    stdout_txt="$(PATH="${prefix}/bin:${PATH}" "${prefix}/bin/${binary}" \
                  "${version_argv[@]:1}" 2> "${stderr_txt}" || true)"
    if [ -z "${stdout_txt//[[:space:]]/}" ]; then
        stdout_txt="$(cat -- "${stderr_txt}")"
    fi
    rm -f -- "${stderr_txt}"

    reported="$(first_version_token "${stdout_txt}")"
    if [ "${reported}" = "${expected}" ]; then
        printf '  ver   %s\n' "${reported}"
    else
        note_failure "${tool}: ${binary} reports '${reported:-<nothing>}' but the pin is '${expected}'"
    fi
}

# --------------------------------------------------------------------------- #
#  Drive it                                                                    #
# --------------------------------------------------------------------------- #
printf 'pins:     %s\n' "src/curation-triage/tool_pins.tsv"
printf 'env-dir:  %s\n' "${ENV_DIR}"
printf 'mode:     %s\n\n' "$([ "${VERIFY_ONLY}" = "yes" ] && echo "verify only" || echo "create and verify")"

SEEN=0
while IFS=$'\t' read -r tool env_name conda_package conda_build binary version_command \
        expected_version lockfile lockfile_sha256 database_kind database_version \
        database_digest notes; do
    # shellcheck disable=SC2034  # conda_package and notes are documentation columns
    : "${conda_package}" "${notes}"
    [ -n "${tool}" ] || continue
    selected "${env_name}" || continue
    SEEN=$((SEEN + 1))
    handle_env "${tool}" "${env_name}" "${conda_build}" "${binary}" "${version_command}" \
        "${expected_version}" "${lockfile}" "${lockfile_sha256}" "${database_kind}" \
        "${database_version}" "${database_digest}"
done < <(grep -v '^#' -- "${PINS}" | tail -n +2)

if [ "${SEEN}" -eq 0 ]; then
    die "--only '${ONLY}' matched no environment in tool_pins.tsv"
fi

printf '\n%d environment(s) checked, %d created, %d problem(s)\n' \
    "${CHECKED}" "${CREATED}" "${#FAILURES[@]}"

if [ "${#FAILURES[@]}" -gt 0 ]; then
    printf '\nThe pinned environments are NOT usable as they stand:\n'
    for f in "${FAILURES[@]}"; do
        printf '  - %s\n' "${f}"
    done
    exit 1
fi

printf 'All pinned tools resolve and report the version recorded in tool_pins.tsv.\n'
if [ "${DB_DIGEST_VERIFIED}" = "yes" ]; then
    printf "AMRFinderPlus's frozen database matches %s file for file. That proves it is\n" \
        "${DB_MANIFEST_REL}"
    printf 'the SAME database, not that it can be downloaded again: NCBI only ever serves\n'
    printf 'the latest one, so keep an archive of the dated directory.\n'
fi
printf "ECTyper's 943 MiB Zenodo MASH sketch is deliberately not provisioned: ECTyper is\n"
printf 'context-only, so the pinned pipeline never plans it and never needs the sketch.\n'
