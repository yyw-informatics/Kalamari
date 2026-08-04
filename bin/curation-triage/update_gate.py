#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
update_gate.py -- the four-point gate a pin change has to pass.

A Tier-2 verdict is only worth reading if a reader can tell WHY it says what it
says.  When a pinned tool, a database, a threshold or one of our own reading rules
moves, some calls move with it -- and the danger is not the movement, it is a moved
call that nobody noticed and nobody can explain.  This module is the check that
makes that impossible to do quietly.  It implements four points, each as a real
check rather than a promise in a document:

  1. **Immutable digest.**  Every pinned artifact has an identity that can be
     recomputed.  For a conda environment that is the sha256 of its explicit
     lockfile, recomputed here from the file on disk and compared with the value
     recorded in ``tool_pins.tsv``.  One environment is different on purpose:
     AMRFinderPlus keeps its database OUTSIDE the conda package, so no lockfile can
     capture it and its identity is the dated version string ``2026-05-15.1``.
     That class (``database_kind = post-install``) is handled explicitly, and so are
     the other three kinds -- see ``digest_problems``.

  2. **Validation panel.**  The panel is not a synthetic fixture, it is the
     committed Tier-2 evidence replayed offline: the 41 units of ``PANEL_SCOPE`` and
     the Tier-2 leads they produce.  ``--panel`` replays it and writes a CALL TABLE --
     unit, check, verdict, call, flags.  Deliberately no ANI number: skani moves the
     last decimal on almost any version bump, so a table containing 96.83 would report
     a difference on every refresh and teach everyone to ignore it.  What has to stay
     put is the CALL, not the measurement behind it.

     The scope is FROZEN in this file rather than derived from the ledger, because
     curation moves the ledger every month: an SME resolving a Tier-1 lead would
     shrink the panel, move its digest and fail a gate that has nothing to do with
     pins.  A pinned tool must be validated against a fixed body of evidence, and
     ``panel:scope`` is itself a pinned component so that changing that body is a
     recorded, reviewed decision (a ``panel`` class changelog row) instead of a
     side effect.  Point 2 also checks PROVENANCE: the committed caches say which
     tool versions produced them, and those have to be the pinned ones -- a panel
     replayed from the previous pin's caches signs off on nothing.

  3. **Classification.**  ``classify_difference()`` is a pure function that maps each
     changed call to ``software`` / ``database`` / ``threshold`` / ``rule``, decided
     by which pinned component the call depends on actually moved between the two pin
     states.  A fifth class, ``panel``, covers a call that appeared or disappeared
     because the panel's scope moved -- a curation decision, not a tool difference.
     ``unattributable`` -- a call that changed while nothing it depends on moved -- is
     a FAILURE, not a class.  It means the pipeline is not reproducible, and no amount
     of reviewing makes that acceptable.

  4. **Changelog.**  ``src/curation-triage/pin_changelog.tsv`` is the record of every
     pin identity this repo has ever had.  Each component's rows form a chain: the
     first row starts from ``-`` and every later row starts where the previous one
     ended.  A pin that moved without a new row fails the gate, and the failure says
     exactly which row to add.

CLI
---
    # what CI runs: offline, no tools, no network.  Must pass on the committed tree.
    python3 bin/curation-triage/update_gate.py --check

    # the call table for one pin state (write it before AND after a pin change)
    python3 bin/curation-triage/update_gate.py --panel --out before.tsv

    # classify what moved between two call tables
    python3 bin/curation-triage/update_gate.py --compare before.tsv after.tsv

    # the current pin state, one component per line (the changelog's raw material)
    python3 bin/curation-triage/update_gate.py --print-state

The full loop for a pin change is three steps, because the "before" call table is
not a committed file -- it is a property of the tree you are leaving:

    1. on the CURRENT tree:   --panel --out before.tsv
    2. change the pin, then:  --panel --out after.tsv
    3.                        --compare before.tsv after.tsv
       then add the pin_changelog.tsv row it prints, and re-run --check.

``--compare`` runs at step 3, BEFORE the new row exists, so the changelog's last row
per component still holds the identity the tree had at step 1.  That is what
``previous_state`` reconstructs the before-state from: each component's last recorded
``to``.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import os
import re
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import markers as markers_mod  # noqa: E402
import skani_runner  # noqa: E402
import tier2  # noqa: E402
import tier2_confirm  # noqa: E402
import typers as typers_mod  # noqa: E402

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
SRC = os.path.join(REPO_ROOT, "src", "curation-triage")
DEFAULT_PINS = os.path.join(SRC, "tool_pins.tsv")
DEFAULT_CHANGELOG = os.path.join(SRC, "pin_changelog.tsv")
DEFAULT_MANIFEST = os.path.join(SRC, "marker_manifest.tsv")

# The five classes a real difference can have, plus the one that is a failure.
CLASS_SOFTWARE = "software"
CLASS_DATABASE = "database"
CLASS_THRESHOLD = "threshold"
CLASS_RULE = "rule"
# `panel` is not a pin that produced a different answer -- it is a decision about WHICH
# evidence the gate replays.  It exists so that changing that decision has a legal way
# into the record: without it the only class left for a call that came or went with the
# scope would be `unattributable`, which would file a curation outcome as a
# reproducibility failure and leave the gate with no row anyone could add.
CLASS_PANEL = "panel"
CLASS_UNATTRIBUTABLE = "unattributable"
CHANGE_CLASSES = (CLASS_SOFTWARE, CLASS_DATABASE, CLASS_THRESHOLD, CLASS_RULE,
                  CLASS_PANEL)

# When several kinds of component moved at once the class reported is the first of
# these that moved.  Neither half of the order is arbitrary:
#   * `rule` and `threshold` come first because they are OUR OWN edits.  A reviewer
#     has to be shown those before a third party's version bump: they are the changes
#     we can revert and re-test in minutes.
#   * `software` comes before `database` because a BUNDLED database moves only
#     because its conda package moved -- the package is the root cause and the one
#     thing to revert.  A `database` verdict therefore means the data moved on its
#     own, which today only AMRFinderPlus's dated database can do.
# Every moved component is still listed on the difference row, so where two things
# really did move at once the ambiguity is visible rather than hidden.
CLASS_PRECEDENCE = (CLASS_RULE, CLASS_THRESHOLD, CLASS_SOFTWARE, CLASS_DATABASE)

# Component id prefixes.  The prefix carries the class, so a changelog row can be
# classified without consulting any other file.
PREFIX_CLASS = {
    "tool": CLASS_SOFTWARE,
    "db": CLASS_DATABASE,
    "threshold": CLASS_THRESHOLD,
    "rule": CLASS_RULE,
    "panel": CLASS_PANEL,
}

# The tier2 evaluation thresholds live in code, not in a table, so they get one
# synthetic component of their own.  Without it a changed species boundary would be
# an unattributable difference, which is exactly the failure this gate exists for.
ANI_THRESHOLD_COMPONENT = "threshold:tier2_ani"

# The panel's scope is a pinned component too: it decides which evidence the gate
# replays, so moving it moves the digest that every other pin is signed off with.
PANEL_SCOPE_COMPONENT = "panel:scope"

# skani answers no marker question, so it is absent from the typer cache; its version
# is recorded in the ANI cache under its own key.
SKANI_TOOL_KEY = "skani"
# The one database whose version the caches record on its own, because it is the one
# database that can move without its conda package moving.
AMRFINDER_DB_COMPONENT = "db:amrfinderplus"
# What a version probe writes when the tool refuses to say (`shigeifinder --version`
# exits 2 and prints usage).  It is a genuine recorded identity, not a failure.
UNKNOWN_VERSION = "unknown"

DATABASE_KINDS = ("bundled", "post-install", "runtime-download", "none")

# --------------------------------------------------------------------------- #
#  The frozen panel scope                                                      #
# --------------------------------------------------------------------------- #
# WHY THIS LIST IS WRITTEN OUT AND NOT COMPUTED
# The obvious alternative -- whatever `trigger_units` returns for the CURRENT ledger
# and policy table -- is wrong.  Neither of those is a pin: the monthly collate appends
# to the ledger every time an SME resolves a lead, and a resolved lead leaves the
# trigger set.  The panel would then shrink, its digest would move, and `--check` would
# fail -- with no legal row to add, because the pins had not moved and a changelog row
# may not record from == to.  Correct curation would jam the gate, and `--compare`
# would call the vanished call `unattributable`: a reproducibility failure that never
# happened.
#
# So the scope is frozen here, in the gate itself.  The gate asks "do the pinned tools
# still produce the same calls on the same evidence?", and that question only means
# something if the evidence is held still.  Adding or removing a unit is a real
# decision about what the gate validates, and it goes through the changelog like any
# other pinned thing: a `panel:scope` row whose from/to are the scope identities, plus
# the new panel digest in its validation column.
#
# The 21 confounded units carry the regimes Tier 2 exists for (the sub-species panels
# and the marker-decided groups); the 20 below them are the Tier-1 lead units the scope
# was frozen with.  Both halves are spelled out, so nothing outside this file can
# change what the gate replays.
PANEL_CONFOUNDED_UNITS = (
    "Bacillus_cereus~GCF_000283675.1",
    "Bacillus_thuringiensis~GCF_000008505.1",
    "Clostridium_botulinum_groupI~GCF_000020285.1",
    "Clostridium_botulinum_groupII~GCF_000019305.1",
    "Escherichia_coli~GCF_003018575.1",
    "Listeria_monocytogenes_I~GCF_013282705.1",
    "Listeria_monocytogenes_II~GCF_013282685.1",
    "Listeria_monocytogenes_III~GCF_013282645.1",
    "Listeria_monocytogenes_IV~GCF_013282665.1",
    "Salmonella_enterica_I~GCF_000006945.2",
    "Salmonella_enterica_IIa~GCF_006051055.1",
    "Salmonella_enterica_IIb~GCF_900635535.1",
    "Salmonella_enterica_IIIa~GCF_000018625.1",
    "Salmonella_enterica_IIIb~GCA_013162305.1",
    "Salmonella_enterica_IV~GCF_013161805.1",
    "Salmonella_enterica_VI~GCF_006051135.1",
    "Salmonella_enterica_VII~GCF_013162165.1",
    "Salmonella_enterica_VIII~GCA_013138165.1",
    "Salmonella_enterica_IX~GCF_013348865.1",
    "Salmonella_enterica_X~GCF_013162025.1",
    "Yersinia_pseudotuberculosis~GCF_000834295.1",
)

PANEL_LEAD_UNITS = (
    "Buchnera_aphidicola~GCF_000009605.1",
    "Campylobacter_fetus~GCF_000495505.1",
    "Candidatus_Korarchaeum_cryptofilum~GCF_000019605.1",
    "Citrobacter_freundii~GCF_000648515.1",
    "Clavibacter_michiganensis_sepedonicus~GCF_000069225.1",
    "Clostridium_argentinense~GCF_002074155.1",
    "Finegoldia_magna~GCF_000010185.1",
    "Geobacter_sulfurreducens~GCF_000007985.2",
    "Ketogulonicigenium_vulgare~GCF_000223375.1",
    "Leptothrix_cholodnii~GCF_000019785.1",
    "Lysinibacillus_sphaericus~GCF_000017965.1",
    "Mycobacterium_leprae~GCF_000195855.1",
    "Mycoplasma_mycoides~GCF_000011445.1",
    "Polynucleobacter_necessarius~GCF_900096765.1",
    "Prochlorococcus_marinus~GCF_000007925.1",
    "Sinorhizobium_fredii~GCF_000018545.1",
    "Streptococcus_mitis~GCF_000027165.1",
    "Thermoanaerobacter_pseudethanolicus~GCF_000019085.1",
    "Yersinia_frederiksenii~GCF_002591095.1",
    "Yersinia_massiliensis~GCF_013282765.1",
)

PANEL_SCOPE = tuple(sorted(set(PANEL_CONFOUNDED_UNITS) | set(PANEL_LEAD_UNITS)))

CHANGELOG_COLUMNS = ("date", "component", "class", "from", "to",
                     "panel_calls_total", "panel_calls_changed", "validation",
                     "reviewer", "note")

# A changelog row records WHICH panel state was accepted, not merely that a panel
# was run.  The digest is what makes point 2 checkable offline: no baseline call
# table is committed, so the accepted panel lives in this column.
VALIDATION_PREFIX = "panel-replay:"
DIGEST_LENGTH = 12

# `from`/`to` value meaning "there was no pin" -- used by a baseline row (nothing
# before it) and by a removal row (nothing after it).
NO_PIN = "-"

PANEL_COLUMNS = ("unit_id", "check", "key", "verdict", "call", "flags", "components")
# The fields the digest covers.  `components` is a property of the pin state rather
# than of the call, so it stays out: re-pointing a call at a renamed component must
# not read as the call having changed.
DIGEST_COLUMNS = ("unit_id", "check", "key", "verdict", "call", "flags")

PANEL_CMD = "python3 bin/curation-triage/update_gate.py --panel --out <path>"
CHECK_CMD = "python3 bin/curation-triage/update_gate.py --check"
COMPARE_CMD = "python3 bin/curation-triage/update_gate.py --compare <before> <after>"


def log(msg: str) -> None:
    print(msg, file=sys.stderr)


def rel(path: str) -> str:
    """Repo-relative when inside the repo, so a message never prints a local path."""
    full = os.path.abspath(path)
    root = REPO_ROOT + os.sep
    return full[len(root):] if full.startswith(root) else path


def sha256_of(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def _cell(value: str) -> str:
    """One TSV cell: a tab or newline inside a tool's label would split the row."""
    return (value or "").replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()


def read_tsv(path: str, key_column: str) -> List[Dict[str, str]]:
    """Header-keyed TSV rows; leading ``#`` lines are documentation, as everywhere."""
    with open(path, encoding="utf-8", newline="") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return [r for r in csv.DictReader(lines, delimiter="\t")
            if (r.get(key_column) or "").strip()]


# --------------------------------------------------------------------------- #
#  The pin state: every component whose movement can move a call               #
# --------------------------------------------------------------------------- #
def load_pins(path: str = DEFAULT_PINS) -> List[Dict[str, str]]:
    return read_tsv(path, "tool")


def component_class(component: str) -> str:
    """The class of a component id, from its prefix; '' when the prefix is unknown."""
    return PREFIX_CLASS.get(component.split(":", 1)[0], "")


def software_identity(row: Dict[str, str]) -> str:
    """A conda environment's identity: the build plus a short lockfile digest.

    Both are needed.  The build alone misses a dependency bump that keeps the same
    primary version but reinstalls half the prefix; the digest alone is unreadable in
    a changelog.  The digest is taken from ``tool_pins.tsv`` rather than from the file
    on disk deliberately -- point 1 is what checks the table against the disk, and
    keeping the two checks independent means a tampered lockfile produces one clear
    failure instead of two confusing ones.
    """
    return "%s+lock:%s" % (row.get("conda_build", ""),
                           (row.get("lockfile_sha256") or "")[:DIGEST_LENGTH])


def database_identity(row: Dict[str, str]) -> str:
    """A database's identity: its version string plus a short content digest.

    The version string alone is a *label*, and a label is exactly what an unattended
    update can move without telling anyone -- the AMRFinderPlus database is fetched
    outside the conda package, so nothing about the package changes when it does.
    Appending the digest recorded in ``tool_pins.tsv`` makes the identity name the
    bytes, so a database swapped under an unchanged dated string is a gate failure
    rather than only a ``setup_envs.sh`` failure.

    A row with no digest keeps the bare version, which is right for the bundled
    databases: their bytes are already pinned by the lockfile, and point 1 checks that.

    Deliberately used ONLY to build the pin state.  The two comparisons that check the
    pins against other files -- ``marker_manifest.tsv`` in point 1 and the caches'
    ``amrfinder_database_version`` in point 2 -- keep reading the plain version column,
    because those files record what a tool printed and must not be asked to carry our
    digest.
    """
    version = (row.get("database_version") or "").strip()
    digest = (row.get("database_digest") or "").strip()
    if not digest or not digest.startswith("sha256:"):
        return version
    return "%s+db:%s" % (version, digest.split(":", 1)[1][:DIGEST_LENGTH])


def threshold_identity() -> str:
    """The Tier-2 numeric thresholds that decide an ANI or panel call.

    skani's own screen floors are in here too: they decide whether a pair is reported
    at all, so raising them turns a measured comparison into ``screened_out``, which
    is a different call.
    """
    values = {
        "species_ani": tier2.DEFAULT_SPECIES_ANI,
        "af_gate": tier2.DEFAULT_AF_GATE,
        "derep": tier2.DEFAULT_DEREP,
        "epsilon": tier2.DEFAULT_EPSILON,
        "screen": skani_runner.DEFAULT_SCREEN,
        "min_af": skani_runner.DEFAULT_MIN_AF,
    }
    return ";".join("%s=%s" % (k, values[k]) for k in sorted(values))


def panel_scope_identity(scope: Sequence[str] = PANEL_SCOPE) -> str:
    """The frozen panel scope's identity: how many units, and which ones.

    A digest rather than the list itself, because a changelog cell has to stay
    readable; the count is in front of it so a shrinking panel is visible at a glance,
    which is the failure this whole component exists to make loud.
    """
    body = "\n".join(sorted(set(scope)))
    return "%d units+scope:%s" % (
        len(set(scope)),
        hashlib.sha256(body.encode("utf-8")).hexdigest()[:DIGEST_LENGTH])


def task_threshold_identity(rows: Sequence["markers_mod.MarkerRow"]) -> str:
    """One task's per-component marker thresholds, in a stable spelling."""
    parts = []
    for row in sorted(rows, key=lambda r: r.component):
        if not row.thresholds:
            continue
        inner = ",".join("%s=%s" % (k, row.thresholds[k]) for k in sorted(row.thresholds))
        parts.append("%s[%s]" % (row.component, inner))
    return " | ".join(parts) or NO_PIN


def pin_state(pins_path: str = DEFAULT_PINS,
              manifest_path: str = DEFAULT_MANIFEST) -> Dict[str, str]:
    """Every pinned component -> its current identity.

    Four kinds of thing can move a Tier-2 call, so all four are components here:
    the software (one conda environment per row of ``tool_pins.tsv``), the databases
    those environments read, the thresholds, and our own interpretation rules.  A
    difference that no component explains is unattributable, so a component missing
    from this map turns a real cause into a gate failure -- which is why the map is
    built from the committed tables rather than hand-listed.

    The panel's scope is the fifth: it moves no call's content, but it decides which
    calls exist at all, so a call that came or went with it must be explainable too.
    """
    state: Dict[str, str] = {}
    for row in load_pins(pins_path):
        env = (row.get("env") or "").strip()
        state["tool:%s" % env] = software_identity(row)
        if (row.get("database_kind") or "").strip() != "none":
            state["db:%s" % env] = database_identity(row)

    manifest = markers_mod.load_manifest(manifest_path)
    for task_id in sorted({r.task_id for r in manifest}):
        rows = [r for r in manifest if r.task_id == task_id]
        state["rule:%s" % task_id] = rows[0].interpretation_rule_version or NO_PIN
        state["threshold:%s" % task_id] = task_threshold_identity(rows)

    state[ANI_THRESHOLD_COMPONENT] = threshold_identity()
    state[PANEL_SCOPE_COMPONENT] = panel_scope_identity()
    return dict(sorted(state.items()))


def call_components(check: str, task_id: str = "",
                    state: Optional[Dict[str, str]] = None) -> List[str]:
    """The pinned components one call depends on.

    A marker call depends on the environments of its task's tools, on those
    environments' databases, on the task's thresholds and on the task's rule.  An ANI
    or panel call depends on skani and on the Tier-2 thresholds.  Nothing else in the
    tree can change either kind of call without being one of these -- which is the
    assumption ``unattributable`` exists to test.
    """
    state = state if state is not None else {}
    if check == tier2.CHECK_MARKER:
        components = {"rule:%s" % task_id, "threshold:%s" % task_id}
        for tool in typers_mod.TASK_TOOLS.get(task_id, []):
            env = typers_mod.TOOL_ENVS.get(tool, tool)
            components.add("tool:%s" % env)
            if "db:%s" % env in state:
                components.add("db:%s" % env)
        return sorted(components)
    return sorted({"tool:skani", ANI_THRESHOLD_COMPONENT})


# --------------------------------------------------------------------------- #
#  Point 1 -- immutable digest                                                 #
# --------------------------------------------------------------------------- #
def regen_hint(env: str) -> str:
    return ("conda list --explicit --md5 -p <env-dir>/%s > src/curation-triage/"
            "envs/%s.linux-64.lock (keep the header), then update lockfile_sha256 "
            "in src/curation-triage/tool_pins.tsv" % (env, env))


def digest_problems(pins: Sequence[Dict[str, str]],
                    repo_root: str = REPO_ROOT,
                    panel: Optional["Panel"] = None) -> List[str]:
    """Point 1: recompute every pinned artifact's identity and report every mismatch.

    Each ``database_kind`` gets the check its physical situation allows, because
    pretending they are all the same is how an unpinned database slips through:

      ``bundled``          the database ships inside the conda package, so the
                           lockfile digest already pins it; only its recorded
                           identity has to be present.
      ``post-install``     AMRFinderPlus: the database is fetched after install and
                           lives outside the package.  NO lockfile can capture it, so
                           its digest IS the dated version string -- which therefore
                           has to look like a date and has to agree with the second
                           place it is recorded, ``marker_manifest.tsv``.
      ``runtime-download`` the tool fetches its data when it runs, so nothing pins
                           it.  A call that depends on such a tool is not
                           reproducible, so the check is that no panel call does.
      ``none``             no database; the lockfile is the whole identity.
    """
    problems: List[str] = []
    manifest_db = {}
    manifest_path = os.path.join(repo_root, "src", "curation-triage", "marker_manifest.tsv")
    if os.path.exists(manifest_path):
        for row in markers_mod.load_manifest(manifest_path):
            if row.conda_package and row.database_version:
                manifest_db.setdefault(row.conda_package, set()).add(row.database_version)

    for row in pins:
        tool = (row.get("tool") or "").strip()
        env = (row.get("env") or "").strip()
        kind = (row.get("database_kind") or "").strip()
        lockfile = (row.get("lockfile") or "").strip()
        recorded = (row.get("lockfile_sha256") or "").strip()
        path = os.path.join(repo_root, lockfile)

        if kind not in DATABASE_KINDS:
            problems.append(
                "%s: database_kind %r is not one of %s. Fix the row in "
                "src/curation-triage/tool_pins.tsv."
                % (tool, kind, "/".join(DATABASE_KINDS)))

        if not os.path.exists(path):
            problems.append(
                "%s: %s names lockfile %s, which does not exist. Restore it or rerun: %s"
                % (tool, rel(DEFAULT_PINS), lockfile, regen_hint(env)))
        else:
            actual = sha256_of(path)
            if actual != recorded:
                problems.append(
                    "%s: %s has sha256 %s on disk but src/curation-triage/tool_pins.tsv "
                    "records %s. Either the lockfile was edited by hand (restore it) or "
                    "the environment was rebuilt -- in that case rerun: %s, then add a "
                    "src/curation-triage/pin_changelog.tsv row for tool:%s."
                    % (tool, lockfile, actual, recorded or "(nothing)",
                       regen_hint(env), env))
            # The recorded build is half of this component's changelog identity, so it
            # has to be the build the lockfile actually installs.  A build edited in
            # the table alone would put an identity in the changelog that no artifact
            # can reproduce -- the opposite of an immutable digest.
            with open(path, encoding="utf-8") as fh:
                text = fh.read()
            stamp = "/%s-%s." % ((row.get("conda_package") or "").strip(),
                                 (row.get("conda_build") or "").strip())
            if stamp not in text:
                problems.append(
                    "%s: src/curation-triage/tool_pins.tsv records build %s but %s "
                    "installs no such package. Record the build the lockfile installs, "
                    "or rerun: %s"
                    % (tool, stamp.strip("/."), lockfile, regen_hint(env)))

        database_version = (row.get("database_version") or "").strip()
        if kind == "post-install":
            # The one identity in the whole pin table that no committed file can
            # reproduce.  It is checked against its shape and against its second
            # recording instead, so a silent database bump still has to show up here.
            if not re.match(r"^\d{4}-\d{2}-\d{2}", database_version):
                problems.append(
                    "%s: database_kind is post-install, so its digest is the dated "
                    "database version, but database_version is %r. Record the exact "
                    "string `%s --database_version` prints in "
                    "src/curation-triage/tool_pins.tsv."
                    % (tool, database_version, row.get("binary") or tool))
            seen = manifest_db.get((row.get("conda_package") or "").strip(), set())
            if seen and database_version not in seen:
                problems.append(
                    "%s: src/curation-triage/tool_pins.tsv records database %s but "
                    "src/curation-triage/marker_manifest.tsv records %s. A dated "
                    "database that is pinned in one table and not the other makes "
                    "every call it produced unattributable; update both."
                    % (tool, database_version or "(nothing)", "/".join(sorted(seen))))
        elif kind in ("bundled", "runtime-download"):
            if not database_version:
                problems.append(
                    "%s: database_kind is %s but database_version is empty. Record what "
                    "the installed database calls itself in "
                    "src/curation-triage/tool_pins.tsv." % (tool, kind))

    if panel is not None:
        unpinnable = {"tool:%s" % (r.get("env") or "").strip()
                      for r in pins
                      if (r.get("database_kind") or "").strip() == "runtime-download"}
        for call in panel.calls:
            for component in call.components:
                if component in unpinnable:
                    problems.append(
                        "%s depends on %s, whose data is fetched at run time and is "
                        "pinned by nothing. Either provision and pin that data, or keep "
                        "the tool context-only so it decides no call."
                        % (call.key, component))
    return problems


# --------------------------------------------------------------------------- #
#  Point 2 -- the validation panel                                             #
# --------------------------------------------------------------------------- #
@dataclass
class Call:
    """One Tier-2 call, in the only form that is stable across a tool bump."""

    unit_id: str = ""
    check: str = ""
    key: str = ""            # the Tier-2 lead key: unique, and it names the task
    verdict: str = ""
    call: str = ""           # the label the tool produced ('' where a check has none)
    flags: str = ""          # the pipeline's own classification of the verdict
    components: List[str] = field(default_factory=list)

    def digest_row(self) -> str:
        return "\t".join(_cell(getattr(self, c)) for c in DIGEST_COLUMNS)


@dataclass
class Panel:
    """A replayed validation panel: its calls, how completely it ran, and who made it.

    ``missing`` names scope units the committed inputs no longer supply, and
    ``provenance`` is what the caches say about the tools that produced them.  Both
    are carried on the panel rather than checked inside the replay, so the replay
    stays a plain "run it and report", and every verdict lives in one of the four
    numbered points.
    """

    calls: List[Call] = field(default_factory=list)
    units: List[str] = field(default_factory=list)
    unevaluated: List[str] = field(default_factory=list)
    missing: List[str] = field(default_factory=list)
    provenance: Dict[str, str] = field(default_factory=dict)

    @property
    def digest(self) -> str:
        return panel_digest(self.calls)


def panel_digest(calls: Sequence[Call]) -> str:
    """A short, stable digest of a whole call table.

    Short on purpose: it goes in a changelog column a human reads, and 12 hex digits
    is far beyond any chance of two different panels colliding by accident.
    """
    body = "\n".join(sorted(c.digest_row() for c in calls))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()[:DIGEST_LENGTH]


def task_of(key: str) -> str:
    """The marker task inside a Tier-2 lead key (``t2:<unit>#marker:<task>``)."""
    suffix = key.split("#", 1)[1] if "#" in key else ""
    return suffix.split("marker:", 1)[1] if suffix.startswith("marker:") else ""


def cache_provenance(ani_meta: Dict[str, object],
                     typer_meta: Dict[str, object]) -> Dict[str, str]:
    """What the committed caches say about the tools that produced them.

    Keyed the way the caches key it -- by ``typers.py`` tool key, not by conda
    environment -- plus the one dated database that lives outside its package.  The
    two spellings are reconciled in ``provenance_problems``; keeping the raw strings
    here means a mismatch is reported as the two files actually spell it.
    """
    out: Dict[str, str] = {}
    skani = str(ani_meta.get("skani_version") or "").strip()
    if skani:
        out[SKANI_TOOL_KEY] = skani
    for tool, version in (typer_meta.get("tool_versions") or {}).items():
        text = str(version or "").strip()
        if text:
            out[str(tool)] = text
    database = str(typer_meta.get("amrfinder_database_version") or "").strip()
    if database:
        out[AMRFINDER_DB_COMPONENT] = database
    return out


def replay_panel(src: str = SRC, chromosomes: Optional[str] = None,
                 cache_dir: Optional[str] = None,
                 state: Optional[Dict[str, str]] = None,
                 scope: Optional[Sequence[str]] = None) -> Panel:
    """Point 2: replay the committed Tier-2 evidence offline and collect the calls.

    This calls the SAME functions the monthly run calls -- ``build_plan``,
    ``evaluate`` -- so the panel cannot drift away from the thing it is validating.
    It reads only committed caches: no tool, no network.

    The one deliberate difference from the monthly run is WHICH units it runs over.
    The monthly run asks the ledger ("what is still open?"); the panel does not, and
    must not: the ledger changes every time an SME resolves a lead, and a gate whose
    evidence moves under it validates nothing.  The scope is ``PANEL_SCOPE``, frozen
    and pinned, and a scope unit the policy table no longer supplies is reported in
    ``missing`` rather than quietly dropped.  For the same reason the replay passes no
    survivor set to ``evaluate``: that argument only fills in the "confirms this
    Tier-1 lead" cross-link, which is a property of the ledger, not of a call.
    """
    chromosomes = chromosomes or os.path.join(REPO_ROOT, "src", "chromosomes.tsv")
    cache_dir = cache_dir or src
    policy = os.path.join(src, "taxon_policy.tsv")
    wanted = list(PANEL_SCOPE if scope is None else scope)

    import coverage  # noqa: E402  (imported here: only the panel needs it)
    import genomes as genomes_mod  # noqa: E402
    import policy as policy_mod  # noqa: E402

    source_lanes = coverage.load_source_lanes(chromosomes)
    picks = coverage.load_picks(policy, source_lanes)
    unit_meta = tier2_confirm.load_unit_meta(policy)
    confounded = policy_mod.load_confounded_manifest(
        os.path.join(src, "confounded_manifest.tsv"))
    manifest = markers_mod.load_manifest(os.path.join(src, "marker_manifest.tsv"))
    by_unit = {pick.unit_id: pick for pick in picks}
    triggers = [(by_unit[unit], ["frozen panel scope"], unit_meta.get(unit, {}))
                for unit in sorted(set(wanted)) if unit in by_unit]
    missing = sorted(set(wanted) - set(by_unit))

    ts_cache = genomes_mod.TypeStrainCache(
        cache_path=os.path.join(cache_dir, tier2_confirm.TYPE_STRAIN_CACHE_NAME))
    plan = tier2_confirm.build_plan(triggers, manifest, confounded, ts_cache,
                                    os.path.join(src, "panels"))
    ani_cache = skani_runner.AniCache(
        cache_path=os.path.join(cache_dir, tier2_confirm.ANI_CACHE_NAME))
    typer_cache = typers_mod.TyperCache(
        cache_path=os.path.join(cache_dir, tier2_confirm.TYPER_CACHE_NAME))

    # The thresholds are the pinned defaults, never a caller's override: a panel run
    # with different numbers would compare two different questions.
    args = argparse.Namespace(species_ani=tier2.DEFAULT_SPECIES_ANI,
                              af_gate=tier2.DEFAULT_AF_GATE,
                              derep=tier2.DEFAULT_DEREP,
                              epsilon=tier2.DEFAULT_EPSILON,
                              emit_unevaluated=False)
    leads, _resolutions, unevaluated = tier2_confirm.evaluate(
        triggers, plan, ani_cache, typer_cache, manifest, args)

    state = state if state is not None else pin_state(
        manifest_path=os.path.join(src, "marker_manifest.tsv"))
    calls = []
    for lead in leads:
        evidence = lead.evidence or {}
        check = str(evidence.get("t2_check") or "")
        calls.append(Call(
            unit_id=lead.unit_id,
            check=check,
            key=lead.lead_key,
            verdict=str(evidence.get("verdict") or ""),
            call=str(evidence.get("call") or ""),
            flags=",".join(lead.actionable_flags or []),
            components=call_components(check, task_of(lead.lead_key), state),
        ))
    calls.sort(key=lambda c: c.key)
    return Panel(calls=calls,
                 units=sorted({p.unit_id for p, _r, _m in triggers}),
                 unevaluated=sorted(unevaluated),
                 missing=missing,
                 provenance=cache_provenance(ani_cache.meta, typer_cache.meta))


def write_panel(panel: Panel, path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("# Tier-2 validation panel: %d call(s) over %d unit(s), "
                 "%s%s\n" % (len(panel.calls), len(panel.units),
                             VALIDATION_PREFIX, panel.digest))
        fh.write("# Written by %s. No ANI number is recorded here on purpose: skani "
                 "moves the last decimal on almost any bump, and a table that reports "
                 "a difference every time teaches everyone to ignore it.\n" % PANEL_CMD)
        fh.write("\t".join(PANEL_COLUMNS) + "\n")
        for call in panel.calls:
            fh.write("\t".join([
                _cell(call.unit_id), _cell(call.check), _cell(call.key),
                _cell(call.verdict), _cell(call.call), _cell(call.flags),
                ",".join(call.components)]) + "\n")


def read_panel(path: str) -> List[Call]:
    """Read a call table written by ``--panel`` (its digest is reproducible)."""
    calls = []
    for row in read_tsv(path, "key"):
        calls.append(Call(
            unit_id=(row.get("unit_id") or "").strip(),
            check=(row.get("check") or "").strip(),
            key=(row.get("key") or "").strip(),
            verdict=(row.get("verdict") or "").strip(),
            call=(row.get("call") or "").strip(),
            flags=(row.get("flags") or "").strip(),
            components=[c for c in (row.get("components") or "").split(",") if c],
        ))
    calls.sort(key=lambda c: c.key)
    return calls


# --------------------------------------------------------------------------- #
#  Point 3 -- classification                                                   #
# --------------------------------------------------------------------------- #
def moved_components(components: Sequence[str], before: Dict[str, str],
                     after: Dict[str, str]) -> List[str]:
    """Which of a call's components have a different identity in the two states."""
    return sorted(c for c in set(components) if before.get(c) != after.get(c))


def classify_difference(before: Optional[Call], after: Optional[Call],
                        before_state: Dict[str, str],
                        after_state: Dict[str, str]) -> str:
    """Point 3 (pure): why did this call change?

    Either side may be ``None`` -- a call that appeared or disappeared is a change
    too, and it is classified the same way, by what moved underneath it.

    ``unattributable`` is returned when nothing the call depends on moved.  That is
    not a mild finding to be recorded and waved through: it means the same pins
    produced two different answers, so the pipeline is not reproducible and the gate
    must fail.  Chasing it usually ends at a cache that was refreshed without its
    tool version changing, or at an input file edited outside the pin tables.
    """
    if ((before is None) != (after is None)
            and before_state.get(PANEL_SCOPE_COMPONENT)
            != after_state.get(PANEL_SCOPE_COMPONENT)):
        # A unit that entered or left the frozen scope fully explains a call that
        # appeared or disappeared, and no pin has to have moved for it.  It is decided
        # before the pin classes for the same reason `rule` outranks a vendor bump: it
        # is OUR OWN edit, and it is the one a reviewer can revert in seconds.
        # A call present on BOTH sides is deliberately excluded: the scope decides
        # which units are replayed, never what a call says, so a verdict that moved
        # during a scope change is still unattributable and still a failure.
        return CLASS_PANEL
    components = set()
    for call in (before, after):
        if call is not None:
            components |= set(call.components)
    classes = {component_class(c) for c in moved_components(components, before_state,
                                                            after_state)}
    for class_name in CLASS_PRECEDENCE:
        if class_name in classes:
            return class_name
    return CLASS_UNATTRIBUTABLE


def describe(call: Optional[Call]) -> str:
    """One side of a difference, compact enough to read in a terminal."""
    if call is None:
        return "(absent)"
    return "%s/%s/%s" % (call.verdict, call.call or "-", call.flags or "-")


@dataclass
class Difference:
    """One changed call and the verdict on why it changed."""

    key: str = ""
    change: str = ""                # changed | added | removed
    before: Optional[Call] = None
    after: Optional[Call] = None
    classification: str = ""
    moved: List[str] = field(default_factory=list)


def diff_panels(before: Sequence[Call], after: Sequence[Call],
                before_state: Dict[str, str],
                after_state: Dict[str, str]) -> List[Difference]:
    """Every call that changed between two panels, each one classified."""
    by_key_before = {c.key: c for c in before}
    by_key_after = {c.key: c for c in after}
    out = []
    for key in sorted(set(by_key_before) | set(by_key_after)):
        old, new = by_key_before.get(key), by_key_after.get(key)
        if old is not None and new is not None and old.digest_row() == new.digest_row():
            continue
        change = "changed" if (old and new) else ("added" if new else "removed")
        components = set((old.components if old else []) + (new.components if new else []))
        if change != "changed":
            # Whether a call exists at all depends on the scope, so the scope belongs
            # on the row: without it a call the scope removed would print `-` under
            # `moved` and read like an unexplained disappearance.
            components.add(PANEL_SCOPE_COMPONENT)
        out.append(Difference(
            key=key, change=change, before=old, after=new,
            classification=classify_difference(old, new, before_state, after_state),
            moved=moved_components(components, before_state, after_state)))
    return out


def classification_problems(rows: Sequence["ChangeRow"]) -> List[str]:
    """Point 3 over the record: no recorded change may be left unexplained.

    The live classification runs in ``--compare``, over two call tables.  What
    ``--check`` can assert about a single committed tree is the residue that
    classification leaves behind: a row that says calls moved has to say why.
    """
    problems = []
    for row in rows:
        where = "%s:%d (%s)" % (rel(DEFAULT_CHANGELOG), row.line, row.component)
        if (row.panel_calls_changed.isdigit() and int(row.panel_calls_changed) > 0
                and not row.note):
            problems.append(
                "%s: %s call(s) changed but the note is empty. A classified difference "
                "with no explanation is a number nobody can act on."
                % (where, row.panel_calls_changed))
    return problems


# --------------------------------------------------------------------------- #
#  Point 4 -- the changelog                                                    #
# --------------------------------------------------------------------------- #
@dataclass
class ChangeRow:
    """One row of ``pin_changelog.tsv``."""

    date: str = ""
    component: str = ""
    class_name: str = ""
    from_value: str = ""
    to_value: str = ""
    panel_calls_total: str = ""
    panel_calls_changed: str = ""
    validation: str = ""
    reviewer: str = ""
    note: str = ""
    line: int = 0

    @property
    def validation_digest(self) -> str:
        text = (self.validation or "").strip()
        return text[len(VALIDATION_PREFIX):].split()[0] if text.startswith(
            VALIDATION_PREFIX) else ""


def load_changelog(path: str = DEFAULT_CHANGELOG) -> List[ChangeRow]:
    with open(path, encoding="utf-8", newline="") as fh:
        raw = [(number, line) for number, line in enumerate(fh, start=1)
               if line.strip() and not line.startswith("#")]
    reader = csv.DictReader([line for _n, line in raw], delimiter="\t")
    numbers = [n for n, line in raw][1:]
    rows = []
    for offset, record in enumerate(reader):
        if not (record.get("component") or "").strip():
            continue
        rows.append(ChangeRow(
            date=(record.get("date") or "").strip(),
            component=(record.get("component") or "").strip(),
            class_name=(record.get("class") or "").strip(),
            from_value=(record.get("from") or "").strip(),
            to_value=(record.get("to") or "").strip(),
            panel_calls_total=(record.get("panel_calls_total") or "").strip(),
            panel_calls_changed=(record.get("panel_calls_changed") or "").strip(),
            validation=(record.get("validation") or "").strip(),
            reviewer=(record.get("reviewer") or "").strip(),
            note=(record.get("note") or "").strip(),
            line=numbers[offset] if offset < len(numbers) else 0))
    return rows


def changelog_header(path: str = DEFAULT_CHANGELOG) -> List[str]:
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.startswith("#") and line.strip():
                return [c.strip() for c in line.rstrip("\n").split("\t")]
    return []


def changelog_problems(rows: Sequence[ChangeRow],
                       path: str = DEFAULT_CHANGELOG) -> List[str]:
    """The changelog must be well formed before it can be trusted as a record."""
    problems = []
    header = changelog_header(path)
    if tuple(header) != CHANGELOG_COLUMNS:
        problems.append("%s: header is %s but must be exactly %s."
                        % (rel(path), header or "(missing)", list(CHANGELOG_COLUMNS)))
    if not rows:
        problems.append("%s: no rows. A gate with an empty record cannot tell a first "
                        "pin from a silent change." % rel(path))

    previous_date = ""
    for row in rows:
        where = "%s:%d (%s)" % (rel(path), row.line, row.component)
        if row.class_name == CLASS_UNATTRIBUTABLE:
            problems.append(
                "%s: class is `unattributable`. That is a FAILURE, not a category: a "
                "call that moved while nothing pinned moved means the pipeline is not "
                "reproducible. Find the real cause -- usually a cache refreshed "
                "without its tool version changing -- and record that instead." % where)
        elif row.class_name not in CHANGE_CLASSES:
            problems.append("%s: class is %r, must be one of %s."
                            % (where, row.class_name, "/".join(CHANGE_CLASSES)))
        expected = component_class(row.component)
        if expected and row.class_name in CHANGE_CLASSES and row.class_name != expected:
            problems.append("%s: component %s is a %s pin but the row says %s."
                            % (where, row.component, expected, row.class_name))
        if not re.match(r"^\d{4}-\d{2}-\d{2}$", row.date):
            problems.append("%s: date %r is not YYYY-MM-DD." % (where, row.date))
        elif row.date < previous_date:
            problems.append("%s: date %s is earlier than the row above it (%s). The "
                            "changelog is append-ordered; add new rows at the end."
                            % (where, row.date, previous_date))
        else:
            previous_date = row.date
        if not row.reviewer:
            problems.append("%s: reviewer is empty. A pin change nobody signed for is "
                            "the thing this file exists to prevent." % where)
        if not row.to_value or not row.from_value:
            problems.append("%s: from/to must both be filled (%s means 'no pin')."
                            % (where, NO_PIN))
        if row.from_value == row.to_value:
            problems.append("%s: from and to are identical (%s), so the row records no "
                            "change." % (where, row.to_value))
        for name, value in (("panel_calls_total", row.panel_calls_total),
                            ("panel_calls_changed", row.panel_calls_changed)):
            if not value.isdigit():
                problems.append("%s: %s is %r, must be a whole number."
                                % (where, name, value))
        if (row.panel_calls_total.isdigit() and row.panel_calls_changed.isdigit()
                and int(row.panel_calls_changed) > int(row.panel_calls_total)):
            problems.append("%s: panel_calls_changed (%s) exceeds panel_calls_total (%s)."
                            % (where, row.panel_calls_changed, row.panel_calls_total))
        if not re.match(r"^[0-9a-f]{%d}$" % DIGEST_LENGTH, row.validation_digest):
            problems.append(
                "%s: validation is %r, must be `%s<%d hex>` -- the digest of the call "
                "table this pin state was accepted with. Get it from: %s"
                % (where, row.validation, VALIDATION_PREFIX, DIGEST_LENGTH, PANEL_CMD))
    return problems


def chain_problems(rows: Sequence[ChangeRow], path: str = DEFAULT_CHANGELOG) -> List[str]:
    """Each component's rows must form an unbroken chain of identities.

    A gap in the chain -- a row that starts somewhere the previous row did not end --
    means a pin state existed that nobody recorded, so a call produced in between can
    never be explained.
    """
    problems = []
    by_component: Dict[str, List[ChangeRow]] = {}
    for row in rows:
        by_component.setdefault(row.component, []).append(row)
    for component, component_rows in sorted(by_component.items()):
        if component_rows[0].from_value != NO_PIN:
            problems.append(
                "%s:%d: the first row for %s starts from %r; the first row for a "
                "component must start from %s, because there was no pin before it."
                % (rel(path), component_rows[0].line, component,
                   component_rows[0].from_value, NO_PIN))
        for previous, row in zip(component_rows, component_rows[1:]):
            if row.from_value != previous.to_value:
                problems.append(
                    "%s:%d: %s starts from %r but the row above it (line %d) ended at "
                    "%r. Every pin state has to be recorded, so the chain cannot skip."
                    % (rel(path), row.line, component, row.from_value,
                       previous.line, previous.to_value))
    return problems


def latest_by_component(rows: Sequence[ChangeRow]) -> Dict[str, ChangeRow]:
    """The last recorded row per component (the changelog is append-ordered)."""
    out: Dict[str, ChangeRow] = {}
    for row in rows:
        out[row.component] = row
    return out


def pin_change_problems(state: Dict[str, str], rows: Sequence[ChangeRow],
                        panel: Optional[Panel] = None,
                        path: str = DEFAULT_CHANGELOG) -> List[str]:
    """Point 4: every pin the tree holds must be the one the changelog last recorded.

    This is the whole gate in one comparison.  Move a pin without adding a row and
    the recorded ``to`` no longer equals what the tree says, so the check fails and
    names the row to add.
    """
    problems = []
    latest = latest_by_component(rows)
    total = str(len(panel.calls)) if panel is not None else "<n>"
    digest = "%s%s" % (VALIDATION_PREFIX, panel.digest) if panel is not None \
        else "%s<digest>" % VALIDATION_PREFIX
    today = datetime.date.today().isoformat()

    for component in sorted(state):
        identity = state[component]
        row = latest.get(component)
        if row is None:
            problems.append(
                "the tree pins %s = %r but %s has no row for it. Every pinned component "
                "needs a first row: date=%s component=%s class=%s from=%s to=%r "
                "panel_calls_total=%s panel_calls_changed=0 validation=%s "
                "reviewer=<you> note=<why it is pinned>."
                % (component, identity, rel(path), today, component,
                   component_class(component) or "<class>", NO_PIN, identity, total,
                   digest))
        elif row.to_value != identity:
            problems.append(
                "%s changed: the tree says %r but %s last records %r (line %d). Add a "
                "row: date=%s component=%s class=%s from=%r to=%r "
                "panel_calls_total=%s panel_calls_changed=<n> validation=%s "
                "reviewer=<you> note=<why>. Get <n> and the digest from `%s` before and "
                "after the change, then `%s`."
                % (component, identity, rel(path), row.to_value, row.line, today,
                   component, component_class(component) or "<class>", row.to_value,
                   identity, total, digest, PANEL_CMD, COMPARE_CMD))

    for component, row in sorted(latest.items()):
        if component not in state and row.to_value != NO_PIN:
            problems.append(
                "%s:%d records %s = %r but nothing in the tree pins it any more. If the "
                "pin was removed on purpose, add a row with to=%s; otherwise restore it."
                % (rel(path), row.line, component, row.to_value, NO_PIN))
    return problems


def panel_problems(panel: Panel, rows: Sequence[ChangeRow],
                   path: str = DEFAULT_CHANGELOG) -> List[str]:
    """Point 2's verdict: does the replayed panel match the accepted one?"""
    problems = []
    if not panel.calls:
        problems.append(
            "the validation panel replayed 0 calls. The committed Tier-2 caches no "
            "longer answer the plan, so the gate is validating nothing. Rebuild them "
            "with bin/curation-triage/refresh_tier2.py and rerun %s." % CHECK_CMD)
    if panel.unevaluated:
        problems.append(
            "%d of %d panel unit(s) produced no Tier-2 evidence (%s). A unit that "
            "silently drops out of the panel stops validating anything; rebuild the "
            "caches rather than accepting a smaller panel."
            % (len(panel.unevaluated), len(panel.units),
               ", ".join(panel.unevaluated[:5])))
    if panel.missing:
        problems.append(
            "%d unit(s) of the frozen panel scope are not in the policy table any more "
            "(%s). Either restore them, or -- if they left on purpose -- take them out "
            "of PANEL_SCOPE in bin/curation-triage/update_gate.py and record that as a "
            "`%s` row for %s, with the new panel digest in its validation column."
            % (len(panel.missing), ", ".join(panel.missing[:5]), CLASS_PANEL,
               PANEL_SCOPE_COMPONENT))
    if not rows:
        return problems

    last = rows[-1]
    if last.validation_digest and last.validation_digest != panel.digest:
        problems.append(
            "the validation panel moved: replaying the committed Tier-2 evidence gives "
            "%d call(s) with digest %s%s, but the last row of %s (line %d) accepted "
            "%s%s. Reproduce it with three steps: (1) on the previous revision run `%s` "
            "into before.tsv, (2) here run `%s` into after.tsv, (3) run `%s` and record "
            "the classified differences in a new %s row -- a `%s` row for %s if what "
            "moved was the panel's scope rather than a pin."
            % (len(panel.calls), VALIDATION_PREFIX, panel.digest, rel(path), last.line,
               VALIDATION_PREFIX, last.validation_digest, PANEL_CMD, PANEL_CMD,
               COMPARE_CMD, rel(path), CLASS_PANEL, PANEL_SCOPE_COMPONENT))
    if last.panel_calls_total.isdigit() and int(last.panel_calls_total) != len(panel.calls):
        problems.append(
            "the validation panel has %d call(s) but %s:%d records panel_calls_total=%s. "
            "Re-run `%s` and record the new total."
            % (len(panel.calls), rel(path), last.line, last.panel_calls_total, PANEL_CMD))
    return problems


def provenance_problems(panel: Panel, pins: Sequence[Dict[str, str]],
                        path: str = DEFAULT_PINS) -> List[str]:
    """Point 2's second half: was this panel produced BY the pinned tools?

    A digest only certifies a pin state if the evidence behind it came from that pin
    state.  Without this check a pin can move in the lockfile, the pin table and the
    manifest, get a changelog row, and be signed off by 75 calls that the PREVIOUS pin
    produced -- and the changelog would then record, permanently, a digest belonging to
    a tool that no longer exists in the tree.  The caches record their own producer, so
    the check is simply that the two agree.
    """
    problems: List[str] = []
    recorded = panel.provenance
    if not recorded:
        return ["the committed Tier-2 caches record no tool versions, so the panel "
                "cannot say which tools produced it. Rebuild them with "
                "bin/curation-triage/refresh_tier2.py, which stamps the versions it "
                "ran into each cache's `meta`."]

    by_env = {(row.get("env") or "").strip(): row for row in pins}
    for key in sorted(recorded):
        if key == AMRFINDER_DB_COMPONENT:
            amrfinder = by_env.get("amrfinderplus", {})
            expected = (amrfinder.get("database_version") or "").strip()
            if recorded[key] != expected:
                problems.append(
                    "the Tier-2 caches were built with AMRFinderPlus database %r but %s "
                    "pins %r. The panel then validates the old database, so its digest "
                    "certifies nothing about the pinned one: rebuild the caches with "
                    "bin/curation-triage/refresh_tier2.py."
                    % (recorded[key], rel(path), expected or "(nothing)"))
            continue
        env = typers_mod.TOOL_ENVS.get(key, key)
        row = by_env.get(env)
        if row is None:
            problems.append(
                "the Tier-2 caches were built with %s (recorded %s), which %s does not "
                "pin. Every tool that produced panel evidence needs a row there, or its "
                "calls can never be attributed."
                % (key, recorded[key], rel(path)))
            continue
        # The recorded exception: `shigeifinder --version` exits 2 and prints usage, so
        # the probe can only ever write `shigeifinder@unknown`.  That IS its genuine
        # identity (the conda build is what point 1 checks), so `unknown` is accepted
        # for a tool whose version_command is `none` -- and only for such a tool.  A
        # version string turning up there would mean the cache came from somewhere else.
        cannot_report = (row.get("version_command") or "").strip() == "none"
        expected_version = UNKNOWN_VERSION if cannot_report \
            else (row.get("expected_version") or "").strip()
        expected = "%s@%s" % (key, expected_version)
        aside = (" (it cannot report its own version, so `@%s` is the only value the "
                 "cache may carry)" % UNKNOWN_VERSION) if cannot_report else ""
        if recorded[key] != expected:
            problems.append(
                "the Tier-2 caches were built with %s but %s pins %s%s. The panel that "
                "signs off a pin state has to have been produced by it: rebuild the "
                "caches with bin/curation-triage/refresh_tier2.py, or put the pin back."
                % (recorded[key], rel(path), expected, aside))

    covered = {typers_mod.TOOL_ENVS.get(key, key) for key in recorded
               if key != AMRFINDER_DB_COMPONENT}
    needed = {c.split(":", 1)[1] for call in panel.calls for c in call.components
              if c.startswith("tool:")}
    for env in sorted(needed - covered):
        problems.append(
            "panel calls depend on tool:%s but the committed caches record no version "
            "for it, so nothing says the pinned build produced them. Rebuild the caches "
            "with bin/curation-triage/refresh_tier2.py." % env)
    if AMRFINDER_DB_COMPONENT not in recorded and any(
            AMRFINDER_DB_COMPONENT in call.components for call in panel.calls):
        problems.append(
            "panel calls depend on %s but the caches record no "
            "amrfinder_database_version. That database moves without its conda package, "
            "so an unrecorded one makes every call it produced unattributable."
            % AMRFINDER_DB_COMPONENT)
    return problems


def previous_state(state: Dict[str, str], rows: Sequence[ChangeRow]) -> Dict[str, str]:
    """The last RECORDED pin state: what the tree held before the change being compared.

    Reconstructed from the changelog so ``--compare`` needs no second copy of
    ``tool_pins.tsv``.  The rollback target is each component's last recorded ``to``,
    not its ``from``, because of WHEN ``--compare`` runs: step 3 of the documented
    loop, after the pin moved and before the new row exists.  At that moment the last
    row still records the identity the "before" call table was produced with, and its
    ``from`` is one state older still -- rolling back to it would put the before-state
    a whole change too early, which reads as "nothing moved" for every component that
    has never moved twice.  That is not a cosmetic error: with nothing moved, every
    real pin bump comes out ``unattributable`` and the documented recipe hard-stops.
    """
    out = dict(state)
    for component, row in latest_by_component(rows).items():
        if row.to_value == NO_PIN:
            # The last record says there was no such pin.  Dropping it (rather than
            # writing an identity of `-`) is what makes a re-added pin read as `added`
            # instead of inventing a component that existed with no identity.
            out.pop(component, None)
        else:
            out[component] = row.to_value
    return out


# --------------------------------------------------------------------------- #
#  Commands                                                                    #
# --------------------------------------------------------------------------- #
def require(path: str, what: str) -> None:
    """A missing input is a gate failure with a sentence, not a traceback."""
    if not os.path.exists(path):
        raise SystemExit("%s is missing, so the gate has nothing to check. %s"
                         % (rel(path), what))


def run_check(args) -> int:
    """The CI gate: all four points, offline, over the committed tree."""
    require(args.pins, "It is the table of pinned tools; restore it from the "
                       "repository before running the gate.")
    require(args.changelog,
            "It is the record every pin state is checked against. Seed it with one "
            "row per component from `python3 bin/curation-triage/update_gate.py "
            "--print-state`, each with from=%s and the digest that `%s` reports."
            % (NO_PIN, PANEL_CMD))
    pins = load_pins(args.pins)
    state = pin_state(args.pins, args.manifest)
    panel = replay_panel(src=args.src, cache_dir=args.cache_dir, state=state)
    rows = load_changelog(args.changelog)

    results = [
        ("1. immutable digest", digest_problems(pins, REPO_ROOT, panel)),
        ("2. validation panel", panel_problems(panel, rows, args.changelog)
         + provenance_problems(panel, pins, args.pins)),
        ("3. classification", classification_problems(rows)),
        ("4. changelog", changelog_problems(rows, args.changelog)
         + chain_problems(rows, args.changelog)
         + pin_change_problems(state, rows, panel, args.changelog)),
    ]

    print("=" * 78)
    print("Update gate -- %d pinned component(s), %d panel call(s) over %d unit(s)"
          % (len(state), len(panel.calls), len(panel.units)))
    print("=" * 78)
    failed = 0
    for name, problems in results:
        print("%-22s %s" % (name, "PASS" if not problems else
                            "FAIL (%d)" % len(problems)))
        for problem in problems:
            failed += 1
            print("    - %s" % problem)
    print("-" * 78)
    if failed:
        print("%d problem(s). Fix them and rerun: %s" % (failed, CHECK_CMD))
        return 1
    print("panel digest %s%s, accepted by %s:%d"
          % (VALIDATION_PREFIX, panel.digest, rel(args.changelog),
             rows[-1].line if rows else 0))
    print("all four points pass.")
    return 0


def run_panel(args) -> int:
    state = pin_state(args.pins, args.manifest)
    panel = replay_panel(src=args.src, cache_dir=args.cache_dir, state=state)
    log("panel: %d call(s) over %d unit(s); %d unit(s) had no evidence, %d scope "
        "unit(s) missing from the policy table"
        % (len(panel.calls), len(panel.units), len(panel.unevaluated),
           len(panel.missing)))
    if args.out:
        write_panel(panel, args.out)
        log("wrote %s" % args.out)
    else:
        for call in panel.calls:
            print("\t".join([call.unit_id, call.check, call.key, call.verdict,
                             _cell(call.call), call.flags]))
    print("%s%s\t%d" % (VALIDATION_PREFIX, panel.digest, len(panel.calls)))
    return 0


def run_compare(args) -> int:
    before_calls = read_panel(args.compare[0])
    after_calls = read_panel(args.compare[1])
    after_state = pin_state(args.pins_after or args.pins, args.manifest)
    if args.pins_before:
        before_state = pin_state(args.pins_before, args.manifest)
    else:
        before_state = previous_state(after_state, load_changelog(args.changelog))

    differences = diff_panels(before_calls, after_calls, before_state, after_state)
    moved = moved_components(sorted(set(before_state) | set(after_state)),
                             before_state, after_state)
    print("pinned components moved: %d%s"
          % (len(moved), (" (%s)" % ", ".join(moved)) if moved else ""))
    print("calls: %d before, %d after, %d changed"
          % (len(before_calls), len(after_calls), len(differences)))
    print("\t".join(("key", "change", "class", "before", "after", "moved")))
    counts: Dict[str, int] = {}
    for difference in differences:
        counts[difference.classification] = counts.get(difference.classification, 0) + 1
        print("\t".join((difference.key, difference.change, difference.classification,
                         describe(difference.before), describe(difference.after),
                         ",".join(difference.moved) or "-")))
    print("-" * 78)
    for class_name in CHANGE_CLASSES + (CLASS_UNATTRIBUTABLE,):
        if counts.get(class_name):
            print("%-16s %d" % (class_name, counts[class_name]))

    if before_state.get(PANEL_SCOPE_COMPONENT) != after_state.get(PANEL_SCOPE_COMPONENT):
        # Said plainly and separately, because the failure this replaces was a false
        # alarm: a call that came or went with the scope is a curation decision, and
        # calling it a reproducibility failure trains reviewers to ignore the gate.
        print("-" * 78)
        print("the panel's scope changed: %s -> %s. The %d call(s) classified `%s` "
              "appeared or disappeared because of that, not because a pinned tool "
              "answered differently. Record it as a %s row for %s in %s."
              % (before_state.get(PANEL_SCOPE_COMPONENT) or "(unrecorded)",
                 after_state.get(PANEL_SCOPE_COMPONENT) or "(unrecorded)",
                 counts.get(CLASS_PANEL, 0), CLASS_PANEL, CLASS_PANEL,
                 PANEL_SCOPE_COMPONENT, rel(DEFAULT_CHANGELOG)))

    unattributable = counts.get(CLASS_UNATTRIBUTABLE, 0)
    if unattributable:
        print("-" * 78)
        print("FAIL: %d call(s) changed with no pinned component to explain them. That "
              "is not a category, it is a reproducibility failure: the same pins "
              "produced two different answers. Look first at a cache refreshed without "
              "its tool version changing, then at an input edited outside %s."
              % (unattributable, rel(DEFAULT_PINS)))
        return 1
    if differences:
        print("-" * 78)
        print("Record these in %s: one row per moved component, with "
              "panel_calls_total=%d panel_calls_changed=%d."
              % (rel(DEFAULT_CHANGELOG), len(after_calls), len(differences)))
    return 0


def run_print_state(args) -> int:
    for component, identity in sorted(pin_state(args.pins, args.manifest).items()):
        print("\t".join((component, component_class(component), identity)))
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="The four-point gate a pin change has to pass.")
    p.add_argument("--check", action="store_true",
                   help="Run all four points over the committed tree (what CI runs). "
                        "Offline: no tool is executed and nothing is fetched.")
    p.add_argument("--panel", action="store_true",
                   help="Replay the Tier-2 validation panel and write its call table.")
    p.add_argument("--compare", nargs=2, metavar=("BEFORE", "AFTER"), default=None,
                   help="Classify every call that differs between two call tables.")
    p.add_argument("--print-state", action="store_true",
                   help="Print every pinned component and its current identity.")
    p.add_argument("--out", default=None, help="Where --panel writes its call table.")
    p.add_argument("--pins", default=DEFAULT_PINS)
    p.add_argument("--pins-before", default=None,
                   help="tool_pins.tsv from the earlier state (default: reconstruct it "
                        "from the changelog's last recorded change per component).")
    p.add_argument("--pins-after", default=None,
                   help="tool_pins.tsv from the later state (default: --pins).")
    p.add_argument("--changelog", default=DEFAULT_CHANGELOG)
    p.add_argument("--manifest", default=DEFAULT_MANIFEST)
    p.add_argument("--src", default=SRC,
                   help="Directory holding the committed policy/manifest/panel inputs.")
    p.add_argument("--cache-dir", default=None,
                   help="Directory holding the committed Tier-2 caches (default: --src).")
    args = p.parse_args(argv)

    modes = [args.check, args.panel, bool(args.compare), args.print_state]
    if sum(1 for m in modes if m) != 1:
        p.error("pick exactly one mode: --check, --panel, --compare or --print-state")

    if args.check:
        return run_check(args)
    if args.panel:
        return run_panel(args)
    if args.compare:
        return run_compare(args)
    return run_print_state(args)


if __name__ == "__main__":
    raise SystemExit(main())
