#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tier 2 core: the ANI aligned-fraction gate, sub-species panel placement, the
regime-A marker override, and the Tier-2 leads that attach to Tier 1's ledger
(pure & offline).

Tier 2 is **confirmation and sub-species residue, never coverage.**  Tier 1 sees
every pick -- the coverage probe measures 100% ANI-report coverage -- so nothing
here exists to find leads that Tier 1 missed.  It runs on two trigger groups only:

  1. **Tier-1 survivors** -- the leads already on the worklist.  Tier 2 attaches a
     sequence-level number to NCBI's metadata verdict: our pick against the type
     strain of its declared species, with the aligned fraction as a QC gate.
  2. **Confounded units** -- mandatory, whether or not Tier 1 flagged them, because
     for these the metadata cannot answer the question at all.

Three checks, one per lead
--------------------------
* ``ani``    -- pick vs type strain (Tier-1 survivors that are not regime A).
* ``panel``  -- regime B: place the unit against its sibling panel and test that
                Kalamari's own lineage references are still distinct genomes.
* ``marker`` -- the confounded taxa: run the pinned typer and apply our versioned
                interpretation rule (``markers.py``).

**The regime-A override is structural, not advisory.**  For Shigella/E. coli,
anthracis/cereus and pestis/pseudotuberculosis, ANI is the wrong ruler, so a
regime-A unit gets **no independent ANI verdict at all** -- the ANI number is
folded into the marker lead as context (``ani_context``) and the marker call is
the verdict.  There is no code path in which an ANI number can contradict a
regime-A marker call, because only one of them is ever a verdict.

Attaching to the existing loop
------------------------------
A Tier-2 result is a **new signal on the same ledger**, not a new ledger: leads
carry ``signal='T2'``, a ``t2:`` key namespace, and flow through the same
``ledger.partition`` (delta + reviews suppression) Tier 1 uses.  A confirming
Tier-2 verdict **annotates** the Tier-1 row (via ``confirms_lead``) and never
removes it -- the tool is informational, so the SME still decides.

No network and no subprocess calls here: ``skani_runner.py`` computes the ANI,
``typers.py`` runs the typers, and both hand this module plain values.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import markers as markers_mod
from leads import GENOME_LINK_BASE, Lead, SIGNAL_T2

# Config defaults (design section 10; all overridable from the CLI).
DEFAULT_SPECIES_ANI = 95.0   # the species boundary
DEFAULT_AF_GATE = 0.5        # aligned fraction below this makes an ANI uninterpretable
DEFAULT_DEREP = 99.5         # at or above this, two genomes are the same thing
DEFAULT_EPSILON = 0.2        # ANI band around a boundary that is too close to call

CHECK_ANI = "ani"
CHECK_PANEL = "panel"
CHECK_MARKER = "marker"

VERDICT_CONFIRMED = markers_mod.VERDICT_CONFIRMED
VERDICT_REFUTED = markers_mod.VERDICT_REFUTED
VERDICT_INDETERMINATE = markers_mod.VERDICT_INDETERMINATE

# Aligned-fraction / ANI gate outcomes.
GATE_UNUSABLE = "unusable_af"        # too little alignment for the ANI to mean anything
GATE_SAME = "same_species"
GATE_BORDERLINE = "borderline"       # within epsilon of the boundary
GATE_DIFFERENT = "different_species"
# Two states that are answers rather than comparisons, both read straight from NCBI:
GATE_SELF_TYPE = "is_type_strain"      # the pick IS the declared type strain
GATE_NO_TYPE = "no_type_material"      # NCBI holds no type material for the species
GATE_SCREENED = "screened_out"         # skani declined the pair: too distant to report

# Worklist flag vocabulary for Tier 2 (the "why is this row here" column).  Kept
# separate from Tier 1's frozen six so the two signals never collide.
T2_FLAGS = {
    "t2_marker_refutes_label": "marker call contradicts the declared label",
    "t2_marker_confirms_label": "marker call agrees with the declared label",
    "t2_marker_indeterminate": "marker evidence does not decide",
    "t2_ani_refutes_species": "ANI to the type strain is below the species boundary",
    "t2_ani_confirms_species": "ANI to the type strain confirms the species",
    "t2_ani_unusable": "aligned fraction too low for the ANI to be interpretable",
    "t2_ani_borderline": "ANI sits within epsilon of the species boundary",
    "t2_panel_lineage_mismatch": "nearest lineage reference is not the declared lineage",
    "t2_panel_not_distinct": "another Kalamari lineage reference is not distinct from this one",
    "t2_panel_ambiguous": "two lineage references are equally close",
    "t2_panel_consistent": "placement against the lineage panel is consistent",
    "t2_ani_screened_out": "skani found the pair too distant to report at all",
    "t2_is_type_strain": "this pick IS the declared type strain of its species",
    "t2_no_type_material": "NCBI holds no type material for this species, so no "
                           "type-strain comparison exists",
    "t2_no_evidence": "no Tier-2 evidence available for this unit",
}

# An ANI lead's flag comes from the gate, not just the verdict: "no type material
# exists" and "the aligned fraction was too low" are both indeterminate, but they tell
# an SME completely different things and must not share a label.
_ANI_GATE_FLAG = {
    GATE_SCREENED: "t2_ani_screened_out",
    GATE_SELF_TYPE: "t2_is_type_strain",
    GATE_NO_TYPE: "t2_no_type_material",
    GATE_UNUSABLE: "t2_ani_unusable",
    GATE_BORDERLINE: "t2_ani_borderline",
}

# Verdicts that put a Tier-2 row on the SME's "needs a look" list.  A confirmation
# is recorded and rendered as an annotation, never as new work.
ACTIONABLE_T2_VERDICTS = frozenset({VERDICT_REFUTED, VERDICT_INDETERMINATE})


# --------------------------------------------------------------------------- #
#  ANI + the aligned-fraction gate                                             #
# --------------------------------------------------------------------------- #
@dataclass
class AniResult:
    """One skani comparison (query = our unit, reference = the comparison genome)."""

    query: str = ""                     # accession of our unit's assembly
    reference: str = ""                 # accession compared against
    ani: Optional[float] = None         # percent
    af_query: Optional[float] = None    # aligned fraction of the query, 0-1
    af_reference: Optional[float] = None
    reference_label: str = ""           # 'type strain of X' / a lineage label
    # skani returns NO row for a pair below its identity screen or aligned-fraction
    # floor.  That is "too distant to report", which is the opposite of "not measured",
    # so it is carried as an explicit state rather than as a missing result.
    screened_out: bool = False
    screen_note: str = ""

    @property
    def af(self) -> Optional[float]:
        """The conservative aligned fraction: the smaller of the two sides.

        Taking the minimum is deliberate.  A small genome fully contained in a much
        larger one has a high aligned fraction on its own side while covering little
        of the other, and that asymmetry is exactly the case the gate exists to
        catch, so the weaker side must decide.
        """
        vals = [v for v in (self.af_query, self.af_reference) if v is not None]
        return min(vals) if vals else None

    @property
    def usable(self) -> bool:
        return self.ani is not None and self.af is not None


def ani_gate(result: AniResult, species_ani: float = DEFAULT_SPECIES_ANI,
             af_gate: float = DEFAULT_AF_GATE,
             epsilon: float = DEFAULT_EPSILON) -> str:
    """Classify one ANI comparison against the species boundary.

    The aligned fraction is checked **first** and vetoes everything: an ANI computed
    over a sliver of the genome is not a weak signal, it is not a signal, and
    reporting it as "different species" would be a false flag on an SME's worklist.
    """
    if not result.usable:
        return GATE_UNUSABLE
    if result.af < af_gate:
        return GATE_UNUSABLE
    if result.ani >= species_ani + epsilon:
        return GATE_SAME
    if result.ani <= species_ani - epsilon:
        return GATE_DIFFERENT
    return GATE_BORDERLINE


@dataclass
class AniConfirmation:
    """Tier-2's sequence-level answer for one Tier-1 survivor."""

    unit_id: str = ""
    declared_species: str = ""
    result: Optional[AniResult] = None
    gate: str = GATE_UNUSABLE
    verdict: str = VERDICT_INDETERMINATE
    reasons: List[str] = field(default_factory=list)

    def as_evidence(self) -> Dict[str, str]:
        r = self.result
        return {
            "t2_check": CHECK_ANI,
            "verdict": self.verdict,
            "call": self.declared_species if self.verdict == VERDICT_CONFIRMED else "",
            "declared_species": self.declared_species,
            "type_strain_assembly": (r.reference if r else "") or "unresolved",
            "ani": ("%.2f" % r.ani) if (r and r.ani is not None) else "",
            "aligned_fraction": ("%.3f" % r.af) if (r and r.af is not None) else "",
            "gate": self.gate,
            "reasons": "; ".join(self.reasons),
        }


def confirm_self_type(unit_id: str, declared_species: str, accession: str) -> AniConfirmation:
    """The pick IS its species' type strain -- the strongest identity there is.

    NCBI writes ``same`` in the declared-type-assembly column for these.  Comparing
    such a genome to itself would produce a meaningless 100%, so Tier 2 skips the
    computation and records the fact, which is the better evidence anyway.
    """
    conf = AniConfirmation(unit_id=unit_id, declared_species=declared_species,
                           gate=GATE_SELF_TYPE, verdict=VERDICT_CONFIRMED)
    conf.reasons.append("this pick (%s) IS the declared type strain of %s, per NCBI's "
                        "ANI report -- no comparison needed"
                        % (accession, declared_species or "its species"))
    return conf


def confirm_no_type_material(unit_id: str, declared_species: str) -> AniConfirmation:
    """No type strain exists for the declared species, so no ANI check is possible.

    This is an explanation, not a gap in our tooling: it is also why NCBI's own
    taxonomy check returns Inconclusive for these picks.  Reporting it converts an
    unexplained Tier-1 flag into an explained one.
    """
    conf = AniConfirmation(unit_id=unit_id, declared_species=declared_species,
                           gate=GATE_NO_TYPE)
    conf.reasons.append("NCBI holds no type material for %s, so there is no type strain "
                        "to compare against -- which is also why its taxonomy check "
                        "cannot conclude" % (declared_species or "this species"))
    return conf


def confirm_species(unit_id: str, declared_species: str, result: Optional[AniResult],
                    species_ani: float = DEFAULT_SPECIES_ANI,
                    af_gate: float = DEFAULT_AF_GATE,
                    epsilon: float = DEFAULT_EPSILON) -> AniConfirmation:
    """Our pick vs the type strain of its declared species -> confirmed / refuted.

    This attaches a precise number to NCBI's metadata verdict; it does not replace
    it.  A pick NCBI called "Inconclusive" that sits at 98.7% ANI over 82% of its
    genome against its own type strain is almost certainly fine, and saying so is
    the single most useful thing Tier 2 can tell an SME.
    """
    conf = AniConfirmation(unit_id=unit_id, declared_species=declared_species, result=result)
    if result is None:
        conf.reasons.append("no ANI result: the type strain for this species was not resolved")
        return conf
    if result.screened_out:
        # skani declined to report this pair at all.  Against a *type strain* that is
        # itself the answer: the two genomes are either below the identity screen or
        # share almost none of their sequence, and neither is compatible with being
        # the same species.  Reporting it as "no data" would lose a real finding.
        conf.gate = GATE_SCREENED
        conf.verdict = VERDICT_REFUTED
        conf.reasons.append(
            "skani reported no comparison against the type strain %s (%s): the two "
            "genomes are either below the identity screen or share too little sequence "
            "to align, and neither is compatible with the declared species"
            % (result.reference or "(unknown)", result.screen_note or "screened out"))
        return conf
    conf.gate = ani_gate(result, species_ani, af_gate, epsilon)
    if conf.gate == GATE_UNUSABLE:
        af = result.af
        conf.reasons.append(
            "aligned fraction %s below the gate %.2f: ANI not interpretable"
            % (("%.3f" % af) if af is not None else "unavailable", af_gate)
            if af is not None or result.ani is None else
            "skani returned no aligned fraction")
    elif conf.gate == GATE_SAME:
        conf.verdict = VERDICT_CONFIRMED
        conf.reasons.append("%.2f%% ANI over aligned fraction %.3f against the type strain %s"
                            % (result.ani, result.af, result.reference or "(unknown)"))
    elif conf.gate == GATE_DIFFERENT:
        conf.verdict = VERDICT_REFUTED
        conf.reasons.append("%.2f%% ANI to the type strain of %s is below the %.1f%% species "
                            "boundary (aligned fraction %.3f)"
                            % (result.ani, declared_species, species_ani, result.af))
    else:
        conf.reasons.append("%.2f%% ANI is within %.1f of the %.1f%% boundary: too close to call"
                            % (result.ani, epsilon, species_ani))
    return conf


# --------------------------------------------------------------------------- #
#  Sub-species panels                                                          #
# --------------------------------------------------------------------------- #
@dataclass
class PanelMember:
    accession: str = ""
    lineage_label: str = ""
    unit_id: str = ""
    role: str = "sibling"


def panel_path(panel_dir: str, kalamari_taxid) -> str:
    return os.path.join(panel_dir, "%s.acclist" % kalamari_taxid)


def load_panel(path: str) -> List[PanelMember]:
    """Read one ``panels/<taxid>.acclist`` (``#`` lines are documentation)."""
    members: List[PanelMember] = []
    if not path or not os.path.exists(path):
        return members
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            cells = line.split("\t")
            if cells[0].strip() == "accession":       # header
                continue
            members.append(PanelMember(
                accession=cells[0].strip(),
                lineage_label=cells[1].strip() if len(cells) > 1 else "",
                unit_id=cells[2].strip() if len(cells) > 2 else "",
                role=cells[3].strip() if len(cells) > 3 else "sibling",
            ))
    return members


@dataclass
class PanelPlacement:
    """Where a genome lands against its lineage panel."""

    unit_id: str = ""
    declared_label: str = ""
    nearest_label: str = ""
    nearest_accession: str = ""
    nearest_ani: Optional[float] = None
    margin: Optional[float] = None            # ANI gap to the runner-up
    not_distinct: List[str] = field(default_factory=list)   # members at/above derep
    verdict: str = VERDICT_INDETERMINATE
    ambiguous: bool = False
    reasons: List[str] = field(default_factory=list)

    def as_evidence(self) -> Dict[str, str]:
        return {
            "t2_check": CHECK_PANEL,
            "verdict": self.verdict,
            "call": self.nearest_label,
            "declared_lineage": self.declared_label,
            "nearest_label": self.nearest_label,
            "nearest_accession": self.nearest_accession,
            "nearest_ani": ("%.2f" % self.nearest_ani) if self.nearest_ani is not None else "",
            "margin_to_runner_up": ("%.2f" % self.margin) if self.margin is not None else "",
            "redundant_with": ",".join(sorted(self.not_distinct)),
            "reasons": "; ".join(self.reasons),
        }


def place_against_panel(unit_id: str, declared_label: str, results: Sequence[AniResult],
                        derep: float = DEFAULT_DEREP, epsilon: float = DEFAULT_EPSILON,
                        af_gate: float = DEFAULT_AF_GATE,
                        self_panel: bool = True) -> PanelPlacement:
    """Place a genome against its lineage panel and check the panel's own health.

    Two different questions, deliberately answered by one pass:

    * **Is the panel still a panel?**  A sibling at or above ``derep`` means two of
      Kalamari's own lineage references are effectively the same genome.  That is a
      curation defect nothing else in the tool can see, and with a self-panel (the
      members are Kalamari's other units) it is the primary finding.
    * **Where does this genome sit?**  The nearest member and the margin to the
      runner-up.  With a self-panel the unit's own lineage is *not* in the panel, so
      a nearest-neighbour label is a neighbourhood, not a lineage assignment -- the
      verdict stays indeterminate unless the panel carries external references.
    """
    place = PanelPlacement(unit_id=unit_id, declared_label=declared_label)
    screened = [r for r in results if r.screened_out]
    usable = [r for r in results
              if not r.screened_out and r.usable and r.af is not None
              and r.af >= af_gate and r.ani is not None]
    dropped = len(results) - len(usable) - len(screened)
    if dropped:
        place.reasons.append("%d panel comparison(s) dropped: aligned fraction below %.2f"
                             % (dropped, af_gate))
    if screened:
        # Distant siblings are the NORMAL case for a lineage panel, so this is context,
        # not a defect -- unlike the same state against a type strain.
        place.reasons.append("%d sibling(s) too distant for skani to report" % len(screened))
    if not usable:
        place.reasons.append("no usable panel comparison")
        return place

    ranked = sorted(usable, key=lambda r: r.ani, reverse=True)
    top = ranked[0]
    place.nearest_label = top.reference_label
    place.nearest_accession = top.reference
    place.nearest_ani = top.ani
    if len(ranked) > 1:
        place.margin = round(top.ani - ranked[1].ani, 4)
        place.ambiguous = place.margin < epsilon

    place.not_distinct = [r.reference_label or r.reference for r in usable if r.ani >= derep]
    if place.not_distinct:
        place.verdict = VERDICT_REFUTED
        place.reasons.append(
            "%.2f%% ANI to %s is at or above the derep threshold %.1f%%: these two Kalamari "
            "lineage references are not distinct genomes"
            % (max(r.ani for r in usable), ", ".join(place.not_distinct), derep))
        return place

    if self_panel:
        place.verdict = VERDICT_CONFIRMED
        place.reasons.append(
            "distinct from all %d sibling reference(s); nearest is %s at %.2f%%%s"
            % (len(usable), place.nearest_label or place.nearest_accession, top.ani,
               (" (margin %.2f to the runner-up)" % place.margin)
               if place.margin is not None else ""))
        if place.ambiguous:
            place.reasons.append("two siblings are within %.1f ANI of each other, so the "
                                 "nearest-neighbour label is not meaningful here" % epsilon)
        return place

    # External-reference panel: the nearest member IS a lineage assignment.
    if place.ambiguous:
        place.reasons.append("nearest two lineage references differ by only %.2f ANI (< %.1f): "
                             "placement is ambiguous" % (place.margin, epsilon))
        return place
    place.verdict = (VERDICT_CONFIRMED
                     if _same_label(place.nearest_label, declared_label) else VERDICT_REFUTED)
    place.reasons.append("nearest lineage reference is %s at %.2f%% ANI; declared %s"
                         % (place.nearest_label, top.ani, declared_label))
    return place


def _same_label(a: str, b: str) -> bool:
    return (a or "").strip().lower() == (b or "").strip().lower()


# --------------------------------------------------------------------------- #
#  Per-unit resolution (this is where the regime-A override happens)           #
# --------------------------------------------------------------------------- #
@dataclass
class UnitResolution:
    """Everything Tier 2 concluded about one unit, and which check decided it."""

    unit_id: str
    regime: str = ""
    ani: Optional[AniConfirmation] = None
    panel: Optional[PanelPlacement] = None
    marker_calls: List["markers_mod.MarkerCall"] = field(default_factory=list)
    deciding_check: str = ""
    verdict: str = VERDICT_INDETERMINATE
    overrode_ani: bool = False

    @property
    def is_actionable(self) -> bool:
        return self.verdict in ACTIONABLE_T2_VERDICTS


def resolve_unit(unit_id: str, regime: str,
                 ani: Optional[AniConfirmation] = None,
                 panel: Optional[PanelPlacement] = None,
                 marker_calls: Sequence["markers_mod.MarkerCall"] = ()) -> UnitResolution:
    """Combine the checks into one verdict, with the regime-A marker override.

    Precedence:
      * **regime A** -- a marker call that decides (confirmed or refuted) wins, and
        ``overrode_ani`` records that it did.  An indeterminate marker call does NOT
        hand the decision back to ANI: for these taxa the ANI answer was never
        admissible, so "the markers could not tell" is the honest verdict.
      * **regime B** -- the panel decides membership; a refuting marker label is
        reported alongside and downgrades a confirmed panel result to refuted,
        because a wrong label on a correctly-placed genome is still a defect.
      * **neither** -- the ANI confirmation is the verdict.
    """
    res = UnitResolution(unit_id=unit_id, regime=regime, ani=ani, panel=panel,
                         marker_calls=list(marker_calls))
    deciding = [c for c in res.marker_calls if c.verdict != VERDICT_INDETERMINATE]

    if regime == markers_mod.REGIME_A:
        res.deciding_check = CHECK_MARKER
        res.overrode_ani = ani is not None
        if deciding:
            res.verdict = (VERDICT_REFUTED
                           if any(c.verdict == VERDICT_REFUTED for c in deciding)
                           else VERDICT_CONFIRMED)
        return res

    if regime == markers_mod.REGIME_B:
        res.deciding_check = CHECK_PANEL if panel is not None else CHECK_MARKER
        base = panel.verdict if panel is not None else VERDICT_INDETERMINATE
        if any(c.verdict == VERDICT_REFUTED for c in deciding):
            res.verdict = VERDICT_REFUTED
            res.deciding_check = CHECK_MARKER
        else:
            res.verdict = base
        return res

    res.deciding_check = CHECK_ANI
    res.verdict = ani.verdict if ani is not None else VERDICT_INDETERMINATE
    return res


# --------------------------------------------------------------------------- #
#  Tier-2 leads (the attachment to Tier 1's ledger)                            #
# --------------------------------------------------------------------------- #
def marker_epoch(pins: dict) -> str:
    """Tooling epoch for a marker task: software + database + our rule.

    All three go in because all three change the call independently.  A database
    refresh is a tooling epoch by design: it re-opens dismissed Tier-2 leads so a
    changed call is re-reviewed as a DB effect rather than read as new biology.
    """
    comps = pins.get("components") or []
    classifier = next((c for c in comps if c.get("status") != "gap"), comps[0] if comps else {})
    method = classifier.get("method", {})
    db = classifier.get("database", {})
    return "%s@%s+db:%s+rule:%s" % (
        (method.get("conda_package") or method.get("name") or "unset").replace(" ", "-"),
        method.get("software_version") or "unset",
        db.get("version_or_scheme_id") or "unset",
        pins.get("interpretation_rule_version") or "unset")


def _flags_for(check: str, verdict: str, extra: Sequence[str] = ()) -> List[str]:
    base = {
        (CHECK_MARKER, VERDICT_REFUTED): "t2_marker_refutes_label",
        (CHECK_MARKER, VERDICT_CONFIRMED): "t2_marker_confirms_label",
        (CHECK_MARKER, VERDICT_INDETERMINATE): "t2_marker_indeterminate",
        (CHECK_ANI, VERDICT_REFUTED): "t2_ani_refutes_species",
        (CHECK_ANI, VERDICT_CONFIRMED): "t2_ani_confirms_species",
        (CHECK_ANI, VERDICT_INDETERMINATE): "t2_ani_unusable",
        (CHECK_PANEL, VERDICT_REFUTED): "t2_panel_not_distinct",
        (CHECK_PANEL, VERDICT_CONFIRMED): "t2_panel_consistent",
        (CHECK_PANEL, VERDICT_INDETERMINATE): "t2_panel_ambiguous",
    }[(check, verdict)]
    return [base] + [e for e in extra if e]


def _lead(pick, check: str, key_suffix: str, epoch: str, evidence: Dict[str, str],
          flags: Sequence[str], summary: str, confirms_lead: str = "") -> Lead:
    ev = dict(evidence)
    if confirms_lead:
        ev["confirms_lead"] = confirms_lead
    return Lead(
        signal=SIGNAL_T2,
        lead_key="t2:%s#%s" % (pick.unit_id, key_suffix),
        unit_id=pick.unit_id,
        scientific_name=pick.scientific_name,
        source_lane=pick.source_lane,
        epoch=epoch,
        actionable_flags=list(flags),
        flags=list(flags),
        candidate_accession=pick.parent_assembly,
        evidence=ev,
        link="%s%s/" % (GENOME_LINK_BASE, pick.parent_assembly),
        summary=summary,
    )


def ani_lead(pick, conf: AniConfirmation, epoch: str, confirms_lead: str = "") -> Lead:
    if conf.gate in (GATE_SELF_TYPE, GATE_NO_TYPE):
        summary = "%s (%s): %s" % (pick.scientific_name, pick.parent_assembly,
                                   "; ".join(conf.reasons) or "no detail")
    else:
        verb = {VERDICT_CONFIRMED: "confirms", VERDICT_REFUTED: "contradicts",
                VERDICT_INDETERMINATE: "cannot decide"}[conf.verdict]
        summary = "%s (%s): Tier-2 ANI %s the declared species -- %s" % (
            pick.scientific_name, pick.parent_assembly, verb,
            "; ".join(conf.reasons) or "no detail")
    override = _ANI_GATE_FLAG.get(conf.gate)
    flags = [override] if override else _flags_for(CHECK_ANI, conf.verdict)
    return _lead(pick, CHECK_ANI, CHECK_ANI, epoch, conf.as_evidence(),
                 flags, summary, confirms_lead)


def panel_lead(pick, place: PanelPlacement, epoch: str, confirms_lead: str = "") -> Lead:
    summary = "%s (%s): sub-species panel -- %s" % (
        pick.scientific_name, pick.parent_assembly, "; ".join(place.reasons) or "no detail")
    extra = ["t2_panel_ambiguous"] if (place.ambiguous and place.verdict == VERDICT_CONFIRMED) else []
    return _lead(pick, CHECK_PANEL, CHECK_PANEL, epoch, place.as_evidence(),
                 _flags_for(CHECK_PANEL, place.verdict, extra), summary, confirms_lead)


def marker_lead(pick, call: "markers_mod.MarkerCall", epoch: str, pins: dict,
                ani: Optional[AniConfirmation] = None, confirms_lead: str = "",
                context: Optional[Dict[str, object]] = None) -> Lead:
    """A marker lead.  For regime A the ANI number rides along as *context only*.

    ``context`` carries the fields of a tool that informs the task without deciding it
    (today: ECTyper's serotype/pathotype on the E. coli pick).  It is written first so
    that a context field can never overwrite a verdict field, and its keys are tool-
    prefixed and absent from ``T2_FINGERPRINT_KEYS``, so a changed serotype never
    re-opens a dismissed lead on its own.
    """
    evidence: Dict[str, object] = dict(context or {})
    evidence["t2_check"] = CHECK_MARKER
    evidence.update(call.as_evidence())
    evidence["verdict"] = call.verdict
    evidence["call"] = call.call
    evidence["marker_tools"] = ",".join(
        "%s %s" % (c["method"].get("name", "?"), c["method"].get("software_version", "?"))
        for c in pins.get("components", []) if c.get("status") != "gap")
    if ani is not None and ani.result is not None and ani.result.ani is not None:
        # Context, never a verdict: for regime A the marker call is the only verdict.
        evidence["ani_context"] = "%.2f%% ANI to %s (af %.3f) -- context only, the marker " \
            "call decides" % (ani.result.ani, ani.result.reference or "type strain",
                              ani.result.af if ani.result.af is not None else float("nan"))
    verb = {VERDICT_CONFIRMED: "agree with", VERDICT_REFUTED: "contradict",
            VERDICT_INDETERMINATE: "do not decide"}[call.verdict]
    summary = "%s (%s): markers %s the declared label%s -- %s" % (
        pick.scientific_name, pick.parent_assembly, verb,
        " (overrides ANI)" if call.overrides_ani else "",
        "; ".join(call.reasons) or "no detail")
    return _lead(pick, CHECK_MARKER, "marker:%s" % call.task_id, epoch, evidence,
                 _flags_for(CHECK_MARKER, call.verdict), summary, confirms_lead)


def split_by_verdict(leads: Sequence[Lead]):
    """Split Tier-2 leads into (needs review, confirmations).

    A confirmation is real evidence and is recorded in the ledger, but it is not
    work: it renders as an annotation on the Tier-1 row it confirms.
    """
    needs, confirms = [], []
    for lead in leads:
        verdict = (lead.evidence or {}).get("verdict", VERDICT_INDETERMINATE)
        (needs if verdict in ACTIONABLE_T2_VERDICTS else confirms).append(lead)
    return needs, confirms


# --------------------------------------------------------------------------- #
#  Render                                                                      #
# --------------------------------------------------------------------------- #
def _cell(s: str) -> str:
    return (s or "").replace("|", "\\|").replace("\n", " ")


def _table(leads: Sequence[Lead]) -> List[str]:
    out = ["| unit | check | verdict | detail | link |", "| --- | --- | --- | --- | --- |"]
    for lead in leads:
        ev = lead.evidence or {}
        out.append("| %s | %s | %s | %s | %s |" % (
            _cell(lead.scientific_name), ev.get("t2_check", ""), ev.get("verdict", ""),
            _cell(lead.summary), ("[genome](%s)" % lead.link) if lead.link else ""))
    return out


def render_markdown(part, run_id: str) -> str:
    """Tier-2 worklist markdown: findings first, confirmations as annotations."""
    surfaced = list(part.new) + list(part.accepted) + list(part.standing)
    needs, confirms = split_by_verdict(surfaced)
    new_needs, _ = split_by_verdict(part.new)

    lines = ["## Kalamari curation-triage — Tier 2 confirmation (run %s)" % run_id, ""]
    lines.append("**%d Tier-2 finding(s) need review** (%d of them new this run) · "
                 "%d confirmation(s) recorded · %d dismissed"
                 % (len(needs), len(new_needs), len(confirms), len(part.dismissed)))
    lines.append("")
    lines.append("_Informational triage, never a pass/fail gate. For the confounded taxa "
                 "(regime A) the marker call overrides ANI; a confirming Tier-2 verdict "
                 "annotates the Tier-1 row and never removes it._")
    lines.append("")

    lines.append("### Findings needing review")
    lines.extend(_table(needs) if needs else ["_None this run._"])
    lines.append("")

    if confirms:
        lines.append("<details><summary>Confirmations (%d) — Tier-2 agrees with the current "
                     "label</summary>" % len(confirms))
        lines.append("")
        lines.extend(_table(confirms))
        lines.append("")
        lines.append("</details>")
        lines.append("")
    return "\n".join(lines)
