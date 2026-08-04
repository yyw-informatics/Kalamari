#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tier 2 marker resolvers: the committed marker manifest + our own versioned rules
for reading a typer's output (pure & offline).

For the confounded taxa, genome-wide similarity is the wrong ruler: the meaningful
distinction lives in a small set of genes (a virulence plasmid, a toxin gene, a
chromosomal target).  So Tier 2 asks a validated typer, and this module turns the
typer's raw fields into a call.

The split this module defends
-----------------------------
A marker call is the product of three things that move independently:

    software  x  database  x  the rule we apply to the output

The first two are pinned in ``src/curation-triage/marker_manifest.tsv``.  The
third is *ours*, lives here as code, and carries its own
``interpretation_rule_version``.  Recording only the tool and the database makes a
changed call unattributable -- you cannot tell a new BTyper3 build from a new
reading of the same output.  Every rule below therefore states its version, and
every call records all three.

Regime A vs regime B (design section 5)
---------------------------------------
* **Regime A** -- Shigella/EIEC vs E. coli, B. anthracis vs the B. cereus group,
  Y. pestis vs Y. pseudotuberculosis.  ANI cannot split these, so the marker call
  **overrides** the ANI verdict (``MarkerCall.overrides_ani``).
* **Regime B** -- Listeria lineages, Salmonella subspecies/serovars, C. botulinum
  groups.  ANI/sketch already resolves membership; the typer supplies the *label*
  only, and a label mismatch is reported without overriding anything.

What the rules deliberately refuse to do
----------------------------------------
* A plasmid marker never makes a Yersinia species call on its own -- plasmid-cured
  *Y. pestis* and plasmid-bearing *Y. pseudotuberculosis* both exist.
* A Listeria serogroup is never reported as a lineage; serogroup, MLST ST/CC and
  cgMLST/LIN stay separate fields because they answer different questions.
* The *C. botulinum* genomic group is never inferred from the toxin type: no
  genome-wide group I-IV classifier is pinned yet (the manifest marks that row
  ``gap``), so the group verdict comes from the ANI panel or not at all.

This module holds **no network and no subprocess calls**.  ``typers.py`` runs the
tools and normalises their output; everything here is a pure function over that
normalised dict -- the same separation Tier 0/1 kept between ``policy.py`` and
``ncbi.py``.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

# Verdicts a marker call can return about the *declared* label of a Kalamari unit.
VERDICT_CONFIRMED = "confirmed"          # the markers agree with what we call it
VERDICT_REFUTED = "refuted"              # the markers say it is something else
VERDICT_INDETERMINATE = "indeterminate"  # the evidence does not decide

REGIME_A = "A"
REGIME_B = "B"

# Role of a manifest row within its task.
ROLE_CLASSIFIER = "classifier"
ROLE_PRE_SCREEN = "pre-screen"
ROLE_LABEL = "label"
ROLE_CONTEXT = "context"
ROLE_GAP = "gap"

_ROLE_ORDER = {ROLE_CLASSIFIER: 0, ROLE_LABEL: 1, ROLE_PRE_SCREEN: 2,
               ROLE_CONTEXT: 3, ROLE_GAP: 4}


# --------------------------------------------------------------------------- #
#  The committed manifest                                                      #
# --------------------------------------------------------------------------- #
@dataclass
class MarkerRow:
    """One (task, component) row of ``marker_manifest.tsv``."""

    task_id: str = ""
    component: str = ""
    trigger_taxids: List[int] = field(default_factory=list)
    regime: str = ""
    role: str = ""
    tool: str = ""
    conda_package: str = ""
    conda_build: str = ""
    database: str = ""
    database_version: str = ""
    db_pin_class: str = ""
    reproducibility_anchor: str = ""
    markers: List[str] = field(default_factory=list)        # 'gene:context' strings
    corroborating: List[str] = field(default_factory=list)
    thresholds: Dict[str, str] = field(default_factory=dict)
    interpretation_rule_version: str = ""
    status: str = ""
    note: str = ""

    @property
    def is_gap(self) -> bool:
        return self.status == "gap" or self.role == ROLE_GAP

    @property
    def marker_genes(self) -> List[str]:
        return [m.split(":", 1)[0] for m in self.markers]

    def pins(self) -> dict:
        """The provenance block stamped into the ledger for this component.

        Shape follows the implementation contract: what ran, what data it read, and
        whose reading of it produced the call.
        """
        return {
            "component": self.component,
            "method": {
                "name": self.tool,
                "conda_package": self.conda_package,
                "conda_build": self.conda_build,
                "software_version": self.conda_build or "unset",
            },
            "database": {
                "name": self.database,
                "version_or_scheme_id": self.database_version,
                "pin_class": self.db_pin_class,
                "anchor": self.reproducibility_anchor,
            },
            "analysis": {
                "thresholds": dict(self.thresholds),
                "interpretation_rule_version": self.interpretation_rule_version,
            },
            "status": self.status,
        }


def _ints(cell: str) -> List[int]:
    return [int(x) for x in (cell or "").split(",") if x.strip().isdigit()]


def _list(cell: str) -> List[str]:
    return [x.strip() for x in (cell or "").split(";") if x.strip()]


def _thresholds(cell: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for item in _list(cell):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def load_manifest(path: str) -> List[MarkerRow]:
    """Read ``marker_manifest.tsv`` (leading ``#`` lines are documentation)."""
    with open(path, encoding="utf-8", newline="") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    rows: List[MarkerRow] = []
    for r in csv.DictReader(lines, delimiter="\t"):
        if not (r.get("task_id") or "").strip():
            continue
        rows.append(MarkerRow(
            task_id=r["task_id"].strip(),
            component=(r.get("component") or "").strip(),
            trigger_taxids=_ints(r.get("trigger_taxids", "")),
            regime=(r.get("regime") or "").strip(),
            role=(r.get("role") or "").strip(),
            tool=(r.get("tool") or "").strip(),
            conda_package=(r.get("conda_package") or "").strip(),
            conda_build=(r.get("conda_build") or "").strip(),
            database=(r.get("database") or "").strip(),
            database_version=(r.get("database_version") or "").strip(),
            db_pin_class=(r.get("db_pin_class") or "").strip(),
            reproducibility_anchor=(r.get("reproducibility_anchor") or "").strip(),
            markers=_list(r.get("markers", "")),
            corroborating=_list(r.get("corroborating", "")),
            thresholds=_thresholds(r.get("thresholds", "")),
            interpretation_rule_version=(r.get("interpretation_rule_version") or "").strip(),
            status=(r.get("status") or "").strip(),
            note=(r.get("note") or "").strip(),
        ))
    return rows


def rows_for_task(manifest: Sequence[MarkerRow], task_id: str) -> List[MarkerRow]:
    """Every component of one task, classifier first."""
    rows = [r for r in manifest if r.task_id == task_id]
    return sorted(rows, key=lambda r: (_ROLE_ORDER.get(r.role, 9), r.component))


def tasks_for_taxids(manifest: Sequence[MarkerRow], taxids: Sequence[int],
                     decisive_only: bool = True) -> List[str]:
    """Task ids triggered by any of ``taxids`` (a unit's species + group triggers).

    Taxid ``0`` in the manifest means "generic fallback, never auto-triggered", so
    it is excluded here and can only be selected explicitly.

    ``decisive_only`` keeps just the tasks that can produce a VERDICT -- those with a
    regime A or B row.  Context-only tasks (ECTyper's serotype, the generic screens)
    are deliberately left out: they are useful background but never a judgement, and
    routing them through the lead path would put a verdict-less row on the worklist
    for every E. coli pick.  Pass ``False`` to list them too.
    """
    want = {int(t) for t in taxids if t}
    seen: List[str] = []
    for r in manifest:
        if r.task_id in seen:
            continue
        if not (want & {t for t in r.trigger_taxids if t}):
            continue
        if decisive_only and task_regime(manifest, r.task_id) not in (REGIME_A, REGIME_B):
            continue
        seen.append(r.task_id)
    return seen


def task_regime(manifest: Sequence[MarkerRow], task_id: str) -> str:
    """The regime of a task (its classifier row decides; '' if context-only)."""
    for r in rows_for_task(manifest, task_id):
        if r.regime in (REGIME_A, REGIME_B):
            return r.regime
    return ""


def task_pins(manifest: Sequence[MarkerRow], task_id: str) -> dict:
    """Tool + database + rule pins for every component of a task."""
    rows = rows_for_task(manifest, task_id)
    return {
        "task_id": task_id,
        "regime": task_regime(manifest, task_id),
        "interpretation_rule_version": rows[0].interpretation_rule_version if rows else "",
        "components": [r.pins() for r in rows],
    }


# --------------------------------------------------------------------------- #
#  The call                                                                    #
# --------------------------------------------------------------------------- #
@dataclass
class MarkerCall:
    """What the markers say, and what that means for the unit's declared label."""

    task_id: str
    rule_version: str
    call: str = ""                # organism / label the evidence supports
    verdict: str = VERDICT_INDETERMINATE
    confidence: str = "low"       # high | medium | low
    overrides_ani: bool = False   # regime A only
    fields: Dict[str, str] = field(default_factory=dict)   # un-merged output fields
    reasons: List[str] = field(default_factory=list)

    def as_evidence(self) -> Dict[str, str]:
        """Flat, human-readable evidence for the ledger record."""
        out = {
            "marker_task": self.task_id,
            "marker_call": self.call,
            "marker_verdict": self.verdict,
            "marker_confidence": self.confidence,
            "interpretation_rule_version": self.rule_version,
            "overrides_ani": "yes" if self.overrides_ani else "no",
            "marker_reasons": "; ".join(self.reasons),
        }
        out.update({("marker_%s" % k): v for k, v in sorted(self.fields.items())})
        return out


def _present(evidence: Dict, key: str) -> List[str]:
    """A marker-presence list from the normalised typer output."""
    val = evidence.get(key)
    if isinstance(val, (list, tuple)):
        return [str(v) for v in val if str(v).strip()]
    if isinstance(val, str) and val.strip():
        return [v.strip() for v in val.split(";") if v.strip()]
    return []


def _norm(s: Optional[str]) -> str:
    return (s or "").strip().lower()


def _same_epithet(a: str, b: str) -> bool:
    """Compare two organism names on their species epithet alone.

    A tool may report a bare epithet ('paranthracis') where we hold a binomial
    ('Bacillus cereus'), so the last word is the only reliably comparable part.
    """
    ea, eb = _norm(a).split(), _norm(b).split()
    return bool(ea and eb and ea[-1] == eb[-1])


# --------------------------------------------------------------------------- #
#  Rule: Shigella / EIEC vs Escherichia coli      (regime A)                    #
# --------------------------------------------------------------------------- #
RULE_SHIGELLA = "shigella_eiec_vs_e_coli@1"


def call_shigella_eiec(evidence: Dict, declared_species: str) -> MarkerCall:
    """ShigEiFinder output -> is this pick really E. coli, or Shigella/EIEC?

    ``ipaH`` is the entry marker for both Shigella and enteroinvasive E. coli, so
    ipaH alone separates "invasive" from "ordinary E. coli" but not Shigella from
    EIEC -- that needs the cluster-specific loci.  An ipaH hit with no cluster
    assignment is therefore reported as undecided rather than guessed.
    """
    c = MarkerCall(task_id="shigella_eiec_vs_escherichia_coli", rule_version=RULE_SHIGELLA,
                   overrides_ani=True)
    ipah = _norm(evidence.get("ipaH"))
    organism = (evidence.get("organism") or "").strip()
    cluster = (evidence.get("cluster") or "").strip()
    c.fields = {"ipaH": ipah or "unknown", "shigeifinder_organism": organism,
                "cluster": cluster, "serotype": (evidence.get("serotype") or "").strip()}

    if ipah not in ("present", "absent"):
        c.reasons.append("ShigEiFinder reported no ipaH result")
        return c
    if ipah == "absent":
        c.call = "Escherichia coli"
        c.verdict = (VERDICT_CONFIRMED if "escherichia coli" in _norm(declared_species)
                     else VERDICT_REFUTED)
        c.confidence = "high"
        c.reasons.append("ipaH absent: not Shigella and not EIEC")
        return c
    # ipaH present.
    if not cluster and not organism:
        c.reasons.append("ipaH present but no cluster assignment: Shigella vs EIEC undecided")
        c.confidence = "medium"
        return c
    c.call = organism or ("Shigella (cluster %s)" % cluster)
    c.confidence = "high" if cluster else "medium"
    c.verdict = (VERDICT_REFUTED if "escherichia coli" in _norm(declared_species)
                 and "escherichia coli" not in _norm(c.call) else VERDICT_CONFIRMED)
    c.reasons.append("ipaH present with cluster assignment %s" % (cluster or "(none)"))
    return c


# --------------------------------------------------------------------------- #
#  Rule: Bacillus anthracis vs the B. cereus group      (regime A)              #
# --------------------------------------------------------------------------- #
RULE_BACILLUS = "bacillus_anthracis_vs_b_cereus_group@1"

# The four states the chromosome + plasmid evidence can combine into.  Keeping them
# named makes the select-agent-relevant cases explicit rather than implied.
ANTHRACIS = "Bacillus anthracis"
ANTHRACIS_CURED = "Bacillus anthracis (plasmid-cured)"
CEREUS_ANTHRAX_TOXIN = "Bacillus cereus group, carrying anthrax toxin plasmid(s)"


def _is_anthracis(name: str) -> bool:
    """True only for *B. anthracis* itself.

    Compare the species epithet, never a substring: *Bacillus paranthracis* contains
    the string "anthracis" and is a different organism.  Reading it as anthracis
    would turn an ordinary *B. cereus*-group genome into a false select-agent call --
    and *B. paranthracis* is exactly what NCBI's best ANI match reports for one of
    Kalamari's own picks, so this is a live case, not a hypothetical one.
    """
    parts = _norm(name).split()
    return len(parts) >= 2 and parts[0] == "bacillus" and parts[1] == "anthracis"


def call_bacillus_anthracis(evidence: Dict, declared_species: str) -> MarkerCall:
    """BTyper3 output -> anthracis, or another member of the B. cereus group?

    BTyper3 supplies group taxonomy plus a marker profile; it does not make a
    standalone anthracis call, so the combination rule below is ours and is
    versioned as ours.  Two cases must not collapse into each other:

      * **plasmid-cured anthracis** -- anthracis chromosome, pXO1/pXO2 lost.  Still
        *B. anthracis*, and still a select agent.
      * **B. cereus carrying anthrax plasmids** (the "biovar anthracis" case) --
        not an anthracis chromosome, but carrying pXO1/pXO2.  Clinically important
        and never to be silently labelled *B. cereus*.
    """
    c = MarkerCall(task_id="bacillus_anthracis_vs_b_cereus_group", rule_version=RULE_BACILLUS,
                   overrides_ani=True)
    genomospecies = (evidence.get("genomospecies") or "").strip()
    chrom = _norm(evidence.get("anthracis_chromosome"))     # yes / no / unknown
    pxo1 = _present(evidence, "pxo1_genes")
    pxo2 = _present(evidence, "pxo2_genes")
    subspecies = (evidence.get("btyper3_subspecies") or "").strip()
    bt_genes = _present(evidence, "bt_genes")
    closest = (evidence.get("closest_type_strain") or "").strip()
    closest_ani = (evidence.get("closest_type_strain_ani") or "").strip()
    c.fields = {
        "genomospecies": genomospecies,
        "btyper3_subspecies": subspecies or "none",
        # The nearest type strain by ANI is a different question from the taxon call
        # and the two can disagree.  Keeping it visible is what explains an unhappy
        # NCBI taxonomy check to the reviewer.
        "closest_type_strain": closest or "none",
        "closest_type_strain_ani": closest_ani or "",
        "anthracis_chromosome": chrom or "unknown",
        "pxo1_markers": ",".join(pxo1) or "none",
        "pxo2_markers": ",".join(pxo2) or "none",
        "bt_toxin_genes": ",".join(bt_genes) or "none",
    }
    if chrom not in ("yes", "no"):
        c.reasons.append("BTyper3 gave no usable chromosome evidence")
        if genomospecies:
            c.call = genomospecies
            c.fields["genomospecies"] = genomospecies
        return c

    has_plasmid = bool(pxo1 or pxo2)
    if chrom == "yes":
        c.call = ANTHRACIS if has_plasmid else ANTHRACIS_CURED
        c.confidence = "high" if has_plasmid else "medium"
        c.reasons.append("anthracis chromosome signature present; pXO1=%s pXO2=%s"
                         % (",".join(pxo1) or "none", ",".join(pxo2) or "none"))
    elif has_plasmid:
        c.call = CEREUS_ANTHRAX_TOXIN
        c.confidence = "high"
        c.reasons.append("no anthracis chromosome signature, but anthrax plasmid marker(s) "
                         "present (B. cereus biovar anthracis pattern)")
    else:
        c.call = genomospecies or "Bacillus cereus group (non-anthracis)"
        c.confidence = "high" if genomospecies else "medium"
        c.reasons.append("no anthracis chromosome signature and no anthrax plasmid markers")

    declared = _norm(declared_species)
    if c.call in (ANTHRACIS, ANTHRACIS_CURED, CEREUS_ANTHRAX_TOXIN):
        c.verdict = (VERDICT_CONFIRMED
                     if _is_anthracis(declared) and c.call != CEREUS_ANTHRAX_TOXIN
                     else VERDICT_REFUTED)
    elif subspecies:
        # BTyper3 3.4 reports the revised nomenclature, where the familiar species
        # name is the SUBSPECIES epithet: our 'Bacillus cereus' pick comes back as
        # B. mosaicus subsp. cereus.  So the declared epithet is matched against the
        # subspecies and the full taxon-name list, not against the genomospecies --
        # comparing 'cereus' with 'mosaicus' would refute every correct label.
        epithet = declared.split()[-1] if declared else ""
        haystack = "%s %s" % (_norm(subspecies), _norm(genomospecies))
        c.verdict = (VERDICT_CONFIRMED
                     if epithet and re.search(r"\b%s\b" % re.escape(epithet), haystack)
                     else VERDICT_REFUTED)
        if c.verdict == VERDICT_REFUTED:
            c.reasons.append("BTyper3 calls this %s, which does not include the "
                             "declared %s" % (genomospecies or subspecies, declared_species))
    else:
        # BTyper3 assigned no subspecies and no biovar, so it cannot speak to the
        # declared NAME.  It has still answered the regime-A question -- this is not
        # anthracis -- and saying "cannot confirm the label" is the honest verdict.
        # Kalamari's B. thuringiensis reference is the real case: strain Al Hakam
        # carries no Bt toxin genes, and the biovar is defined by those genes.
        c.reasons.append("BTyper3 assigned no subspecies or biovar (Bt toxin genes: %s), "
                         "so it cannot confirm the declared %s; it does rule out "
                         "anthracis" % (",".join(bt_genes) or "none", declared_species))
    if closest and not _same_epithet(closest, declared):
        c.reasons.append("nearest type strain by ANI is %s%s -- a different organism "
                         "from the declared %s, which is what an unhappy NCBI "
                         "taxonomy check is reacting to"
                         % (closest, " (%s%%)" % closest_ani if closest_ani else "",
                            declared_species))
    return c


# --------------------------------------------------------------------------- #
#  Rule: Yersinia pestis vs Y. pseudotuberculosis      (regime A)               #
# --------------------------------------------------------------------------- #
RULE_YERSINIA = "yersinia_pestis_vs_pseudotuberculosis@1"

PESTIS = "Yersinia pestis"
PSEUDOTB = "Yersinia pseudotuberculosis"


def call_yersinia_pestis(evidence: Dict, declared_species: str,
                         min_loci_called: int = 450) -> MarkerCall:
    """cgMLST assignment (classifier) with chromosome markers as the fallback.

    The plasmid markers (caf1, pst, pla) are corroboration only, never the call:
    plasmid-cured *Y. pestis* would be called *pseudotuberculosis* and a
    plasmid-bearing *Y. pseudotuberculosis* would be called *pestis*.  That is the
    mistake the 2026-08-02 pinning pass removed by dropping abricate as the
    classifier, so the rule encodes it explicitly.
    """
    c = MarkerCall(task_id="yersinia_pestis_vs_y_pseudotuberculosis", rule_version=RULE_YERSINIA,
                   overrides_ani=True)
    assignment = (evidence.get("cgmlst_species_assignment") or "").strip()
    try:
        loci = int(evidence.get("loci_called") or 0)
    except (TypeError, ValueError):
        loci = 0
    chrom = _present(evidence, "chromosome_markers")     # ypo2088 / yihN / opgG
    plasmid = _present(evidence, "plasmid_markers")      # caf1 / pst / pla
    c.fields = {
        "cgmlst_species_assignment": assignment or "none",
        "cgmlst_loci_called": str(loci),
        "chromosome_markers": ",".join(chrom) or "none",
        "plasmid_markers": ",".join(plasmid) or "none",
    }

    if assignment and loci >= min_loci_called:
        c.call = assignment
        c.confidence = "high"
        c.reasons.append("cgMLST scheme 1 assignment over %d/%d loci" % (loci, min_loci_called))
    elif assignment:
        c.reasons.append("cgMLST assignment present but only %d loci called (< %d): not used"
                         % (loci, min_loci_called))

    if not c.call:
        pestis_chrom = [m for m in chrom if _norm(m) in ("ypo2088", "yihn")]
        ptb_chrom = [m for m in chrom if _norm(m) == "opgg"]
        if pestis_chrom and not ptb_chrom:
            c.call, c.confidence = PESTIS, "medium"
            c.reasons.append("no usable cgMLST; pestis-specific chromosome target(s) %s present"
                             % ",".join(pestis_chrom))
        elif ptb_chrom and not pestis_chrom:
            c.call, c.confidence = PSEUDOTB, "medium"
            c.reasons.append("no usable cgMLST; pseudotuberculosis-specific target opgG present")
        elif plasmid:
            c.reasons.append("plasmid markers %s only: never a species call on their own "
                             "(plasmid-cured pestis and plasmid-bearing pseudotuberculosis "
                             "both exist)" % ",".join(plasmid))
            return c
        else:
            c.reasons.append("no usable cgMLST and no discriminating chromosome target")
            return c

    if plasmid:
        c.reasons.append("plasmid corroboration: %s" % ",".join(plasmid))
    declared = _norm(declared_species)
    c.verdict = VERDICT_CONFIRMED if _norm(c.call) == declared else VERDICT_REFUTED
    return c


# --------------------------------------------------------------------------- #
#  Rule: Listeria monocytogenes serogroup / lineage / CC      (regime B)        #
# --------------------------------------------------------------------------- #
RULE_LISTERIA = "listeria_serogroup_lineage_cc@1"


def _declared_lineage(scientific_name: str) -> str:
    """'Listeria_monocytogenes_III' -> 'III' (the Kalamari lineage suffix)."""
    tail = (scientific_name or "").rsplit("_", 1)[-1].strip().upper()
    return tail if tail in ("I", "II", "III", "IV") else ""


def call_listeria_lineage(evidence: Dict, scientific_name: str) -> MarkerCall:
    """LisSero + BIGSdb-Lm -> the four un-merged fields, and a lineage check.

    Serogroup, classical MLST ST/CC and cgMLST/LIN answer different questions, so
    they are reported side by side and never merged into a single "lineage" field.
    Only the cgMLST/MLST-derived broad lineage is allowed to decide the verdict --
    a serogroup does not determine a lineage.
    """
    c = MarkerCall(task_id="listeria_monocytogenes_serogroup_lineage_cc",
                   rule_version=RULE_LISTERIA, overrides_ani=False)
    serogroup = (evidence.get("serogroup") or "").strip()
    st = (evidence.get("mlst_st") or "").strip()
    cc = (evidence.get("mlst_cc") or "").strip()
    lin = (evidence.get("cgmlst_lin") or "").strip()
    broad = (evidence.get("broad_lineage") or "").strip().upper()
    c.fields = {
        "predicted_molecular_serogroup": serogroup or "none",
        "mlst_st": st or "none",
        "mlst_clonal_complex": cc or "none",
        "lin_code": lin or "none",
        "broad_lineage": broad or "none",
    }
    declared = _declared_lineage(scientific_name)
    c.call = ("Listeria monocytogenes lineage %s" % broad) if broad else "Listeria monocytogenes"
    if not declared:
        c.reasons.append("unit is not one of the four Kalamari lineage references; "
                         "labels recorded, no lineage verdict")
        return c
    if not broad:
        c.reasons.append("no cgMLST/MLST-derived lineage: serogroup alone cannot decide a lineage")
        return c
    c.confidence = "high" if lin or cc else "medium"
    c.verdict = VERDICT_CONFIRMED if broad == declared else VERDICT_REFUTED
    c.reasons.append("declared lineage %s, cgMLST/MLST lineage %s (serogroup %s reported "
                     "separately)" % (declared, broad, serogroup or "none"))
    return c


# --------------------------------------------------------------------------- #
#  Rule: Salmonella subspecies / serovar      (regime B)                        #
# --------------------------------------------------------------------------- #
RULE_SALMONELLA = "salmonella_subspecies_serovar@1"

# Kalamari's subspecies suffix -> the subspecies name SISTR reports.  The five
# Kalamari units with no accepted subspecies name (IIa, IIb, VIII, IX, X) are
# deliberately absent: SISTR has no name to return for them, so the rule reports
# the labels and abstains rather than forcing a match.
_SALMONELLA_SUBSP = {
    "I": "enterica", "II": "salamae", "IIIA": "arizonae", "IIIB": "diarizonae",
    "IV": "houtenae", "VI": "indica", "VII": "VII",
}


def _declared_subspecies(scientific_name: str) -> str:
    return (scientific_name or "").rsplit("_", 1)[-1].strip().upper()


def call_salmonella(evidence: Dict, scientific_name: str) -> MarkerCall:
    """SISTR (+ SeqSero2 as an independent antigen call) -> subspecies and serovar.

    Two typers are run on purpose: a SISTR/SeqSero2 serovar disagreement is itself
    worth an SME's attention, so it is recorded and lowers confidence, but it never
    refutes the subspecies on its own.
    """
    c = MarkerCall(task_id="salmonella_subspecies_and_serovar", rule_version=RULE_SALMONELLA,
                   overrides_ani=False)
    subsp = (evidence.get("subspecies") or "").strip()
    sistr_serovar = (evidence.get("sistr_serovar") or "").strip()
    seqsero_serovar = (evidence.get("seqsero2_serovar") or "").strip()
    c.fields = {
        "subspecies": subsp or "none",
        "sistr_serovar": sistr_serovar or "none",
        "seqsero2_serovar": seqsero_serovar or "none",
        "serovar_agreement": ("agree" if sistr_serovar and seqsero_serovar
                              and _norm(sistr_serovar) == _norm(seqsero_serovar)
                              else "differ" if sistr_serovar and seqsero_serovar else "single"),
    }
    c.call = ("Salmonella enterica subsp. %s" % subsp) if subsp else "Salmonella enterica"
    if c.fields["serovar_agreement"] == "differ":
        c.reasons.append("SISTR serovar %s vs SeqSero2 serovar %s: the two typers disagree"
                         % (sistr_serovar, seqsero_serovar))

    declared = _declared_subspecies(scientific_name)
    expected = _SALMONELLA_SUBSP.get(declared)
    if not expected:
        c.reasons.append("Kalamari subspecies %s has no accepted subspecies name for SISTR to "
                         "return; labels recorded, no verdict" % (declared or "(none)"))
        return c
    if not subsp:
        c.reasons.append("SISTR returned no subspecies")
        return c
    c.verdict = VERDICT_CONFIRMED if _norm(subsp) == _norm(expected) else VERDICT_REFUTED
    c.confidence = "medium" if c.fields["serovar_agreement"] == "differ" else "high"
    c.reasons.append("declared subspecies %s (expects %s), SISTR reported %s"
                     % (declared, expected, subsp))
    return c


# --------------------------------------------------------------------------- #
#  Rule: C. botulinum genomic group + bont type      (regime B, with a gap)     #
# --------------------------------------------------------------------------- #
RULE_BOTULINUM = "botulinum_group_and_bont@1"


def _declared_group(scientific_name: str) -> str:
    tail = (scientific_name or "").rsplit("group", 1)[-1].strip().upper()
    return tail if tail in ("I", "II", "III", "IV") else ""


def call_botulinum(evidence: Dict, scientific_name: str) -> MarkerCall:
    """AMRFinderPlus --plus -> the bont type; the genomic group stays unresolved.

    The toxin type and the genomic group are different facts -- group I and group II
    strains can both carry bont/B, and group III carries bont/C and bont/D.  No
    genome-wide group I-IV classifier is pinned yet (the manifest marks that row
    ``gap``), so this rule records the toxin type and refuses to infer the group.
    The group check for the two Kalamari group units comes from the ANI panel.
    """
    c = MarkerCall(task_id="clostridium_botulinum_group_and_bont_type",
                   rule_version=RULE_BOTULINUM, overrides_ani=False)
    types = _present(evidence, "bont_types")
    subtype = (evidence.get("bont_subtype") or "").strip()
    c.fields = {
        "bont_gene_type": ",".join(types) or "none",
        "bont_gene_subtype": subtype or "none",
        "bont_detection_method": (evidence.get("detection_method")
                                  or "AMRFinderPlus --plus").strip(),
        "clostridium_genome_group": "not called (no pinned group I-IV classifier)",
    }
    c.call = ("Clostridium botulinum bont/%s" % ",".join(types)) if types \
        else "Clostridium botulinum (no bont gene detected)"
    c.confidence = "high" if types else "medium"
    declared = _declared_group(scientific_name)
    c.reasons.append("bont type recorded; genomic group %s NOT verified by markers -- no "
                     "pinned genome-wide group classifier (manifest row status=gap)"
                     % (declared or "(none declared)"))
    if not types:
        c.reasons.append("no bont gene detected: for a C. botulinum reference this is itself "
                         "worth an SME look")
    return c


# --------------------------------------------------------------------------- #
#  Dispatch                                                                    #
# --------------------------------------------------------------------------- #
# task_id -> the rule that reads that task's normalised typer output.
RULES = {
    "shigella_eiec_vs_escherichia_coli": call_shigella_eiec,
    "bacillus_anthracis_vs_b_cereus_group": call_bacillus_anthracis,
    "yersinia_pestis_vs_y_pseudotuberculosis": call_yersinia_pestis,
    "listeria_monocytogenes_serogroup_lineage_cc": call_listeria_lineage,
    "salmonella_subspecies_and_serovar": call_salmonella,
    "clostridium_botulinum_group_and_bont_type": call_botulinum,
}

# Rules keyed on the unit's Kalamari name (a lineage/subspecies/group suffix) rather
# than on its species name.
_NAME_KEYED = {
    "listeria_monocytogenes_serogroup_lineage_cc",
    "salmonella_subspecies_and_serovar",
    "clostridium_botulinum_group_and_bont_type",
}


def resolve(task_id: str, evidence: Dict, declared_species: str,
            scientific_name: str) -> MarkerCall:
    """Run the versioned interpretation rule for ``task_id`` over ``evidence``.

    ``evidence`` is the normalised typer output from ``typers.py`` (never raw tool
    text).  An unknown task returns an indeterminate call rather than raising, so a
    manifest row added ahead of its rule degrades to "no verdict" instead of
    breaking the monthly run.
    """
    rule = RULES.get(task_id)
    if rule is None:
        return MarkerCall(task_id=task_id, rule_version="none",
                          reasons=["no interpretation rule implemented for this task"])
    if task_id in _NAME_KEYED:
        return rule(evidence, scientific_name)
    return rule(evidence, declared_species)
