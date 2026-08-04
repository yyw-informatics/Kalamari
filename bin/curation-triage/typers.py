#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cached, mockable typer layer for Tier 2: run a pinned typer over one genome and
normalise its output into the flat evidence dict ``markers.py`` reads.

Why normalise at all.  Each typer prints its own table with its own column names,
and none of them answer our question directly -- BTyper3 reports group taxonomy and
a marker profile, not "is this anthracis".  Letting the interpretation rules read
raw tool tables would weld our judgement to one tool's output format, so a tool
swap (three of the four happened on 2026-08-02) would rewrite the rules.  Instead
each tool has a small adapter here, the rules read normalised fields, and swapping
a tool means writing one adapter.

Same offline contract as the rest of the tool: an **injectable runner** per tool
plus a **committed cache of normalised results**, keyed ``tool|accession``.  The
raw tool output is not committed -- only the handful of fields we interpret --
which keeps the committed artifact small and readable in a diff.

**Every command line and column name below was captured from a real run** of the
pinned build on a real Kalamari genome (2026-08-04), not read from documentation.
That mattered: several documented guesses were wrong, and each would have produced
a confident wrong answer rather than an obvious failure:

* AMRFinderPlus names the column ``Element symbol``, not ``Gene symbol``, and
  writes ``bont_E3`` with an underscore, not ``bont/E3``.
* BTyper3 3.4 reports the *revised* B. cereus group nomenclature -- ``species(ANI)``
  is ``mosaicus`` and the familiar name appears as ``subspecies(ANI)`` -- and writes
  its results to a file, not to stdout.
* Four of the nine tools write a file rather than printing, each with its own path
  convention; reading stdout for those yields nothing at all.
* ``mlst`` emits **no header line**; its columns are positional.
* LisSero reports gene presence as ``FULL``/``PARTIAL``/``NONE``, not yes/no.
* ShigEiFinder's ``ipaH`` is a bare ``+``/``-``, and it has no organism column: the
  species call is folded into ``CLUSTER`` and ``SEROTYPE``.

Two of the conda environments ship without ``setuptools``, so their tools die on
``import pkg_resources`` until it is installed (``shigeifinder``, ``sistr``).  That
is an environment defect rather than a code one; it is recorded in the Tier-2 doc.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
from typing import Callable, Dict, Iterable, List, Optional, Sequence

# --------------------------------------------------------------------------- #
#  Tool identifiers                                                            #
# --------------------------------------------------------------------------- #
SHIGEIFINDER = "shigeifinder"
BTYPER3 = "btyper3"
LISSERO = "lissero"
SISTR = "sistr"
SEQSERO2 = "seqsero2"
ECTYPER = "ectyper"
AMRFINDERPLUS = "amrfinderplus"
ABRICATE = "abricate"
# Classical 7-locus MLST via tseemann/mlst, named for what it IS.  This is NOT the
# 500-locus BIGSdb cgMLST scheme the Yersinia task needs, and a verified fact from the
# probe run is that the classical ST alone cannot separate Y. pestis from
# Y. pseudotuberculosis.  It is context; the cgMLST classifier remains a gap.
MLST_LISTERIA = "mlst_listeria"
MLST_YERSINIA = "mlst_yersinia"

# Which tools feed which task, in the order their fields are merged.
TASK_TOOLS: Dict[str, List[str]] = {
    "shigella_eiec_vs_escherichia_coli": [SHIGEIFINDER],
    "bacillus_anthracis_vs_b_cereus_group": [BTYPER3],
    "yersinia_pestis_vs_y_pseudotuberculosis": [MLST_YERSINIA, ABRICATE],
    "listeria_monocytogenes_serogroup_lineage_cc": [LISSERO, MLST_LISTERIA],
    "salmonella_subspecies_and_serovar": [SISTR, SEQSERO2],
    "clostridium_botulinum_group_and_bont_type": [AMRFINDERPLUS],
    "escherichia_coli_serotype_pathotype": [ECTYPER],
}

# The executable each tool installs, where it differs from the tool key.
TOOL_BINARIES = {
    AMRFINDERPLUS: "amrfinder",
    SEQSERO2: "SeqSero2_package.py",
    MLST_LISTERIA: "mlst",
    MLST_YERSINIA: "mlst",
}

# The conda env each tool lives in, where it differs from the tool key.
TOOL_ENVS = {
    MLST_LISTERIA: "mlst",
    MLST_YERSINIA: "mlst",
}

# mlst scheme names as the installed database actually spells them (`mlst --list`).
MLST_SCHEMES = {
    MLST_LISTERIA: "listeria_2",
    MLST_YERSINIA: "ypseudotuberculosis_achtman_3",
}


# --------------------------------------------------------------------------- #
#  Generic parsing helpers                                                     #
# --------------------------------------------------------------------------- #
def parse_tsv(text: str) -> List[Dict[str, str]]:
    """Header-keyed TSV rows; a leading ``#`` on the header is stripped."""
    lines = [ln for ln in (text or "").splitlines() if ln.strip()]
    if not lines:
        return []
    header = [h.strip().lstrip("#").strip() for h in lines[0].split("\t")]
    rows: List[Dict[str, str]] = []
    for line in lines[1:]:
        if line.startswith("#"):
            continue
        cells = line.split("\t")
        rows.append({h: (cells[i].strip() if i < len(cells) else "")
                     for i, h in enumerate(header)})
    return rows


def _first(rows: Sequence[Dict[str, str]]) -> Dict[str, str]:
    return dict(rows[0]) if rows else {}


def _get(row: Dict[str, str], *names: str, default: str = "") -> str:
    """First present, non-empty value among ``names`` (case-insensitive)."""
    lowered = {k.strip().lower(): v for k, v in row.items()}
    for name in names:
        val = lowered.get(name.strip().lower())
        if val:
            return val.strip()
    return default


_NAME_VALUE = re.compile(r"^\s*([^()]*?)\s*\(([^()]*)\)\s*$")


def _split_name_value(cell: str):
    """``mosaicus(95.19)`` -> ``('mosaicus', 95.19)``; tolerant of odd shapes."""
    m = _NAME_VALUE.match(cell or "")
    if not m:
        return (cell or "").strip(), None
    name, val = m.group(1).strip(), m.group(2).strip()
    try:
        return name, float(val)
    except ValueError:
        return name, None


_FRACTION_GENES = re.compile(r"^\s*(\d+)\s*/\s*(\d+)\s*\((.*)\)\s*$", re.DOTALL)


def _split_gene_hits(cell: str) -> List[str]:
    """BTyper3 gene cell ``4/4(cesA;cesB;cesC;cesD)`` -> the gene list (``0/3()`` -> [])."""
    m = _FRACTION_GENES.match(cell or "")
    if not m:
        return []
    return [g.strip() for g in m.group(3).split(";") if g.strip()]


# LisSero reports gene presence as one of these words, not as yes/no.  PARTIAL is
# deliberately NOT counted as present: a partial hit is what the tool reports when it
# is unsure, and treating it as a hit would invent a Doumith pattern.
_LISSERO_PRESENT = {"full"}


# --------------------------------------------------------------------------- #
#  Per-tool adapters: raw text -> normalised evidence fields                    #
# --------------------------------------------------------------------------- #
# Verified 2026-08-04 against shigeifinder 1.3.4-pyhdfd78af_0.
SHIGEIFINDER_COLUMNS = ("#SAMPLE", "ipaH", "VIRULENCE_PLASMID", "CLUSTER", "SEROTYPE",
                        "O_ANTIGEN", "H_ANTIGEN", "NOTES")
NOT_SHIGELLA = "not shigella/eiec"


def normalize_shigeifinder(text: str) -> Dict[str, object]:
    """ShigEiFinder has no organism column: the call lives in CLUSTER/SEROTYPE.

    On a negative genome BOTH hold the exact string ``Not Shigella/EIEC``; on a
    positive one CLUSTER holds a cluster code and SEROTYPE the serotype.  Reading
    both is safer than reading either.
    """
    row = _first(parse_tsv(text))
    ipah_raw = _get(row, "ipaH")
    cluster = _get(row, "CLUSTER")
    serotype = _get(row, "SEROTYPE")
    negative = (cluster.strip().lower() == NOT_SHIGELLA
                or serotype.strip().lower() == NOT_SHIGELLA)
    if ipah_raw in ("+", "1"):
        ipah = "present"
    elif ipah_raw in ("-", "0"):
        ipah = "absent"
    else:
        ipah = ""
    return {
        "ipaH": ipah,
        "organism": "Escherichia coli" if negative else (serotype or cluster),
        "cluster": "" if negative else cluster,
        "serotype": "" if negative else serotype,
        "virulence_plasmid_gene_count": _get(row, "VIRULENCE_PLASMID"),
        "shigeifinder_notes": _get(row, "NOTES"),
    }


# Verified against btyper3 3.4.0-pyhdfd78af_0.  BTyper3 3.4 uses the revised B. cereus
# group nomenclature: species(ANI) is the genomospecies (e.g. 'mosaicus') and the
# familiar name appears as subspecies(ANI) (e.g. 'cereus').
BTYPER3_COLUMNS = ("#filename", "prefix", "species(ANI)", "subspecies(ANI)",
                   "Pseudo_Gene_Flow_Unit(ANI)", "Closest_Type_Strain(ANI)",
                   "anthrax_toxin(genes)", "capsule_Cap(genes)", "capsule_Has(genes)",
                   "capsule_Bps(genes)", "Bt(genes)",
                   "PubMLST_ST[clonal_complex](perfect_matches)",
                   "Adjusted_panC_Group(predicted_species)", "final_taxon_names")


def normalize_btyper3(text: str) -> Dict[str, object]:
    """BTyper3 -> genomospecies, chromosome evidence, and the pXO1/pXO2 profiles.

    ``anthracis`` is matched as a whole taxon name, never as a substring:
    *B. paranthracis* contains the string and is a different organism.

    The taxon call and the nearest type strain are kept as separate fields because
    they answer different questions and can disagree -- for Kalamari's own
    *B. cereus* pick BTyper3 calls it *B. mosaicus* subsp. *cereus* while its closest
    type strain by ANI is *paranthracis*, which is exactly why NCBI's taxonomy check
    is unhappy about it.  Collapsing the two would destroy that explanation.
    """
    row = _first(parse_tsv(text))
    species, species_ani = _split_name_value(_get(row, "species(ANI)"))
    subspecies, subspecies_ani = _split_name_value(_get(row, "subspecies(ANI)"))
    closest, closest_ani = _split_name_value(_get(row, "Closest_Type_Strain(ANI)"))
    pxo1 = _split_gene_hits(_get(row, "anthrax_toxin(genes)"))
    pxo2 = _split_gene_hits(_get(row, "capsule_Cap(genes)"))
    taxa = _get(row, "final_taxon_names")

    # 'No subspecies' is BTyper3's sentinel for "not assigned", not a taxon name.
    if subspecies.strip().lower() in ("no subspecies", "na", "none", ""):
        subspecies = ""
    # The Bt toxin genes are what define the thuringiensis biovar.  Kalamari's own
    # B. thuringiensis reference (Al Hakam) is Bt-toxin-negative, so BTyper3 assigns
    # it no biovar -- a fact the rule needs in order to say "cannot confirm" instead
    # of "wrong".
    bt_genes = _split_gene_hits(_get(row, "Bt(genes)")) or []
    is_anthracis = "anthracis" in (subspecies.lower(), species.lower())
    return {
        "genomespecies_taxa": taxa,
        "genomospecies": taxa or ("%s subsp. %s" % (species, subspecies)).strip(),
        "btyper3_species": species,
        "btyper3_subspecies": subspecies,
        "btyper3_subspecies_ani": "" if subspecies_ani is None else "%.2f" % subspecies_ani,
        "btyper3_species_ani": "" if species_ani is None else "%.2f" % species_ani,
        "closest_type_strain": closest,
        "closest_type_strain_ani": "" if closest_ani is None else "%.2f" % closest_ani,
        "anthracis_chromosome": "yes" if is_anthracis else "no",
        "pxo1_genes": pxo1,
        "pxo2_genes": pxo2,
        "bt_genes": bt_genes,
        "panc_group": _get(row, "Adjusted_panC_Group(predicted_species)"),
        "pubmlst_st": _get(row, "PubMLST_ST[clonal_complex](perfect_matches)"),
    }


# Verified against lissero 0.4.10-pyhdfd78af_0.
LISSERO_COLUMNS = ("ID", "SEROTYPE", "PRS", "LMO0737", "LMO1118", "ORF2110",
                   "ORF2819", "COMMENT")


def normalize_lissero(text: str) -> Dict[str, object]:
    row = _first(parse_tsv(text))
    present = [g for g in ("PRS", "LMO0737", "LMO1118", "ORF2110", "ORF2819")
               if _get(row, g).strip().lower() in _LISSERO_PRESENT]
    return {
        "serogroup": _get(row, "SEROTYPE"),
        "doumith_pattern": ",".join(present),
        "lissero_comment": _get(row, "COMMENT"),
    }


# `mlst` prints NO header.  Columns are positional: file, scheme, ST, then one
# locus(allele) per scheme locus.
MLST_POSITIONAL = ("FILE", "SCHEME", "ST", "<locus(allele)>...")


def parse_mlst(text: str) -> Dict[str, str]:
    """Positional parse of one ``mlst`` row (there is no header line to key on)."""
    for line in (text or "").splitlines():
        if not line.strip():
            continue
        cells = line.split("\t")
        if len(cells) < 3:
            continue
        return {"file": cells[0], "scheme": cells[1], "st": cells[2],
                "loci": ";".join(cells[3:])}
    return {}


def normalize_mlst_listeria(text: str) -> Dict[str, object]:
    row = parse_mlst(text)
    return {
        "mlst_st": row.get("st", ""),
        "mlst_scheme": row.get("scheme", ""),
        "mlst_loci": row.get("loci", ""),
        # tseemann/mlst reports no clonal complex and no lineage.  The BIGSdb-Lm
        # CC->lineage mapping is a separate pin that does not exist yet, so both stay
        # empty rather than being guessed from the ST.
        "mlst_cc": "",
        "broad_lineage": "",
    }


def normalize_mlst_yersinia(text: str) -> Dict[str, object]:
    row = parse_mlst(text)
    return {
        "mlst_st": row.get("st", ""),
        "mlst_scheme": row.get("scheme", ""),
        # Deliberately NOT mapped to cgmlst_species_assignment.  Verified in the probe
        # run: the classical 7-locus ST cannot separate Y. pestis from
        # Y. pseudotuberculosis -- pestis is a clone nested inside pseudotuberculosis --
        # so feeding this in as a species call would manufacture a verdict.
        "classical_mlst_only": "yes",
    }


# Verified against sistr_cmd 1.1.3-pyhdc42f0e_2.
SISTR_COLUMNS = ("antigenic_formula", "cgmlst_ST", "cgmlst_distance", "cgmlst_found_loci",
                 "cgmlst_genome_match", "cgmlst_matching_alleles", "cgmlst_subspecies",
                 "genome", "h1", "h2", "o_antigen", "qc_status", "serogroup", "serovar",
                 "serovar_antigen", "serovar_cgmlst")


def normalize_sistr(text: str) -> Dict[str, object]:
    """SISTR emits JSON or a tab table depending on ``-f``; accept both."""
    text = (text or "").strip()
    if text.startswith(("[", "{")):
        blob = json.loads(text)
        raw = (blob[0] if isinstance(blob, list) and blob else blob) or {}
        row = {str(k): ("" if v is None else str(v)) for k, v in raw.items()}
    else:
        row = _first(parse_tsv(text))
    return {
        "subspecies": _get(row, "cgmlst_subspecies"),
        "sistr_serovar": _get(row, "serovar"),
        "sistr_serogroup": _get(row, "serogroup"),
        "cgmlst_matching_alleles": _get(row, "cgmlst_matching_alleles"),
        "cgmlst_found_loci": _get(row, "cgmlst_found_loci"),
        "sistr_qc_status": _get(row, "qc_status"),
    }


# Verified against seqsero2 1.3.2-pyhdfd78af_0.
SEQSERO2_COLUMNS = ("Sample name", "Output directory", "Input files",
                    "O antigen prediction", "H1 antigen prediction(fliC)",
                    "H2 antigen prediction(fljB)", "Predicted identification",
                    "Predicted antigenic profile", "Predicted serotype", "Note")


def normalize_seqsero2(text: str) -> Dict[str, object]:
    row = _first(parse_tsv(text))
    return {
        "seqsero2_serovar": _get(row, "Predicted serotype"),
        "seqsero2_antigenic_profile": _get(row, "Predicted antigenic profile"),
        # SeqSero2 also states the subspecies in prose, e.g. "Salmonella enterica
        # subspecies arizonae (subspecies IIIa)" -- a second, independent read on the
        # field SISTR supplies.
        "seqsero2_identification": _get(row, "Predicted identification"),
        "seqsero2_note": _get(row, "Note"),
    }


# Verified against ectyper 2.0.0-pyhdfd78af_4.
ECTYPER_COLUMNS = ("Name", "Species", "O-type", "H-type", "Serotype", "QC",
                   "DatabaseVer", "Pathotype", "StxSubtypes")


def normalize_ectyper(text: str) -> Dict[str, object]:
    row = _first(parse_tsv(text))
    return {
        "ectyper_species": _get(row, "Species"),
        "serotype": _get(row, "Serotype"),
        "o_type": _get(row, "O-type"),
        "h_type": _get(row, "H-type"),
        "ectyper_qc": _get(row, "QC"),
        "ectyper_db_version": _get(row, "DatabaseVer"),
        "pathotype": _get(row, "Pathotype"),
        "stx_subtypes": _get(row, "StxSubtypes"),
    }


# Verified against ncbi-amrfinderplus 4.2.7-hf69ffd2_0, database 2026-05-15.1.
# The column is 'Element symbol' (NOT 'Gene symbol'), and symbols use an underscore:
# 'bont_E3', not 'bont/E3'.
AMRFINDER_COLUMNS = ("Name", "Contig id", "Element symbol", "Element name", "Scope",
                     "Type", "Subtype", "Method", "% Coverage of reference",
                     "% Identity to reference")
_BONT = re.compile(r"^bont[_/]?([A-G])(\S*)$", re.IGNORECASE)


def normalize_amrfinderplus(text: str) -> Dict[str, object]:
    types, subtypes, names = [], [], []
    for row in parse_tsv(text):
        symbol = _get(row, "Element symbol", "Gene symbol")
        m = _BONT.match(symbol)
        if not m:
            continue
        letter, tail = m.group(1).upper(), m.group(2)
        if letter not in types:
            types.append(letter)
        if tail:
            subtypes.append("%s%s" % (letter, tail))
        names.append(_get(row, "Element name"))
    return {
        "bont_types": types,
        "bont_subtype": ",".join(sorted(set(subtypes))),
        "bont_element_names": "; ".join(sorted(set(n for n in names if n))),
        "detection_method": "AMRFinderPlus --plus",
    }


# Verified against abricate 1.4.0-h05cac1d_0 over vfdb + plasmidfinder.
ABRICATE_COLUMNS = ("#FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "COVERAGE",
                    "COVERAGE_MAP", "GAPS", "%COVERAGE", "%IDENTITY", "DATABASE",
                    "ACCESSION", "PRODUCT", "RESISTANCE")
# The Yersinia markers the design names, split by what they mean.  Verified fact from
# the probe run: vfdb + plasmidfinder contain NONE of pst, ypo2088, yihN or opgG, so
# the chromosome discriminators are not detectable with these databases at all.  The
# absence is recorded rather than papered over -- it is precisely why abricate is a
# pre-screen here and not the classifier.
_YERSINIA_PLASMID = ("caf1", "pst", "pla", "ymt")
_YERSINIA_CHROMOSOME = ("ypo2088", "yihn", "opgg")
_YERSINIA_CONTEXT = ("lcrv", "lcrg")


def normalize_abricate(text: str) -> Dict[str, object]:
    plasmid, chromosome, context, hits = [], [], [], []
    for row in parse_tsv(text):
        gene = _get(row, "GENE")
        if not gene:
            continue
        hits.append(gene)
        low = gene.lower()
        if low in _YERSINIA_PLASMID:
            plasmid.append(gene)
        elif low in _YERSINIA_CHROMOSOME:
            chromosome.append(gene)
        elif low in _YERSINIA_CONTEXT:
            context.append(gene)
    return {
        "plasmid_markers": plasmid,
        "chromosome_markers": chromosome,
        "virulence_context_markers": context,
        "abricate_hit_count": str(len(hits)),
    }


NORMALIZERS: Dict[str, Callable[[str], Dict[str, object]]] = {
    SHIGEIFINDER: normalize_shigeifinder,
    BTYPER3: normalize_btyper3,
    LISSERO: normalize_lissero,
    SISTR: normalize_sistr,
    SEQSERO2: normalize_seqsero2,
    ECTYPER: normalize_ectyper,
    AMRFINDERPLUS: normalize_amrfinderplus,
    ABRICATE: normalize_abricate,
    MLST_LISTERIA: normalize_mlst_listeria,
    MLST_YERSINIA: normalize_mlst_yersinia,
}


# --------------------------------------------------------------------------- #
#  Commands + where each tool leaves its answer                                 #
# --------------------------------------------------------------------------- #
STDOUT = "stdout"


def tool_binary(tool: str, env_dir: Optional[str] = None) -> str:
    """The executable for a tool, from its own conda env when one is given.

    Every typer gets its own environment because bioconda pins them against
    conflicting python and BLAST versions -- installing them together does not solve.
    ``env_dir/<env>/bin/<binary>`` is that layout; without it the tool must be on PATH.
    """
    binary = TOOL_BINARIES.get(tool, tool)
    if env_dir:
        candidate = os.path.join(env_dir, TOOL_ENVS.get(tool, tool), "bin", binary)
        if os.path.exists(candidate):
            return candidate
    return binary


def _stderr_tail(text: str, limit: int = 400) -> str:
    """The LAST part of a tool's stderr, which is where the real error is.

    Taking the head instead is how a failure gets misdiagnosed: several of these
    tools open with a deprecation warning, so the first 400 characters can be pure
    noise while the actual traceback sits at the end.
    """
    lines = [ln for ln in (text or "").strip().splitlines()
             if ln.strip() and "UserWarning" not in ln]
    return " | ".join(lines[-4:])[-limit:]


def tool_prefix(tool: str, env_dir: Optional[str] = None) -> Optional[str]:
    """The conda prefix a tool lives in (``<env_dir>/<env>``), or None."""
    if not env_dir:
        return None
    return os.path.join(env_dir, TOOL_ENVS.get(tool, tool))


def amrfinder_db_dir(env_dir: Optional[str] = None) -> str:
    """Where AMRFinderPlus keeps its database INSIDE its own pinned environment."""
    prefix = tool_prefix(AMRFINDERPLUS, env_dir)
    return os.path.join(prefix, "share", "amrfinderplus", "data", "latest") if prefix else ""


def tool_env(tool: str, env_dir: Optional[str] = None) -> Optional[Dict[str, str]]:
    """The process environment a tool needs: its own conda env, not the ambient one.

    Putting the env's ``bin`` first on PATH is the part that is obvious in hindsight,
    and the two failures it caused here are worth remembering:

      * LisSero shells out to ``makeblastdb``, which lives in the same env's ``bin``.
        Off PATH it aborts with "Could not find executable makeblastdb".
      * ``mlst`` is a Perl script.  Without the env's ``perl`` first on PATH the system
        perl runs it and cannot find ``List::MoreUtils``.

    ``CONDA_PREFIX`` is the part that is not obvious, and getting it wrong produced a
    wrong *answer* rather than a failure.  A conda tool that keeps a database beside
    itself finds it through ``CONDA_PREFIX``; inheriting the value the interactive
    shell happened to export points the tool at a **different installation**.  That is
    exactly what happened to AMRFinderPlus: the recorded database version
    ``2026-05-15.1`` was read from the base mambaforge prefix, not from the pinned env,
    and on a machine without that base prefix the same pin would resolve to something
    else or to nothing.  So the prefix is pinned to the tool's own env here, and
    ``AMRFINDER_DB`` is set explicitly when the in-env database exists -- when it does
    not, amrfinder now says so loudly instead of quietly reading a stranger's database.
    """
    prefix = tool_prefix(tool, env_dir)
    if not prefix:
        return None
    bin_dir = os.path.join(prefix, "bin")
    if not os.path.isdir(bin_dir):
        return None
    env = dict(os.environ)
    env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
    env["CONDA_PREFIX"] = prefix
    if tool == AMRFINDERPLUS:
        db = amrfinder_db_dir(env_dir)
        if os.path.isdir(db):
            env["AMRFINDER_DB"] = db
        else:
            # Never inherit a database path from the calling shell: it is how the
            # pinned tool ends up reading an unpinned database.
            env.pop("AMRFINDER_DB", None)
    return env


# Every tool key this module can run.  ``mlst_listeria``/``mlst_yersinia`` share one
# binary, so the nine distinct executables sit behind ten tool keys.
ALL_TOOLS = tuple(sorted(NORMALIZERS))


def preflight(tools: Iterable[str] = (), env_dir: Optional[str] = None) -> List[str]:
    """Resolve every tool's binary up front; return one problem string per failure.

    ``tool_binary`` falls back to the bare name when the pinned env has no such file,
    which means a missing environment does not fail -- it runs whatever the runner's
    PATH happens to offer, or dies one job at a time inside ``TyperCache.build``,
    where each failure is recorded as a miss and the run still exits 0.  For a refresh
    that is the worst outcome available: a cache committed with silently absent calls.
    This check is the loud alternative, and it is deliberately a plain list of strings
    so the caller decides whether to warn or to stop.
    """
    problems: List[str] = []
    wanted = sorted(set(tools)) or list(ALL_TOOLS)
    for tool in wanted:
        binary = TOOL_BINARIES.get(tool, tool)
        prefix = tool_prefix(tool, env_dir)
        path = tool_binary(tool, env_dir)
        if prefix and not os.path.isdir(prefix):
            problems.append("%s: no environment at %s" % (tool, prefix))
        elif path == binary:
            problems.append("%s: %s is not in the pinned environment%s (a bare-name "
                            "PATH lookup would run an unpinned build)"
                            % (tool, binary,
                               " %s" % prefix if prefix else " -- no --typer-env-dir given"))
        elif not os.access(path, os.X_OK):
            problems.append("%s: %s exists but is not executable" % (tool, path))
    return problems


def typer_command(tool: str, fasta: str, outdir: str,
                  env_dir: Optional[str] = None, threads: int = 4) -> List[str]:
    """The verified invocation for one tool over one genome."""
    b = tool_binary(tool, env_dir)
    sample = os.path.basename(fasta).rsplit(".fna", 1)[0]
    cmds = {
        # --tmpdir is effectively mandatory: the default is a RELATIVE path whose
        # parent is never created, so the run dies -- but only AFTER writing the header
        # line, leaving a header-only file that looks like a clean "no results".
        SHIGEIFINDER: [b, "-i", fasta, "-t", str(threads),
                       "--tmpdir", os.path.join(outdir, "tmp")],
        BTYPER3: [b, "-i", fasta, "-o", outdir],
        LISSERO: [b, fasta],
        SISTR: [b, "-f", "tab", "-o", os.path.join(outdir, "out"), "-t", str(threads),
                "-T", os.path.join(outdir, "tmp"), "--qc", fasta],
        SEQSERO2: [b, "-m", "k", "-t", "4", "-p", str(threads), "-i", fasta,
                   "-d", os.path.join(outdir, "seqsero")],
        ECTYPER: [b, "-i", fasta, "-o", outdir, "--verify", "-c", str(threads)],
        AMRFINDERPLUS: [b, "--nucleotide", fasta, "--plus", "--threads", str(threads),
                        "--name", sample,
                        "--output", os.path.join(outdir, "amrfinder.tsv")],
        ABRICATE: [b, "--db", "vfdb", "--threads", str(threads), fasta],
        MLST_LISTERIA: [b, "--quiet", "--threads", str(threads),
                        "--scheme", MLST_SCHEMES[MLST_LISTERIA], fasta],
        MLST_YERSINIA: [b, "--quiet", "--threads", str(threads),
                        "--scheme", MLST_SCHEMES[MLST_YERSINIA], fasta],
    }
    if tool not in cmds:
        raise KeyError("no pinned command for tool %r" % tool)
    return cmds[tool]


def output_location(tool: str, fasta: str, outdir: str) -> str:
    """Where the tool leaves its table: ``STDOUT`` or a path.

    Four of the nine write a file instead of printing, each with its own convention.
    Reading stdout for those returns an empty string, which the adapter would then
    normalise into an empty evidence dict -- a silent "no result" for a run that
    actually succeeded.
    """
    prefix = os.path.basename(fasta).rsplit(".fna", 1)[0]
    return {
        BTYPER3: os.path.join(outdir, "btyper3_final_results",
                              "%s_final_results.txt" % prefix),
        SISTR: os.path.join(outdir, "out.tab"),
        SEQSERO2: os.path.join(outdir, "seqsero", "SeqSero_result.tsv"),
        ECTYPER: os.path.join(outdir, "output.tsv"),
        AMRFINDERPLUS: os.path.join(outdir, "amrfinder.tsv"),
    }.get(tool, STDOUT)


# abricate runs one database per invocation, and the Yersinia pre-screen needs both;
# the DATABASE column keeps the merged rows distinguishable.
ABRICATE_DBS = ("vfdb", "plasmidfinder")


# Tools that need their scratch directory to exist ALREADY.  ShigEiFinder does a bare
# os.mkdir() of a uuid folder inside --tmpdir, so a missing parent kills the run --
# and it kills it AFTER printing the header line, leaving a header-only output file
# that reads as a clean "no results".
_NEEDS_TMPDIR = (SHIGEIFINDER, SISTR)


def run_typer(tool: str, fasta: str, outdir: str, env_dir: Optional[str] = None) -> str:
    """Run one typer and return its result table, wherever the tool put it."""
    os.makedirs(outdir, exist_ok=True)
    if tool in _NEEDS_TMPDIR:
        os.makedirs(os.path.join(outdir, "tmp"), exist_ok=True)
    if tool == ABRICATE:
        chunks: List[str] = []
        for db in ABRICATE_DBS:
            cmd = typer_command(tool, fasta, outdir, env_dir)
            cmd[cmd.index("--db") + 1] = db
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False,
                                  env=tool_env(tool, env_dir))  # noqa: S603
            if proc.returncode != 0:
                raise RuntimeError("abricate/%s failed (rc=%s): %s"
                                   % (db, proc.returncode, _stderr_tail(proc.stderr)))
            # Keep the header once, then append only data rows.
            chunks.append(proc.stdout if not chunks
                          else "\n".join(proc.stdout.splitlines()[1:]))
        return "\n".join(c for c in chunks if c.strip())

    cmd = typer_command(tool, fasta, outdir, env_dir)
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False,
                          env=tool_env(tool, env_dir))  # noqa: S603
    if proc.returncode != 0:
        raise RuntimeError("%s failed (rc=%s): %s" % (tool, proc.returncode,
                                                      _stderr_tail(proc.stderr)))
    dest = output_location(tool, fasta, outdir)
    if dest == STDOUT:
        return proc.stdout
    if not os.path.exists(dest):
        raise RuntimeError("%s exited 0 but wrote no result at %s" % (tool, dest))
    with open(dest, encoding="utf-8", errors="replace") as fh:
        return fh.read()


def tool_version(tool: str, env_dir: Optional[str] = None) -> str:
    """``<tool>@<version>``; ``@unknown`` when the tool reports none."""
    try:
        proc = subprocess.run([tool_binary(tool, env_dir), "--version"],
                              capture_output=True, text=True, check=False,
                              env=tool_env(tool, env_dir))  # noqa: S603
        out = (proc.stdout or proc.stderr).strip().splitlines()
        # Take the FIRST version-looking token, not the last word: ectyper prints
        # "ectyper 2.0.0 running database version 1.0", where the last word is the
        # database version and recording it as the software version would be wrong.
        ver = "unknown"
        if proc.returncode == 0 and out:
            m = re.search(r"\b(\d+\.\d+(?:\.\d+)*)\b", out[0])
            ver = m.group(1) if m else out[0].split()[-1]
    except (OSError, IndexError, KeyError):
        ver = "unknown"
    return "%s@%s" % (tool, ver)


def amrfinder_database_version(env_dir: Optional[str] = None) -> str:
    """AMRFinderPlus's dated DB version -- the best reproducibility anchor we have."""
    try:
        proc = subprocess.run([tool_binary(AMRFINDERPLUS, env_dir), "--database_version"],
                              capture_output=True, text=True, check=False,
                              env=tool_env(AMRFINDERPLUS, env_dir))  # noqa: S603
        for line in (proc.stdout or proc.stderr).splitlines():
            if "database version" in line.lower():
                return line.split(":")[-1].strip()
    except OSError:
        pass
    return "unknown"


# --------------------------------------------------------------------------- #
#  Committed cache of NORMALISED results                                       #
# --------------------------------------------------------------------------- #
class TyperCache:
    """``tool|accession`` -> normalised evidence fields, committed as JSON."""

    def __init__(self, cache_path: Optional[str] = None,
                 runner: Callable[[str, str, str], str] = None,
                 env_dir: Optional[str] = None) -> None:
        self.cache_path = cache_path
        self.env_dir = env_dir
        self._runner = runner or (
            lambda tool, fasta, outdir: run_typer(tool, fasta, outdir, env_dir))
        self.results: Dict[str, dict] = {}
        self.meta: Dict[str, object] = {}
        self.misses: List[str] = []
        if cache_path and os.path.exists(cache_path):
            with open(cache_path, encoding="utf-8") as fh:
                blob = json.load(fh)
            self.results = blob.get("results", {})
            self.meta = blob.get("meta", {})

    @staticmethod
    def key(tool: str, accession: str) -> str:
        return "%s|%s" % (tool, accession)

    def get(self, tool: str, accession: str) -> Optional[dict]:
        return self.results.get(self.key(tool, accession))

    def covers(self, tool: str, accessions: Iterable[str]) -> bool:
        return all(self.key(tool, a) in self.results for a in accessions)

    def build(self, jobs: Sequence["tuple"], store, workdir: str,
              log=lambda m: None, only_tools: Optional[Iterable[str]] = None) -> None:
        """Run each (tool, accession) job that is not already cached.

        A tool that is missing or fails is recorded as a miss and skipped: with nine
        possible typers, one broken environment must not take down the other eight.
        An adapter that returns nothing is treated as a failure too -- that means the
        output shape changed, and a silent empty result would read as "no findings".

        ``only_tools`` runs one **tool shard** of the job list.  The eleven pinned
        environments total about 13 GB, so no single runner holds them all and a
        refresh has to be split; splitting by *unit* instead would not work, because
        the sub-species panels are cliques and a unit split duplicates whole cliques.

        **Provenance is only ever refreshed by work that actually happened.**  A tool
        version IS the epoch, so stamping ``<tool>@unknown`` over a recorded version is
        not a cosmetic error: it re-baselines every Tier-2 lead and voids every
        ``reviews.tsv`` dismissal.  Two rules keep that from happening, and both are
        needed because they fail in different ways:

          * only the tools that produced a NEW result in this call are re-probed.  A
            job list is not evidence that anything ran -- on a full refresh where every
            job is a cache hit, nothing runs at all, and re-probing the whole list
            there (or the other nine environments on a shard runner) invents a tooling
            change out of an idle run.
          * a probe that comes back ``@unknown`` never replaces a version already on
            record.  ``@unknown`` is what the probe returns when the binary is simply
            not there, which is a missing tool, not a new one.

        ``shigeifinder`` is the deliberate exception to keep in mind: its
        ``--version`` exits 2, so ``shigeifinder@unknown`` is its GENUINE recorded
        value.  The guard preserves it -- it refuses to *replace* a recorded version
        with unknown, and on a first build there is nothing to replace.
        """
        if only_tools is not None:
            keep = set(only_tools)
            jobs = [(t, a) for t, a in jobs if t in keep]
        ran: List[str] = []
        for tool, accession in jobs:
            if self.key(tool, accession) in self.results:
                continue
            try:
                fasta = store.ensure(accession)
                outdir = os.path.join(workdir, "typer-run", tool, accession)
                raw = self._runner(tool, fasta, outdir)
                normalizer = NORMALIZERS.get(tool)
                if normalizer is None:
                    raise KeyError("no adapter for tool %r" % tool)
                fields = dict(normalizer(raw))
                if not any(v for v in fields.values()):
                    raise RuntimeError("adapter produced no populated field -- the "
                                       "output shape probably changed; re-check the "
                                       "column names against a real run")
                self.results[self.key(tool, accession)] = fields
                if tool not in ran:
                    ran.append(tool)
                log("%s: typed %s" % (tool, accession))
            except Exception as exc:            # noqa: BLE001 - degrade, never crash
                self.misses.append(self.key(tool, accession))
                log("%s: skipped %s (%s)" % (tool, accession, exc))
        meta = dict(self.meta)
        versions = dict(meta.get("tool_versions") or {})
        for tool in sorted(ran):
            probed = tool_version(tool, self.env_dir)
            previous = str(versions.get(tool) or "")
            if probed.endswith("@unknown") and previous:
                log("%s: version probe failed; keeping the recorded %s" % (tool, previous))
                continue
            versions[tool] = probed
        meta["tool_versions"] = versions
        meta["n_results"] = len(self.results)
        if AMRFINDERPLUS in ran:
            # Same rule for the database pin: a dated version like '2026-05-15.1' is
            # the only reproducibility anchor AMRFinderPlus gives us, and "unknown"
            # means the probe could not reach the tool, not that the pin changed.
            database = amrfinder_database_version(self.env_dir)
            recorded = str(meta.get("amrfinder_database_version") or "")
            if database == "unknown" and recorded:
                log("%s: database probe failed; keeping the recorded %s"
                    % (AMRFINDERPLUS, recorded))
            else:
                meta["amrfinder_database_version"] = database
        self.meta = meta

    def save(self) -> None:
        if not self.cache_path:
            return
        os.makedirs(os.path.dirname(os.path.abspath(self.cache_path)), exist_ok=True)
        payload = {"schema": 1, "meta": self.meta,
                   "results": dict(sorted(self.results.items()))}
        tmp = "%s.tmp.%d" % (self.cache_path, os.getpid())
        try:
            with open(tmp, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=1, sort_keys=True)
                fh.write("\n")
            os.replace(tmp, self.cache_path)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)


def evidence_for_task(task_id: str, accession: str, cache: TyperCache) -> Dict[str, object]:
    """Merge every tool's normalised fields for one task into one evidence dict.

    Merge order follows ``TASK_TOOLS``, so a later tool can fill a field the earlier
    one left blank but never overwrites a populated one -- two typers disagreeing is
    something the rules must SEE (Salmonella records it explicitly), not something
    the merge quietly resolves.
    """
    evidence: Dict[str, object] = {}
    for tool in TASK_TOOLS.get(task_id, []):
        fields = cache.get(tool, accession)
        if not fields:
            continue
        for key, val in fields.items():
            if key not in evidence or evidence[key] in ("", [], None):
                evidence[key] = val
    return evidence


def jobs_for_task(task_id: str, accession: str) -> List["tuple"]:
    return [(tool, accession) for tool in TASK_TOOLS.get(task_id, [])]


# --------------------------------------------------------------------------- #
#  Background context: tools that inform a task without deciding it            #
# --------------------------------------------------------------------------- #
# A decisive task -> the tools whose output rides along as CONTEXT for it.
#
# ECTyper is the only entry, and it is context by design, not by omission.  It types
# an E. coli serotype and pathotype, which is genuinely useful background for the
# Shigella/EIEC question, but it answers a different question and ``markers.RULES``
# has no rule for it -- routing it as its own task would return
# ``rule_version='none'`` and put a verdict-less row on the worklist for every E. coli
# pick.  So its fields are attached to the marker lead the SME is already reading.
TASK_CONTEXT_TOOLS: Dict[str, List[str]] = {
    "shigella_eiec_vs_escherichia_coli": [ECTYPER],
}

# Why the fields get a tool prefix: ``normalize_shigeifinder`` and
# ``normalize_ectyper`` BOTH emit a bare ``serotype`` key, and ``evidence_for_task``
# fills only missing keys.  Merging them naively would drop ECTyper's serotype into
# ShigEiFinder's slot, where it would read as a Shigella serotype -- a context field
# silently impersonating a verdict field.  Prefixing removes the collision entirely.
CONTEXT_ROLE_NOTE = ("background context only -- this tool decides no Tier-2 verdict "
                     "and is not part of the lead's evidence fingerprint")


def context_for_task(task_id: str, accession: str, cache: "TyperCache") -> Dict[str, object]:
    """Cached context-tool fields for one task, namespaced by the tool that made them.

    Returns ``{}`` when nothing is cached, so a unit whose context tool never ran is
    simply a lead without the extra fields -- context is never a reason to hold up a
    verdict.
    """
    out: Dict[str, object] = {}
    for tool in TASK_CONTEXT_TOOLS.get(task_id, []):
        fields = cache.get(tool, accession)
        if not fields:
            continue
        for key, val in sorted(fields.items()):
            out[key if key.startswith("%s_" % tool) else "%s_%s" % (tool, key)] = val
        out["%s_role" % tool] = CONTEXT_ROLE_NOTE
    return out
