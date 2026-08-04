# -*- coding: utf-8 -*-
"""Shared fixtures: a tiny offline taxonomy + nuccore-record builders.

Everything here is synthetic and self-contained so the Tier 0 logic is tested
with zero network and zero dependence on the multi-hundred-MB real taxdump.
"""

import os
import sys

import pytest

# Make "import taxonomy/ncbi/policy" work from the tests dir.
BIN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BIN_DIR)

from taxonomy import TaxonomyResolver  # noqa: E402
from ncbi import NuccoreRecord  # noqa: E402
import policy  # noqa: E402


# (taxid, parent, rank, scientific_name) for the base fixture taxonomy.
# Simplified lineages: genera hang directly under Bacteria/Eukaryota, which is
# all the climb-to-species / superkingdom logic needs.
_BASE_NODES = [
    (1, 1, "no rank", "root"),
    (131567, 1, "no rank", "cellular organisms"),
    (2, 131567, "superkingdom", "Bacteria"),
    (2759, 131567, "superkingdom", "Eukaryota"),
    (10239, 1, "superkingdom", "Viruses"),
    # Listeria
    (1637, 2, "genus", "Listeria"),
    (1639, 1637, "species", "Listeria monocytogenes"),
    # Clostridium
    (1485, 2, "genus", "Clostridium"),
    (1491, 1485, "species", "Clostridium botulinum"),
    # Salmonella
    (590, 2, "genus", "Salmonella"),
    (28901, 590, "species", "Salmonella enterica"),
    (54736, 590, "species", "Salmonella bongori"),
    # Yersinia
    (629, 2, "genus", "Yersinia"),
    (630, 629, "species", "Yersinia enterocolitica"),
    (633, 629, "species", "Yersinia pseudotuberculosis"),
    (632, 629, "species", "Yersinia pestis"),
    (502800, 633, "strain", "Yersinia pseudotuberculosis IP 32953"),
    # Aliivibrio
    (511678, 2, "genus", "Aliivibrio"),
    (668, 511678, "species", "Aliivibrio fischeri"),
    # Campylobacter (hyointestinalis subspecies share one accession)
    (194, 2, "genus", "Campylobacter"),
    (9999, 194, "species", "Campylobacter hyointestinalis"),
    (91352, 9999, "subspecies", "Campylobacter hyointestinalis subsp. hyointestinalis"),
    (91353, 9999, "subspecies", "Campylobacter hyointestinalis subsp. lawsonii"),
    # Escherichia (Shigella confounding)
    (561, 2, "genus", "Escherichia"),
    (562, 561, "species", "Escherichia coli"),
    (564, 561, "species", "Escherichia fergusonii"),
    # Bacillus cereus group (anthracis confounding)
    (1386, 2, "genus", "Bacillus"),
    (86661, 1386, "species group", "Bacillus cereus group"),
    (1396, 86661, "species", "Bacillus cereus"),
    (1428, 86661, "species", "Bacillus thuringiensis"),
    (1423, 1386, "species", "Bacillus subtilis"),
    # Staphylococcus (for merged.dmp test: 46170 -> 1280)
    (1279, 2, "genus", "Staphylococcus"),
    (1280, 1279, "species", "Staphylococcus aureus"),
    # Ruminococcus (genus-level Kalamari entry; organism fallback -> species)
    (216572, 2, "family", "Oscillospiraceae"),
    (1263, 216572, "genus", "Ruminococcus"),
    (2564099, 1263, "species", "Ruminococcus bovis"),
    # Arripis (fish; Eukaryota) -- caught by the organelle guard, not superkingdom
    (270544, 2759, "species", "Arripis trutta"),
    # Deliberately unresolvable BACTERIAL leaves with no species ancestor
    # (models a partial taxdump that dropped the species node). One is
    # 'strain'-rank, one is the ambiguous 'no rank' NCBI uses for many leaves --
    # both must HARD-FAIL, not soft-demote to genus-level.
    (7000000, 2, "genus", "Faketobacter"),
    (7000001, 7000000, "strain", "Faketobacter unresolvabilis strain X"),
    (7000002, 7000000, "no rank", "Faketobacter norankensis leaf Y"),
]

# merged.dmp: old taxid -> current taxid
_MERGED = {46170: 1280}

# Kalamari custom synthetic overlay (build/nodes.dmp equivalent).
_OVERLAY_NODES = [
    (9000000, 1639, "subspecies", "Listeria monocytogenes lineage I"),
    (9000001, 1639, "subspecies", "Listeria monocytogenes lineage II"),
    (9000002, 1639, "subspecies", "Listeria monocytogenes lineage III"),
    (9000003, 1639, "subspecies", "Listeria monocytogenes lineage IV"),
    (9000004, 1491, "subspecies", "Clostridium botulinum group I"),
    (9000005, 1491, "subspecies", "Clostridium botulinum group II"),
    (9000009, 28901, "subspecies", "Salmonella enterica subsp. VIII"),
    (9000010, 28901, "subspecies", "Salmonella enterica subsp. IIa"),
    (9000011, 28901, "subspecies", "Salmonella enterica subsp. IIb"),
    (9000014, 28901, "subspecies", "Salmonella enterica subsp. IIIa"),
    (9000015, 28901, "subspecies", "Salmonella enterica subsp. IIIb"),
    (9000016, 28901, "subspecies", "Salmonella enterica subsp. IX"),
    (9000017, 28901, "subspecies", "Salmonella enterica subsp. X"),
]


def _build_resolver(nodes, merged, overlay):
    parents, ranks, names = {}, {}, {}
    for taxid, parent, rank, name in nodes:
        parents[taxid] = parent
        ranks[taxid] = rank
        names[taxid] = name
    for taxid, parent, rank, name in overlay:
        parents[taxid] = parent
        ranks[taxid] = rank
        names[taxid] = name
    return TaxonomyResolver(parents=parents, ranks=ranks, names=names, merged=dict(merged))


@pytest.fixture
def resolver():
    return _build_resolver(_BASE_NODES, _MERGED, _OVERLAY_NODES)


@pytest.fixture
def manifest():
    repo_root = os.path.dirname(os.path.dirname(BIN_DIR))  # bin/curation-triage -> repo
    path = os.path.join(repo_root, "src", "curation-triage", "confounded_manifest.tsv")
    return policy.load_confounded_manifest(path)


def make_record(accession, strain_taxid, length=3_000_000, genome="chromosome",
                moltype="dna", organism=None, assembly_acc=None):
    return NuccoreRecord(
        accession=accession,
        length=length,
        genome=genome,
        moltype=moltype,
        strain_taxid=strain_taxid,
        organism=organism,
        assembly_acc=assembly_acc,
    )


@pytest.fixture
def make_records():
    """Return a helper that builds an {accession: NuccoreRecord} dict."""
    def _make(*recs):
        return {r.accession: r for r in recs}
    return _make


def row(scientific_name, accession, taxid, parent=1, source="TEST"):
    return policy.ChromosomeRow(
        scientific_name=scientific_name,
        accession=accession,
        taxid=taxid,
        parent=parent,
        source=source,
    )
