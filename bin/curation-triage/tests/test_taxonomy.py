# -*- coding: utf-8 -*-
"""Taxonomy resolver: synthetic taxids, merged taxids, superkingdom, climbing."""

import pytest


# --- MANDATED CASE: synthetic taxids resolve to the real parent species ------
@pytest.mark.parametrize(
    "synthetic,expected_species",
    [
        (9000000, 1639),  # Listeria monocytogenes lineage I  -> L. monocytogenes
        (9000001, 1639),
        (9000002, 1639),
        (9000003, 1639),  # lineage IV
        (9000004, 1491),  # Clostridium botulinum group I     -> C. botulinum
        (9000005, 1491),  # group II
        (9000009, 28901),  # Salmonella enterica subsp. VIII  -> S. enterica
        (9000010, 28901),
        (9000011, 28901),
        (9000014, 28901),
        (9000015, 28901),
        (9000016, 28901),
        (9000017, 28901),
    ],
)
def test_synthetic_taxids_resolve_to_parent_species(resolver, synthetic, expected_species):
    res = resolver.resolve(synthetic)
    assert res.resolves_to_species
    assert res.species_taxid == expected_species
    assert res.superkingdom == "Bacteria"


def test_merged_taxid_is_followed(resolver):
    # 46170 (S. aureus subsp. aureus) was merged into 1280 (S. aureus).
    res = resolver.resolve(46170)
    assert res.species_taxid == 1280
    assert res.species_name == "Staphylococcus aureus"
    assert res.resolved_from == "merged"


def test_strain_climbs_to_species(resolver):
    res = resolver.resolve(502800)  # Y. pseudotuberculosis IP 32953 (strain)
    assert res.species_taxid == 633
    assert res.input_rank == "strain"
    assert res.is_species_or_below


def test_species_resolves_to_itself(resolver):
    res = resolver.resolve(668)  # Aliivibrio fischeri (species)
    assert res.species_taxid == 668
    assert res.resolved_from == "self"


def test_genus_has_no_species(resolver):
    res = resolver.resolve(1263)  # Ruminococcus (genus)
    assert res.species_taxid is None
    assert res.input_rank == "genus"
    assert not res.is_species_or_below  # genus is ABOVE species


def test_superkingdom_detection(resolver):
    assert resolver.superkingdom(1639) == "Bacteria"
    assert resolver.superkingdom(270544) == "Eukaryota"  # Arripis (fish)
    assert resolver.superkingdom(10239) == "Viruses"


def test_unknown_taxid_is_unresolved(resolver):
    res = resolver.resolve(424242)
    assert res.species_taxid is None
    assert res.resolved_from == "unknown"
    assert res.superkingdom is None


def test_from_dumps_roundtrip_with_overlay(tmp_path):
    """The disk-loading path merges a base dump with the synthetic overlay."""
    from taxonomy import TaxonomyResolver

    base = tmp_path / "base"
    overlay = tmp_path / "overlay"
    base.mkdir()
    overlay.mkdir()
    (base / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n"
        "2\t|\t1\t|\tsuperkingdom\t|\n"
        "1637\t|\t2\t|\tgenus\t|\n"
        "1639\t|\t1637\t|\tspecies\t|\n"
    )
    (base / "names.dmp").write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
        "1639\t|\tListeria monocytogenes\t|\t\t|\tscientific name\t|\n"
    )
    (base / "merged.dmp").write_text("999\t|\t1639\t|\n")
    (overlay / "nodes.dmp").write_text("9000000\t|\t1639\t|\tsubspecies\t|\n")
    (overlay / "names.dmp").write_text(
        "9000000\t|\tListeria monocytogenes lineage I\t|\t\t|\tscientific name\t|\n"
    )

    r = TaxonomyResolver.from_dumps(str(base), overlay_dirs=[str(overlay)])
    assert r.resolve(9000000).species_taxid == 1639  # overlay node climbs into base
    assert r.resolve(999).species_taxid == 1639       # merged taxid followed
    assert r.name(9000000) == "Listeria monocytogenes lineage I"
