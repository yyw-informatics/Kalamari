# -*- coding: utf-8 -*-
"""Programmatic eligibility guard (organelle / virus / eukaryote / eligible)."""

from conftest import make_record, row
import policy


def _one(rows, records, resolver, manifest):
    units = policy.build_units(rows, resolver, records, manifest)
    assert len(units) == 1
    return units[0]


# --- MANDATED CASE: organelle guard catches a fish mitogenome -----------------
def test_organelle_guard_catches_fish_mitogenome(resolver, manifest, make_records):
    # Arripis_trutta AP006810: a fish mitogenome a hand list already missed.
    # The guard must catch it programmatically from the 'genome' field.
    rows = [row("Arripis_trutta", "AP006810", 270544)]
    records = make_records(
        make_record("AP006810", strain_taxid=270544, length=17244, genome="mitochondrion"),
    )
    u = _one(rows, records, resolver, manifest)
    assert u.eligibility == policy.ELIG_ORGANELLE


def test_host_mitogenome_excluded_over_eukaryote(resolver, manifest, make_records):
    # A eukaryote whose *molecule* is a mitochondrion is excluded as organelle,
    # not classified 'eukaryote' -- the guard order matters.
    rows = [row("Homo_sapiens", "J01415", 270544)]  # reuse a Eukaryota taxid
    records = make_records(
        make_record("J01415", strain_taxid=270544, length=16569, genome="mitochondrion"),
    )
    u = _one(rows, records, resolver, manifest)
    assert u.eligibility == policy.ELIG_ORGANELLE


def test_oversized_assembly_excluded(resolver, manifest, make_records):
    # A grouped assembly summing above the size cap is treated as organelle/
    # eukaryote even without a 'genome' organelle flag.
    rows = [row("Big_eukaryote", "XX000001", 270544)]
    records = make_records(
        make_record("XX000001", strain_taxid=270544, length=60_000_000, genome="chromosome"),
    )
    u = _one(rows, records, resolver, manifest)
    assert u.eligibility == policy.ELIG_ORGANELLE


def test_bacteria_is_eligible(resolver, manifest, make_records):
    rows = [row("Yersinia_enterocolitica", "AM286415", 630)]
    records = make_records(make_record("AM286415", strain_taxid=630, length=4_600_000))
    u = _one(rows, records, resolver, manifest)
    assert u.eligibility == policy.ELIG_ELIGIBLE
    assert u.superkingdom == "Bacteria"


def test_eukaryote_pathogen_deferred(resolver, manifest, make_records):
    # A eukaryote nuclear genome under the size cap -> 'eukaryote' (deferred),
    # not organelle.
    rows = [row("Some_eukaryote", "YY000001", 270544)]
    records = make_records(
        make_record("YY000001", strain_taxid=270544, length=9_000_000, genome="chromosome"),
    )
    u = _one(rows, records, resolver, manifest)
    assert u.eligibility == policy.ELIG_EUKARYOTE


def test_nuclear_assembly_with_bundled_mito_is_eukaryote_not_organelle(
    resolver, manifest, make_records
):
    # A eukaryote nuclear assembly (many chromosomes) that also bundles its
    # mitochondrion must be classified 'eukaryote', NOT excluded as organelle --
    # the guard fires only when EVERY replicon is non-nuclear.
    rows = [
        row("Cryptococcus_neoformans", "CP003820", 270544),
        row("Cryptococcus_neoformans", "CP003834", 270544),  # the mito replicon
    ]
    records = make_records(
        make_record("CP003820", strain_taxid=270544, length=2_291_499, genome="chromosome"),
        make_record("CP003834", strain_taxid=270544, length=24_919, genome="mitochondrion"),
    )
    u = _one(rows, records, resolver, manifest)
    assert u.n_replicons == 2
    assert u.eligibility == policy.ELIG_EUKARYOTE


def test_all_organelle_multireplicon_is_excluded(resolver, manifest, make_records):
    # If EVERY replicon is an organelle, it is an organelle submission.
    rows = [
        row("Some_plant", "PL000001", 270544),
        row("Some_plant", "PL000002", 270544),
    ]
    records = make_records(
        make_record("PL000001", strain_taxid=270544, length=300_000, genome="mitochondrion"),
        make_record("PL000002", strain_taxid=270544, length=150_000, genome="chloroplast"),
    )
    u = _one(rows, records, resolver, manifest)
    assert u.eligibility == policy.ELIG_ORGANELLE


def test_genus_level_falls_back_to_organism_species(resolver, manifest, make_records):
    # Ruminococcus_sp. declares genus taxid 1263 (no species). The genome's own
    # organism taxid (R. bovis) resolves to a species -> eligible, flagged.
    rows = [row("Ruminococcus_sp.", "CP039381", 1263)]
    records = make_records(
        make_record("CP039381", strain_taxid=2564099, length=2_440_231, genome="chromosome"),
    )
    u = _one(rows, records, resolver, manifest)
    assert u.eligibility == policy.ELIG_ELIGIBLE
    assert u.species_taxid == 2564099
    assert "species_from_organism" in u.flags
