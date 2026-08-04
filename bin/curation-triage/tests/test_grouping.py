# -*- coding: utf-8 -*-
"""Grouping rows into genome-units by assembly, never by the Kalamari taxid."""

from conftest import make_record, row
import policy


def _build(rows, records, resolver, manifest):
    return policy.build_units(rows, resolver, records, manifest)


# --- MANDATED CASE: two strains under one taxid stay SEPARATE -----------------
def test_two_strains_one_taxid_stay_separate(resolver, manifest, make_records):
    # Yersinia enterocolitica taxid 630 has two rows that are two different
    # strains (distinct strain-level organism taxids) -> must be TWO units.
    rows = [
        row("Yersinia_enterocolitica", "AM286415", 630),
        row("Yersinia_enterocolitica", "CP002246", 630),
    ]
    records = make_records(
        make_record("AM286415", strain_taxid=393305, length=4_615_899, genome=None),
        make_record("CP002246", strain_taxid=994476, length=4_552_107),
    )
    units = _build(rows, records, resolver, manifest)
    assert len(units) == 2
    accs = sorted(u.accession_field() for u in units)
    assert accs == ["AM286415", "CP002246"]
    # Both correctly resolve to the same species, but remain distinct units.
    assert all(u.species_taxid == 630 for u in units)


# --- MANDATED CASE: multi-replicon genome concatenates to ONE unit ------------
def test_multireplicon_concatenates_to_one_unit(resolver, manifest, make_records):
    # Aliivibrio fischeri: two chromosomes of one assembly (same strain taxid).
    rows = [
        row("Aliivibrio_fischeri", "CP000020", 668),
        row("Aliivibrio_fischeri", "CP000021", 668),
    ]
    records = make_records(
        make_record("CP000020", strain_taxid=312309, length=2_897_536),
        make_record("CP000021", strain_taxid=312309, length=1_330_333),
    )
    units = _build(rows, records, resolver, manifest)
    assert len(units) == 1
    u = units[0]
    assert u.accession_set == ["CP000020", "CP000021"]
    assert u.n_replicons == 2
    assert u.total_length == 2_897_536 + 1_330_333
    assert u.species_taxid == 668
    assert u.eligibility == policy.ELIG_ELIGIBLE


# --- MANDATED CASE: shared-accession quirk is handled/flagged -----------------
def test_shared_accession_is_flagged_not_double_counted(resolver, manifest, make_records):
    # CP015575 appears on TWO rows (two Campylobacter subspecies) -- one physical
    # sequence labelled two ways. Must become ONE unit and be flagged.
    rows = [
        row("Campylobacter_hyointestinalis_hyointestinalis", "CP015575", 91352),
        row("Campylobacter_hyointestinalis_lawsonii", "CP015575", 91353),
    ]
    records = make_records(
        make_record("CP015575", strain_taxid=1031746, length=1_753_385),
    )
    units = _build(rows, records, resolver, manifest)
    assert len(units) == 1  # not double-counted
    u = units[0]
    assert "shared_accession" in u.flags
    assert "label_conflict" in u.flags
    # both labels + both Kalamari taxids surfaced
    assert set(u.kalamari_taxids) == {91352, 91353}
    assert "hyointestinalis_hyointestinalis" in u.scientific_name_field()
    assert "hyointestinalis_lawsonii" in u.scientific_name_field()


def test_distinct_species_are_distinct_units(resolver, manifest, make_records):
    # Sanity: two genuinely different genomes never merge.
    rows = [
        row("Escherichia_coli", "CP027582", 562),
        row("Salmonella_enterica_I", "AE006468", 28901),
    ]
    records = make_records(
        make_record("CP027582", strain_taxid=562, length=5_000_000),
        make_record("AE006468", strain_taxid=28901, length=4_800_000),
    )
    units = _build(rows, records, resolver, manifest)
    assert len(units) == 2
