# -*- coding: utf-8 -*-
"""Preflight: HARD-FAIL when an eligible unit has no real species taxid."""

import pytest

from conftest import make_record, row
import policy


# --- MANDATED CASE: preflight fails on a deliberately unresolvable taxon ------
def test_preflight_fails_on_unresolvable_eligible_taxon(resolver, manifest, make_records):
    # 7000001 is a bacterial *strain*-level node with NO species ancestor (a
    # broken/partial-taxonomy situation). It classifies as eligible (Bacteria)
    # but has no species taxid, so the preflight must catch it.
    rows = [row("Faketobacter_unresolvabilis", "ZZ000001", 7000001)]
    records = make_records(
        # organism taxid is also unresolvable, so the fallback cannot save it
        make_record("ZZ000001", strain_taxid=7000001, length=3_000_000),
    )
    units = policy.build_units(rows, resolver, records, manifest)
    assert len(units) == 1
    u = units[0]
    assert u.eligibility == policy.ELIG_ELIGIBLE
    assert u.species_taxid is None

    offenders = policy.run_preflight(units)
    assert offenders == [u]
    with pytest.raises(policy.PreflightError):
        policy.assert_preflight(units)


def test_preflight_fails_on_no_rank_bacterial_leaf(resolver, manifest, make_records):
    # A 'no rank' bacterial leaf with no species ancestor must HARD-FAIL, not
    # soft-demote -- "no rank" is the ambiguous rank NCBI gives most strains, so
    # treating it as genus-level would let a strain key silently pass.
    rows = [row("Faketobacter_norankensis", "ZZ000002", 7000002)]
    records = make_records(make_record("ZZ000002", strain_taxid=7000002, length=3_000_000))
    units = policy.build_units(rows, resolver, records, manifest)
    assert units[0].eligibility == policy.ELIG_ELIGIBLE
    assert units[0].species_taxid is None
    with pytest.raises(policy.PreflightError):
        policy.assert_preflight(units)


def test_preflight_passes_on_clean_bacteria(resolver, manifest, make_records):
    rows = [row("Aliivibrio_fischeri", "CP000020", 668)]
    records = make_records(make_record("CP000020", strain_taxid=312309, length=2_897_536))
    units = policy.build_units(rows, resolver, records, manifest)
    assert policy.run_preflight(units) == []
    policy.assert_preflight(units)  # does not raise


def test_genus_level_entry_does_not_fail_preflight(resolver, manifest, make_records):
    # A genuine genus-level 'sp.' entry where even the organism taxid has no
    # species is demoted to review:genus-level (surfaced), NOT eligible -- so it
    # must not crash the whole build.
    rows = [row("Ruminococcus_sp.", "CP039381", 1263)]
    records = make_records(
        # organism taxid is the genus itself -> still no species
        make_record("CP039381", strain_taxid=1263, length=2_440_231),
    )
    units = policy.build_units(rows, resolver, records, manifest)
    u = units[0]
    assert u.eligibility == policy.ELIG_GENUS
    assert u.species_taxid is None
    assert policy.run_preflight(units) == []  # not an offender
    policy.assert_preflight(units)  # does not raise


def test_excluded_units_never_trigger_preflight(resolver, manifest, make_records):
    # An organelle/virus/eukaryote unit with no species taxid is fine.
    rows = [row("Arripis_trutta", "AP006810", 270544)]
    records = make_records(
        make_record("AP006810", strain_taxid=270544, length=17244, genome="mitochondrion"),
    )
    units = policy.build_units(rows, resolver, records, manifest)
    assert units[0].eligibility == policy.ELIG_ORGANELLE
    policy.assert_preflight(units)  # does not raise
