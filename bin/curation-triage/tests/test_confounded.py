# -*- coding: utf-8 -*-
"""Confounded-taxa tagging (regime A/marker, regime B/ANI) from the manifest."""

from conftest import make_record, row
import policy


def _one(rows, records, resolver, manifest):
    units = policy.build_units(rows, resolver, records, manifest)
    assert len(units) == 1
    return units[0]


def test_salmonella_enterica_is_regime_B_ani(resolver, manifest, make_records):
    rows = [row("Salmonella_enterica_VIII", "CP053318", 9000009)]  # synthetic subsp.
    records = make_records(make_record("CP053318", strain_taxid=28901, length=4_800_000))
    u = _one(rows, records, resolver, manifest)
    assert u.species_taxid == 28901
    assert u.regime == "B"
    assert u.resolver == "ANI"
    assert "confounded" in u.flags


def test_listeria_monocytogenes_lineage_is_regime_B(resolver, manifest, make_records):
    rows = [row("Listeria_monocytogenes_I", "CP054040", 9000000)]
    records = make_records(make_record("CP054040", strain_taxid=1639, length=3_000_000))
    u = _one(rows, records, resolver, manifest)
    assert u.regime == "B"
    assert u.resolver == "ANI"


def test_clostridium_botulinum_group_is_regime_B(resolver, manifest, make_records):
    rows = [row("Clostridium_botulinum_groupI", "CP001078", 9000005)]
    records = make_records(make_record("CP001078", strain_taxid=1491, length=3_900_000))
    u = _one(rows, records, resolver, manifest)
    assert u.regime == "B"
    assert u.resolver == "ANI"


def test_ecoli_is_regime_A_marker(resolver, manifest, make_records):
    rows = [row("Escherichia_coli", "CP027582", 562)]
    records = make_records(make_record("CP027582", strain_taxid=562, length=5_000_000))
    u = _one(rows, records, resolver, manifest)
    assert u.regime == "A"
    assert u.resolver == "marker"


def test_bacillus_cereus_group_is_regime_A_via_lineage(resolver, manifest, make_records):
    # Trigger 86661 (Bacillus cereus group) tags member species by lineage.
    rows = [row("Bacillus_cereus", "AP007209", 1396)]
    records = make_records(make_record("AP007209", strain_taxid=1396, length=5_400_000))
    u = _one(rows, records, resolver, manifest)
    assert u.regime == "A"
    assert u.resolver == "marker"


def test_yersinia_pseudotuberculosis_is_regime_A(resolver, manifest, make_records):
    rows = [row("Yersinia_pseudotuberculosis", "CP009712", 502800)]
    records = make_records(make_record("CP009712", strain_taxid=502800, length=4_700_000))
    u = _one(rows, records, resolver, manifest)
    assert u.regime == "A"
    assert u.resolver == "marker"


def test_non_confounded_bacterium_has_blank_regime(resolver, manifest, make_records):
    # Bacillus subtilis is NOT in the cereus group -> not confounded.
    rows = [row("Bacillus_subtilis", "AL009126", 1423)]
    records = make_records(make_record("AL009126", strain_taxid=1423, length=4_200_000))
    u = _one(rows, records, resolver, manifest)
    assert u.regime == ""
    assert u.resolver == ""
    assert "confounded" not in u.flags


def test_salmonella_bongori_not_confounded(resolver, manifest, make_records):
    # S. bongori is a different species from S. enterica -> not tagged.
    rows = [row("Salmonella_bongori", "FR877557", 54736)]
    records = make_records(make_record("FR877557", strain_taxid=54736, length=4_500_000))
    u = _one(rows, records, resolver, manifest)
    assert u.regime == ""
