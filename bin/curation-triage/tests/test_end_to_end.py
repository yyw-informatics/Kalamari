# -*- coding: utf-8 -*-
"""End-to-end: mocked fetch -> build_units -> write_table -> read back."""

import csv

from conftest import make_record, row
from ncbi import NuccoreResolver
import policy


def test_full_pipeline_with_mocked_network(resolver, manifest, tmp_path):
    rows = [
        # multi-replicon Aliivibrio -> one unit
        row("Aliivibrio_fischeri", "CP000020", 668),
        row("Aliivibrio_fischeri", "CP000021", 668),
        # two Yersinia strains -> two units
        row("Yersinia_enterocolitica", "AM286415", 630),
        row("Yersinia_enterocolitica", "CP002246", 630),
        # organelle
        row("Arripis_trutta", "AP006810", 270544),
        # confounded regime B
        row("Salmonella_enterica_VIII", "CP053318", 9000009),
        # shared accession
        row("Campylobacter_hyointestinalis_hyointestinalis", "CP015575", 91352),
        row("Campylobacter_hyointestinalis_lawsonii", "CP015575", 91353),
    ]
    fixtures = {
        "CP000020": make_record("CP000020", 312309, 2_897_536),
        "CP000021": make_record("CP000021", 312309, 1_330_333),
        "AM286415": make_record("AM286415", 393305, 4_615_899, genome=None),
        "CP002246": make_record("CP002246", 994476, 4_552_107),
        "AP006810": make_record("AP006810", 270544, 17244, genome="mitochondrion"),
        "CP053318": make_record("CP053318", 28901, 4_800_000),
        "CP015575": make_record("CP015575", 1031746, 1_753_385),
    }

    def fake_fetch(accs):
        return {a: fixtures[a] for a in accs if a in fixtures}

    resolver_nc = NuccoreResolver(cache_path=None, fetch_fn=fake_fetch)
    records = resolver_nc.resolve_all([r.accession for r in rows])

    units = policy.build_units(rows, resolver, records, manifest)
    # 2 (Aliivibrio=1, Yersinia=2) + Arripis + Salmonella + Campylobacter(shared=1)
    assert len(units) == 6
    policy.assert_preflight(units)

    out = tmp_path / "taxon_policy.tsv"
    policy.write_table(units, str(out))

    with open(out, newline="") as fh:
        table = list(csv.DictReader(fh, delimiter="\t"))
    # required columns present and in the right leading order
    assert list(table[0].keys())[:9] == [
        "unit_id", "scientificName", "kalamari_taxid", "species_taxid",
        "accession_set", "parent_assembly", "eligibility", "regime", "resolver",
    ]
    assert len(table) == 6

    by_name = {r["scientificName"]: r for r in table}
    ali = by_name["Aliivibrio_fischeri"]
    assert ali["accession_set"] == "CP000020,CP000021"
    assert ali["eligibility"] == "eligible"

    arripis = by_name["Arripis_trutta"]
    assert arripis["eligibility"] == "exclude:organelle"

    salm = by_name["Salmonella_enterica_VIII"]
    assert salm["regime"] == "B" and salm["resolver"] == "ANI"

    # shared accession: one row, both labels joined, flagged
    shared = [r for r in table if r["accession_set"] == "CP015575"]
    assert len(shared) == 1
    assert "shared_accession" in shared[0]["flags"]
