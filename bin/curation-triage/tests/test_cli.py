# -*- coding: utf-8 -*-
"""build_policy.main() CLI orchestration + preflight-driven exit code.

Self-contained tiny fixtures (a 2-row chromosomes.tsv, a tiny taxdump, a tiny
cache) so these run deterministically offline without the real data or network.
"""

import json
import os

import build_policy


def _write_taxdump(d):
    d.mkdir()
    (d / "nodes.dmp").write_text(
        "1\t|\t1\t|\tno rank\t|\n"
        "131567\t|\t1\t|\tno rank\t|\n"
        "2\t|\t131567\t|\tsuperkingdom\t|\n"
        "1350\t|\t2\t|\tgenus\t|\n"
        "1351\t|\t1350\t|\tspecies\t|\n"       # Enterococcus faecalis (clean)
        "7000000\t|\t2\t|\tgenus\t|\n"
        "7000001\t|\t7000000\t|\tstrain\t|\n"  # unresolvable strain (no species)
    )
    (d / "names.dmp").write_text(
        "1\t|\troot\t|\t\t|\tscientific name\t|\n"
        "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
        "1351\t|\tEnterococcus faecalis\t|\t\t|\tscientific name\t|\n"
        "7000001\t|\tFaketobacter unresolvabilis\t|\t\t|\tscientific name\t|\n"
    )


def _write_cache(path, records):
    path.write_text(json.dumps({"schema": 1, "records": records}))


# .../Kalamari-master/bin/curation-triage/tests/test_cli.py -> repo = Kalamari-master
_BIN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_REPO = os.path.dirname(os.path.dirname(_BIN))
MANIFEST = os.path.join(_REPO, "src", "curation-triage", "confounded_manifest.tsv")


def _args(chrom, taxdump, cache, out):
    return [
        "--offline",
        "--chromosomes", str(chrom),
        "--taxdump", str(taxdump),
        "--cache", str(cache),
        "--manifest", MANIFEST,
        "--out", str(out),
    ]


def test_cli_exits_zero_and_writes_table(tmp_path, capsys):
    chrom = tmp_path / "chromosomes.tsv"
    chrom.write_text(
        "scientificName\tnuccoreAcc\ttaxid\tparent\tsource\n"
        "Enterococcus_faecalis\tAE016830\t1351\t1350\tNCBI-REF\n"
    )
    taxdump = tmp_path / "taxdump"
    _write_taxdump(taxdump)
    cache = tmp_path / "cache.json"
    _write_cache(cache, {"AE016830": {"accession": "AE016830", "strain_taxid": 1351,
                                      "length": 3_000_000, "genome": "chromosome"}})
    out = tmp_path / "policy.tsv"

    rc = build_policy.main(_args(chrom, taxdump, cache, out))
    assert rc == 0
    assert out.exists()
    header = out.read_text().splitlines()[0].split("\t")
    assert header[:9] == [
        "unit_id", "scientificName", "kalamari_taxid", "species_taxid",
        "accession_set", "parent_assembly", "eligibility", "regime", "resolver",
    ]
    assert "Preflight: passed" in capsys.readouterr().out


def test_cli_exits_nonzero_on_preflight_failure(tmp_path, capsys):
    chrom = tmp_path / "chromosomes.tsv"
    chrom.write_text(
        "scientificName\tnuccoreAcc\ttaxid\tparent\tsource\n"
        "Faketobacter_x\tZZ000001\t7000001\t7000000\tTEST\n"
    )
    taxdump = tmp_path / "taxdump"
    _write_taxdump(taxdump)
    cache = tmp_path / "cache.json"
    _write_cache(cache, {"ZZ000001": {"accession": "ZZ000001", "strain_taxid": 7000001,
                                      "length": 3_000_000, "genome": "chromosome"}})
    out = tmp_path / "policy.tsv"

    rc = build_policy.main(_args(chrom, taxdump, cache, out))
    assert rc == 1
    assert "FAILED" in capsys.readouterr().out
