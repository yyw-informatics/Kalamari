# -*- coding: utf-8 -*-
"""coverage_probe.main() end-to-end, fully offline.

Exercises the WHOLE path -- stream + distil the three NCBI files (via ``--*-file``
local copies, so no network), join, flag, write the report + cache -- then re-runs
``--offline`` from the committed cache.  Tiny self-contained fixtures.
"""

import os

import coverage_probe
from test_ani_reports import ani_line, summary_line, ANI_HEADER, SUMMARY_HEADER_1, SUMMARY_HEADER_2


POLICY_HEADER = ("unit_id\tscientificName\tkalamari_taxid\tspecies_taxid\taccession_set\t"
                 "parent_assembly\teligibility\tregime\tresolver\tspecies_name\tsuperkingdom\t"
                 "total_length_bp\tn_replicons\tflags\tnote")


def _policy_row(unit, name, species_taxid, accs, parent, elig="eligible", regime="", resolver="",
                species_name="Sp."):
    return (f"{unit}\t{name}\t{species_taxid}\t{species_taxid}\t{accs}\t{parent}\t{elig}\t"
            f"{regime}\t{resolver}\t{species_name}\tBacteria\t3000000\t1\t\t")


def _write_inputs(tmp_path):
    policy = tmp_path / "taxon_policy.tsv"
    policy.write_text("\n".join([
        POLICY_HEADER,
        _policy_row("Foo~GCF_1.1", "Foo", 100, "CP1", "GCF_1.1"),
        _policy_row("Bar~GCF_2.1", "Bar", 200, "CP2", "GCF_2.1"),
        _policy_row("Baz~GCF_3.1", "Baz", 300, "CP3", "GCF_3.1", regime="A", resolver="marker"),
        _policy_row("Org~strain-9", "Org", 400, "CP4", "strain:9", elig="exclude:organelle"),
    ]) + "\n")

    chrom = tmp_path / "chromosomes.tsv"
    chrom.write_text(
        "scientificName\tnuccoreAcc\ttaxid\tparent\tsource\n"
        "Foo\tCP1\t100\t1\tNCBI-REF\n"
        "Bar\tCP2\t200\t1\tCDC SME\n"
        "Baz\tCP3\t300\t1\tNCBI-REF\n"
        "Org\tCP4\t400\t1\tNCBI-REF\n"
    )

    ani = tmp_path / "ani.txt"
    ani.write_text("".join([
        ANI_HEADER + "\n",
        ani_line(genbank="GCA_1.1", refseq="GCF_1.1", species_taxid=100, best_taxid=100),
        ani_line(genbank="GCA_2.1", refseq="GCF_2.1", species_taxid=200, best_taxid=999,
                 best_name="Wrong sp.", best_status="mismatch", tax_check="Failed"),
        ani_line(genbank="GCA_3.1", refseq="GCF_3.1", species_taxid=300, best_taxid=888,
                 best_name="Sister group sp.", best_status="mismatch", tax_check="Failed"),
    ]))

    refseq = tmp_path / "refseq.txt"
    refseq.write_text("".join([
        SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
        summary_line(acc="GCF_1.1", species_taxid=100, paired="GCA_1.1", relation="assembly from type material"),
        summary_line(acc="GCF_2.1", species_taxid=200, paired="GCA_2.1"),
        summary_line(acc="GCF_3.1", species_taxid=300, paired="GCA_3.1"),
    ]))

    genbank = tmp_path / "genbank.txt"
    genbank.write_text("".join([
        SUMMARY_HEADER_1 + "\n", SUMMARY_HEADER_2 + "\n",
        summary_line(acc="GCA_1.1", species_taxid=100, paired="GCF_1.1", source_db="genbank"),
        summary_line(acc="GCA_2.1", species_taxid=200, paired="GCF_2.1", source_db="genbank"),
        summary_line(acc="GCA_3.1", species_taxid=300, paired="GCF_3.1", source_db="genbank"),
    ]))
    return policy, chrom, ani, refseq, genbank


def _read_report(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = {}
        for line in fh:
            vals = line.rstrip("\n").split("\t")
            d = dict(zip(header, vals))
            rows[d["scientificName"]] = d
    return rows


def test_cli_builds_report_and_flags(tmp_path, capsys):
    policy, chrom, ani, refseq, genbank = _write_inputs(tmp_path)
    cache = tmp_path / "cache.json"
    out = tmp_path / "coverage_report.tsv"

    rc = coverage_probe.main([
        "--policy", str(policy), "--chromosomes", str(chrom),
        "--cache", str(cache), "--out", str(out),
        "--ani-file", str(ani), "--refseq-file", str(refseq), "--genbank-file", str(genbank),
    ])
    assert rc == 0
    assert cache.exists() and out.exists()

    rows = _read_report(out)
    assert set(rows) == {"Foo", "Bar", "Baz"}          # organelle excluded from the join
    assert rows["Foo"]["ani_match"] == "exact"
    assert rows["Foo"]["flags"] == ""                   # clean pick
    assert "tax_check_failed" in rows["Bar"]["flags"]   # actionable
    assert rows["Bar"]["source_lane"] == "SME"
    assert "confounded_ani_ambiguity" in rows["Baz"]["flags"]
    assert "tax_check_failed" not in rows["Baz"]["flags"]

    out_txt = capsys.readouterr().out
    assert "ANI report        : 3/3  (100.0%)" in out_txt
    assert "Units with >=1 actionable flag: 1/3" in out_txt   # only Bar


def test_cli_offline_reuses_cache(tmp_path):
    policy, chrom, ani, refseq, genbank = _write_inputs(tmp_path)
    cache = tmp_path / "cache.json"
    out = tmp_path / "coverage_report.tsv"
    args = ["--policy", str(policy), "--chromosomes", str(chrom),
            "--cache", str(cache), "--out", str(out)]

    # Build once (local files), then run offline from the committed cache.
    coverage_probe.main(args + ["--ani-file", str(ani), "--refseq-file", str(refseq),
                                "--genbank-file", str(genbank)])
    first = out.read_text()
    rc = coverage_probe.main(args + ["--offline"])
    assert rc == 0
    assert out.read_text() == first          # byte-identical from cache


def test_cli_offline_refuses_stale_cache_after_retaxonomy(tmp_path):
    # Same assembly bases but a changed species_taxid (a Phase-A retaxonomy) makes the
    # cached per-species presence counts stale; offline must refuse, not silently reuse.
    policy, chrom, ani, refseq, genbank = _write_inputs(tmp_path)
    cache = tmp_path / "cache.json"
    out = tmp_path / "coverage_report.tsv"
    coverage_probe.main(["--policy", str(policy), "--chromosomes", str(chrom),
                         "--cache", str(cache), "--out", str(out),
                         "--ani-file", str(ani), "--refseq-file", str(refseq),
                         "--genbank-file", str(genbank)])
    policy2 = tmp_path / "policy2.tsv"          # Foo's species 100 -> 555, same GCF_1 base
    policy2.write_text(policy.read_text().replace("\t100\t100\t", "\t555\t555\t"))
    try:
        coverage_probe.main(["--policy", str(policy2), "--chromosomes", str(chrom),
                             "--cache", str(cache), "--out", str(out), "--offline"])
        assert False, "expected SystemExit on species-stale cache"
    except SystemExit as e:
        assert "does not cover" in str(e)


def test_cli_offline_incomplete_cache_fails(tmp_path):
    policy, chrom, ani, refseq, genbank = _write_inputs(tmp_path)
    cache = tmp_path / "missing.json"        # no cache built
    out = tmp_path / "coverage_report.tsv"
    try:
        coverage_probe.main(["--policy", str(policy), "--chromosomes", str(chrom),
                             "--cache", str(cache), "--out", str(out), "--offline"])
        assert False, "expected SystemExit on incomplete offline cache"
    except SystemExit as e:
        assert "does not cover" in str(e)
