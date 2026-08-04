# -*- coding: utf-8 -*-
"""Regression tests for issues found by the adversarial review of Tier 0."""

import json

import pytest

from conftest import make_record, row
from ncbi import NuccoreRecord, NuccoreResolver, is_real_assembly
import policy


# --- real GCF/GCA assembly grouping (the PRODUCTION path, not the proxy) ------
def test_grouping_by_real_assembly_merges_and_separates(resolver, manifest, make_records):
    # Two replicons that share a real assembly -> one unit; a third genome of the
    # SAME species with a DIFFERENT assembly -> its own unit. This exercises the
    # ground-truth path (records carry assembly_acc), unlike the proxy fallback.
    rows = [
        row("Aliivibrio_fischeri", "CP000020", 668),
        row("Aliivibrio_fischeri", "CP000021", 668),
        row("Aliivibrio_fischeri", "ZZ999999", 668),  # different strain, diff assembly
    ]
    records = make_records(
        make_record("CP000020", 312309, 2_897_536, assembly_acc="GCF_000011805.1"),
        make_record("CP000021", 312309, 1_330_333, assembly_acc="GCF_000011805.1"),
        make_record("ZZ999999", 668, 4_200_000, assembly_acc="GCF_000099999.1"),
    )
    units = policy.build_units(rows, resolver, records, manifest)
    assert len(units) == 2
    merged = [u for u in units if u.parent_assembly == "GCF_000011805.1"][0]
    assert merged.accession_set == ["CP000020", "CP000021"]
    assert "grouped_by_proxy" not in merged.flags


# --- fallback must NEVER merge distinct genomes sharing a species-level taxid --
def test_proxy_fallback_does_not_merge_distinct_lineages(resolver, manifest, make_records):
    # Two Salmonella lineages, DIFFERENT Kalamari names/taxids, both reporting the
    # SPECIES organism taxid 28901, and NO real assembly resolved. They must stay
    # separate (the strain-taxid proxy is scoped to the Kalamari identity).
    rows = [
        row("Salmonella_enterica_IIa", "CP053411", 9000010),
        row("Salmonella_enterica_IV", "CP053579", 59205),
    ]
    records = make_records(
        make_record("CP053411", 28901, 4_800_000, assembly_acc=None),
        make_record("CP053579", 28901, 4_900_000, assembly_acc=None),
    )
    units = policy.build_units(rows, resolver, records, manifest)
    assert len(units) == 2  # NOT merged into a chimera


# --- virus + segment eligibility (previously untested) -----------------------
def test_virus_is_excluded(resolver, manifest, make_records):
    rows = [row("Some_virus", "NC_000001", 10239)]  # taxid under Viruses
    records = make_records(make_record("NC_000001", 10239, 30000, genome="genomic"))
    units = policy.build_units(rows, resolver, records, manifest)
    assert units[0].eligibility == policy.ELIG_VIRUS


def test_segmented_virus_stays_one_virus_unit(resolver, manifest, make_records):
    # A multi-segment virus sharing one assembly is one unit, classified virus
    # (not organelle) -- virus classification wins over the segment rule.
    rows = [
        row("Segmented_virus", "NC_00A", 10239),
        row("Segmented_virus", "NC_00B", 10239),
    ]
    records = make_records(
        make_record("NC_00A", 10239, 1871, genome="genomic", moltype="segment", assembly_acc="GCF_000000011.1"),
        make_record("NC_00B", 10239, 6562, genome="genomic", moltype="segment", assembly_acc="GCF_000000011.1"),
    )
    units = policy.build_units(rows, resolver, records, manifest)
    assert len(units) == 1
    assert units[0].eligibility == policy.ELIG_VIRUS


# --- is_real_assembly rejects source INSDC accession from esummary -----------
def test_is_real_assembly_rejects_insdc_accession():
    assert is_real_assembly("GCF_000011805.1")
    assert is_real_assembly("GCA_013138165.1")
    assert not is_real_assembly("AF291702")   # source INSDC, not an assembly
    assert not is_real_assembly("PP230457")
    assert not is_real_assembly(None)


def test_insdc_assemblyacc_does_not_split_segments(resolver, manifest, make_records):
    # Emulate the Andes hantavirus quirk: each segment's esummary assemblyacc is a
    # source INSDC accession (AF2917xx), which must be ignored so they group by
    # the shared organism taxid instead of splitting into one-unit-per-segment.
    rows = [
        row("Orthohantavirus", "NC_003466", 10239),
        row("Orthohantavirus", "NC_003467", 10239),
        row("Orthohantavirus", "NC_003468", 10239),
    ]
    # is_real_assembly filtering happens at fetch time; here records already have
    # assembly_acc=None (filtered) and share the organism taxid.
    records = make_records(
        make_record("NC_003466", 1980456, 1871, genome="genomic", assembly_acc=None),
        make_record("NC_003467", 1980456, 3671, genome="genomic", assembly_acc=None),
        make_record("NC_003468", 1980456, 6562, genome="genomic", assembly_acc=None),
    )
    units = policy.build_units(rows, resolver, records, manifest)
    assert len(units) == 1
    assert units[0].n_replicons == 3


# --- injectable assembly resolver (symmetry with fetch_fn) -------------------
def test_resolve_assemblies_is_injectable(tmp_path):
    calls = {}

    def fake_assembly_fn(recs, log=lambda m: None):
        calls["n"] = len(recs)
        for a, r in recs.items():
            r.assembly_acc = "GCF_000000001.1"

    def fake_fetch(accs):
        return {a: NuccoreRecord(a, strain_taxid=42) for a in accs}

    r = NuccoreResolver(cache_path=None, fetch_fn=fake_fetch, assembly_fn=fake_assembly_fn)
    records = r.resolve_all(["CP000001"])
    r.resolve_assemblies(records)
    assert calls["n"] == 1
    assert records["CP000001"].assembly_acc == "GCF_000000001.1"


# --- taxid -> WRONG species is detected and flagged (mislabel surfacing) ------
def test_taxid_organism_mismatch_is_flagged(resolver, manifest, make_records):
    # Kalamari row declares taxid 630 (Yersinia enterocolitica) but the deposited
    # genome's NCBI organism is taxid 632 (Yersinia pestis): a mislabel. It must
    # be flagged (declared species kept, organism surfaced in the note).
    rows = [row("Yersinia_enterocolitica", "CP999999", 630)]
    records = make_records(make_record("CP999999", strain_taxid=632, length=4_600_000))
    units = policy.build_units(rows, resolver, records, manifest)
    u = units[0]
    assert u.species_taxid == 630  # declared species kept (faithful snapshot)
    assert "taxid_organism_mismatch" in u.flags
    assert "632" in u.note


def test_no_mismatch_flag_when_organism_agrees(resolver, manifest, make_records):
    rows = [row("Yersinia_enterocolitica", "CP888888", 630)]
    records = make_records(make_record("CP888888", strain_taxid=630, length=4_600_000))
    u = policy.build_units(rows, resolver, records, manifest)[0]
    assert "taxid_organism_mismatch" not in u.flags


# --- organelle guard must not fire when a replicon record is MISSING ----------
def test_all_present_guard_blocks_organelle_when_a_record_is_missing(resolver, make_records):
    # Directly exercise classify_eligibility on a unit that MIXES a present mito
    # record with a missing nuclear-chromosome record. The all_present guard must
    # stop it being called organelle-only (the missing record is unknown, not
    # organelle) -- removing the guard would make this fail.
    unit = policy.GenomeUnit(
        unit_id="u", group_key="GCF_x", scientific_names=["Some_eukaryote"],
        kalamari_taxids=[270544], accession_set=["MITO1", "NUCLEAR1"],
        parent_assembly="GCF_x", superkingdom="Eukaryota", total_length=2_000_000,
    )
    records = make_records(
        make_record("MITO1", strain_taxid=270544, length=20000, genome="mitochondrion"),
        # NUCLEAR1 deliberately absent
    )
    policy.classify_eligibility(unit, records)
    assert unit.eligibility != policy.ELIG_ORGANELLE  # not organelle-only: a record is missing
    assert unit.eligibility == policy.ELIG_EUKARYOTE

    # sanity: when ALL replicons are present and organelle, it IS organelle-only
    unit2 = policy.GenomeUnit(
        unit_id="u2", group_key="GCF_y", scientific_names=["Some_eukaryote"],
        kalamari_taxids=[270544], accession_set=["MITO1"], parent_assembly="GCF_y",
        superkingdom="Eukaryota", total_length=20000,
    )
    policy.classify_eligibility(unit2, records)
    assert unit2.eligibility == policy.ELIG_ORGANELLE


# --- shared accession that is part of a real multi-replicon assembly must NOT split
def test_shared_accession_with_real_assembly_does_not_split(resolver, manifest, make_records):
    # A(chromosome) is shared across two subspecies labels; B(plasmid) is a sibling
    # replicon of the same real assembly. They must form ONE 2-replicon unit, not a
    # split with an orphaned B -- the sharedacc override must not clobber the assembly.
    rows = [
        row("Species_subspA", "ACC_A", 630),
        row("Species_subspB", "ACC_A", 630),   # ACC_A shared across two labels
        row("Species_subspA", "ACC_B", 630),   # sibling replicon, same assembly
    ]
    records = make_records(
        make_record("ACC_A", strain_taxid=630, length=4_000_000, assembly_acc="GCF_000000001.1"),
        make_record("ACC_B", strain_taxid=630, length=100_000, genome="plasmid",
                    assembly_acc="GCF_000000001.1"),
    )
    units = policy.build_units(rows, resolver, records, manifest)
    assert len(units) == 1
    u = units[0]
    assert u.accession_set == ["ACC_A", "ACC_B"]
    assert u.n_replicons == 2
    assert "shared_accession" in u.flags
    assert "label_conflict" in u.flags


def test_genome_group_ignores_insdc_assemblyacc():
    from ncbi import NuccoreRecord
    rec = NuccoreRecord("NC_003466", strain_taxid=1980456, assembly_acc="AF291702")
    assert rec.genome_group() == "strain:1980456"  # INSDC accession ignored


# --- elink assembly-resolution HTTP path (mocked _http_json) -----------------
def test_entrez_resolve_assemblies_maps_and_splits(monkeypatch):
    import ncbi

    # canned: nuccore uid "10" and "11" -> assembly uid "900"; uid "12" -> none.
    # A chunk containing the poison uid "99" fails unless it is alone (forces the
    # adaptive split to exercise).
    def fake_http_json(endpoint, params, timeout=60, retries=4):
        ids = [v for k, v in params if k == "id"]
        if "elink" in endpoint:
            if "77" in ids:
                raise ConnectionError("permanent failure even when alone")
            if "99" in ids and len(ids) > 1:
                raise ConnectionError("simulated truncation")
            linksets = []
            for i in ids:
                dbs = []
                if i in ("10", "11"):
                    dbs = [{"dbto": "assembly", "links": ["900"]}]
                elif i == "99":
                    dbs = [{"dbto": "assembly", "links": ["901"]}]
                linksets.append({"ids": [i], "linksetdbs": dbs})
            return {"linksets": linksets}
        # assembly esummary
        asm_ids = ",".join(ids).split(",")
        res = {"uids": asm_ids}
        for a in asm_ids:
            res[a] = {"assemblyaccession": "GCF_000000900.1" if a == "900" else "GCF_000000901.1"}
        return {"result": res}

    monkeypatch.setattr(ncbi, "_http_json", fake_http_json)
    monkeypatch.setattr(ncbi.time, "sleep", lambda *a, **k: None)

    recs = {
        "A": ncbi.NuccoreRecord("A", strain_taxid=1, uid="10"),
        "B": ncbi.NuccoreRecord("B", strain_taxid=1, uid="11"),
        "C": ncbi.NuccoreRecord("C", strain_taxid=2, uid="12"),
        "D": ncbi.NuccoreRecord("D", strain_taxid=3, uid="99"),
        "E": ncbi.NuccoreRecord("E", strain_taxid=4, uid="77"),  # fails even alone
    }
    ncbi.entrez_resolve_assemblies(recs, chunk_size=4)  # must NOT raise
    assert recs["A"].assembly_acc == "GCF_000000900.1"
    assert recs["B"].assembly_acc == "GCF_000000900.1"
    assert recs["C"].assembly_acc is None          # no assembly link
    assert recs["D"].assembly_acc == "GCF_000000901.1"  # recovered via adaptive split
    assert recs["E"].assembly_acc is None          # permanent failure -> left unlinked, non-fatal


# --- cache robustness: tolerate unknown keys, atomic write -------------------
def test_cache_tolerates_unknown_keys(tmp_path):
    cache = tmp_path / "c.json"
    cache.write_text(json.dumps({
        "schema": 99,
        "records": {"CP1": {"accession": "CP1", "strain_taxid": 5,
                            "future_field": "ignored", "another": 1}},
    }))
    r = NuccoreResolver(cache_path=str(cache))
    assert r.records["CP1"].strain_taxid == 5  # loaded, unknown keys dropped


def test_cache_write_is_atomic_and_roundtrips(tmp_path):
    cache = tmp_path / "c.json"
    r = NuccoreResolver(cache_path=str(cache), fetch_fn=lambda a: {x: NuccoreRecord(x, strain_taxid=7) for x in a})
    r.resolve_all(["CP9"])
    r.save()
    # no leftover temp files, valid JSON, reloads
    assert not list(tmp_path.glob("*.tmp*"))
    blob = json.loads(cache.read_text())
    assert blob["schema"] == 1 and "CP9" in blob["records"]


def test_cache_load_backfills_missing_accession_key(tmp_path):
    cache = tmp_path / "c.json"
    cache.write_text(json.dumps({"schema": 1,
                                 "records": {"CP7": {"strain_taxid": 9}}}))  # no 'accession'
    r = NuccoreResolver(cache_path=str(cache))
    assert r.records["CP7"].accession == "CP7"  # backfilled from the dict key
    assert r.records["CP7"].strain_taxid == 9


# --- grouped_by_proxy flag coverage ------------------------------------------
def test_grouped_by_proxy_flag_set_and_unset(resolver, manifest, make_records):
    # Two replicons of one genome grouped WITHOUT a real assembly (strain proxy)
    # -> grouped_by_proxy. The same pair WITH a real assembly -> no such flag.
    rows = [row("Aliivibrio_fischeri", "R1", 668), row("Aliivibrio_fischeri", "R2", 668)]
    proxy = make_records(
        make_record("R1", strain_taxid=668, length=3_000_000, assembly_acc=None),
        make_record("R2", strain_taxid=668, length=1_000_000, assembly_acc=None),
    )
    u = policy.build_units(rows, resolver, proxy, manifest)[0]
    assert u.n_replicons == 2 and "grouped_by_proxy" in u.flags

    real = make_records(
        make_record("R1", strain_taxid=668, length=3_000_000, assembly_acc="GCF_000011805.1"),
        make_record("R2", strain_taxid=668, length=1_000_000, assembly_acc="GCF_000011805.1"),
    )
    u2 = policy.build_units(rows, resolver, real, manifest)[0]
    assert u2.n_replicons == 2 and "grouped_by_proxy" not in u2.flags


# --- unit_id disambiguation (#2 suffix on collision) -------------------------
def test_unit_id_disambiguation_suffix(resolver, manifest, make_records):
    # Two distinct genomes: same scientificName + same strain_taxid proxy but
    # different Kalamari taxids and no real assembly -> two units whose base
    # unit_id collides, so the second must get a '#2' suffix (ids stay unique).
    rows = [row("Yersinia_x", "AAA", 630), row("Yersinia_x", "BBB", 632)]
    records = make_records(
        make_record("AAA", strain_taxid=629, length=4_000_000, assembly_acc=None),
        make_record("BBB", strain_taxid=629, length=4_100_000, assembly_acc=None),
    )
    units = policy.build_units(rows, resolver, records, manifest)
    assert len(units) == 2
    ids = sorted(u.unit_id for u in units)
    assert len(set(ids)) == 2                 # unique
    assert any(i.endswith("#2") for i in ids)  # collision disambiguated


# --- mismatch: chimera / label-conflict variants -----------------------------
def test_mismatch_flags_multiple_distinct_organisms(resolver, manifest, make_records):
    # One assembly whose two replicons report DIFFERENT organism species is a
    # chimera hint -> flagged even though the declared species matches one of them.
    rows = [row("Yersinia_enterocolitica", "M1", 630), row("Yersinia_enterocolitica", "M2", 630)]
    records = make_records(
        make_record("M1", strain_taxid=630, length=4_000_000, assembly_acc="GCF_000000050.1"),
        make_record("M2", strain_taxid=632, length=100_000, assembly_acc="GCF_000000050.1"),
    )
    u = policy.build_units(rows, resolver, records, manifest)[0]
    assert "taxid_organism_mismatch" in u.flags


def test_mismatch_not_flagged_on_organelle_unit(resolver, manifest, make_records):
    # A Eukaryota organelle whose declared vs organism species differ is NOT a
    # curation signal -> no mismatch flag (only Bacteria/Archaea are flagged).
    rows = [row("Some_host", "H1", 270544)]
    records = make_records(
        make_record("H1", strain_taxid=668, length=17000, genome="mitochondrion"),
    )
    u = policy.build_units(rows, resolver, records, manifest)[0]
    assert u.eligibility == policy.ELIG_ORGANELLE
    assert "taxid_organism_mismatch" not in u.flags
