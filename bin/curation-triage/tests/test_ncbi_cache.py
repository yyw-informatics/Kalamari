# -*- coding: utf-8 -*-
"""The NCBI layer: genome-group keys, cache roundtrip, mockable fetching."""

import pytest

from ncbi import NuccoreRecord, NuccoreResolver, strip_version


def test_strip_version():
    assert strip_version("CP000020.2") == "CP000020"
    assert strip_version("J01415") == "J01415"


def test_genome_group_prefers_assembly_then_strain_then_acc():
    assert NuccoreRecord("A", assembly_acc="GCF_1", strain_taxid=5).genome_group() == "GCF_1"
    assert NuccoreRecord("A", strain_taxid=312309).genome_group() == "strain:312309"
    assert NuccoreRecord("A").genome_group() == "acc:A"


def test_is_organelle():
    assert NuccoreRecord("A", genome="mitochondrion").is_organelle
    assert NuccoreRecord("A", genome="chloroplast").is_organelle
    assert not NuccoreRecord("A", genome="chromosome").is_organelle
    assert not NuccoreRecord("A", genome=None).is_organelle


def test_resolver_uses_injected_fetch_and_caches(tmp_path):
    calls = []

    def fake_fetch(accs):
        calls.append(list(accs))
        return {a: NuccoreRecord(a, strain_taxid=100 + i) for i, a in enumerate(accs)}

    cache = tmp_path / "cache.json"
    r = NuccoreResolver(cache_path=str(cache), fetch_fn=fake_fetch)
    recs = r.resolve_all(["CP000020.2", "CP000021"])
    assert set(recs) == {"CP000020", "CP000021"}
    assert calls == [["CP000020", "CP000021"]]
    r.save()
    assert cache.exists()

    # A fresh resolver loads the cache and does NOT fetch again.
    calls.clear()
    r2 = NuccoreResolver(cache_path=str(cache), fetch_fn=fake_fetch)
    recs2 = r2.resolve_all(["CP000020", "CP000021"])
    assert calls == []
    assert recs2["CP000020"].strain_taxid == 100


def test_offline_raises_on_cache_miss(tmp_path):
    def boom(accs):  # pragma: no cover - must never be called
        raise AssertionError("network used in offline mode")

    r = NuccoreResolver(cache_path=None, fetch_fn=boom)
    with pytest.raises(LookupError):
        r.resolve_all(["NOPE"], allow_fetch=False)
