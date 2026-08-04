# -*- coding: utf-8 -*-
"""Golden test: rebuild the whole table from the REAL chromosomes.tsv + committed
cache (offline) and assert the key invariants of the committed policy table.

Needs a full NCBI taxdump for species resolution; skips cleanly if none is
available (so the offline unit suite still runs everywhere, while CI / a dev box
with a real dump gets full end-to-end coverage)."""

import csv
import os

import pytest

from taxonomy import TaxonomyResolver
from ncbi import NuccoreResolver
import policy

BIN = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(os.path.dirname(BIN))
CHROMOSOMES = os.path.join(REPO, "src", "chromosomes.tsv")
CACHE = os.path.join(REPO, "src", "curation-triage", "nuccore_cache.json")
MANIFEST = os.path.join(REPO, "src", "curation-triage", "confounded_manifest.tsv")
OVERLAY = os.path.join(REPO, "src", "taxonomy", "build")


def _full_taxdump():
    """A full dump (not the bundled partial one), or None."""
    return TaxonomyResolver.find_base_taxdump(None)


requires_data = pytest.mark.skipif(
    not (os.path.exists(CACHE) and _full_taxdump()),
    reason="needs the committed nuccore cache + a full NCBI taxdump",
)


@requires_data
def test_golden_policy_table_invariants():
    rows = policy.load_chromosome_rows(CHROMOSOMES)
    manifest = policy.load_confounded_manifest(MANIFEST)
    resolver = TaxonomyResolver.from_dumps(_full_taxdump(), overlay_dirs=[OVERLAY])
    nuccore = NuccoreResolver(cache_path=CACHE)
    records = nuccore.resolve_all([r.accession for r in rows], allow_fetch=False)

    units = policy.build_units(rows, resolver, records, manifest)

    # every chromosome accession is accounted for exactly once
    seen = [a for u in units for a in u.accession_set]
    assert len(seen) == len(set(seen))
    assert set(seen) == {r.accession for r in rows}

    # preflight must pass on the real data
    policy.assert_preflight(units)

    # the confounded lineages must be SEPARATE units (no chimera merge)
    def named(prefix):
        return [u for u in units if any(n.startswith(prefix) for n in u.scientific_names)]

    salmonella_enterica = [u for u in units
                           if any(n.startswith("Salmonella_enterica") for n in u.scientific_names)]
    assert len(salmonella_enterica) == 11  # 11 lineages, each its own assembly
    assert all(u.species_taxid == 28901 for u in salmonella_enterica)

    listeria = named("Listeria_monocytogenes_")
    assert len(listeria) == 4  # lineages I-IV, distinct units
    assert all(u.regime == "B" and u.resolver == "ANI" for u in listeria)

    # the two Yersinia enterocolitica strains stay separate
    yent = [u for u in units if u.scientific_names == ["Yersinia_enterocolitica"]]
    assert len(yent) == 2

    # Arripis fish mitogenome excluded as organelle
    arripis = named("Arripis_trutta")[0]
    assert arripis.eligibility == policy.ELIG_ORGANELLE

    # CP015575 shared accession -> exactly one unit, flagged
    shared = [u for u in units if "CP015575" in u.accession_set]
    assert len(shared) == 1 and "shared_accession" in shared[0].flags


@requires_data
def test_golden_matches_committed_table():
    """The freshly built table must equal the committed taxon_policy.tsv."""
    committed = os.path.join(REPO, "src", "curation-triage", "taxon_policy.tsv")
    if not os.path.exists(committed):
        pytest.skip("no committed table")

    rows = policy.load_chromosome_rows(CHROMOSOMES)
    manifest = policy.load_confounded_manifest(MANIFEST)
    resolver = TaxonomyResolver.from_dumps(_full_taxdump(), overlay_dirs=[OVERLAY])
    records = NuccoreResolver(cache_path=CACHE).resolve_all(
        [r.accession for r in rows], allow_fetch=False
    )
    units = policy.build_units(rows, resolver, records, manifest)

    with open(committed, newline="") as fh:
        committed_rows = list(csv.DictReader(fh, delimiter="\t"))
    built_rows = [u.as_row() for u in units]
    assert len(built_rows) == len(committed_rows)

    # 1) row ORDER must match exactly (build is deterministically sorted)
    assert [b["unit_id"] for b in built_rows] == [c["unit_id"] for c in committed_rows]

    # 2) EVERY column must match per unit_id -- protects flags/note/superkingdom
    #    (the whole review surface), not just the required 9.
    by_id_committed = {r["unit_id"]: r for r in committed_rows}
    for b in built_rows:
        c = by_id_committed[b["unit_id"]]
        for col in policy.TABLE_COLUMNS:
            assert str(b[col]) == c[col], f"{b['unit_id']} col {col}: built={b[col]!r} committed={c[col]!r}"
