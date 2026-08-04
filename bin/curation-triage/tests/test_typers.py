# -*- coding: utf-8 -*-
"""Tier 2: the typer adapters, tested against REAL captured tool output.

Every fixture under ``fixtures/tier2/raw/`` is verbatim output from the pinned build
of that tool, run on a real Kalamari genome on 2026-08-04.  That matters more than it
sounds: the first version of these adapters was written from each tool's
documentation, and running the tools for real showed several column names and value
formats were wrong.  Each would have produced a confident wrong answer rather than an
obvious failure, so the tests below pin the shapes that were actually observed.
"""

import argparse
import copy
import os

import pytest

import markers
import tier2_confirm
import typers

RAW = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "tier2", "raw")


def raw(name: str) -> str:
    with open(os.path.join(RAW, name), encoding="utf-8") as fh:
        return fh.read()


# --------------------------------------------------------------------------- #
#  ShigEiFinder -- the E. coli pick                                            #
# --------------------------------------------------------------------------- #
def test_shigeifinder_negative_genome_reads_as_e_coli():
    out = typers.normalize_shigeifinder(raw("shigeifinder.tsv"))
    assert out["ipaH"] == "absent"                     # the raw cell is a bare '-'
    assert out["organism"] == "Escherichia coli"
    assert out["cluster"] == "" and out["serotype"] == ""


def test_the_e_coli_pick_is_confirmed_end_to_end():
    call = markers.resolve("shigella_eiec_vs_escherichia_coli",
                           typers.normalize_shigeifinder(raw("shigeifinder.tsv")),
                           "Escherichia coli", "Escherichia_coli")
    assert call.verdict == markers.VERDICT_CONFIRMED
    assert call.overrides_ani is True


def test_shigeifinder_has_no_organism_column_so_the_call_comes_from_cluster():
    """Both CLUSTER and SEROTYPE carry the literal string 'Not Shigella/EIEC'."""
    assert "Not Shigella/EIEC" in raw("shigeifinder.tsv")
    assert not any(c.lower() in ("organism", "species", "prediction")
                   for c in raw("shigeifinder.tsv").splitlines()[0].split("\t"))


# --------------------------------------------------------------------------- #
#  BTyper3 -- the Bacillus units                                               #
# --------------------------------------------------------------------------- #
def test_btyper3_revised_nomenclature_is_read_correctly():
    """3.4 reports B. mosaicus subsp. cereus, not 'Bacillus cereus'."""
    out = typers.normalize_btyper3(raw("btyper3.tsv"))
    assert out["btyper3_species"] == "mosaicus"
    assert out["btyper3_subspecies"] == "cereus"
    assert out["closest_type_strain"] == "paranthracis"
    assert out["closest_type_strain_ani"] == "97.58"


def test_btyper3_gene_cells_are_parsed_out_of_the_fraction_format():
    out = typers.normalize_btyper3(raw("btyper3.tsv"))
    assert out["pxo1_genes"] == []           # the raw cell is '0/3()'
    assert out["pxo2_genes"] == []           # the raw cell is '0/5()'


def test_the_real_bacillus_pick_is_not_anthracis():
    out = typers.normalize_btyper3(raw("btyper3.tsv"))
    assert out["anthracis_chromosome"] == "no"
    call = markers.resolve("bacillus_anthracis_vs_b_cereus_group", out,
                           "Bacillus cereus", "Bacillus_cereus")
    # BTyper3's taxon call includes 'cereus', so the declared label stands...
    assert call.verdict == markers.VERDICT_CONFIRMED
    # ...but the nearest type strain being a different organism is recorded, because
    # that is what NCBI's taxonomy check is reacting to.
    assert "paranthracis" in "; ".join(call.reasons)


def test_paranthracis_is_never_read_as_anthracis():
    """The substring trap, on the exact string the real tool emits."""
    out = typers.normalize_btyper3(
        raw("btyper3.tsv").replace("cereus(99.96074676513672)",
                                   "paranthracis(99.96074676513672)"))
    assert out["anthracis_chromosome"] == "no"


# --------------------------------------------------------------------------- #
#  LisSero + mlst -- the Listeria units                                        #
# --------------------------------------------------------------------------- #
def test_lissero_presence_words_are_not_booleans():
    out = typers.normalize_lissero(raw("lissero.tsv"))
    assert out["serogroup"] == "Nontypeable"
    assert out["doumith_pattern"] == "PRS"        # only PRS is FULL; the rest are NONE


def test_mlst_is_positional_because_it_prints_no_header():
    out = typers.normalize_mlst_listeria(raw("mlst_listeria.tsv"))
    assert out["mlst_scheme"] == "listeria_2"
    assert out["mlst_st"] == "1779"


def test_mlst_supplies_no_lineage_so_none_is_invented():
    out = typers.normalize_mlst_listeria(raw("mlst_listeria.tsv"))
    assert out["mlst_cc"] == "" and out["broad_lineage"] == ""
    call = markers.resolve("listeria_monocytogenes_serogroup_lineage_cc",
                           out, "Listeria monocytogenes", "Listeria_monocytogenes_III")
    assert call.verdict == markers.VERDICT_INDETERMINATE
    assert "serogroup alone cannot decide" in "; ".join(call.reasons)


# --------------------------------------------------------------------------- #
#  mlst + abricate -- the Yersinia unit                                        #
# --------------------------------------------------------------------------- #
def test_the_classical_yersinia_st_is_never_used_as_a_species_call():
    """pestis is a clone inside pseudotuberculosis; the 7-locus ST cannot split them."""
    out = typers.normalize_mlst_yersinia(raw("mlst_yersinia.tsv"))
    assert out["mlst_st"] == "42"
    assert "cgmlst_species_assignment" not in out
    assert out["classical_mlst_only"] == "yes"


def test_abricate_finds_none_of_the_yersinia_discriminators():
    """VFDB + PlasmidFinder contain no pst / ypo2088 / yihN / opgG at all."""
    out = typers.normalize_abricate(raw("abricate.tsv"))
    assert out["chromosome_markers"] == []
    assert out["plasmid_markers"] == []
    # What they DO return is virulence-plasmid context the design excludes as a
    # species discriminator.
    assert sorted(out["virulence_context_markers"]) == ["lcrG", "lcrV"]


def test_the_yersinia_task_returns_undecided_rather_than_guessing():
    evidence = dict(typers.normalize_mlst_yersinia(raw("mlst_yersinia.tsv")))
    evidence.update(typers.normalize_abricate(raw("abricate.tsv")))
    call = markers.resolve("yersinia_pestis_vs_y_pseudotuberculosis", evidence,
                           "Yersinia pseudotuberculosis", "Yersinia_pseudotuberculosis")
    assert call.verdict == markers.VERDICT_INDETERMINATE


# --------------------------------------------------------------------------- #
#  SISTR + SeqSero2 -- the Salmonella units                                    #
# --------------------------------------------------------------------------- #
def test_sistr_supplies_the_subspecies():
    out = typers.normalize_sistr(raw("sistr.tab"))
    assert out["subspecies"] == "arizonae"
    assert out["sistr_serovar"] == "IIIa 41:z4,z23:-"
    assert out["sistr_qc_status"] == "PASS"


def test_seqsero2_states_the_subspecies_in_prose_as_a_cross_check():
    out = typers.normalize_seqsero2(raw("seqsero2.tsv"))
    assert out["seqsero2_serovar"] == "IIIa 41:z4,z23:-"
    assert "arizonae" in out["seqsero2_identification"]


def test_the_real_salmonella_iiia_unit_is_confirmed():
    evidence = dict(typers.normalize_sistr(raw("sistr.tab")))
    evidence.update(typers.normalize_seqsero2(raw("seqsero2.tsv")))
    call = markers.resolve("salmonella_subspecies_and_serovar", evidence,
                           "Salmonella enterica", "Salmonella_enterica_IIIa")
    assert call.verdict == markers.VERDICT_CONFIRMED
    assert call.fields["serovar_agreement"] == "agree"


# --------------------------------------------------------------------------- #
#  AMRFinderPlus -- the C. botulinum units                                     #
# --------------------------------------------------------------------------- #
def test_amrfinder_column_is_element_symbol_and_symbols_use_an_underscore():
    out = typers.normalize_amrfinderplus(raw("amrfinderplus.tsv"))
    assert out["bont_types"] == ["E"]
    assert out["bont_subtype"] == "E3"          # raw symbol is 'bont_E3'


def test_the_group_is_still_not_inferred_from_the_real_toxin_call():
    out = typers.normalize_amrfinderplus(raw("amrfinderplus.tsv"))
    call = markers.resolve("clostridium_botulinum_group_and_bont_type", out,
                           "Clostridium botulinum", "Clostridium_botulinum_groupI")
    assert call.fields["bont_gene_type"] == "E"
    assert call.verdict == markers.VERDICT_INDETERMINATE


# --------------------------------------------------------------------------- #
#  ECTyper                                                                     #
# --------------------------------------------------------------------------- #
def test_ectyper_reports_its_own_database_version_per_row():
    out = typers.normalize_ectyper(raw("ectyper.tsv"))
    assert out["ectyper_species"] == "Escherichia coli"
    assert out["serotype"] == "O118/O151:H16"
    assert out["ectyper_db_version"].startswith("v1.0")


# --------------------------------------------------------------------------- #
#  Robustness of the adapter layer                                             #
# --------------------------------------------------------------------------- #
def test_missing_output_is_empty_not_an_exception():
    for tool, fn in typers.NORMALIZERS.items():
        assert isinstance(fn(""), dict), tool


def test_a_renamed_column_yields_no_field_rather_than_a_wrong_one():
    out = typers.normalize_lissero("ID\tSOMETHING_ELSE\nx\t4b\n")
    assert out["serogroup"] == ""


def test_every_task_maps_to_installed_tools():
    for task, tools in typers.TASK_TOOLS.items():
        for tool in tools:
            assert tool in typers.NORMALIZERS, "%s: no adapter for %s" % (task, tool)
            assert typers.typer_command(tool, "g.fna", "out"), tool


def test_each_tool_gets_its_own_env_on_path():
    """LisSero needs makeblastdb and mlst needs its Perl libs, both from their env."""
    env = typers.tool_env(typers.MLST_LISTERIA, "/envs")
    assert env is None or env["PATH"].startswith("/envs/mlst/bin")


def test_tools_that_write_a_file_are_not_read_from_stdout():
    for tool in (typers.BTYPER3, typers.SISTR, typers.SEQSERO2, typers.ECTYPER,
                 typers.AMRFINDERPLUS):
        assert typers.output_location(tool, "/g/GCF_1.1.fna", "/out") != typers.STDOUT
    for tool in (typers.LISSERO, typers.SHIGEIFINDER, typers.ABRICATE,
                 typers.MLST_LISTERIA):
        assert typers.output_location(tool, "/g/GCF_1.1.fna", "/out") == typers.STDOUT


def test_a_failing_tool_is_recorded_as_a_miss_and_the_others_still_run():
    class Store:
        def ensure(self, acc):
            return "/tmp/%s.fna" % acc

    def runner(tool, fasta, outdir):
        if tool == typers.LISSERO:
            raise RuntimeError("lissero: could not find makeblastdb")
        return raw("mlst_listeria.tsv")

    cache = typers.TyperCache(runner=runner)
    cache.build([(typers.LISSERO, "GCF_1.1"), (typers.MLST_LISTERIA, "GCF_1.1")],
                Store(), "/tmp")
    assert cache.misses == ["lissero|GCF_1.1"]
    assert cache.get(typers.MLST_LISTERIA, "GCF_1.1")["mlst_st"] == "1779"


def test_an_adapter_that_yields_nothing_counts_as_a_failure():
    """An empty result must not be cached as a clean 'no findings'."""
    class Store:
        def ensure(self, acc):
            return "/tmp/%s.fna" % acc

    cache = typers.TyperCache(runner=lambda t, f, o: "totally\tunexpected\nshape\there\n")
    cache.build([(typers.LISSERO, "GCF_1.1")], Store(), "/tmp")
    assert cache.misses == ["lissero|GCF_1.1"]
    assert cache.get(typers.LISSERO, "GCF_1.1") is None


# --------------------------------------------------------------------------- #
#  Sharding: one runner, one environment                                       #
# --------------------------------------------------------------------------- #
class FakeStore:
    def ensure(self, accession):
        return "/tmp/%s.fna" % accession


def test_a_shard_runs_only_the_tools_it_was_given():
    """The eleven pinned environments are ~13 GB, so one runner holds one of them."""
    cache = typers.TyperCache(runner=lambda t, f, o: raw("mlst_listeria.tsv"))
    cache.build([(typers.MLST_LISTERIA, "GCF_1.1"), (typers.SISTR, "GCF_2.1")],
                FakeStore(), "/tmp", only_tools=[typers.MLST_LISTERIA])
    assert cache.get(typers.MLST_LISTERIA, "GCF_1.1")["mlst_st"] == "1779"
    assert cache.get(typers.SISTR, "GCF_2.1") is None
    assert cache.misses == []


def test_a_shard_probes_a_version_only_for_the_tools_it_ran():
    """A tool version IS the epoch: stamping '@unknown' re-baselines every T2 lead.

    A shard runner has one environment installed, so probing the other nine would
    record them all as unknown -- which reads as a tooling change and voids every
    dismissal in reviews.tsv.
    """
    cache = typers.TyperCache(runner=lambda t, f, o: raw("mlst_listeria.tsv"))
    cache.meta = {"tool_versions": {"sistr": "sistr@1.1.3", "btyper3": "btyper3@3.4.0"},
                  "amrfinder_database_version": "2026-05-15.1"}
    cache.build([(typers.MLST_LISTERIA, "GCF_1.1"), (typers.SISTR, "GCF_2.1")],
                FakeStore(), "/tmp", only_tools=[typers.MLST_LISTERIA])
    versions = cache.meta["tool_versions"]
    assert versions["sistr"] == "sistr@1.1.3"          # untouched by this shard
    assert versions["btyper3"] == "btyper3@3.4.0"
    assert typers.MLST_LISTERIA in versions            # the tool it did run
    # An AMRFinder-free shard must not drop the database pin either.
    assert cache.meta["amrfinder_database_version"] == "2026-05-15.1"


# --------------------------------------------------------------------------- #
#  Provenance: a version probe must never invent a tooling change              #
# --------------------------------------------------------------------------- #
# The committed meta, as ``src/curation-triage/tier2_typer_cache.json`` records it.
# ``shigeifinder@unknown`` is not a hole: ``shigeifinder --version`` exits 2, so
# unknown is that tool's genuine recorded value and it has to survive every guard.
COMMITTED_META = {
    "tool_versions": {"amrfinderplus": "amrfinderplus@4.2.7",
                      "mlst_listeria": "mlst_listeria@2.23.0",
                      "shigeifinder": "shigeifinder@unknown",
                      "sistr": "sistr@1.1.3"},
    "amrfinder_database_version": "2026-05-15.1",
}


def probes_that_reach_no_tool(monkeypatch):
    """Stand in for a runner whose pinned environments are absent.

    Returns the list every probe appends to, so a test can assert a probe that must
    not happen did not happen -- 'meta came out right' alone would also pass if the
    probes ran and their answers were merely discarded.
    """
    probed = []

    def tool_version(tool, env_dir=None):
        probed.append(tool)
        return "%s@unknown" % tool

    def amrfinder_database_version(env_dir=None):
        probed.append("amrfinderplus/database")
        return "unknown"

    monkeypatch.setattr(typers, "tool_version", tool_version)
    monkeypatch.setattr(typers, "amrfinder_database_version", amrfinder_database_version)
    return probed


class RefuseStore:
    """Proves no genome was fetched: every job in these tests is a cache hit."""

    def ensure(self, accession):
        raise AssertionError("no genome should be fetched: %s" % accession)


def refuse_to_run(tool, fasta, outdir):
    raise AssertionError("no typer should run: %s" % tool)


def test_a_build_in_which_nothing_ran_leaves_the_committed_provenance_untouched(monkeypatch):
    """The documented refresh, on a machine whose typer environments are absent.

    Every planned job is already cached, so nothing runs and ``misses`` stays empty.
    ``build`` used to re-probe the whole job list regardless, stamping ten tools
    '@unknown' and database 'unknown' over the committed pins, then exiting 0 with no
    warning at all.  A tool version IS the epoch, so that silently re-baselined every
    Tier-2 lead, voided every reviews.tsv dismissal, and wrote the false provenance
    into append-only ledger records that cannot afterwards be corrected.
    """
    probed = probes_that_reach_no_tool(monkeypatch)
    cache = typers.TyperCache(runner=refuse_to_run, env_dir="no-such-typer-env-dir")
    cache.meta = copy.deepcopy(COMMITTED_META)
    cache.results = {typers.TyperCache.key(tool, "GCF_1.1"): {"mlst_st": "1779"}
                     for tool in typers.ALL_TOOLS}

    cache.build([(tool, "GCF_1.1") for tool in typers.ALL_TOOLS], RefuseStore(), "/tmp")

    assert probed == []                 # nothing ran, so nothing was re-probed
    assert cache.misses == []
    assert cache.meta["tool_versions"] == COMMITTED_META["tool_versions"]
    assert cache.meta["amrfinder_database_version"] == "2026-05-15.1"


def test_a_failed_probe_never_overwrites_a_version_that_is_already_on_record(monkeypatch):
    """'@unknown' means the binary was not found, which is a missing tool, not a new one."""
    probes_that_reach_no_tool(monkeypatch)
    cache = typers.TyperCache(runner=lambda t, f, o: raw("mlst_listeria.tsv"))
    cache.meta = copy.deepcopy(COMMITTED_META)
    cache.build([(typers.MLST_LISTERIA, "GCF_9.9")], FakeStore(), "/tmp")
    # The tool really did run and really did produce a new result ...
    assert cache.get(typers.MLST_LISTERIA, "GCF_9.9")["mlst_st"] == "1779"
    # ... and its recorded epoch still survived the failed probe.
    assert cache.meta["tool_versions"][typers.MLST_LISTERIA] == "mlst_listeria@2.23.0"


def test_the_amrfinder_database_pin_survives_a_probe_that_cannot_reach_the_tool(monkeypatch):
    """The dated DB version is the only reproducibility anchor AMRFinderPlus gives us."""
    probes_that_reach_no_tool(monkeypatch)
    cache = typers.TyperCache(runner=lambda t, f, o: raw("amrfinderplus.tsv"))
    cache.meta = copy.deepcopy(COMMITTED_META)
    cache.build([(typers.AMRFINDERPLUS, "GCF_9.9")], FakeStore(), "/tmp")
    assert cache.get(typers.AMRFINDERPLUS, "GCF_9.9")["bont_types"]
    assert cache.meta["amrfinder_database_version"] == "2026-05-15.1"


def test_shigeifinder_reports_no_version_so_unknown_is_recorded_as_its_real_one(monkeypatch):
    """The deliberate exception: ``shigeifinder --version`` exits 2.

    The guard refuses to REPLACE a recorded version with unknown.  It must not refuse
    to record one, or the one tool that legitimately reports no version would never
    get an epoch at all.
    """
    probes_that_reach_no_tool(monkeypatch)
    cache = typers.TyperCache(runner=lambda t, f, o: raw("shigeifinder.tsv"))
    cache.build([(typers.SHIGEIFINDER, "GCF_1.1")], FakeStore(), "/tmp")
    assert cache.meta["tool_versions"][typers.SHIGEIFINDER] == "shigeifinder@unknown"


def test_a_refresh_that_cannot_find_its_typer_environments_stops_before_it_runs(tmp_path):
    """The loud alternative to a cache quietly refreshed with absent calls.

    This lives with the typers rather than with the Tier-2 CLI tests because what it
    pins is the typer preflight wiring: ``tool_binary`` falls back to a bare-name PATH
    lookup, so without this gate a mistyped ``--typer-env-dir`` records every job as a
    miss and still exits 0.
    """
    plan = tier2_confirm.Plan()
    plan.typer_jobs = [(typers.SISTR, "GCF_1.1")]
    cache = typers.TyperCache()
    args = argparse.Namespace(typer_env_dir=str(tmp_path / "no-envs"),
                              allow_unpinned=False)

    with pytest.raises(SystemExit) as exc:
        # ani_cache is unused: the plan has no ANI pair to look up.
        tier2_confirm.check_typer_environments(plan, cache, None, args)
    assert "preflight failed" in str(exc.value)
    assert "sistr" in str(exc.value)

    # A tool whose jobs are ALL cached runs nothing, so its 3-GB environment is not
    # needed to replay the committed calls and the gate must stay out of the way.
    cache.results = {typers.TyperCache.key(typers.SISTR, "GCF_1.1"): {"subspecies": "IIIa"}}
    assert tier2_confirm.check_typer_environments(plan, cache, None, args) is None


def test_the_preflight_names_every_tool_that_is_not_in_a_pinned_environment(tmp_path):
    """tool_binary() falls back to the bare name, which is a silent unpinned run."""
    problems = typers.preflight([typers.ABRICATE, typers.SISTR], str(tmp_path))
    assert len(problems) == 2
    assert all("abricate" in p or "sistr" in p for p in problems)
    (tmp_path / "abricate" / "bin").mkdir(parents=True)
    binary = tmp_path / "abricate" / "bin" / "abricate"
    binary.write_text("#!/bin/sh\n", encoding="utf-8")
    assert typers.preflight([typers.ABRICATE], str(tmp_path)) == \
        ["abricate: %s exists but is not executable" % binary]
    binary.chmod(0o755)
    assert typers.preflight([typers.ABRICATE], str(tmp_path)) == []


def test_the_amrfinder_database_comes_from_the_pinned_env_not_the_calling_shell(tmp_path,
                                                                                monkeypatch):
    """The committed 2026-05-15.1 pin was read from a database OUTSIDE the pinned env.

    ``tool_env`` used to prepend PATH and nothing else, so amrfinder resolved its
    database through whatever CONDA_PREFIX the interactive shell had exported -- a
    different installation, and a pin nobody could reproduce.
    """
    monkeypatch.setenv("CONDA_PREFIX", "/somewhere/else")
    monkeypatch.setenv("AMRFINDER_DB", "/somewhere/else/share/amrfinderplus/data/latest")
    (tmp_path / "amrfinderplus" / "bin").mkdir(parents=True)

    env = typers.tool_env(typers.AMRFINDERPLUS, str(tmp_path))
    assert env["CONDA_PREFIX"] == str(tmp_path / "amrfinderplus")
    assert "AMRFINDER_DB" not in env      # no in-env database: fail loudly, never guess

    db = tmp_path / "amrfinderplus" / "share" / "amrfinderplus" / "data" / "latest"
    db.mkdir(parents=True)
    assert typers.tool_env(typers.AMRFINDERPLUS, str(tmp_path))["AMRFINDER_DB"] == str(db)


def test_two_mlst_tasks_still_share_one_environment(tmp_path):
    (tmp_path / "mlst" / "bin").mkdir(parents=True)
    for tool in (typers.MLST_LISTERIA, typers.MLST_YERSINIA):
        env = typers.tool_env(tool, str(tmp_path))
        assert env["CONDA_PREFIX"] == str(tmp_path / "mlst")


# --------------------------------------------------------------------------- #
#  ECTyper as background context                                               #
# --------------------------------------------------------------------------- #
def test_ectyper_fields_are_attached_under_their_own_names():
    """ECTyper informs the Shigella/EIEC question but never answers it."""
    cache = typers.TyperCache()
    cache.results = {"ectyper|GCF_1.1": typers.normalize_ectyper(raw("ectyper.tsv"))}
    context = typers.context_for_task("shigella_eiec_vs_escherichia_coli", "GCF_1.1", cache)
    assert context["ectyper_serotype"] == "O118/O151:H16"
    assert context["ectyper_species"] == "Escherichia coli"   # already prefixed, not twice
    assert "ectyper_ectyper_species" not in context
    assert "context" in context["ectyper_role"]


def test_the_context_prefix_is_what_keeps_two_serotypes_apart():
    """normalize_shigeifinder and normalize_ectyper both emit a bare 'serotype'."""
    assert "serotype" in typers.normalize_shigeifinder(raw("shigeifinder.tsv"))
    assert "serotype" in typers.normalize_ectyper(raw("ectyper.tsv"))
    cache = typers.TyperCache()
    cache.results = {"ectyper|GCF_1.1": typers.normalize_ectyper(raw("ectyper.tsv"))}
    context = typers.context_for_task("shigella_eiec_vs_escherichia_coli", "GCF_1.1", cache)
    assert "serotype" not in context


def test_a_task_with_no_cached_context_simply_has_none():
    assert typers.context_for_task("shigella_eiec_vs_escherichia_coli", "GCF_9.9",
                                   typers.TyperCache()) == {}
    assert typers.context_for_task("salmonella_subspecies_and_serovar", "GCF_1.1",
                                   typers.TyperCache()) == {}


def test_ectyper_is_never_routed_as_a_task_of_its_own():
    """markers.RULES has no ectyper rule, so a lead from it would carry no rule."""
    import markers as markers_mod
    for task_id in typers.TASK_CONTEXT_TOOLS:
        assert task_id in markers_mod.RULES
    assert "escherichia_coli_serotype_pathotype" not in markers_mod.RULES


def test_merging_two_tools_never_lets_the_later_one_hide_a_disagreement():
    cache = typers.TyperCache()
    cache.results = {
        "sistr|GCF_1.1": {"subspecies": "arizonae", "sistr_serovar": "IIIa"},
        "seqsero2|GCF_1.1": {"seqsero2_serovar": "IIIb", "subspecies": "diarizonae"},
    }
    ev = typers.evidence_for_task("salmonella_subspecies_and_serovar", "GCF_1.1", cache)
    assert ev["subspecies"] == "arizonae"      # first tool wins a populated field
    assert ev["sistr_serovar"] == "IIIa" and ev["seqsero2_serovar"] == "IIIb"
