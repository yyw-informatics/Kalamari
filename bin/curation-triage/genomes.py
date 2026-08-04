#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cached, mockable genome-FASTA layer for Tier 2 -- the first tier that needs
sequence.

Tiers 0 and 1 are metadata-only by design: they read NCBI's own published verdicts
and never download a genome.  Tier 2 cannot be: skani and every typer need the
actual sequence.  This module quarantines that new capability the same way
``ani_reports.py`` quarantined the ASSEMBLY_REPORTS download -- an **injectable
fetcher** plus a cache -- so the pure logic in ``tier2.py``/``markers.py`` and all
of its tests stay offline.

**What is and is not committed.**  A genome FASTA is megabytes and must never enter
the repo, so the FASTA cache lives in a gitignored directory keyed by accession.
What *is* committed is the small derived result -- the ANI numbers
(``skani_runner.py``) and the normalised typer output (``typers.py``).  That keeps
the offline-replay property Tiers 0 and 1 have, without committing sequence.

The type-strain lookup has the same shape.  Knowing *which* assembly is the type
strain of a species is a metadata question, so its answer is a small committed JSON
map (species taxid -> assembly accession) with an injectable resolver behind it.
When the map has no entry, Tier 2 reports "type strain unresolved" and returns an
indeterminate verdict -- it never guesses an accession.
"""

from __future__ import annotations

import json
import os
import subprocess
from typing import Callable, Dict, Iterable, List, Optional, Sequence

DEFAULT_DATASETS_BIN = "datasets"


class GenomeUnavailable(RuntimeError):
    """A genome is needed but is not cached and cannot be fetched (offline)."""


# --------------------------------------------------------------------------- #
#  FASTA cache                                                                 #
# --------------------------------------------------------------------------- #
def datasets_download_command(accession: str, dest_zip: str,
                              datasets_bin: str = DEFAULT_DATASETS_BIN) -> List[str]:
    """The ncbi-datasets-cli genome download (flags centralised in one place)."""
    return [datasets_bin, "download", "genome", "accession", accession,
            "--include", "genome", "--filename", dest_zip]


def fetch_with_datasets(accession: str, dest_fasta: str,
                        datasets_bin: str = DEFAULT_DATASETS_BIN) -> str:
    """Download one assembly's FASTA to ``dest_fasta`` via ncbi-datasets-cli.

    Kept deliberately thin: it shells out, unpacks the one file we want, and returns
    the path.  Everything interesting happens in the pure layer.
    """
    import zipfile

    dest_zip = dest_fasta + ".zip"
    cmd = datasets_download_command(accession, dest_zip, datasets_bin)
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)  # noqa: S603
    if proc.returncode != 0:
        raise GenomeUnavailable("datasets download failed for %s (rc=%s): %s"
                                % (accession, proc.returncode, proc.stderr.strip()[:400]))
    try:
        with zipfile.ZipFile(dest_zip) as zf:
            names = [n for n in zf.namelist() if n.endswith((".fna", ".fa", ".fasta"))]
            if not names:
                raise GenomeUnavailable("no FASTA inside the download for %s" % accession)
            with zf.open(names[0]) as src, open(dest_fasta, "wb") as out:
                out.write(src.read())
    finally:
        if os.path.exists(dest_zip):
            os.remove(dest_zip)
    return dest_fasta


def datasets_download_many_command(accessions: Sequence[str], dest_zip: str,
                                   datasets_bin: str = DEFAULT_DATASETS_BIN) -> List[str]:
    """One download call for many accessions (flags centralised in one place)."""
    return [datasets_bin, "download", "genome", "accession", *accessions,
            "--include", "genome", "--filename", dest_zip, "--no-progressbar"]


def fetch_many_with_datasets(accessions: Sequence[str], cache_dir: str,
                             datasets_bin: str = DEFAULT_DATASETS_BIN,
                             log=lambda m: None) -> List[str]:
    """Download many assemblies in ONE datasets call and unpack them into the cache.

    A monthly Tier-2 run touches on the order of a hundred genomes.  Asking NCBI for
    them one at a time is both slower and ruder than asking once, so the bulk path is
    the default and the per-accession fetcher stays as the fallback for a single miss.

    Returns the accessions actually written.  A requested accession NCBI did not
    return is simply absent -- the caller reports it as unresolved rather than
    substituting something else.
    """
    import zipfile

    if not accessions:
        return []
    os.makedirs(cache_dir, exist_ok=True)
    dest_zip = os.path.join(cache_dir, "_bulk.%d.zip" % os.getpid())
    cmd = datasets_download_many_command(sorted(set(accessions)), dest_zip, datasets_bin)
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)  # noqa: S603
    if proc.returncode != 0:
        if os.path.exists(dest_zip):
            os.remove(dest_zip)
        raise GenomeUnavailable("datasets bulk download failed (rc=%s): %s"
                                % (proc.returncode, proc.stderr.strip()[:400]))
    written: List[str] = []
    try:
        with zipfile.ZipFile(dest_zip) as zf:
            by_acc: Dict[str, List[str]] = {}
            for name in zf.namelist():
                if not name.endswith((".fna", ".fa", ".fasta")):
                    continue
                parts = name.split("/")
                # ncbi_dataset/data/<ACCESSION>/<file>.fna
                if "data" in parts:
                    idx = parts.index("data")
                    if idx + 2 < len(parts):
                        by_acc.setdefault(parts[idx + 1], []).append(name)
            for acc, members in sorted(by_acc.items()):
                # An assembly can ship several FASTA files (one per replicon set);
                # concatenating them keeps the unit whole, which is what ANI needs.
                dest = os.path.join(cache_dir, "%s.fna" % acc)
                tmp = dest + ".tmp.%d" % os.getpid()
                with open(tmp, "wb") as out:
                    for member in sorted(members):
                        with zf.open(member) as src:
                            out.write(src.read())
                os.replace(tmp, dest)
                written.append(acc)
    finally:
        if os.path.exists(dest_zip):
            os.remove(dest_zip)
    log("datasets: downloaded %d of %d requested genome(s)"
        % (len(written), len(set(accessions))))
    return written


class GenomeStore:
    """Accession -> local FASTA path, fetched once and reused.

    ``offline=True`` never fetches: a missing genome raises ``GenomeUnavailable``,
    which the driver turns into an "unresolved" verdict rather than a crash, so a
    monthly run degrades gracefully instead of failing on one bad accession.
    """

    def __init__(self, cache_dir: str, fetcher: Callable[[str, str], str] = None,
                 offline: bool = False, datasets_bin: str = DEFAULT_DATASETS_BIN) -> None:
        self.cache_dir = cache_dir
        self.offline = offline
        self.datasets_bin = datasets_bin
        self._fetcher = fetcher or (
            lambda acc, dest: fetch_with_datasets(acc, dest, datasets_bin))
        self.fetched: List[str] = []

    def path_for(self, accession: str) -> str:
        return os.path.join(self.cache_dir, "%s.fna" % accession)

    def have(self, accession: str) -> bool:
        return os.path.exists(self.path_for(accession))

    def ensure(self, accession: str) -> str:
        """Local FASTA path for ``accession``, fetching it if that is allowed."""
        path = self.path_for(accession)
        if os.path.exists(path):
            return path
        if self.offline:
            raise GenomeUnavailable(
                "genome %s is not in the local cache %s and --offline forbids fetching"
                % (accession, self.cache_dir))
        os.makedirs(self.cache_dir, exist_ok=True)
        self._fetcher(accession, path)
        if not os.path.exists(path):
            raise GenomeUnavailable("fetcher returned without producing %s" % path)
        self.fetched.append(accession)
        return path

    def missing(self, accessions: Iterable[str]) -> List[str]:
        return sorted({a for a in accessions if a and not self.have(a)})

    def prefetch(self, accessions: Iterable[str], log=lambda m: None) -> List[str]:
        """Fetch every not-yet-cached accession in one bulk call.

        Best-effort by design: a bulk failure is logged and the run continues, because
        each genome is retried individually by ``ensure`` and a unit whose genome never
        arrives is reported as unevaluated rather than taking down the run.
        """
        want = self.missing(accessions)
        if not want or self.offline:
            return []
        try:
            got = fetch_many_with_datasets(want, self.cache_dir, self.datasets_bin, log)
        except GenomeUnavailable as exc:
            log("bulk genome download failed, falling back per accession (%s)" % exc)
            return []
        self.fetched.extend(got)
        still = self.missing(want)
        if still:
            log("%d genome(s) were not returned by datasets: %s"
                % (len(still), ", ".join(still[:5])))
        return got


# --------------------------------------------------------------------------- #
#  Type-strain map (small, committed)                                          #
# --------------------------------------------------------------------------- #
def type_strain_query(taxid: int, datasets_bin: str = DEFAULT_DATASETS_BIN) -> List[str]:
    """The intended lookup for a species' type-strain assembly.

    ``--reference`` narrows to NCBI's designated reference/representative genome;
    the caller then keeps the record whose ``type_material`` block is populated.
    Centralised here so the flags are pinned in one place.
    """
    return [datasets_bin, "summary", "genome", "taxon", str(taxid),
            "--reference", "--as-json-lines"]


class TypeStrainCache:
    """The type-strain assembly to compare a pick against.

    Keyed by the pick's **own assembly**, with the species taxid only as a fallback.
    That ordering is not a detail: NCBI records a declared type strain per *assembly*,
    and for a species with subspecies it differs between them.  *Salmonella enterica*
    is the case that matters here -- all eleven Kalamari units declare species 28901,
    but NCBI names a different, subspecies-appropriate type assembly for each, and
    those subspecies units are precisely what Tier 2 exists to check.  Keying on the
    species taxid would compare ten of the eleven against the wrong type strain.

    Committed as a small JSON map.  An unpopulated map is a valid state: Tier 2 then
    reports the ANI check as unresolved instead of inventing a comparison genome.
    """

    def __init__(self, cache_path: Optional[str] = None,
                 resolver: Callable[[int], Optional[dict]] = None) -> None:
        self.cache_path = cache_path
        self._resolver = resolver
        self.by_assembly: Dict[str, dict] = {}
        self.by_species: Dict[int, dict] = {}
        self.meta: Dict[str, object] = {}
        if cache_path and os.path.exists(cache_path):
            with open(cache_path, encoding="utf-8") as fh:
                blob = json.load(fh)
            self.by_assembly = dict(blob.get("by_assembly") or {})
            self.by_species = {int(k): v for k, v in (blob.get("by_species") or {}).items()}
            self.meta = blob.get("meta", {})

    def get(self, assembly_base: Optional[str] = None,
            taxid: Optional[int] = None) -> Optional[dict]:
        if assembly_base and assembly_base in self.by_assembly:
            return self.by_assembly[assembly_base]
        if taxid is None:
            return None
        entry = self.by_species.get(int(taxid))
        if entry or self._resolver is None:
            return entry
        entry = self._resolver(int(taxid))
        if entry:
            self.by_species[int(taxid)] = entry
        return entry

    def accession(self, assembly_base: Optional[str] = None,
                  taxid: Optional[int] = None) -> str:
        return (self.get(assembly_base, taxid) or {}).get("accession", "")

    def covers(self, taxids: Iterable[int]) -> bool:
        return {int(t) for t in taxids if t}.issubset(self.by_species)

    def save(self) -> None:
        if not self.cache_path:
            return
        os.makedirs(os.path.dirname(os.path.abspath(self.cache_path)), exist_ok=True)
        payload = {
            "schema": 1,
            "meta": self.meta or {"status": "unpopulated",
                                  "how": " ".join(type_strain_query(0)).replace(" 0 ", " <taxid> ")},
            "by_species": {str(k): v for k, v in sorted(self.by_species.items())},
        }
        tmp = "%s.tmp.%d" % (self.cache_path, os.getpid())
        try:
            with open(tmp, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=1, sort_keys=True)
                fh.write("\n")
            os.replace(tmp, self.cache_path)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)
