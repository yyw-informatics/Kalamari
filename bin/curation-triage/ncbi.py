#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cached, mockable NCBI network layer for Tier 0.

This is the *only* part of ``build_policy.py`` that needs the network, and it is
deliberately quarantined here so the rest of the pipeline (taxonomy resolution,
grouping, eligibility, confounded tagging, preflight) stays pure and offline.

For each nuccore accession in ``chromosomes.tsv`` we need three facts:

    * **the parent genome** — to group replicons into one genome-unit *by
      assembly, never by the Kalamari taxid column*;
    * **the organelle / molecule type** — so the eligibility guard can catch a
      mitochondrion/plastid programmatically (a hand list already missed the
      Arripis fish mitogenome);
    * **the replicon length** — for the ">50 Mbp = organelle/eukaryote" guard.

All three come from one fast, batched Entrez ``esummary`` call.  A crucial
detail (verified 2026-07): NCBI's nuccore ``esummary`` returns the record's own
**strain-level ``taxid``**, which cleanly *merges* the replicons of one assembly
(both Aliivibrio fischeri chromosomes report strain taxid 312309) and *splits*
distinct strains sharing a Kalamari taxid (the two Yersinia enterocolitica rows
under taxid 630 report 393305 vs 994476).  So the genome-group key is:

    real assembly accession (GCF/GCA, if resolved)  ->  ``strain:<taxid>``  ->  ``acc:<accession>``

``AssemblyAcc`` is empty for most of Kalamari's old direct-submission records, so
the real GCF/GCA is fetched via ``elink`` (slow, rate-limited) and is the DEFAULT
grouping key.  The ``strain:<taxid>`` proxy is only a fallback for when no
assembly resolves: it CANNOT by itself separate distinct genomes that report the
same species-level organism taxid (the 12 Salmonella / 4 Listeria lineages all
report their species taxid), so the fallback key in ``policy.py`` is additionally
scoped to the Kalamari identity to keep those apart.

The fetch function is injectable (``fetch_fn``) and results are persisted to a
committed JSON cache, so a normal run and every unit test are fully offline.
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, asdict, fields
from typing import Callable, Dict, Iterable, List, Optional

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Values of the esummary ``genome`` field that mean "not a cellular chromosome":
# an organelle or an extrachromosomal element.  Used by the eligibility guard.
ORGANELLE_GENOME_VALUES = frozenset(
    {"mitochondrion", "chloroplast", "plastid", "apicoplast", "kinetoplast", "cyanelle"}
)


@dataclass
class NuccoreRecord:
    """The cached NCBI facts we need about one nuccore accession."""

    accession: str  # version-stripped, e.g. "CP000020"
    length: Optional[int] = None  # Slen (bp of this replicon)
    genome: Optional[str] = None  # esummary 'genome': chromosome/mitochondrion/plasmid/...
    moltype: Optional[str] = None  # esummary 'moltype'
    strain_taxid: Optional[int] = None  # the record's OWN organism taxid (strain level)
    organism: Optional[str] = None
    assembly_acc: Optional[str] = None  # real GCF_/GCA_ (filled by elink resolution)
    title: Optional[str] = None
    uid: Optional[str] = None  # nuccore UID (used to map elink linksets robustly)

    def genome_group(self) -> str:
        """Assembly-level grouping key (never the Kalamari taxid column).

        The real GCF/GCA assembly is the ground truth: it merges the replicons of
        one genome and keeps distinct genomes apart even when they share an
        organism taxid (the 5 Salmonella lineages all report species taxid 28901
        but are 5 different assemblies).  Only when no assembly is resolved (e.g.
        an organelle with no genome assembly) do we fall back to the record's own
        strain taxid, then the accession.
        """
        if is_real_assembly(self.assembly_acc):
            return self.assembly_acc  # type: ignore[return-value]
        if self.strain_taxid:
            return f"strain:{self.strain_taxid}"
        return f"acc:{self.accession}"

    @property
    def is_organelle(self) -> bool:
        return (self.genome or "").lower() in ORGANELLE_GENOME_VALUES


def strip_version(accession: str) -> str:
    """``CP000020.2`` -> ``CP000020`` (Kalamari lists unversioned accessions)."""
    return accession.split(".")[0].strip()


def is_real_assembly(acc: Optional[str]) -> bool:
    """True only for genome-assembly accessions (``GCF_``/``GCA_``).

    Guards against a nuccore ``esummary`` quirk: for RefSeq records the
    ``assemblyacc`` field can hold the *source INSDC accession* (e.g. the Andes
    hantavirus RefSeq segments report ``AF291702``/``AF291703``/``AF291704``),
    which is not a genome assembly and would wrongly split a multi-segment
    genome into one unit per segment.
    """
    return bool(acc) and acc.startswith(("GCF_", "GCA_"))


# --------------------------------------------------------------------------- #
#  Default Entrez fetcher (network).  Injectable so tests never hit the wire.  #
# --------------------------------------------------------------------------- #
def _entrez_common_params() -> Dict[str, str]:
    params: Dict[str, str] = {"tool": "kalamari-curation-triage"}
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]
    if os.environ.get("EMAIL"):
        params["email"] = os.environ["EMAIL"]
    return params


def _http_json(endpoint: str, params, timeout: int = 60, retries: int = 4) -> dict:
    """POST an E-utilities request and parse JSON, with retry/backoff.

    POST (not GET) keeps large id lists off the URL, avoiding the truncated
    responses that a multi-hundred-id GET triggers.  ``params`` is a list of
    (key, value) pairs so repeated ``id`` keys are preserved.
    """
    body = urllib.parse.urlencode(params).encode()
    last_exc = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(endpoint, data=body)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode())
        except Exception as exc:  # noqa: BLE001 - transient network/HTTP errors
            last_exc = exc
            time.sleep(min(2 ** attempt, 8))
    raise last_exc  # type: ignore[misc]


def entrez_esummary(
    accessions: List[str],
    chunk_size: int = 180,
    pause: float = 0.34,
    timeout: int = 60,
) -> Dict[str, NuccoreRecord]:
    """Batched nuccore ``esummary`` (v2.0 JSON).  Network."""
    out: Dict[str, NuccoreRecord] = {}
    for i in range(0, len(accessions), chunk_size):
        chunk = accessions[i : i + chunk_size]
        params = list(_entrez_common_params().items())
        params += [("db", "nuccore"), ("id", ",".join(chunk)),
                   ("version", "2.0"), ("retmode", "json")]
        data = _http_json(f"{EUTILS}/esummary.fcgi", params, timeout=timeout)
        result = data.get("result", {})
        for uid in result.get("uids", []):
            d = result[uid]
            acc = strip_version(d.get("accessionversion") or d.get("caption") or "")
            if not acc:
                continue
            taxid = d.get("taxid")
            asm = d.get("assemblyacc") or None
            out[acc] = NuccoreRecord(
                accession=acc,
                length=int(d["slen"]) if d.get("slen") not in (None, "") else None,
                genome=(d.get("genome") or None),
                moltype=(d.get("moltype") or None),
                strain_taxid=int(taxid) if taxid not in (None, "", 0) else None,
                organism=(d.get("organism") or None),
                assembly_acc=asm if is_real_assembly(asm) else None,
                title=(d.get("title") or None),
                uid=str(uid),
            )
        if i + chunk_size < len(accessions):
            time.sleep(pause)
    return out


def _pick_assembly(accessions: List[str]) -> Optional[str]:
    """Choose one assembly accession deterministically (prefer GCF over GCA)."""
    real = sorted(a for a in accessions if is_real_assembly(a))
    if not real:
        return None
    gcf = [a for a in real if a.startswith("GCF_")]
    return (gcf or real)[0]


def _elink_batch(uids: List[str], common, timeout: int) -> Dict[str, List[str]]:
    """One elink call: {nuccore_uid: [assembly_uid,...]} keyed by the query UID."""
    q = [("dbfrom", "nuccore"), ("db", "assembly"), ("retmode", "json")]
    q += [("id", u) for u in uids]
    q += list(common.items())
    data = _http_json(f"{EUTILS}/elink.fcgi", q, timeout=timeout)
    out: Dict[str, List[str]] = {}
    for ls in data.get("linksets", []):
        src = ls.get("ids") or []
        if not src:
            continue
        asm_uids: List[str] = []
        for db in ls.get("linksetdbs", []):
            if db.get("dbto") == "assembly":
                asm_uids.extend(str(x) for x in db.get("links", []))
        if asm_uids:
            out[str(src[0])] = asm_uids
    return out


def entrez_resolve_assemblies(
    records: Dict[str, NuccoreRecord],
    chunk_size: int = 30,
    pause: float = 0.34,
    timeout: int = 90,
    log=lambda msg: None,
) -> None:
    """Fill each record's real GCF/GCA ``assembly_acc`` via batched ``elink``.

    This is what makes grouping-by-assembly correct rather than a heuristic.  It
    maps results through the nuccore UID (not response order), so it is robust
    even if NCBI reorders or omits linksets.  NCBI truncates large elink
    responses, so a chunk that fails is split in half and retried down to a
    single id; a chunk that still fails is skipped non-fatally (those records
    stay unlinked and fall back to the strain/accession grouping key).
    """
    common = _entrez_common_params()
    uid_records = {r.uid: r for r in records.values() if r.uid}
    uids = sorted(uid_records)
    if not uids:
        return

    nuc_to_asm: Dict[str, List[str]] = {}

    def fetch(chunk: List[str]) -> None:
        if not chunk:
            return
        try:
            nuc_to_asm.update(_elink_batch(chunk, common, timeout))
            time.sleep(pause)
        except Exception as exc:  # noqa: BLE001 - truncation/throttle; split & retry
            if len(chunk) == 1:
                log(f"  elink failed for nuccore uid {chunk[0]} (left unlinked): {exc}")
                return
            mid = len(chunk) // 2
            fetch(chunk[:mid])
            fetch(chunk[mid:])

    for i in range(0, len(uids), chunk_size):
        fetch(uids[i : i + chunk_size])

    # assembly UID -> assembly accession (GCF_/GCA_)
    all_asm_uids = sorted({u for lst in nuc_to_asm.values() for u in lst})
    asm_acc: Dict[str, str] = {}
    for i in range(0, len(all_asm_uids), chunk_size):
        chunk = all_asm_uids[i : i + chunk_size]
        p = list(common.items())
        p += [("db", "assembly"), ("id", ",".join(chunk)),
              ("version", "2.0"), ("retmode", "json")]
        try:
            data = _http_json(f"{EUTILS}/esummary.fcgi", p, timeout=timeout)
        except Exception as exc:  # noqa: BLE001 - non-fatal; those stay unlinked
            log(f"  assembly esummary failed for a chunk (left unlinked): {exc}")
            continue
        res = data.get("result", {})
        for u in res.get("uids", []):
            acc = res[u].get("assemblyaccession")
            if acc:
                asm_acc[str(u)] = acc
        if i + chunk_size < len(all_asm_uids):
            time.sleep(pause)

    linked = 0
    for nuc_uid, asm_uids in nuc_to_asm.items():
        rec = uid_records.get(nuc_uid)
        if not rec:
            continue
        chosen = _pick_assembly([asm_acc.get(u, "") for u in asm_uids])
        if chosen:
            rec.assembly_acc = chosen
            linked += 1
    log(f"Resolved real assemblies for {linked}/{len(records)} accessions via elink")


# --------------------------------------------------------------------------- #
#  Cache-backed resolver                                                       #
# --------------------------------------------------------------------------- #
class NuccoreResolver:
    """Resolve accessions to :class:`NuccoreRecord`, backed by a JSON cache.

    Parameters
    ----------
    cache_path:
        JSON file to read/write (committed for reproducibility).  ``None`` = no
        persistence (used in tests).
    fetch_fn:
        ``fetch_fn(list_of_accessions) -> {accession: NuccoreRecord}`` for the
        accessions not already cached.  Defaults to :func:`entrez_esummary`
        (network).  Tests inject a dict-backed stub so nothing hits the wire.
    assembly_fn:
        ``assembly_fn(records_dict, log) -> None`` that fills each record's real
        ``assembly_acc`` in place.  Defaults to :func:`entrez_resolve_assemblies`
        (network).  Injectable so the assembly-grouping path is unit-testable
        offline, symmetric with ``fetch_fn``.
    """

    _RECORD_FIELDS = frozenset(f.name for f in fields(NuccoreRecord))

    def __init__(
        self,
        cache_path: Optional[str] = None,
        fetch_fn: Optional[Callable[[List[str]], Dict[str, NuccoreRecord]]] = None,
        assembly_fn: Optional[Callable[..., None]] = None,
    ) -> None:
        self.cache_path = cache_path
        self.fetch_fn = fetch_fn or entrez_esummary
        self.assembly_fn = assembly_fn or entrez_resolve_assemblies
        self.records: Dict[str, NuccoreRecord] = {}
        if cache_path and os.path.exists(cache_path):
            self._load_cache(cache_path)

    def _load_cache(self, path: str) -> None:
        with open(path) as fh:
            blob = json.load(fh)
        for acc, d in blob.get("records", blob).items():
            # Tolerate cache written by a newer/older schema: keep only fields
            # the current NuccoreRecord understands, and backfill the required
            # ``accession`` from the dict key if the record body omits it.
            known = {k: v for k, v in d.items() if k in self._RECORD_FIELDS}
            known.setdefault("accession", acc)
            self.records[acc] = NuccoreRecord(**known)

    def save(self) -> None:
        if not self.cache_path:
            return
        os.makedirs(os.path.dirname(os.path.abspath(self.cache_path)), exist_ok=True)
        payload = {
            "schema": 1,
            "records": {acc: asdict(rec) for acc, rec in sorted(self.records.items())},
        }
        # Atomic write: a crash mid-write must never truncate the committed cache,
        # and a failed serialization must not leave an orphaned temp file behind.
        tmp = f"{self.cache_path}.tmp.{os.getpid()}"
        try:
            with open(tmp, "w") as fh:
                json.dump(payload, fh, indent=1, sort_keys=True)
                fh.write("\n")
            os.replace(tmp, self.cache_path)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)

    def resolve_assemblies(self, records: Dict[str, NuccoreRecord], log=lambda m: None) -> None:
        """Fill real GCF/GCA assembly accessions on the given records in place."""
        unlinked = {a: r for a, r in records.items() if not is_real_assembly(r.assembly_acc)}
        if unlinked:
            log(f"Resolving real assemblies via elink for {len(unlinked)} accession(s) ...")
            self.assembly_fn(unlinked, log=log)

    def resolve_all(
        self, accessions: Iterable[str], allow_fetch: bool = True, log=lambda m: None
    ) -> Dict[str, NuccoreRecord]:
        """Return records for ``accessions``, fetching any cache misses.

        Raises ``LookupError`` if fetching is disabled and something is missing.
        """
        wanted = [strip_version(a) for a in accessions]
        missing = sorted({a for a in wanted if a not in self.records})
        if missing:
            if not allow_fetch:
                raise LookupError(
                    f"{len(missing)} accession(s) not in cache and fetching disabled: "
                    f"{missing[:5]}{'...' if len(missing) > 5 else ''}"
                )
            log(f"Fetching {len(missing)} uncached accession(s) from NCBI ...")
            fetched = self.fetch_fn(missing)
            self.records.update(fetched)
            still_missing = [a for a in missing if a not in self.records]
            if still_missing:
                log(f"WARNING: NCBI returned no record for: {still_missing}")
        return {a: self.records[a] for a in wanted if a in self.records}
