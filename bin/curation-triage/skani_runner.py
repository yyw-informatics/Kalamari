#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cached, mockable skani layer for Tier 2.

skani computes ANI **and** the aligned fraction, which is why the design picked it:
the aligned fraction is the QC gate that stops a high ANI over a sliver of the
genome from being read as "same species".  The gate itself lives in ``tier2.py``
(pure); this module only produces the numbers.

Same pattern as every other external step in the tool: an **injectable runner**
(default = the real subprocess) plus a **committed JSON cache** of the results,
keyed by ``query|reference``.  The FASTA files stay in the gitignored genome cache
(``genomes.py``); only the handful of numbers per comparison is committed, so a
monthly run and every unit test replay offline without shipping sequence.

The cache ``meta`` records the resolved skani version.  A version bump is a tooling
epoch: it re-opens dismissed Tier-2 ANI leads, so a moved number is re-reviewed as a
tooling change rather than read as biology.
"""

from __future__ import annotations

import json
import os
import subprocess
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

from tier2 import AniResult

DEFAULT_SKANI_BIN = "skani"

# skani dist column headers.  Named rather than positional so a column added
# upstream cannot silently shift the parse.
COL_REF_FILE = "Ref_file"
COL_QUERY_FILE = "Query_file"
COL_ANI = "ANI"
COL_AF_REF = "Align_fraction_ref"
COL_AF_QUERY = "Align_fraction_query"


# skani silently drops pairs it considers too distant, and BOTH defaults are traps
# for a curation tool: `-s 80` screens out pairs below ~80% identity and `--min-af 15`
# drops pairs where neither genome aligns over 15% of its length.  A dropped pair
# produces no output row at all, which reads as "not measured" when it actually means
# "very distant" -- the opposite conclusion.  We widen both so the interesting range is
# visible, and record any pair that is still screened out as a RESULT (see AniCache).
DEFAULT_SCREEN = 70.0     # skani -s
DEFAULT_MIN_AF = 5.0      # skani --min-af


def skani_command(query_fasta: str, ref_fastas: Sequence[str],
                  skani_bin: str = DEFAULT_SKANI_BIN,
                  screen: float = DEFAULT_SCREEN,
                  min_af: float = DEFAULT_MIN_AF) -> List[str]:
    """``skani dist`` of one query against N references (flags pinned in one place).

    No ``-o``: skani writes the table to stdout by default, and ``-o -`` creates a
    file literally named ``-`` instead of meaning stdout.
    """
    return [skani_bin, "dist", "-q", query_fasta, "-r", *list(ref_fastas),
            "-s", str(screen), "--min-af", str(min_af)]


def skani_version(skani_bin: str = DEFAULT_SKANI_BIN) -> str:
    """The pinned tool epoch, e.g. ``skani@0.2.2`` (``skani@unknown`` if absent)."""
    try:
        proc = subprocess.run([skani_bin, "--version"], capture_output=True,
                              text=True, check=False)  # noqa: S603
        out = (proc.stdout or proc.stderr).strip()
        ver = out.split()[-1] if (proc.returncode == 0 and out) else "unknown"
    except (OSError, IndexError):
        ver = "unknown"
    return "skani@%s" % ver


def run_skani(query_fasta: str, ref_fastas: Sequence[str],
              skani_bin: str = DEFAULT_SKANI_BIN,
              screen: float = DEFAULT_SCREEN, min_af: float = DEFAULT_MIN_AF) -> str:
    """Call ``skani dist`` and return its raw TSV stdout."""
    cmd = skani_command(query_fasta, ref_fastas, skani_bin, screen, min_af)
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)  # noqa: S603
    if proc.returncode != 0:
        raise RuntimeError("skani failed (rc=%s): %s" % (proc.returncode,
                                                         proc.stderr.strip()[:400]))
    return proc.stdout


def _accession_of(path: str) -> str:
    """``.../GCF_000009045.1.fna`` -> ``GCF_000009045.1`` (the cache key we use)."""
    return os.path.basename(path).rsplit(".fna", 1)[0].rsplit(".fasta", 1)[0]


def _float(cell: str) -> Optional[float]:
    try:
        return float(cell)
    except (TypeError, ValueError):
        return None


def parse_skani(text: str) -> List[AniResult]:
    """Pure: ``skani dist`` TSV -> AniResult list.

    skani reports the aligned fraction as a percentage; the gate works in 0-1, so it
    is converted here, once, at the boundary -- keeping two conventions alive inside
    the pure logic is how off-by-100 bugs get written.
    """
    lines = [ln for ln in (text or "").splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    if COL_ANI not in header:
        raise ValueError("unexpected skani output: no %r column in %r" % (COL_ANI, header[:8]))
    idx = {name: i for i, name in enumerate(header)}

    def cell(cells, name):
        i = idx.get(name)
        return cells[i] if (i is not None and i < len(cells)) else ""

    out: List[AniResult] = []
    for line in lines[1:]:
        cells = line.split("\t")
        af_ref, af_q = _float(cell(cells, COL_AF_REF)), _float(cell(cells, COL_AF_QUERY))
        out.append(AniResult(
            query=_accession_of(cell(cells, COL_QUERY_FILE)),
            reference=_accession_of(cell(cells, COL_REF_FILE)),
            ani=_float(cell(cells, COL_ANI)),
            af_query=(af_q / 100.0) if af_q is not None else None,
            af_reference=(af_ref / 100.0) if af_ref is not None else None,
        ))
    return out


class AniCache:
    """Committed ANI results, keyed ``query|reference``.

    ``build`` resolves each genome through the ``GenomeStore`` and runs the injected
    runner once per query (skani takes N references in one call).  A genome that
    cannot be fetched is recorded as a miss rather than raising, so one unavailable
    accession degrades that comparison instead of the whole run.
    """

    def __init__(self, cache_path: Optional[str] = None,
                 runner: Callable[[str, Sequence[str]], str] = None,
                 skani_bin: str = DEFAULT_SKANI_BIN,
                 screen: float = DEFAULT_SCREEN,
                 min_af: float = DEFAULT_MIN_AF) -> None:
        self.cache_path = cache_path
        self.skani_bin = skani_bin
        self.screen, self.min_af = screen, min_af
        self._runner = runner or (
            lambda q, refs: run_skani(q, refs, skani_bin, screen, min_af))
        self.results: Dict[str, dict] = {}
        self.meta: Dict[str, object] = {}
        self.misses: List[str] = []
        if cache_path and os.path.exists(cache_path):
            with open(cache_path, encoding="utf-8") as fh:
                blob = json.load(fh)
            self.results = blob.get("results", {})
            self.meta = blob.get("meta", {})

    # -- lookup ------------------------------------------------------------- #
    @staticmethod
    def key(query: str, reference: str) -> str:
        return "%s|%s" % (query, reference)

    def get(self, query: str, reference: str, reference_label: str = "") -> Optional[AniResult]:
        """Cached ANI for a pair, falling back to the reverse comparison.

        An ANI comparison is symmetric, so B-vs-A answers A-vs-B and computing both
        would double the work -- which matters most for the sibling panels, where
        every unit is compared against every other.  The two aligned fractions are
        *not* symmetric, so they are swapped when the reverse entry is used; getting
        that backwards would apply the gate to the wrong side of the comparison.
        """
        rec = self.results.get(self.key(query, reference))
        if rec is not None:
            af_q, af_r = rec.get("af_query"), rec.get("af_reference")
        else:
            rec = self.results.get(self.key(reference, query))
            if rec is None:
                return None
            af_q, af_r = rec.get("af_reference"), rec.get("af_query")
        return AniResult(query=query, reference=reference, ani=rec.get("ani"),
                         af_query=af_q, af_reference=af_r,
                         screened_out=bool(rec.get("screened_out")),
                         screen_note=rec.get("screen_note", ""),
                         reference_label=reference_label or rec.get("reference_label", ""))

    def has(self, query: str, reference: str) -> bool:
        return (self.key(query, reference) in self.results
                or self.key(reference, query) in self.results)

    def covers(self, pairs: Iterable[Tuple[str, str]]) -> bool:
        return all(self.has(q, r) for q, r in pairs)

    def epoch(self) -> str:
        return str(self.meta.get("skani_version") or "skani@unset")

    # -- build -------------------------------------------------------------- #
    def build(self, pairs: Sequence[Tuple[str, str]], store, labels: Dict[str, str] = None,
              version_epoch: Optional[str] = None, log=lambda m: None) -> None:
        """Compute every (query, reference) pair not already cached."""
        labels = labels or {}
        todo: Dict[str, List[str]] = {}
        for query, reference in pairs:
            if not self.has(query, reference):     # the reverse comparison also counts
                todo.setdefault(query, []).append(reference)
        for query, refs in sorted(todo.items()):
            try:
                qpath = store.ensure(query)
            except Exception as exc:            # noqa: BLE001 - degrade, never crash the run
                self.misses.append(query)
                log("skani: cannot resolve query genome %s (%s)" % (query, exc))
                continue
            rpaths, kept = [], []
            for ref in sorted(set(refs)):
                try:
                    rpaths.append(store.ensure(ref))
                    kept.append(ref)
                except Exception as exc:        # noqa: BLE001
                    self.misses.append(ref)
                    log("skani: cannot resolve reference genome %s (%s)" % (ref, exc))
            if not rpaths:
                continue
            returned = set()
            for res in parse_skani(self._runner(qpath, rpaths)):
                returned.add(res.reference)
                self.results[self.key(res.query, res.reference)] = {
                    "ani": res.ani, "af_query": res.af_query,
                    "af_reference": res.af_reference,
                    "screened_out": False,
                    "reference_label": labels.get(res.reference, ""),
                }
            # A pair skani did not return was screened out, and that is a FINDING, not a
            # gap: the two genomes are either below the identity screen or share almost
            # no sequence.  Recording it keeps "too distant to report" distinct from
            # "never computed", which the caller must not confuse.
            for ref in kept:
                if ref in returned or self.has(query, ref):
                    continue
                self.results[self.key(query, ref)] = {
                    "ani": None, "af_query": None, "af_reference": None,
                    "screened_out": True,
                    "screen_note": "skani returned no row at -s %s --min-af %s"
                                   % (self.screen, self.min_af),
                    "reference_label": labels.get(ref, ""),
                }
            log("skani: %s vs %d reference(s) (%d returned, %d screened out)"
                % (query, len(kept), len(returned), len(kept) - len(returned)))
        # An epoch change re-baselines every Tier-2 ANI lead and voids every dismissal
        # in reviews.tsv, so it must record a real version change and nothing else.  A
        # probe that cannot find skani returns '@unknown', which is not a version
        # change -- it is a missing tool -- and overwriting a known epoch with it would
        # bust the whole review history on a runner with a broken install.
        epoch = version_epoch or skani_version(self.skani_bin)
        previous = str(self.meta.get("skani_version") or "")
        if epoch.endswith("@unknown") and previous:
            log("skani: version probe failed; keeping the recorded epoch %s" % previous)
            epoch = previous
        self.meta = {
            "skani_version": epoch,
            "screen": self.screen,
            "min_af": self.min_af,
            "n_results": len(self.results),
            "n_screened_out": sum(1 for r in self.results.values() if r.get("screened_out")),
        }

    def save(self) -> None:
        if not self.cache_path:
            return
        os.makedirs(os.path.dirname(os.path.abspath(self.cache_path)), exist_ok=True)
        payload = {"schema": 1, "meta": self.meta,
                   "results": dict(sorted(self.results.items()))}
        tmp = "%s.tmp.%d" % (self.cache_path, os.getpid())
        try:
            with open(tmp, "w", encoding="utf-8") as fh:
                json.dump(payload, fh, indent=1, sort_keys=True)
                fh.write("\n")
            os.replace(tmp, self.cache_path)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)
