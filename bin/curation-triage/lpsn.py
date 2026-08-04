#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cached, mockable LPSN layer for the Tier 1 J2 "coverage feeder".

J2's second source is *newly validly-published species* (LPSN = the List of
Prokaryotic names with Standing in Nomenclature) within the pathogen genera we track:
a validly-published species in a tracked genus that Kalamari does not carry is a
coverage gap worth an SME's attention.

LPSN publishes a downloadable complete-list CSV.  Like ``ani_reports.py`` /
``datasets_summary.py``, reading it is quarantined behind an **injectable opener**
(an operator's downloaded file, or the link it came from; tests pass fixture lines)
and distilled into a **committed JSON cache**, so a run and every test are offline.
The pure genus-gating lives in ``leads.j2_leads``; this module only reads + distils.

Epoch pinning: LPSN has real dated releases, so the cache ``meta`` records the LPSN
``release`` (a date/DOI) -- an LPSN update is a tooling epoch that re-opens dismissed
J2 leads for re-review.  The exact CSV column names are pinned in ``LPSN_COLS`` (one
place, header-validated) so a schema change fails loudly instead of mis-reading.

How the cache is refreshed
--------------------------
An LPSN account exists and the committed ``src/curation-triage/lpsn_cache.json`` was
built from a real complete-list download (release 2026-08-04).  The refresh is
**operator-driven, not unattended**: the public ``/downloads`` address is a landing
page, and the complete-list export sits behind the account, so a human signs in,
downloads the CSV and re-runs the distiller (``LPSN_REFRESH_HELP`` has the exact
command).  The cache IS committed, so ``--with-lpsn`` works fully offline afterwards.

The intended unattended path is LPSN's REST API at ``https://api.lpsn.dsmz.de/``
(``LPSN_API_ROOT``).  It is deliberately left as **future work for a developer**: it
needs the same account's credentials held as a CI secret, a login/token exchange, and
paging over ~34 000 records, none of which the CSV path needs.  When someone builds
it, it becomes a third opener next to ``open_csv_file`` / ``open_url_lines`` -- the
distiller, the cache format and the epoch all stay as they are.
"""

from __future__ import annotations

import csv
import json
import os
import re
import urllib.request
from dataclasses import asdict, dataclass, fields
from typing import Callable, Dict, Iterable, Iterator, List, Optional, Sequence

# LPSN's download LANDING PAGE -- not a CSV.  The complete-list export sits behind a
# registered account, so ``LpsnCache`` has no working default opener (see below).
LPSN_DOWNLOAD_PAGE = "https://lpsn.dsmz.de/downloads"

# LPSN's REST API: the intended unattended refresh path, NOT built yet (see the module
# docstring).  Recorded here so the next developer starts from the right address.
LPSN_API_ROOT = "https://api.lpsn.dsmz.de/"

# What an operator with LPSN access has to do to refresh the cache.  Printed verbatim
# by every failure path, so nobody has to reverse-engineer it from the source.
LPSN_REFRESH_HELP = (
    "LPSN's complete-list CSV is behind a registered account, so the cache is built by\n"
    "hand from a download; the committed src/curation-triage/lpsn_cache.json came from\n"
    "exactly this. To refresh it:\n"
    "  1. sign in at %s and download the complete list as CSV;\n"
    "  2. python3 bin/curation-triage/tier1_metadata.py --offline --no-write \\\n"
    "         --with-j2 --with-lpsn \\\n"
    "         --lpsn-csv <downloaded.csv> --lpsn-release <YYYY-MM-DD> \\\n"
    "         --lpsn-cache src/curation-triage/lpsn_cache.json\n"
    "     (--offline is right here: distilling a downloaded CSV needs no network, and\n"
    "      without it the run would also rebuild the 2.7 GB assembly-reports cache);\n"
    "  3. if step 2 fails on a missing column, pin the real header name in\n"
    "     LPSN_COLS (bin/curation-triage/lpsn.py) and bump --lpsn-release.\n"
    "An unattended refresh would go through LPSN's REST API (%s) with the same\n"
    "account's credentials as a CI secret -- future work, not built."
    % (LPSN_DOWNLOAD_PAGE, LPSN_API_ROOT)
)

# Column-name contract, PINNED against the real complete-list export (release
# 2026-08-04), whose full header is:
#
#   genus_name, sp_epithet, subsp_epithet, reference, status, authors, address,
#   risk_grp, nomenclatural_type, record_no, record_lnk
#
# Kept in ONE place so a schema change is a loud epoch bump; adjust here + the release
# string together when LPSN changes format.  Each entry is a TUPLE of accepted header
# names, tried in order, so re-pinning a renamed column stays a one-line edit here.
#
# Note what is NOT in the header: there is **no date column at all**.  The date of
# valid publication has to be derived from the authority string -- see
# ``valid_publication_year``.
LPSN_COLS = {
    "genus": ("genus_name",),
    "species": ("sp_epithet",),
    "subspecies": ("subsp_epithet",),
    "status": ("status",),          # 'VL; sp. nov.; ...; correct name' / '...; synonym'
    "authority": ("authors",),      # 'Tahon et al. 2018' -- the ONLY source of the year
    "reference": ("reference",),    # the citation; year fallback when 'authors' has none
}

# The logical columns whose absence is a schema change we must NOT paper over: exactly
# what the J2 gate depends on.  ``authority`` is in here because it now carries the
# valid-publication year, and ``--lpsn-since`` is a *gate* -- fed by a silently-blank
# column it would quietly drop every gap it was meant to surface.  ``reference`` is a
# fallback only (no species row in the 2026-08-04 release needed it), and
# ``subspecies`` only decides whether a row is collapsed, so neither is required.
LPSN_REQUIRED_COLS = ("genus", "species", "status", "authority")

# A four-digit year in an authority string.  Bounded to 1500-2099 so a page range or a
# journal volume ('68:3379-3393') can never be mistaken for a year.
_YEAR_RE = re.compile(r"\b(1[5-9]\d\d|20\d\d)\b")

# Everything after 'emend.' / 'corrig.' is a LATER amendment of an existing name -- an
# emended description or a corrected spelling -- not the act that validly published it.
# Cutting here is what keeps 'Anandham et al. 2011 emend. Radolfova-Krizova et al. 2026'
# at 2011.  ('corrig.' does not occur in the 2026-08-04 release, but it is the same
# class of post-hoc annotation, so it is cut too.)
_AMENDMENT_RE = re.compile(r"\b(?:emend|corrig)\.", re.IGNORECASE)


def valid_publication_year(authority: str, reference: str = "") -> str:
    """The year the name became validly published, derived from LPSN's authority string.

    LPSN's complete-list CSV has no date column, so the year is read out of ``authors``.
    The rule is: **cut at ``emend.`` / ``corrig.``, then take the LAST four-digit year
    in what remains** (falling back to the whole string, then to ``reference``).  That
    single rule gives the right answer for all four shapes LPSN uses:

      ``Tahon et al. 2018``                              -> 2018  plain valid publication
      ``(Bouvet et al. 1989) Kawamura et al. 1995``      -> 1995  a new combination cites
                                                                 the BASIONYM year first;
                                                                 the current name dates
                                                                 from the combination
      ``Anandham et al. 2011 emend. ... et al. 2026``    -> 2011  2026 is an emendation
      ``(Murray et al. 1926) Pirie 1940
        (Approved Lists 1980)``                          -> 1980  the Approved Lists are
                                                                 when the name gained
                                                                 standing under the code

    Returns a bare YEAR, not a full date -- that is all LPSN gives, and it is the right
    granularity anyway: "validly published" is an event dated by its publication, and
    the J2 ``--lpsn-since`` cutoff compares left-anchored text, so a year works against
    a ``YYYY`` or ``YYYY-MM-DD`` cutoff alike.  Empty string when no year can be found.
    """
    for text in (authority, reference):
        text = (text or "").strip()
        if not text:
            continue
        head = _AMENDMENT_RE.split(text)[0]
        for candidate in (head, text):       # head first; whole string as a fallback
            years = _YEAR_RE.findall(candidate)
            if years:
                return years[-1]
    return ""


def taxonomic_standing(status: str) -> str:
    """The taxonomic standing LPSN assigns a name: the LAST semicolon-clause of ``status``.

    ``status`` is a semicolon-separated trail, e.g.
    ``VL; sp. nov.; validly published under the ICNP; correct name``.  Only the last
    clause is the standing (``correct name`` / ``synonym`` / ``misspelling`` /
    ``orphaned species`` / ...); its trailing ``, recommended for medical use`` rider is
    dropped so ``correct name, recommended for medical use`` reads as ``correct name``.

    Reading the last clause -- rather than searching the whole field -- matters: 12 rows
    in the 2026-08-04 release read ``... validly published under the ICNP, inappropriate
    correction; misspelling``.  A substring test for "correct" matches those, so it
    would hand an SME 12 misspellings to add as new species.
    """
    last_clause = (status or "").split(";")[-1].strip().lower()
    return last_clause.split(",")[0].strip()


@dataclass
class LpsnRecord:
    """One validly-published prokaryotic name (species granularity)."""

    genus: str = ""
    species: str = ""               # the full binomial, e.g. 'Escherichia coli'
    status: str = ""
    authority: str = ""
    # Year of valid publication, DERIVED from ``authority`` (LPSN ships no date column).
    # Named ``valid_date`` because the J2 cutoff treats it as a left-anchored date.
    valid_date: str = ""

    @property
    def is_correct_name(self) -> bool:
        return taxonomic_standing(self.status) == "correct name"


def _binomial(genus: str, epithet: str) -> str:
    genus, epithet = (genus or "").strip(), (epithet or "").strip()
    return ("%s %s" % (genus, epithet)).strip()


def column_candidates(logical: str) -> tuple:
    """The accepted header spellings for one logical column, as a tuple.

    A single string is accepted as well as a tuple, so re-pinning a renamed column
    really is the one-line edit the comment above promises --
    ``"authority": "authors"`` works exactly like the tuple form.
    """
    value = LPSN_COLS.get(logical, ())
    return (value,) if isinstance(value, str) else tuple(value)


def resolve_columns(fieldnames: Sequence[str]) -> Dict[str, str]:
    """Map each logical LPSN column to the real header name this CSV uses.

    Every ``LPSN_COLS`` entry is a list of accepted spellings; the first one present in
    the header wins.  A required column (``LPSN_REQUIRED_COLS``) with no match raises,
    quoting the header we actually got -- that is how an LPSN rename of, say, ``authors``
    surfaces as a loud, fixable failure instead of a silently blank year.
    """
    present = set(fieldnames or ())
    resolved: Dict[str, str] = {}
    missing: List[str] = []
    for logical in LPSN_COLS:
        candidates = column_candidates(logical)
        hit = next((c for c in candidates if c in present), "")
        if hit:
            resolved[logical] = hit
        elif logical in LPSN_REQUIRED_COLS:
            missing.append("%s (tried %s)" % (logical, ", ".join(candidates)))
    if missing:
        raise ValueError(
            "LPSN CSV: missing expected column(s): %s. Header was: %s. This is a schema "
            "change -> pin the real name in LPSN_COLS and bump the release (epoch)."
            % ("; ".join(missing), ", ".join(sorted(present)) or "(empty)"))
    return resolved


def distill_lpsn(lines: Iterable[str]) -> List[LpsnRecord]:
    """Pure: LPSN CSV lines -> species-level records, deduped, status/year preserved.

    Header-validated via ``resolve_columns``; a missing required column raises loudly.
    Subspecies rows (a non-empty subspecies epithet) are collapsed to their species.
    The valid-publication year is derived per row by ``valid_publication_year`` -- LPSN
    ships no date column, so this is the only place a date can come from.

    Deliberately NOT filtered: synonyms and old names stay in the distilled cache with
    their ``status`` and ``valid_date`` intact.  The cache is a faithful copy of the
    source; deciding which names count as a coverage gap is policy, and policy lives in
    the pure lead layer (``leads.j2_leads``).  Filtering here would bake today's gate
    into the committed artifact and force a re-download to change our minds.
    """
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        return []
    cols = resolve_columns(reader.fieldnames)
    # Dedup by collapsed binomial, but STATUS-AWARE: a correct-name row must win over a
    # synonym/subspecies row for the same binomial regardless of CSV row order (else a
    # subspecies 'synonym' row listed first would shadow the real species and drop a
    # genuine coverage gap).
    by_binom: Dict[str, int] = {}
    out: List[LpsnRecord] = []
    for r in reader:
        genus = (r.get(cols["genus"]) or "").strip()
        sp = (r.get(cols["species"]) or "").strip()
        if not genus or not sp:
            continue
        binom = _binomial(genus, sp)
        authority = (r.get(cols["authority"]) or "").strip()
        reference = (r.get(cols.get("reference", ""), "") or "").strip()
        rec = LpsnRecord(
            genus=genus, species=binom,
            status=(r.get(cols["status"]) or "").strip(),
            authority=authority,
            valid_date=valid_publication_year(authority, reference),
        )
        if binom in by_binom:
            stored = out[by_binom[binom]]
            if rec.is_correct_name and not stored.is_correct_name:
                out[by_binom[binom]] = rec       # correct name supersedes a synonym/subsp
            continue
        by_binom[binom] = len(out)
        out.append(rec)
    return out


# --------------------------------------------------------------------------- #
#  Injectable opener + cache                                                    #
# --------------------------------------------------------------------------- #
def open_url_lines(url: str, timeout: int = 300) -> Iterator[str]:
    req = urllib.request.Request(url, headers={"User-Agent": "kalamari-curation-triage"})
    resp = urllib.request.urlopen(req, timeout=timeout)  # noqa: S310 - fixed LPSN host
    for raw in resp:
        yield raw.decode("utf-8", "replace")


def open_csv_file(path: str) -> Iterator[str]:
    """Read an operator-downloaded LPSN complete-list CSV (the only refresh path)."""
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for line in fh:
            yield line


def is_remote_source(location: str) -> bool:
    """True when ``--lpsn-csv`` names a URL rather than a file on disk.

    The caller needs this to keep ``--offline`` honest: distilling a downloaded CSV is
    a plain file read and stays offline, but a short-lived per-account download LINK is
    a network fetch and must be refused under ``--offline``.
    """
    loc = (location or "").strip().lower()
    return loc.startswith("http://") or loc.startswith("https://")


def open_lpsn_source(location: str) -> Iterator[str]:
    """Open a complete-list CSV given a local path OR a direct http(s) link.

    A URL is accepted because LPSN's per-account download links are short-lived: an
    operator holding one should not have to save the file first just to feed it in.
    """
    loc = (location or "").strip()
    return open_url_lines(loc) if is_remote_source(loc) else open_csv_file(loc)


def _no_opener() -> Iterator[str]:
    """The default 'opener': there isn't one yet, and pretending otherwise wastes a run.

    Fetching ``LPSN_DOWNLOAD_PAGE`` returns HTML, which csv.DictReader would report as
    a missing-column schema error -- a confusing way to say "you need an account".  So
    say it directly instead.  An API client (``LPSN_API_ROOT``) would go here.
    """
    raise RuntimeError("LPSN cannot be downloaded unattended.\n" + LPSN_REFRESH_HELP)


def _record_from_dict(d: dict) -> LpsnRecord:
    known = {f.name for f in fields(LpsnRecord)}
    return LpsnRecord(**{k: v for k, v in d.items() if k in known})


class LpsnCache:
    """Distilled LPSN species list + its release epoch; injectable opener, committed
    JSON cache -- mirrors ``ani_reports.AssemblyReportsCache``."""

    def __init__(self, cache_path: Optional[str] = None,
                 opener: Callable[[], Iterable[str]] = None) -> None:
        self.cache_path = cache_path
        self._opener = opener or _no_opener
        self.records: List[LpsnRecord] = []
        self.release: str = ""
        self.meta: Dict[str, object] = {}
        if cache_path and os.path.exists(cache_path):
            self._load_cache(cache_path)

    def _load_cache(self, path: str) -> None:
        with open(path, encoding="utf-8") as fh:
            blob = json.load(fh)
        self.records = [_record_from_dict(d) for d in blob.get("records", [])]
        self.release = blob.get("release", "")
        self.meta = blob.get("meta", {})

    def save(self) -> None:
        if not self.cache_path:
            return
        os.makedirs(os.path.dirname(os.path.abspath(self.cache_path)), exist_ok=True)
        payload = {
            "schema": 1,
            "meta": self.meta,
            "release": self.release,
            "records": [asdict(r) for r in self.records],
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

    def build(self, release: str, log=lambda m: None) -> None:
        self.records = distill_lpsn(self._opener())
        self.release = release
        n_correct = sum(1 for r in self.records if r.is_correct_name)
        n_dated = sum(1 for r in self.records if r.valid_date)
        self.meta = {"lpsn_release": release, "n_records": len(self.records),
                     "n_correct_names": n_correct, "n_with_valid_date": n_dated,
                     # Recorded so a reader of the committed cache knows which refresh
                     # path produced it (the API path would say so here instead).
                     "source": "LPSN complete-list CSV"}
        log("LPSN: distilled %d name(s) -- %d correct name(s), %d with a valid-publication "
            "year (release %s)" % (len(self.records), n_correct, n_dated, release))

    def in_genera(self, genera: Iterable[str]) -> List[LpsnRecord]:
        gset = {g.strip().lower() for g in genera if g}
        return [r for r in self.records
                if r.is_correct_name and r.genus.strip().lower() in gset]
