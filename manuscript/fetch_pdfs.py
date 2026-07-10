#!/usr/bin/env python3
"""
Fetch PDFs for a manifest of articles (pmid, doi, ...), OA-first.

Order per article:  Europe PMC OA  ->  Unpaywall OA  ->  Sci-Hub (fallback).
Writes <outdir>/<pmid>.pdf and a provenance log <outdir>/_fetch_log.csv.

  python manuscript/fetch_pdfs.py <manifest.csv> <outdir> [--email you@x.com] [--delay 4]

manifest.csv must have columns: pmid, doi  (extra columns ignored).
"""
import csv, sys, time, argparse
from pathlib import Path
import requests, urllib3

# reuse the vendored Sci-Hub client from the literature_review project
sys.path.insert(0, r"C:\Users\jonghoon.kim\Workspace\literature_review\scripts")
try:
    from pubmed_search.scihub_client import SciHubClient
except Exception:
    SciHubClient = None

urllib3.disable_warnings()
HEADERS = {"User-Agent": ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
                          "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36")}


def is_pdf(b): return b[:4] == b"%PDF"


def get(url, **kw):
    return requests.get(url, headers=HEADERS, timeout=90, verify=False, allow_redirects=True, **kw)


def pmcid_for(pmid, email):
    try:
        r = get("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/",
                params={"ids": pmid, "format": "json", "tool": "typhoid-ori", "email": email})
        recs = r.json().get("records", [])
        return recs[0].get("pmcid") if recs and recs[0].get("pmcid") else None
    except Exception:
        return None


def try_europepmc(pmcid, dest):
    url = f"https://europepmc.org/backend/ptpmcrender.fcgi?accid={pmcid}&blobtype=pdf"
    try:
        r = get(url)
        if r.ok and is_pdf(r.content):
            dest.write_bytes(r.content); return "europepmc"
    except Exception:
        pass
    return None


def try_unpaywall(doi, email, dest):
    if not doi: return None
    try:
        r = get(f"https://api.unpaywall.org/v2/{doi}", params={"email": email})
        loc = (r.json() or {}).get("best_oa_location") or {}
        purl = loc.get("url_for_pdf")
        if purl:
            rr = get(purl)
            if rr.ok and is_pdf(rr.content):
                dest.write_bytes(rr.content); return "unpaywall"
    except Exception:
        pass
    return None


def try_scihub(doi, client, dest):
    if not doi or client is None: return None
    try:
        res = client.fetch(doi)  # returns dict with 'pdf' bytes or url on this vendored client
        data = res.get("pdf") if isinstance(res, dict) else None
        if isinstance(data, (bytes, bytearray)) and is_pdf(data):
            dest.write_bytes(data); return "scihub"
        url = res.get("url") if isinstance(res, dict) else None
        if url:
            rr = get(url)
            if rr.ok and is_pdf(rr.content):
                dest.write_bytes(rr.content); return "scihub"
    except Exception:
        pass
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest"); ap.add_argument("outdir")
    ap.add_argument("--email", default="kimfinale@gmail.com")
    ap.add_argument("--delay", type=float, default=4.0)
    ap.add_argument("--no-scihub", action="store_true")
    a = ap.parse_args()

    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(a.manifest, encoding="utf-8")))
    client = None if a.no_scihub or SciHubClient is None else SciHubClient()

    log = []
    for i, row in enumerate(rows, 1):
        pmid = (row.get("pmid") or "").strip()
        doi = (row.get("doi") or "").strip() or None
        dest = outdir / f"{pmid}.pdf"
        if dest.exists() and dest.stat().st_size > 5000:
            print(f"[{i}/{len(rows)}] {pmid} SKIP (exists)"); log.append((pmid, doi, "exists", dest.stat().st_size)); continue

        src = None
        pmcid = pmcid_for(pmid, a.email)
        if pmcid:
            src = try_europepmc(pmcid, dest)
        if not src:
            src = try_unpaywall(doi, a.email, dest)
        if not src and not a.no_scihub:
            src = try_scihub(doi, client, dest)

        sz = dest.stat().st_size if dest.exists() else 0
        status = src or "FAILED"
        print(f"[{i}/{len(rows)}] {pmid} pmc={pmcid or '-'} -> {status} ({sz:,}B)")
        log.append((pmid, doi, status, sz))
        time.sleep(a.delay)

    with open(outdir / "_fetch_log.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["pmid", "doi", "source", "bytes"]); w.writerows(log)
    got = sum(1 for r in log if r[2] not in ("FAILED",))
    print(f"\nDONE: {got}/{len(rows)} obtained. Log -> {outdir/'_fetch_log.csv'}")


if __name__ == "__main__":
    main()
