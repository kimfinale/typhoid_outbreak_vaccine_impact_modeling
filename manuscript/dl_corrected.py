import requests, urllib3
from pathlib import Path
urllib3.disable_warnings()
H = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/120 Safari/537.36"}
out = Path("manuscript/data/pdfs/corrected"); out.mkdir(parents=True, exist_ok=True)
targets = {"Davis2018_Harare_PMC5868204": "PMC5868204", "Poncin2022_HarareTCV_PMC9379662": "PMC9379662"}
for name, acc in targets.items():
    dest = out / f"{name}.pdf"
    try:
        r = requests.get(f"https://europepmc.org/backend/ptpmcrender.fcgi?accid={acc}&blobtype=pdf",
                         headers=H, timeout=90, verify=False)
        if r.ok and r.content[:4] == b"%PDF":
            dest.write_bytes(r.content); print(f"{name}: OK {len(r.content):,} bytes")
        else:
            print(f"{name}: FAILED ({r.status_code}, {r.content[:30]!r})")
    except Exception as e:
        print(f"{name}: ERROR {e}")
