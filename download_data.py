"""Download Dimuon_DoubleMu.csv from CERN Open Data record 5201 (only this file)."""
import argparse
import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

import config

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
CSV = ROOT / config.DATA_FILE
REQUIRED = ["Run", "Event", "type1", "E1", "px1", "py1", "pz1", "pt1", "eta1", "phi1", "Q1",
            "type2", "E2", "px2", "py2", "pz2", "pt2", "eta2", "phi2", "Q2", "M"]
HEADERS = {"User-Agent": "open-data-educational-analysis"}


def log(msg: str) -> None:
    with open(ROOT / "RUN_LOG.md", "a") as f:
        f.write(f"- {time.strftime('%F %T')} {msg}\n")


def find_link_on_record_page() -> str:
    r = requests.get(config.RECORD_URL, headers=HEADERS, timeout=60)
    r.raise_for_status()
    m = re.findall(r'["\'(]([^"\'()\s]*Dimuon_DoubleMu\.csv[^"\'()\s]*)', r.text)
    if not m:
        raise RuntimeError("no link to Dimuon_DoubleMu.csv on the record page")
    link = m[0]
    if link.startswith("/"):
        link = "https://opendata.cern.ch" + link
    return link


def fetch(url: str) -> bytes:
    r = requests.get(url, headers=HEADERS, timeout=300)
    r.raise_for_status()
    return r.content


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()
    if CSV.exists() and not args.refresh:
        print("data file exists, nothing to do (use --refresh to download again)")
        return
    DATA.mkdir(exist_ok=True)
    url = config.DATA_URL
    content = None
    for attempt in range(1, config.RETRY_COUNT + 1):
        try:
            content = fetch(url)
            break
        except Exception as exc:  # noqa: BLE001
            log(f"download attempt {attempt} with {url} failed: {exc}")
            if attempt == 1:
                try:
                    url = find_link_on_record_page()
                    log(f"link found on the record page: {url}")
                except Exception as exc2:  # noqa: BLE001
                    log(f"could not read the link from the record page: {exc2}")
            if attempt < config.RETRY_COUNT:
                time.sleep(config.RETRY_WAIT_S)
    if content is None:
        log("HARD STOP B1: the data file could not be downloaded.")
        raise SystemExit("hard stop B1: download failed")
    CSV.write_bytes(content)
    df = pd.read_csv(CSV)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        log(f"HARD STOP B1: missing required columns {missing}")
        raise SystemExit(f"hard stop B1: missing columns {missing}")
    meta = {"download_time_utc": datetime.now(timezone.utc).isoformat(), "url": url,
            "file_size_bytes": len(content), "sha256": hashlib.sha256(content).hexdigest(),
            "rows": int(len(df)), "columns": list(df.columns)}
    (DATA / "metadata.json").write_text(json.dumps(meta, indent=2))
    log(f"downloaded {url}: {len(content)} bytes, {len(df)} rows, sha256 {meta['sha256']}")
    print(meta)


if __name__ == "__main__":
    main()
