#!/usr/bin/env python3
"""
scripts/data_ingest_pubchem.py

Fetch compound data from PubChem via pubchempy and write to data/raw/pubchem.json.
Looks first for a CID list file; if missing or empty, falls back to an embedded list.
"""

import os
import sys
import json
import time
import traceback

import pubchempy as pcp

# --- Configuration ---
RAW_DIR     = os.path.abspath(os.path.join(__file__, '..', '../data/raw'))
CID_FILE    = os.path.join(RAW_DIR, 'CID.txt')
OUTPUT_FILE = os.path.join(RAW_DIR, 'pubchem.json')

BATCH_SIZE  = 100       # PubChem allows max 100 per batch
BATCH_DELAY = 0.1       # seconds between batches

# Fallback manual list of CIDs (edit as needed)
MANUAL_CIDS = [2244, 702]

# Properties to fetch
PROPERTIES = [
    "MolecularFormula",
    "IsomericSMILES",
    "InChIKey",
    "IUPACName",
    "MolecularWeight",
    "XLogP"
]

# --- 1) Determine CIDs to fetch ---
cids = []
source_desc = None

if os.path.exists(CID_FILE):
    with open(CID_FILE, 'r') as f:
        cids = [int(line.strip()) for line in f if line.strip().isdigit()]
    if cids:
        source_desc = f"from {os.path.basename(CID_FILE)}"
    else:
        print(f"CID file is empty; falling back to MANUAL_CIDS")

if not cids and MANUAL_CIDS:
    cids = MANUAL_CIDS.copy()
    source_desc = "from hardcoded MANUAL_CIDS"

if not cids:
    print("No CIDs provided (file missing/empty and MANUAL_CIDS empty). Exiting.")
    sys.exit(1)

total = len(cids)
print(f"PubChem ingest: {total} CIDs {source_desc}")

# --- 2) Fetch in batches ---
out = []
fetched = 0
start_time = time.time()

for batch_idx in range(0, total, BATCH_SIZE):
    batch = cids[batch_idx: batch_idx + BATCH_SIZE]
    print(f"  Batch {batch_idx//BATCH_SIZE + 1} / {(total-1)//BATCH_SIZE + 1} (size {len(batch)})…")

    try:
        props_list = pcp.get_properties(PROPERTIES, batch, namespace='cid')
    except Exception as e:
        print(f"    ERROR fetching batch {batch_idx//BATCH_SIZE + 1}: {e}", file=sys.stderr)
        traceback.print_exc()
        continue

    for entry in props_list:
        rec = {
            "id":         str(entry["CID"]),
            "source":     "pubchem",
            "source_id":  str(entry["CID"]),
            "name":       entry.get("IUPACName") or entry.get("MolecularFormula") or f"CID_{entry['CID']}",
            "type":       "molecule",
            "formula":    entry.get("MolecularFormula"),
            "smiles":     entry.get("IsomericSMILES"),
            "inchikey":   entry.get("InChIKey"),
            "properties": {
                "mw":   entry.get("MolecularWeight"),
                "logP": entry.get("XLogP")
            }
        }
        out.append(rec)
        fetched += 1

    time.sleep(BATCH_DELAY)

# --- 3) Summary and write ---
elapsed = time.time() - start_time
print(f"Fetched {fetched}/{total} compounds in {elapsed:.1f}s")

os.makedirs(RAW_DIR, exist_ok=True)
with open(OUTPUT_FILE, 'w') as f:
    json.dump(out, f, indent=2)

print(f"Wrote {len(out)} records to {OUTPUT_FILE}")
