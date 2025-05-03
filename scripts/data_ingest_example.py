#!/usr/bin/env python3
"""
scripts/data_ingest_example.py

Template for a universal data ingestion script. Follow this pattern to create source‐specific ingestors:

1. Set SOURCE_NAME to a unique identifier for your data source.
2. If your source requires an API key, place it in scripts/getenv.py and import it.
3. Implement `fetch_records()` to yield raw entries from your API or file.
4. Implement `transform_record(raw)` to map source‐specific fields to a unified record schema.
5. Adjust CHUNK_SIZE, BATCH_DELAY, and OUTPUT_DIR as needed.

COMMON RECORD FIELDS:
  - id:         unique record ID (string)
  - source:     SOURCE_NAME
  - source_id:  original ID from the source
  - name:       human‐readable name or formula
  - type:       e.g. 'molecule' or 'crystal'
  - formula:    chemical formula
  - smiles:     SMILES string (if any)
  - inchikey:   InChIKey (if any)
  - composition:dict of element→count
  - structure:  raw structure data (if any)
  - sites:      site list (if any)
  - properties: dict of source‐computed properties

Usage:
    python scripts/data_ingest_example.py
"""

import os
import json
import time

# --- Optional API key import (if needed) ---
try:
    from scripts.getenv import API_KEY
except ImportError:
    API_KEY = None

# --- Configuration ---
SOURCE_NAME = 'example_source'
CHUNK_SIZE = 100  # number of records per batch in fetch_records
BATCH_DELAY = 0.1  # seconds to sleep between batches
OUTPUT_DIR = os.path.abspath(
    os.path.join(__file__, '..', '../data/raw')
)
OUTPUT_FILE = os.path.join(OUTPUT_DIR, f"{SOURCE_NAME}.json")


# --- 1) Fetch raw records from the data source ---
def fetch_records():
    """
    Yield raw entries from the source. Replace this stub with
    API calls or file reads for your specific source.
    Each yielded item `raw` can be a dict, object, etc.
    """
    # Example (pseudo-code):
    # client = SomeApiClient(api_key=API_KEY)
    # total = client.count_items()
    # for offset in range(0, total, CHUNK_SIZE):
    #     batch = client.get_items(offset=offset, limit=CHUNK_SIZE)
    #     for item in batch:
    #         yield item

    raise NotImplementedError("fetch_records() must be implemented per source")


# --- 2) Transform raw entry into unified schema ---
def transform_record(raw):
    """
    Map a source-specific record `raw` to the common fields:
      id, source, source_id, name, type, formula, smiles,
      inchikey, composition, structure, sites, properties
    """
    rec = {
        'id': None,
        'source': SOURCE_NAME,
        'source_id': None,
        'name': None,
        'type': None,
        'formula': None,
        'smiles': None,
        'inchikey': None,
        'composition': {},
        'structure': None,
        'sites': None,
        'properties': {},
    }
    # Example mappings (override for your source):
    # rec['source_id'] = raw.get('material_id') or raw.get('CID')
    # rec['id']        = rec['source_id']
    # rec['name']      = raw.get('formula') or raw.get('name')
    # rec['type']      = 'molecule' or 'crystal'
    # rec['formula']   = raw.get('formula')
    # rec['smiles']    = raw.get('smiles')
    # rec['inchikey']  = raw.get('inchikey')
    # rec['composition'] = raw.get('composition_reduced') or {}
    # rec['properties']  = raw.get('properties', {})
    return rec


# --- 3) Main ingestion logic with progress reporting ---
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    out = []
    count = 0
    start = time.time()

    for raw in fetch_records():
        rec = transform_record(raw)
        out.append(rec)
        count += 1

        # batch progress report
        if count % CHUNK_SIZE == 0:
            elapsed = time.time() - start
            print(f"Processed {count} records ({elapsed:.1f}s elapsed)", flush=True)
            time.sleep(BATCH_DELAY)

    # write output
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"Wrote {count} records to {OUTPUT_FILE} (total {time.time() - start:.1f}s)")


if __name__ == '__main__':
    main()
