# scripts/data_ingest_mp.py

import json, math
from scripts.getenv import MAPI_KEY
from mp_api.client import MPRester

CHUNK_SIZE = 500
criteria = {"deprecated": False}

with MPRester(api_key=MAPI_KEY) as m:
    # 1) count how many non-deprecated materials
    total = m.materials.summary.count(criteria)
    num_chunks = math.ceil(total / CHUNK_SIZE)
    print(f"Fetching {total} MP entries over {num_chunks}×{CHUNK_SIZE}…")

    # 2) fetch every field for each entry
    entries = m.materials.search(
        deprecated=False,
        chunk_size=CHUNK_SIZE,
        num_chunks=num_chunks,
        all_fields=True
    )

out = []
for doc in entries:
    rec = doc.dict()                 # full MaterialsDoc → dict
    rec["source"]    = "materialsproject"
    rec["source_id"] = doc.material_id
    # give it a human‐readable name too
    rec["name"]      = doc.formula_pretty or doc.material_id
    out.append(rec)

with open("../data/raw/mp_full.json","w") as f:
    json.dump(out, f, indent=2)

print(f"Wrote {len(out)} records to data/raw/mp_full.json")
