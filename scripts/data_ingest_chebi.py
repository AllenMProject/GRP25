#!/usr/bin/env python3
"""
scripts/data_ingest_chebi.py

Fetch the complete ChEBI SDF (all small molecules) and write to data/raw/chebi.json
in our unified record format.
"""

import os
import json
import time
import gzip
from urllib.request import urlretrieve

from rdkit import RDLogger, Chem

# -----------------------------------------------------------------------------
# 1) Disable RDKit logging so we don’t flood the console with sanitization errors
# -----------------------------------------------------------------------------
RDLogger.DisableLog('rdApp.*')


# -----------------------------------------------------------------------------
# 2) Configuration
# -----------------------------------------------------------------------------
RAW_DIR     = os.path.abspath(os.path.join(__file__, "..", "../data/raw"))
os.makedirs(RAW_DIR, exist_ok=True)

SDF_URL     = "ftp://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete.sdf.gz"
LOCAL_GZ    = os.path.join(RAW_DIR, "ChEBI_complete.sdf.gz")
OUTPUT_JSON = os.path.join(RAW_DIR, "chebi.json")


# -----------------------------------------------------------------------------
# 3) Download the SDF if it's not already present
# -----------------------------------------------------------------------------
if not os.path.exists(LOCAL_GZ):
    print("Downloading ChEBI_complete.sdf.gz (≈136 MB)...")
    urlretrieve(SDF_URL, LOCAL_GZ)
    print("Download complete.")


# -----------------------------------------------------------------------------
# 4) Stream through the gzipped SDF & build unified records
# -----------------------------------------------------------------------------
print("Opening gzipped SDF via ForwardSDMolSupplier (this may take a minute)...")
start_time = time.time()
records = []

with gzip.open(LOCAL_GZ, "rb") as gz_f:
    # disable internal sanitization, we'll do it per-molecule
    suppl = Chem.ForwardSDMolSupplier(gz_f, sanitize=False, removeHs=False)

    for idx, mol in enumerate(suppl):
        if mol is None:
            continue

        # try to sanitize this molecule; skip on failure
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            continue

        # helper to pick the first available property
        def get_prop(*keys):
            for k in keys:
                if mol.HasProp(k):
                    return mol.GetProp(k)
            return None

        cid      = get_prop("ChEBI ID", "CHEBI ID", "ChEBI_ID")
        name     = mol.GetProp("_Name") if mol.HasProp("_Name") else None

        # safe SMILES conversion
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception:
            smiles = None

        formula  = get_prop("Formulae", "Molecular Formula", "formula")
        inchi    = get_prop("InChI")
        inchikey = get_prop("InChIKey", "InChI_KEY")
        mass     = get_prop("Mass")

        rec = {
            "id":         cid or f"chebi_{idx}",
            "source":     "chebi",
            "source_id":  cid,
            "name":       name,
            "type":       "molecule",
            "formula":    formula,
            "smiles":     smiles,
            "inchikey":   inchikey,
            "properties": {
                "inchi": inchi,
                "mass":  float(mass) if mass and mass.replace('.','',1).isdigit() else None
            }
        }
        # after building `rec` but before appending:
        # 1) try synonyms if no name
        if not rec['name']:
            syn = get_prop("ChEBI Synonyms")
            if syn:
                # first synonym is usually the best trivial name
                rec['name'] = syn.split("|")[0]
        # 2) fallback to formula or ID
        if not rec['name']:
            rec['name'] = rec['formula'] or rec['source_id']

        records.append(rec)

        # log progress every 100 000 molecules
        if (idx + 1) % 100_000 == 0:
            elapsed = time.time() - start_time
            pct = (idx + 1) / len(suppl) * 100 if hasattr(suppl, "__len__") else ""
            print(f"[{elapsed:.0f}s] Parsed {idx+1} molecules {pct and f'({pct:.1f}% )'}...")

# -----------------------------------------------------------------------------
# 5) Write out JSON
# -----------------------------------------------------------------------------
duration = time.time() - start_time
print(f"Writing {len(records)} records to {OUTPUT_JSON} ({duration:.0f}s)...")
with open(OUTPUT_JSON, "w") as out_f:
    json.dump(records, out_f, indent=2)

print("Done.")
