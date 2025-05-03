#!/usr/bin/env python3
"""
scripts/standardize.py

Standardize raw material data and categorize entries into multiple JSON buckets,
including sanitization of SMILES strings.

USAGE:
    python scripts/standardize.py

CUSTOM CATEGORIES:
- Categories are defined in the `categories` dict at the top.
- To add or modify a category:
    1. Add a new key to the `categories` dict with an empty list as its value.
    2. Under the "# Classification logic" section, insert logic to append `rec` to your new category when appropriate.
- Common record fields available for classification:
    - `rec['source']`: original source name (e.g., 'pubchem', 'materialsproject')
    - `rec['composition']`: dict of element symbols → stoichiometric counts
    - `rec['nelements']`: number of unique elements
    - `rec['type']`: 'molecule', 'crystal', or 'unknown'
    - `rec['formula']`: formula string
    - `rec['properties']`: computed properties (e.g., 'band_gap')
    - `rec['smiles']`: SMILES string (now sanitized)
"""

import os
import glob
import json
import time

from rdkit import RDLogger, Chem
from pymatgen.core.periodic_table import Element

# -----------------------------------------------------------------------------
# Disable RDKit logs so sanitization failures don’t spam the console
# -----------------------------------------------------------------------------
RDLogger.DisableLog('rdApp.*')

# --- Configurable reporting settings ---
REPORT_INTERVAL = 30          # seconds between time-based updates
REPORT_EVERY    = 1000        # records between count-based updates



# --- 1) Prepare folders ---
RAW_DIR   = os.path.abspath(os.path.join(__file__, '..', '../data/raw'))
CLEAN_DIR = os.path.abspath(os.path.join(__file__, '..', '../data/cleaned'))
os.makedirs(CLEAN_DIR, exist_ok=True)

# --- 2) Load & unify records from every JSON in raw/ ---
records = []
for fp in glob.glob(os.path.join(RAW_DIR, '*.json')):

    with open(fp, 'r') as f:
        data = json.load(f)
    entries = data if isinstance(data, list) else [data]
    for idx, e in enumerate(entries):
        rec = {
            'source':      e.get('source') or os.path.splitext(os.path.basename(fp))[0],
            'source_file': os.path.basename(fp),
            'id':          None,
            'formula':     None,
            'type':        None,
            'composition': {},
            'nelements':   None,
            'smiles':      None,
            'inchikey':    None,
            'structure':   None,
            'sites':       None,
            'builder_meta':None,
            'properties':  None,
        }
        # Assign unique ID
        if 'material_id' in e:
            rec['id'] = e['material_id']
        elif 'id' in e:
            rec['id'] = str(e['id'])
        elif 'inchikey' in e:
            rec['id'] = e['inchikey']
        else:
            rec['id'] = f"{rec['source_file']}_{idx}"

        # Formula
        rec['formula'] = e.get('formula_pretty') or e.get('formula') or None

        # Type inference
        if 'smiles' in e or 'inchikey' in e:
            rec['type'] = 'molecule'
        elif 'structure' in e or 'lattice' in e:
            rec['type'] = 'crystal'
        else:
            rec['type'] = 'unknown'

        # Composition & element count
        comp = e.get('composition_reduced') or e.get('composition') or {}
        rec['composition'] = comp
        rec['nelements']   = len(comp) if isinstance(comp, dict) else None

        # Other fields
        rec['smiles']       = e.get('smiles')
        rec['inchikey']     = e.get('inchikey')
        rec['structure']    = e.get('structure') or e.get('lattice')
        rec['sites']        = e.get('sites')
        rec['builder_meta'] = e.get('builder_meta')
        rec['properties']   = e.get('properties')

        # -----------------------------------------------------------------------------
        # SANITIZATION: validate SMILES strings via RDKit; drop if invalid
        # -----------------------------------------------------------------------------
        if rec['type'] == 'molecule' and rec['smiles']:
            try:
                m = Chem.MolFromSmiles(rec['smiles'], sanitize=True)
                if m is None:
                    rec['smiles'] = None
                    rec['inchikey'] = None
                else:
                    # Optionally, rewrite to a canonical SMILES
                    rec['smiles'] = Chem.MolToSmiles(m, canonical=True)
            except Exception:
                rec['smiles'] = None
                rec['inchikey'] = None

        records.append(rec)


# --- 3) Define all desired categories ---
categories = {
    # By number of elements
    'simple_elements': [], 'binary_compounds': [], 'ternary_compounds': [], 'quaternary_compounds': [],
    # Broad chemical classes
    'metals_intermetallics': [], 'ceramics_oxides': [], 'ceramics_nitrides': [], 'ceramics_carbides': [],
    'ceramics_borides': [], 'ceramics_sulfides': [], 'ceramics_selenides': [], 'ceramics_tellurides': [],
    'ceramics_halides': [], 'covalent_network_solids': [], 'molecular_crystals': [],
    'metal_organic_frameworks': [], 'two_dimensional_materials': [],
    # Electronic / functional classes
    'semiconductors': [], 'insulators_dielectrics': [], 'conductors_metals': [], 'superconductors': [],
    'magnetic_materials': [], 'piezoelectrics_ferroelectrics': [], 'thermoelectrics': [],
    'catalysts': [], 'battery_materials': [], 'battery_cathodes': [], 'battery_anodes': [],
    'battery_solid_electrolytes': [],
    # Structural families (stubs)
    'perovskites': [], 'spinels': [], 'wurtzite_zincblende': [], 'rocksalt': [], 'fluorites': [],
    'heusler_laves': [], 'zeolites': [],
    # Application / use-case (stubs)
    'aerospace_alloys': [], 'biomedical_materials': [], 'polymers': [], 'glasses': [],
}

# --- Progress tracking setup ---
total = len(records)
processed = 0
start_time = time.time()
last_time_report = start_time

# --- 4) Classification logic (unchanged but now working on sanitized records) ---
for rec in records:
    comp_keys = set(rec['composition'].keys())
    ne = rec['nelements'] or 0
    props = rec.get('properties') or {}

    # --- Number of elements ---
    if ne == 1:   categories['simple_elements'].append(rec)
    elif ne == 2: categories['binary_compounds'].append(rec)
    elif ne == 3: categories['ternary_compounds'].append(rec)
    elif ne >= 4: categories['quaternary_compounds'].append(rec)

    # --- Broad chemical classes ---
    if any(Element(el).is_metal for el in comp_keys):
        categories['metals_intermetallics'].append(rec)
    if rec['type'] == 'crystal':
        if 'O'  in comp_keys: categories['ceramics_oxides'].append(rec)
        if 'N'  in comp_keys: categories['ceramics_nitrides'].append(rec)
        if 'C' in comp_keys and 'H' not in comp_keys:
            categories['ceramics_carbides'].append(rec)
        if 'B'  in comp_keys: categories['ceramics_borides'].append(rec)
        if 'S'  in comp_keys: categories['ceramics_sulfides'].append(rec)
        if 'Se' in comp_keys: categories['ceramics_selenides'].append(rec)
        if 'Te' in comp_keys: categories['ceramics_tellurides'].append(rec)
        if comp_keys & {'F','Cl','Br','I'}:
            categories['ceramics_halides'].append(rec)
        if ne == 1 and rec['formula'] in {'C','Si'} or comp_keys == {'Si','C'}:
            categories['covalent_network_solids'].append(rec)
        if ne > 1 and all(not Element(el).is_metal and not Element(el).is_metalloid
                         for el in comp_keys):
            categories['molecular_crystals'].append(rec)
    if rec['type'] == 'crystal' and any(Element(el).is_metal for el in comp_keys)\
                                  and 'C' in comp_keys and 'O' in comp_keys:
        categories['metal_organic_frameworks'].append(rec)
    # 2D, structural families, application stubs omitted for brevity…

    # --- Electronic / functional classes ---
    if 'band_gap' in props and props['band_gap'] is not None:
        bg = props['band_gap']
        if 0 < bg < 3.5: categories['semiconductors'].append(rec)
        elif bg >= 3.5: categories['insulators_dielectrics'].append(rec)
    if props.get('superconductor'):     categories['superconductors'].append(rec)
    if props.get('magnetic_type'):      categories['magnetic_materials'].append(rec)
    if props.get('piezoelectric_tensor') or props.get('polarization'):
        categories['piezoelectrics_ferroelectrics'].append(rec)
    if props.get('seebeck_coefficient'): categories['thermoelectrics'].append(rec)
    if props.get('catalytic_activity'):  categories['catalysts'].append(rec)
    if 'Li' in comp_keys or 'Na' in comp_keys:
        categories['battery_materials'].append(rec)
        if 'Li' in comp_keys and 'O' in comp_keys and any(Element(el).is_metal and el!='Li'
                                                           for el in comp_keys):
            categories['battery_cathodes'].append(rec)
        if 'C' in comp_keys and 'Li' in comp_keys:
            categories['battery_anodes'].append(rec)
        if rec['type']=='crystal' and 'Li' in comp_keys and comp_keys&{'S','P','Cl','Br','I'}:
            categories['battery_solid_electrolytes'].append(rec)

    # --- Progress reporting ---
    processed += 1
    now = time.time()
    if processed % REPORT_EVERY == 0:
        pct = (processed/total)*100
        left = total-processed
        print(f"Processed {processed}/{total} ({pct:.1f}%), {left} remaining", flush=True)
    if now - last_time_report >= REPORT_INTERVAL:
        elapsed = now - start_time
        mins, secs = divmod(int(elapsed), 60)
        pct = (processed/total)*100
        left = total-processed
        print(f"[{mins}m{secs}s] Processed {processed}/{total} ({pct:.1f}%), {left} remaining",
              flush=True)
        last_time_report = now

# --- 5) Dump one JSON per category ---
for cat, recs in categories.items():
    out_fp = os.path.join(CLEAN_DIR, f"{cat}.json")
    with open(out_fp, 'w') as fw:
        json.dump(recs, fw, indent=2)

print(f"Done: standardized {total} records into {len(categories)} categories.")
