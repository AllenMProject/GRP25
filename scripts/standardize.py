#!/usr/bin/env python3
"""
scripts/standardize.py

Standardize raw material data and categorize entries into multiple JSON buckets,
including sanitization of SMILES strings.

USAGE:
    python scripts/standardize.py
"""

import os
import glob
import json
import time
import ijson
from decimal import Decimal
from rdkit import RDLogger, Chem
from pymatgen.core.periodic_table import Element

# Disable RDKit logs
RDLogger.DisableLog('rdApp.*')

# Configurable reporting settings
REPORT_INTERVAL = 30  # seconds between time-based updates
REPORT_EVERY = 1000   # records between count-based updates

# Prepare folders
RAW_DIR = os.path.abspath(os.path.join(__file__, '..', '../data/raw'))
CLEAN_DIR = os.path.abspath(os.path.join(__file__, '..', '../data/cleaned'))
print("Cleaning into:", CLEAN_DIR)
os.makedirs(CLEAN_DIR, exist_ok=True)

# Custom JSON encoder to handle Decimal
class DecimalEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Decimal):
            return float(obj)  # Convert Decimal to float
        return super().default(obj)

# Function to process a single record
def process_record(entry, source_file, idx):
    rec = {
        'source': entry.get('source') or os.path.splitext(os.path.basename(source_file))[0],
        'source_file': os.path.basename(source_file),
        'id': None,
        'formula': None,
        'type': None,
        'composition': {},
        'nelements': None,
        'smiles': None,
        'inchikey': None,
        'structure': None,
        'sites': None,
        'builder_meta': None,
        'properties': None,
    }
    # Assign unique ID
    if 'material_id' in entry:
        rec['id'] = entry['material_id']
    elif 'id' in entry:
        rec['id'] = str(entry['id'])
    elif 'inchikey' in entry:
        rec['id'] = entry['inchikey']
    else:
        rec['id'] = f"{rec['source_file']}_{idx}"

    # Formula
    rec['formula'] = entry.get('formula_pretty') or entry.get('formula') or None

    # Type inference
    if 'smiles' in entry or 'inchikey' in entry:
        rec['type'] = 'molecule'
    elif 'structure' in entry or 'lattice' in entry:
        rec['type'] = 'crystal'
    else:
        rec['type'] = 'unknown'

    # Composition & element count
    comp = entry.get('composition_reduced') or entry.get('composition') or {}
    rec['composition'] = comp
    rec['nelements'] = len(comp) if isinstance(comp, dict) else None

    # Other fields
    rec['smiles'] = entry.get('smiles')
    rec['inchikey'] = entry.get('inchikey')
    rec['structure'] = entry.get('structure') or entry.get('lattice')
    rec['sites'] = entry.get('sites')
    rec['builder_meta'] = entry.get('builder_meta')
    rec['properties'] = entry.get('properties')

    # Sanitize SMILES
    if rec['type'] == 'molecule' and rec['smiles']:
        try:
            m = Chem.MolFromSmiles(rec['smiles'], sanitize=True)
            if m is None:
                rec['smiles'] = rec['inchikey'] = None
            else:
                rec['smiles'] = Chem.MolToSmiles(m, canonical=True)
        except Exception:
            rec['smiles'] = rec['inchikey'] = None

    return rec

# Define categories
categories = {
    'simple_elements': [],
    'binary_compounds': [],
    'ternary_compounds': [],
    'quaternary_compounds': [],
    'metals_intermetallics': [],
    'ceramics_oxides': [],
    'ceramics_nitrides': [],
    'ceramics_carbides': [],
    'ceramics_borides': [],
    'ceramics_sulfides': [],
    'ceramics_selenides': [],
    'ceramics_tellurides': [],
    'ceramics_halides': [],
    'covalent_network_solids': [],
    'molecular_crystals': [],
    'metal_organic_frameworks': [],
    'two_dimensional_materials': [],
    'semiconductors': [],
    'insulators_dielectrics': [],
    'conductors_metals': [],
    'superconductors': [],
    'magnetic_materials': [],
    'piezoelectrics_ferroelectrics': [],
    'thermoelectrics': [],
    'catalysts': [],
    'battery_materials': [],
    'battery_cathodes': [],
    'battery_anodes': [],
    'battery_solid_electrolytes': [],
    'perovskites': [],
    'spinels': [],
    'wurtzite_zincblende': [],
    'rocksalt': [],
    'fluorites': [],
    'heusler_laves': [],
    'zeolites': [],
    'aerospace_alloys': [],
    'biomedical_materials': [],
    'polymers': [],
    'glasses': [],
    'uncategorized': [],  # Fallback category
}

# List to store all records for master.json
all_records = []

# Track category counts for debugging
category_counts = {cat: 0 for cat in categories}

# Progress tracking setup
processed = 0
start_time = time.time()
last_time_report = start_time

# Count total records for progress reporting
total = 0
for fp in glob.glob(os.path.join(RAW_DIR, '*.json')):
    print(fp, " being initialized...")

    with open(fp, 'rb') as f:
        data = json.load(f)
        total += len(data) if isinstance(data, list) else 1
    print(fp, " initialization Complete...")
if total == 0:
    print("⚠️  No records found in RAW_DIR!")
    total = float('inf')

# Process records incrementally using streaming
for fp in glob.glob(os.path.join(RAW_DIR, '*.json')):
    print(fp, " being loaded...")
    with open(fp, 'rb') as f:
        parser = ijson.items(f, 'item')  # Adjust prefix if needed
        for idx, entry in enumerate(parser):
            rec = process_record(entry, fp, idx)
            comp_keys = set(rec['composition'].keys())
            ne = rec['nelements'] or 0
            props = rec.get('properties') or {}

            # Add to all_records for master.json
            all_records.append(rec)

            # Track whether this record gets categorized
            assigned = False

            # Classification logic
            # Number of elements
            if ne == 1:
                categories['simple_elements'].append(rec)
                category_counts['simple_elements'] += 1
                assigned = True
            elif ne == 2:
                categories['binary_compounds'].append(rec)
                category_counts['binary_compounds'] += 1
                assigned = True
            elif ne == 3:
                categories['ternary_compounds'].append(rec)
                category_counts['ternary_compounds'] += 1
                assigned = True
            elif ne >= 4:
                categories['quaternary_compounds'].append(rec)
                category_counts['quaternary_compounds'] += 1
                assigned = True

            # Broad chemical classes
            if any(Element(el).is_metal for el in comp_keys):
                categories['metals_intermetallics'].append(rec)
                category_counts['metals_intermetallics'] += 1
                assigned = True

            if rec['type'] == 'crystal':
                if 'O' in comp_keys:
                    categories['ceramics_oxides'].append(rec)
                    category_counts['ceramics_oxides'] += 1
                    assigned = True
                if 'N' in comp_keys:
                    categories['ceramics_nitrides'].append(rec)
                    category_counts['ceramics_nitrides'] += 1
                    assigned = True
                if 'C' in comp_keys and 'H' not in comp_keys:
                    categories['ceramics_carbides'].append(rec)
                    category_counts['ceramics_carbides'] += 1
                    assigned = True
                if 'B' in comp_keys:
                    categories['ceramics_borides'].append(rec)
                    category_counts['ceramics_borides'] += 1
                    assigned = True
                if 'S' in comp_keys:
                    categories['ceramics_sulfides'].append(rec)
                    category_counts['ceramics_sulfides'] += 1
                    assigned = True
                if 'Se' in comp_keys:
                    categories['ceramics_selenides'].append(rec)
                    category_counts['ceramics_selenides'] += 1
                    assigned = True
                if 'Te' in comp_keys:
                    categories['ceramics_tellurides'].append(rec)
                    category_counts['ceramics_tellurides'] += 1
                    assigned = True
                if comp_keys & {'F', 'Cl', 'Br', 'I'}:
                    categories['ceramics_halides'].append(rec)
                    category_counts['ceramics_halides'] += 1
                    assigned = True
                if ne == 1 and rec['formula'] in {'C', 'Si'} or comp_keys == {'Si', 'C'}:
                    categories['covalent_network_solids'].append(rec)
                    category_counts['covalent_network_solids'] += 1
                    assigned = True
                if ne > 1 and all(not Element(el).is_metal and not Element(el).is_metalloid for el in comp_keys):
                    categories['molecular_crystals'].append(rec)
                    category_counts['molecular_crystals'] += 1
                    assigned = True

            if rec['type'] == 'crystal' and any(Element(el).is_metal for el in comp_keys) and 'C' in comp_keys and 'O' in comp_keys:
                categories['metal_organic_frameworks'].append(rec)
                category_counts['metal_organic_frameworks'] += 1
                assigned = True

            # Electronic / functional classes
            if 'band_gap' in props and props['band_gap'] is not None:
                bg = props['band_gap']
                if 0 < bg < 3.5:
                    categories['semiconductors'].append(rec)
                    category_counts['semiconductors'] += 1
                    assigned = True
                elif bg >= 3.5:
                    categories['insulators_dielectrics'].append(rec)
                    category_counts['insulators_dielectrics'] += 1
                    assigned = True
            else:
                if processed < 10:  # Log for first few records to avoid spam
                    print(f"Record {rec['id']} missing 'band_gap' in properties")

            if props.get('superconductor'):
                categories['superconductors'].append(rec)
                category_counts['superconductors'] += 1
                assigned = True
            else:
                if processed < 10:
                    print(f"Record {rec['id']} missing or false 'superconductor' in properties: {props.get('superconductor')}")

            if props.get('magnetic_type'):
                categories['magnetic_materials'].append(rec)
                category_counts['magnetic_materials'] += 1
                assigned = True
            else:
                if processed < 10:
                    print(f"Record {rec['id']} missing 'magnetic_type' in properties")

            if props.get('piezoelectric_tensor') or props.get('polarization'):
                categories['piezoelectrics_ferroelectrics'].append(rec)
                category_counts['piezoelectrics_ferroelectrics'] += 1
                assigned = True

            if props.get('seebeck_coefficient'):
                categories['thermoelectrics'].append(rec)
                category_counts['thermoelectrics'] += 1
                assigned = True

            if props.get('catalytic_activity'):
                categories['catalysts'].append(rec)
                category_counts['catalysts'] += 1
                assigned = True

            if 'Li' in comp_keys or 'Na' in comp_keys:
                categories['battery_materials'].append(rec)
                category_counts['battery_materials'] += 1
                assigned = True
                if 'Li' in comp_keys and 'O' in comp_keys and any(Element(el).is_metal and el != 'Li' for el in comp_keys):
                    categories['battery_cathodes'].append(rec)
                    category_counts['battery_cathodes'] += 1
                    assigned = True
                if 'C' in comp_keys and 'Li' in comp_keys:
                    categories['battery_anodes'].append(rec)
                    category_counts['battery_anodes'] += 1
                    assigned = True
                if rec['type'] == 'crystal' and 'Li' in comp_keys and comp_keys & {'S', 'P', 'Cl', 'Br', 'I'}:
                    categories['battery_solid_electrolytes'].append(rec)
                    category_counts['battery_solid_electrolytes'] += 1
                    assigned = True

            # Structural families and application categories
            # Simplified logic for previously empty stub categories
            if rec['type'] == 'crystal':
                # Perovskites: General formula ABO3
                if ne == 3 and 'O' in comp_keys and len([el for el in comp_keys if Element(el).is_metal]) == 2:
                    categories['perovskites'].append(rec)
                    category_counts['perovskites'] += 1
                    assigned = True
                # Spinels: General formula AB2O4
                if ne == 3 and 'O' in comp_keys and rec['composition'].get('O', 0) == 4:
                    categories['spinels'].append(rec)
                    category_counts['spinels'] += 1
                    assigned = True
                # Wurtzite/Zincblende: e.g., ZnS, GaN
                if ne == 2 and comp_keys & {'Zn', 'Ga', 'Al'} and comp_keys & {'S', 'N'}:
                    categories['wurtzite_zincblende'].append(rec)
                    category_counts['wurtzite_zincblende'] += 1
                    assigned = True
                # Rocksalt: e.g., NaCl, MgO
                if ne == 2 and comp_keys & {'Na', 'Mg', 'Ca'} and comp_keys & {'Cl', 'O'}:
                    categories['rocksalt'].append(rec)
                    category_counts['rocksalt'] += 1
                    assigned = True
                # Fluorites: e.g., CaF2
                if ne == 2 and 'F' in comp_keys and comp_keys & {'Ca', 'Zr'}:
                    categories['fluorites'].append(rec)
                    category_counts['fluorites'] += 1
                    assigned = True
                # Heusler/Laves: e.g., Cu2MnAl
                if ne >= 3 and any(Element(el).is_transition_metal for el in comp_keys):
                    categories['heusler_laves'].append(rec)
                    category_counts['heusler_laves'] += 1
                    assigned = True
                # Zeolites: Si/Al with O
                if comp_keys & {'Si', 'Al'} and 'O' in comp_keys:
                    categories['zeolites'].append(rec)
                    category_counts['zeolites'] += 1
                    assigned = True
                # Two-dimensional materials: e.g., Graphene, MoS2
                if ne <= 2 and rec['structure'] and 'layered' in str(rec['structure']).lower():
                    categories['two_dimensional_materials'].append(rec)
                    category_counts['two_dimensional_materials'] += 1
                    assigned = True

            # Application/use-case categories
            # Aerospace alloys: High metal content, e.g., Ti, Al
            if any(el in comp_keys for el in {'Ti', 'Al', 'Ni'}) and all(Element(el).is_metal for el in comp_keys):
                categories['aerospace_alloys'].append(rec)
                category_counts['aerospace_alloys'] += 1
                assigned = True
            # Biomedical materials: e.g., Ti for implants
            if 'Ti' in comp_keys:
                categories['biomedical_materials'].append(rec)
                category_counts['biomedical_materials'] += 1
                assigned = True
            # Polymers: Organic molecules with C, H
            if 'C' in comp_keys and 'H' in comp_keys and rec['type'] == 'molecule':
                categories['polymers'].append(rec)
                category_counts['polymers'] += 1
                assigned = True
            # Glasses: Amorphous SiO2 or similar
            if comp_keys & {'Si', 'B'} and 'O' in comp_keys and not rec['structure']:
                categories['glasses'].append(rec)
                category_counts['glasses'] += 1
                assigned = True

            # Fallback for uncategorized records
            if not assigned:
                categories['uncategorized'].append(rec)
                category_counts['uncategorized'] += 1

            # Progress reporting
            processed += 1
            now = time.time()
            if processed % REPORT_EVERY == 0:
                pct = (processed / total) * 100 if total != float('inf') else 0
                left = total - processed if total != float('inf') else 'unknown'
                print(f"Processed {processed}/{total} records ({pct:.1f}%), {left} remaining", flush=True)
            if now - last_time_report >= REPORT_INTERVAL:
                elapsed = now - start_time
                mins, secs = divmod(int(elapsed), 60)
                print(f"[{mins}m{secs}s] Processed {processed} records", flush=True)
                last_time_report = now
    print(fp, " Load Complete...")
# Dump one JSON per category
for cat, recs in categories.items():
    out_fp = os.path.join(CLEAN_DIR, f"{cat}.json")
    with open(out_fp, 'w') as fw:
        json.dump(recs, fw, indent=2, cls=DecimalEncoder)
    print(f"  • {cat}.json → {len(recs)} records")

# Dump master.json with all records
print(f"About to write master.json with {len(all_records)} total records…")
master_fp = os.path.join(CLEAN_DIR, "master.json")
with open(master_fp, 'w') as fw:
    json.dump(all_records, fw, indent=2, cls=DecimalEncoder)
print(f"✅ Wrote master.json → {master_fp}")

# Print category distribution for debugging
print("\nCategory Distribution:")
for cat, count in category_counts.items():
    print(f"{cat}: {count} records")

print(f"\nDone: standardized {processed} records into {len(categories)} categories.")