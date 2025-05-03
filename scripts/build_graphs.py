#!/usr/bin/env python3
"""
scripts/build_graphs.py

Stage 1 of pipeline: Stream each record from cleaned/master.json,
build its NetworkX graph (with explicit dummy‑atom attachment points),
and dump three pickle files:

  data/graphs/mol_graphs.pkl
  data/graphs/crystal_graphs.pkl
  data/graphs/with_dummies.pkl   ← materials that had at least one dummy atom
"""

import sys, time, pickle
from pathlib import Path
from decimal import Decimal

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')  # silence RDKit warnings

import ijson
import networkx as nx
import rdkit.Chem as Chem
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure

# —————————————————————————————————————————————————————————————
# Helper to coerce Decimal→float, lists, dicts recursively
# —————————————————————————————————————————————————————————————
def to_native(x):
    if isinstance(x, Decimal):
        return float(x)
    if isinstance(x, dict):
        return {k: to_native(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [to_native(v) for v in x]
    return x

# —————————————————————————————————————————————————————————————
# Paths & Constants
# —————————————————————————————————————————————————————————————
ROOT         = Path(__file__).resolve().parent.parent
PULLED_JSON  = ROOT / "data" / "cleaned" / "master.json" # Pulls all the materials from the master.json file
# PULLED_JSON  = ROOT / "data" / "cleaned" / "rocksalt.json" # This can be changed to any specific JSON from cleaned
OUT_DIR      = ROOT / "data" / "graphs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
REPORT_EVERY = 1000
MAX_RECORDS  = None  # set to int for testing

# —————————————————————————————————————————————————————————————
# Build a single NX graph from one record
# —————————————————————————————————————————————————————————————
def build_graph(rec):
    G   = nx.Graph()
    src = rec.get("source", "unknown")

    if rec["type"] == "molecule" and rec.get("smiles"):
        mol = Chem.MolFromSmiles(rec["smiles"])
        if mol is None:
            raise ValueError(f"Invalid SMILES for {rec['id']}: {rec['smiles']!r}")
        for atom in mol.GetAtoms():
            idx      = atom.GetIdx()
            atomic_n = atom.GetAtomicNum()
            is_dummy = (atomic_n == 0)
            G.add_node(
                idx,
                element     = atom.GetSymbol() or "*",
                atomic_num  = atomic_n,
                is_dummy    = is_dummy,
                aromatic    = atom.GetIsAromatic(),
                source      = src
            )
        for b in mol.GetBonds():
            u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            G.add_edge(u, v,
                       order      = str(b.GetBondType()),
                       source     = src)

    elif rec["type"] == "crystal":
        struct_raw   = rec.get("structure", {})
        lattice_dict = struct_raw.get("lattice", struct_raw)
        sites_raw    = struct_raw.get("sites", [])

        lattice_dict = to_native(lattice_dict)
        try:
            lat = Lattice.from_dict(lattice_dict)
        except Exception as e:
            print(f"\n⚠️  Lattice.from_dict failed for {rec['id']}: {e}")
            for k, v in list(lattice_dict.items())[:5]:
                print(f"    - {k}: {type(v).__name__} -> {repr(v)[:50]}")
            raise

        species, fracs = [], []
        for s in sites_raw:
            sp = s["species"]
            if isinstance(sp, list):
                sp0 = sp[0]
                sp  = sp0["element"] if isinstance(sp0, dict) else sp0
            el = sp if isinstance(sp, str) else sp.get("element","?")
            species.append(el)

            raw = s.get("frac_coords") or s.get("abc") or s.get("xyz")
            if raw is None:
                raise ValueError(f"No coords for site in {rec['id']}: {s}")
            fracs.append(to_native(raw))

        struct = Structure(lat, species, fracs)
        for i, site in enumerate(struct.sites):
            G.add_node(i,
                       element     = site.species_string,
                       frac_coords = list(site.frac_coords),
                       is_dummy    = False,
                       source      = src)
        n = len(struct.sites)
        for i in range(n):
            for j in range(i+1, n):
                d = struct.get_distance(i, j)
                if d <= 3.0:
                    G.add_edge(i, j,
                               distance = round(d,4),
                               source   = src)
    else:
        raise ValueError(f"Unknown record type for {rec['id']}: {rec['type']}")

    return G

# —————————————————————————————————————————————————————————————
# Main: stream JSON, build graphs, detect dummies, dump pickles
# —————————————————————————————————————————————————————————————
def main():
    mol_graphs      = []  # list of (id, G, props)
    crystal_graphs  = []  # list of (id, G, props)
    with_dummies    = []  # subset that had dummy atoms

    total_dummy_atoms        = 0
    materials_with_dummy_cnt = 0
    count = 0
    t0    = time.time()

    print(f"⏳ Streaming {'all' if MAX_RECORDS is None else MAX_RECORDS} records from {PULLED_JSON}")
    with open(PULLED_JSON, "rb") as f:
        for rec in ijson.items(f, "item"):
            if MAX_RECORDS and count >= MAX_RECORDS:
                break
            mid = rec.get("id", f"unknown_{count}")

            # Build the graph
            try:
                G = build_graph(rec)
            except Exception as e:
                print(f"⚠️  Skipping {mid}: {e}")
                continue

            # Prepare material‐level properties
            props = {
                "formula":   rec.get("formula"),
                "source":    rec.get("source"),
                **to_native(rec.get("properties") or {})
            }

            # Count dummy atoms
            dummy_nodes = [n for n, attrs in G.nodes(data=True) if attrs.get("is_dummy")]
            if dummy_nodes:
                materials_with_dummy_cnt += 1
                total_dummy_atoms        += len(dummy_nodes)
                with_dummies.append((mid, G, props))

            # Append to the right list
            if rec["type"] == "molecule":
                mol_graphs.append((mid, G, props))
            else:
                crystal_graphs.append((mid, G, props))

            count += 1
            if count % REPORT_EVERY == 0:
                print(f"  >> Built {count} graphs so far…")

    dt = time.time() - t0
    print(f"\n✅ Built {count} graphs in {dt:.1f}s")
    print(f"   ↳ molecules     : {len(mol_graphs):,}")
    print(f"   ↳ crystals      : {len(crystal_graphs):,}")
    print(f"   ↳ with dummies  : {materials_with_dummy_cnt} materials, total {total_dummy_atoms} dummy atoms")

    # Dump pickles
    with open(OUT_DIR/"mol_graphs.pkl",     "wb") as fp: pickle.dump(mol_graphs, fp)
    with open(OUT_DIR/"crystal_graphs.pkl", "wb") as fp: pickle.dump(crystal_graphs, fp)
    with open(OUT_DIR/"with_dummies.pkl",   "wb") as fp: pickle.dump(with_dummies, fp)

    print("🎉 All graphs written, including `with_dummies.pkl`.")

if __name__ == "__main__":
    main()
