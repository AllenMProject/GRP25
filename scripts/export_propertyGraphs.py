#!/usr/bin/env python3
"""
scripts/export_propertyGraphs.py

Load your pickled graphs (with props), then export:

  • JSON with material‑level props + node/edge data
  • CSV for nodes & edges (including props)
  • plain edgelists for NetworkX
  • PyTorch‑Geometric `.pt` ready files for GNN training
"""

import json
import csv
import pickle
from pathlib import Path

import networkx as nx
import numpy as np
import torch
from torch_geometric.data import Data

# ─── Paths ────────────────────────────────────────────────────────────
ROOT       = Path(__file__).resolve().parent.parent
GRAPHS_DIR = ROOT / "data" / "graphs"
EXPORT_DIR = ROOT / "exports" / "PropertyGraphs"
EXPORT_DIR.mkdir(parents=True, exist_ok=True)

# ─── Load pickled graphs ─────────────────────────────────────────────
mol_file    = GRAPHS_DIR / "mol_graphs.pkl"
crys_file   = GRAPHS_DIR / "crystal_graphs.pkl"
mol_graphs  = pickle.load(open(mol_file, "rb"))    # list of (mid,G,props)
crys_graphs = pickle.load(open(crys_file, "rb"))

all_graphs = mol_graphs + crys_graphs

# ─── JSON export ─────────────────────────────────────────────────────
json_out = EXPORT_DIR / "materials_graphs.json"
print(f"⏳ Writing JSON → {json_out}")
data = {}
for mid, G, props in all_graphs:
    data[mid] = {
        "props": props,
        "nodes": [{"id": n, **attrs} for n, attrs in G.nodes(data=True)],
        "edges": [{"u": u, "v": v, **attrs} for u, v, attrs in G.edges(data=True)]
    }
with open(json_out, "w") as fp:
    json.dump(data, fp, indent=2)

# ─── CSV nodes ───────────────────────────────────────────────────────
nodes_csv = EXPORT_DIR / "nodes.csv"
print(f"⏳ Writing nodes CSV → {nodes_csv}")
with open(nodes_csv, "w", newline="") as fp:
    w = csv.writer(fp)
    w.writerow(["material_id", "node_id", "element", "props", "node_attrs"])
    for mid, G, props in all_graphs:
        props_json = json.dumps(props)
        for n, attrs in G.nodes(data=True):
            w.writerow([mid, n, attrs.get("element",""), props_json, json.dumps(attrs)])

# ─── CSV edges ───────────────────────────────────────────────────────
edges_csv = EXPORT_DIR / "edges.csv"
print(f"⏳ Writing edges CSV → {edges_csv}")
with open(edges_csv, "w", newline="") as fp:
    w = csv.writer(fp)
    w.writerow(["material_id", "u", "v", "props", "edge_attrs"])
    for mid, G, props in all_graphs:
        props_json = json.dumps(props)
        for u, v, attrs in G.edges(data=True):
            w.writerow([mid, u, v, props_json, json.dumps(attrs)])

# ─── Edge‑list dumps ─────────────────────────────────────────────────
print("⏳ Writing edgelists…")
for mid, G, _ in all_graphs:
    out = EXPORT_DIR / f"{mid}_edgelist.txt"
    nx.write_edgelist(G, str(out), data=False)

# ─── PyTorch‑Geometric export ────────────────────────────────────────

# define a fixed element→index mapping
ELEMENTS = {
    "H":0, "C":1, "N":2, "O":3, "F":4,
    "Si":5, "P":6, "S":7, "Cl":8, "Dummy":9
}
NUM_ELEMENTS = len(ELEMENTS)

def atom_features(attrs):
    # one‑hot element + atomic_num + aromatic + is_dummy
    onehot = np.zeros(NUM_ELEMENTS, dtype=np.float32)
    el = attrs.get("element","")
    idx = ELEMENTS.get(el, ELEMENTS["Dummy"])
    onehot[idx] = 1.0
    extras = np.array([
        attrs.get("atomic_num", 0),
        1.0 if attrs.get("aromatic") else 0.0,
        1.0 if attrs.get("is_dummy") else 0.0
    ], dtype=np.float32)
    return np.concatenate([onehot, extras])

def graph_to_data(mid, G, props=None):
    # 1) Node features
    feats = [atom_features(attrs) for _, attrs in G.nodes(data=True)]
    if not feats:
        return None
    x = torch.from_numpy(np.stack(feats, axis=0)).float()

    # 2) Edges
    edge_index = []
    edge_attr  = []
    for u, v, attrs in G.edges(data=True):
        # bidirectional
        edge_index += [[u, v], [v, u]]
        order = attrs.get("order","")
        bond_type = {"SINGLE":1, "DOUBLE":2, "TRIPLE":3}.get(order, 0)
        edge_attr += [[bond_type], [bond_type]]

    if edge_index:
        ei = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        ea = torch.tensor(edge_attr,  dtype=torch.float)
    else:
        ei = torch.empty((2,0), dtype=torch.long)
        ea = torch.empty((0,1), dtype=torch.float)

    data = Data(x=x, edge_index=ei, edge_attr=ea)

    # 3) Graph‑level props
    if props:
        for k, v in props.items():
            if isinstance(v, (int, float, bool)):
                data[k] = torch.tensor([v])

    data.mid = mid
    return data

print("⏳ Exporting PyG .pt files…")
mol_data = [d for d in (graph_to_data(mid, G, props) for mid, G, props in mol_graphs) if d]
crys_data = [d for d in (graph_to_data(mid, G, props) for mid, G, props in crys_graphs) if d]

pt_mol = EXPORT_DIR / "mol_graphs.pt"
pt_crys = EXPORT_DIR / "crystal_graphs.pt"
torch.save(mol_data,  str(pt_mol))
torch.save(crys_data, str(pt_crys))

print(f"✅ Saved {len(mol_data)} molecule graphs → {pt_mol}")
print(f"✅ Saved {len(crys_data)} crystal  graphs → {pt_crys}")
print("✅ export_propertyGraphs complete.")
