#!/usr/bin/env python3
# scripts/export_knowledgeGraphs.py

"""
Export a file‐based RDF/OWL knowledge graph from your pickled graphs,
including material‐level props and marking dummy atoms as mt:AttachmentPoint,
and also produce PyTorch‐Geometric `.pt` files for GNN consumption.

Outputs:
  • Turtle TTL → exports/KnowledgeGraphs/materials.ttl
  • PyG `.pt` for molecules → exports/KnowledgeGraphs/kg_mol_graphs.pt
  • PyG `.pt` for crystals  → exports/KnowledgeGraphs/kg_crystal_graphs.pt
"""

import pickle
from pathlib import Path
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF, XSD, OWL

import numpy as np
import torch
from torch_geometric.data import Data

# ─── Configuration ───────────────────────────────────────────
ROOT       = Path(__file__).resolve().parent.parent
GRAPHS_DIR = ROOT / "data" / "graphs"
MOL_PKL    = GRAPHS_DIR / "mol_graphs.pkl"      # list of (mid, G, props)
CRYS_PKL   = GRAPHS_DIR / "crystal_graphs.pkl"
OUT_KG_DIR = ROOT / "exports" / "KnowledgeGraphs"
OUT_TTL    = OUT_KG_DIR / "materials.ttl"
OUT_PT_M   = OUT_KG_DIR / "kg_mol_graphs.pt"
OUT_PT_C   = OUT_KG_DIR / "kg_crystal_graphs.pt"

MT = Namespace("http://materialsdb.org/schema#")
MB = Namespace("http://materialsdb.org/data/")

# ─── Helpers for RDF build ────────────────────────────────────
def material_uri(mid): return MB[f"material/{mid}"]
def atom_uri(mid, idx): return MB[f"atom/{mid}_{idx}"]
def bond_uri(mid, u, v):  return MB[f"bond/{mid}_{u}_{v}"]

def add_props(kg, subj, props):
    for k, v in props.items():
        if isinstance(v, bool):
            lit = Literal(v, datatype=XSD.boolean)
        elif isinstance(v, (int, float)):
            lit = Literal(v, datatype=XSD.double)
        else:
            lit = Literal(v, datatype=XSD.string)
        kg.add((subj, MT[k], lit))

# ─── Helpers for GNN export ──────────────────────────────────
ELEMENTS = {
    "H":0, "C":1, "N":2, "O":3, "F":4,
    "Si":5, "P":6, "S":7, "Cl":8, "*":9
}
NUM_ELEMENTS = len(ELEMENTS)

def atom_feat(attrs):
    el = attrs.get("element", "*")
    idx = ELEMENTS.get(el, ELEMENTS["*"])
    onehot = np.zeros(NUM_ELEMENTS, dtype=np.float32)
    onehot[idx] = 1.0
    extras = np.array([
        attrs.get("atomic_num", 0),
        1.0 if attrs.get("is_dummy") else 0.0
    ], dtype=np.float32)
    return np.hstack([onehot, extras])

def graph_to_data(mid, G, props):
    # node features
    feat_list = [atom_feat(attrs) for _, attrs in G.nodes(data=True)]
    if not feat_list:
        return None
    x = torch.from_numpy(np.stack(feat_list, axis=0)).float()

    # edges
    edge_index = []
    edge_attr  = []
    for u, v, attrs in G.edges(data=True):
        order = attrs.get("order", "")
        bond_type = {"SINGLE":1, "DOUBLE":2, "TRIPLE":3}.get(order, 0)
        edge_index += [[u, v], [v, u]]
        edge_attr  += [[bond_type], [bond_type]]

    if edge_index:
        ei = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        ea = torch.tensor(edge_attr, dtype=torch.float)
    else:
        ei = torch.empty((2,0), dtype=torch.long)
        ea = torch.empty((0,1), dtype=torch.float)

    data = Data(x=x, edge_index=ei, edge_attr=ea)
    # attach material props as scalar tensors
    for k, v in (props or {}).items():
        if isinstance(v, (int, float, bool)):
            setattr(data, k, torch.tensor([v]))
    data.mid = mid
    return data

# ─── Main export function ─────────────────────────────────────
def export_kg():
    # ensure output directories exist
    OUT_KG_DIR.mkdir(parents=True, exist_ok=True)

    # --- build RDF/OWL TTL ---
    kg = Graph()
    kg.bind("mt", MT)
    kg.bind("owl", OWL)

    mols = pickle.load(open(MOL_PKL, "rb"))   # each item: (mid, G, props)
    crys = pickle.load(open(CRYS_PKL, "rb"))

    # molecules
    for mid, G, props in mols:
        m = URIRef(material_uri(mid))
        kg.add((m, RDF.type, MT.Material))
        kg.add((m, RDF.type, MT.Molecule))
        add_props(kg, m, props)
        for n, attrs in G.nodes(data=True):
            a = URIRef(atom_uri(mid, n))
            cls = MT.AttachmentPoint if attrs.get("is_dummy") else MT.Atom
            kg.add((a, RDF.type, cls))
            kg.add((m, MT.hasAtom, a))
            kg.add((a, MT.element, Literal(attrs.get("element",""), datatype=XSD.string)))
            if attrs.get("is_dummy"):
                kg.add((a, MT.isDummy, Literal(True, datatype=XSD.boolean)))
        for u, v, attrs in G.edges(data=True):
            b = URIRef(bond_uri(mid, u, v))
            kg.add((b, RDF.type, MT.Bond))
            kg.add((b, MT.connects, URIRef(atom_uri(mid, u))))
            kg.add((b, MT.connects, URIRef(atom_uri(mid, v))))
            kg.add((b, MT.order, Literal(attrs.get("order",""), datatype=XSD.string)))
            kg.add((m, MT.hasBond, b))

    # crystals
    for mid, G, props in crys:
        m = URIRef(material_uri(mid))
        kg.add((m, RDF.type, MT.Material))
        kg.add((m, RDF.type, MT.Crystal))
        add_props(kg, m, props)
        for n, attrs in G.nodes(data=True):
            a = URIRef(atom_uri(mid, n))
            kg.add((a, RDF.type, MT.Site))
            kg.add((m, MT.hasSite, a))
            kg.add((a, MT.element, Literal(attrs.get("element",""), datatype=XSD.string)))
            coord = ",".join(map(str, attrs.get("frac_coords", [])))
            kg.add((a, MT.fracCoords, Literal(coord, datatype=XSD.string)))
        for u, v, attrs in G.edges(data=True):
            b = URIRef(bond_uri(mid, u, v))
            kg.add((b, RDF.type, MT.Contact))
            kg.add((b, MT.connects, URIRef(atom_uri(mid, u))))
            kg.add((b, MT.connects, URIRef(atom_uri(mid, v))))
            kg.add((b, MT.distance, Literal(attrs.get("distance",0.0), datatype=XSD.float)))
            kg.add((m, MT.hasContact, b))

    # write TTL
    print(f"⏳ Writing TTL → {OUT_TTL}")
    kg.serialize(str(OUT_TTL), format="turtle")
    print("✅ TTL export complete.")

    # --- export PyG .pt files ---
    print("⏳ Exporting PyG .pt for molecules…")
    mol_data = [graph_to_data(mid, G, props) for mid, G, props in mols]
    mol_data = [d for d in mol_data if d is not None]
    torch.save(mol_data, str(OUT_PT_M))
    print(f"✅ Saved {len(mol_data)} molecule graphs → {OUT_PT_M}")

    print("⏳ Exporting PyG .pt for crystals…")
    crys_data = [graph_to_data(mid, G, props) for mid, G, props in crys]
    crys_data = [d for d in crys_data if d is not None]
    torch.save(crys_data, str(OUT_PT_C))
    print(f"✅ Saved {len(crys_data)} crystal graphs → {OUT_PT_C}")

if __name__ == "__main__":
    export_kg()
