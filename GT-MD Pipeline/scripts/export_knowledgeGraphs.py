#!/usr/bin/env python3
# scripts/export_knowledgeGraphs.py

"""
Export a file‑based RDF/OWL knowledge graph from your pickled graphs,
including material‑level props and marking dummy atoms as mt:AttachmentPoint,
and also produce PyTorch‑Geometric `.pt` files for GNN consumption.

Outputs:
  • Turtle TTL → exports/KnowledgeGraphs/materials.ttl
  • PyG `.pt` for molecules → exports/KnowledgeGraphs/kg_mol_graphs.pt
  • PyG `.pt` for crystals  → exports/KnowledgeGraphs/kg_crystal_graphs.pt
"""

import pickle
from pathlib import Path
import time
import re

import numpy as np
import torch
from torch_geometric.data import Data
from rdflib import Namespace, URIRef, Literal
from rdflib.namespace import RDF, XSD, OWL

t0 = time.time()

# ─── Configuration ───────────────────────────────────────────
ROOT         = Path(__file__).resolve().parent.parent
GRAPHS_DIR   = ROOT / "data" / "graphs"
MOL_PKL      = GRAPHS_DIR / "mol_graphs.pkl"
CRYS_PKL     = GRAPHS_DIR / "crystal_graphs.pkl"
OUT_KG_DIR   = ROOT / "exports" / "KnowledgeGraphs"
OUT_KG_DIR.mkdir(parents=True, exist_ok=True)
OUT_TTL      = OUT_KG_DIR / "materials.ttl"
OUT_PT_M     = OUT_KG_DIR / "kg_mol_graphs.pt"
OUT_PT_C     = OUT_KG_DIR / "kg_crystal_graphs.pt"

MT = Namespace("http://materialsdb.org/schema#")
MB = Namespace("http://materialsdb.org/data/")

# Sanitize IDs to avoid invalid IRIs
_id_clean = lambda s: re.sub(r'[^0-9A-Za-z_]', '_', str(s))
material_uri = lambda mid: MB[f"material/{_id_clean(mid)}"]
atom_uri     = lambda mid, idx: MB[f"atom/{_id_clean(mid)}_{idx}"]
bond_uri     = lambda mid, u, v: MB[f"bond/{_id_clean(mid)}_{u}_{v}"]

# Lightweight TTL emitter
class TTLWriter:
    def __init__(self, path):
        self.f = open(path, 'w', encoding='utf-8')
        # prefixes
        self.f.write(f"@prefix mt: <{MT}> .\n")
        self.f.write(f"@prefix owl: <{OWL}> .\n")
        self.f.write(f"@prefix rdf: <{RDF}> .\n\n")
    def write(self, s, p, o):
        self.f.write(f"{s.n3()} {p.n3()} {o.n3()} .\n")
    def close(self):
        self.f.close()

# Convert props (assume no Decimal left)
def format_literal(v):
    if isinstance(v, bool):
        return Literal(v, datatype=XSD.boolean)
    if isinstance(v, (int, float)):
        return Literal(v, datatype=XSD.double)
    return Literal(str(v), datatype=XSD.string)

# GNN export helpers
ELEMENTS = {"H":0,"C":1,"N":2,"O":3,"F":4,"Si":5,"P":6,"S":7,"Cl":8,"*":9}
NUM_ELEMENTS = len(ELEMENTS)
def atom_feat(attrs):
    idx = ELEMENTS.get(attrs.get("element","*"), ELEMENTS['*'])
    onehot = np.zeros(NUM_ELEMENTS, np.float32); onehot[idx]=1
    extras = np.array([attrs.get('atomic_num',0), 1.0 if attrs.get('is_dummy') else 0.0], np.float32)
    return np.hstack([onehot, extras])

def graph_to_data(mid, G, props):
    feats = [atom_feat(a) for _, a in G.nodes(data=True)]
    if not feats: return None
    x = torch.tensor(np.stack(feats), dtype=torch.float)
    ei, ea = [], []
    for u,v,a in G.edges(data=True):
        bt = {"SINGLE":1,"DOUBLE":2,"TRIPLE":3}.get(a.get('order',''),0)
        ei += [[u,v],[v,u]]; ea += [[bt],[bt]]
    if ei:
        edge_index = torch.tensor(ei, dtype=torch.long).t().contiguous()
        edge_attr  = torch.tensor(ea, dtype=torch.float)
    else:
        edge_index = torch.empty((2,0),dtype=torch.long)
        edge_attr  = torch.empty((0,1),dtype=torch.float)
    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    for k,v in props.items():
        if isinstance(v,(int,float,bool)): setattr(data, k, torch.tensor([v]))
    data.mid = mid; return data

# ─── Main export ─────────────────────────────────────────────
with open(MOL_PKL,'rb') as f: mols = pickle.load(f)
with open(CRYS_PKL,'rb') as f: crys = pickle.load(f)

# 1) TTL streaming
ttl = TTLWriter(OUT_TTL)
count=0
for mid, G, props in mols:
    m = URIRef(material_uri(mid))
    ttl.write(m, RDF.type, MT.Material)
    ttl.write(m, RDF.type, MT.Molecule)
    for k,v in props.items(): ttl.write(m, MT[k], format_literal(v))
    for n,attrs in G.nodes(data=True):
        a=URIRef(atom_uri(mid,n));
        cls = MT.AttachmentPoint if attrs.get('is_dummy') else MT.Atom
        ttl.write(a,RDF.type,cls)
        ttl.write(m,MT.hasAtom,a)
        ttl.write(a,MT.element,Literal(attrs.get('element',''),datatype=XSD.string))
        if attrs.get('is_dummy'): ttl.write(a,MT.isDummy,Literal(True,datatype=XSD.boolean))
    for u,v,attrs in G.edges(data=True):
        b=URIRef(bond_uri(mid,u,v))
        ttl.write(b,RDF.type,MT.Bond)
        ttl.write(b,MT.connects,URIRef(atom_uri(mid,u)))
        ttl.write(b,MT.connects,URIRef(atom_uri(mid,v)))
        ttl.write(b,MT.order,Literal(attrs.get('order',''),datatype=XSD.string))
        ttl.write(m,MT.hasBond,b)
    count+=1
    if count%1000==0: print(f"{count} molecules done")
count=0
for mid,G,props in crys:
    m = URIRef(material_uri(mid))
    ttl.write(m,RDF.type,MT.Material)
    ttl.write(m,RDF.type,MT.Crystal)
    for k,v in props.items(): ttl.write(m,MT[k],format_literal(v))
    for n,attrs in G.nodes(data=True):
        a=URIRef(atom_uri(mid,n))
        ttl.write(a,RDF.type,MT.Site)
        ttl.write(m,MT.hasSite,a)
        ttl.write(a,MT.element,Literal(attrs.get('element',''),datatype=XSD.string))
        coord=','.join(map(str,attrs.get('frac_coords',[])))
        ttl.write(a,MT.fracCoords,Literal(coord,datatype=XSD.string))
    for u,v,attrs in G.edges(data=True):
        b=URIRef(bond_uri(mid,u,v))
        ttl.write(b,RDF.type,MT.Contact)
        ttl.write(b,MT.connects,URIRef(atom_uri(mid,u)))
        ttl.write(b,MT.connects,URIRef(atom_uri(mid,v)))
        ttl.write(b,MT.distance,Literal(attrs.get('distance',0.0),datatype=XSD.float))
        ttl.write(m,MT.hasContact,b)
    count+=1
    if count%1000==0: print(f"{count} crystals done")
ttl.close()
print(f"✅ TTL written in {(time.time()-t0):.1f}s")

# 2) PyG export
print("⏳ Exporting PyG .pt…")
mol_data=[graph_to_data(mid,G,props) for mid,G,props in mols if graph_to_data(mid,G,props)]
torch.save(mol_data,str(OUT_PT_M))
crys_data=[graph_to_data(mid,G,props) for mid,G,props in crys if graph_to_data(mid,G,props)]
torch.save(crys_data,str(OUT_PT_C))
print(f"✅ PyG saved in {(time.time()-t0):.1f}s")
