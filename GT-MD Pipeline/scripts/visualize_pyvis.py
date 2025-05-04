#!/usr/bin/env python3
# scripts/visualize_pyvis.py

"""
Pick and render molecule/crystal graphs from your pickles as interactive HTML via PyVis.

Usage:
  # visualize all molecules
  python scripts/visualize_pyvis.py --type molecules

  # visualize the 3rd crystal (0‑based index)
  python scripts/visualize_pyvis.py --type crystals --index 2

  # visualize both sets
  python scripts/visualize_pyvis.py --type all
"""

import argparse
import pickle
from pathlib import Path
from pyvis.network import Network

ROOT = Path(__file__).resolve().parent.parent

def load_graph_list(pkl_path):
    if not pkl_path.exists():
        return []
    return pickle.load(open(pkl_path, "rb"))

def normalize(graph_items):
    """
    Each item may be (mid, G) or (mid, G, props); we only need mid & G.
    """
    return [(item[0], item[1]) for item in graph_items]

def visualize_one(mid, G, prefix):
    out_dir = ROOT / "html_maps"
    out_dir.mkdir(exist_ok=True)
    net = Network(notebook=False)
    net.force_atlas_2based()

    for n, attrs in G.nodes(data=True):
        label = attrs.get("element", str(n))
        color = "#FFCCCC" if attrs.get("is_dummy") else "#CCCCFF"
        net.add_node(n, label=label, title=str(attrs), color=color)

    for u, v, attrs in G.edges(data=True):
        net.add_edge(u, v, title=str(attrs))

    out_file = out_dir / f"{prefix}_{mid}.html"
    # pass a string, not a Path
    net.write_html(str(out_file), open_browser=False, notebook=False)
    print(f"🔍 Written {out_file}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--type", choices=["molecules","crystals","all"], default="molecules")
    p.add_argument("--index", type=int, default=None)
    args = p.parse_args()

    mol_pkl  = ROOT / "data" / "graphs" / "mol_graphs.pkl"
    crys_pkl = ROOT / "data" / "graphs" / "crystal_graphs.pkl"

    mols = normalize(load_graph_list(mol_pkl))
    crys = normalize(load_graph_list(crys_pkl))

    if args.type in ("molecules","all"):
        if not mols:
            print("⚠️  No molecule graphs found.")
        else:
            sel = mols if args.index is None else [mols[args.index]] if 0 <= args.index < len(mols) else []
            for mid, G in sel:
                visualize_one(mid, G, prefix="molecule")

    if args.type in ("crystals","all"):
        if not crys:
            print("⚠️  No crystal graphs found.")
        else:
            sel = crys if args.index is None else [crys[args.index]] if 0 <= args.index < len(crys) else []
            for mid, G in sel:
                visualize_one(mid, G, prefix="crystal")

if __name__ == "__main__":
    main()
