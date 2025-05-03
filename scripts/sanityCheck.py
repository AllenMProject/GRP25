#!/usr/bin/env python3
# scripts/sanity_check.py

"""
Sanity‑check your pickled graphs before downstream steps.

Checks both:
  data/graphs/mol_graphs.pkl
  data/graphs/crystal_graphs.pkl

and for each nonempty file, prints min/avg/max of node‑counts and edge‑counts.
"""

import pickle
from pathlib import Path

def stats(graph_items):
    """
    Given a list of tuples where the second element is a NetworkX graph,
    return (nodemin, nodeavg, nodemax, edgemin, edgeavg, edgemax).
    """
    # extract the graphs (always second element)
    graphs = [item[1] for item in graph_items]
    node_counts = [g.number_of_nodes() for g in graphs]
    edge_counts = [g.number_of_edges()  for g in graphs]
    return (
        min(node_counts),
        sum(node_counts) // len(node_counts),
        max(node_counts),
        min(edge_counts),
        sum(edge_counts) // len(edge_counts),
        max(edge_counts),
    )

def print_stats(graph_items, label):
    if not graph_items:
        print(f"{label}:  ✗ no entries, skipping.")
    else:
        nmin, navg, nmax, emin, eavg, emax = stats(graph_items)
        print(f"{label}:")
        print(f"  nodes: min={nmin}, avg={navg}, max={nmax}")
        print(f"  edges: min={emin}, avg={eavg}, max={emax}")

def main():
    root      = Path(__file__).resolve().parent.parent
    graphs_dir = root / "data" / "graphs"
    targets = [
        ("Molecules",   graphs_dir / "mol_graphs.pkl"),
        ("Crystals",    graphs_dir / "crystal_graphs.pkl"),
    ]

    for label, pkl_path in targets:
        if not pkl_path.exists():
            print(f"{label}:  ✗ file not found at {pkl_path}, skipping.")
            continue

        try:
            items = pickle.load(open(pkl_path, "rb"))
        except Exception as e:
            print(f"{label}:  ✗ error loading {pkl_path}: {e}")
            continue

        print_stats(items, label)

if __name__ == "__main__":
    main()
