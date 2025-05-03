#!/usr/bin/env python3
# scripts/neo4j_load_pg.py

"""
Load material property‑graphs into Neo4j in bulk.

Reads:
  data/graphs/mol_graphs.pkl      # list of (mid, G) or (mid, G, props)
  data/graphs/crystal_graphs.pkl  # same

For each material it creates:
  • (m:Material:Molecule {…props})
  • UNWIND nodes → (a:Atom {mat_id, idx, …})
     + MERGE (m)-[:HAS_ATOM]->(a)
  • UNWIND edges → MATCH (a1:Atom …),(a2:Atom …) MERGE (a1)-[:BOND {…}]->(a2)

And similarly for `Crystal` / `Site`.

Usage:
    python scripts/neo4j_load_pg.py
"""

import os, pickle
from pathlib import Path
from neo4j import GraphDatabase

# ─── Credentials ─────────────────────────────────────────────────────────
try:
    from scripts.getenv import NEO4J_URI, NEO4J_USER, NEO4J_PW
except ImportError:
    NEO4J_URI  = os.getenv('NEO4J_URI', 'bolt://localhost:7687')
    NEO4J_USER = os.getenv('NEO4J_USER', 'neo4j')
    NEO4J_PW   = os.getenv('NEO4J_PW', 'neo4j')

# ─── Paths ────────────────────────────────────────────────────────────────
ROOT     = Path(__file__).resolve().parent.parent
GRAPHS   = ROOT / "data" / "graphs"
MOL_PKL  = GRAPHS / "mol_graphs.pkl"
CRYS_PKL = GRAPHS / "crystal_graphs.pkl"

# ─── Driver ───────────────────────────────────────────────────────────────
driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PW))


def ingest_material(tx, mid, kind, props):
    """
    Create the Material node with label = kind (Molecule or Crystal).
    """
    tx.run(
        f"""
        MERGE (m:Material:{kind} {{id:$mid}})
        SET m += $props
        """,
        mid=mid, props=props or {}
    )


def ingest_nodes(tx, mid, nodes, label, rel):
    """
    Bulk ingest all nodes (Atom or Site) and attach them to the material.
    nodes = [ {idx: int, attrs: dict}, ... ]
    """
    tx.run(
        f"""
        MATCH (m:Material {{id:$mid}})
        UNWIND $nodes AS n
          MERGE (a:{label} {{ mat_id:$mid, idx:n.idx }})
          SET a += n.attrs
          MERGE (m)-[:{rel}]->(a)
        """,
        mid=mid, nodes=nodes
    )


def ingest_edges(tx, mid, edges):
    """
    Bulk ingest all bonds (between Atom or Site) for this material.
    edges = [ {u:int, v:int, props:dict}, ... ]
    """
    tx.run(
        """
        UNWIND $edges AS e
          MATCH (a1 {mat_id:$mid, idx:e.u}), (a2 {mat_id:$mid, idx:e.v})
          MERGE (a1)-[r:BOND]->(a2)
          SET r += e.props
        """,
        mid=mid, edges=edges
    )


with driver.session() as session:
    # --- Molecules ---
    mols = pickle.load(open(MOL_PKL, 'rb'))
    for item in mols:
        if len(item) == 2:
            mid, G = item
            props   = {}
        else:
            mid, G, props = item

        # 1) Material node
        session.execute_write(ingest_material, mid, "Molecule", props)

        # 2) Bulk Atom nodes
        atom_list = [
            {"idx": n, "attrs": attrs}
            for n, attrs in G.nodes(data=True)
        ]
        session.execute_write(ingest_nodes, mid, atom_list, "Atom", "HAS_ATOM")

        # 3) Bulk BOND relationships
        edge_list = [
            {"u": u, "v": v, "props": attrs}
            for u, v, attrs in G.edges(data=True)
        ]
        session.execute_write(ingest_edges, mid, edge_list)

    # --- Crystals ---
    crys = pickle.load(open(CRYS_PKL, 'rb'))
    for item in crys:
        if len(item) == 2:
            mid, G = item
            props   = {}
        else:
            mid, G, props = item

        session.execute_write(ingest_material, mid, "Crystal", props)

        site_list = [
            {"idx": n, "attrs": attrs}
            for n, attrs in G.nodes(data=True)
        ]
        session.execute_write(ingest_nodes, mid, site_list, "Site", "HAS_SITE")

        contact_list = [
            {"u": u, "v": v, "props": attrs}
            for u, v, attrs in G.edges(data=True)
        ]
        session.execute_write(ingest_edges, mid, contact_list)

    print(f"✅ Ingested {len(mols)} Molecules and {len(crys)} Crystals into Neo4j.")

# swallow that destructor bug gracefully
try:
    driver.close()
except Exception:
    pass
