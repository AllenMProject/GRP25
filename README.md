Open‑Access Graph‑Theoretic Materials Database (GTDB)
===========================================

This project ingests molecular & crystalline materials data,
converts them into both property‐graphs (NetworkX) and
knowledge‐graphs (RDF/OWL), exports multiple formats
(JSON, CSV, edgelists, .ttl), builds GNN‐ready .pt files
(PyTorch‐Geometric), serves SPARQL + GraphQL endpoints,
and provides sanity checks & interactive visualizations.

Overview
--------
This repository contains an end‑to‑end pipeline for:
  • Ingesting heterogeneous materials data (PubChem, MaterialsProject, etc.)
  • Converting each material into a graph representation (molecules and crystals)
  • Exporting property graphs (NetworkX pickles) and a knowledge graph (RDF/Turtle)
  • Serving the knowledge graph via SPARQL and GraphQL APIs
  • Loading both property and knowledge graphs into Neo4j
  • Sanity checking and visualizing graphs via PyVis
  • Tracking pipeline runs with MLflow

Key components:
  • build_graphs.py       → Stage 1: JSON → NetworkX graphs
  • export_propertyGraphs → Stage 2: export JSON/CSV/edgelist + .pt
  • export_knowledgeGraphs→ Stage 3: graphs → Turtle + .pt
  • kg_service.py         → SPARQL & GraphQL server
  • sanity_check.py       → Quick stats on graph sizes
  • visualize_pyvis.py    → PyVis HTML visualizations
  • neo4j_load_pg.py      → Load property‐graphs into Neo4j
  • neo4j_load_kg.py      → Load knowledge‐graphs into Neo4j
  • mlflow tracking       → track runs & artifacts
  • JupyterLab (optional) → interactive exploration

Repository Layout
-----------------
README.txt                   # this overview
Requirements.txt             # Python dependencies
Getting_Started.txt          # initial setup instructions
Pipeline_Explanation.txt     # full pipeline flow
How_to_Check.txt             # verification & sanity checking
Getting_Results.txt          # how to view outputs

/data
  /cleaned/        	     # containing standardize data
  /graphs/         	     # containing all built graphs
/exports		     # containing exported graphs (.csv, .txt, .json, .pt)
/scripts
  build_graphs.py             # produce property graphs
  export_knowledgeGraphs.py   # produce TTL knowledge graph
  kg_service.py               # SPARQL/GraphQL server
  neo4j_load_pg.py            # load property graphs into Neo4j
  neo4j_load_kg.py            # load knowledge graph into Neo4j via n10s
  sanity_check.py             # statistics sanity check
  visualize_pyvis.py          # interactive HTML graph visualizations
  ...                         # helper scripts





