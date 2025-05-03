#!/usr/bin/env python3
# scripts/neo4j_load_kg.py
"""

REQUIRES Neo4j Neosemantics plugin

"""

import os, sys
from pathlib import Path
from neo4j import GraphDatabase

# credentials
try:
    from scripts.getenv import NEO4J_URI, NEO4J_USER, NEO4J_PW
except ImportError:
    NEO4J_URI  = os.getenv('NEO4J_URI','bolt://localhost:7687')
    NEO4J_USER = os.getenv('NEO4J_USER','neo4j')
    NEO4J_PW   = os.getenv('NEO4J_PW','neo4j')

ROOT    = Path(__file__).resolve().parent.parent
TTL     = ROOT / "exports" / "KnowledgeGraphs" / "materials.ttl"

driver  = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PW))

with driver.session() as session:
    # initialize n10s
    session.run("""
      CALL n10s.graphconfig.init({
        handleRDF:'OVERRIDE',
        handleVocabUris:'SHORTEN',
        keepLangTag:true,
        baseUri:'http://materialsdb.org/data/',
        namespace:{mt:'http://materialsdb.org/schema#'}
      })
    """)
    # import TTL
    ttl_data = TTL.read_text()
    session.run("""
      CALL n10s.rdf.import.inline($ttl, 'Turtle')
    """, ttl=ttl_data)

driver.close()
print("✅ Knowledge graph imported into Neo4j.")
