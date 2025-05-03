#!/usr/bin/env python3
# scripts/kg_service.py

"""
1) --export-kg : regenerate TTL from pickles
2) --serve     : serve SPARQL + GraphQL on :5000

Serve Instructions:
1. Open Powershell Prompt
2. Go to the project anaconda environment and directory
3. Run python scripts/kg_service.py --export-kg IF AND ONLY IF you don't already have a .ttl file available
4. Run python scripts/kg_service.py --serve # This will start the server
5. Open browser and go to http://127.0.0.1:5000/graphql

GraphQL Example Query:

query {
  materials(limit:3) {
    uri
    formula
    source
    atoms {
      uri
      element
      isDummy
    }
    bonds {
      uri
      order
    }
  }
}



SPARQL Example Query:
6. Open browser and go to http://127.0.0.1:5000/sparql

SPARQL Example Query:
http://127.0.0.1:5000/sparql?query=PREFIX%20mt%3A%20%3Chttp%3A//materialsdb.org/schema%23%3E%20SELECT%20%3Fm%20%3Ff%20%3Fs%20WHERE%20%7B%20%3Fm%20a%20mt%3AMaterial%20%3B%20mt%3Aformula%20%3Ff%20%3B%20mt%3Asource%20%3Fs%20.%20%7D%20LIMIT%203

curling Example:

curl -G http://127.0.0.1:5000/sparql \
     --data-urlencode "query=PREFIX mt: <http://materialsdb.org/schema#> \
SELECT ?m ?f ?s WHERE { \
  ?m a mt:Material ; \
     mt:formula ?f ; \
     mt:source  ?s . \
} LIMIT 3"

should return something akin to
[
  {
    "m": "http://materialsdb.org/data/material/mp-1009127",
    "f": "MgO",
    "s": "materialsproject"
  },
  {
    "m": "http://materialsdb.org/data/material/mp-1009129",
    "f": "MgO",
    "s": "materialsproject"
  },
  {
    "m": "http://materialsdb.org/data/material/mp-1065577",
    "f": "MgO",
    "s": "materialsproject"
  }
]
"""

import sys
import argparse
from pathlib import Path

from flask import Flask, request, jsonify, redirect
from rdflib import Graph, Namespace, Literal
from rdflib.namespace import RDF, XSD
from flask_graphql import GraphQLView
import graphene

# ─── Configuration ────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
KG_FILE = ROOT / "exports" / "KnowledgeGraphs" / "materials.ttl"

MT = Namespace("http://materialsdb.org/schema#")
MB = Namespace("http://materialsdb.org/data/")


# ─── Helpers ──────────────────────────────────────────────────────────
def load_exporter():
    """
    Dynamically load scripts/export_knowledgeGraphs.py as a module
    and return its export_kg() function.
    """
    from importlib.machinery import SourceFileLoader
    loader = SourceFileLoader(
        "export_kg_module",
        str(SCRIPTS / "export_knowledgeGraphs.py")
    )
    mod = loader.load_module()
    return mod.export_kg


def build_sparql_graph():
    """Load the TTL file into an rdflib.Graph for querying."""
    g = Graph()
    try:
        g.parse(KG_FILE.as_uri(), format="turtle")
    except Exception:
        g.parse(data=KG_FILE.read_bytes(), format="turtle")
    return g


# ─── CLI parsing ──────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--export-kg", action="store_true",
                    help="(Re)export the TTL knowledge graph from pickles")
parser.add_argument("--serve", action="store_true",
                    help="Launch SPARQL+GraphQL server on :5000")
args = parser.parse_args()

# ─── Handle export first ─────────────────────────────────────────────
if args.export_kg:
    export_kg = load_exporter()
    print("⏳ Exporting knowledge graph…")
    export_kg()
    print("✅ TTL regenerated at", KG_FILE)
    if not args.serve:
        sys.exit(0)

# ─── If serving, verify TTL and spin up endpoints ────────────────────
if args.serve:
    if not KG_FILE.exists():
        print(f"‼️  Knowledge graph not found at {KG_FILE!r}.")
        print("    First run with --export-kg to generate it, then rerun --serve.")
        sys.exit(1)

    SPARQL_KG = build_sparql_graph()
    app = Flask(__name__)


    # → Redirect root to GraphiQL UI
    @app.route("/")
    def index():
        return redirect("/graphql")


    # → SPARQL endpoint
    @app.route("/sparql", methods=["GET", "POST"])
    def sparql():
        q = request.values.get("query")
        if not q:
            return "Provide SPARQL as `?query=` or POST form field.", 400
        res = SPARQL_KG.query(q)
        out = [
            {str(var): row[i] for i, var in enumerate(res.vars)}
            for row in res
        ]
        return jsonify(out)


    # → GraphQL schema
    class Atom(graphene.ObjectType):
        uri = graphene.String()
        element = graphene.String()
        isDummy = graphene.Boolean()


    class Bond(graphene.ObjectType):
        uri = graphene.String()
        order = graphene.String()


    class Material(graphene.ObjectType):
        uri = graphene.String()
        formula = graphene.String()
        source = graphene.String()
        atoms = graphene.List(Atom)
        bonds = graphene.List(Bond)


    class Query(graphene.ObjectType):
        materials = graphene.List(Material, limit=graphene.Int())

        def resolve_materials(self, info, limit=20):
            prefix = f"PREFIX mt: <{MT}>"
            q = (
                f"{prefix}\n"
                f"SELECT ?m WHERE {{\n"
                f"  ?m a mt:Material .\n"
                f"}} LIMIT {limit}"
            )
            mids = [str(r["m"]) for r in SPARQL_KG.query(q)]
            results = []
            for m in mids:
                # material-level props
                f_q = (
                    f"{prefix}\n"
                    f"SELECT ?f ?s WHERE {{\n"
                    f"  <{m}> mt:formula ?f ;\n"
                    f"         mt:source  ?s .\n"
                    f"}}"
                )
                formula = ""
                source = ""
                for r in SPARQL_KG.query(f_q):
                    formula = str(r[0])
                    source = str(r[1])
                    break

                # atoms (with optional isDummy)
                a_q = (
                    f"{prefix}\n"
                    f"SELECT ?a ?el ?d WHERE {{\n"
                    f"  <{m}> mt:hasAtom ?a .\n"
                    f"  ?a mt:element ?el .\n"
                    f"  OPTIONAL {{ ?a mt:isDummy ?d . }}\n"
                    f"}}"
                )
                atoms = [
                    Atom(
                        uri=str(r["a"]),
                        element=str(r["el"]),
                        isDummy=bool(r.get("d", False))
                    )
                    for r in SPARQL_KG.query(a_q)
                ]

                # bonds
                b_q = (
                    f"{prefix}\n"
                    f"SELECT ?b ?ord WHERE {{\n"
                    f"  <{m}> mt:hasBond ?b .\n"
                    f"  ?b mt:order ?ord .\n"
                    f"}}"
                )
                bonds = [
                    Bond(uri=str(r["b"]), order=str(r["ord"]))
                    for r in SPARQL_KG.query(b_q)
                ]

                results.append(Material(
                    uri=m,
                    formula=formula,
                    source=source,
                    atoms=atoms,
                    bonds=bonds
                ))

            return results


    schema = graphene.Schema(query=Query)
    app.add_url_rule(
        "/graphql",
        view_func=GraphQLView.as_view("graphql", schema=schema, graphiql=True)
    )


# → Start server
    print("🚀 Serving SPARQL & GraphQL at http://0.0.0.0:5000")
    app.run(host="0.0.0.0", port=5000)

else:
    parser.print_help()

