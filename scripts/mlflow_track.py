#!/usr/bin/env python3
# scripts/mlflow_track.py

import time
import pickle
import mlflow
from pathlib import Path

ROOT     = Path(__file__).resolve().parent.parent
GRAPHS   = ROOT / "data" / "graphs"
MOL_PKL  = GRAPHS / "mol_graphs.pkl"
CRYS_PKL = GRAPHS / "crystal_graphs.pkl"

def main():
    mlflow.set_experiment("materials-graph-pipeline")
    with mlflow.start_run():
        t0 = time.time()
        mols = pickle.load(open(MOL_PKL, "rb"))
        crys = pickle.load(open(CRYS_PKL, "rb"))
        dt = time.time() - t0

        # log metrics
        mlflow.log_metric("build_time_s", dt)
        mlflow.log_metric("n_molecules", len(mols))
        mlflow.log_metric("n_crystals", len(crys))

        # log artifacts
        mlflow.log_artifact(str(MOL_PKL), artifact_path="graphs")
        mlflow.log_artifact(str(CRYS_PKL), artifact_path="graphs")

        print(f"Logged run in {dt:.1f}s with {len(mols)} mols, {len(crys)} crys.")

if __name__ == "__main__":
    main()
