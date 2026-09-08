"""Run the prescribed-contraction benchmark across independent geometries."""

from __future__ import annotations

import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from .model import G4V3Config, run_prescribed_contraction


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "seed_check_summary.json"


def _worker(seed: int) -> dict:
    result = run_prescribed_contraction(G4V3Config(seed=seed), contraction=0.20)
    initial = result["initial_radial_order_by_shell"]
    final = result["final_radial_order_by_shell"]
    return {
        "seed": seed,
        "initial_radial_order": initial.tolist(),
        "delta_radial_order": (final - initial).tolist(),
        "stretch_to_bend_ratio": result["final_stretch_to_bend_ratio"],
        "required_inward_force": result["frames"][-1]["required_inward_boundary_force"],
        "natural_contact_fragments": result["report"]["contact_fragments"],
    }


def run(seeds=(41, 42, 43, 44, 45)) -> dict:
    with Pool(processes=min(len(seeds), 5)) as pool:
        cases = pool.map(_worker, seeds)
    delta = np.asarray([case["delta_radial_order"] for case in cases])
    summary = {
        "protocol": "20% prescribed contraction; full permanent intersection linking",
        "seeds": list(seeds),
        "cases": cases,
        "delta_radial_order_mean": np.mean(delta, axis=0).tolist(),
        "delta_radial_order_median": np.median(delta, axis=0).tolist(),
        "near_cell_positive_fraction": float(np.mean(delta[:, 0] > 0.0)),
        "near_cell_above_0_05_fraction": float(np.mean(delta[:, 0] > 0.05)),
    }
    OUTPUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
