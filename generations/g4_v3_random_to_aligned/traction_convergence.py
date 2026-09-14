"""Numerical checks for the weak v3-C deterministic-traction response."""

from __future__ import annotations

from dataclasses import replace
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from .model import G4V3Config, make_random_void_spec, run_fixed_material_traction


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "traction_convergence_summary.json"


def _worker(task) -> dict:
    label, cfg, spec = task
    result = run_fixed_material_traction(
        cfg,
        spec=spec,
        total_force=24.0,
        crosslink_probability=0.35,
    )
    frames = result["frames"]
    last_change = (
        frames[-1]["delta_radial_order_by_shell"]
        - frames[-2]["delta_radial_order_by_shell"]
    )
    return {
        "label": label,
        "dt": cfg.dt,
        "duration": cfg.traction_duration,
        "final_delta_radial_order": result["final_delta_radial_order_by_shell"].tolist(),
        "last_sample_delta_change": last_change.tolist(),
        "max_displacement": result["max_displacement"],
    }


def run() -> dict:
    base = G4V3Config()
    spec = make_random_void_spec(base)
    tasks = [
        ("half_dt", replace(base, dt=0.025), spec),
        ("double_duration", replace(base, traction_duration=960.0), spec),
    ]
    with Pool(processes=2) as pool:
        cases = pool.map(_worker, tasks)
    baseline = {
        "dt": 0.05,
        "duration": 480.0,
        "final_delta_radial_order": [0.0102137090, 0.0000174134, 0.0000197437],
        "max_displacement": 0.8497749128,
    }
    summary = {
        "protocol": "seed 41; same geometry; 24 nN; p_x=0.35; fixed material contacts",
        "baseline": baseline,
        "cases": cases,
        "near_cell_relative_dt_difference": abs(
            cases[0]["final_delta_radial_order"][0]
            - baseline["final_delta_radial_order"][0]
        ) / max(abs(baseline["final_delta_radial_order"][0]), 1e-12),
        "near_cell_change_when_duration_doubled": (
            cases[1]["final_delta_radial_order"][0]
            - baseline["final_delta_radial_order"][0]
        ),
    }
    OUTPUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
