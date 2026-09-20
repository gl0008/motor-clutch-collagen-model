"""Short diagnostic separating front-selection bias from long-time motion."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import json
import math
from pathlib import Path

import numpy as np

from .model import G4V4FConfig, make_matched_specs, run_v4f_condition


OUTPUT = Path(__file__).resolve().parent / "front_bias_diagnostic.json"


def _one(seed: int) -> dict:
    cfg = replace(
        G4V4FConfig(seed=seed),
        migration_sample_interval=30.0,
        telemetry_interval=30.0,
        event_duration=10.0,
    )
    matched = make_matched_specs(cfg)
    result = run_v4f_condition(
        cfg,
        "protrusion_guidance",
        seed,
        matched=matched,
        duration=240.0,
        store_positions=False,
    )
    event = result["front_events"][0] if result["front_events"] else None
    return {
        "seed": seed,
        "event": event,
        "final_direction_angle": result["metrics"]["final_direction_angle"],
        "net_displacement": result["metrics"]["net_displacement"],
    }


def main() -> None:
    with ProcessPoolExecutor(max_workers=10) as pool:
        rows = list(pool.map(_one, range(41, 61)))
    front_angles = np.asarray([
        math.atan2(row["event"]["direction"][1], row["event"]["direction"][0])
        for row in rows if row["event"] is not None
    ])
    summary = {
        "rows": rows,
        "front_resultant": float(abs(np.mean(np.exp(1j * front_angles)))),
        "front_mean_vector": np.mean(
            np.column_stack((np.cos(front_angles), np.sin(front_angles))), axis=0
        ).tolist(),
    }
    OUTPUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
