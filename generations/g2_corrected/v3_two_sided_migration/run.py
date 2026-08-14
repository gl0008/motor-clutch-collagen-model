"""Generate corrected V3 fixed/moving animation data and calibration summary."""

from pathlib import Path
import argparse
import json

import numpy as np

from model import (
    MigrationConfig,
    _random_stream,
    make_network_spec,
    run_fixed_moving_pair,
    run_migration,
    run_speed_ensemble,
)

HERE = Path(__file__).resolve().parent
import sys
sys.path.insert(0, str(HERE.parent))
from common.serialization import write_data_js  # noqa: E402


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reuse-fixed",
        action="store_true",
        help="reuse an existing fixed-cell encoded result and recompute moving only",
    )
    args = parser.parse_args()
    cfg = MigrationConfig()
    data_path = HERE / "demo" / "data.js"
    reused_fixed = False
    if args.reuse_fixed and data_path.exists():
        text = data_path.read_text(encoding="utf-8")
        old = json.loads(text[len("window.MODEL_DATA=") : -2])
        old["fixed"]["config"]["cell_drag"] = cfg.cell_drag
        spec = make_network_spec(cfg)
        stream = _random_stream(cfg, cfg.seed + 3001)
        pair = {
            "fixed": old["fixed"],
            "moving": run_migration(
                cfg, spec=spec, moving=True, random_stream=stream
            ),
            "config": cfg.__dict__,
            "comparison": "same ECM and same random stream; x-translation constraint only",
        }
        reused_fixed = True
    else:
        pair = run_fixed_moving_pair(cfg)
    speeds = run_speed_ensemble(cfg, trials=20)
    pair["calibration"] = {
        "cell_line": "MDA-MB-231",
        "trials": 20,
        "speeds_um_per_min": np.round(speeds, 5).tolist(),
        "median_um_per_min": float(np.median(speeds)),
        "target_um_per_min": [0.2, 0.4],
        "note": "lumped clutch ensemble; full visual trajectory includes ECM deformation",
    }
    write_data_js(
        data_path,
        pair,
        quantize_keys=("moving",) if reused_fixed else ("fixed", "moving"),
    )
    summary = {
        "median_speed_um_per_min": float(np.median(speeds)),
        "full_moving_path_um": float(
            np.sum(np.linalg.norm(np.diff(pair["moving"]["cell_center"], axis=0), axis=1))
        ),
        "full_net_displacement_um": float(
            np.linalg.norm(pair["moving"]["cell_center"][-1])
        ),
        "fixed_net_displacement_um": 0.0,
    }
    (HERE / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
