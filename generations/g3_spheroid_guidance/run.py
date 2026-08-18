"""Run a G3 spheroid-guidance simulation and save frames + traces to disk.

Usage:
    python run.py --out <dir> [--fibers N] [--duration S] [--dt DT] [--seed K] [--fixed]
"""

from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np

from model import SpheroidConfig, run_spheroid


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--fibers", type=int, default=120)
    ap.add_argument("--domain", type=float, default=220.0)
    ap.add_argument("--radius", type=float, default=20.0)
    ap.add_argument("--gap", type=float, default=9.0)
    ap.add_argument("--duration", type=float, default=5400.0)
    ap.add_argument("--dt", type=float, default=0.1)
    ap.add_argument("--sample", type=float, default=60.0)
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--aligned", action="store_true")
    ap.add_argument("--director", type=float, default=0.0)
    ap.add_argument("--fixed", action="store_true", help="clamp the cell (no migration)")
    ap.add_argument("--attempts", type=int, default=40)
    args = ap.parse_args()

    cfg = SpheroidConfig(
        n_fibers=args.fibers, domain_size=args.domain, cell_radius=args.radius,
        gap=args.gap, duration=args.duration, dt=args.dt, sample_interval=args.sample,
        seed=args.seed, aligned=args.aligned, director_angle=args.director,
        generation_attempts=args.attempts,
    )
    res = run_spheroid(cfg, seed=args.seed, moving=not args.fixed)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    meta = {
        "edges": res["edges"],
        "edge_fiber": res["edge_fiber"],
        "fixed": res["fixed"],
        "initial_positions": res["initial_positions"],
        "crosslinks": [(c.edge_a, c.alpha_a, c.edge_b, c.alpha_b) for c in res["crosslinks"]],
        "cell_radius": cfg.cell_radius,
        "domain_size": cfg.domain_size,
        "gap": cfg.gap,
        "n_sectors": cfg.n_sectors,
    }
    with open(out / "frames.pkl", "wb") as fh:
        pickle.dump(res["frames"], fh)
    with open(out / "meta.pkl", "wb") as fh:
        pickle.dump(meta, fh)
    np.savez(out / "traces.npz", **res["traces"])
    tr = res["traces"]
    summary = {
        "fibers": len(res["network"].fibers),
        "beads": int(len(res["initial_positions"])),
        "crosslinks": len(res["crosslinks"]),
        "connected_fraction": res["report"]["connected_fraction"],
        "duration_s": cfg.duration,
        "frames": len(res["frames"]),
        "final_engaged": int(tr["engaged_count"][-1]),
        "peak_traction_nN": float(np.max(tr["total_traction"])),
        "path_length_um": float(tr["path_length"][-1]),
        "net_displacement_um": float(np.hypot(tr["cell_x"][-1], tr["cell_y"][-1])),
        "speed_um_per_min": float(tr["path_length"][-1] / (cfg.duration / 60.0)),
        "order_near_start": float(tr["order_near"][0]),
        "order_near_end": float(tr["order_near"][-1]),
        "max_bead_disp_um": float(np.max(tr["max_disp"])),
        "polarity_angle_end": float(tr["polarity_angle"][-1]),
    }
    with open(out / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
