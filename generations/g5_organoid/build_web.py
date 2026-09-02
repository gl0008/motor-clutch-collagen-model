"""Export a compact G5 animation for the GitHub Pages site (docs/g5-data.js).

Runs two short G5 simulations on a web-sized network and writes their frames as a
single ``window.G5_DATA = {...}`` script (same lazy-load-friendly pattern as
g4-v2-manifest.js; no fetch, so it also renders from a local file).  Positions are
stored as int16 offsets from the initial geometry to keep the payload small.

  contract : Stage B — fixed organoid, collagen turns radial (radial-order readout)
  invade   : Stage D — released cells invade outward (cell trails)

    python generations/g5_organoid/build_web.py
"""

from __future__ import annotations

import base64
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.model import (  # noqa: E402
    OrganoidConfig, run_organoid_pull, run_organoid_invasion, radial_alignment_profile)

SCALE = 0.01  # µm per int16 step (±327 µm range, 0.01 µm resolution)


def _enc(frames: np.ndarray, initial: np.ndarray) -> dict:
    """frames (F,N,2), initial (N,2) -> base64 int16 offsets."""
    q = np.round((frames - initial[None]) / SCALE).astype(np.int16)
    return {"b64": base64.b64encode(np.ascontiguousarray(q).tobytes()).decode(),
            "F": int(q.shape[0]), "N": int(q.shape[1])}


# Web-sized network: small enough to animate smoothly, with the A+B connected corona.
WEB = OrganoidConfig(
    n_fibers=120, n_corona_fibers=48, corona_band=10.0, gap=1.0,
    organoid_radius=40.0, cell_spacing=18.0, domain_size=300.0, boundary_width=6.0,
    bead_spacing=1.0, total_pull_force=20.0, generation_attempts=15,
)


def main(out: str = "docs/g5-data.js") -> None:
    # --- Stage B: contract -> radial ---
    cB = OrganoidConfig(**{**WEB.__dict__, "duration": 220.0, "sample_interval": 8.0})
    b = run_organoid_pull(cB, seed=23, snapshots=True)
    initial = b["initial_positions"]
    edges = b["edges"]
    snapsB = b["snapshots"]
    orderB = [fr["global_radial_order"] for fr in b["frames"]]
    timesB = [fr["time"] for fr in b["frames"]]
    centers = b["centers"]

    # --- Stage D: released cells invade ---
    cD = OrganoidConfig(**{**WEB.__dict__, "duration": 420.0, "sample_interval": 14.0,
                           "cc_adhesion": 2.0, "max_cell_speed": 0.03,
                           "contact_update_interval": 2.0})
    d = run_organoid_invasion(cD, seed=23, snapshots=True)
    snapsD = d["bead_snapshots"]
    cellsD = d["cell_snapshots"]
    initialD = d["initial_positions"]
    edgesD = d["edges"]
    invaded = [fr["mean_cell_radial_disp"] for fr in d["frames"]]
    timesD = [fr["time"] for fr in d["frames"]]

    data = {
        "domain": WEB.domain_size,
        "cellRadius": WEB.cell_radius,
        "scale": SCALE,
        "contract": {
            "initial": [round(float(v), 3) for v in initial.reshape(-1)],
            "edges": [int(v) for v in edges.reshape(-1)],
            "cells": [round(float(v), 2) for v in centers.reshape(-1)],
            "times": [round(float(t), 1) for t in timesB],
            "order": [round(float(o), 4) for o in orderB],
            "pos": _enc(snapsB, initial),
            "beads": int(len(initial)),
        },
        "invade": {
            "initial": [round(float(v), 3) for v in initialD.reshape(-1)],
            "edges": [int(v) for v in edgesD.reshape(-1)],
            "times": [round(float(t), 1) for t in timesD],
            "invaded": [round(float(v), 2) for v in invaded],
            "pos": _enc(snapsD, initialD),
            "cells": _enc(cellsD, cellsD[0]),
            "cells0": [round(float(v), 2) for v in cellsD[0].reshape(-1)],
            "beads": int(len(initialD)),
        },
    }
    text = "window.G5_DATA=" + json.dumps(data, separators=(",", ":")) + ";\n"
    Path(REPO / out).write_text(text, encoding="utf-8")
    kb = len(text) / 1024
    print("wrote %s (%.0f KB) | contract %d beads %d frames | invade %d beads %d frames" % (
        out, kb, data["contract"]["beads"], data["contract"]["pos"]["F"],
        data["invade"]["beads"], data["invade"]["pos"]["F"]))


if __name__ == "__main__":
    main()
