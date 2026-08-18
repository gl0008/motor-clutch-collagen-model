"""Export a saved G3 spheroid run to a compact JS payload for the web notebook.

Positions are stored as base64 int16 displacements from the initial geometry (same
scheme as Gloria's G2 demos). Per-frame protrusion/polarity data are small plain
arrays. Writes ``docs/g3_spheroid_demo_data.js`` defining ``window.G3_DATA``.
"""

from __future__ import annotations

import argparse
import base64
import json
import pickle
from pathlib import Path

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="run dir with frames.pkl + meta.pkl")
    ap.add_argument("--out", required=True, help="output .js path")
    ap.add_argument("--max-frames", type=int, default=80)
    ap.add_argument("--quantum", type=float, default=0.006)
    args = ap.parse_args()

    run = Path(args.run)
    frames = pickle.load(open(run / "frames.pkl", "rb"))
    meta = pickle.load(open(run / "meta.pkl", "rb"))

    if len(frames) > args.max_frames:
        keep = np.unique(np.linspace(0, len(frames) - 1, args.max_frames).round().astype(int))
        frames = [frames[i] for i in keep]

    r0 = np.asarray(meta["initial_positions"], dtype=float)
    edges = np.asarray(meta["edges"], dtype=int)
    l0 = np.linalg.norm(r0[edges[:, 1]] - r0[edges[:, 0]], axis=1)

    pos = np.stack([f["positions"] for f in frames], axis=0)          # F x M x 2
    delta = np.rint((pos - r0[None, :, :]) / args.quantum)
    if np.max(np.abs(delta)) > np.iinfo(np.int16).max:
        raise ValueError("quantisation overflow; raise --quantum")
    packed = np.asarray(delta, dtype="<i2").tobytes(order="C")

    S = meta["n_sectors"]

    def eng_pts(f):
        ep = np.asarray(f["engaged_points"], dtype=float)
        out = []
        for k in range(S):
            if np.isfinite(ep[k, 0]):
                out.append([round(float(ep[k, 0]), 2), round(float(ep[k, 1]), 2)])
            else:
                out.append(None)
        return out

    act_vmax = max(float(np.max([f["activity"].max() for f in frames])), 1e-6)

    data = {
        "meta": {
            "domain_size": meta["domain_size"],
            "cell_radius": meta["cell_radius"],
            "gap": meta.get("gap", 0.0),
            "n_sectors": S,
            "edges": edges.tolist(),
            "fixed": np.asarray(meta["fixed"], dtype=int).tolist(),
            "initial_positions": np.round(r0, 3).tolist(),
            "l0": np.round(l0, 4).tolist(),
            "crosslinks": [[int(a), round(float(aa), 4), int(b), round(float(ab), 4)]
                           for (a, aa, b, ab) in meta["crosslinks"]],
        },
        "encoding": {
            "quantum": args.quantum,
            "shape": list(delta.shape),
            "base64": base64.b64encode(packed).decode("ascii"),
        },
        "time": [round(float(f["time"]), 1) for f in frames],
        "center": [[round(float(f["center"][0]), 3), round(float(f["center"][1]), 3)] for f in frames],
        "activity": [np.round(f["activity"], 3).tolist() for f in frames],
        "length": [np.round(f["length"], 2).tolist() for f in frames],
        "engaged": [np.asarray(f["engaged"], dtype=int).tolist() for f in frames],
        "engaged_points": [eng_pts(f) for f in frames],
        "traction": [round(float(np.sum(f["traction"])), 2) for f in frames],
        "act_vmax": round(act_vmax, 4),
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    text = "window.G3_DATA=" + json.dumps(data, separators=(",", ":")) + ";\n"
    out.write_text(text, encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size/1e6:.2f} MB, {len(frames)} frames, "
          f"{len(r0)} beads, {len(edges)} bonds)")


if __name__ == "__main__":
    main()
