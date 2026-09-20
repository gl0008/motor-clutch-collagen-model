"""Export the G5 organoid gv_*_full.npz runs to a Gloria-style interactive web viewer dataset.

Mirrors gl0008's g4-v4e packing: bead positions are stored as base64-packed int16 DELTAS from the
frame-0 base geometry (decoded in JS as base + scale*value), so the browser redraws every frame on a
<canvas> exactly like her viewer.  One JS file per condition (loaded on demand) + a shared manifest.

    python .../export_web_viewer.py            # export all gv_*_full.npz found in output/
Writes: docs/g5-data/<tag>.js  and  docs/g5-organoid-manifest.js
Personal testing visuals (CLAUDE.md 7.5).  No swirling claim implied.
"""

from __future__ import annotations

import base64
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

DOCS = REPO / "docs"
DATA = DOCS / "g5-data"
OUT = REPO / "output"

# condition tag -> label.  NARRATIVE = the honest PHASE story we found: adhesion robustly sets
# collective-vs-escape, but the leader-led strand + its R4 guidance rescue are CROSSLINK-DENSITY
# GATED (work only in dense ECM where leaders detach; fail in sparse G4D-parity ECM where the
# organoid holds together on its own).  Redundant R0 and R2 no-leader are left out.
CASES = [
    ("r1_cohesive", "① Cohesive organoid — stays collective (baseline)"),
    ("r1_escape", "② Full partial-EMT → single-cell escape (ROBUST)"),
    ("r2_leader_front", "③ Strong leaders, cohesive matrix → no strand (leaders alone aren't enough)"),
    ("r4_dense_unguided", "④ DENSE ECM + leaders, NO guidance → leaders DETACH from the bulk"),
    ("r4_dense_guided", "⑤ DENSE ECM + follower guidance → followers FOLLOW (R4 rescue works)"),
    ("r4_unguided", "⑥ SPARSE ECM (G4D-parity) + leaders → organoid HOLDS (nothing to rescue)"),
    ("r4_guided", "⑦ SPARSE ECM + guidance → over-driven → SCATTERS (R4 fails here)"),
]

FRAME_STRIDE = 2   # decimate 121 -> ~61 frames to keep the JS payload light (still smooth at 10 fps)


def pack_i16(states: np.ndarray, base: np.ndarray) -> dict:
    """base64-packed int16 deltas from base (Gloria's format)."""
    delta = states - base[None, :, :]
    m = float(np.max(np.abs(delta))) if delta.size else 0.0
    scale = max(1e-6, m / 32700.0)
    vals = np.rint(delta / scale).astype("<i2")
    return {"scale": scale, "shape": list(vals.shape),
            "base64": base64.b64encode(vals.tobytes()).decode("ascii")}


def net_traction(sf, snrm, n_sec, n_cells):
    """Per-frame per-cell net traction vector (outward = invasion dir): -sum_s f_s * inward_s."""
    F = sf.shape[0]
    out = np.zeros((F, n_cells, 2), dtype=np.float32)
    for f in range(F):
        live = np.isfinite(snrm[f, :, 0]) & (sf[f] > 1e-9)
        for s in np.flatnonzero(live):
            out[f, s // n_sec] -= sf[f, s] * snrm[f, s]
    return out


def export_case(tag: str, label: str) -> dict | None:
    path = OUT / f"gv_{tag}_full.npz"
    if not path.exists():
        print("skip (missing)", path.name); return None
    z = np.load(path, allow_pickle=True)
    beads = np.asarray(z["bead_snapshots"], dtype=np.float64)[::FRAME_STRIDE]
    cells = np.asarray(z["cell_snapshots"], dtype=np.float64)[::FRAME_STRIDE]
    times = np.asarray(z["times"], dtype=np.float64)[::FRAME_STRIDE]
    edges = np.asarray(z["edges"], dtype=np.int32)
    sf = np.asarray(z["site_force_snapshots"], dtype=np.float64)[::FRAME_STRIDE]
    snrm = np.asarray(z["site_normal_snapshots"], dtype=np.float64)[::FRAME_STRIDE]
    n_sec = int(z["n_contact_sectors"])
    xl_edge = np.asarray(z["crosslink_edge"], dtype=np.int32)
    xl_alpha = np.asarray(z["crosslink_alpha"], dtype=np.float64)
    red_ids = [int(x) for x in z["red_ids"]] if "red_ids" in z.files else []
    red_pf = ([[int(c) for c in r] for r in z["red_per_frame"]][::FRAME_STRIDE]
              if "red_per_frame" in z.files else None)
    base = beads[0]
    net = net_traction(sf, snrm, n_sec, cells.shape[1])
    case = {
        "tag": tag, "label": label,
        "base": np.round(base, 3).tolist(),
        "positions": pack_i16(beads, base),
        "edges": edges.tolist(),
        "cells": np.round(cells, 3).tolist(),
        "net": np.round(net, 3).tolist(),
        "xlEdge": xl_edge.tolist(), "xlAlpha": np.round(xl_alpha, 4).tolist(),
        "times": np.round(times, 1).tolist(),
        "redIds": red_ids, "redPerFrame": red_pf,
        "cellRadius": float(z["cell_radius"]) if "cell_radius" in z.files else 9.0,
        "span": float(np.max(np.abs(base)) * 1.04),
        "nCells": int(cells.shape[1]), "nFrames": int(beads.shape[0]),
        "summary": {k: (float(z[k]) if z[k].dtype.kind == "f" else int(z[k]))
                    for k in ("mean_inv", "max_inv", "detached", "lcc", "fpair", "grip_fails")
                    if k in z.files},
    }
    DATA.mkdir(parents=True, exist_ok=True)
    p = DATA / f"{tag}.js"
    p.write_text("window.G5_CASES=window.G5_CASES||{};window.G5_CASES[%s]=%s;\n"
                 % (json.dumps(tag), json.dumps(case, separators=(",", ":"))), encoding="utf-8")
    kb = p.stat().st_size / 1024
    print("wrote %-22s %6.0f KB | %d frames %d beads %d cells" % (p.name, kb, case["nFrames"], len(base), case["nCells"]))
    return {"tag": tag, "label": label, "file": f"g5-data/{tag}.js",
            "nFrames": case["nFrames"], "summary": case["summary"]}


def main():
    cases_meta = [m for m in (export_case(t, l) for t, l in CASES) if m]
    manifest = {"cases": cases_meta, "units": {"length": "um", "time": "s"}}
    (DOCS / "g5-organoid-manifest.js").write_text(
        "window.G5_MANIFEST=%s;\n" % json.dumps(manifest, separators=(",", ":")), encoding="utf-8")
    print("wrote docs/g5-organoid-manifest.js with", len(cases_meta), "cases")


if __name__ == "__main__":
    main()
