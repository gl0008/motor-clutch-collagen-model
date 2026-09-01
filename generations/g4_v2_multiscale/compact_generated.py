"""Deterministically thin already-generated overview geometry frames.

The solver and dense metric series are unchanged.  This utility selects every
second saved geometry frame from the first full G4 v2 build so the committed
site data match the production sampling intervals now declared in
``build_demo.py``: 60 s for G4C and 240 s for G4D.
"""

from __future__ import annotations

import base64
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
DATA = REPO / "docs" / "g4-v2-data"
MANIFEST = REPO / "docs" / "g4-v2-manifest.js"
SUMMARY = HERE / "generated_summary.json"


def decimate(path: Path, stride: int = 2) -> int:
    payload = json.loads(path.read_text(encoding="utf-8"))
    packed = payload["positions"]
    shape = tuple(int(x) for x in packed["shape"])
    values = np.frombuffer(base64.b64decode(packed["base64"]), dtype="<i2").reshape(shape)
    indices = list(range(0, shape[0], stride))
    if indices[-1] != shape[0] - 1:
        indices.append(shape[0] - 1)
    selected = np.ascontiguousarray(values[indices], dtype="<i2")
    packed["shape"] = list(selected.shape)
    packed["base64"] = base64.b64encode(selected.tobytes()).decode("ascii")
    payload["frames"] = [payload["frames"][i] for i in indices]
    payload["config"]["sample_interval"] = float(payload["config"]["sample_interval"]) * stride
    path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    return path.stat().st_size


def main() -> None:
    case_sizes = {}
    for case_id in ("c_independent", "c_shared", "d_fixed", "d_moving", "d_mobile_ecm"):
        case_sizes[case_id] = decimate(DATA / f"{case_id}.json")

    prefix = "window.G4V2_MANIFEST="
    raw = MANIFEST.read_text(encoding="utf-8").strip()
    manifest = json.loads(raw[len(prefix):-1])
    for case_id, size in case_sizes.items():
        manifest["cases"][case_id]["bytes"] = size
    MANIFEST.write_text(prefix + json.dumps(manifest, separators=(",", ":")) + ";\n", encoding="utf-8")

    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    for case_id, size in case_sizes.items():
        summary["cases"][case_id]["bytes"] = size
    SUMMARY.write_text(json.dumps(summary, separators=(",", ":")), encoding="utf-8")

    print(json.dumps(case_sizes, indent=2))


if __name__ == "__main__":
    main()
