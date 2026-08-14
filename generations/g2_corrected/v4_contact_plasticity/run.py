"""Generate V4 only after the saved baseline-validation gates pass."""

from pathlib import Path
import json
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from common.serialization import write_data_js  # noqa: E402
from model import PlasticityConfig, run_v4_pair  # noqa: E402


if __name__ == "__main__":
    gate_path = ROOT / "validation_summary.json"
    if not gate_path.exists():
        raise SystemExit("Run generations/g2_corrected/validate.py before V4")
    gates = json.loads(gate_path.read_text())
    required = ("v2_network", "v2_mechanics", "v3_mobility")
    passed = all(gates.get("gates", {}).get(name, False) for name in required)
    pair = run_v4_pair(PlasticityConfig(), gates_passed=passed)
    pair["validation_gates"] = gates
    write_data_js(
        HERE / "demo" / "data.js",
        pair,
        quantize_keys=("elastic", "new_links"),
    )
    elastic = pair["elastic"]
    plastic = pair["new_links"]
    summary = {
        "new_weak_links": int(plastic["new_link_count"][-1]),
        "elastic_residual_rms_um": float(elastic["rms_displacement"][-1]),
        "plastic_residual_rms_um": float(plastic["rms_displacement"][-1]),
        "elastic_residual_near_alignment": float(elastic["alignment_by_shell"][-1, 0] - elastic["alignment_by_shell"][0, 0]),
        "plastic_residual_near_alignment": float(plastic["alignment_by_shell"][-1, 0] - plastic["alignment_by_shell"][0, 0]),
    }
    (HERE / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
