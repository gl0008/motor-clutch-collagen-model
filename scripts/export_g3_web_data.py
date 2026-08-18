"""Export solved G3 snapshots for the Gloria G1/G2-style SVG renderer."""

from __future__ import annotations

import base64
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from g3.run import load_saved_run  # noqa: E402
from g3.visualization import (  # noqa: E402
    _contact_fibre_ids,
    _fibre_reorientation_deg,
    _mechanism_phase,
    _protrusion_switches,
)


def _round(value, digits=7):
    return np.round(np.asarray(value, dtype=float), digits).tolist()


def _encode_positions(result):
    initial_um = result.initial_positions * 1.0e6
    positions_um = np.asarray([frame.positions for frame in result.snapshots]) * 1.0e6
    delta = positions_um - initial_um[None, :, :]
    maximum = float(np.abs(delta).max(initial=0.0))
    quantum = max(1.0e-7, maximum / 30000.0)
    quantized = np.rint(delta / quantum).astype("<i2")
    return {
        "quantum_um": quantum,
        "shape": list(quantized.shape),
        "base64": base64.b64encode(quantized.tobytes()).decode("ascii"),
    }


def _serialize(result):
    network = result.fixture.network
    snapshots = result.snapshots
    fibers = [
        np.flatnonzero(network.fiber_id == fiber).astype(int).tolist()
        for fiber in range(network.n_fibers)
    ]
    edges = network.fiber_bonds.astype(int)
    crosslinks = []
    if network.xl_edge_a is not None:
        crosslinks = [
            {
                "edge_a": int(a), "edge_b": int(b),
                "alpha_a": float(aa), "alpha_b": float(ab),
            }
            for a, b, aa, ab in zip(
                network.xl_edge_a, network.xl_edge_b,
                network.xl_alpha_a, network.xl_alpha_b,
            )
        ]
    contact_ids = sorted(_contact_fibre_ids(result))
    metrics = {
        "phase": [], "load_nN": [], "angle_mdeg": [], "max_bead_dr_um": [],
        "switches": [], "cell_dr_um": [], "cell_rotation_deg": [],
    }
    start_center = snapshots[0].cell_center
    start_angle = snapshots[0].cell_angle
    for index, frame in enumerate(snapshots):
        attached = set(map(int, frame.bound_fibre_ids))
        metrics["phase"].append(_mechanism_phase(frame))
        metrics["load_nN"].append(float(np.linalg.norm(frame.clutch_forces, axis=1).sum() * 1.0e9))
        metrics["angle_mdeg"].append(1000.0 * _fibre_reorientation_deg(
            result, frame, attached if attached else contact_ids,
        ))
        metrics["max_bead_dr_um"].append(float(np.linalg.norm(
            frame.positions - result.initial_positions, axis=1,
        ).max(initial=0.0) * 1.0e6))
        metrics["switches"].append(_protrusion_switches(result, index))
        metrics["cell_dr_um"].append(float(np.linalg.norm(frame.cell_center - start_center) * 1.0e6))
        metrics["cell_rotation_deg"].append(float(np.degrees(frame.cell_angle - start_angle)))
    return {
        "stage": result.stage,
        "fixture": result.fixture_name,
        "config": {
            "domain_size": result.config.scaled_domain_size * 1.0e6,
            "cell_radius": result.config.cell_radius * 1.0e6,
            "n_sectors": result.config.n_sectors,
        },
        "initial_positions": _round(result.initial_positions * 1.0e6, 6),
        "position_encoding": _encode_positions(result),
        "fibers": fibers,
        "bead_fiber": network.fiber_id.astype(int).tolist(),
        "edges": edges.tolist(),
        "edge_fiber": network.fiber_id[edges[:, 0]].astype(int).tolist(),
        "fixed": result.fixture.fixed_mask.astype(bool).tolist(),
        "crosslinks": crosslinks,
        "contact_fibers": contact_ids,
        "time": _round([frame.time for frame in snapshots], 4),
        "cell_center": _round([frame.cell_center * 1.0e6 for frame in snapshots], 7),
        "cell_angle": _round([frame.cell_angle for frame in snapshots], 9),
        "active_sectors": [list(map(int, frame.active_sectors)) for frame in snapshots],
        "activity": _round([frame.polarity_activity for frame in snapshots], 7),
        "protrusion_lengths": _round([frame.protrusion_lengths * 1.0e6 for frame in snapshots], 7),
        "protrusion_tips": _round([frame.protrusion_tips * 1.0e6 for frame in snapshots], 7),
        "bound_points": [_round(frame.bound_points * 1.0e6, 7) for frame in snapshots],
        "motor_points": [_round(frame.motor_points * 1.0e6, 7) for frame in snapshots],
        "clutch_forces_nN": [_round(frame.clutch_forces * 1.0e9, 7) for frame in snapshots],
        "bound_fibre_ids": [list(map(int, frame.bound_fibre_ids)) for frame in snapshots],
        "bound_sector_ids": [list(map(int, frame.bound_sector_ids)) for frame in snapshots],
        "metrics": {key: _round(value, 7) if key != "phase" else value
                    for key, value in metrics.items()},
    }


def main():
    results = ROOT / "results" / "g3_revision"
    runs = {
        "g3a": load_saved_run(results / "g3a_multifibre"),
        "g3b_off": load_saved_run(results / "g3b_feedback_off"),
        "g3b_on": load_saved_run(results / "g3b_feedback_on"),
        "g3c_fixed": load_saved_run(results / "g3c_fixed_control"),
        "g3c_released": load_saved_run(results / "g3c_released"),
    }
    payload = {name: _serialize(run) for name, run in runs.items()}
    target = ROOT / "docs" / "g3-web-data.js"
    target.write_text(
        "window.G3_WEB_DATA=" + json.dumps(payload, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": str(target),
        "bytes": target.stat().st_size,
        "runs": {name: len(run.snapshots) for name, run in runs.items()},
    }, indent=2))


if __name__ == "__main__":
    main()
