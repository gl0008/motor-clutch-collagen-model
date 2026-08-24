"""Precompute every exact trajectory used by the static G4 web laboratory."""

from __future__ import annotations

import argparse
import base64
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from model import G4Config, _assemble_g4_spec, run_clutch_pair, run_elastic

HERE = Path(__file__).resolve().parent


def _rounded(values, digits=5):
    return np.round(np.asarray(values, dtype=float), digits).tolist()


def _pack_positions(frames: list[dict], initial: np.ndarray) -> dict:
    scale = 1.0e-4
    positions = np.asarray([frame["positions"] for frame in frames])
    quantized = np.rint((positions - initial[None, :, :]) / scale)
    if np.max(np.abs(quantized)) > np.iinfo(np.int16).max:
        raise ValueError("G4 displacement exceeds int16 web packing range")
    raw = np.asarray(quantized, dtype="<i2").tobytes()
    return {
        "scale": scale,
        "shape": list(positions.shape),
        "base64": base64.b64encode(raw).decode("ascii"),
    }


def _links(result: dict) -> list[list[float]]:
    return [
        [x.edge_a, round(x.alpha_a, 5), x.edge_b, round(x.alpha_b, 5)]
        for x in result["network"].crosslinks
    ]


def _case(result: dict, initial: np.ndarray, label: str, stage: str) -> dict:
    frames = result["frames"]
    return {
        "label": label,
        "stage": stage,
        "config": result["config"],
        "positions": _pack_positions(frames, initial),
        "crosslinks": _links(result),
        "graphDistance": np.asarray(result["graph_distance"], dtype=int).tolist(),
        "directFibers": [int(x) for x in result["direct_fibers"]],
        "time": _rounded([frame["time"] for frame in frames], 3),
        "centers": _rounded([frame["center"] for frame in frames], 5),
        "theta": _rounded([frame["theta"] for frame in frames], 6),
        "contacts": [_rounded(frame["contact_points"], 4) for frame in frames],
        "vectors": [_rounded(frame["contact_vectors"], 4) for frame in frames],
        "traction": _rounded([frame["traction"] for frame in frames], 4),
        "bound": [int(frame["bound"]) for frame in frames],
        "slips": [int(frame["slips"]) for frame in frames],
        "tractionDrops": [int(frame["traction_drops"]) for frame in frames],
        "recoilEvents": [int(frame["recoil_events"]) for frame in frames],
        "radialOrder": _rounded([frame["radial_order"] for frame in frames], 6),
        "graphDisplacement": [frame["graph_displacement"] for frame in frames],
        "report": result["report"],
        "pathLength": round(float(result.get("path_length", 0.0)), 5),
        "netDisplacement": round(float(result.get("net_displacement", 0.0)), 5),
        "slipCount": int(len(result.get("slip_times", []))),
        "tractionDropCount": int(result.get("traction_drop_events", 0)),
        "recoilCount": int(result.get("recoil_events", 0)),
        "maxDisplacement": round(float(result.get("max_displacement", 0.0)), 5),
    }


def build(output: Path) -> None:
    baseline = G4Config(duration=120.0, sample_interval=4.0)
    spec = _assemble_g4_spec(baseline, baseline.seed)
    initial = spec.positions
    cases: dict[str, dict] = {}

    axes = {
        "pull": ("total_pull_force", [12.0, 24.0, 48.0]),
        "bend": ("bending_multiplier", [0.10, 0.25, 0.50, 1.00]),
        "drag": ("bead_drag", [90.0, 180.0, 360.0]),
        "probability": ("crosslink_probability", [0.00, 0.20, 0.35, 0.60]),
        "link": ("crosslink_stiffness", [5.0, 10.0, 25.0, 75.0]),
    }
    seen: dict[tuple, str] = {}
    for short, (field, values) in axes.items():
        for value in values:
            cfg = replace(baseline, **{field: value})
            signature = (
                cfg.total_pull_force, cfg.bending_multiplier, cfg.bead_drag,
                cfg.crosslink_probability, cfg.crosslink_stiffness,
            )
            key = f"a_{short}_{str(value).replace('.', 'p')}"
            if signature in seen:
                cases[key] = {"alias": seen[signature], "label": f"{field} = {value}", "stage": "G4A"}
                continue
            result = run_elastic(cfg, spec=spec)
            cases[key] = _case(result, initial, f"{field} = {value}", "G4A")
            seen[signature] = key

    # G4B intentionally combines the accepted more-mobile settings so indirect
    # transmission is measurable, then changes the link or mobility condition.
    b_base = replace(
        baseline, duration=240.0, sample_interval=4.0,
        total_pull_force=48.0, crosslink_probability=0.60,
        crosslink_stiffness=25.0, bead_drag=180.0,
    )
    b_cases = {
        "b_no_links": replace(b_base, crosslink_probability=0.0),
        "b_mobile": replace(b_base, bead_drag=90.0),
        "b_baseline": b_base,
        "b_constrained": replace(b_base, bead_drag=360.0),
    }
    for key, cfg in b_cases.items():
        result = run_elastic(cfg, spec=spec)
        label = {
            "b_no_links": "no crosslinks (negative control)",
            "b_mobile": "linked + low drag",
            "b_baseline": "linked baseline",
            "b_constrained": "linked + high drag",
        }[key]
        cases[key] = _case(result, initial, label, "G4B")

    # G4C/D use one identical geometry, link set and counter-addressed stream.
    clutch_cfg = replace(
        baseline, duration=300.0, sample_interval=5.0,
        total_pull_force=24.0, crosslink_probability=0.60,
        crosslink_stiffness=25.0, bead_drag=180.0,
    )
    pair = run_clutch_pair(clutch_cfg, spec=spec)
    cases["c_fixed"] = _case(pair["fixed"], initial, "G4C · fixed cell", "G4C")
    cases["d_moving"] = _case(pair["moving"], initial, "G4D · moving cell", "G4D")

    data = {
        "units": {"length": "µm", "force": "nN", "time": "s"},
        "domain": baseline.domain_size,
        "cellRadius": baseline.cell_radius,
        "initial": _rounded(initial, 4),
        "edges": np.asarray(pair["fixed"]["network"].edges, dtype=int).tolist(),
        "edgeFiber": np.asarray(pair["fixed"]["network"].edge_fiber, dtype=int).tolist(),
        "fixed": np.flatnonzero(pair["fixed"]["network"].fixed).astype(int).tolist(),
        "fibers": [np.asarray(ids, dtype=int).tolist() for ids in pair["fixed"]["network"].fibers],
        "axes": {
            short: {"field": field, "values": values}
            for short, (field, values) in axes.items()
        },
        "baselineKeys": {
            "pull": "a_pull_24p0", "bend": "a_bend_0p25",
            "drag": "a_drag_180p0", "probability": "a_probability_0p35",
            "link": "a_link_10p0",
        },
        "bell": {
            "kClutch": clutch_cfg.clutch_stiffness,
            "kOff0": clutch_cfg.clutch_off_rate0,
            "bellForce": clutch_cfg.bell_force,
            "nClutches": clutch_cfg.n_clutches_per_site,
            "zeroMedian": round(np.log(2) / clutch_cfg.clutch_off_rate0, 2),
            "oneMedian": round(np.log(2) / (clutch_cfg.clutch_off_rate0 * np.e), 2),
            "twoMedian": round(np.log(2) / (clutch_cfg.clutch_off_rate0 * np.e**2), 2),
        },
        "cases": cases,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("window.G4_DATA=" + json.dumps(data, separators=(",", ":")) + ";\n", encoding="utf-8")
    print(f"wrote {output} ({output.stat().st_size / 1e6:.2f} MB, {len(cases)} cases)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=HERE.parent.parent / "docs" / "g4-demo-data.js")
    args = parser.parse_args()
    build(args.out)


if __name__ == "__main__":
    main()
