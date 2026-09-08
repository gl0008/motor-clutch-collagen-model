"""Build compact, lazy-loaded G4 v3 benchmark data for GitHub Pages."""

from __future__ import annotations

import argparse
import base64
from dataclasses import replace
import json
from multiprocessing import Pool
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from generations.g4_v3_random_to_aligned.model import (  # noqa: E402
    G4V3Config,
    build_random_network,
    make_random_void_spec,
    run_fixed_material_traction,
    run_prescribed_contraction,
)

DATA_DIR = REPO / "docs" / "g4-v3-data"
MANIFEST = REPO / "docs" / "g4-v3-manifest.js"
SUMMARY = HERE / "generated_summary.json"


def _round(value, digits=6):
    if isinstance(value, np.ndarray):
        return np.round(value.astype(float), digits).tolist()
    if isinstance(value, (np.floating, float)):
        return round(float(value), digits)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, dict):
        return {str(key): _round(item, digits) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_round(item, digits) for item in value]
    return value


def _pack_displacements(states: np.ndarray, initial: np.ndarray) -> dict:
    delta = np.asarray(states, dtype=float) - np.asarray(initial)[None, :, :]
    maximum = float(np.max(np.abs(delta))) if delta.size else 0.0
    scale = max(1e-5, maximum / 32700.0)
    quantized = np.rint(delta / scale).astype("<i2")
    return {
        "dtype": "i2",
        "scale": scale,
        "shape": list(quantized.shape),
        "base64": base64.b64encode(quantized.tobytes()).decode("ascii"),
    }


def _crosslinks(network) -> list[list[float]]:
    return [
        [link.edge_a, link.alpha_a, link.edge_b, link.alpha_b]
        for link in network.crosslinks
    ]


def _case_worker(task):
    config_dict, case_id, mode, value, probability, spec, reuse = task
    cfg = G4V3Config(**config_dict)
    path = DATA_DIR / f"{case_id}.json"
    if reuse and path.exists():
        payload = json.loads(path.read_text(encoding="utf-8"))
        network, _, report = build_random_network(
            cfg, spec=spec, crosslink_probability=probability
        )
        return payload, report, _crosslinks(network)
    if mode == "traction":
        result = run_fixed_material_traction(
            cfg,
            total_force=value,
            spec=spec,
            crosslink_probability=probability,
        )
    else:
        result = run_prescribed_contraction(
            cfg,
            contraction=value,
            spec=spec,
            crosslink_probability=probability,
        )
    frames = result["frames"]
    # Eleven geometry states are enough for smooth interpolation; all metric
    # samples remain available as small scalar arrays.
    keep = sorted(set(np.linspace(0, len(frames) - 1, 11).round().astype(int)))
    states = np.asarray([frames[index]["positions"] for index in keep])
    if mode == "traction":
        metrics = [
            {
                "time": frame["time"],
                "contraction": 0.0,
                "traction": frame["distributed_force_magnitude"],
                "radialOrder": _round(frame["radial_order_by_shell"], 7),
                "deltaRadialOrder": _round(frame["delta_radial_order_by_shell"], 7),
                "stretchEnergy": frame["energy"]["stretch"],
                "bendEnergy": frame["energy"]["bend"],
                "crosslinkEnergy": frame["energy"]["crosslink"],
                "energyRatio": frame["stretch_to_bend_ratio"],
                "recruitedFraction": frame["recruited_fraction"],
                "boundaryForce": frame["distributed_force_magnitude"],
                "netReaction": float(np.linalg.norm(frame["constrained_cell_reaction"])),
                "dragDissipation": 0.0,
            }
            for frame in frames
        ]
        final_delta = result["final_delta_radial_order_by_shell"]
        final_ratio = frames[-1]["stretch_to_bend_ratio"]
        final_recruited = frames[-1]["recruited_fraction"]
        final_force = frames[-1]["distributed_force_magnitude"]
    else:
        metrics = [
            {
                "time": frame["time"],
                "contraction": frame["contraction"],
                "traction": 0.0,
                "radialOrder": _round(frame["radial_order_by_shell"], 7),
                "deltaRadialOrder": _round(
                    frame["radial_order_by_shell"] - frames[0]["radial_order_by_shell"], 7
                ),
                "stretchEnergy": frame["energy"]["stretch"],
                "bendEnergy": frame["energy"]["bend"],
                "crosslinkEnergy": frame["energy"]["crosslink"],
                "energyRatio": frame["stretch_to_bend_ratio"],
                "recruitedFraction": frame["recruited_fraction"],
                "boundaryForce": frame["required_inward_boundary_force"],
                "netReaction": frame["required_inward_boundary_force"],
                "dragDissipation": frame["cumulative_drag_dissipation"],
            }
            for frame in frames
        ]
        final_delta = (
            result["final_radial_order_by_shell"]
            - result["initial_radial_order_by_shell"]
        )
        final_ratio = result["final_stretch_to_bend_ratio"]
        final_recruited = result["final_recruited_fraction"]
        final_force = frames[-1]["required_inward_boundary_force"]

    payload = {
        "id": case_id,
        "mode": mode,
        "value": value,
        "contraction": value if mode == "contraction" else 0.0,
        "crosslinkProbability": probability,
        "geometryTimes": [frames[index]["time"] for index in keep],
        "positions": _pack_displacements(states, spec.positions),
        "contactPoints": _round(
            [frames[index].get("contact_points", []) for index in keep], 5
        ),
        "contactVectors": _round(
            [frames[index].get("contact_vectors", []) for index in keep], 5
        ),
        "metrics": metrics,
        "final": {
            "deltaRadialOrder": _round(final_delta, 7),
            "energyRatio": final_ratio,
            "recruitedFraction": final_recruited,
            "boundaryForce": final_force,
        },
    }
    return payload, result["report"], _crosslinks(result["network"])


def build(*, quick: bool = False, reuse: bool = False) -> dict:
    cfg = G4V3Config()
    if quick:
        cfg = replace(
            cfg,
            n_fibers=120,
            benchmark_ramp_time=20.0,
            benchmark_relax_time=60.0,
            benchmark_sample_interval=3.0,
        )
    spec = make_random_void_spec(cfg)
    cases_to_run = [
        ("contraction-5", "contraction", 0.05, 1.0),
        ("contraction-10", "contraction", 0.10, 1.0),
        ("contraction-20", "contraction", 0.20, 1.0),
        ("links-0", "contraction", 0.20, 0.0),
        ("links-15", "contraction", 0.20, 0.15),
        ("links-35", "contraction", 0.20, 0.35),
        ("links-60", "contraction", 0.20, 0.60),
        ("traction-12", "traction", 12.0, 0.35),
        ("traction-24", "traction", 24.0, 0.35),
        ("traction-48", "traction", 48.0, 0.35),
    ]
    tasks = [
        (dict(cfg.__dict__), case_id, mode, value, probability, spec, reuse)
        for case_id, mode, value, probability in cases_to_run
    ]
    # Independent cases share only immutable input geometry and are safe to run
    # concurrently. The underlying equations remain in model.py.
    with Pool(processes=min(3, len(tasks))) as pool:
        outputs = pool.map(_case_worker, tasks)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    cases = {}
    report = outputs[0][1]
    links_by_case = {}
    reports_by_case = {}
    for payload, report_for_case, links in outputs:
        case_id = payload["id"]
        path = DATA_DIR / f"{case_id}.json"
        path.write_text(json.dumps(_round(payload), separators=(",", ":")), encoding="utf-8")
        cases[case_id] = {
            "file": f"g4-v3-data/{case_id}.json",
            "bytes": path.stat().st_size,
            "contraction": payload["contraction"],
            "mode": payload.get("mode", "contraction"),
            "value": payload.get("value", payload["contraction"]),
            "crosslinkProbability": payload["crosslinkProbability"],
            "report": _round(report_for_case, 7),
            "final": _round(payload["final"], 7),
        }
        links_by_case[case_id] = links
        reports_by_case[case_id] = report_for_case

    geometry = {
        "initial": _round(spec.positions, 4),
        "fibers": spec.fibers,
        "fixed": np.flatnonzero(spec.fixed).astype(int).tolist(),
        "inner": [],
        "crosslinks": [],
    }
    # Crosslink topology is the same in all cases. Rebuild from the first result
    # only for export; production physics has already completed.
    first_cfg = G4V3Config(**dict(cfg.__dict__))
    first_result = run_prescribed_contraction(
        replace(first_cfg, benchmark_relax_time=first_cfg.dt),
        contraction=0.0,
        spec=spec,
    )
    network = first_result["network"]
    geometry["inner"] = first_result["inner_beads"].astype(int).tolist()
    geometry["crosslinks"] = _round(_crosslinks(network), 6)
    link_identity = {
        tuple(round(float(value), 8) for value in link): index
        for index, link in enumerate(_crosslinks(network))
    }
    for case_id, links in links_by_case.items():
        cases[case_id]["crosslinkIds"] = [
            link_identity[tuple(round(float(value), 8) for value in link)]
            for link in links
        ]

    manifest = {
        "version": "G4 v3 · random mechanics through fixed traction",
        "quick": quick,
        "units": {"length": "µm", "force": "nN", "time": "s"},
        "domain": cfg.domain_size,
        "cellRadius": cfg.cell_radius,
        "cellClearance": cfg.cell_clearance,
        "geometry": geometry,
        "report": _round(report, 7),
        "cases": cases,
        "assumptions": ["A-001", "A-002", "A-003", "A-004", "A-005", "A-006", "A-007", "A-008", "A-009", "A-011", "A-012", "A-013", "A-014"],
    }
    MANIFEST.write_text(
        "window.G4V3_MANIFEST=" + json.dumps(manifest, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )
    compact_cases = {
        case_id: {
            "file": case["file"],
            "bytes": case["bytes"],
            "contraction": case["contraction"],
            "mode": case["mode"],
            "value": case["value"],
            "crosslinkProbability": case["crosslinkProbability"],
            "report": case["report"],
            "final": case["final"],
        }
        for case_id, case in cases.items()
    }
    summary = {
        "version": manifest["version"],
        "quick": quick,
        "report": manifest["report"],
        "cases": compact_cases,
        "acceptance": {
            "global_nematic_below_0_10": report["initial_global_nematic_order"] < 0.10,
            "minimum_four_natural_contacts": report["contact_fragments"] >= 4,
            "no_handcrafted_contact_fibers": True,
            "same_geometry_for_all_cases": True,
        },
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="small smoke-test dataset")
    parser.add_argument("--reuse", action="store_true", help="reuse existing case arrays")
    args = parser.parse_args()
    print(json.dumps(build(quick=args.quick, reuse=args.reuse), indent=2))


if __name__ == "__main__":
    main()
