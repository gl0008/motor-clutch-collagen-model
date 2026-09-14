"""Generate a self-contained G4 v4 GitHub Pages dataset.

The website loads JavaScript data directly instead of fetching JSON.  This
keeps the historical GitHub Pages behavior and also lets Safari open a local
copy without failing on ``file://`` fetch restrictions.
"""

from __future__ import annotations

import argparse
import base64
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from generations.g4_v4_single_cell_force_alignment.model import (  # noqa: E402
    G4V4Config,
    default_spec,
    run_alignment_ensemble,
    run_cavity_benchmark,
    run_fixed_cell_passive,
    run_fixed_force_direction,
    run_motor_clutch,
)


DATA_JS = REPO / "docs" / "g4-v4-data.js"  # obsolete monolith; removed after a build
DATA_DIR = REPO / "docs" / "g4-v4-data"
MANIFEST_JS = REPO / "docs" / "g4-v4-manifest.js"
SUMMARY = HERE / "generated_summary.json"


def _round(value, digits=6):
    if isinstance(value, np.ndarray):
        return np.round(value.astype(float), digits).tolist()
    if isinstance(value, (np.floating, float)):
        return round(float(value), digits)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, dict):
        return {str(key): _round(item, digits) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_round(item, digits) for item in value]
    return value


def _pack(states: np.ndarray, base: np.ndarray, *, byte: bool = False) -> dict:
    states = np.asarray(states, dtype=float)
    base = np.asarray(base, dtype=float)
    delta = states - base[None, :, :]
    maximum = float(np.max(np.abs(delta))) if delta.size else 0.0
    ceiling = 120.0 if byte else 32700.0
    scale = max(1e-6, maximum / ceiling)
    values = np.rint(delta / scale).astype("<i1" if byte else "<i2")
    return {
        "dtype": "i1" if byte else "i2",
        "scale": scale,
        "shape": list(values.shape),
        "base64": base64.b64encode(values.tobytes()).decode("ascii"),
    }


def _links(network):
    return [
        [link.edge_a, link.alpha_a, link.edge_b, link.alpha_b]
        for link in network.crosslinks
    ]


def _scalar_frame(frame):
    graph = frame.get("graph_class_metrics", {})
    bound = frame.get("bound")
    site_force = frame.get("site_force")
    actin = frame.get("actin_speed")
    substrate = frame.get("substrate_speed")
    return {
        "time": frame.get("time", 0.0),
        "cellRadius": frame.get("cell_radius"),
        "requestedForce": frame.get("requested_force", 0.0),
        "distributedForce": frame.get("distributed_force", 0.0),
        "parallelMagnitude": frame.get("parallel_magnitude", 0.0),
        "perpendicularMagnitude": frame.get("perpendicular_magnitude", 0.0),
        "deltaRadialOrder": frame.get("delta_radial_order_by_shell", [0.0, 0.0, 0.0]),
        "radialOrder": frame.get("radial_order_by_shell", [0.0, 0.0, 0.0]),
        "meanRadialAngle": frame.get("mean_radial_angle_deg"),
        "meanBoundaryAngle": frame.get("mean_boundary_angle_deg"),
        "energyRatio": frame.get("stretch_to_bend_ratio", 0.0),
        "recruitedFraction": frame.get("recruited_fraction", 0.0),
        "meanAbsAxialStrain": float(np.mean(np.abs(frame.get("bond_strain", [0.0])))),
        "meanBendingChange": float(np.mean(frame.get("bending_change", [0.0]))),
        "forceBalanceError": frame.get("force_balance_error", 0.0),
        "graph": graph,
        "boundCount": [] if bound is None else np.asarray(bound).sum(axis=1).tolist(),
        "siteForce": [] if site_force is None else np.asarray(site_force).tolist(),
        "meanActinSpeed": 0.0 if actin is None else float(np.mean(actin)),
        "meanSubstrateSpeed": 0.0 if substrate is None else float(np.mean(substrate)),
        "ruptures": frame.get("cumulative_ruptures", 0),
        "siteFailures": frame.get("cumulative_site_failures", 0),
    }


def _serialize_v4_case(result, *, case_id, title, central_object, base):
    frames = result["frames"]
    return {
        "id": case_id,
        "title": title,
        "centralObject": central_object,
        "mode": result.get("mode", case_id),
        "positions": _pack(np.asarray([frame["positions"] for frame in frames]), base),
        "times": [frame["time"] for frame in frames],
        "centers": _round([frame.get("center", [0.0, 0.0]) for frame in frames], 5),
        "contacts": _round([frame.get("contact_points", []) for frame in frames], 5),
        "forces": _round([frame.get("contact_vectors", []) for frame in frames], 5),
        "parallel": _round([frame.get("force_parallel", []) for frame in frames], 5),
        "perpendicular": _round([frame.get("force_perpendicular", []) for frame in frames], 5),
        "guidance": _round([frame.get("local_guidance_direction", []) for frame in frames], 5),
        "metrics": _round([_scalar_frame(frame) for frame in frames], 7),
        "crosslinks": _round(_links(result["network"]), 6),
        "fiberGraphDistance": _round(result.get("graph_distance", []), 0),
        "cellRadiusConstant": result.get("cell_radius_constant", True),
        "eventTime": result.get("event_time"),
        # Overview files retain site-level events only. Individual bind/rupture
        # records are kept in the high-frequency event-window file below.
        "eventLog": [
            event for event in result.get("event_log", [])
            if event.get("kind") in {"site_failure", "relocate"}
        ],
        "netDisplacement": result.get("net_displacement", 0.0),
        "pathLength": result.get("path_length", 0.0),
        "siteFailures": result.get("site_failures", 0),
        "ruptureCount": result.get("ruptures", 0),
    }


def _serialize_cavity(result, base):
    frames = result["frames"]
    return {
        "id": "cavity",
        "title": "Contractile-cavity mechanics benchmark",
        "centralObject": "contractile cavity/inclusion — not a tumor cell",
        "mode": "cavity",
        "positions": _pack(np.asarray([frame["positions"] for frame in frames]), base),
        "times": [frame["time"] for frame in frames],
        "centers": [[0.0, 0.0] for _ in frames],
        "contacts": [[] for _ in frames],
        "forces": [[] for _ in frames],
        "parallel": [[] for _ in frames],
        "perpendicular": [[] for _ in frames],
        "guidance": [[] for _ in frames],
        "metrics": _round([
            {
                "time": frame["time"],
                "cellRadius": None,
                "cavityRadius": result["config"]["cell_radius"] * (1.0 - frame["contraction"]),
                "contraction": frame["contraction"],
                "requestedForce": frame["required_inward_boundary_force"],
                "distributedForce": 0.0,
                "parallelMagnitude": 0.0,
                "perpendicularMagnitude": 0.0,
                "deltaRadialOrder": frame["radial_order_by_shell"] - frames[0]["radial_order_by_shell"],
                "radialOrder": frame["radial_order_by_shell"],
                "energyRatio": frame["stretch_to_bend_ratio"],
                "recruitedFraction": frame["recruited_fraction"],
                "meanAbsAxialStrain": float(np.mean(np.abs(frame["bond_strain"]))),
                "meanBendingChange": float(np.mean(frame.get("bending_change", [0.0]))),
                "forceBalanceError": 0.0,
                "graph": {},
                "boundCount": [],
                "siteForce": [],
                "meanActinSpeed": 0.0,
                "meanSubstrateSpeed": 0.0,
                "ruptures": 0,
                "siteFailures": 0,
            }
            for frame in frames
        ], 7),
        "crosslinks": _round(_links(result["network"]), 6),
        "fiberGraphDistance": [],
        "cellRadiusConstant": None,
        "eventLog": [],
        "netDisplacement": 0.0,
        "pathLength": 0.0,
        "siteFailures": 0,
        "ruptureCount": 0,
    }


def _serialize_event(result, base):
    frames = result.get("event_frames", [])
    if not frames:
        return None
    local = np.asarray(result["local_beads"], dtype=int)
    return {
        "time": result.get("event_time"),
        "localBeads": local.tolist(),
        "positions": _pack(
            np.asarray([frame["positions"] for frame in frames]), base[local], byte=True
        ),
        "times": [frame["time"] for frame in frames],
        "centers": _round([frame["center"] for frame in frames], 5),
        "contacts": _round([frame["contact_points"] for frame in frames], 5),
        "bound": [np.asarray(frame["bound"], dtype=np.uint8).tolist() for frame in frames],
        "boundCount": [np.asarray(frame["bound"]).sum(axis=1).tolist() for frame in frames],
        "siteForce": _round([frame["site_force"] for frame in frames], 6),
        "ruptures": [np.asarray(frame["ruptures"], dtype=np.uint8).tolist() for frame in frames],
        "binds": [np.asarray(frame["binds"], dtype=np.uint8).tolist() for frame in frames],
        "eventLog": [
            event for event in result.get("event_log", [])
            if frames[0]["time"] - 1e-9 <= event["time"] <= frames[-1]["time"] + 1e-9
        ],
        "ruptureCount": [int(np.asarray(frame["ruptures"]).sum()) for frame in frames],
        "bindCount": [int(np.asarray(frame["binds"]).sum()) for frame in frames],
    }


def _passive_case(result, base):
    network = result["network"]
    return {
        "id": "passive",
        "title": "Passive fixed-size cell",
        "centralObject": "fixed biological cell",
        "mode": "passive",
        "positions": _pack(base[None, :, :], base),
        "times": [0.0],
        "centers": [[0.0, 0.0]],
        "contacts": [[]],
        "forces": [[]],
        "parallel": [[]],
        "perpendicular": [[]],
        "guidance": [[]],
        "metrics": [{
            "time": 0.0,
            "cellRadius": result["config"]["cell_radius"],
            "requestedForce": 0.0,
            "distributedForce": 0.0,
            "parallelMagnitude": 0.0,
            "perpendicularMagnitude": 0.0,
            "deltaRadialOrder": _round(result["delta_radial_order_by_shell"], 7),
            "radialOrder": [0.0, 0.0, 0.0],
            "energyRatio": 0.0,
            "recruitedFraction": 0.0,
            "meanAbsAxialStrain": 0.0,
            "meanBendingChange": 0.0,
            "forceBalanceError": 0.0,
            "graph": {},
            "boundCount": [],
            "siteForce": [],
            "meanActinSpeed": 0.0,
            "meanSubstrateSpeed": 0.0,
            "ruptures": 0,
            "siteFailures": 0,
        }],
        "crosslinks": _round(_links(network), 6),
        "fiberGraphDistance": [],
        "cellRadiusConstant": True,
        "eventLog": [],
        "netDisplacement": 0.0,
        "pathLength": 0.0,
        "siteFailures": 0,
        "ruptureCount": 0,
    }


def build(*, quick=False):
    cfg = G4V4Config()
    if quick:
        cfg = replace(
            cfg,
            n_fibers=120,
            geometry_attempt_factor=80,
            benchmark_ramp_time=20.0,
            benchmark_relax_time=60.0,
            benchmark_sample_interval=5.0,
            fixed_duration=120.0,
            fixed_sample_interval=10.0,
            migration_duration=240.0,
            migration_sample_interval=20.0,
            minimum_event_time=5.0,
        )
    spec = default_spec(cfg)
    passive = run_fixed_cell_passive(cfg, spec=spec)
    cavity = run_cavity_benchmark(cfg, spec=spec)
    normal = run_fixed_force_direction(cfg, geometry="normal", spec=spec)
    dipole = run_fixed_force_direction(cfg, geometry="dipole", spec=spec)
    independent_clutch = run_motor_clutch(
        replace(cfg, clutch_mode="independent"),
        moving=False,
        spec=spec,
        capture_event=False,
    )
    fixed_clutch = run_motor_clutch(cfg, moving=False, spec=spec, capture_event=True)

    normal_frame = normal["frames"][-1]
    mechanics_gate = (
        normal["cell_radius_constant"]
        and normal_frame["force_balance_error"] < 1e-9
        and np.all(np.isfinite(normal_frame["delta_radial_order_by_shell"]))
    )
    clutch_gate = (
        fixed_clutch["cell_radius_constant"]
        and fixed_clutch["ruptures"] > 0
        and fixed_clutch["site_failures"] > 0
        and fixed_clutch["event_time"] is not None
    )
    moving = run_motor_clutch(cfg, moving=True, spec=spec, capture_event=False)

    # The quick ensemble checks the pipeline only. The full builder runs the
    # predeclared five independent geometries over the actual two-hour clock.
    ensemble_seeds = (41, 42, 43) if quick else (41, 42, 43, 44, 45)
    ensemble = run_alignment_ensemble(
        cfg,
        seeds=ensemble_seeds,
        duration=cfg.fixed_duration,
        sample_interval=cfg.fixed_duration,
    )

    base = np.asarray(spec.positions, dtype=float)
    cases = {
        "passive": _passive_case(passive, base),
        "cavity": _serialize_cavity(cavity, base),
        "normal": _serialize_v4_case(normal, case_id="normal", title="Fixed cell · normal surface traction", central_object="fixed biological cell", base=base),
        "dipole": _serialize_v4_case(dipole, case_id="dipole", title="Fixed cell · dipole control", central_object="fixed biological cell", base=base),
        "clutch-independent": _serialize_v4_case(independent_clutch, case_id="clutch-independent", title="Fixed cell · independent-clutch control", central_object="fixed biological cell", base=base),
        "clutch": _serialize_v4_case(fixed_clutch, case_id="clutch", title="Fixed cell · shared-load clutch cycles", central_object="fixed biological cell", base=base),
        "moving": _serialize_v4_case(moving, case_id="moving", title="Released cell · traction imbalance", central_object="moving biological cell", base=base),
    }
    cases["clutch"]["event"] = _serialize_event(fixed_clutch, base)

    data = {"cases": cases}
    manifest = {
        "version": "G4 v4 · literature-grounded single-cell force–collagen alignment",
        "quick": bool(quick),
        "actualClocks": {
            "fixedSeconds": cfg.fixed_duration,
            "movingSeconds": cfg.migration_duration,
            "fixedObservationSeconds": cfg.fixed_sample_interval,
            "movingObservationSeconds": cfg.migration_sample_interval,
        },
        "units": {"length": "µm", "force": "nN", "time": "s"},
        "domain": cfg.domain_size,
        "cellRadius": cfg.cell_radius,
        "geometry": {
            "initial": _round(spec.positions, 4),
            "fibers": spec.fibers,
            "edges": normal["network"].edges.astype(int).tolist(),
            "fixed": np.flatnonzero(spec.fixed).astype(int).tolist(),
        },
        "config": _round({
            "sourceSystem": cfg.source_system,
            "nFibers": cfg.n_fibers,
            "beadSpacing": cfg.bead_spacing,
            "crosslinkProbability": cfg.crosslink_probability,
            "crosslinkStiffness": cfg.crosslink_stiffness,
            "contactWidth": cfg.contact_width,
            "gaussianSigma": cfg.gaussian_sigma,
            "totalPullForce": cfg.total_pull_force,
            "nClutchesPerSite": cfg.n_clutches_per_site,
            "primaryClutchMode": cfg.clutch_mode,
            "clutchStiffness": cfg.clutch_stiffness,
            "clutchOnRate": cfg.clutch_on_rate,
            "clutchOffRate0": cfg.clutch_off_rate0,
            "bellForce": cfg.bell_force,
            "unloadedActinSpeed": cfg.unloaded_actin_speed,
            "motorStallPerSite": cfg.motor_stall_per_site,
        }),
        "ensemble": _round(ensemble, 7),
        "gates": {
            "passiveZero": passive["max_displacement"] < 1e-10,
            "biologicalCellRadiusConstant": all(case.get("cellRadiusConstant") is True for key, case in cases.items() if key != "cavity"),
            "mechanicsImplementation": bool(mechanics_gate),
            "clutchRuptureObserved": bool(clutch_gate),
            "completeSiteFailureObserved": bool(
                fixed_clutch["site_failures"] > 0 and fixed_clutch["event_time"] is not None
            ),
            "alignmentExcludesPassive": bool(ensemble["near_alignment_excludes_passive_zero"]),
            "oneHopExceedsUnconnected": bool(ensemble["one_hop_exceeds_unconnected"]),
            "movingScientificStatus": "diagnostic only until all upstream biological gates pass",
        },
        "assumptions": [f"V4-A{index:03d}" for index in range(1, 18)],
        "caseFiles": {
            case_id: f"g4-v4-data/{case_id}.js" for case_id in cases
        },
        "eventFiles": {"clutch": "g4-v4-data/clutch-event.js"},
    }
    MANIFEST_JS.write_text(
        "window.G4V4_MANIFEST=" + json.dumps(_round(manifest), separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    event = data["cases"]["clutch"].pop("event", None)
    for case_id, case in data["cases"].items():
        payload = (
            "window.G4V4_CASES=window.G4V4_CASES||{};"
            f"window.G4V4_CASES[{json.dumps(case_id)}]="
            + json.dumps(_round(case), separators=(",", ":"))
            + ";\n"
        )
        (DATA_DIR / f"{case_id}.js").write_text(payload, encoding="utf-8")
    (DATA_DIR / "clutch-event.js").write_text(
        "window.G4V4_EVENTS=window.G4V4_EVENTS||{};window.G4V4_EVENTS.clutch="
        + json.dumps(_round(event), separators=(",", ":"))
        + ";\n",
        encoding="utf-8",
    )
    DATA_JS.unlink(missing_ok=True)
    summary = {
        "version": manifest["version"],
        "quick": bool(quick),
        "actual_clocks": manifest["actualClocks"],
        "gates": manifest["gates"],
        "ensemble": manifest["ensemble"],
        "cases": {
            key: {
                "frames": len(case["times"]),
                "duration": case["times"][-1],
                "ruptures": case.get("ruptureCount", 0),
                "site_failures": case.get("siteFailures", 0),
                "net_displacement": case.get("netDisplacement", 0.0),
            }
            for key, case in cases.items()
        },
    }
    SUMMARY.write_text(json.dumps(_round(summary), indent=2) + "\n", encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="short, explicitly labeled pipeline check")
    args = parser.parse_args()
    print(json.dumps(build(quick=args.quick), indent=2))


if __name__ == "__main__":
    main()
