"""Build lazy-loaded static datasets for the G4 v2 web laboratory.

The Python model remains the only physics implementation.  The browser reads a
small manifest first, then fetches one exact precomputed case at a time.
"""

from __future__ import annotations

import argparse
import base64
from dataclasses import asdict, replace
import json
from multiprocessing import Pool
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from generations.g4_interactive_calibration.model import _assemble_g4_spec  # noqa: E402
from generations.g4_v2_multiscale.model import (  # noqa: E402
    G4V2Config,
    build_v2_network,
    run_clutch,
    run_elastic,
    shared_cluster_ensemble,
    with_counter_seed,
)


DEFAULT_DATA = REPO / "docs" / "g4-v2-data"
DEFAULT_MANIFEST = REPO / "docs" / "g4-v2-manifest.js"
DEFAULT_SUMMARY = HERE / "generated_summary.json"


def _round(value, digits=6):
    if isinstance(value, np.ndarray):
        return np.round(value.astype(float), digits).tolist()
    if isinstance(value, (np.floating, float)):
        return round(float(value), digits)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, dict):
        return {str(k): _round(v, digits) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_round(v, digits) for v in value]
    return value


def _pack_signed(values: np.ndarray, base: np.ndarray, minimum_scale=1e-4) -> dict:
    values = np.asarray(values, dtype=float)
    delta = values - np.asarray(base, dtype=float)[None, :, :]
    maximum = float(np.max(np.abs(delta))) if delta.size else 0.0
    scale = max(float(minimum_scale), maximum / 32700.0)
    quantized = np.rint(delta / scale).astype("<i2")
    return {
        "dtype": "i2", "scale": scale, "shape": list(values.shape),
        "base64": base64.b64encode(quantized.tobytes()).decode("ascii"),
    }


def _links(result: dict) -> list[list[float]]:
    return [
        [x.edge_a, round(x.alpha_a, 5), x.edge_b, round(x.alpha_b, 5)]
        for x in result["network"].crosslinks
    ]


def _record(frame: dict, *, geometry=False, compact_metric=False) -> dict:
    bound = frame.get("bound")
    site_force = frame.get("site_force")
    out = {
        "time": round(float(frame["time"]), 4),
        "center": _round(frame["center"], 5),
        "theta": round(float(frame["theta"]), 7),
        "traction": round(float(frame.get("traction", 0.0)), 5),
        "boundBySite": [] if bound is None else np.asarray(bound).sum(axis=1).astype(int).tolist(),
        "slips": int(frame.get("cumulative_slips", 0)),
        "siteFailures": int(frame.get("site_failures", 0)),
        "radialOrder": round(float(frame.get("radial_order", 0.0)), 7),
    }
    if not compact_metric:
        out["siteForce"] = [] if site_force is None else _round(site_force, 5)
        out["graphStats"] = _round(frame.get("graph_stats", {}), 7)
    if geometry:
        out.update({
            "contacts": _round(frame.get("contact_points", []), 4),
            "vectors": _round(frame.get("contact_vectors", []), 4),
            "fiberDeltaTheta": _round(frame.get("fiber_delta_theta", []), 6),
        })
    return out


def _case_payload(result: dict, initial: np.ndarray, case_id: str, label: str) -> dict:
    frames = result["frames"]
    compact_metrics = result.get("mode", "elastic") != "elastic"
    payload = {
        "id": case_id, "label": label, "config": result["config"],
        "positions": _pack_signed(
            np.asarray([f["positions"] for f in frames]), initial,
        ),
        "frames": [_record(f, geometry=True) for f in frames],
        # Clutch metrics are sampled much more densely than geometry.  Keep only
        # the fields used by the traction/cell-path plots; graph summaries and
        # individual site forces remain in sparse geometry and event frames.
        "metrics": [_record(f, compact_metric=compact_metrics) for f in result["metrics"]],
        "crosslinks": _links(result),
        "fixed": np.flatnonzero(result["network"].fixed).astype(int).tolist(),
        "graphDistance": np.asarray(result["graph_distance"], dtype=int).tolist(),
        "directFibers": [int(x) for x in result["direct_fibers"]],
        "report": _round(result["report"], 6),
        "finalGraphStats": _round(result["final_graph_stats"], 7),
        "maxDisplacement": round(float(result["max_displacement"]), 7),
        "mode": result.get("mode", "elastic"),
        "moving": bool(result.get("moving", False)),
        "pathLength": round(float(result.get("path_length", 0.0)), 7),
        "netDisplacement": round(float(result.get("net_displacement", 0.0)), 7),
        "siteFailureTimes": _round(result.get("site_failure_times", []), 4),
    }
    return payload


def _event_payload(result: dict, initial: np.ndarray, case_id: str, label: str) -> dict | None:
    frames = result.get("event_frames", [])
    if not frames:
        return None
    local = np.asarray(result["local_beads"], dtype=int)
    positions = np.asarray([f["positions"] for f in frames])
    return {
        "id": case_id, "label": label,
        "eventTime": None if result.get("event_time") is None else round(float(result["event_time"]), 4),
        "eventSite": int(result.get("event_site", -1)),
        "positions": _pack_signed(positions, initial[local]),
        "localInitial": _round(initial[local], 4),
        "localBeads": local.tolist(),
        "localEdges": np.asarray(result["local_edges"], dtype=int).tolist(),
        "localEdgeFiber": np.asarray(result["network"].edge_fiber[result["local_edge_ids"]], dtype=int).tolist(),
        "graphDistance": np.asarray(result["graph_distance"], dtype=int).tolist(),
        "frames": [{
            "time": round(float(f["time"]), 4),
            "center": _round(f["center"], 5),
            "theta": round(float(f["theta"]), 7),
            "contacts": _round(f["contact_points"], 4),
            "bound": np.asarray(f["bound"], dtype=np.uint8).tolist(),
            "boundBySite": np.asarray(f["bound"], dtype=np.uint8).sum(axis=1).astype(int).tolist(),
            "siteForce": _round(f["site_force"], 5),
            "breaks": np.asarray(f["breaks"], dtype=np.uint8).tolist(),
            "binds": np.asarray(f["binds"], dtype=np.uint8).tolist(),
            "substrateSpeed": _round(f["substrate_speed"], 7),
        } for f in frames],
    }


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")


def _elastic_worker(task):
    case_id, label, config_dict, data_dir, reuse = task
    path = Path(data_dir) / f"{case_id}.json"
    if reuse and path.exists():
        payload = json.loads(path.read_text(encoding="utf-8"))
        return case_id, label, {
            "file": f"g4-v2-data/{case_id}.json",
            "bytes": path.stat().st_size,
            "maxDisplacement": payload["maxDisplacement"],
            "finalGraphStats": payload["finalGraphStats"],
        }
    cfg = G4V2Config(**config_dict)
    spec = _assemble_g4_spec(cfg, cfg.seed)
    result = run_elastic(cfg, spec=spec)
    payload = _case_payload(result, spec.positions, case_id, label)
    _write_json(path, payload)
    return case_id, label, {
        "file": f"g4-v2-data/{case_id}.json",
        "bytes": path.stat().st_size,
        "maxDisplacement": payload["maxDisplacement"],
        "finalGraphStats": payload["finalGraphStats"],
    }


def _clutch_worker(task):
    case_id, label, config_dict, mode, moving, data_dir = task
    cfg = G4V2Config(**config_dict)
    spec = _assemble_g4_spec(cfg, cfg.seed)
    result = run_clutch(cfg, spec=spec, mode=mode, moving=moving)
    payload = _case_payload(result, spec.positions, case_id, label)
    path = Path(data_dir) / f"{case_id}.json"
    _write_json(path, payload)
    return case_id, label, {
        "file": f"g4-v2-data/{case_id}.json",
        "bytes": path.stat().st_size,
        "netDisplacement": payload["netDisplacement"],
        "pathLength": payload["pathLength"],
        "siteFailures": int(payload["metrics"][-1]["siteFailures"]),
        "slips": int(payload["metrics"][-1]["slips"]),
    }


def _axis_cases(base: G4V2Config):
    axes = {
        "pull": ("total_pull_force", [12.0, 24.0, 48.0], "nN"),
        "bend": ("bending_multiplier", [0.10, 0.25, 1.00], "× EI₀"),
        "stretch": ("collagen_modulus_mpa", [1.5, 3.0, 6.0], "MPa"),
        "drag": ("bead_drag", [90.0, 180.0, 360.0], "nN·s/µm"),
        "probability": ("crosslink_probability", [0.0, 0.35, 0.60], "fraction"),
        "link": ("crosslink_stiffness", [5.0, 10.0, 25.0], "nN/µm"),
        "sigma": ("gaussian_sigma", [0.75, 1.50, 2.25], "µm"),
    }
    cases = {"a_baseline": ("baseline", base)}
    manifest_axes = {}
    for short, (field, values, unit) in axes.items():
        entries = []
        baseline_value = getattr(base, field)
        for value in values:
            if value == baseline_value:
                case_id = "a_baseline"
            else:
                case_id = f"a_{short}_{str(value).replace('.', 'p')}"
                cases[case_id] = (f"{field} = {value} {unit}", replace(base, **{field: value}))
            entries.append({"value": value, "case": case_id})
        manifest_axes[short] = {
            "field": field, "unit": unit, "baseline": baseline_value, "entries": entries,
        }
    cases["a_mobile_boundary"] = ("mobile outer boundary control", replace(base, boundary_mode="mobile"))
    return axes, cases, manifest_axes


def build(*, data_dir: Path, manifest_path: Path, summary_path: Path,
          quick: bool = False, workers: int = 1, reuse_elastic: bool = False) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    if quick:
        elastic_duration, elastic_sample = 60.0, 5.0
        c_duration, c_geometry, c_metric = 90.0, 3.0, 0.5
        d_duration, d_geometry, d_metric = 120.0, 4.0, 1.0
        ensemble_trials, ensemble_duration = 20, 60.0
    else:
        elastic_duration, elastic_sample = 7200.0, 240.0
        c_duration, c_geometry, c_metric = 1800.0, 60.0, 1.0
        d_duration, d_geometry, d_metric = 7200.0, 240.0, 10.0
        ensemble_trials, ensemble_duration = 200, 1800.0

    a_base = G4V2Config(
        duration=elastic_duration, sample_interval=elastic_sample,
        metric_sample_interval=elastic_sample, boundary_mode="anchored",
    )
    spec = _assemble_g4_spec(a_base, a_base.seed)
    baseline_network, _, _ = build_v2_network(a_base, spec=spec)
    _, a_cases, axes = _axis_cases(a_base)

    b_base = replace(
        a_base, total_pull_force=48.0, crosslink_probability=0.60,
        crosslink_stiffness=25.0,
    )
    b_cases = {
        "b_no_links": ("same geometry · no crosslinks", replace(b_base, crosslink_probability=0.0)),
        "b_linked": ("linked · anchored ECM", b_base),
        "b_mobile_boundary": ("linked · mobile-boundary control", replace(b_base, boundary_mode="mobile")),
    }

    elastic_tasks = []
    for case_id, (label, cfg) in {**a_cases, **b_cases}.items():
        elastic_tasks.append((case_id, label, asdict(cfg), str(data_dir), reuse_elastic))
    if workers > 1:
        with Pool(processes=workers) as pool:
            elastic_results = pool.map(_elastic_worker, elastic_tasks)
    else:
        elastic_results = [_elastic_worker(task) for task in elastic_tasks]
    case_index = {case_id: {"label": label, **summary}
                  for case_id, label, summary in elastic_results}

    clutch_base = G4V2Config(
        duration=c_duration, sample_interval=c_geometry,
        metric_sample_interval=c_metric, crosslink_probability=0.60,
        crosslink_stiffness=25.0, bead_drag=180.0,
    )
    ensemble = shared_cluster_ensemble(
        clutch_base, trials=ensemble_trials, duration=ensemble_duration, seed0=9000
    )
    clutch_base = with_counter_seed(clutch_base, ensemble["representative_seed"])
    c_spec = _assemble_g4_spec(clutch_base, clutch_base.seed)
    shared = run_clutch(
        clutch_base, spec=c_spec, mode="shared", moving=False, capture_event=True
    )
    shared_event_time = shared.get("event_time")
    independent = run_clutch(
        clutch_base, spec=c_spec, mode="independent", moving=False,
        capture_event=True, event_center_time=shared_event_time,
    )
    independent["event_site"] = shared.get("event_site", -1)
    for case_id, label, result in (
        ("c_independent", "independent Bell clutches", independent),
        ("c_shared", "shared-load adhesion cluster", shared),
    ):
        payload = _case_payload(result, c_spec.positions, case_id, label)
        path = data_dir / f"{case_id}.json"
        _write_json(path, payload)
        case_index[case_id] = {
            "label": label, "file": f"g4-v2-data/{case_id}.json",
            "bytes": path.stat().st_size,
            "siteFailures": int(payload["metrics"][-1]["siteFailures"]),
            "slips": int(payload["metrics"][-1]["slips"]),
        }
        event = _event_payload(result, c_spec.positions, case_id + "_event", label)
        if event is not None:
            event_path = data_dir / f"{case_id}_event.json"
            _write_json(event_path, event)
            case_index[case_id]["eventFile"] = f"g4-v2-data/{case_id}_event.json"
            case_index[case_id]["eventBytes"] = event_path.stat().st_size

    d_base = replace(
        clutch_base, duration=d_duration, sample_interval=d_geometry,
        metric_sample_interval=d_metric, event_half_window=5.0,
    )
    d_tasks = [
        ("d_fixed", "shared-load · fixed cell", asdict(d_base), "shared", False, str(data_dir)),
        ("d_moving", "shared-load · released cell", asdict(d_base), "shared", True, str(data_dir)),
        ("d_mobile_ecm", "released cell · mobile ECM control",
         asdict(replace(d_base, boundary_mode="mobile")), "shared", True, str(data_dir)),
    ]
    if workers > 1:
        with Pool(processes=min(workers, len(d_tasks))) as pool:
            d_results = pool.map(_clutch_worker, d_tasks)
    else:
        d_results = [_clutch_worker(task) for task in d_tasks]
    for case_id, label, summary in d_results:
        case_index[case_id] = {"label": label, **summary}

    lifetimes = np.asarray(ensemble["lifetimes"], dtype=float)
    if len(lifetimes):
        upper = max(float(np.percentile(lifetimes, 99)), 1.0)
        hist, bins = np.histogram(lifetimes, bins=24, range=(0.0, upper))
    else:
        hist, bins = np.zeros(24, dtype=int), np.linspace(0.0, 1.0, 25)
    ensemble_web = {
        "trials": ensemble["trials"], "duration": ensemble["duration"],
        "episodes": int(len(lifetimes)), "medianLifetime": _round(ensemble["median_lifetime"], 4),
        "medianFirstFailure": _round(ensemble["median_first_failure"], 4),
        "representativeSeed": ensemble["representative_seed"],
        "histogram": {"bins": _round(bins, 4), "counts": hist.astype(int).tolist()},
        "protocol": "isolated shared-load site, zero substrate speed; same motor/Bell/rebinding equations",
    }
    manifest = {
        "version": "G4 v2", "quick": bool(quick),
        "units": {"length": "µm", "force": "nN", "time": "s"},
        "domain": a_base.domain_size, "cellRadius": a_base.cell_radius,
        "initial": _round(spec.positions, 4),
        "edges": np.asarray(baseline_network.edges, dtype=int).tolist(),
        "edgeFiber": np.asarray(baseline_network.edge_fiber, dtype=int).tolist(),
        "fibers": [np.asarray(ids, dtype=int).tolist() for ids in baseline_network.fibers],
        "axes": axes,
        "boundaryCases": {"anchored": "a_baseline", "mobile": "a_mobile_boundary"},
        "stageCases": {
            "a": list(a_cases) + ["a_mobile_boundary"],
            "b": list(b_cases),
            "c": ["c_independent", "c_shared"],
            "d": ["d_fixed", "d_moving", "d_mobile_ecm"],
        },
        "ensemble": ensemble_web, "cases": case_index,
        "clutchDefaults": {
            "N": clutch_base.n_clutches_per_site,
            "kClutch": clutch_base.clutch_stiffness,
            "kOn": clutch_base.clutch_on_rate,
            "kOff0": clutch_base.clutch_off_rate0,
            "bellForce": clutch_base.bell_force,
            "v0": clutch_base.unloaded_actin_speed,
            "stallPerSite": clutch_base.motor_stall_per_site,
        },
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        "window.G4V2_MANIFEST=" + json.dumps(manifest, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )
    compact_summary = {
        "version": manifest["version"], "quick": manifest["quick"],
        "durations": {
            "elastic": elastic_duration, "clutch": c_duration, "movement": d_duration,
        },
        "ensemble": ensemble_web,
        "cases": {key: {k: v for k, v in value.items() if k not in ("file", "eventFile")}
                  for key, value in case_index.items()},
        "interpretation": {
            "visual_scaling": "geometry remains true-scale; displacement arrows may be 1x/10x/50x",
            "shared_load": "coarse-grained equal-load-sharing hypothesis, compared with independent baseline",
            "plasticity": "absent; crosslinks remain permanent elastic links",
        },
    }
    _write_json(summary_path, compact_summary)
    print(f"wrote {len(case_index)} lazy cases; manifest {manifest_path.stat().st_size/1e3:.1f} kB")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--quick", action="store_true", help="short smoke-test datasets")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--reuse-elastic", action="store_true",
        help="reuse existing A/B JSON while rebuilding clutch and movement cases",
    )
    args = parser.parse_args()
    build(
        data_dir=args.data_dir, manifest_path=args.manifest,
        summary_path=args.summary, quick=args.quick,
        workers=max(1, args.workers), reuse_elastic=args.reuse_elastic,
    )


if __name__ == "__main__":
    main()
