"""Generate compact, split JavaScript datasets for the G4 v4E website."""

from __future__ import annotations

import argparse
import base64
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from generations.g4_v4e_contact_guided_single_cell.model import (  # noqa: E402
    G4V4EConfig,
    MODES,
    make_matched_specs,
    run_v4e_condition,
    run_v4e_panel,
)


DATA_DIR = REPO / "docs" / "g4-v4e-data"
MANIFEST = REPO / "docs" / "g4-v4e-manifest.js"
SUMMARY = HERE / "generated_summary.json"


def _plain(value, digits=6):
    if isinstance(value, np.ndarray):
        return np.round(value.astype(float), digits).tolist()
    if isinstance(value, (np.floating, float)):
        return round(float(value), digits)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, complex):
        return [round(value.real, digits), round(value.imag, digits)]
    if isinstance(value, dict):
        return {str(key): _plain(item, digits) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item, digits) for item in value]
    return value


def _pack(states, base, byte=False):
    states = np.asarray(states, dtype=float)
    base = np.asarray(base, dtype=float)
    delta = states - base[None, :, :]
    maximum = float(np.max(np.abs(delta))) if delta.size else 0.0
    scale = max(1e-6, maximum / (120.0 if byte else 32700.0))
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


def _force_components(result):
    parallel, perpendicular, tangent = [], [], []
    fibers = result["network"].fibers
    for positions, points, vectors, fiber_ids in zip(
        result["positions"], result["contact_points"],
        result["contact_vectors"], result["patch_fibers"],
    ):
        frame_parallel, frame_perpendicular, frame_tangent = [], [], []
        for point, vector, fid in zip(points, vectors, fiber_ids):
            ids = fibers[int(fid)]
            a = positions[ids[:-1]]
            b = positions[ids[1:]]
            delta = b - a
            alpha = np.clip(
                np.sum((point - a) * delta, axis=1)
                / np.maximum(np.sum(delta * delta, axis=1), 1e-12), 0.0, 1.0
            )
            nearest = int(np.argmin(np.linalg.norm(a + alpha[:, None] * delta - point, axis=1)))
            local = delta[nearest] / max(float(np.linalg.norm(delta[nearest])), 1e-12)
            component = float(vector @ local) * local
            frame_parallel.append(component)
            frame_perpendicular.append(vector - component)
            frame_tangent.append(local)
        parallel.append(frame_parallel)
        perpendicular.append(frame_perpendicular)
        tangent.append(frame_tangent)
    return _plain(parallel, 5), _plain(perpendicular, 5), _plain(tangent, 5)


def _case_payload(result):
    base = np.asarray(result["spec"].positions)
    parallel, perpendicular, tangents = _force_components(result)
    telemetry = result["telemetry"]
    event = result["event_clip"]
    return {
        "mode": result["mode"],
        "seed": result["seed"],
        "geometry": "cue" if result["mode"] in ("matrix_cue", "combined") else "random",
        "positions": _pack(result["positions"], base),
        "times": _plain(result["overview_times"], 3),
        "centers": _plain(result["centers"], 5),
        "contacts": _plain(result["contact_points"], 5),
        "forces": _plain(result["contact_vectors"], 5),
        "parallel": parallel,
        "perpendicular": perpendicular,
        "tangents": tangents,
        "siteForce": _plain(result["site_force"], 5),
        "boundCount": _plain(result["bound_count"], 0),
        "fronts": _plain(result["fronts"], 5),
        "velocity": _plain(result["cell_velocity"], 6),
        "patchFibers": _plain(result["patch_fibers"], 0),
        "radialOrder": _plain(result["radial_order_by_shell"], 6),
        "deltaRadialOrder": _plain(result["delta_radial_order_by_shell"], 6),
        "crosslinks": _plain(_links(result["network"]), 6),
        "tractFibers": result["tract_fibers"],
        "frontEvents": _plain(result["front_events"], 6),
        "metrics": _plain(result["metrics"], 7),
        "telemetry": {
            "times": _plain(result["telemetry_times"], 3),
            "center": _plain(telemetry["center"], 5),
            "velocity": _plain(telemetry["velocity"], 6),
            "front": _plain(telemetry["front"], 5),
            "siteForce": _plain(telemetry["site_force"], 5),
            "boundCount": _plain(telemetry["bound_count"], 0),
            "cellFiberAngle": _plain(np.degrees(telemetry["cell_fiber_angle"]), 4),
            "Ct": _plain(telemetry["contact_guidance_Ct"], 5),
        },
        "event": {
            "times": _plain(event["times"], 3),
            "center": _plain(event["center"], 5),
            "front": _plain(event["front"], 5),
            "contacts": _plain(event["contact_points"], 5),
            "siteForce": _plain(event["site_force"], 5),
            "boundCount": _plain(event["bound_count"], 0),
        },
    }


def _representative_seed(panel):
    rows = panel["cases"]["combined"]
    keys = ("D_cue", "persistence", "contact_guidance_Ct")
    values = np.asarray([[row[key] for key in keys] for row in rows], dtype=float)
    center = np.median(values, axis=0)
    scale = np.maximum(np.std(values, axis=0), 1e-9)
    distance = np.sum(((values - center) / scale) ** 2, axis=1)
    return int(rows[int(np.argmin(distance))]["seed"])


def _one_seed(seed, cfg, duration):
    """Worker: construct one matched ECM pair, then run all four additions."""

    case_cfg = replace(cfg, seed=int(seed))
    matched = make_matched_specs(case_cfg)
    return {
        mode: dict(
            seed=int(seed),
            **run_v4e_condition(
                case_cfg, mode, int(seed), matched=matched,
                duration=duration, store_positions=False,
            )["metrics"],
        )
        for mode in MODES
    }


def _ci95(values):
    values = np.asarray(values, dtype=float)
    mean = float(np.mean(values))
    if len(values) < 2:
        return mean, mean
    half = 1.96 * float(np.std(values, ddof=1)) / np.sqrt(len(values))
    return mean - half, mean + half


def _aggregate_panel(seed_rows):
    cases = {mode: [] for mode in MODES}
    for row in sorted(seed_rows, key=lambda item: item["d_control"]["seed"]):
        for mode in MODES:
            cases[mode].append(row[mode])
    summary = {}
    for mode, rows in cases.items():
        summary[mode] = {}
        for metric in ("D_cue", "persistence", "contact_guidance_Ct"):
            values = [row[metric] for row in rows]
            summary[mode][metric] = {"mean": float(np.mean(values)), "ci95": _ci95(values)}
        summary[mode]["D_cue_positive_ci"] = _ci95([r["D_cue"] for r in rows])[0] > 0.0
    directions = np.asarray([row["final_direction_angle"] for row in cases["protrusion_guidance"]])
    resultant = abs(np.mean(np.exp(1j * directions)))
    control_ci = summary["d_control"]["D_cue"]["ci95"]
    return {
        "seeds": [row["seed"] for row in cases["d_control"]],
        "cases": cases,
        "summary": summary,
        "gates": {
            "E0_projected_ci_contains_zero": control_ci[0] <= 0.0 <= control_ci[1],
            "E2_no_fixed_world_axis": resultant < 0.5,
            "E3_positive_cue_displacement": summary["combined"]["D_cue_positive_ci"],
            "negative_results_are_retained": True,
        },
    }


def _parallel_panel(cfg, count, duration):
    seeds = list(range(cfg.seed, cfg.seed + count))
    rows = []
    # Each worker keeps all four matched conditions on the same process and
    # geometry. Parallelism changes wall time only, never a model equation.
    with ProcessPoolExecutor(max_workers=min(6, count)) as pool:
        futures = [pool.submit(_one_seed, seed, cfg, duration) for seed in seeds]
        for future in as_completed(futures):
            rows.append(future.result())
    return _aggregate_panel(rows)


def build(*, quick=False):
    cfg = G4V4EConfig()
    if quick:
        cfg = replace(
            cfg,
            migration_duration=600.0,
            migration_sample_interval=30.0,
            telemetry_interval=5.0,
            event_duration=60.0,
            front_maturation=60.0,
            front_loss_time=60.0,
        )
    panel_count = 3 if quick else cfg.panel_seeds
    panel = (
        run_v4e_panel(cfg, panel_count, duration=cfg.migration_duration)
        if quick else _parallel_panel(cfg, panel_count, cfg.migration_duration)
    )
    representative_seed = _representative_seed(panel)
    representative_cfg = replace(cfg, seed=representative_seed)
    matched = make_matched_specs(representative_cfg)
    results = {
        mode: run_v4e_condition(
            representative_cfg, mode, representative_seed,
            matched=matched, duration=cfg.migration_duration,
        )
        for mode in MODES
    }

    random_spec = matched["random_spec"]
    cue_spec = matched["cue_spec"]
    geometry = {
        "fibers": random_spec.fibers,
        "edges": results["d_control"]["network"].edges.astype(int).tolist(),
        "fixed": np.flatnonzero(random_spec.fixed).astype(int).tolist(),
        "randomInitial": _plain(random_spec.positions, 4),
        "cueInitial": _plain(cue_spec.positions, 4),
        "tractFibers": matched["tract_fibers"],
    }
    manifest = {
        "version": "G4 v4E · contact-guided persistent single-cell migration",
        "quick": bool(quick),
        "scientificStatus": "quick pipeline check" if quick else "predeclared full-clock experiment",
        "domain": cfg.domain_size,
        "cellRadius": cfg.cell_radius,
        "units": {"length": "µm", "force": "nN", "time": "s"},
        "actualClock": {
            "duration": cfg.migration_duration,
            "overviewInterval": cfg.migration_sample_interval,
            "telemetryInterval": cfg.telemetry_interval,
            "eventInterval": cfg.event_interval,
        },
        "representativeSeed": representative_seed,
        "selectionRule": "closest-to-multimetric E3 ensemble median; no extreme seed",
        "geometry": geometry,
        "topology": _plain(matched["tract_report"], 7),
        "config": _plain({
            "dt": cfg.dt,
            "cellRadius": cfg.cell_radius,
            "probeReach": cfg.probing_reach,
            "cueHalfWidthDeg": np.degrees(cfg.radial_cue_half_width),
            "tractAngularSdDeg": np.degrees(cfg.radial_tract_angular_sd),
            "frontMaturation": cfg.front_maturation,
            "frontLossTime": cfg.front_loss_time,
            "frontConeHalfWidthDeg": np.degrees(cfg.front_cone_half_width),
            "betaAlignment": cfg.beta_alignment,
            "betaMemory": cfg.beta_memory,
            "nClutchesPerSite": cfg.n_clutches_per_site,
            "clutchStiffness": cfg.clutch_stiffness,
            "onRate": cfg.clutch_on_rate,
            "offRate0": cfg.clutch_off_rate0,
            "bellForce": cfg.bell_force,
            "motorStallPerSite": cfg.motor_stall_per_site,
            "cellDrag": cfg.cell_drag,
        }),
        "panel": _plain(panel, 7),
        "caseFiles": {mode: f"g4-v4e-data/{mode}.js" for mode in MODES},
        "assumptions": [f"V4E-A{index:03d}" for index in range(1, 9)],
        "references": [
            {"label": "Riching et al. 2014", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4255204/"},
            {"label": "Carey et al. 2016", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4980151/"},
            {"label": "Ray et al. 2017", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5394287/"},
        ],
    }

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for mode, result in results.items():
        payload = _case_payload(result)
        (DATA_DIR / f"{mode}.js").write_text(
            "window.G4V4E_CASES=window.G4V4E_CASES||{};"
            f"window.G4V4E_CASES[{json.dumps(mode)}]="
            + json.dumps(_plain(payload), separators=(",", ":")) + ";\n",
            encoding="utf-8",
        )
    MANIFEST.write_text(
        "window.G4V4E_MANIFEST="
        + json.dumps(_plain(manifest), separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )
    summary = {
        "version": manifest["version"],
        "quick": bool(quick),
        "duration": cfg.migration_duration,
        "representative_seed": representative_seed,
        "topology": manifest["topology"],
        "ensemble": panel["summary"],
        "gates": panel["gates"],
        "representative_metrics": {
            mode: result["metrics"] for mode, result in results.items()
        },
    }
    SUMMARY.write_text(json.dumps(_plain(summary), indent=2) + "\n", encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    print(json.dumps(_plain(build(quick=args.quick)), indent=2))


if __name__ == "__main__":
    main()
