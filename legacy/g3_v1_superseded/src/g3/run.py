"""Command-line entry point for reproducible G3A/G3B/G3C runs."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import yaml

from .config import DEFAULT_CONFIG_PATH, G3Config
from .fixtures import build_fixture
from .simulation import G3RunResult, run_g3, run_load_unload
from .state import ClutchState, G3Snapshot, ProtrusionState, RigidCellState
from .visualization import make_stage_animation, make_summary_figure


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _git_commit():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def save_run(result, output: Path, make_gif=True, make_mp4=False):
    output.mkdir(parents=True, exist_ok=True)
    with open(output / "resolved_config.yaml", "w", encoding="utf-8") as handle:
        yaml.safe_dump({"g3": result.config.to_dict()}, handle, sort_keys=False)
    with open(output / "metrics.json", "w", encoding="utf-8") as handle:
        json.dump(_json_ready(result.summary), handle, indent=2, ensure_ascii=False)
    manifest = {
        "stage": result.stage,
        "fixture": result.fixture_name,
        "seed": result.seed,
        "git_commit": _git_commit(),
        "status": result.status,
        "interpretation": result.summary["interpretation"],
    }
    with open(output / "manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
    np.savez_compressed(output / "traces.npz", **result.traces)
    n_frames = len(result.snapshots)
    n_clutches = result.config.n_clutches
    bound_mask = np.zeros((n_frames, n_clutches), dtype=bool)
    bound_points = np.full((n_frames, n_clutches, 2), np.nan)
    motor_points = np.full((n_frames, n_clutches, 2), np.nan)
    clutch_forces = np.zeros((n_frames, n_clutches, 2))
    bound_fibre_ids = np.full((n_frames, n_clutches), -1, dtype=int)
    bound_sector_ids = np.full((n_frames, n_clutches), -1, dtype=int)
    for frame, snapshot in enumerate(result.snapshots):
        count = snapshot.bound_points.shape[0]
        bound_mask[frame, :count] = True
        bound_points[frame, :count] = snapshot.bound_points
        motor_points[frame, :count] = snapshot.motor_points
        clutch_forces[frame, :count] = snapshot.clutch_forces
        bound_fibre_ids[frame, :count] = snapshot.bound_fibre_ids
        bound_sector_ids[frame, :count] = snapshot.bound_sector_ids
    np.savez_compressed(
        output / "frames.npz",
        time=np.asarray([snapshot.time for snapshot in result.snapshots]),
        positions=np.asarray([snapshot.positions for snapshot in result.snapshots]),
        cell_center=np.asarray([snapshot.cell_center for snapshot in result.snapshots]),
        cell_angle=np.asarray([snapshot.cell_angle for snapshot in result.snapshots]),
        active_sectors=np.asarray([
            np.isin(np.arange(result.config.n_sectors), snapshot.active_sectors)
            for snapshot in result.snapshots
        ]),
        active_sector_age=np.asarray([
            snapshot.active_sector_age for snapshot in result.snapshots
        ]),
        protrusion_tips=np.asarray([snapshot.protrusion_tips for snapshot in result.snapshots]),
        protrusion_lengths=np.asarray([
            snapshot.protrusion_lengths for snapshot in result.snapshots
        ]),
        polarity_activity=np.asarray([
            snapshot.polarity_activity for snapshot in result.snapshots
        ]),
        bound_mask=bound_mask,
        bound_points=bound_points,
        motor_points=motor_points,
        clutch_forces=clutch_forces,
        bound_fibre_ids=bound_fibre_ids,
        bound_sector_ids=bound_sector_ids,
    )
    make_summary_figure(result, output / f"{result.stage}_summary.png")
    if make_gif:
        make_stage_animation(result, output / f"{result.stage}_{result.fixture_name}.gif")
    if make_mp4:
        make_stage_animation(result, output / f"{result.stage}_{result.fixture_name}.mp4")


def load_saved_run(output: str | Path) -> G3RunResult:
    """Reconstruct a renderable run from saved numerical frames without rerunning physics."""
    output = Path(output)
    with open(output / "resolved_config.yaml", encoding="utf-8") as handle:
        config_values = (yaml.safe_load(handle) or {}).get("g3", {})
    # ``to_dict`` records these useful derived values for provenance, but they
    # are properties rather than constructor inputs.
    config_values.pop("rotational_drag", None)
    config_values.pop("beads_per_fibre", None)
    cfg = G3Config.from_dict(config_values)
    with open(output / "manifest.json", encoding="utf-8") as handle:
        manifest = json.load(handle)
    with open(output / "metrics.json", encoding="utf-8") as handle:
        summary = json.load(handle)
    seed = int(manifest["seed"])
    fixture = build_fixture(manifest["fixture"], cfg, np.random.default_rng(seed))
    frames = np.load(output / "frames.npz")
    traces_npz = np.load(output / "traces.npz")
    traces = {name: traces_npz[name] for name in traces_npz.files}
    n_active = 1 if manifest["stage"] == "g3a" else cfg.n_active_protrusions
    protrusions = ProtrusionState.initialize(
        cfg.n_sectors, n_active, np.random.default_rng(seed),
        prescribed=[0] if manifest["stage"] == "g3a" else None,
    )

    snapshots = []
    trace_time = traces["time"]
    for frame, time in enumerate(frames["time"]):
        bound = frames["bound_mask"][frame]
        trace_index = int(np.argmin(np.abs(trace_time - time)))
        active = np.flatnonzero(frames["active_sectors"][frame])
        snapshots.append(G3Snapshot(
            time=float(time),
            positions=frames["positions"][frame].copy(),
            cell_center=frames["cell_center"][frame].copy(),
            cell_angle=float(frames["cell_angle"][frame]),
            bound_points=frames["bound_points"][frame, bound].copy(),
            bound_fibre_ids=frames["bound_fibre_ids"][frame, bound].copy(),
            bound_sector_ids=frames["bound_sector_ids"][frame, bound].copy(),
            motor_points=frames["motor_points"][frame, bound].copy(),
            clutch_forces=frames["clutch_forces"][frame, bound].copy(),
            active_sectors=active,
            active_sector_age=frames["active_sector_age"][frame].copy(),
            protrusion_tips=frames["protrusion_tips"][frame].copy(),
            protrusion_lengths=frames["protrusion_lengths"][frame].copy(),
            polarity_activity=frames["polarity_activity"][frame].copy(),
            geometry_scores=np.zeros(cfg.n_sectors),
            traction_scores=np.zeros(cfg.n_sectors),
            foi=float(traces["foi"][trace_index]),
            cell_force=np.array([
                traces["cell_force_x"][trace_index], traces["cell_force_y"][trace_index]
            ]),
            cell_torque=float(traces["cell_torque"][trace_index]),
        ))

    final = snapshots[-1]
    protrusions.active[:] = False
    protrusions.active[final.active_sectors] = True
    protrusions.active_age[:] = final.active_sector_age
    protrusions.length[:] = final.protrusion_lengths
    protrusions.activity[:] = final.polarity_activity
    cell = RigidCellState.at_origin(cfg.cell_radius)
    cell.center[:] = final.cell_center
    cell.body_angle = final.cell_angle
    return G3RunResult(
        stage=manifest["stage"], fixture_name=manifest["fixture"], seed=seed,
        status=manifest["status"], config=cfg, fixture=fixture,
        initial_positions=fixture.initial_positions.copy(),
        final_positions=final.positions.copy(), cell=cell,
        clutches=ClutchState.empty(cfg.n_clutches), protrusions=protrusions,
        traces=traces, snapshots=snapshots, summary=summary,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description="G3 emergent cell-collagen guidance")
    parser.add_argument("--stage", choices=("g3a", "g3b", "g3c"), required=True)
    parser.add_argument("--fixture", default="single_fibre")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--duration", type=float, default=None, help="override simulated seconds")
    parser.add_argument("--output", required=True)
    parser.add_argument("--feedback-off", action="store_true")
    parser.add_argument("--no-gif", action="store_true")
    parser.add_argument("--mp4", action="store_true")
    parser.add_argument("--load-unload", action="store_true",
                        help="also write the G3A elastic FOI/kappa recovery diagnostic")
    args = parser.parse_args(argv)

    cfg = G3Config.from_yaml(args.config)
    result = run_g3(args.stage, args.fixture, cfg, args.seed, args.duration,
                    feedback_enabled=not args.feedback_off)
    output = Path(args.output)
    save_run(result, output, make_gif=not args.no_gif, make_mp4=args.mp4)
    if args.load_unload:
        diagnostic = run_load_unload(cfg, args.fixture, args.seed)
        with open(output / "load_unload_kappa.json", "w", encoding="utf-8") as handle:
            json.dump(_json_ready(diagnostic), handle, indent=2, ensure_ascii=False)
    print(json.dumps(_json_ready(result.summary), indent=2, ensure_ascii=False))
    return 0 if result.status == "complete" else 2


if __name__ == "__main__":
    raise SystemExit(main())
