"""Checkpointed preregistered G3B/G3C ensemble-validation campaign."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
import yaml

from .config import DEFAULT_CONFIG_PATH, G3Config
from .fixtures import build_fixture, mirrored_fixture
from .simulation import run_g3
from .validation import CALIBRATION_SEEDS, VALIDATION_SEEDS


@dataclass(frozen=True)
class CampaignCondition:
    name: str
    stage: str
    fixture: str
    director_rad: float | None = None
    feedback_enabled: bool = True
    config_updates: tuple[tuple[str, float], ...] = ()
    mirror_axis: str | None = None

    def to_dict(self):
        value = asdict(self)
        value["config_updates"] = dict(self.config_updates)
        return value


G3B_CONDITIONS = (
    CampaignCondition("balanced", "g3b", "balanced_8"),
    CampaignCondition("isotropic", "g3b", "isotropic_random_8"),
    CampaignCondition("aligned", "g3b", "aligned_8", 0.0),
    CampaignCondition("aligned_rotated_30", "g3b", "aligned_8_rotated_30", np.pi / 6.0),
    CampaignCondition("aligned_feedback_off", "g3b", "aligned_8", 0.0, False),
    CampaignCondition("no_fibre", "g3b", "empty"),
)

G3C_CONDITIONS = (
    CampaignCondition("no_fibre", "g3c", "empty"),
    CampaignCondition("balanced", "g3c", "balanced_8"),
    CampaignCondition("isotropic", "g3c", "isotropic_random_8"),
    CampaignCondition("aligned", "g3c", "aligned_8", 0.0),
    CampaignCondition("aligned_rotated_30", "g3c", "aligned_8_rotated_30", np.pi / 6.0),
    CampaignCondition("asymmetric_torque", "g3c", "asymmetric_torque"),
    CampaignCondition(
        "asymmetric_torque_mirror_x", "g3c", "asymmetric_torque", mirror_axis="x"
    ),
    CampaignCondition("aligned_drag_150", "g3c", "aligned_8", 0.0,
                      config_updates=(("cell_drag", 0.15),)),
    CampaignCondition("aligned_drag_600", "g3c", "aligned_8", 0.0,
                      config_updates=(("cell_drag", 0.6),)),
    CampaignCondition("aligned_rot_drag_0_5", "g3c", "aligned_8", 0.0,
                      config_updates=(("rotational_drag_factor", 0.5),)),
    CampaignCondition("aligned_rot_drag_2", "g3c", "aligned_8", 0.0,
                      config_updates=(("rotational_drag_factor", 2.0),)),
)

ALL_CONDITIONS = G3B_CONDITIONS + G3C_CONDITIONS


def _git_commit():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _config_payload(cfg: G3Config):
    return asdict(cfg)


def campaign_fingerprint(cfg: G3Config, duration: float) -> str:
    payload = {
        "config": _config_payload(cfg),
        "duration_s": float(duration),
        "conditions": [condition.to_dict() for condition in ALL_CONDITIONS],
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _condition_config(base: G3Config, condition: CampaignCondition) -> G3Config:
    return replace(base, **dict(condition.config_updates))


def _axis_metrics(result):
    """Return one independent run's axial and polar summaries from its second half."""
    start = len(result.snapshots) // 2
    snapshots = result.snapshots[start:]
    angles = []
    for snapshot in snapshots:
        angles.extend(
            result.protrusions.sector_angles[snapshot.active_sectors] + snapshot.cell_angle
        )
    angles = np.asarray(angles, dtype=float)
    if angles.size:
        axis_cos2 = float(np.mean(np.cos(2.0 * angles)))
        axis_sin2 = float(np.mean(np.sin(2.0 * angles)))
    else:
        axis_cos2 = axis_sin2 = float("nan")

    if result.stage == "g3b":
        midpoint = len(result.traces["protrusion_x"]) // 2
        polar_x = float(np.mean(result.traces["protrusion_x"][midpoint:]))
        polar_y = float(np.mean(result.traces["protrusion_y"][midpoint:]))
    else:
        polar_x = float(result.cell.center[0])
        polar_y = float(result.cell.center[1])
    magnitude = math.hypot(polar_x, polar_y)
    direction = math.atan2(polar_y, polar_x) if magnitude > np.finfo(float).eps else None
    return axis_cos2, axis_sin2, polar_x, polar_y, direction


def _run_one(condition_data, seed: int, config_data, duration: float):
    condition = CampaignCondition(
        name=condition_data["name"],
        stage=condition_data["stage"],
        fixture=condition_data["fixture"],
        director_rad=condition_data.get("director_rad"),
        feedback_enabled=condition_data.get("feedback_enabled", True),
        config_updates=tuple(sorted(condition_data.get("config_updates", {}).items())),
        mirror_axis=condition_data.get("mirror_axis"),
    )
    base = G3Config.from_dict(config_data)
    cfg = _condition_config(base, condition)
    fixture = None
    if condition.mirror_axis:
        fixture_rng = np.random.default_rng(seed)
        fixture = mirrored_fixture(
            build_fixture(condition.fixture, cfg, fixture_rng), condition.mirror_axis
        )
    started = time.perf_counter()
    result = run_g3(
        condition.stage,
        condition.fixture,
        cfg,
        seed,
        duration,
        feedback_enabled=condition.feedback_enabled,
        fixture=fixture,
    )
    axis_cos2, axis_sin2, polar_x, polar_y, direction = _axis_metrics(result)
    force = np.hypot(result.traces["cell_force_x"], result.traces["cell_force_y"])
    record = {
        "condition": condition.name,
        "stage": condition.stage,
        "fixture": result.fixture_name,
        "seed": int(seed),
        "status": result.status,
        "duration_s": float(result.summary["simulated_time_s"]),
        "director_rad": condition.director_rad,
        "feedback_enabled": condition.feedback_enabled,
        "config_updates": dict(condition.config_updates),
        "axis_cos2": axis_cos2,
        "axis_sin2": axis_sin2,
        "polar_x": polar_x,
        "polar_y": polar_y,
        "direction_rad": direction,
        "cell_x_m": float(result.cell.center[0]),
        "cell_y_m": float(result.cell.center[1]),
        "net_displacement_m": float(result.summary["cell_net_displacement_m"]),
        "path_length_m": float(result.summary["cell_path_length_m"]),
        "final_cell_angle_rad": float(result.cell.body_angle),
        "max_abs_cell_torque_N_m": float(np.max(np.abs(result.traces["cell_torque"]))),
        "max_cell_force_N": float(np.max(force)),
        "max_bound_clutches": int(result.summary["max_bound_clutches"]),
        "initial_foi": float(result.summary["initial_foi"]),
        "final_foi": float(result.summary["final_foi"]),
        "max_force_error": float(result.summary["max_force_error"]),
        "max_torque_error": float(result.summary["max_torque_error"]),
        "wall_time_s": float(time.perf_counter() - started),
    }
    return record


def _atomic_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
    temporary.replace(path)


def _record_path(root: Path, phase: str, condition: CampaignCondition, seed: int):
    return root / phase / "records" / condition.stage / condition.name / f"seed_{seed}.json"


def _checkpoint_is_terminal(path: Path) -> bool:
    """Return true for a completed or scientific-negative run, but retry infrastructure errors."""
    if not path.exists():
        return False
    try:
        status = json.loads(path.read_text(encoding="utf-8")).get("status")
    except (OSError, json.JSONDecodeError):
        return False
    return status not in (None, "worker_error")


def selected_conditions(stages: set[str], names: set[str] | None = None):
    return [condition for condition in ALL_CONDITIONS
            if condition.stage in stages and (not names or condition.name in names)]


def run_campaign(root: Path, phase: str, cfg: G3Config, duration: float, workers: int,
                 conditions):
    seeds = CALIBRATION_SEEDS if phase == "calibration" else VALIDATION_SEEDS
    fingerprint = campaign_fingerprint(cfg, duration)
    if phase == "validation":
        lock_path = root / "config_lock.json"
        if not lock_path.exists():
            raise RuntimeError("validation requires config_lock.json; run the freeze command first")
        locked = json.loads(lock_path.read_text(encoding="utf-8"))
        if locked["campaign_fingerprint"] != fingerprint:
            raise RuntimeError("active config/conditions differ from the frozen calibration lock")

    manifest = {
        "phase": phase,
        "seeds": list(seeds),
        "duration_s": duration,
        "workers": workers,
        "git_commit_at_launch": _git_commit(),
        "campaign_fingerprint": fingerprint,
        "conditions": [condition.to_dict() for condition in conditions],
    }
    _atomic_json(root / phase / "manifest.json", manifest)
    with open(root / phase / "resolved_config.yaml", "w", encoding="utf-8") as handle:
        yaml.safe_dump({"g3": _config_payload(cfg)}, handle, sort_keys=False)

    pending = []
    for condition in conditions:
        for seed in seeds:
            path = _record_path(root, phase, condition, seed)
            if not _checkpoint_is_terminal(path):
                pending.append((condition, int(seed), path))
    print(f"phase={phase} total={len(conditions) * len(seeds)} pending={len(pending)} workers={workers}",
          flush=True)
    if not pending:
        return

    completed = 0
    failed = 0
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(
                _run_one, condition.to_dict(), seed, _config_payload(cfg), duration
            ): (condition, seed, path)
            for condition, seed, path in pending
        }
        for future in as_completed(futures):
            condition, seed, path = futures[future]
            try:
                record = future.result()
            except Exception as error:  # preserve failures as resumable evidence
                failed += 1
                record = {
                    "condition": condition.name,
                    "stage": condition.stage,
                    "seed": seed,
                    "status": "worker_error",
                    "error_type": type(error).__name__,
                    "error": str(error),
                }
            _atomic_json(path, record)
            completed += 1
            elapsed = time.perf_counter() - started
            rate = completed / elapsed if elapsed > 0.0 else 0.0
            remaining = (len(pending) - completed) / rate if rate > 0.0 else float("nan")
            print(
                f"[{completed}/{len(pending)}] {condition.stage}/{condition.name} seed={seed} "
                f"status={record['status']} failed={failed} eta_s={remaining:.0f}",
                flush=True,
            )


def _load_records(root: Path, phase: str, condition: CampaignCondition):
    folder = root / phase / "records" / condition.stage / condition.name
    return [json.loads(path.read_text(encoding="utf-8"))
            for path in sorted(folder.glob("seed_*.json"))] if folder.exists() else []


def _bootstrap_mean_ci(values, rng, repetitions=5000):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return [None, None]
    draws = rng.choice(values, size=(repetitions, values.size), replace=True).mean(axis=1)
    return [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))]


def _axial_difference_degrees(first, second):
    delta = (second - first + 0.5 * np.pi) % np.pi - 0.5 * np.pi
    return abs(float(np.degrees(delta)))


def summarize_condition(records, director):
    valid = [record for record in records if record.get("status") == "complete"]
    rng = np.random.default_rng(20260817)
    c2 = np.asarray([record["axis_cos2"] for record in valid], dtype=float)
    s2 = np.asarray([record["axis_sin2"] for record in valid], dtype=float)
    axis_angle = (0.5 * math.atan2(float(np.nanmean(s2)), float(np.nanmean(c2))) % np.pi
                  if valid else None)
    reference = 0.0 if director is None else director
    alignment = c2 * math.cos(2.0 * reference) + s2 * math.sin(2.0 * reference)
    px = np.asarray([record["polar_x"] for record in valid], dtype=float)
    py = np.asarray([record["polar_y"] for record in valid], dtype=float)
    individual = np.hypot(px, py)
    polar_ratio = (float(math.hypot(float(px.mean()), float(py.mean())) / individual.mean())
                   if individual.size and individual.mean() > np.finfo(float).eps else 0.0)
    direction = np.asarray([
        record["direction_rad"] for record in valid if record["direction_rad"] is not None
    ], dtype=float)
    positive_fraction = (float(np.mean(np.cos(direction - reference) > 0.0))
                         if direction.size else None)
    displacement = np.asarray([record["net_displacement_m"] for record in valid])
    path = np.asarray([record["path_length_m"] for record in valid])
    rotation = np.asarray([record["final_cell_angle_rad"] for record in valid])
    return {
        "n_records": len(records),
        "n_valid": len(valid),
        "n_invalid_overlap": sum(record.get("status") == "invalid_geometry_overlap"
                                 for record in records),
        "n_worker_error": sum(record.get("status") == "worker_error" for record in records),
        "director_rad": director,
        "estimated_axis_rad": axis_angle,
        "estimated_axis_deg": None if axis_angle is None else float(np.degrees(axis_angle)),
        "mean_nematic_alignment": float(np.nanmean(alignment)) if alignment.size else None,
        "nematic_alignment_ci95": _bootstrap_mean_ci(alignment, rng),
        "positive_fraction": positive_fraction,
        "ensemble_to_individual_polar_ratio": polar_ratio,
        "mean_net_displacement_m": float(displacement.mean()) if displacement.size else None,
        "mean_path_length_m": float(path.mean()) if path.size else None,
        "mean_speed_um_per_min": (float(path.mean() / np.mean([
            record["duration_s"] for record in valid]) * 60.0e6) if path.size else None),
        "mean_final_rotation_rad": float(rotation.mean()) if rotation.size else None,
        "mean_abs_final_rotation_rad": float(np.abs(rotation).mean()) if rotation.size else None,
        "mean_cell_x_m": float(np.mean([record["cell_x_m"] for record in valid])) if valid else None,
        "mean_cell_y_m": float(np.mean([record["cell_y_m"] for record in valid])) if valid else None,
        "max_bound_clutches": max((record["max_bound_clutches"] for record in valid), default=None),
        "max_cell_force_N": max((record["max_cell_force_N"] for record in valid), default=None),
        "max_abs_cell_torque_N_m": max(
            (record["max_abs_cell_torque_N_m"] for record in valid), default=None
        ),
        "max_force_error": max((record["max_force_error"] for record in valid), default=None),
        "max_torque_error": max((record["max_torque_error"] for record in valid), default=None),
        "mean_wall_time_s": float(np.mean([record["wall_time_s"] for record in valid]))
        if valid else None,
    }


def summarize_campaign(root: Path, phase: str, conditions):
    summaries = {}
    for condition in conditions:
        records = _load_records(root, phase, condition)
        summaries[f"{condition.stage}/{condition.name}"] = summarize_condition(
            records, condition.director_rad
        )
    _atomic_json(root / phase / "summary.json", summaries)
    rows = []
    for name, values in summaries.items():
        row = {"condition": name}
        row.update({key: value for key, value in values.items() if not isinstance(value, list)})
        rows.append(row)
    if rows:
        with open(root / phase / "summary.csv", "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    gates = evaluate_gates(summaries)
    _atomic_json(root / phase / "gates.json", gates)
    print(json.dumps(summaries, indent=2), flush=True)
    print(json.dumps({"gates": gates}, indent=2), flush=True)
    return summaries


def evaluate_gates(summaries):
    """Evaluate only gates whose required conditions are present; never tune from this output."""
    gates = {}

    def available(*names):
        return all(name in summaries and summaries[name]["n_records"] > 0 for name in names)

    if available("g3b/balanced"):
        gates["g3b_balanced_no_fixed_polarity"] = (
            summaries["g3b/balanced"]["ensemble_to_individual_polar_ratio"] < 0.1
        )
    if available("g3b/isotropic"):
        gates["g3b_isotropic_polar_ratio_lt_0_1"] = (
            summaries["g3b/isotropic"]["ensemble_to_individual_polar_ratio"] < 0.1
        )
    if available("g3b/aligned"):
        aligned = summaries["g3b/aligned"]
        gates["g3b_aligned_guidance_ci_above_zero"] = (
            aligned["nematic_alignment_ci95"][0] is not None
            and aligned["nematic_alignment_ci95"][0] > 0.0
        )
        gates["g3b_aligned_plus_minus_40_60"] = (
            aligned["positive_fraction"] is not None
            and 0.4 <= aligned["positive_fraction"] <= 0.6
        )
    if available("g3b/aligned", "g3b/aligned_rotated_30"):
        base = summaries["g3b/aligned"]["estimated_axis_rad"]
        rotated = summaries["g3b/aligned_rotated_30"]["estimated_axis_rad"]
        error = _axial_difference_degrees(base + np.pi / 6.0, rotated)
        gates["g3b_rotational_covariance_error_deg"] = error
        gates["g3b_rotational_covariance_within_5deg"] = error <= 5.0
    if available("g3b/aligned", "g3b/aligned_feedback_off"):
        aligned = summaries["g3b/aligned"]["mean_nematic_alignment"]
        ablated = summaries["g3b/aligned_feedback_off"]["mean_nematic_alignment"]
        gates["g3b_feedback_reduction_fraction"] = (
            None if aligned is None or aligned <= 0.0 else 1.0 - ablated / aligned
        )
        gates["g3b_feedback_reduces_guidance_50pct"] = (
            aligned is not None and aligned > 0.0 and ablated <= 0.5 * aligned
        )
    if available("g3b/no_fibre"):
        no_fibre = summaries["g3b/no_fibre"]
        gates["g3b_no_fibre_has_zero_traction"] = (
            no_fibre["max_bound_clutches"] == 0 and no_fibre["max_cell_force_N"] == 0.0
        )

    if available("g3c/no_fibre"):
        no_fibre = summaries["g3c/no_fibre"]
        gates["g3c_no_hidden_self_propulsion"] = (
            no_fibre["mean_net_displacement_m"] == 0.0
            and no_fibre["mean_abs_final_rotation_rad"] == 0.0
        )
    if available("g3c/isotropic"):
        gates["g3c_isotropic_polar_ratio_lt_0_1"] = (
            summaries["g3c/isotropic"]["ensemble_to_individual_polar_ratio"] < 0.1
        )
    if available("g3c/aligned"):
        aligned = summaries["g3c/aligned"]
        gates["g3c_aligned_guidance_ci_above_zero"] = (
            aligned["nematic_alignment_ci95"][0] is not None
            and aligned["nematic_alignment_ci95"][0] > 0.0
        )
        gates["g3c_aligned_plus_minus_40_60"] = (
            aligned["positive_fraction"] is not None
            and 0.4 <= aligned["positive_fraction"] <= 0.6
        )
    if available("g3c/aligned", "g3c/aligned_rotated_30"):
        base = summaries["g3c/aligned"]["estimated_axis_rad"]
        rotated = summaries["g3c/aligned_rotated_30"]["estimated_axis_rad"]
        error = _axial_difference_degrees(base + np.pi / 6.0, rotated)
        gates["g3c_rotational_covariance_error_deg"] = error
        gates["g3c_rotational_covariance_within_5deg"] = error <= 5.0
    if available("g3c/asymmetric_torque"):
        asymmetric = summaries["g3c/asymmetric_torque"]
        gates["g3c_asymmetric_torque_nonzero_rotation"] = (
            asymmetric["mean_abs_final_rotation_rad"] is not None
            and asymmetric["mean_abs_final_rotation_rad"] > 1.0e-12
            and asymmetric["max_abs_cell_torque_N_m"] > 0.0
        )
    if available("g3c/asymmetric_torque", "g3c/asymmetric_torque_mirror_x"):
        base = summaries["g3c/asymmetric_torque"]
        mirror = summaries["g3c/asymmetric_torque_mirror_x"]
        gates["g3c_mirror_reverses_mean_rotation"] = (
            base["mean_final_rotation_rad"] * mirror["mean_final_rotation_rad"] < 0.0
        )
        gates["g3c_mirror_reverses_mean_lateral_trajectory"] = (
            base["mean_cell_y_m"] * mirror["mean_cell_y_m"] < 0.0
        )
    if available("g3c/aligned_drag_150", "g3c/aligned", "g3c/aligned_drag_600"):
        fast = summaries["g3c/aligned_drag_150"]["mean_speed_um_per_min"]
        baseline = summaries["g3c/aligned"]["mean_speed_um_per_min"]
        slow = summaries["g3c/aligned_drag_600"]["mean_speed_um_per_min"]
        gates["g3c_drag_speed_monotonic"] = fast > baseline > slow

    for name, summary in summaries.items():
        gates[f"{name}_all_runs_valid"] = (
            summary["n_records"] == summary["n_valid"]
            and summary["n_invalid_overlap"] == 0
            and summary["n_worker_error"] == 0
        )
    return gates


def freeze_campaign(root: Path, cfg: G3Config, duration: float):
    calibration = root / "calibration" / "summary.json"
    if not calibration.exists():
        raise RuntimeError("summarize calibration before freezing the campaign")
    lock = {
        "campaign_fingerprint": campaign_fingerprint(cfg, duration),
        "duration_s": duration,
        "config": _config_payload(cfg),
        "conditions": [condition.to_dict() for condition in ALL_CONDITIONS],
        "calibration_summary_sha256": hashlib.sha256(calibration.read_bytes()).hexdigest(),
        "git_commit_at_freeze": _git_commit(),
        "rule": "Validation seeds 1000-1099 must not be used for parameter tuning.",
    }
    _atomic_json(root / "config_lock.json", lock)
    print(json.dumps(lock, indent=2), flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description="G3 checkpointed ensemble campaign")
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("run", "summarize"):
        sub = subparsers.add_parser(command)
        sub.add_argument("--phase", choices=("calibration", "validation"), required=True)
        sub.add_argument("--stage", choices=("g3b", "g3c", "both"), default="both")
        sub.add_argument("--condition", action="append", default=[])
        sub.add_argument("--output", type=Path, required=True)
        sub.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH)
        sub.add_argument("--duration", type=float, default=600.0)
        if command == "run":
            sub.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    freeze = subparsers.add_parser("freeze")
    freeze.add_argument("--output", type=Path, required=True)
    freeze.add_argument("--config", type=Path, default=DEFAULT_CONFIG_PATH)
    freeze.add_argument("--duration", type=float, default=600.0)
    args = parser.parse_args(argv)

    cfg = G3Config.from_yaml(args.config)
    if args.command == "freeze":
        freeze_campaign(args.output, cfg, args.duration)
        return 0
    stages = {"g3b", "g3c"} if args.stage == "both" else {args.stage}
    conditions = selected_conditions(stages, set(args.condition) if args.condition else None)
    if not conditions:
        raise ValueError("no campaign conditions selected")
    if args.command == "run":
        run_campaign(args.output, args.phase, cfg, args.duration, args.workers, conditions)
    else:
        summarize_campaign(args.output, args.phase, conditions)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
