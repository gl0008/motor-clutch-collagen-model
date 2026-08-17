"""Preregistered ensemble controls for G3B/G3C directionality claims."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .analysis import nematic_order, polar_order
from .config import G3Config
from .simulation import G3RunResult, run_g3


CALIBRATION_SEEDS = tuple(range(20))
VALIDATION_SEEDS = tuple(range(1000, 1100))


@dataclass(frozen=True)
class EnsembleMetrics:
    n_runs: int
    n_valid: int
    nematic_order: float
    polar_order: float
    positive_fraction: float
    ensemble_polar_ratio: float
    mean_net_displacement_m: float
    mean_path_length_m: float

    def to_dict(self):
        return {
            "n_runs": self.n_runs,
            "n_valid": self.n_valid,
            "nematic_order": self.nematic_order,
            "polar_order": self.polar_order,
            "positive_fraction": self.positive_fraction,
            "ensemble_polar_ratio": self.ensemble_polar_ratio,
            "mean_net_displacement_m": self.mean_net_displacement_m,
            "mean_path_length_m": self.mean_path_length_m,
        }


def _g3b_angles(result: G3RunResult):
    snapshots = result.snapshots[len(result.snapshots) // 2:]
    angles = []
    for snapshot in snapshots:
        angles.extend(snapshot.cell_angle + result.protrusions.sector_angles[snapshot.active_sectors])
    return np.asarray(angles, dtype=float)


def _g3c_angle(result: G3RunResult):
    displacement = result.cell.center
    if np.linalg.norm(displacement) <= np.finfo(float).eps:
        return np.empty(0)
    return np.array([np.arctan2(displacement[1], displacement[0])])


def summarize_ensemble(results: list[G3RunResult], director: float = 0.0) -> EnsembleMetrics:
    valid = [result for result in results if result.status == "complete"]
    angle_sets = [_g3b_angles(result) if result.stage == "g3b" else _g3c_angle(result)
                  for result in valid]
    angles = np.concatenate([value for value in angle_sets if value.size]) if any(
        value.size for value in angle_sets) else np.empty(0)
    vectors = np.column_stack((np.cos(angles), np.sin(angles))) if angles.size else np.empty((0, 2))
    polar_ratio = (float(np.linalg.norm(vectors.mean(axis=0))) if vectors.size else float("nan"))
    positive = float(np.mean(np.cos(angles - director) > 0.0)) if angles.size else float("nan")
    return EnsembleMetrics(
        n_runs=len(results),
        n_valid=len(valid),
        nematic_order=nematic_order(angles, director),
        polar_order=polar_order(angles, director),
        positive_fraction=positive,
        ensemble_polar_ratio=polar_ratio,
        mean_net_displacement_m=float(np.mean([
            result.summary["cell_net_displacement_m"] for result in valid])) if valid else float("nan"),
        mean_path_length_m=float(np.mean([
            result.summary["cell_path_length_m"] for result in valid])) if valid else float("nan"),
    )


def run_ensemble(stage: str, fixture: str, seeds, cfg: G3Config | None = None,
                 duration: float | None = None, feedback_enabled: bool = True):
    cfg = cfg or G3Config()
    return [run_g3(stage, fixture, cfg, int(seed), duration, feedback_enabled)
            for seed in seeds]


def validate_g3b(cfg: G3Config | None = None, seeds=VALIDATION_SEEDS,
                 duration: float = 600.0):
    """Run the aligned/rotated/isotropic/ablation G3B control suite."""
    cfg = cfg or G3Config()
    conditions = {
        "isotropic": ("isotropic_random_8", True, None),
        "aligned": ("aligned_8", True, 0.0),
        "aligned_rotated_30": ("aligned_8_rotated_30", True, np.radians(30.0)),
        "aligned_feedback_off": ("aligned_8", False, 0.0),
        "no_fibre": ("empty", True, None),
    }
    output = {}
    for label, (fixture, feedback, director) in conditions.items():
        results = run_ensemble("g3b", fixture, seeds, cfg, duration, feedback)
        output[label] = summarize_ensemble(results, 0.0 if director is None else director).to_dict()
    aligned = output["aligned"]
    ablated = output["aligned_feedback_off"]
    output["gates"] = {
        "isotropic_polar_ratio_lt_0_1": output["isotropic"]["ensemble_polar_ratio"] < 0.1,
        "aligned_nematic_positive": aligned["nematic_order"] > 0.0,
        "aligned_plus_minus_40_60": 0.4 <= aligned["positive_fraction"] <= 0.6,
        "feedback_reduces_guidance_50pct": (
            aligned["nematic_order"] > 0.0
            and ablated["nematic_order"] <= 0.5 * aligned["nematic_order"]
        ),
    }
    return output


def validate_g3c(cfg: G3Config | None = None, seeds=VALIDATION_SEEDS,
                 duration: float = 600.0):
    """Run G3C hidden-propulsion, aligned-axis, and drag-sensitivity controls."""
    cfg = cfg or G3Config()
    output = {}
    for label, fixture, director in (
        ("no_fibre", "empty", 0.0),
        ("isotropic", "isotropic_random_8", 0.0),
        ("aligned", "aligned_8", 0.0),
        ("aligned_rotated_30", "aligned_8_rotated_30", np.radians(30.0)),
        ("asymmetric_torque", "asymmetric_torque", 0.0),
    ):
        output[label] = summarize_ensemble(
            run_ensemble("g3c", fixture, seeds, cfg, duration), director).to_dict()
    drag = {}
    for value in (0.15, 0.3, 0.6):
        varied = G3Config.from_dict({**{key: value for key, value in cfg.to_dict().items()
                                       if key not in ("rotational_drag", "beads_per_fibre")},
                                     "cell_drag": value})
        drag[str(value)] = summarize_ensemble(
            run_ensemble("g3c", "aligned_8", seeds, varied, duration), 0.0).to_dict()
    output["drag_sweep_N_s_per_m"] = drag
    output["gates"] = {
        "no_hidden_self_propulsion": output["no_fibre"]["mean_net_displacement_m"] == 0.0,
        "aligned_axis_positive": output["aligned"]["nematic_order"] > 0.0,
        "aligned_plus_minus_40_60": 0.4 <= output["aligned"]["positive_fraction"] <= 0.6,
    }
    return output
