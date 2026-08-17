"""Gated G3A/G3B/G3C simulation loop."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .analysis import displacement_statistics, elastic_energy, fibre_orientation_index, nam_plasticity_index
from .config import G3Config
from .fixtures import G3Fixture, build_fixture
from .mechanics import (
    advance_cell,
    advance_ecm,
    conservation_errors,
    project_clutch_forces,
    reaction_force_and_torque,
    update_spatial_clutches,
)
from .protrusions import geometry_scores, step_protrusion_turnover, update_traction_scores
from .state import ClutchState, G3Snapshot, ProtrusionState, RigidCellState


@dataclass
class G3RunResult:
    stage: str
    fixture_name: str
    seed: int
    status: str
    config: G3Config
    fixture: G3Fixture
    initial_positions: np.ndarray
    final_positions: np.ndarray
    cell: RigidCellState
    clutches: ClutchState
    protrusions: ProtrusionState
    traces: dict[str, np.ndarray]
    snapshots: list[G3Snapshot]
    summary: dict


def _duration_for_stage(stage: str, cfg: G3Config) -> float:
    return {"g3a": cfg.duration_g3a, "g3b": cfg.duration_g3b,
            "g3c": cfg.duration_g3c}[stage]


def _snapshot(time, positions, cell, clutches, protrusions, material, motor,
              force, torque, foi):
    bound = clutches.bound
    return G3Snapshot(
        time=float(time),
        positions=positions.copy(),
        cell_center=cell.center.copy(),
        cell_angle=float(cell.body_angle),
        bound_points=material[bound].copy(),
        motor_points=motor[bound].copy(),
        clutch_forces=clutches.force_vector[bound].copy(),
        active_sectors=np.flatnonzero(protrusions.active).copy(),
        geometry_scores=protrusions.geometry_score.copy(),
        traction_scores=protrusions.traction_score.copy(),
        foi=float(foi),
        cell_force=force.copy(),
        cell_torque=float(torque),
    )


def run_g3(
    stage: str,
    fixture_name: str,
    cfg: G3Config | None = None,
    seed: int = 0,
    duration: float | None = None,
    feedback_enabled: bool = True,
    fixture: G3Fixture | None = None,
) -> G3RunResult:
    """Run one G3 trajectory. G3A/G3B fix the cell; only G3C releases it."""
    stage = stage.lower()
    if stage not in ("g3a", "g3b", "g3c"):
        raise ValueError("stage must be g3a, g3b, or g3c")
    cfg = cfg or G3Config()
    cfg.validate()
    rng = np.random.default_rng(seed)
    fixture = fixture or build_fixture(fixture_name, cfg, rng)
    positions = fixture.initial_positions.copy()
    cell = RigidCellState.at_origin(cfg.cell_radius)
    n_active = 1 if stage == "g3a" else cfg.n_active_protrusions
    prescribed = [0] if stage == "g3a" else None
    protrusions = ProtrusionState.initialize(cfg.n_sectors, n_active, rng, prescribed)
    clutches = ClutchState.empty(cfg.n_clutches)

    duration = _duration_for_stage(stage, cfg) if duration is None else float(duration)
    if duration <= 0.0:
        raise ValueError("duration must be positive")
    n_steps = int(round(duration / cfg.dt))
    geometry_every = max(1, int(round(cfg.geometry_update_interval / cfg.dt)))
    metrics_every = max(1, int(round(cfg.metrics_interval / cfg.dt)))
    frame_every = max(1, int(round(cfg.frame_interval / cfg.dt)))

    trace = {name: [] for name in (
        "time", "cell_x", "cell_y", "cell_angle", "bound_count", "foi", "cell_force_x",
        "cell_force_y", "cell_torque", "force_error", "torque_error", "elastic_energy",
        "protrusion_x", "protrusion_y",
    )}
    snapshots = []
    status = "complete"
    material = np.full((cfg.n_clutches, 2), np.nan)
    motor = np.zeros((cfg.n_clutches, 2))
    cell_force = np.zeros(2)
    cell_torque = 0.0
    force_error = 0.0
    torque_error = 0.0

    geometry_scores(protrusions, cell, positions, fixture, cfg)

    for step in range(n_steps + 1):
        time = step * cfg.dt
        if step % geometry_every == 0:
            geometry_scores(protrusions, cell, positions, fixture, cfg)

        if stage in ("g3b", "g3c"):
            update_traction_scores(protrusions, clutches, cfg)
            step_protrusion_turnover(protrusions, cfg, rng, feedback_enabled)

        material, _, motor, _ = update_spatial_clutches(
            clutches, protrusions, cell, positions, fixture, cfg, time, rng)
        bead_traction = project_clutch_forces(clutches, material, positions, fixture, cfg)
        cell_force, cell_torque = reaction_force_and_torque(clutches, motor, cell)
        force_error, torque_error = conservation_errors(
            clutches, material, motor, bead_traction, positions, cell)

        if step % metrics_every == 0 or step == n_steps:
            foi = fibre_orientation_index(positions, fixture, cell.center)
            active = np.flatnonzero(protrusions.active)
            angles = cell.body_angle + protrusions.sector_angles[active]
            weights = 1.0 + protrusions.traction_score[active]
            pvec = (weights[:, None] * np.column_stack((np.cos(angles), np.sin(angles)))).sum(axis=0)
            values = (
                time, cell.center[0], cell.center[1], cell.body_angle, int(clutches.bound.sum()),
                foi, cell_force[0], cell_force[1], cell_torque, force_error, torque_error,
                elastic_energy(positions, fixture, cfg), pvec[0], pvec[1],
            )
            for key, value in zip(trace, values):
                trace[key].append(value)

        if step % frame_every == 0 or step == n_steps:
            snapshots.append(_snapshot(
                time, positions, cell, clutches, protrusions, material, motor,
                cell_force, cell_torque, fibre_orientation_index(positions, fixture, cell.center)))

        if step == n_steps:
            break

        positions = advance_ecm(positions, bead_traction, fixture, cfg, rng)
        if stage == "g3c":
            advance_cell(cell, cell_force, cell_torque, cfg)

        if positions.size:
            distance_to_cell = np.linalg.norm(positions - cell.center[None, :], axis=1)
            if np.any(distance_to_cell < cell.radius - cfg.overlap_tolerance):
                status = "invalid_geometry_overlap"
                break

    traces = {key: np.asarray(value) for key, value in trace.items()}
    displacement = displacement_statistics(fixture.initial_positions, positions)
    path = np.column_stack((traces["cell_x"], traces["cell_y"]))
    path_length = float(np.linalg.norm(np.diff(path, axis=0), axis=1).sum()) if path.shape[0] > 1 else 0.0
    net_displacement = float(np.linalg.norm(cell.center))
    summary = {
        "status": status,
        "stage": stage,
        "fixture": fixture.name,
        "seed": int(seed),
        "simulated_time_s": float(traces["time"][-1]) if traces["time"].size else 0.0,
        "n_fibres": int(fixture.network.n_fibers),
        "n_beads": int(fixture.network.n_beads),
        "final_bound_clutches": int(clutches.bound.sum()),
        "max_bound_clutches": int(traces["bound_count"].max(initial=0)),
        "initial_foi": float(traces["foi"][0]) if traces["foi"].size else float("nan"),
        "final_foi": float(traces["foi"][-1]) if traces["foi"].size else float("nan"),
        "max_force_error": float(traces["force_error"].max(initial=0.0)),
        "max_torque_error": float(traces["torque_error"].max(initial=0.0)),
        "cell_net_displacement_m": net_displacement,
        "cell_path_length_m": path_length,
        "cell_final_angle_rad": float(cell.body_angle),
        "bead_displacement": displacement,
        "interpretation": "Personal simulation / mechanism validation; not yet a realistic 3D tumor-migration prediction.",
    }
    return G3RunResult(stage, fixture.name, seed, status, cfg, fixture,
                       fixture.initial_positions.copy(), positions.copy(), cell, clutches,
                       protrusions, traces, snapshots, summary)


def run_load_unload(
    cfg: G3Config | None = None,
    fixture_name: str = "single_fibre",
    seed: int = 0,
    pull_duration: float = 15.0,
    recovery_duration: float = 600.0,
):
    """G3 elastic load--unload protocol for FOI recovery and Nam's kappa diagnostic."""
    cfg = cfg or G3Config()
    # Fixtures are constructed exactly at their elastic rest lengths, so the planned 60 s
    # no-load equilibration is an analytical no-op (zero force at every step).
    loaded = run_g3("g3a", fixture_name, cfg, seed, duration=pull_duration)
    positions = loaded.final_positions.copy()
    foi_initial = fibre_orientation_index(
        loaded.initial_positions, loaded.fixture, loaded.cell.center)
    foi_pre = fibre_orientation_index(positions, loaded.fixture, loaded.cell.center)
    peak_energy = elastic_energy(positions, loaded.fixture, cfg)
    signal = abs(foi_pre - foi_initial)
    rng = np.random.default_rng(seed + 1_000_000)
    max_steps = int(round(recovery_duration / cfg.dt))
    recovery_time = 0.0
    resolved = peak_energy <= np.finfo(float).tiny
    for step in range(max_steps):
        positions = advance_ecm(positions, np.zeros_like(positions), loaded.fixture, cfg, rng)
        recovery_time = (step + 1) * cfg.dt
        energy = elastic_energy(positions, loaded.fixture, cfg)
        if peak_energy > 0.0 and energy <= 0.01 * peak_energy:
            current_foi = fibre_orientation_index(positions, loaded.fixture, loaded.cell.center)
            current_kappa = nam_plasticity_index(foi_pre, current_foi, reference=foi_initial)
            if signal < 1.0e-4 or (np.isfinite(current_kappa) and current_kappa < 0.1):
                resolved = True
                break
    foi_post = fibre_orientation_index(positions, loaded.fixture, loaded.cell.center)
    kappa = nam_plasticity_index(foi_pre, foi_post, reference=foi_initial)
    if signal < 1.0e-4:
        status = "insufficient_foi_signal"
    else:
        status = "complete" if resolved else "unresolved_recovery"
    return {
        "foi_initial_measured": foi_initial,
        "foi_pre_unload": foi_pre,
        "foi_post_recovery": foi_post,
        "foi_pull_signal": signal,
        "kappa": kappa,
        "kappa_reference": "measured initial FOI (equals Nam's random baseline only for an ideal random start)",
        "peak_elastic_energy_J": peak_energy,
        "final_elastic_energy_J": elastic_energy(positions, loaded.fixture, cfg),
        "recovery_time_s": recovery_time,
        "recovery_resolved": resolved,
        "status": status,
    }
