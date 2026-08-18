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
    cell_contact_forces,
    conservation_errors,
    project_clutch_forces,
    reaction_force_and_torque,
    update_spatial_clutches,
)
from .protrusions import (
    geometry_scores,
    protrusion_tips,
    step_intrinsic_polarity,
    step_protrusion_lengths,
    update_traction_scores,
)
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
            "g3c": cfg.duration_g3c, "g3s": cfg.duration_g3s}[stage]


def _snapshot(time, positions, cell, clutches, protrusions, material, motor,
              force, torque, foi):
    bound = clutches.bound
    _, tips = protrusion_tips(protrusions, cell)
    return G3Snapshot(
        time=float(time),
        positions=positions.copy(),
        cell_center=cell.center.copy(),
        cell_angle=float(cell.body_angle),
        bound_points=material[bound].copy(),
        bound_fibre_ids=clutches.fiber_id[bound].copy(),
        bound_sector_ids=clutches.sector_id[bound].copy(),
        motor_points=motor[bound].copy(),
        clutch_forces=clutches.force_vector[bound].copy(),
        active_sectors=np.flatnonzero(protrusions.active).copy(),
        active_sector_age=protrusions.active_age.copy(),
        protrusion_tips=tips.copy(),
        protrusion_lengths=protrusions.length.copy(),
        polarity_activity=protrusions.activity.copy(),
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
    if stage not in ("g3a", "g3b", "g3c", "g3s"):
        raise ValueError("stage must be g3a, g3b, g3c, or g3s")
    cfg = cfg or G3Config()
    cfg.validate()
    rng = np.random.default_rng(seed)
    fixture = fixture or build_fixture(fixture_name, cfg, rng)
    positions = fixture.initial_positions.copy()
    cell = RigidCellState.at_origin(cfg.cell_radius)
    n_active = 1 if stage == "g3a" else cfg.n_active_protrusions
    prescribed = [0] if stage == "g3a" else None
    if stage == "g3s":
        n_active = min(cfg.spheroid_surface_patches, cfg.n_sectors)
        prescribed = np.linspace(0, cfg.n_sectors - 1, n_active, dtype=int)
    protrusions = ProtrusionState.initialize(cfg.n_sectors, n_active, rng, prescribed)
    clutches = ClutchState.empty(cfg.n_clutches)

    duration = _duration_for_stage(stage, cfg) if duration is None else float(duration)
    if duration <= 0.0:
        raise ValueError("duration must be positive")
    n_steps = int(round(duration / cfg.dt))
    geometry_every = max(1, int(round(cfg.geometry_update_interval / cfg.dt)))
    polarity_every = max(1, int(round(cfg.polarity_update_interval / cfg.dt)))
    metrics_every = max(1, int(round(cfg.metrics_interval / cfg.dt)))
    frame_every = max(1, int(round(cfg.frame_interval / cfg.dt)))

    trace = {name: [] for name in (
        "time", "cell_x", "cell_y", "cell_angle", "bound_count", "foi", "cell_force_x",
        "cell_force_y", "cell_torque", "force_error", "torque_error", "elastic_energy",
        "contact_count", "contact_penetration", "contact_energy", "protrusion_x",
        "protrusion_y", "polarity_magnitude", "max_protrusion_length",
        "attached_protrusions",
    )}
    snapshots = []
    status = "complete"
    material = np.full((cfg.n_clutches, 2), np.nan)
    motor = np.zeros((cfg.n_clutches, 2))
    cell_force = np.zeros(2)
    cell_torque = 0.0
    force_error = 0.0
    torque_error = 0.0
    contact_count = 0
    contact_penetration = 0.0
    contact_energy = 0.0

    geometry_scores(protrusions, cell, positions, fixture, cfg)

    for step in range(n_steps + 1):
        time = step * cfg.dt
        diagnostic_step = step % metrics_every == 0 or step == n_steps
        if step % geometry_every == 0:
            geometry_scores(protrusions, cell, positions, fixture, cfg)

        if stage in ("g3b", "g3c"):
            update_traction_scores(protrusions, clutches, cfg)
            if step % polarity_every == 0:
                step_intrinsic_polarity(
                    protrusions,
                    cfg,
                    rng,
                    dt=polarity_every * cfg.dt,
                    feedback_enabled=feedback_enabled,
                )

        step_protrusion_lengths(protrusions, cfg, clutches=clutches)

        material, _, motor, _ = update_spatial_clutches(
            clutches, protrusions, cell, positions, fixture, cfg, time, rng)
        bead_traction = project_clutch_forces(clutches, material, positions, fixture, cfg)
        (
            bead_contact, contact_cell_force, contact_cell_torque,
            contact_penetration, contact_count, contact_energy,
        ) = cell_contact_forces(positions, cell, cfg)
        # Reaction drives G3C every step. Fixed-cell stages need it only when a diagnostic is
        # recorded. Conservation is an identity checked by unit tests and sampled here at the
        # registered metrics cadence; recomputing it on every unrecorded step changes no state.
        if stage == "g3c" or diagnostic_step:
            clutch_cell_force, clutch_cell_torque = reaction_force_and_torque(
                clutches, motor, cell)
            cell_force = clutch_cell_force + contact_cell_force
            cell_torque = clutch_cell_torque + contact_cell_torque
        if diagnostic_step:
            force_error, torque_error = conservation_errors(
                clutches, material, motor, bead_traction, positions, cell,
                contact_bead_forces=bead_contact,
                contact_cell_force=contact_cell_force,
                contact_cell_torque=contact_cell_torque,
            )

        if diagnostic_step:
            foi = fibre_orientation_index(positions, fixture, cell.center)
            active = np.flatnonzero(protrusions.active)
            angles = cell.body_angle + protrusions.sector_angles[active]
            weights = 1.0 + protrusions.traction_score[active]
            pvec = (weights[:, None] * np.column_stack((np.cos(angles), np.sin(angles)))).sum(axis=0)
            activity_angles = cell.body_angle + protrusions.sector_angles
            activity_vector = np.sum(
                protrusions.activity[:, None]
                * np.column_stack((np.cos(activity_angles), np.sin(activity_angles))),
                axis=0,
            )
            attached = np.unique(clutches.sector_id[clutches.bound & (clutches.sector_id >= 0)])
            values = (
                time, cell.center[0], cell.center[1], cell.body_angle, int(clutches.bound.sum()),
                foi, cell_force[0], cell_force[1], cell_torque, force_error, torque_error,
                elastic_energy(positions, fixture, cfg), contact_count,
                contact_penetration, contact_energy, activity_vector[0], activity_vector[1],
                np.linalg.norm(activity_vector), protrusions.length.max(initial=0.0),
                attached.size,
            )
            for key, value in zip(trace, values):
                trace[key].append(value)

        if step % frame_every == 0 or step == n_steps:
            snapshots.append(_snapshot(
                time, positions, cell, clutches, protrusions, material, motor,
                cell_force, cell_torque, fibre_orientation_index(positions, fixture, cell.center)))

        if step == n_steps:
            break

        positions = advance_ecm(positions, bead_traction + bead_contact, fixture, cfg, rng)
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
        "max_contact_count": int(traces["contact_count"].max(initial=0)),
        "max_contact_penetration_m": float(
            traces["contact_penetration"].max(initial=0.0)
        ),
        "max_contact_energy_J": float(traces["contact_energy"].max(initial=0.0)),
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
    contact_force, _, _, _, _, peak_contact_energy = cell_contact_forces(
        positions, loaded.cell, cfg)
    peak_elastic_energy = elastic_energy(positions, loaded.fixture, cfg)
    peak_energy = peak_elastic_energy + peak_contact_energy
    signal = abs(foi_pre - foi_initial)
    rng = np.random.default_rng(seed + 1_000_000)
    max_steps = int(round(recovery_duration / cfg.dt))
    recovery_time = 0.0
    resolved = peak_energy <= np.finfo(float).tiny
    for step in range(max_steps):
        contact_force, _, _, _, _, _ = cell_contact_forces(positions, loaded.cell, cfg)
        positions = advance_ecm(positions, contact_force, loaded.fixture, cfg, rng)
        recovery_time = (step + 1) * cfg.dt
        _, _, _, _, _, contact_energy = cell_contact_forces(positions, loaded.cell, cfg)
        energy = elastic_energy(positions, loaded.fixture, cfg) + contact_energy
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
        "peak_elastic_energy_J": peak_elastic_energy,
        "peak_contact_energy_J": peak_contact_energy,
        "peak_recoverable_energy_J": peak_energy,
        "final_elastic_energy_J": elastic_energy(positions, loaded.fixture, cfg),
        "final_contact_energy_J": cell_contact_forces(
            positions, loaded.cell, cfg
        )[-1],
        "recovery_time_s": recovery_time,
        "recovery_resolved": resolved,
        "status": status,
    }
