"""Generation 2 / V3: two-sided motor clutches and rigid-cell translation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from common.model import (  # noqa: E402
    CollagenConfig,
    Network,
    NetworkSpec,
    _empty_output,
    _finalise_output,
    _radial_metrics,
    contact_patches,
    forces_from_patches,
    make_network_spec,
)


@dataclass(frozen=True)
class MigrationConfig(CollagenConfig):
    """30-minute MDA-MB-231-scale migration experiment."""

    bead_drag: float = 360.0
    dt: float = 0.05
    duration: float = 1800.0
    sample_interval: float = 6.0
    contact_update_interval: float = 0.50
    force_ramp_time: float = 0.0
    n_clutches_per_side: int = 12
    clutch_stiffness: float = 2.0  # nN / um
    clutch_on_rate: float = 0.055  # 1 / s
    clutch_off_rate0: float = 0.018  # 1 / s
    bell_force: float = 1.5  # nN
    unloaded_actin_speed: float = 0.025  # um / s
    motor_stall_force: float = 8.0  # nN per side
    polarity_probability: float = 0.65
    cell_drag: float = 300.0  # nN s / um; calibrated to the speed ensemble below
    max_cell_speed: float = 0.012  # um/s numerical/biological guard


def _random_stream(cfg: MigrationConfig, seed: int):
    nsteps = int(round(cfg.duration / cfg.dt))
    rng = np.random.default_rng(seed)
    shape = (nsteps, 2, cfg.n_clutches_per_side)
    return rng.random(shape), rng.random(shape)


def _material_velocity(network: Network, velocity: np.ndarray, patches) -> float:
    if not patches:
        return 0.0
    values = []
    weights = []
    for patch in patches:
        i, j = network.edges[patch.edge]
        v = (1.0 - patch.alpha) * velocity[i] + patch.alpha * velocity[j]
        values.append(float(v @ patch.normal_in))
        weights.append(patch.weight)
    return float(np.dot(values, weights))


def _clutch_step(
    cfg: MigrationConfig,
    bound: np.ndarray,
    extension: np.ndarray,
    substrate_speed: np.ndarray,
    u_on: np.ndarray,
    u_off: np.ndarray,
):
    force_before = cfg.clutch_stiffness * extension * bound
    traction_before = force_before.sum(axis=1)
    actin_speed = cfg.unloaded_actin_speed * np.maximum(
        0.0, 1.0 - traction_before / cfg.motor_stall_force
    )
    extension[bound] += cfg.dt * np.maximum(
        0.0,
        np.repeat(actin_speed[:, None], cfg.n_clutches_per_side, axis=1)[bound]
        - np.repeat(substrate_speed[:, None], cfg.n_clutches_per_side, axis=1)[bound],
    )
    force = cfg.clutch_stiffness * extension
    off_rate = cfg.clutch_off_rate0 * np.exp(
        np.minimum(force / cfg.bell_force, 12.0)
    )
    breaking = bound & (u_off < 1.0 - np.exp(-off_rate * cfg.dt))
    bound[breaking] = False
    extension[breaking] = 0.0

    side_factor = np.asarray(
        [2.0 * (1.0 - cfg.polarity_probability), 2.0 * cfg.polarity_probability]
    )
    on_probability = 1.0 - np.exp(
        -cfg.clutch_on_rate * side_factor[:, None] * cfg.dt
    )
    binding = (~bound) & (u_on < on_probability)
    bound[binding] = True
    extension[binding] = 0.0

    force = cfg.clutch_stiffness * extension * bound
    traction = force.sum(axis=1)
    actin_speed = cfg.unloaded_actin_speed * np.maximum(
        0.0, 1.0 - traction / cfg.motor_stall_force
    )
    return force, traction, actin_speed


def _record_migration_frame(
    out: dict,
    network: Network,
    center: np.ndarray,
    time: float,
    patches,
    contact_vectors,
    velocity: np.ndarray,
    energy: dict,
    strain: np.ndarray,
    link_force: np.ndarray,
    active: np.ndarray,
    *,
    bound: np.ndarray,
    clutch_force: np.ndarray,
    traction: np.ndarray,
    actin_speed: np.ndarray,
    reaction: np.ndarray,
    cell_velocity: np.ndarray,
) -> None:
    align, shell_disp = _radial_metrics(network, center)
    out["time"].append(float(time))
    out["positions"].append(network.r.copy())
    out["bond_strain"].append(strain.copy())
    out["crosslink_force"].append(link_force.copy())
    out["contact_points"].append([x.point.tolist() for side in patches for x in side])
    out["contact_weights"].append([x.weight for side in patches for x in side])
    out["contact_forces"].append([v.tolist() for side in contact_vectors for v in side])
    out["contact_fibers"].append([x.fiber for side in patches for x in side])
    out["alignment_by_shell"].append(align)
    out["displacement_by_shell"].append(shell_disp)
    disp = network.r - network.r0
    mobile = ~network.fixed
    out["rms_displacement"].append(
        float(np.sqrt(np.mean(np.sum(disp[mobile] ** 2, axis=1))))
    )
    out["stretch_energy"].append(energy["stretch"])
    out["bend_energy"].append(energy["bend"])
    out["crosslink_energy"].append(energy["crosslink"])
    out["drag_power"].append(float(network.cfg.bead_drag * np.sum(velocity**2)))
    out["active_work_rate"].append(float(np.sum(active * velocity)))
    out["min_cell_gap"].append(
        float(np.min(np.linalg.norm(network.r - center, axis=1)) - network.cfg.cell_radius)
    )
    out["cell_center"].append(center.copy())
    out["cell_velocity"].append(cell_velocity.copy())
    out["reaction"].append(reaction.copy())
    out["bound"].append(bound.copy())
    out["clutch_force"].append(clutch_force.copy())
    out["traction"].append(traction.copy())
    out["actin_speed"].append(actin_speed.copy())


def run_migration(
    cfg: MigrationConfig,
    *,
    spec: NetworkSpec,
    moving: bool,
    random_stream,
) -> dict:
    """Run one fixed or moving cell with the same counter-addressed stream."""

    network = Network(spec, cfg)
    center = np.zeros(2)
    bound = np.zeros((2, cfg.n_clutches_per_side), dtype=bool)
    extension = np.zeros_like(bound, dtype=float)
    clutch_force = np.zeros_like(extension)
    traction = np.zeros(2)
    actin_speed = np.full(2, cfg.unloaded_actin_speed)
    substrate_speed = np.zeros(2)
    velocity = np.zeros_like(network.r)
    energy = {"stretch": 0.0, "bend": 0.0, "crosslink": 0.0}
    strain = np.zeros(len(network.edges))
    link_force = np.zeros(len(network.crosslinks))
    reaction = np.zeros(2)
    cell_velocity = np.zeros(2)
    out = _empty_output()
    for name in (
        "cell_center",
        "cell_velocity",
        "reaction",
        "bound",
        "clutch_force",
        "traction",
        "actin_speed",
    ):
        out[name] = []

    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    contact_every = max(1, int(round(cfg.contact_update_interval / cfg.dt)))
    patches = [
        contact_patches(network, center, 180.0),
        contact_patches(network, center, 0.0),
    ]
    u_on, u_off = random_stream

    for step in range(nsteps + 1):
        time = step * cfg.dt
        if step % every == 0:
            vectors = [
                forces_from_patches(network, patches[side], traction[side])[1]
                for side in range(2)
            ]
            active = sum(
                (
                    forces_from_patches(network, patches[side], traction[side])[0]
                    for side in range(2)
                ),
                np.zeros_like(network.r),
            )
            _record_migration_frame(
                out,
                network,
                center,
                time,
                patches,
                vectors,
                velocity,
                energy,
                strain,
                link_force,
                active,
                bound=bound,
                clutch_force=clutch_force,
                traction=traction,
                actin_speed=actin_speed,
                reaction=reaction,
                cell_velocity=cell_velocity,
            )
        if step == nsteps:
            break
        if step and step % contact_every == 0:
            patches = [
                contact_patches(network, center, 180.0),
                contact_patches(network, center, 0.0),
            ]
        clutch_force, traction, actin_speed = _clutch_step(
            cfg,
            bound,
            extension,
            substrate_speed,
            u_on[step],
            u_off[step],
        )
        active = np.zeros_like(network.r)
        for side in range(2):
            nodal, _ = forces_from_patches(network, patches[side], traction[side])
            active += nodal
        reaction = -active.sum(axis=0)
        if moving:
            # V3 intentionally has one translational degree of freedom: the
            # left/right protrusion axis.  Letting incidental y-asymmetry drive
            # a full 2-D trajectory would add a new mechanism before the
            # two-sided baseline is understood.
            cell_velocity = np.array([reaction[0] / cfg.cell_drag, 0.0])
            speed = np.linalg.norm(cell_velocity)
            if speed > cfg.max_cell_speed:
                cell_velocity *= cfg.max_cell_speed / speed
            center += cfg.dt * cell_velocity
        else:
            cell_velocity[:] = 0.0
            center[:] = 0.0
        velocity, energy, strain, link_force, _ = network.advance(
            active, center, cfg.dt, include_crosslinks=True
        )
        substrate_speed = np.asarray(
            [_material_velocity(network, velocity, patches[side]) for side in range(2)]
        )

    out = _finalise_output(out, network, True)
    for name in (
        "cell_center",
        "cell_velocity",
        "reaction",
        "bound",
        "clutch_force",
        "traction",
        "actin_speed",
    ):
        out[name] = np.asarray(out[name])
    out["moving"] = moving
    return out


def run_fixed_moving_pair(cfg: MigrationConfig = MigrationConfig()) -> dict:
    spec = make_network_spec(cfg)
    stream = _random_stream(cfg, cfg.seed + 3001)
    return {
        "fixed": run_migration(
            cfg, spec=spec, moving=False, random_stream=stream
        ),
        "moving": run_migration(
            cfg, spec=spec, moving=True, random_stream=stream
        ),
        "config": asdict(cfg),
        "comparison": "same ECM and same random stream; translation constraint only",
    }


def lumped_speed_trial(cfg: MigrationConfig, seed: int) -> float:
    """Fast mobility calibration using the same clutch law on a rigid ECM."""

    u_on, u_off = _random_stream(cfg, seed)
    bound = np.zeros((2, cfg.n_clutches_per_side), dtype=bool)
    extension = np.zeros_like(bound, dtype=float)
    substrate_speed = np.zeros(2)
    path = 0.0
    nsteps = len(u_on)
    for step in range(nsteps):
        _, traction, _ = _clutch_step(
            cfg, bound, extension, substrate_speed, u_on[step], u_off[step]
        )
        # Right reaction is +x; left reaction is -x.
        velocity_x = np.clip(
            (traction[1] - traction[0]) / cfg.cell_drag,
            -cfg.max_cell_speed,
            cfg.max_cell_speed,
        )
        path += cfg.dt * abs(velocity_x)
    # The experimental comparison reports instantaneous/path speed, not net
    # displacement divided by time (which is reduced by stochastic reversals).
    return float(path / cfg.duration * 60.0)


def run_speed_ensemble(cfg: MigrationConfig, trials: int = 20) -> np.ndarray:
    return np.asarray(
        [lumped_speed_trial(cfg, cfg.seed + 10_000 + i) for i in range(trials)]
    )
