"""Simulation driver and ensemble aggregation."""

from dataclasses import asdict, dataclass
from typing import Dict, List

import numpy as np

from .clutch import MotorClutchCell, MotorClutchParameters
from .network import make_diluted_triangular_network


@dataclass
class SimulationConfig:
    duration: float = 90.0
    dt: float = 0.02
    sample_interval: float = 0.5
    nx: int = 16
    ny: int = 12
    spacing: float = 2.5
    bond_dilution: float = 0.15
    k0: float = 18.0
    kinf: float = 4.5
    bend_stiffness: float = 0.35
    bead_drag: float = 45.0
    boundary_margin: float = 4.0


def run_condition(
    tau: float,
    seed: int,
    config: SimulationConfig,
    clutch_parameters: MotorClutchParameters = None,
) -> Dict[str, np.ndarray]:
    """Run one coupled cell/network realization."""

    if clutch_parameters is None:
        clutch_parameters = MotorClutchParameters()
    rng = np.random.default_rng(seed)
    network = make_diluted_triangular_network(
        nx=config.nx,
        ny=config.ny,
        spacing=config.spacing,
        bond_dilution=config.bond_dilution,
        seed=seed,
        k0=config.k0,
        kinf=config.kinf,
        tau=tau,
        bend_stiffness=config.bend_stiffness,
        drag=config.bead_drag,
    )
    low, high = network.bounds
    cell = MotorClutchCell((low + high) / 2.0, clutch_parameters)
    initial_cell_position = cell.position.copy()

    n_steps = int(round(config.duration / config.dt))
    sample_every = max(1, int(round(config.sample_interval / config.dt)))
    records: List[List[float]] = []
    for step in range(n_steps + 1):
        if step % sample_every == 0:
            displacement = cell.position - initial_cell_position
            records.append(
                [
                    step * config.dt,
                    cell.position[0],
                    cell.position[1],
                    float(np.dot(displacement, displacement)),
                    float(cell.bound_count),
                    cell.mean_bound_force,
                    float(np.mean(cell.completed_lifetimes)) if cell.completed_lifetimes else 0.0,
                    float(np.mean(network.axial_force_magnitudes())),
                ]
            )
        if step == n_steps:
            break
        bead_forces, _ = cell.step(network.positions, config.dt, rng)
        network.advance(bead_forces, config.dt)

        # Keep the point cell inside the calibrated interior; this is a finite
        # simulation box, not a biological confinement force.
        cell.position = np.minimum(
            np.maximum(cell.position, low + config.boundary_margin),
            high - config.boundary_margin,
        )

        if not np.all(np.isfinite(network.positions)) or not np.all(np.isfinite(cell.position)):
            raise FloatingPointError("simulation became non-finite; reduce dt or stiffness")

    data = np.asarray(records, dtype=float)
    return {
        "time": data[:, 0],
        "x": data[:, 1],
        "y": data[:, 2],
        "squared_displacement": data[:, 3],
        "bound_clutches": data[:, 4],
        "mean_clutch_force": data[:, 5],
        "completed_lifetime_mean": data[:, 6],
        "mean_network_axial_force": data[:, 7],
        "final_network_positions": network.positions.copy(),
        "reference_network_positions": network.reference_positions.copy(),
        "edges": network.edges.copy(),
        "config": asdict(config),
    }

