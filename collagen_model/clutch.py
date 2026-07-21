"""Coarse stochastic motor–clutch coupling between a cell and network beads."""

from dataclasses import dataclass
from math import exp, pi
from typing import List, Tuple

import numpy as np


@dataclass
class MotorClutchParameters:
    n_modules: int = 8
    clutches_per_module: int = 5
    cell_radius: float = 5.5
    capture_radius: float = 3.2
    clutch_stiffness: float = 2.0
    on_rate: float = 0.7
    unloaded_off_rate: float = 0.12
    rupture_force: float = 2.0
    unloaded_actin_speed: float = 0.08
    motors_per_module: int = 6
    motor_stall_force: float = 2.0
    cell_drag: float = 35.0


class MotorClutchCell:
    def __init__(self, position: np.ndarray, parameters: MotorClutchParameters) -> None:
        self.position = np.asarray(position, dtype=float).copy()
        self.parameters = parameters
        angles = np.linspace(0.0, 2.0 * pi, parameters.n_modules, endpoint=False)
        self.directions = np.column_stack((np.cos(angles), np.sin(angles)))
        shape = (parameters.n_modules, parameters.clutches_per_module)
        self.bead_index = np.full(shape, -1, dtype=int)
        self.rest_offset = np.zeros(shape + (2,), dtype=float)
        self.bind_time = np.zeros(shape, dtype=float)
        self.force_magnitude = np.zeros(shape, dtype=float)
        self.actin_slip = np.zeros(parameters.n_modules, dtype=float)
        self.completed_lifetimes: List[float] = []
        self.time = 0.0

    @property
    def bound_count(self) -> int:
        return int(np.count_nonzero(self.bead_index >= 0))

    @property
    def mean_bound_force(self) -> float:
        mask = self.bead_index >= 0
        return float(self.force_magnitude[mask].mean()) if np.any(mask) else 0.0

    def _nearest_bead(self, target: np.ndarray, positions: np.ndarray) -> int:
        delta = positions - target
        distances_sq = np.einsum("ij,ij->i", delta, delta)
        candidate = int(np.argmin(distances_sq))
        if distances_sq[candidate] <= self.parameters.capture_radius ** 2:
            return candidate
        return -1

    def step(
        self,
        network_positions: np.ndarray,
        dt: float,
        rng: np.random.Generator,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Advance clutch states and return bead forces and cell reaction force."""

        p = self.parameters
        bead_forces = np.zeros_like(network_positions)
        cell_force = np.zeros(2, dtype=float)

        # Force–velocity relation for the module's myosin-driven actin flow.
        previous_module_force = self.force_magnitude.sum(axis=1)
        stall = max(p.motors_per_module * p.motor_stall_force, 1.0e-12)
        actin_speed = p.unloaded_actin_speed * np.clip(1.0 - previous_module_force / stall, 0.0, 1.0)
        self.actin_slip += actin_speed * dt
        protrusion_tips = self.position[None, :] + p.cell_radius * self.directions
        actin_anchors = protrusion_tips - self.actin_slip[:, None] * self.directions

        self.force_magnitude.fill(0.0)
        for module in range(p.n_modules):
            anchor = actin_anchors[module]
            for clutch in range(p.clutches_per_module):
                bead = int(self.bead_index[module, clutch])
                if bead >= 0:
                    extension = anchor - network_positions[bead] - self.rest_offset[module, clutch]
                    force_on_bead = p.clutch_stiffness * extension
                    magnitude = float(np.linalg.norm(force_on_bead))
                    off_rate = p.unloaded_off_rate * exp(min(magnitude / p.rupture_force, 30.0))
                    if rng.random() < 1.0 - exp(-off_rate * dt):
                        self.completed_lifetimes.append(self.time - self.bind_time[module, clutch])
                        self.bead_index[module, clutch] = -1
                        continue
                    bead_forces[bead] += force_on_bead
                    cell_force -= force_on_bead
                    self.force_magnitude[module, clutch] = magnitude
                elif rng.random() < 1.0 - exp(-p.on_rate * dt):
                    bead = self._nearest_bead(protrusion_tips[module], network_positions)
                    if bead >= 0:
                        self.bead_index[module, clutch] = bead
                        self.rest_offset[module, clutch] = anchor - network_positions[bead]
                        self.bind_time[module, clutch] = self.time

        self.position += dt * cell_force / p.cell_drag
        self.time += dt
        return bead_forces, cell_force

