"""State containers for spatial clutches, protrusions, and the rigid G3 cell."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class RigidCellState:
    center: np.ndarray
    body_angle: float
    radius: float
    velocity: np.ndarray
    angular_velocity: float = 0.0

    @classmethod
    def at_origin(cls, radius: float) -> "RigidCellState":
        return cls(np.zeros(2, dtype=float), 0.0, radius, np.zeros(2, dtype=float), 0.0)


@dataclass
class ClutchState:
    """Vectorized clutch states; material coordinates persist until Bell unbinding."""

    bound: np.ndarray
    fiber_id: np.ndarray
    segment_id: np.ndarray
    alpha: np.ndarray
    body_anchor_angle: np.ndarray
    rest_length: np.ndarray
    actin_displacement: np.ndarray
    force_vector: np.ndarray
    binding_time: np.ndarray
    sector_id: np.ndarray

    @classmethod
    def empty(cls, n_clutches: int) -> "ClutchState":
        return cls(
            bound=np.zeros(n_clutches, dtype=bool),
            fiber_id=np.full(n_clutches, -1, dtype=int),
            segment_id=np.full(n_clutches, -1, dtype=int),
            alpha=np.zeros(n_clutches, dtype=float),
            body_anchor_angle=np.zeros(n_clutches, dtype=float),
            rest_length=np.zeros(n_clutches, dtype=float),
            actin_displacement=np.zeros(n_clutches, dtype=float),
            force_vector=np.zeros((n_clutches, 2), dtype=float),
            binding_time=np.full(n_clutches, np.nan, dtype=float),
            sector_id=np.full(n_clutches, -1, dtype=int),
        )

    @property
    def n_clutches(self) -> int:
        return int(self.bound.size)

    def detach(self, indices: np.ndarray) -> None:
        indices = np.asarray(indices, dtype=int)
        if indices.size == 0:
            return
        self.bound[indices] = False
        self.fiber_id[indices] = -1
        self.segment_id[indices] = -1
        self.alpha[indices] = 0.0
        self.rest_length[indices] = 0.0
        self.actin_displacement[indices] = 0.0
        self.force_vector[indices] = 0.0
        self.binding_time[indices] = np.nan


@dataclass
class ProtrusionState:
    sector_angles: np.ndarray
    active: np.ndarray
    geometry_score: np.ndarray
    traction_score: np.ndarray

    @classmethod
    def initialize(
        cls, n_sectors: int, n_active: int, rng: np.random.Generator, prescribed=None
    ) -> "ProtrusionState":
        angles = 2.0 * np.pi * np.arange(n_sectors) / n_sectors
        active = np.zeros(n_sectors, dtype=bool)
        if prescribed is None:
            chosen = rng.choice(n_sectors, size=n_active, replace=False)
        else:
            chosen = np.atleast_1d(prescribed).astype(int)
        active[chosen] = True
        return cls(angles, active, np.zeros(n_sectors), np.zeros(n_sectors))


@dataclass
class G3Snapshot:
    time: float
    positions: np.ndarray
    cell_center: np.ndarray
    cell_angle: float
    bound_points: np.ndarray
    motor_points: np.ndarray
    clutch_forces: np.ndarray
    active_sectors: np.ndarray
    geometry_scores: np.ndarray
    traction_scores: np.ndarray
    foi: float
    cell_force: np.ndarray
    cell_torque: float

