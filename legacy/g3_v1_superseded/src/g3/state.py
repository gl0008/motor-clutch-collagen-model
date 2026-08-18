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
    # The controller start is persistent across individual clutch ruptures so a
    # declared protocol ramps from the first *network attachment*, not from the
    # most recently rebound molecular clutch.
    traction_start_time: float

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
            traction_start_time=float("nan"),
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
    activity: np.ndarray
    length: np.ndarray
    active_age: np.ndarray

    @classmethod
    def initialize(
        cls, n_sectors: int, n_active: int, rng: np.random.Generator, prescribed=None
    ) -> "ProtrusionState":
        angles = 2.0 * np.pi * np.arange(n_sectors) / n_sectors
        active = np.zeros(n_sectors, dtype=bool)
        activity = np.full(n_sectors, 1.0 / n_sectors, dtype=float)
        if prescribed is None:
            activity += rng.normal(0.0, 1.0e-3 / n_sectors, size=n_sectors)
            activity = np.maximum(activity, 0.0)
            activity /= activity.sum()
            chosen = np.argpartition(activity, -n_active)[-n_active:]
        else:
            chosen = np.atleast_1d(prescribed).astype(int)
            activity[:] = 0.0
            activity[chosen] = 1.0 / max(chosen.size, 1)
        active[chosen] = True
        return cls(
            angles,
            active,
            np.zeros(n_sectors),
            np.zeros(n_sectors),
            activity,
            np.zeros(n_sectors),
            np.zeros(n_sectors),
        )


@dataclass
class G3Snapshot:
    time: float
    positions: np.ndarray
    cell_center: np.ndarray
    cell_angle: float
    bound_points: np.ndarray
    bound_fibre_ids: np.ndarray
    bound_sector_ids: np.ndarray
    motor_points: np.ndarray
    clutch_forces: np.ndarray
    active_sectors: np.ndarray
    active_sector_age: np.ndarray
    protrusion_tips: np.ndarray
    protrusion_lengths: np.ndarray
    polarity_activity: np.ndarray
    geometry_scores: np.ndarray
    traction_scores: np.ndarray
    foi: float
    cell_force: np.ndarray
    cell_torque: float
