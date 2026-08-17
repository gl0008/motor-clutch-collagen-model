"""Coarse-grained collagen sensing and traction-reinforced protrusion turnover."""

from __future__ import annotations

import numpy as np

from .config import G3Config
from .fixtures import G3Fixture
from .state import ClutchState, ProtrusionState, RigidCellState


def geometry_scores(protrusions: ProtrusionState, cell: RigidCellState,
                    positions: np.ndarray, fixture: G3Fixture, cfg: G3Config):
    """Compute collagen availability, local nematic alignment, and their combined score."""
    bonds = fixture.network.fiber_bonds
    n_sectors = protrusions.sector_angles.size
    if bonds.size == 0:
        zeros = np.zeros(n_sectors, dtype=float)
        protrusions.geometry_score[:] = 0.0
        return zeros, zeros.copy(), zeros.copy()

    a = positions[bonds[:, 0]]
    b = positions[bonds[:, 1]]
    midpoint = 0.5 * (a + b)
    tangent = b - a
    length = np.linalg.norm(tangent, axis=1)
    tangent = np.divide(tangent, length[:, None], out=np.zeros_like(tangent),
                        where=length[:, None] > 0.0)

    lab_angles = cell.body_angle + protrusions.sector_angles
    normal = np.column_stack((np.cos(lab_angles), np.sin(lab_angles)))
    anchors = cell.center[None, :] + cell.radius * normal
    delta = midpoint[None, :, :] - anchors[:, None, :]
    distance = np.linalg.norm(delta, axis=2)
    kernel = np.exp(-0.5 * (distance / cfg.sensing_sigma) ** 2)
    availability = kernel.sum(axis=1)
    orientation = (normal @ tangent.T) ** 2
    denom = kernel.sum(axis=1)
    alignment = np.divide((kernel * orientation).sum(axis=1), denom,
                          out=np.zeros(n_sectors), where=denom > 0.0)
    normalized_availability = availability / max(float(availability.max()), np.finfo(float).eps)
    combined = normalized_availability * alignment
    protrusions.geometry_score[:] = combined
    return availability, alignment, combined


def update_traction_scores(protrusions: ProtrusionState, clutches: ClutchState,
                           cfg: G3Config):
    active = np.flatnonzero(protrusions.active)
    alpha = min(cfg.dt / cfg.feedback_time, 1.0)
    instant = np.zeros_like(protrusions.traction_score)
    motors_per_sector = cfg.n_motors / max(active.size, 1)
    stall = motors_per_sector * cfg.motor_force
    for sector in active:
        assigned = clutches.sector_id == sector
        n_assigned = int(assigned.sum())
        if n_assigned == 0:
            continue
        bound = assigned & clutches.bound
        bound_fraction = float(bound.sum() / n_assigned)
        traction = float(np.linalg.norm(clutches.force_vector[bound], axis=1).sum())
        instant[sector] = bound_fraction * min(traction / max(stall, np.finfo(float).eps), 1.0)
    protrusions.traction_score += alpha * (instant - protrusions.traction_score)
    protrusions.traction_score[~protrusions.active] = 0.0
    return instant


def step_protrusion_turnover(protrusions: ProtrusionState, cfg: G3Config,
                             rng: np.random.Generator, feedback_enabled: bool = True):
    """Apply persistence hazard and replace lost protrusions with no global polar bias."""
    active_before = np.flatnonzero(protrusions.active)
    if active_before.size == 0:
        return False
    beta_g = cfg.beta_geometry if feedback_enabled else 0.0
    beta_q = cfg.beta_traction if feedback_enabled else 0.0
    hazard = (1.0 / cfg.protrusion_lifetime) * np.exp(
        -beta_g * protrusions.geometry_score[active_before]
        -beta_q * protrusions.traction_score[active_before]
    )
    probability = 1.0 - np.exp(-hazard * cfg.dt)
    lost = active_before[rng.random(active_before.size) < probability]
    if lost.size == 0:
        return False

    changed = False
    for old_sector in lost:
        protrusions.active[old_sector] = False
        protrusions.traction_score[old_sector] = 0.0
        candidates = np.flatnonzero(~protrusions.active)
        if candidates.size == 0:
            protrusions.active[old_sector] = True
            continue
        logits = beta_g * protrusions.geometry_score[candidates]
        logits -= logits.max()
        weights = np.exp(logits)
        weights /= weights.sum()
        replacement = int(rng.choice(candidates, p=weights))
        protrusions.active[replacement] = True
        changed = changed or replacement != old_sector
    return changed


def protrusion_axis(protrusions: ProtrusionState, cell: RigidCellState,
                    use_traction: bool = True):
    """Polar resultant of active protrusions, used only as an analysis observable."""
    active = np.flatnonzero(protrusions.active)
    if active.size == 0:
        return np.zeros(2)
    angles = cell.body_angle + protrusions.sector_angles[active]
    vectors = np.column_stack((np.cos(angles), np.sin(angles)))
    if use_traction:
        weights = 1.0 + protrusions.traction_score[active]
    else:
        weights = np.ones(active.size)
    return np.sum(weights[:, None] * vectors, axis=0)
