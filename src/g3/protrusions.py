"""Coarse-grained collagen sensing and traction-reinforced protrusion turnover."""

from __future__ import annotations

import numpy as np

from .config import G3Config
from .fixtures import G3Fixture
from .state import ClutchState, ProtrusionState, RigidCellState


def protrusion_tips(protrusions: ProtrusionState, cell: RigidCellState):
    """Return explicit protrusion bases and tips in the laboratory frame."""
    angles = cell.body_angle + protrusions.sector_angles
    normal = np.column_stack((np.cos(angles), np.sin(angles)))
    base = cell.center[None, :] + cell.radius * normal
    tip = base + protrusions.length[:, None] * normal
    return base, tip


def step_protrusion_lengths(
    protrusions: ProtrusionState,
    cfg: G3Config,
    dt: float | None = None,
    clutches: ClutchState | None = None,
):
    """Grow active protrusions and retract inactive probes at a finite speed."""
    dt = cfg.dt if dt is None else float(dt)
    maximum = max(float(protrusions.activity.max(initial=0.0)), np.finfo(float).eps)
    activity = protrusions.activity / maximum
    activity = np.where(protrusions.active, np.maximum(activity, 0.35), activity)
    growth = cfg.protrusion_growth_speed * activity * (
        1.0 - protrusions.length / cfg.protrusion_max_length
    )
    retraction = cfg.protrusion_retraction_speed * (~protrusions.active)
    if clutches is not None and np.any(clutches.bound):
        attached = np.unique(clutches.sector_id[clutches.bound & (clutches.sector_id >= 0)])
        growth[attached] = 0.0
        retraction[attached] = 0.0
    protrusions.length += dt * (growth - retraction)
    protrusions.length[:] = np.clip(protrusions.length, 0.0, cfg.protrusion_max_length)


def step_intrinsic_polarity(
    protrusions: ProtrusionState,
    cfg: G3Config,
    rng: np.random.Generator,
    dt: float,
    feedback_enabled: bool = True,
):
    """Mass-conserving stochastic cell polarity with adhesion reinforcement.

    Collagen geometry is deliberately absent from this update. Matrix structure
    affects polarity only after a protrusion physically encounters a fibre and
    develops a nonzero traction-success signal.
    """
    activity = protrusions.activity
    adhesion = protrusions.traction_score if feedback_enabled else np.zeros_like(activity)
    laplacian = np.roll(activity, 1) + np.roll(activity, -1) - 2.0 * activity
    logits = (
        cfg.polarity_self_gain * activity
        + cfg.polarity_adhesion_gain * adhesion
        + cfg.polarity_diffusion * laplacian
    )
    failed_probe = (
        (protrusions.length >= 0.80 * cfg.protrusion_max_length)
        & (protrusions.traction_score < 0.02)
    )
    logits -= cfg.polarity_failed_probe_penalty * failed_probe
    logits -= float(logits.max(initial=0.0))
    target = np.exp(np.clip(logits, -40.0, 0.0))
    target *= cfg.polarity_total_activity / max(float(target.sum()), np.finfo(float).eps)
    alpha = min(float(dt) / cfg.polarity_time, 1.0)
    noise = cfg.polarity_noise * np.sqrt(max(alpha, 0.0)) * rng.normal(size=activity.size)
    activity += alpha * (target - activity) + noise / activity.size
    activity[:] = np.maximum(activity, 0.0)
    total = float(activity.sum())
    if total <= np.finfo(float).eps:
        activity[:] = cfg.polarity_total_activity / activity.size
    else:
        activity *= cfg.polarity_total_activity / total
    # Activity determines where a protrusion is *eligible* to form, but a
    # protrusion is a physical object and must not teleport whenever two noisy
    # sector values cross. Keep expressed sectors for a minimum lifetime and
    # require both a relative and absolute activity advantage before replacing
    # the weakest mature protrusion.
    n_active = min(cfg.n_active_protrusions, activity.size)
    protrusions.active_age[protrusions.active] += dt
    protrusions.active_age[~protrusions.active] = 0.0
    active = np.flatnonzero(protrusions.active)
    if active.size < n_active:
        candidates = np.flatnonzero(~protrusions.active)
        add = candidates[np.argsort(activity[candidates])[-(n_active - active.size):]]
        protrusions.active[add] = True
        protrusions.active_age[add] = 0.0
        active = np.flatnonzero(protrusions.active)
    elif active.size > n_active:
        keep = active[np.argsort(activity[active])[-n_active:]]
        protrusions.active[:] = False
        protrusions.active[keep] = True
        protrusions.active_age[:] = 0.0
        active = keep

    inactive = np.flatnonzero(~protrusions.active)
    mature = active[protrusions.active_age[active] >= cfg.protrusion_min_lifetime]
    if mature.size and inactive.size:
        weakest = int(mature[np.argmin(activity[mature])])
        challenger = int(inactive[np.argmax(activity[inactive])])
        threshold = max(
            cfg.protrusion_switch_ratio * activity[weakest],
            activity[weakest] + cfg.protrusion_switch_margin,
        )
        if activity[challenger] > threshold:
            protrusions.active[weakest] = False
            protrusions.active[challenger] = True
            protrusions.active_age[weakest] = 0.0
            protrusions.active_age[challenger] = 0.0
    return np.flatnonzero(protrusions.active)


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
        n_bound = int(bound.sum())
        traction = float(np.linalg.norm(clutches.force_vector[bound], axis=1).sum())
        engagement = 1.0 - np.exp(-n_bound / cfg.adhesion_clutch_scale)
        force_scale = cfg.adhesion_force_fraction * stall
        load_success = min(traction / max(force_scale, np.finfo(float).eps), 1.0)
        instant[sector] = engagement * load_success
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
