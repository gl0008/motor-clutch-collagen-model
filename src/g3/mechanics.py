"""Spatial motor-clutch mechanics and moment-preserving local force projection."""

from __future__ import annotations

import numpy as np

from .config import G3Config
from .elastic import NUMBA_ENABLED, bending_forces, extensional_forces, njit, overdamped_step
from .fixtures import G3Fixture
from .state import ClutchState, ProtrusionState, RigidCellState


_EXPONENT_CLAMP = 40.0


def cross2(a, b):
    """Scalar z-component of a 2D cross product; broadcasts over leading axes."""
    a = np.asarray(a)
    b = np.asarray(b)
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def bell_off_rate(force_magnitude, cfg: G3Config):
    """Bell slip-bond rate, Adebowale et al. 2021 SI Eq. 2."""
    argument = np.clip(np.asarray(force_magnitude) / cfg.bell_force,
                       0.0, _EXPONENT_CLAMP)
    return cfg.unbind_rate * np.exp(argument)


def force_velocity(load: float, motor_count: float, cfg: G3Config) -> float:
    """Hill motor force-velocity closure for one active protrusion module."""
    stall = motor_count * cfg.motor_force
    if stall <= 0.0:
        return 0.0
    return float(np.clip(cfg.unloaded_actin_speed * (1.0 - load / stall),
                         0.0, cfg.unloaded_actin_speed))


def assign_clutch_pool(clutches: ClutchState, protrusions: ProtrusionState) -> None:
    """Assign a fixed total clutch pool to active sectors while preserving surviving sites."""
    active = np.flatnonzero(protrusions.active)
    if active.size == 0:
        clutches.detach(np.flatnonzero(clutches.bound))
        clutches.sector_id[:] = -1
        return

    # In normal operation the assignment changes only when a protrusion turns over. Avoid
    # rebuilding identical masks, pools and anchor offsets on every 0.005 s physics step.
    if np.all(clutches.sector_id >= 0):
        base, remainder = divmod(clutches.n_clutches, active.size)
        expected = np.zeros(protrusions.active.size, dtype=int)
        expected[active] = base
        expected[active[:remainder]] += 1
        observed = np.bincount(
            clutches.sector_id, minlength=protrusions.active.size
        )[:protrusions.active.size]
        if np.array_equal(observed, expected):
            return

    inactive_assignment = ~np.isin(clutches.sector_id, active)
    changed = np.flatnonzero(inactive_assignment & (clutches.sector_id >= 0))
    clutches.detach(changed)
    clutches.sector_id[inactive_assignment] = -1

    base, remainder = divmod(clutches.n_clutches, active.size)
    targets = {int(sector): base + (rank < remainder)
               for rank, sector in enumerate(active)}
    pool = list(np.flatnonzero(clutches.sector_id < 0))

    for sector in active:
        assigned = np.flatnonzero(clutches.sector_id == sector)
        excess = max(0, assigned.size - targets[int(sector)])
        if excess:
            release = assigned[-excess:]
            clutches.detach(release)
            clutches.sector_id[release] = -1
            pool.extend(release.tolist())

    pool = sorted(set(pool))
    cursor = 0
    sector_width = 2.0 * np.pi / protrusions.sector_angles.size
    for sector in active:
        assigned = np.flatnonzero(clutches.sector_id == sector)
        deficit = targets[int(sector)] - assigned.size
        if deficit > 0:
            take = np.asarray(pool[cursor:cursor + deficit], dtype=int)
            cursor += deficit
            clutches.sector_id[take] = sector
            assigned = np.concatenate((assigned, take))
        offsets = np.linspace(-0.45 * sector_width, 0.45 * sector_width,
                              assigned.size, endpoint=True) if assigned.size > 1 else np.array([0.0])
        clutches.body_anchor_angle[assigned] = protrusions.sector_angles[sector] + offsets


def clutch_geometry_numpy(clutches: ClutchState, cell: RigidCellState, positions,
                          fibre_bonds):
    """Return evolving material points, surface anchors, and motor-side endpoints."""
    n = clutches.n_clutches
    material = np.full((n, 2), np.nan, dtype=float)
    bound_idx = np.flatnonzero(clutches.bound)
    if bound_idx.size:
        segments = fibre_bonds[clutches.segment_id[bound_idx]]
        alpha = clutches.alpha[bound_idx, None]
        material[bound_idx] = ((1.0 - alpha) * positions[segments[:, 0]]
                               + alpha * positions[segments[:, 1]])

    lab_angles = cell.body_angle + clutches.body_anchor_angle
    normals = np.column_stack((np.cos(lab_angles), np.sin(lab_angles)))
    anchors = cell.center[None, :] + cell.radius * normals
    motor_points = anchors - clutches.actin_displacement[:, None] * normals
    return material, anchors, motor_points, normals


if njit is not None:
    @njit(cache=True)
    def _clutch_geometry_numba(bound, segment_id, alpha, body_anchor_angle,
                               actin_displacement, center, body_angle, radius,
                               positions, fibre_bonds):
        n_clutches = bound.size
        material = np.empty((n_clutches, 2), dtype=np.float64)
        material[:, :] = np.nan
        anchors = np.empty((n_clutches, 2), dtype=np.float64)
        motor_points = np.empty((n_clutches, 2), dtype=np.float64)
        normals = np.empty((n_clutches, 2), dtype=np.float64)
        for clutch_id in range(n_clutches):
            angle = body_angle + body_anchor_angle[clutch_id]
            normal_x = np.cos(angle)
            normal_y = np.sin(angle)
            normals[clutch_id, 0] = normal_x
            normals[clutch_id, 1] = normal_y
            anchors[clutch_id, 0] = center[0] + radius * normal_x
            anchors[clutch_id, 1] = center[1] + radius * normal_y
            motor_points[clutch_id, 0] = (
                anchors[clutch_id, 0] - actin_displacement[clutch_id] * normal_x
            )
            motor_points[clutch_id, 1] = (
                anchors[clutch_id, 1] - actin_displacement[clutch_id] * normal_y
            )
            if bound[clutch_id]:
                segment = segment_id[clutch_id]
                first = fibre_bonds[segment, 0]
                second = fibre_bonds[segment, 1]
                fraction = alpha[clutch_id]
                material[clutch_id, 0] = (
                    (1.0 - fraction) * positions[first, 0]
                    + fraction * positions[second, 0]
                )
                material[clutch_id, 1] = (
                    (1.0 - fraction) * positions[first, 1]
                    + fraction * positions[second, 1]
                )
        return material, anchors, motor_points, normals


def clutch_geometry(clutches: ClutchState, cell: RigidCellState, positions,
                    fibre_bonds):
    if NUMBA_ENABLED:
        return _clutch_geometry_numba(
            clutches.bound, clutches.segment_id, clutches.alpha,
            clutches.body_anchor_angle, clutches.actin_displacement,
            cell.center, cell.body_angle, cell.radius, positions, fibre_bonds,
        )
    return clutch_geometry_numpy(clutches, cell, positions, fibre_bonds)


def _point_forces_numpy(clutches: ClutchState, material, motor_points, cfg: G3Config):
    forces = np.zeros_like(clutches.force_vector)
    idx = np.flatnonzero(clutches.bound)
    if idx.size == 0:
        return forces
    vector = motor_points[idx] - material[idx]
    distance = np.linalg.norm(vector, axis=1)
    extension = np.maximum(distance - clutches.rest_length[idx], 0.0)
    nonzero = distance > 0.0
    forces[idx[nonzero]] = (cfg.clutch_stiffness * extension[nonzero] / distance[nonzero])[:, None] * vector[nonzero]
    return forces


if njit is not None:
    @njit(cache=True)
    def _point_forces_numba(bound, rest_length, material, motor_points, stiffness):
        forces = np.zeros_like(motor_points)
        for clutch_id in range(bound.size):
            if not bound[clutch_id]:
                continue
            dx = motor_points[clutch_id, 0] - material[clutch_id, 0]
            dy = motor_points[clutch_id, 1] - material[clutch_id, 1]
            distance = np.sqrt(dx * dx + dy * dy)
            extension = max(distance - rest_length[clutch_id], 0.0)
            if distance > 0.0:
                scale = stiffness * extension / distance
                forces[clutch_id, 0] = scale * dx
                forces[clutch_id, 1] = scale * dy
        return forces


def _point_forces(clutches: ClutchState, material, motor_points, cfg: G3Config):
    if NUMBA_ENABLED:
        return _point_forces_numba(
            clutches.bound, clutches.rest_length, material, motor_points,
            cfg.clutch_stiffness,
        )
    return _point_forces_numpy(clutches, material, motor_points, cfg)


def closest_material_point(anchor, positions, fibre_bonds, segment_fiber_id):
    """Nearest point on any fibre segment as (segment index, fibre id, alpha, point, distance)."""
    if fibre_bonds.size == 0:
        return -1, -1, 0.0, np.full(2, np.nan), np.inf
    a = positions[fibre_bonds[:, 0]]
    b = positions[fibre_bonds[:, 1]]
    edge = b - a
    denom = np.einsum("ij,ij->i", edge, edge)
    numerator = np.einsum("ij,ij->i", anchor[None, :] - a, edge)
    alpha = np.divide(numerator, denom, out=np.zeros_like(numerator), where=denom > 0.0)
    alpha = np.clip(alpha, 0.0, 1.0)
    points = a + alpha[:, None] * edge
    distance = np.linalg.norm(points - anchor[None, :], axis=1)
    segment = int(np.argmin(distance))
    return (segment, int(segment_fiber_id[segment]), float(alpha[segment]),
            points[segment], float(distance[segment]))


def update_spatial_clutches(
    clutches: ClutchState,
    protrusions: ProtrusionState,
    cell: RigidCellState,
    positions: np.ndarray,
    fixture: G3Fixture,
    cfg: G3Config,
    time: float,
    rng: np.random.Generator,
):
    """Advance Bell binding/unbinding and motor-side motion by one timestep."""
    assign_clutch_pool(clutches, protrusions)
    material, anchors, motor_points, _ = clutch_geometry(
        clutches, cell, positions, fixture.network.fiber_bonds)
    current = _point_forces(clutches, material, motor_points, cfg)

    # Bell rupture uses the force at the beginning of the step.
    idx = np.flatnonzero(clutches.bound)
    if idx.size:
        rates = bell_off_rate(np.linalg.norm(current[idx], axis=1), cfg)
        probability = 1.0 - np.exp(-rates * cfg.dt)
        broken = idx[rng.random(idx.size) < probability]
        clutches.detach(broken)

    # Binding trials are rare (rate * dt), so perform segment searches only on successes.
    eligible = np.flatnonzero((clutches.sector_id >= 0) & ~clutches.bound)
    if eligible.size:
        p_bind = 1.0 - np.exp(-cfg.bind_rate * cfg.dt)
        trials = eligible[rng.random(eligible.size) < p_bind]
        for clutch_id in trials:
            lab_angle = cell.body_angle + clutches.body_anchor_angle[clutch_id]
            normal = np.array([np.cos(lab_angle), np.sin(lab_angle)])
            anchor = cell.center + cell.radius * normal
            segment, fibre, alpha, point, distance = closest_material_point(
                anchor, positions, fixture.network.fiber_bonds, fixture.segment_fiber_id)
            if distance <= cfg.capture_distance:
                clutches.bound[clutch_id] = True
                clutches.segment_id[clutch_id] = segment
                clutches.fiber_id[clutch_id] = fibre
                clutches.alpha[clutch_id] = alpha
                clutches.actin_displacement[clutch_id] = 0.0
                clutches.rest_length[clutch_id] = distance
                clutches.binding_time[clutch_id] = time
                clutches.force_vector[clutch_id] = 0.0

    material, anchors, motor_points, _ = clutch_geometry(
        clutches, cell, positions, fixture.network.fiber_bonds)
    current = _point_forces(clutches, material, motor_points, cfg)

    active = np.flatnonzero(protrusions.active)
    speeds = np.zeros(protrusions.active.size, dtype=float)
    motors_per_sector = cfg.n_motors / max(active.size, 1)
    for sector in active:
        mask = clutches.bound & (clutches.sector_id == sector)
        load = float(np.linalg.norm(current[mask], axis=1).sum())
        speeds[sector] = force_velocity(load, motors_per_sector, cfg)
        clutches.actin_displacement[mask] += speeds[sector] * cfg.dt

    material, anchors, motor_points, _ = clutch_geometry(
        clutches, cell, positions, fixture.network.fiber_bonds)
    clutches.force_vector[:] = _point_forces(clutches, material, motor_points, cfg)
    return material, anchors, motor_points, speeds


def project_clutch_forces_numpy(
    clutches: ClutchState,
    material_points: np.ndarray,
    positions: np.ndarray,
    fixture: G3Fixture,
    cfg: G3Config,
):
    """Project local point forces to same-fibre beads, preserving force and first moment."""
    bead_forces = np.zeros_like(positions)
    bonds = fixture.network.fiber_bonds
    for clutch_id in np.flatnonzero(clutches.bound):
        force = clutches.force_vector[clutch_id]
        if not np.any(force):
            continue
        point = material_points[clutch_id]
        fibre = clutches.fiber_id[clutch_id]
        bead_ids = fixture.bead_ids_by_fibre[fibre]
        distance = np.linalg.norm(positions[bead_ids] - point[None, :], axis=1)
        support = distance <= cfg.gaussian_support_sigma * cfg.gaussian_sigma
        local_ids = bead_ids[support]
        local_distance = distance[support]
        if local_ids.size == 0:
            local_ids = bonds[clutches.segment_id[clutch_id]]
            local_distance = np.linalg.norm(positions[local_ids] - point[None, :], axis=1)
        weights = np.exp(-0.5 * (local_distance / cfg.gaussian_sigma) ** 2)
        weights /= weights.sum()
        projected = weights[:, None] * force
        np.add.at(bead_forces, local_ids, projected)

        # Correct the Gaussian first moment without changing its net force.
        raw_moment = float(np.sum(cross2(positions[local_ids] - point, projected)))
        i, j = bonds[clutches.segment_id[clutch_id]]
        edge = positions[j] - positions[i]
        edge_sq = float(np.dot(edge, edge))
        if edge_sq > 0.0 and raw_moment != 0.0:
            perpendicular = np.array([-edge[1], edge[0]])
            correction_at_j = (-raw_moment / edge_sq) * perpendicular
            bead_forces[i] -= correction_at_j
            bead_forces[j] += correction_at_j

    return bead_forces


if njit is not None:
    @njit(cache=True)
    def _project_clutch_forces_numba(bound, fibre_ids, segment_ids, force_vectors,
                                     material_points, positions, bonds, beads_per_fibre,
                                     gaussian_sigma, support_sigma):
        bead_forces = np.zeros_like(positions)
        support_distance = support_sigma * gaussian_sigma
        for clutch_id in range(bound.size):
            if not bound[clutch_id]:
                continue
            force_x = force_vectors[clutch_id, 0]
            force_y = force_vectors[clutch_id, 1]
            if force_x == 0.0 and force_y == 0.0:
                continue
            point_x = material_points[clutch_id, 0]
            point_y = material_points[clutch_id, 1]
            first_bead = fibre_ids[clutch_id] * beads_per_fibre
            last_bead = min(first_bead + beads_per_fibre, positions.shape[0])
            weight_sum = 0.0
            for bead_id in range(first_bead, last_bead):
                dx = positions[bead_id, 0] - point_x
                dy = positions[bead_id, 1] - point_y
                distance = np.sqrt(dx * dx + dy * dy)
                if distance <= support_distance:
                    weight_sum += np.exp(-0.5 * (distance / gaussian_sigma) ** 2)
            raw_moment = 0.0
            if weight_sum > 0.0:
                for bead_id in range(first_bead, last_bead):
                    dx = positions[bead_id, 0] - point_x
                    dy = positions[bead_id, 1] - point_y
                    distance = np.sqrt(dx * dx + dy * dy)
                    if distance <= support_distance:
                        weight = np.exp(-0.5 * (distance / gaussian_sigma) ** 2) / weight_sum
                        projected_x = weight * force_x
                        projected_y = weight * force_y
                        bead_forces[bead_id, 0] += projected_x
                        bead_forces[bead_id, 1] += projected_y
                        raw_moment += dx * projected_y - dy * projected_x
            segment = segment_ids[clutch_id]
            first = bonds[segment, 0]
            second = bonds[segment, 1]
            edge_x = positions[second, 0] - positions[first, 0]
            edge_y = positions[second, 1] - positions[first, 1]
            edge_sq = edge_x * edge_x + edge_y * edge_y
            if edge_sq > 0.0 and raw_moment != 0.0:
                correction_x = (-raw_moment / edge_sq) * (-edge_y)
                correction_y = (-raw_moment / edge_sq) * edge_x
                bead_forces[first, 0] -= correction_x
                bead_forces[first, 1] -= correction_y
                bead_forces[second, 0] += correction_x
                bead_forces[second, 1] += correction_y
        return bead_forces


def project_clutch_forces(
    clutches: ClutchState,
    material_points: np.ndarray,
    positions: np.ndarray,
    fixture: G3Fixture,
    cfg: G3Config,
):
    if NUMBA_ENABLED:
        return _project_clutch_forces_numba(
            clutches.bound, clutches.fiber_id, clutches.segment_id,
            clutches.force_vector, material_points, positions,
            fixture.network.fiber_bonds, cfg.beads_per_fibre,
            cfg.gaussian_sigma, cfg.gaussian_support_sigma,
        )
    return project_clutch_forces_numpy(clutches, material_points, positions, fixture, cfg)


def reaction_force_and_torque(clutches: ClutchState, motor_points, cell: RigidCellState):
    reaction = -clutches.force_vector
    force = reaction.sum(axis=0)
    torque = float(np.sum(cross2(motor_points - cell.center[None, :], reaction)))
    return force, torque


def elastic_forces(positions, fixture: G3Fixture, cfg: G3Config):
    network = fixture.network
    return (
        extensional_forces(positions, network.fiber_bonds, cfg.bead_spacing, cfg.kappa_s_f)
        + bending_forces(positions, network.bending_triples, cfg.theta0_f, cfg.kappa_b_f)
    )


def advance_ecm(positions, active_forces, fixture: G3Fixture, cfg: G3Config,
                rng: np.random.Generator):
    updated = positions
    sub_dt = cfg.dt / cfg.ecm_substeps
    for _ in range(cfg.ecm_substeps):
        total = elastic_forces(updated, fixture, cfg) + active_forces
        updated = overdamped_step(updated, total, cfg.bead_drag, sub_dt)
        updated[fixture.fixed_mask] = fixture.initial_positions[fixture.fixed_mask]
    return updated


def advance_cell(cell: RigidCellState, force, torque: float, cfg: G3Config) -> None:
    cell.velocity = np.asarray(force, dtype=float) / cfg.cell_drag
    cell.angular_velocity = float(torque / cfg.rotational_drag)
    cell.center = cell.center + cell.velocity * cfg.dt
    cell.body_angle = float((cell.body_angle + cell.angular_velocity * cfg.dt + np.pi)
                            % (2.0 * np.pi) - np.pi)


def conservation_errors(clutches, material_points, motor_points, bead_forces, positions, cell):
    """Relative internal force and torque residuals for diagnostics/tests."""
    point_force = clutches.force_vector.sum(axis=0)
    reaction_force, reaction_torque = reaction_force_and_torque(clutches, motor_points, cell)
    force_residual = bead_forces.sum(axis=0) + reaction_force
    force_scale = max(float(np.linalg.norm(point_force)), np.finfo(float).eps)
    force_error = float(np.linalg.norm(force_residual) / force_scale)

    ecm_moment = float(np.sum(cross2(positions, bead_forces)))
    cell_moment = float(cross2(cell.center, reaction_force) + reaction_torque)
    bound = clutches.bound
    if np.any(bound):
        point_scale = float(np.sum(
            np.linalg.norm(clutches.force_vector[bound], axis=1)
            * np.maximum(np.linalg.norm(material_points[bound], axis=1), cell.radius)
        ))
    else:
        point_scale = 0.0
    torque_error = abs(ecm_moment + cell_moment) / max(point_scale, np.finfo(float).eps)
    return force_error, float(torque_error)
