"""G3 observables: FOI, Nam plasticity index, displacement, and directional statistics."""

from __future__ import annotations

import numpy as np

from .elastic import bending_energy, extensional_energy

from .config import G3Config
from .fixtures import G3Fixture


FOI_RANDOM = 2.0 / np.pi


def fibre_geometry(positions, fixture: G3Fixture):
    """Return segment midpoints and tangents for local, deformation-sensitive FOI."""
    bonds = fixture.network.fiber_bonds
    if bonds.size == 0:
        return np.empty((0, 2)), np.empty((0, 2))
    a = positions[bonds[:, 0]]
    b = positions[bonds[:, 1]]
    vector = b - a
    length = np.linalg.norm(vector, axis=1)
    valid = length > 0.0
    centroids = 0.5 * (a[valid] + b[valid])
    directions = vector[valid] / length[valid, None]
    return centroids, directions


def fibre_orientation_index(positions, fixture: G3Fixture, cell_center,
                            max_distance: float | None = None) -> float:
    """Nam-style absolute fibre/radial overlap; random 2D baseline is 2/pi."""
    centroids, directions = fibre_geometry(positions, fixture)
    if centroids.size == 0:
        return float("nan")
    radial = centroids - np.asarray(cell_center)[None, :]
    distance = np.linalg.norm(radial, axis=1)
    valid = distance > 0.0
    if max_distance is not None:
        valid &= distance <= max_distance
    if not np.any(valid):
        return float("nan")
    radial_hat = radial[valid] / distance[valid, None]
    return float(np.mean(np.abs(np.sum(directions[valid] * radial_hat, axis=1))))


def foi_profile(positions, fixture: G3Fixture, cell_center, radial_edges):
    centroids, directions = fibre_geometry(positions, fixture)
    radial = centroids - np.asarray(cell_center)[None, :]
    distance = np.linalg.norm(radial, axis=1)
    values = np.zeros_like(distance)
    valid = distance > 0.0
    values[valid] = np.abs(np.sum(directions[valid] * radial[valid] / distance[valid, None], axis=1))
    radial_edges = np.asarray(radial_edges)
    centers = 0.5 * (radial_edges[:-1] + radial_edges[1:])
    output = np.full(centers.size, np.nan)
    for i in range(centers.size):
        mask = (distance >= radial_edges[i]) & (distance < radial_edges[i + 1])
        if np.any(mask):
            output[i] = values[mask].mean()
    return centers, output


def nam_plasticity_index(foi_pre_unload, foi_post, reference=FOI_RANDOM):
    """Nam et al. 2016 Eq. 7, generalized to scalar or spatial FOI arrays."""
    pre = np.asarray(foi_pre_unload, dtype=float) - reference
    post = np.asarray(foi_post, dtype=float) - reference
    denominator = float(np.sqrt(np.nanmean(pre**2)))
    if denominator <= np.finfo(float).eps:
        return float("nan")
    return float(np.sqrt(np.nanmean(post**2)) / denominator)


def displacement_statistics(initial, current):
    displacement = np.linalg.norm(np.asarray(current) - np.asarray(initial), axis=1)
    if displacement.size == 0:
        return {"median": 0.0, "p95": 0.0, "p99": 0.0, "max": 0.0,
                "participation_ratio": 0.0}
    squared = displacement**2
    denominator = displacement.size * np.sum(squared**2)
    participation = float(np.sum(squared) ** 2 / denominator) if denominator > 0.0 else 0.0
    return {
        "median": float(np.median(displacement)),
        "p95": float(np.percentile(displacement, 95.0)),
        "p99": float(np.percentile(displacement, 99.0)),
        "max": float(displacement.max()),
        "participation_ratio": participation,
    }


def elastic_energy(positions, fixture: G3Fixture, cfg: G3Config):
    network = fixture.network
    if network.external_network is not None:
        external = network.external_network
        external.r[:] = np.asarray(positions) * 1.0e6
        _, _, _, energy, _, _ = external.elastic_forces(include_crosslinks=True)
        # nN um -> J
        return float(sum(energy.values()) * 1.0e-15)
    link_energy = (
        extensional_energy(
            positions, network.xl_bonds, network.xl_rest_length, cfg.kappa_s_f
        )
        if network.xl_bonds.size
        else 0.0
    )
    return float(
        extensional_energy(positions, network.fiber_bonds, cfg.bead_spacing, cfg.kappa_s_f)
        + bending_energy(positions, network.bending_triples, cfg.theta0_f, cfg.kappa_b_f)
        + link_energy
    )


def nematic_order(angles, director):
    angles = np.asarray(angles, dtype=float)
    return float(np.mean(np.cos(2.0 * (angles - director)))) if angles.size else float("nan")


def polar_order(angles, director):
    angles = np.asarray(angles, dtype=float)
    return float(np.mean(np.cos(angles - director))) if angles.size else float("nan")
