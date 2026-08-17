"""Small deterministic collagen fixtures for G3 mechanism isolation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ecm.network_init import ECMNetwork
from .config import G3Config


@dataclass
class G3Fixture:
    name: str
    network: ECMNetwork
    initial_positions: np.ndarray
    fixed_mask: np.ndarray
    director_angle: float | None

    @property
    def segment_fiber_id(self) -> np.ndarray:
        if self.network.fiber_bonds.size == 0:
            return np.empty(0, dtype=int)
        return self.network.fiber_id[self.network.fiber_bonds[:, 0]]


def _network_from_rays(starts, directions, cfg: G3Config, name: str, director=None):
    starts = np.asarray(starts, dtype=float).reshape((-1, 2))
    directions = np.asarray(directions, dtype=float).reshape((-1, 2))
    n_fibres = starts.shape[0]
    n_beads = cfg.beads_per_fibre
    if n_fibres == 0:
        positions = np.empty((0, 2), dtype=float)
        fiber_id = np.empty(0, dtype=int)
        bonds = np.empty((0, 2), dtype=int)
        triples = np.empty((0, 3), dtype=int)
        fixed = np.empty(0, dtype=bool)
    else:
        directions = directions / np.linalg.norm(directions, axis=1)[:, None]
        arclength = np.arange(n_beads, dtype=float) * cfg.bead_spacing
        positions = (starts[:, None, :] + arclength[None, :, None] * directions[:, None, :])
        positions = positions.reshape((-1, 2))
        fiber_id = np.repeat(np.arange(n_fibres), n_beads)
        base = (np.arange(n_fibres) * n_beads)[:, None]
        edge = np.arange(n_beads - 1)[None, :]
        first = (base + edge).ravel()
        bonds = np.column_stack((first, first + 1)).astype(int)
        bend = np.arange(n_beads - 2)[None, :]
        first_bend = (base + bend).ravel()
        triples = np.column_stack((first_bend, first_bend + 1, first_bend + 2)).astype(int)
        fixed = np.zeros(positions.shape[0], dtype=bool)
        fixed[np.arange(n_fibres) * n_beads] = True
        fixed[np.arange(n_fibres) * n_beads + (n_beads - 1)] = True

    network = ECMNetwork(
        positions=positions.copy(),
        fiber_bonds=bonds,
        bending_triples=triples,
        xl_bonds=np.empty((0, 2), dtype=int),
        xl_rest_length=np.empty(0, dtype=float),
        fiber_id=fiber_id,
    )
    return G3Fixture(name, network, positions.copy(), fixed, director)


def _surface_start(phi, radius):
    return radius * np.column_stack((np.cos(phi), np.sin(phi)))


def _rotate(points, angle):
    c, s = np.cos(angle), np.sin(angle)
    rotation = np.array([[c, -s], [s, c]])
    return np.asarray(points) @ rotation.T


def build_fixture(name: str, cfg: G3Config, rng: np.random.Generator) -> G3Fixture:
    """Build one of the preregistered few-fibre geometries.

    Fibres begin 0.5 um outside the rigid cell and extend outward. Both endpoints are fixed,
    eliminating the rigid-rotation zero mode while representing continuation into a larger
    matrix without introducing an annular wall.
    """
    contact_radius = cfg.cell_radius + 0.5e-6

    if name == "empty":
        return _network_from_rays([], [], cfg, name)

    if name == "single_fibre":
        # A tangential fibre passes the cell at its midpoint. Pulling that midpoint inward
        # produces a visible V-shaped radial recruitment while both far endpoints stay fixed.
        starts = np.array([[contact_radius, -0.5 * cfg.fibre_length]])
        directions = np.array([[0.0, 1.0]])
        return _network_from_rays(starts, directions, cfg, name, director=0.0)

    if name == "balanced_8":
        phi = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False)
        starts = _surface_start(phi, contact_radius)
        directions = np.column_stack((np.cos(phi), np.sin(phi)))
        return _network_from_rays(starts, directions, cfg, name)

    if name in ("aligned_8", "aligned_8_rotated_30"):
        offsets = np.radians(np.array([-12.0, -4.0, 4.0, 12.0]))
        phi = np.concatenate((offsets, np.pi + offsets))
        starts = _surface_start(phi, contact_radius)
        directions = np.vstack((np.tile([1.0, 0.0], (4, 1)), np.tile([-1.0, 0.0], (4, 1))))
        angle = np.radians(30.0) if name.endswith("rotated_30") else 0.0
        if angle:
            starts = _rotate(starts, angle)
            directions = _rotate(directions, angle)
        return _network_from_rays(starts, directions, cfg, name, director=angle)

    if name == "isotropic_random_8":
        phi = rng.uniform(0.0, 2.0 * np.pi, cfg.fibre_count)
        psi = rng.uniform(0.0, 2.0 * np.pi, cfg.fibre_count)
        radial = np.column_stack((np.cos(phi), np.sin(phi)))
        directions = np.column_stack((np.cos(psi), np.sin(psi)))
        inward = np.sum(directions * radial, axis=1) < 0.0
        directions[inward] *= -1.0
        starts = contact_radius * radial
        return _network_from_rays(starts, directions, cfg, name)

    if name == "asymmetric_torque":
        phi = np.radians(np.array([18.0, 34.0, 52.0]))
        psi = np.radians(np.array([-4.0, 12.0, 28.0]))
        starts = _surface_start(phi, contact_radius)
        directions = np.column_stack((np.cos(psi), np.sin(psi)))
        return _network_from_rays(starts, directions, cfg, name)

    raise ValueError(
        f"unknown fixture {name!r}; expected empty, single_fibre, balanced_8, "
        "isotropic_random_8, aligned_8, aligned_8_rotated_30, or asymmetric_torque"
    )


def rotated_fixture(fixture: G3Fixture, angle: float, name: str | None = None) -> G3Fixture:
    """Rigidly rotate a fixture for covariance tests."""
    positions = _rotate(fixture.network.positions, angle)
    network = ECMNetwork(
        positions=positions.copy(),
        fiber_bonds=fixture.network.fiber_bonds.copy(),
        bending_triples=fixture.network.bending_triples.copy(),
        xl_bonds=fixture.network.xl_bonds.copy(),
        xl_rest_length=fixture.network.xl_rest_length.copy(),
        fiber_id=fixture.network.fiber_id.copy(),
    )
    director = None if fixture.director_angle is None else fixture.director_angle + angle
    return G3Fixture(name or f"{fixture.name}_rotated", network, positions.copy(),
                     fixture.fixed_mask.copy(), director)


def mirrored_fixture(fixture: G3Fixture, axis: str = "x") -> G3Fixture:
    """Mirror a fixture for force/torque sign controls."""
    positions = fixture.network.positions.copy()
    coord = 1 if axis == "x" else 0
    positions[:, coord] *= -1.0
    network = ECMNetwork(
        positions=positions.copy(),
        fiber_bonds=fixture.network.fiber_bonds.copy(),
        bending_triples=fixture.network.bending_triples.copy(),
        xl_bonds=fixture.network.xl_bonds.copy(),
        xl_rest_length=fixture.network.xl_rest_length.copy(),
        fiber_id=fixture.network.fiber_id.copy(),
    )
    director = fixture.director_angle
    if director is not None:
        director = -director if axis == "x" else np.pi - director
    return G3Fixture(f"{fixture.name}_mirror_{axis}", network, positions.copy(),
                     fixture.fixed_mask.copy(), director)
