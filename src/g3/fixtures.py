"""Small deterministic collagen fixtures for G3 mechanism isolation."""

from __future__ import annotations

from dataclasses import dataclass, field
import math

import numpy as np

from .config import G3Config
from .elastic import ECMNetwork

from generations.g2_corrected.common.model import (
    CollagenConfig as G2CollagenConfig,
    Network as G2Network,
    NetworkSpec as G2NetworkSpec,
    connectivity_report as g2_connectivity_report,
)


@dataclass
class G3Fixture:
    name: str
    network: ECMNetwork
    initial_positions: np.ndarray
    fixed_mask: np.ndarray
    director_angle: float | None
    _segment_fiber_id: np.ndarray = field(init=False, repr=False)
    bead_ids_by_fibre: tuple[np.ndarray, ...] = field(init=False, repr=False)

    def __post_init__(self):
        bonds = self.network.fiber_bonds
        self._segment_fiber_id = (
            self.network.fiber_id[bonds[:, 0]] if bonds.size else np.empty(0, dtype=int)
        )
        self.bead_ids_by_fibre = tuple(
            np.flatnonzero(self.network.fiber_id == fibre)
            for fibre in range(self.network.n_fibers)
        )

    @property
    def segment_fiber_id(self) -> np.ndarray:
        return self._segment_fiber_id

    @property
    def has_uniform_bead_count(self) -> bool:
        sizes = [ids.size for ids in self.bead_ids_by_fibre]
        return len(set(sizes)) <= 1


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


def _straight_points(start, end, spacing):
    length = float(np.linalg.norm(end - start))
    n = max(4, int(math.ceil(length / spacing)) + 1)
    return np.linspace(start, end, n)


def _scaled_g2_config(cfg: G3Config) -> G2CollagenConfig:
    """Build the multi-fibre mechanics config in G2's nN--um unit system."""
    return G2CollagenConfig(
        domain_size=cfg.scaled_domain_size * 1.0e6,
        cell_radius=cfg.cell_radius * 1.0e6,
        n_fibers=cfg.scaled_fibre_count,
        min_fiber_length=cfg.scaled_min_fibre_length * 1.0e6,
        max_fiber_length=cfg.scaled_max_fibre_length * 1.0e6,
        bead_spacing=cfg.scaled_bead_spacing * 1.0e6,
        required_connected_fraction=cfg.scaled_required_connected_fraction,
    )


def _sample_scaled_spec(g2cfg, rng, director, seed_used):
    half = g2cfg.domain_size / 2.0
    positions = []
    fibres = []
    fixed = []

    def append(points):
        ids = []
        for point in points:
            ids.append(len(positions))
            positions.append(np.asarray(point, dtype=float))
            fixed.append(bool(np.max(np.abs(point)) >= half - g2cfg.boundary_width))
        fibres.append(ids)

    # Eight near-cell fibres guarantee that every stage has multiple physical
    # candidates. In aligned fixtures they form two nematically equivalent poles.
    if director is None:
        phi = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False)
        directions = np.column_stack((np.cos(phi), np.sin(phi)))
    else:
        offsets = np.radians(np.array([-12.0, -4.0, 4.0, 12.0]))
        phi = np.concatenate((director + offsets, director + np.pi + offsets))
        directions = np.vstack((
            np.tile([np.cos(director), np.sin(director)], (4, 1)),
            np.tile([-np.cos(director), -np.sin(director)], (4, 1)),
        ))
    # Start beyond the clutch capture radius so the movie resolves a real
    # grow-before-contact interval instead of binding at frame zero.
    starts = (g2cfg.cell_radius + 2.5) * np.column_stack((np.cos(phi), np.sin(phi)))
    for start, direction in zip(starts, directions):
        # Reach the outer anchor band while remaining outside the cell.
        candidates = []
        for axis in range(2):
            if abs(direction[axis]) > 1.0e-9:
                target = (math.copysign(half - 0.25, direction[axis]) - start[axis]) / direction[axis]
                if target > 0.0:
                    candidates.append(target)
        length = min(candidates) if candidates else g2cfg.max_fiber_length
        length = min(length, g2cfg.max_fiber_length)
        append(_straight_points(start, start + length * direction, g2cfg.bead_spacing))
    contact_fibres = list(range(8))

    target_boundary = max(24, int(round(0.42 * g2cfg.n_fibers)))
    attempts = 0
    while len(fibres) < g2cfg.n_fibers and attempts < 100_000:
        attempts += 1
        length = rng.uniform(g2cfg.min_fiber_length, g2cfg.max_fiber_length)
        if director is not None and rng.random() < 0.72:
            theta = director + rng.normal(0.0, 0.30)
        else:
            theta = rng.uniform(0.0, np.pi)
        direction = np.array([np.cos(theta), np.sin(theta)])
        boundary_seeded = len(fibres) < target_boundary
        if boundary_seeded:
            side = int(rng.integers(4))
            along = rng.uniform(-0.85 * half, 0.85 * half)
            if side == 0:
                start = np.array([-half + 0.25, along])
            elif side == 1:
                start = np.array([half - 0.25, along])
            elif side == 2:
                start = np.array([along, -half + 0.25])
            else:
                start = np.array([along, half - 0.25])
            if np.dot(direction, -start) < 0.0:
                direction *= -1.0
            end = start + length * direction
        else:
            center = rng.uniform(-0.58 * half, 0.58 * half, size=2)
            start = center - 0.5 * length * direction
            end = center + 0.5 * length * direction
        if np.max(np.abs(start)) > half - 0.15 or np.max(np.abs(end)) > half - 0.15:
            continue
        points = _straight_points(start, end, g2cfg.bead_spacing)
        # Reserve the cell-scale probing shell for the eight declared contact
        # fibres. Background fibres can receive force through crosslinks but
        # cannot become unregistered direct adhesion targets.
        if np.min(np.linalg.norm(points, axis=1)) < g2cfg.cell_radius + 8.5:
            continue
        append(points)
    if len(fibres) != g2cfg.n_fibers:
        raise RuntimeError("could not construct the scaled finite-fibre network")
    return G2NetworkSpec(
        np.asarray(positions), fibres, np.asarray(fixed), contact_fibres, int(seed_used)
    )


def _scaled_fixture(name: str, cfg: G3Config, rng: np.random.Generator, director=None):
    """Generate a 99-fibre boundary-connected network with G2 mechanics."""
    g2cfg = _scaled_g2_config(cfg)
    base_seed = int(rng.integers(0, 2**31 - 1))
    best = None
    best_score = -1.0
    for attempt in range(80):
        seed_used = base_seed + 7919 * attempt
        spec = _sample_scaled_spec(g2cfg, np.random.default_rng(seed_used), director, seed_used)
        external = G2Network(spec, g2cfg)
        report = g2_connectivity_report(external)
        score = float(report["connected_fraction"])
        if score > best_score:
            best = (external, report)
            best_score = score
        if report["contact_fibers_connected"] and score >= cfg.scaled_required_connected_fraction:
            break
    else:
        external, report = best
        if not report["contact_fibers_connected"] or best_score < 0.70:
            raise RuntimeError(
                f"scaled network percolation failed; best connected fraction={best_score:.3f}"
            )

    positions = external.r * 1.0e-6
    fibre_id = np.empty(positions.shape[0], dtype=int)
    for fid, ids in enumerate(external.fibers):
        fibre_id[ids] = fid
    if external.crosslinks:
        edge_a = external.link_edge_a.copy()
        edge_b = external.link_edge_b.copy()
        alpha_a = external.link_alpha_a.copy()
        alpha_b = external.link_alpha_b.copy()
        bead_a = external.edges[edge_a, (alpha_a >= 0.5).astype(int)]
        bead_b = external.edges[edge_b, (alpha_b >= 0.5).astype(int)]
        xl_bonds = np.column_stack((bead_a, bead_b)).astype(int)
        xl_rest = np.linalg.norm(positions[bead_b] - positions[bead_a], axis=1)
    else:
        edge_a = edge_b = np.empty(0, dtype=int)
        alpha_a = alpha_b = np.empty(0, dtype=float)
        xl_bonds = np.empty((0, 2), dtype=int)
        xl_rest = np.empty(0, dtype=float)
    network = ECMNetwork(
        positions=positions.copy(),
        fiber_bonds=external.edges.copy(),
        bending_triples=external.triplets.copy(),
        xl_bonds=xl_bonds,
        xl_rest_length=xl_rest,
        fiber_id=fibre_id,
        external_network=external,
        effective_bead_drag=g2cfg.bead_drag * 1.0e-3,
        xl_edge_a=edge_a,
        xl_alpha_a=alpha_a,
        xl_edge_b=edge_b,
        xl_alpha_b=alpha_b,
    )
    fixture = G3Fixture(name, network, positions.copy(), external.fixed.copy(), director)
    fixture.connectivity_report = report
    return fixture


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

    if name == "scaled_isotropic_99":
        return _scaled_fixture(name, cfg, rng, director=None)

    if name == "scaled_aligned_99":
        return _scaled_fixture(name, cfg, rng, director=0.0)

    if name == "scaled_aligned_99_rotated_30":
        return _scaled_fixture(name, cfg, rng, director=np.radians(30.0))

    raise ValueError(
        f"unknown fixture {name!r}; expected empty, single_fibre, balanced_8, "
        "isotropic_random_8, aligned_8, aligned_8_rotated_30, asymmetric_torque, "
        "scaled_isotropic_99, scaled_aligned_99, or scaled_aligned_99_rotated_30"
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
