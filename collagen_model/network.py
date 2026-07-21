"""Diluted triangular bead network with SLS axial segments and elastic bending."""

from dataclasses import dataclass
from math import exp, sqrt
from typing import Dict, List, Tuple

import numpy as np


@dataclass
class CollagenNetwork:
    positions: np.ndarray
    edges: np.ndarray
    bending_triplets: np.ndarray
    fixed: np.ndarray
    k0: float
    kinf: float
    tau: float
    bend_stiffness: float
    drag: float

    def __post_init__(self) -> None:
        self.positions = np.asarray(self.positions, dtype=float)
        self.edges = np.asarray(self.edges, dtype=int)
        self.bending_triplets = np.asarray(self.bending_triplets, dtype=int).reshape((-1, 3))
        self.fixed = np.asarray(self.fixed, dtype=bool)
        if self.positions.ndim != 2 or self.positions.shape[1] != 2:
            raise ValueError("positions must have shape (n_beads, 2)")
        if self.edges.ndim != 2 or self.edges.shape[1] != 2:
            raise ValueError("edges must have shape (n_edges, 2)")
        if self.k0 <= 0 or not 0 <= self.kinf <= self.k0 or self.tau <= 0 or self.drag <= 0:
            raise ValueError("invalid SLS or drag parameters")
        edge_vectors = self.positions[self.edges[:, 1]] - self.positions[self.edges[:, 0]]
        self.rest_lengths = np.linalg.norm(edge_vectors, axis=1)
        self.previous_lengths = self.rest_lengths.copy()
        self.maxwell_force = np.zeros(len(self.edges), dtype=float)
        self.reference_positions = self.positions.copy()

    @property
    def n_beads(self) -> int:
        return len(self.positions)

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    @property
    def bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        return self.reference_positions.min(axis=0), self.reference_positions.max(axis=0)

    def instantaneous_sync(self) -> None:
        """Treat the current position change as an ideal instantaneous strain step."""

        vectors = self.positions[self.edges[:, 1]] - self.positions[self.edges[:, 0]]
        lengths = np.linalg.norm(vectors, axis=1)
        self.maxwell_force += (self.k0 - self.kinf) * (lengths - self.previous_lengths)
        self.previous_lengths = lengths

    def internal_forces(self, dt: float, zero_fixed: bool = True) -> np.ndarray:
        if dt <= 0:
            raise ValueError("dt must be positive")
        forces = np.zeros_like(self.positions)
        i = self.edges[:, 0]
        j = self.edges[:, 1]
        vectors = self.positions[j] - self.positions[i]
        lengths = np.linalg.norm(vectors, axis=1)
        safe_lengths = np.maximum(lengths, 1.0e-12)
        unit = vectors / safe_lengths[:, None]

        km = self.k0 - self.kinf
        decay = exp(-dt / self.tau)
        rates = (lengths - self.previous_lengths) / dt
        self.maxwell_force = (
            decay * self.maxwell_force + km * self.tau * (1.0 - decay) * rates
        )
        extension = lengths - self.rest_lengths
        axial = self.kinf * extension + self.maxwell_force
        pair_force = axial[:, None] * unit
        np.add.at(forces, i, pair_force)
        np.add.at(forces, j, -pair_force)
        self.previous_lengths = lengths

        # Small-curvature discretization of filament bending. Crosslinks are
        # freely hinged because triplets only join collinear lattice segments.
        if len(self.bending_triplets):
            a = self.bending_triplets[:, 0]
            b = self.bending_triplets[:, 1]
            c = self.bending_triplets[:, 2]
            curvature = self.positions[a] - 2.0 * self.positions[b] + self.positions[c]
            bend_force = self.bend_stiffness * curvature
            np.add.at(forces, a, -bend_force)
            np.add.at(forces, b, 2.0 * bend_force)
            np.add.at(forces, c, -bend_force)

        if zero_fixed:
            forces[self.fixed] = 0.0
        return forces

    def advance(self, external_forces: np.ndarray, dt: float) -> np.ndarray:
        external_forces = np.asarray(external_forces, dtype=float)
        if external_forces.shape != self.positions.shape:
            raise ValueError("external_forces must match positions")
        total = self.internal_forces(dt) + external_forces
        total[self.fixed] = 0.0
        velocity = total / self.drag
        self.positions[~self.fixed] += dt * velocity[~self.fixed]
        self.positions[self.fixed] = self.reference_positions[self.fixed]
        return velocity

    def axial_force_magnitudes(self) -> np.ndarray:
        vectors = self.positions[self.edges[:, 1]] - self.positions[self.edges[:, 0]]
        lengths = np.linalg.norm(vectors, axis=1)
        return np.abs(self.kinf * (lengths - self.rest_lengths) + self.maxwell_force)


def make_diluted_triangular_network(
    nx: int = 16,
    ny: int = 12,
    spacing: float = 2.5,
    bond_dilution: float = 0.15,
    seed: int = 1,
    k0: float = 18.0,
    kinf: float = 4.5,
    tau: float = 5.0,
    bend_stiffness: float = 0.35,
    drag: float = 45.0,
) -> CollagenNetwork:
    """Build a finite 2D triangular fiber network.

    Three families of straight lattice lines act as fibers. Shared vertices are
    freely hinged crosslinks; bending triplets are added only along one line
    family. Perimeter beads are fixed as a finite-domain boundary condition.
    """

    if nx < 4 or ny < 4:
        raise ValueError("nx and ny must be at least 4")
    if not 0.0 <= bond_dilution < 1.0:
        raise ValueError("bond_dilution must be in [0, 1)")
    rng = np.random.default_rng(seed)

    def index(ii: int, jj: int) -> int:
        return jj * nx + ii

    positions = np.zeros((nx * ny, 2), dtype=float)
    for jj in range(ny):
        for ii in range(nx):
            positions[index(ii, jj)] = (
                spacing * (ii + 0.5 * jj),
                spacing * (sqrt(3.0) / 2.0) * jj,
            )

    # direction 0: (1, 0), direction 1: (0, 1), direction 2: (-1, 1)
    candidate_edges: List[Tuple[int, int, int]] = []
    for jj in range(ny):
        for ii in range(nx):
            if ii + 1 < nx:
                candidate_edges.append((index(ii, jj), index(ii + 1, jj), 0))
            if jj + 1 < ny:
                candidate_edges.append((index(ii, jj), index(ii, jj + 1), 1))
            if ii - 1 >= 0 and jj + 1 < ny:
                candidate_edges.append((index(ii, jj), index(ii - 1, jj + 1), 2))

    kept: List[Tuple[int, int, int]] = []
    for edge in candidate_edges:
        if rng.random() >= bond_dilution:
            kept.append(edge)

    # Always keep boundary-parallel edges to provide a robust anchored frame.
    for edge in candidate_edges:
        a, b, direction = edge
        ia, ja = a % nx, a // nx
        ib, jb = b % nx, b // nx
        if ia in (0, nx - 1) or ib in (0, nx - 1) or ja in (0, ny - 1) or jb in (0, ny - 1):
            if edge not in kept:
                kept.append(edge)

    edges = np.asarray([(a, b) for a, b, _ in kept], dtype=int)
    edge_lookup: Dict[Tuple[int, int, int], bool] = {}
    for a, b, direction in kept:
        edge_lookup[(min(a, b), max(a, b), direction)] = True

    triplets: List[Tuple[int, int, int]] = []
    direction_steps = [(1, 0), (0, 1), (-1, 1)]
    for jj in range(ny):
        for ii in range(nx):
            center = index(ii, jj)
            for direction, (di, dj) in enumerate(direction_steps):
                pi, pj = ii - di, jj - dj
                ni, nj = ii + di, jj + dj
                if 0 <= pi < nx and 0 <= pj < ny and 0 <= ni < nx and 0 <= nj < ny:
                    prev_node = index(pi, pj)
                    next_node = index(ni, nj)
                    left_key = (min(prev_node, center), max(prev_node, center), direction)
                    right_key = (min(center, next_node), max(center, next_node), direction)
                    if left_key in edge_lookup and right_key in edge_lookup:
                        triplets.append((prev_node, center, next_node))

    fixed = np.zeros(nx * ny, dtype=bool)
    for jj in range(ny):
        for ii in range(nx):
            if ii in (0, nx - 1) or jj in (0, ny - 1):
                fixed[index(ii, jj)] = True

    return CollagenNetwork(
        positions=positions,
        edges=edges,
        bending_triplets=np.asarray(triplets, dtype=int),
        fixed=fixed,
        k0=k0,
        kinf=kinf,
        tau=tau,
        bend_stiffness=bend_stiffness,
        drag=drag,
    )
