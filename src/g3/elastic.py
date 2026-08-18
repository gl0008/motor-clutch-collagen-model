"""Minimal SI elastic bead-chain core used only by Generation 3 fixtures."""

from __future__ import annotations

from dataclasses import dataclass
import os

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - exercised only in minimal installations
    njit = None


NUMBA_ENABLED = njit is not None and os.environ.get("G3_DISABLE_NUMBA") != "1"


@dataclass
class ECMNetwork:
    """Topology for the inspectable G3 few-fibre network."""

    positions: np.ndarray
    fiber_bonds: np.ndarray
    bending_triples: np.ndarray
    xl_bonds: np.ndarray
    xl_rest_length: np.ndarray
    fiber_id: np.ndarray
    # Optional bridge to the validated G2 network mechanics. Revised 99-fibre
    # fixtures use its variable bond lengths, reference curvature, compression
    # softening and material-point crosslinks, converted to SI at the boundary.
    external_network: object | None = None
    effective_bead_drag: float | None = None
    xl_edge_a: np.ndarray | None = None
    xl_alpha_a: np.ndarray | None = None
    xl_edge_b: np.ndarray | None = None
    xl_alpha_b: np.ndarray | None = None

    @property
    def n_beads(self) -> int:
        return int(self.positions.shape[0])

    @property
    def n_fibers(self) -> int:
        return int(self.fiber_id.max()) + 1 if self.fiber_id.size else 0


def harmonic_bond_energy(positions, bonds, rest_length, stiffness) -> float:
    """Hookean bond energy, 0.5*k*(length-rest_length)^2."""
    if bonds.size == 0:
        return 0.0
    delta = positions[bonds[:, 1]] - positions[bonds[:, 0]]
    length = np.linalg.norm(delta, axis=1)
    return float(0.5 * np.sum(stiffness * (length - rest_length) ** 2))


def harmonic_bond_forces_numpy(positions, bonds, rest_length, stiffness):
    """Equal-and-opposite Hookean forces for an arbitrary bond list."""
    n_beads = positions.shape[0]
    if bonds.size == 0:
        return np.zeros_like(positions)
    first = bonds[:, 0]
    second = bonds[:, 1]
    delta = positions[second] - positions[first]
    length = np.linalg.norm(delta, axis=1)
    direction = np.divide(
        delta,
        length[:, None],
        out=np.zeros_like(delta),
        where=length[:, None] > 0.0,
    )
    vector = (stiffness * (length - rest_length))[:, None] * direction
    fx = (np.bincount(first, vector[:, 0], minlength=n_beads)
          - np.bincount(second, vector[:, 0], minlength=n_beads))
    fy = (np.bincount(first, vector[:, 1], minlength=n_beads)
          - np.bincount(second, vector[:, 1], minlength=n_beads))
    return np.column_stack((fx, fy))


if njit is not None:
    @njit(cache=True)
    def _harmonic_bond_forces_numba(positions, bonds, rest_length, stiffness):
        forces = np.zeros_like(positions)
        for bond_id in range(bonds.shape[0]):
            first = bonds[bond_id, 0]
            second = bonds[bond_id, 1]
            dx = positions[second, 0] - positions[first, 0]
            dy = positions[second, 1] - positions[first, 1]
            length = np.sqrt(dx * dx + dy * dy)
            if length <= 0.0:
                continue
            magnitude = stiffness * (length - rest_length) / length
            fx = magnitude * dx
            fy = magnitude * dy
            forces[first, 0] += fx
            forces[first, 1] += fy
            forces[second, 0] -= fx
            forces[second, 1] -= fy
        return forces


def harmonic_bond_forces(positions, bonds, rest_length, stiffness):
    """Equal-and-opposite Hookean forces, optionally using the verified Numba kernel."""
    if NUMBA_ENABLED:
        return _harmonic_bond_forces_numba(positions, bonds, rest_length, stiffness)
    return harmonic_bond_forces_numpy(positions, bonds, rest_length, stiffness)


def extensional_energy(positions, bonds, rest_length, stiffness) -> float:
    """Fibre extensional energy; Saraswathibhatla et al. 2025 Eq. 7."""
    return harmonic_bond_energy(positions, bonds, rest_length, stiffness)


def extensional_forces(positions, bonds, rest_length, stiffness):
    """Fibre extensional forces; Saraswathibhatla et al. 2025 Eq. 7."""
    return harmonic_bond_forces(positions, bonds, rest_length, stiffness)


def _angles(positions, triples):
    first, middle, last = triples[:, 0], triples[:, 1], triples[:, 2]
    u = positions[middle] - positions[first]
    v = positions[last] - positions[middle]
    lu = np.linalg.norm(u, axis=1)
    lv = np.linalg.norm(v, axis=1)
    cosine = np.divide(
        np.sum(u * v, axis=1),
        lu * lv,
        out=np.ones_like(lu),
        where=(lu * lv) > 0.0,
    )
    return u, v, lu, lv, np.clip(cosine, -1.0, 1.0)


def bending_energy(positions, triples, theta0, stiffness) -> float:
    """Harmonic angle energy; Saraswathibhatla et al. 2025 Eq. 8."""
    if triples.size == 0:
        return 0.0
    _, _, _, _, cosine = _angles(positions, triples)
    theta = np.arccos(cosine)
    return float(0.5 * np.sum(stiffness * (theta - theta0) ** 2))


def bending_forces_numpy(positions, triples, theta0, stiffness):
    """Analytical negative gradient of the harmonic angle energy."""
    n_beads = positions.shape[0]
    if triples.size == 0:
        return np.zeros_like(positions)
    first, middle, last = triples[:, 0], triples[:, 1], triples[:, 2]
    u, v, lu, lv, cosine = _angles(positions, triples)
    theta = np.arccos(cosine)
    sine = np.sqrt(np.maximum(1.0 - cosine**2, 0.0))
    uv = np.sum(u * v, axis=1)
    v_perp_u = v - np.divide(
        uv, lu**2, out=np.zeros_like(uv), where=lu > 0.0
    )[:, None] * u
    u_perp_v = u - np.divide(
        uv, lv**2, out=np.zeros_like(uv), where=lv > 0.0
    )[:, None] * v
    straight = sine <= 1.0e-12
    sine_safe = np.where(straight, 1.0, sine)
    prefactor = np.where(straight, 0.0, stiffness * (theta - theta0) / sine_safe)
    inverse_lengths = np.divide(
        1.0, lu * lv, out=np.zeros_like(lu), where=(lu * lv) > 0.0
    )
    force_first = -(prefactor * inverse_lengths)[:, None] * v_perp_u
    force_last = (prefactor * inverse_lengths)[:, None] * u_perp_v
    force_middle = -(force_first + force_last)

    fx = (
        np.bincount(first, force_first[:, 0], minlength=n_beads)
        + np.bincount(middle, force_middle[:, 0], minlength=n_beads)
        + np.bincount(last, force_last[:, 0], minlength=n_beads)
    )
    fy = (
        np.bincount(first, force_first[:, 1], minlength=n_beads)
        + np.bincount(middle, force_middle[:, 1], minlength=n_beads)
        + np.bincount(last, force_last[:, 1], minlength=n_beads)
    )
    return np.column_stack((fx, fy))


if njit is not None:
    @njit(cache=True)
    def _bending_forces_numba(positions, triples, theta0, stiffness):
        forces = np.zeros_like(positions)
        for triple_id in range(triples.shape[0]):
            first = triples[triple_id, 0]
            middle = triples[triple_id, 1]
            last = triples[triple_id, 2]
            ux = positions[middle, 0] - positions[first, 0]
            uy = positions[middle, 1] - positions[first, 1]
            vx = positions[last, 0] - positions[middle, 0]
            vy = positions[last, 1] - positions[middle, 1]
            lu2 = ux * ux + uy * uy
            lv2 = vx * vx + vy * vy
            if lu2 <= 0.0 or lv2 <= 0.0:
                continue
            lu = np.sqrt(lu2)
            lv = np.sqrt(lv2)
            uv = ux * vx + uy * vy
            cosine = uv / (lu * lv)
            cosine = min(1.0, max(-1.0, cosine))
            sine = np.sqrt(max(1.0 - cosine * cosine, 0.0))
            if sine <= 1.0e-12:
                continue
            theta = np.arccos(cosine)
            prefactor = stiffness * (theta - theta0) / sine
            inverse_lengths = 1.0 / (lu * lv)
            v_perp_ux = vx - (uv / lu2) * ux
            v_perp_uy = vy - (uv / lu2) * uy
            u_perp_vx = ux - (uv / lv2) * vx
            u_perp_vy = uy - (uv / lv2) * vy
            first_x = -prefactor * inverse_lengths * v_perp_ux
            first_y = -prefactor * inverse_lengths * v_perp_uy
            last_x = prefactor * inverse_lengths * u_perp_vx
            last_y = prefactor * inverse_lengths * u_perp_vy
            middle_x = -(first_x + last_x)
            middle_y = -(first_y + last_y)
            forces[first, 0] += first_x
            forces[first, 1] += first_y
            forces[middle, 0] += middle_x
            forces[middle, 1] += middle_y
            forces[last, 0] += last_x
            forces[last, 1] += last_y
        return forces


def bending_forces(positions, triples, theta0, stiffness):
    """Analytical bending forces, optionally using the verified Numba kernel."""
    if NUMBA_ENABLED:
        return _bending_forces_numba(positions, triples, theta0, stiffness)
    return bending_forces_numpy(positions, triples, theta0, stiffness)


def overdamped_step(positions, forces, drag, dt):
    """Deterministic overdamped Euler step: r_next = r + F*dt/drag."""
    drag = np.asarray(drag, dtype=float)
    if drag.ndim == 1:
        drag = drag[:, None]
    return positions + forces * dt / drag
