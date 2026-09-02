"""Generation 5 - multicellular tumor organoid on the G2 collagen engine.

Units: micrometre, nanonewton, second -- identical to the G2 / G3 / G4 engine so
every generation shares one coordinate and force language.

Architecture (docs/G5_organoid_plan.md):

  organoid = N motor-clutch DISKS  +  simplified cell-cell adhesion
             each disk grips the shared G2 bead-spring collagen network via the
             existing Gaussian contact kernel and reels fibres inward (contraction)

Overdamped force balance on every collagen bead (unchanged from G2)::

    zeta * dr_i/dt = F_stretch + F_bend + F_crosslink + F_repulsion + F_active

where ``F_repulsion`` now sums soft no-penetration over ALL cell disks and
``F_active`` sums the inward Gaussian-projected pull of every surface cell.

This module implements Stage A (multi-cell scaffold + resting balance) and
Stage B (fixed contractile organoid -> radial collagen alignment).  Cells do not
translate yet; releasing them is Stage D.  Matrix plasticity is Stage E.

Citations
---------
* G2 engine, Gaussian contact kernel, bead-spring energies: ``common/model.py``
  (Saraswathibhatla 2025 SI Table 2 parameters; Lee 2014 fibre ranges).
* Cell-cell adhesion / tissue surface tension as the collective<->single-cell
  determinant: Ilina & Friedl 2020 Nat Cell Biol (qualitative; magnitudes fit).
* Force-induced radial (TACS-3-like) alignment as the Stage-B target:
  Saraswathibhatla 2025 bioRxiv; Su/Kim 2021 Biomaterials.
* Per-cell contractile pull magnitude (tens of nN): Steinwachs 2016; Mark 2020.
"""

from __future__ import annotations

from dataclasses import dataclass, replace, asdict
from pathlib import Path
import math
import sys

import numpy as np

# --- reuse the frozen Generation-2 collagen engine --------------------------------
HERE = Path(__file__).resolve().parent
G2_ROOT = HERE.parent / "g2_corrected"
if str(G2_ROOT) not in sys.path:
    sys.path.insert(0, str(G2_ROOT))

from common.model import (  # noqa: E402
    CollagenConfig,
    Network,
    NetworkSpec,
    connectivity_report,
    contact_patches,
    forces_from_patches,
    _curved_segment,
)


# =================================================================================
# Configuration
# =================================================================================
@dataclass(frozen=True)
class OrganoidConfig(CollagenConfig):
    """Organoid geometry, cell-cell mechanics, and the shared collagen network.

    Inherits every G2 collagen parameter (fibre stiffness, bending, crosslink
    stiffness, Gaussian sigma, bead drag, ...).  ``cell_radius`` here is the
    radius of a SINGLE cell disk, not the whole organoid.
    """

    # --- domain (large enough for the ~organoid + ~200 um remodelling halo) ---
    domain_size: float = 480.0

    # --- single cell ---
    cell_radius: float = 9.0             # one tumour cell (movies: ~15-20 um dia)

    # --- organoid packing ---
    organoid_radius: float = 63.0        # target radius for cell CENTRES (core)
    cell_spacing: float = 18.0           # hex lattice pitch (= 2*cell_radius touch)
    gap: float = 1.5                     # fibre-free ring beyond the organoid at t=0
    n_corona_fibers: int = 60            # short fibres hugging the organoid (grippable collagen)
    corona_band: float = 12.0            # radial thickness of the near-field corona

    # --- simplified cell-cell interaction (piecewise-linear soft adhesive disk) ---
    # Equilibrium centre separation is ``cell_spacing``; below it disks repel,
    # above it they adhere out to ``cc_adhesion_range``, then feel nothing.
    cc_repulsion: float = 40.0           # nN/um, overlap stiffness (keeps disks apart)
    cc_adhesion: float = 6.0             # nN/um, short-range cohesion (surface tension)
    cc_adhesion_range: float = 6.0       # um, how far past contact adhesion reaches

    # --- collagen network (scaled up from the 99-fibre / 180 um G2 baseline) ---
    # Softened, lightly crosslinked collagen so the pull can visibly reorganise the
    # near-field into a radial aster (same prof-requested softening as G3/G4:
    # 3 MPa modulus, 10 nN/um links).  Not the stiff 32 MPa G2 default.
    collagen_modulus_mpa: float = 3.0
    crosslink_stiffness: float = 10.0
    n_fibers: int = 260
    crosslink_fraction: float = 0.85     # keep most auto-generated links (connectivity)
    boundary_width: float = 5.0          # anchoring band; wider -> better boundary percolation
    contact_width: float = 8.0           # grip shell so surface cells reach the corona
    total_pull_force: float = 12.0       # nN, per cell (tens-of-nN traction scale)

    # --- integration / output ---
    dt: float = 0.02
    duration: float = 120.0
    sample_interval: float = 2.0
    contact_update_interval: float = 0.5
    force_ramp_time: float = 5.0
    required_connected_fraction: float = 0.55
    generation_attempts: int = 25
    seed: int = 23

    def validate(self) -> None:  # noqa: D401 - relax the single-cell G2 gate
        if self.domain_size < 2.2 * (self.organoid_radius + self.cell_radius):
            raise ValueError("domain must comfortably enclose the organoid + halo")
        if self.n_fibers < 16 or self.bead_spacing <= 0:
            raise ValueError("network resolution is too small")
        if self.dt <= 0 or self.sample_interval < self.dt:
            raise ValueError("invalid integration or output interval")

    @property
    def gap_radius_floor(self) -> float:
        """Lower bound for the fibre-free radius (updated once cells are packed)."""
        return self.organoid_radius + self.cell_radius + self.gap


# =================================================================================
# Cell packing
# =================================================================================
def hex_centers(radius: float, spacing: float) -> np.ndarray:
    """Hexagonal-lattice disk centres with ``|centre| <= radius``.

    A hex packing at pitch ``spacing`` puts every interior cell in a locally
    symmetric neighbourhood, so at pitch == cell_spacing the cell-cell forces
    cancel at rest (Stage-A balance).
    """

    dx = spacing
    dy = spacing * math.sqrt(3.0) / 2.0
    n = int(math.ceil(radius / min(dx, dy))) + 2
    pts: list[tuple[float, float]] = []
    for row in range(-n, n + 1):
        y = row * dy
        x_off = 0.5 * dx if (row % 2) else 0.0
        for col in range(-n, n + 1):
            x = col * dx + x_off
            if math.hypot(x, y) <= radius + 1e-9:
                pts.append((x, y))
    centers = np.asarray(sorted(pts, key=lambda p: (math.hypot(*p), math.atan2(p[1], p[0]))))
    return centers.reshape(-1, 2)


# =================================================================================
# Network generation with a fibre-free organoid gap (G3-style, larger gap)
# =================================================================================
def _organoid_fiber(cfg: OrganoidConfig, rng, gap_radius: float, *, boundary_seeded: bool):
    """One curved fibre that avoids the whole fibre-free organoid disk."""

    half = cfg.domain_size / 2.0
    length = rng.uniform(cfg.min_fiber_length, cfg.max_fiber_length)
    if boundary_seeded:
        side = int(rng.integers(4))
        along = rng.uniform(-0.82 * half, 0.82 * half)
        if side == 0:
            start = np.array([-half + 0.3, along]); inward = np.array([1.0, 0.0])
        elif side == 1:
            start = np.array([half - 0.3, along]); inward = np.array([-1.0, 0.0])
        elif side == 2:
            start = np.array([along, -half + 0.3]); inward = np.array([0.0, 1.0])
        else:
            start = np.array([along, half - 0.3]); inward = np.array([0.0, -1.0])
        theta = rng.uniform(-0.75, 0.75)
        rot = np.array([[math.cos(theta), -math.sin(theta)],
                        [math.sin(theta), math.cos(theta)]])
        end = start + length * (rot @ inward)
    else:
        center = rng.uniform(-0.7 * half, 0.7 * half, size=2)
        theta = rng.uniform(0.0, math.pi)
        direction = np.array([math.cos(theta), math.sin(theta)])
        start = center - 0.5 * length * direction
        end = center + 0.5 * length * direction
    if np.max(np.abs(end)) > half - 0.2 or np.max(np.abs(start)) > half - 0.2:
        return None
    points = _curved_segment(start, end, cfg, rng)
    if np.min(np.linalg.norm(points, axis=1)) < gap_radius:
        return None
    return points


def _corona_fiber(cfg: OrganoidConfig, rng, gap_radius: float):
    """A short fibre hugging the organoid, so surface cells have collagen to grip.

    The imaging shows collagen pressed right against the organoid boundary; the
    corona represents that grippable near-field matrix (it is reorganised into the
    radial aster).  Fibres are placed just outside the fibre-free gap with a mild
    tangential-to-random orientation and are NOT boundary-anchored.
    """

    phi = rng.uniform(0.0, 2.0 * math.pi)
    r0 = gap_radius + rng.uniform(0.5, cfg.corona_band)
    radial = np.array([math.cos(phi), math.sin(phi)])
    tangent = np.array([-radial[1], radial[0]])
    mid = r0 * radial
    length = rng.uniform(cfg.min_fiber_length * 0.6, cfg.min_fiber_length * 1.2)
    # orientation: biased tangential but noisy, so t=0 is not already radial
    theta = rng.uniform(0.0, math.pi)
    direction = math.cos(theta) * tangent + math.sin(theta) * radial * 0.6
    direction /= max(np.linalg.norm(direction), 1e-12)
    start = mid - 0.5 * length * direction
    end = mid + 0.5 * length * direction
    half = cfg.domain_size / 2.0
    if np.max(np.abs(end)) > half - 0.3 or np.max(np.abs(start)) > half - 0.3:
        return None
    points = _curved_segment(start, end, cfg, rng)
    if np.min(np.linalg.norm(points, axis=1)) < gap_radius:
        return None
    return points


def _assemble_organoid_spec(cfg: OrganoidConfig, rng, gap_radius: float, seed_used: int) -> NetworkSpec:
    half = cfg.domain_size / 2.0
    positions: list[np.ndarray] = []
    fibers: list[list[int]] = []
    fixed: list[bool] = []

    def append(points: np.ndarray) -> None:
        ids: list[int] = []
        for p in points:
            ids.append(len(positions))
            positions.append(np.asarray(p, dtype=float))
            fixed.append(bool(np.max(np.abs(p)) >= half - cfg.boundary_width))
        fibers.append(ids)

    # 1) near-field corona hugging the organoid (grippable collagen)
    corona = 0
    attempts = 0
    while corona < cfg.n_corona_fibers and attempts < 200_000:
        attempts += 1
        pts = _corona_fiber(cfg, rng, gap_radius)
        if pts is not None:
            append(pts)
            corona += 1
    contact_fibers = list(range(len(fibers)))  # corona fibres are the grippable set

    # 2) bulk network (boundary-seeded first for percolation, then interior)
    target_total = cfg.n_fibers
    target_boundary = corona + max(24, int(round(0.45 * (target_total - corona))))
    attempts = 0
    while len(fibers) < target_total and attempts < 400_000:
        attempts += 1
        pts = _organoid_fiber(
            cfg, rng, gap_radius, boundary_seeded=len(fibers) < target_boundary
        )
        if pts is not None:
            append(pts)
    if len(fibers) != target_total:
        raise RuntimeError("could not construct the requested fibre-free-gap network")
    return NetworkSpec(np.asarray(positions), fibers, np.asarray(fixed), contact_fibers, seed_used)


def make_organoid(cfg: OrganoidConfig = OrganoidConfig(), seed=None):
    """Pack the cells and build a crosslinked, boundary-anchored collagen network.

    Returns ``(network, centers, gap_radius, report)``.  ``gap_radius`` is derived
    from the actual packing so the surface cells' contact shell just reaches the
    innermost fibres.
    """

    cfg.validate()
    centers = hex_centers(cfg.organoid_radius, cfg.cell_spacing)
    organoid_outer = float(np.max(np.linalg.norm(centers, axis=1))) + cfg.cell_radius
    gap_radius = organoid_outer + cfg.gap

    base = cfg.seed if seed is None else int(seed)
    best = (-1.0, None)
    for attempt in range(cfg.generation_attempts):
        seed_used = base + 7919 * attempt
        spec = _assemble_organoid_spec(
            cfg, np.random.default_rng(seed_used), gap_radius, seed_used
        )
        network = Network(spec, cfg)  # crosslinks auto-built in __init__
        if cfg.crosslink_fraction < 1.0 and network.crosslinks:
            rng = np.random.default_rng(seed_used + 101)
            keep = rng.random(len(network.crosslinks)) < cfg.crosslink_fraction
            network.crosslinks = [x for x, k in zip(network.crosslinks, keep) if k]
            network.refresh_crosslink_arrays()
        report = connectivity_report(network)
        score = float(report["connected_fraction"])
        if score > best[0]:
            best = (score, network)
        if score >= cfg.required_connected_fraction:
            return network, centers, gap_radius, report
    if best[1] is not None:
        return best[1], centers, gap_radius, connectivity_report(best[1])
    raise RuntimeError("organoid network percolation gate failed")


# =================================================================================
# Simplified cell-cell interaction (Stage A scaffold; drives motion in Stage D)
# =================================================================================
def cell_cell_forces(centers: np.ndarray, cfg: OrganoidConfig) -> np.ndarray:
    """Pairwise soft adhesive-disk force on each cell centre.

    Piecewise-linear in the centre separation ``d`` with equilibrium at
    ``cell_spacing`` (``d_eq``)::

        d < d_eq                     : repel,  F = cc_repulsion * (d_eq - d)
        d_eq <= d < d_eq + range      : adhere, F = -cc_adhesion * (d - d_eq)
        else                          : 0

    A hex packing at pitch ``d_eq`` is therefore force-balanced at rest, and a
    displaced cell feels a restoring force (cohesion / tissue surface tension).
    Returns an ``(M, 2)`` net force per cell.
    """

    centers = np.asarray(centers, dtype=float)
    m = len(centers)
    forces = np.zeros((m, 2))
    d_eq = cfg.cell_spacing
    cutoff = d_eq + cfg.cc_adhesion_range
    for a in range(m):
        for b in range(a + 1, m):
            delta = centers[b] - centers[a]
            dist = float(np.hypot(delta[0], delta[1]))
            if dist < 1e-9 or dist >= cutoff:
                continue
            unit = delta / dist
            if dist < d_eq:
                mag = cfg.cc_repulsion * (d_eq - dist)      # >0 -> push apart
            else:
                mag = -cfg.cc_adhesion * (dist - d_eq)      # <0 -> pull together
            forces[a] -= mag * unit
            forces[b] += mag * unit
    return forces


# =================================================================================
# Multi-cell coupling to the collagen network
# =================================================================================
def multi_cell_repulsion(network: Network, centers: np.ndarray) -> np.ndarray:
    """Soft no-penetration of collagen beads from EVERY cell disk (summed)."""

    force = np.zeros_like(network.r)
    reach = network.cfg.cell_radius + network.cfg.cell_clearance
    for center in centers:
        radial = network.r - center
        radius = np.linalg.norm(radial, axis=1)
        penetration = np.maximum(0.0, reach - radius)
        active = penetration > 0.0
        if np.any(active):
            force[active] += (
                network.cfg.repulsion_stiffness
                * penetration[active, None]
                * radial[active]
                / np.maximum(radius[active], 1e-12)[:, None]
            )
    return force


def organoid_active_forces(network: Network, centers: np.ndarray, total_force: float):
    """Every cell grips nearby fibres all-around and reels them inward.

    Reuses the G2 Gaussian contact kernel per cell (full 360 deg sector).  Only
    cells whose contact shell overlaps a fibre contribute -- i.e. the organoid
    SURFACE cells -- so the collective effect is an inward radial pull on the
    surrounding matrix.  Returns ``(nodal_force, per_cell_patches)``.
    """

    nodal = np.zeros_like(network.r)
    per_cell: list[list] = []
    for center in centers:
        patches = contact_patches(network, center, angle_deg=0.0, half_width_deg=180.0)
        per_cell.append(patches)
        if not patches:
            continue
        f, _ = forces_from_patches(network, patches, total_force)
        nodal += f
    nodal[network.fixed] = 0.0
    return nodal, per_cell


# =================================================================================
# Metrics
# =================================================================================
def radial_alignment_profile(
    network: Network, organoid_center: np.ndarray, n_shells: int = 8, shell_width: float = 25.0
):
    """Radial fibre-alignment order in distance shells from the organoid centre.

    For each collagen segment, ``order = 2*cos^2(theta) - 1`` where ``theta`` is
    the angle between the segment tangent and the local radial direction:
    ``+1`` fully radial (TACS-3-like), ``-1`` tangential, ``0`` isotropic.
    Also returns the mean bead displacement per shell.  (Lee 2017 alignment
    method; Saraswathibhatla 2025 radial index.)
    """

    center = np.asarray(organoid_center, dtype=float)
    i, j = network.edges.T
    seg = network.r[j] - network.r[i]
    tangent = seg / np.maximum(np.linalg.norm(seg, axis=1), 1e-12)[:, None]
    midpoint = 0.5 * (network.r[i] + network.r[j])
    radial = midpoint - center
    radius = np.linalg.norm(radial, axis=1)
    er = radial / np.maximum(radius, 1e-12)[:, None]
    order = 2.0 * np.square(np.sum(tangent * er, axis=1)) - 1.0

    bead_radius0 = np.linalg.norm(network.r0 - center, axis=1)
    bead_disp = np.linalg.norm(network.r - network.r0, axis=1)

    edges: list[dict] = []
    for s in range(n_shells):
        low = s * shell_width
        high = low + shell_width
        emask = (radius >= low) & (radius < high)
        bmask = (bead_radius0 >= low) & (bead_radius0 < high)
        edges.append(
            {
                "shell": (low, high),
                "n_segments": int(np.count_nonzero(emask)),
                "radial_order": float(np.mean(order[emask])) if np.any(emask) else 0.0,
                "mean_displacement": float(np.mean(bead_disp[bmask])) if np.any(bmask) else 0.0,
            }
        )
    return {
        "global_radial_order": float(np.mean(order)),
        "shells": edges,
    }


# =================================================================================
# Stage B driver: fixed contractile organoid -> radial alignment
# =================================================================================
def run_organoid_pull(cfg: OrganoidConfig = OrganoidConfig(), seed=None) -> dict:
    """Integrate a fixed (non-translating) contractile organoid pulling collagen.

    Cells stay put; each surface cell reels in nearby fibres.  The gate is whether
    the surrounding collagen's radial-alignment order RISES over the run (Hongbo's
    minimal criterion: pull -> alignment).  Cell-cell forces are computed and
    reported but do not move the fixed cells (they act in Stage D).
    """

    network, centers, gap_radius, report = make_organoid(cfg, seed=seed)
    organoid_center = np.zeros(2)
    drag = cfg.bead_drag

    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    contact_every = max(1, int(round(cfg.contact_update_interval / cfg.dt)))

    frames: list[dict] = []
    # ``active_full`` is the nodal pull at full magnitude for the current material
    # contacts; forces_from_patches is linear in the force, so each step just scales
    # it by the ramp -- no per-step per-cell loop.  Refreshed only on contact update.
    active_full, per_cell = organoid_active_forces(network, centers, cfg.total_pull_force)

    for step in range(nsteps + 1):
        time = step * cfg.dt
        ramp = min(1.0, time / cfg.force_ramp_time) if cfg.force_ramp_time else 1.0
        if step and step % contact_every == 0:
            active_full, per_cell = organoid_active_forces(
                network, centers, cfg.total_pull_force
            )
        active = active_full * ramp

        if step % every == 0:
            prof = radial_alignment_profile(network, organoid_center)
            frames.append(
                {
                    "time": time,
                    "global_radial_order": prof["global_radial_order"],
                    "shells": prof["shells"],
                    "n_gripping_cells": int(sum(1 for p in per_cell if p)),
                    "rms_bead_disp": float(
                        np.sqrt(np.mean(np.sum((network.r - network.r0) ** 2, axis=1)))
                    ),
                }
            )

        if step == nsteps:
            break

        fs, fb, fx, energy, strain, link_force = network.elastic_forces(True)
        total = fs + fb + fx + multi_cell_repulsion(network, centers) + active
        total[network.fixed] = 0.0
        network.r[~network.fixed] += cfg.dt * (total[~network.fixed] / drag)
        if not np.all(np.isfinite(network.r)):
            raise FloatingPointError("non-finite bead position; reduce dt")

    return {
        "config": asdict(cfg),
        "centers": centers,
        "gap_radius": gap_radius,
        "organoid_center": organoid_center,
        "connectivity": report,
        "n_cells": len(centers),
        "n_beads": len(network.r),
        "n_fibers": len(network.fibers),
        "n_crosslinks": len(network.crosslinks),
        "frames": frames,
        "final_positions": network.r.copy(),
        "initial_positions": network.r0.copy(),
        "edges": network.edges.copy(),
        "cell_cell_rest_force_max": float(
            np.max(np.linalg.norm(cell_cell_forces(centers, cfg), axis=1))
        ),
    }


def parameter_variant(cfg: OrganoidConfig, **changes) -> OrganoidConfig:
    """Public helper for sweeps / ablations (mirrors the G2 helper)."""

    return replace(cfg, **changes)
