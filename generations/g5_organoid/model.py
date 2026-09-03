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

try:  # Numba accelerates the per-step force balance; NumPy is the fallback.
    from numba import njit
except Exception:  # pragma: no cover
    njit = None

# Reuse the VALIDATED g4 motor-clutch law verbatim (Bell 1978; Chan & Odde 2008;
# Adebowale 2021).  G5 flattens M cells x sectors into the "sites" axis so these
# per-site step functions apply unchanged.
_REPO = Path(__file__).resolve().parents[2]   # repo root (for generations.g4_* imports)
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
from generations.g4_interactive_calibration.model import (  # noqa: E402
    bell_off_rate, patch_point)
from generations.g4_v2_multiscale.model import (  # noqa: E402
    ClutchState, _independent_step, _shared_step, shared_load_hazard,
    _counter_uniform_vector)


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
    domain_size: float = 440.0

    # --- single cell ---
    cell_radius: float = 9.0             # one tumour cell (movies: ~15-20 um dia)

    # --- organoid packing ---
    organoid_radius: float = 63.0        # target radius for cell CENTRES (core)
    cell_spacing: float = 18.0           # hex lattice pitch (= 2*cell_radius touch)
    gap: float = 1.5                     # fibre-free ring beyond the organoid at t=0
    n_corona_fibers: int = 90            # short fibres hugging the organoid (grippable collagen)
    corona_band: float = 12.0            # radial thickness of the near-field corona

    # --- simplified cell-cell interaction (piecewise-linear soft adhesive disk) ---
    # Equilibrium centre separation is ``cell_spacing``; below it disks repel,
    # above it they adhere out to ``cc_adhesion_range``, then feel nothing.
    cc_repulsion: float = 40.0           # nN/um, overlap stiffness (keeps disks apart)
    cc_adhesion: float = 6.0             # nN/um, short-range cohesion (surface tension)
    cc_adhesion_range: float = 6.0       # um, how far past contact adhesion reaches

    # --- Stage D: released rigid cells (translation under overdamped dynamics) ---
    # Each cell feels the reaction of the traction it applies to collagen (grip-and-
    # reel pulls the cell toward the ECM = outward invasion) plus cell-cell adhesion.
    # Adhesion strength sets collective (front stays cohesive) vs single-cell escape.
    cell_drag: float = 400.0             # nN*s/um, cell body drag
    max_cell_speed: float = 0.02         # um/s, biological guard on cell speed

    # --- Stage E: matrix plasticity (Bell crosslink rupture + reformation) ---
    # Loaded crosslinks rupture at Bell rate k_off0*exp(F/F_b); fibres then re-weld at
    # their CURRENT crossings.  Breaking restoring welds + re-welding the deformed
    # config makes strain irreversible (Nam 2016 plasticity index kappa; Wisdom 2018).
    # This is crosslink-topology plasticity, distinct from SLS (still deferred) and
    # from Bell on the CLUTCHES.  Check crosslink-force magnitudes before trusting it.
    plasticity: bool = False
    xl_k_off0: float = 2e-4              # s^-1, zero-force rupture rate (stress-selective:
                                         #   negligible at F=0, accelerates as exp(F/F_b))
    xl_F_b: float = 3.0                  # nN, Bell force scale for crosslink rupture
    xl_reform_prob: float = 0.5          # prob a fresh fibre crossing re-welds per update
    plasticity_interval: float = 5.0     # s between rupture/reform updates
    unload_time: float = 0.0             # s at end with pull ramped to 0 (for kappa)

    # --- cell-fibre MOLECULAR CLUTCH with slippage (opt-in; replaces constant pull) ---
    # A stochastic motor-clutch at every cell grip site: a motor pulls actin at v0 with
    # force-velocity v=v0(1-F/F_stall) (Chan & Odde 2008 Science), clutches load
    # (F=k_c.ext), unbind at Bell rate k_off0.exp(F/F_b) (Bell 1978), rebind at k_on;
    # the cluster load-and-fails and traction resets (slippage).  Distinct from the
    # crosslinker plasticity of Stage E.  Values inherited from G4/Adebowale 2021 Nat
    # Mater SI Table 4 (docs/g4_parameter_provenance.md) -- do NOT invent.
    clutch_dynamics: bool = False        # off -> the validated constant-pull path
    clutch_mode: str = "independent"     # "independent" | "shared" (g4_v2 load sharing)
    n_contact_sectors: int = 12          # grip sites per cell (angular sectors)
    n_clutches_per_site: int = 12        # effective clutch bundle per site (Adebowale 2021)
    clutch_stiffness: float = 2.0        # k_c, nN/um (G2 V3 / g4 provenance)
    clutch_on_rate: float = 0.055        # k_on, 1/s (Adebowale 2021 SI Table 4)
    clutch_off_rate0: float = 0.018      # k_off0, 1/s zero-force (Adebowale 2021)
    bell_force: float = 1.5              # F_b, nN Bell scale (Bell 1978; Adebowale 2021)
    unloaded_actin_speed: float = 0.025  # v0, um/s (~24 nm/s, Adebowale 2021)
    motor_stall_per_site: float = 8.0    # F_stall, nN/site (Chan & Odde 2008 scale)
    clutch_counter_seed: int = 3042      # deterministic counter-RNG seed (g4_v2)

    # --- collagen network (scaled up from the 99-fibre / 180 um G2 baseline) ---
    # Softened, lightly crosslinked collagen so the pull can visibly reorganise the
    # near-field into a radial aster (same prof-requested softening as G3/G4:
    # 3 MPa modulus, 10 nN/um links).  Not the stiff 32 MPa G2 default.
    collagen_modulus_mpa: float = 3.0
    crosslink_stiffness: float = 10.0
    # --- Stage C: nonlinear fibre strain-stiffening (off = Stage B linear) ---
    # Tensile stiffness multiplies by exp(strain / stiffen_strain_ref), capped, so
    # taut (radial) fibres stiffen and transmit the pull further -> sharper, longer-
    # range aster (Steinwachs 2016; Mark 2020 eLife 51912; Shenoy 2014). Compression
    # stays soft (microbuckling), unchanged.  The on/off flag is the Stage-C ablation.
    strain_stiffening: bool = False
    stiffen_strain_ref: float = 0.06     # e-fold tensile stiffening strain
    stiffen_cap: float = 12.0            # max stiffening multiple (numerical guard)
    n_fibers: int = 420
    crosslink_fraction: float = 0.85     # keep most auto-generated links (connectivity)
    boundary_width: float = 6.0          # anchoring band; wider -> better boundary percolation
    contact_width: float = 8.0           # grip shell so surface cells reach the corona
    total_pull_force: float = 12.0       # nN, per cell (tens-of-nN traction scale)

    # --- integration / output ---
    dt: float = 0.05
    duration: float = 120.0
    sample_interval: float = 2.0
    contact_update_interval: float = 4.0  # cells are fixed -> contacts drift slowly
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
# Network generation excluding the UNION of cell disks (collagen reaches every
# perimeter cell, not just the outermost ring of a circular gap)
# =================================================================================
def _inside_any_cell(points: np.ndarray, centers: np.ndarray, clearance: float) -> bool:
    """True if any point lies within (cell_radius + clearance) of any cell centre."""
    # points (P,2), centers (M,2) -> min distance to a centre per point
    diff = points[:, None, :] - centers[None, :, :]
    dmin = np.sqrt(np.min(np.sum(diff * diff, axis=2), axis=1))
    return bool(np.any(dmin < clearance))


def _organoid_fiber(cfg: OrganoidConfig, rng, centers: np.ndarray, exclude: float,
                    organoid_outer: float, *, boundary_seeded: bool):
    """One curved fibre that avoids every cell disk (union exclusion)."""

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
    if _inside_any_cell(points, centers, exclude):
        return None
    return points


def _corona_fiber(cfg: OrganoidConfig, rng, centers: np.ndarray, exclude: float,
                  organoid_outer: float):
    """A grippable near-field fibre: inner end at the organoid, extending OUTWARD.

    Fix (B): instead of short tangential stubs that float, seed a longer fibre whose
    inner end presses against the perimeter (so surface cells grip it) and which
    extends outward with a radial-but-noisy heading (like the radial TACS-3 tracts
    seen around invasive tumours).  Being longer and radial-ish, it crosses many bulk
    fibres -> it crosslinks into the boundary-connected network instead of floating,
    and it carries the pull farther out.  The wide angular noise keeps t=0 only mildly
    ordered so reorganisation is still visible.
    """

    phi = rng.uniform(0.0, 2.0 * math.pi)
    inner_r = organoid_outer + rng.uniform(0.5, cfg.gap + 3.0)
    radial = np.array([math.cos(phi), math.sin(phi)])
    inner = inner_r * radial
    length = rng.uniform(0.6 * cfg.max_fiber_length, cfg.max_fiber_length)
    # heading: outward (+radial) with wide noise (~+/-80 deg) -> bridges to the bulk
    # (connectivity) yet starts near-isotropic so reorganisation is still visible.
    heading = phi + rng.uniform(-1.4, 1.4)
    direction = np.array([math.cos(heading), math.sin(heading)])
    start = inner
    end = inner + length * direction
    half = cfg.domain_size / 2.0
    if np.max(np.abs(end)) > half - 0.3 or np.max(np.abs(start)) > half - 0.3:
        return None
    points = _curved_segment(start, end, cfg, rng)
    if _inside_any_cell(points, centers, exclude):
        return None
    return points


def _assemble_organoid_spec(cfg: OrganoidConfig, rng, centers: np.ndarray,
                            organoid_outer: float, seed_used: int) -> NetworkSpec:
    half = cfg.domain_size / 2.0
    exclude = cfg.cell_radius + cfg.cell_clearance + cfg.gap
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

    # 1) near-field corona hugging the organoid perimeter (grippable collagen)
    corona = 0
    attempts = 0
    while corona < cfg.n_corona_fibers and attempts < 400_000:
        attempts += 1
        pts = _corona_fiber(cfg, rng, centers, exclude, organoid_outer)
        if pts is not None:
            append(pts)
            corona += 1
    contact_fibers = list(range(len(fibers)))  # corona fibres are the grippable set

    # 2) bulk network (boundary-seeded first for percolation, then interior)
    target_total = cfg.n_fibers
    target_boundary = corona + max(24, int(round(0.45 * (target_total - corona))))
    attempts = 0
    while len(fibers) < target_total and attempts < 600_000:
        attempts += 1
        pts = _organoid_fiber(
            cfg, rng, centers, exclude, organoid_outer,
            boundary_seeded=len(fibers) < target_boundary,
        )
        if pts is not None:
            append(pts)
    if len(fibers) != target_total:
        raise RuntimeError("could not construct the requested organoid network")
    return NetworkSpec(np.asarray(positions), fibers, np.asarray(fixed), contact_fibers, seed_used)


def build_crosslinks_grid(network: Network) -> list:
    """O(E) spatial-grid crosslinker: intersections of edges on DIFFERENT fibres.

    Replaces the G2 O(F^2) fibre-pair double loop.  Edges (bead-to-bead segments,
    ~bead_spacing long) are binned by midpoint; only edges in the same / adjacent
    bins are intersection-tested.  Produces the same permanent hinged links joining
    coincident material points on different fibres, deduplicated by proximity.
    """

    from common.model import Crosslink  # local import to avoid cycle at top

    r = network.r
    edges = network.edges
    efib = network.edge_fiber
    p = r[edges[:, 0]]
    d = r[edges[:, 1]] - p
    mid = p + 0.5 * d
    bin_size = max(2.0, 2.0 * float(network.cfg.bead_spacing))
    keys = np.floor(mid / bin_size).astype(np.int64)

    grid: dict = {}
    for e in range(len(edges)):
        grid.setdefault((keys[e, 0], keys[e, 1]), []).append(e)

    stiffness = network.cfg.crosslink_stiffness
    links: list = []
    seen: set = set()
    neigh = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]
    for (bx, by), es in grid.items():
        cand: list = []
        for dx, dy in neigh:
            cand.extend(grid.get((bx + dx, by + dy), ()))
        for e in es:
            pe = p[e]; de = d[e]; fe = efib[e]
            for f in cand:
                if f <= e or efib[f] == fe:
                    continue
                pf = p[f]; df = d[f]
                den = de[0] * df[1] - de[1] * df[0]
                if abs(den) < 1e-10:
                    continue
                qp = pf - pe
                ta = (qp[0] * df[1] - qp[1] * df[0]) / den
                tb = (qp[0] * de[1] - qp[1] * de[0]) / den
                if not (1e-5 < ta < 1.0 - 1e-5 and 1e-5 < tb < 1.0 - 1e-5):
                    continue
                point = pe + ta * de
                fa, fb = (fe, efib[f]) if fe < efib[f] else (efib[f], fe)
                key = (int(fa), int(fb), int(round(point[0] / 0.6)), int(round(point[1] / 0.6)))
                if key in seen:
                    continue
                seen.add(key)
                links.append(Crosslink(int(e), float(ta), int(f), float(tb), np.zeros(2), stiffness))
    return links


def _link_key(network: Network, link) -> tuple:
    """Coarse identity of a crosslink by fibre pair + rounded material point.

    Position-based (uses current geometry) -- used for reform de-duplication, where
    two welds at the same current crossing should count as one.
    """
    fa = int(network.edge_fiber[link.edge_a])
    fb = int(network.edge_fiber[link.edge_b])
    pt = network.material_point(link.edge_a, link.alpha_a)
    lo, hi = (fa, fb) if fa < fb else (fb, fa)
    return (lo, hi, int(round(pt[0] / 0.6)), int(round(pt[1] / 0.6)))


def _link_id(link) -> tuple:
    """Deformation-INVARIANT identity: which edges + where along them (material points).

    Unlike ``_link_key`` this does not move as the network deforms, so it correctly
    tracks whether a specific weld survives -- used for the topological plasticity
    index kappa_topo.
    """
    return (int(link.edge_a), int(link.edge_b),
            round(float(link.alpha_a), 2), round(float(link.alpha_b), 2))


def update_plasticity(network: Network, cfg: OrganoidConfig, dt_eff: float, rng) -> tuple:
    """Bell rupture of loaded crosslinks + re-weld at current fibre crossings.

    Ruptured welds no longer pull the fibres back; new welds (rest = current, i.e.
    zero strain at the deformed crossing) pin the deformed configuration -> the strain
    becomes irreversible (Nam 2016; Wisdom 2018).  Returns ``(ruptured, formed)``.
    """

    ruptured = 0
    if len(network.crosslinks):
        *_, link_force = network.elastic_forces(True)          # |force| per link, nN
        rate = cfg.xl_k_off0 * np.exp(np.minimum(link_force / cfg.xl_F_b, 50.0))
        prob = 1.0 - np.exp(-rate * dt_eff)
        keep = rng.random(len(network.crosslinks)) >= prob
        ruptured = int(np.count_nonzero(~keep))
        if ruptured:
            network.crosslinks = [x for x, k in zip(network.crosslinks, keep) if k]

    # re-weld: current geometric crossings not already linked
    existing = {_link_key(network, x) for x in network.crosslinks}
    formed = 0
    for cand in build_crosslinks_grid(network):
        if _link_key(network, cand) in existing:
            continue
        if rng.random() < cfg.xl_reform_prob:
            network.crosslinks.append(cand)                    # rest_vector = 0 at the
            existing.add(_link_key(network, cand))             # current (deformed) crossing
            formed += 1

    network.refresh_crosslink_arrays()
    return ruptured, formed


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
    best = (-1, -1.0, None)  # (contact_connected, connected_fraction, network)
    for attempt in range(cfg.generation_attempts):
        seed_used = base + 7919 * attempt
        spec = _assemble_organoid_spec(
            cfg, np.random.default_rng(seed_used), centers, organoid_outer, seed_used
        )
        network = Network(spec, cfg, crosslinks=[])  # skip G2 O(F^2) auto-build
        network.crosslinks = build_crosslinks_grid(network)  # O(E) grid crosslinker
        network.refresh_crosslink_arrays()
        if cfg.crosslink_fraction < 1.0 and network.crosslinks:
            rng = np.random.default_rng(seed_used + 101)
            keep = rng.random(len(network.crosslinks)) < cfg.crosslink_fraction
            network.crosslinks = [x for x, k in zip(network.crosslinks, keep) if k]
            network.refresh_crosslink_arrays()
        report = connectivity_report(network)
        score = float(report["connected_fraction"])
        # Fix (A): require the grippable near-field (contact) fibres to be
        # boundary-connected -- the G2 gate that was dropped in the scale-up.  Rank
        # candidates by (contact-connected, then fraction) so the accepted network
        # never has floating grippable fibres.
        cc = 1 if report["contact_fibers_connected"] else 0
        if (cc, score) > (best[0], best[1]):
            best = (cc, score, network)
        if report["contact_fibers_connected"] and score >= cfg.required_connected_fraction:
            return network, centers, gap_radius, report
    if best[2] is not None:
        return best[2], centers, gap_radius, connectivity_report(best[2])
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
    d_eq = cfg.cell_spacing
    cutoff = d_eq + cfg.cc_adhesion_range
    d = centers[None, :, :] - centers[:, None, :]          # d[a,b] = c_b - c_a
    dist = np.linalg.norm(d, axis=2)
    safe = np.maximum(dist, 1e-12)
    unit = d / safe[:, :, None]
    mag = np.where(dist < d_eq, cfg.cc_repulsion * (d_eq - dist),
                   np.where(dist < cutoff, -cfg.cc_adhesion * (dist - d_eq), 0.0))
    mag[dist < 1e-9] = 0.0                                  # self / coincident
    return -np.sum(mag[:, :, None] * unit, axis=1)         # F_a = -sum_b mag*unit


# =================================================================================
# Numba per-step integrator (generalised to M cell centres)
# =================================================================================
if njit is not None:

    @njit(cache=True)
    def _advance_organoid_numba(
        r, r0, fixed, edges, l0, k_tension, k_compression,
        triplets, curvature0, bend_coefficient,
        link_edge_a, link_edge_b, link_alpha_a, link_alpha_b,
        link_rest, link_stiffness, active, centers,
        cell_radius, cell_clearance, repulsion_stiffness,
        bead_drag, dt, stiffen_on, stiffen_ref, stiffen_cap, velocity, total,
    ):
        """Allocation-free G2 force balance with repulsion summed over M disks.

        Same bend / crosslink law as G2 ``elastic_forces`` (lifted from g4_v2
        ``_advance_numba``).  Generalisations: repulsion sums over every cell centre,
        and (Stage C) tensile stiffness stiffens as exp(strain/ref) when stiffen_on.
        """

        total[:, :] = 0.0
        for e in range(edges.shape[0]):
            i = edges[e, 0]
            j = edges[e, 1]
            dx = r[j, 0] - r[i, 0]
            dy = r[j, 1] - r[i, 1]
            length = math.sqrt(dx * dx + dy * dy)
            if length < 1e-12:
                length = 1e-12
            extension = length - l0[e]
            if extension >= 0.0:
                stiffness = k_tension[e]
                if stiffen_on == 1:
                    mult = math.exp((extension / l0[e]) / stiffen_ref)
                    if mult > stiffen_cap:
                        mult = stiffen_cap
                    stiffness = stiffness * mult
            else:
                stiffness = k_compression[e]
            scale = stiffness * extension / length
            fx = scale * dx
            fy = scale * dy
            total[i, 0] += fx
            total[i, 1] += fy
            total[j, 0] -= fx
            total[j, 1] -= fy

        for t in range(triplets.shape[0]):
            a = triplets[t, 0]
            b = triplets[t, 1]
            c = triplets[t, 2]
            qx = r[a, 0] - 2.0 * r[b, 0] + r[c, 0] - curvature0[t, 0]
            qy = r[a, 1] - 2.0 * r[b, 1] + r[c, 1] - curvature0[t, 1]
            fx = bend_coefficient * qx
            fy = bend_coefficient * qy
            total[a, 0] -= fx
            total[a, 1] -= fy
            total[b, 0] += 2.0 * fx
            total[b, 1] += 2.0 * fy
            total[c, 0] -= fx
            total[c, 1] -= fy

        for x in range(link_edge_a.shape[0]):
            ea = link_edge_a[x]
            eb = link_edge_b[x]
            ia = edges[ea, 0]
            ja = edges[ea, 1]
            ib = edges[eb, 0]
            jb = edges[eb, 1]
            aa = link_alpha_a[x]
            ab = link_alpha_b[x]
            pax = (1.0 - aa) * r[ia, 0] + aa * r[ja, 0]
            pay = (1.0 - aa) * r[ia, 1] + aa * r[ja, 1]
            pbx = (1.0 - ab) * r[ib, 0] + ab * r[jb, 0]
            pby = (1.0 - ab) * r[ib, 1] + ab * r[jb, 1]
            fx = link_stiffness[x] * (pbx - pax - link_rest[x, 0])
            fy = link_stiffness[x] * (pby - pay - link_rest[x, 1])
            total[ia, 0] += (1.0 - aa) * fx
            total[ia, 1] += (1.0 - aa) * fy
            total[ja, 0] += aa * fx
            total[ja, 1] += aa * fy
            total[ib, 0] -= (1.0 - ab) * fx
            total[ib, 1] -= (1.0 - ab) * fy
            total[jb, 0] -= ab * fx
            total[jb, 1] -= ab * fy

        min_radius = cell_radius + cell_clearance
        for i in range(r.shape[0]):
            for c in range(centers.shape[0]):
                dx = r[i, 0] - centers[c, 0]
                dy = r[i, 1] - centers[c, 1]
                radius = math.sqrt(dx * dx + dy * dy)
                if radius < min_radius:
                    safe = radius if radius > 1e-12 else 1e-12
                    force = repulsion_stiffness * (min_radius - radius) / safe
                    total[i, 0] += force * dx
                    total[i, 1] += force * dy
            total[i, 0] += active[i, 0]
            total[i, 1] += active[i, 1]
            if fixed[i]:
                velocity[i, 0] = 0.0
                velocity[i, 1] = 0.0
                r[i, 0] = r0[i, 0]
                r[i, 1] = r0[i, 1]
            else:
                velocity[i, 0] = total[i, 0] / bead_drag
                velocity[i, 1] = total[i, 1] / bead_drag
                r[i, 0] += dt * velocity[i, 0]
                r[i, 1] += dt * velocity[i, 1]


class OrganoidStepper:
    """Preallocated overdamped integrator (Numba if available, else NumPy)."""

    def __init__(self, network: Network, centers: np.ndarray):
        self.n = network
        self.centers = np.ascontiguousarray(centers, dtype=float)
        self.velocity = np.zeros_like(network.r)
        self.total = np.zeros_like(network.r)

    def step(self, active: np.ndarray, dt: float) -> None:
        n = self.n
        stiffen_on = 1 if getattr(n.cfg, "strain_stiffening", False) else 0
        stiffen_ref = float(getattr(n.cfg, "stiffen_strain_ref", 0.06))
        stiffen_cap = float(getattr(n.cfg, "stiffen_cap", 12.0))
        if njit is not None:
            _advance_organoid_numba(
                n.r, n.r0, n.fixed, n.edges, n.l0, n.k_tension, n.k_compression,
                n.triplets, n.curvature0, n.bend_coefficient,
                n.link_edge_a, n.link_edge_b, n.link_alpha_a, n.link_alpha_b,
                n.link_rest, n.link_stiffness, active, self.centers,
                n.cfg.cell_radius, n.cfg.cell_clearance, n.cfg.repulsion_stiffness,
                n.cfg.bead_drag, dt, stiffen_on, stiffen_ref, stiffen_cap,
                self.velocity, self.total,
            )
        else:  # pragma: no cover - exercised only without Numba
            _, fb, fx, _, _, _ = n.elastic_forces(True)
            fs = _stiffened_stretch(n, stiffen_on, stiffen_ref, stiffen_cap)
            total = fs + fb + fx + multi_cell_repulsion(n, self.centers) + active
            total[n.fixed] = 0.0
            n.r[~n.fixed] += dt * (total[~n.fixed] / n.cfg.bead_drag)
        if not np.all(np.isfinite(n.r)):
            raise FloatingPointError("non-finite bead position; reduce dt")


def _stiffened_stretch(n: Network, stiffen_on: int, ref: float, cap: float) -> np.ndarray:
    """Vectorised stretch force with optional exp strain-stiffening (NumPy fallback)."""

    i, j = n.edges.T
    d = n.r[j] - n.r[i]
    length = np.linalg.norm(d, axis=1)
    extension = length - n.l0
    taut = extension >= 0.0
    stiffness = np.where(taut, n.k_tension, n.k_compression)
    if stiffen_on == 1:
        mult = np.minimum(cap, np.exp(np.where(taut, extension / n.l0, 0.0) / ref))
        stiffness = np.where(taut, n.k_tension * mult, n.k_compression)
    pair = (stiffness * extension / np.maximum(length, 1e-12))[:, None] * d
    fs = np.zeros_like(n.r)
    np.add.at(fs, i, pair)
    np.add.at(fs, j, -pair)
    return fs


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


def cell_candidate_fibers(network: Network, centers: np.ndarray, reach: float) -> list:
    """Fibres with any bead within ``reach`` of each cell centre (computed once).

    Cells are fixed in Stage B, so this list is stable and lets the contact search
    skip the O(M*F) sweep over every fibre -- only these few near fibres can be
    gripped.
    """

    fiber_of_bead = np.empty(len(network.r), dtype=np.int64)
    for fid, ids in enumerate(network.fibers):
        fiber_of_bead[ids] = fid
    out: list = []
    for center in centers:
        dist = np.linalg.norm(network.r - center, axis=1)
        out.append(np.unique(fiber_of_bead[dist < reach]))
    return out


def _cell_patches(network: Network, center: np.ndarray, fiber_ids) -> list:
    """G2 hard-contact + Gaussian weighting, restricted to candidate fibres."""

    from common.model import ContactPatch  # local import

    cfg = network.cfg
    candidates: list = []
    for fid in fiber_ids:
        edge_ids = np.flatnonzero(network.edge_fiber == fid)
        if not len(edge_ids):
            continue
        pairs = network.edges[edge_ids]
        a = network.r[pairs[:, 0]]
        d = network.r[pairs[:, 1]] - a
        alpha = np.clip(
            np.sum((center - a) * d, axis=1) / np.maximum(np.sum(d * d, axis=1), 1e-12),
            0.0, 1.0,
        )
        point = a + alpha[:, None] * d
        radial = point - center
        radius = np.linalg.norm(radial, axis=1)
        surface = radius - cfg.cell_radius
        eligible = (surface >= 0.0) & (surface <= cfg.contact_width)
        if np.any(eligible):
            local = int(np.flatnonzero(eligible)[np.argmin(surface[eligible])])
            candidates.append(
                ContactPatch(
                    int(fid), int(edge_ids[local]), float(alpha[local]),
                    point[local].copy(), float(surface[local]), 0.0,
                    -radial[local] / max(float(radius[local]), 1e-12),
                )
            )
    if not candidates:
        return []
    raw = np.exp(-np.square([x.surface_distance for x in candidates]) / cfg.gaussian_sigma**2)
    raw /= np.sum(raw)
    for patch, weight in zip(candidates, raw):
        patch.weight = float(weight)
    return candidates


def organoid_active_forces(network: Network, centers: np.ndarray, total_force: float,
                           candidates: list | None = None):
    """Every cell grips nearby fibres all-around and reels them inward.

    Reuses the G2 Gaussian contact kernel per cell.  Only cells whose contact shell
    overlaps a fibre contribute -- the organoid SURFACE cells -- so the collective
    effect is an inward radial pull on the surrounding matrix.  ``candidates`` is a
    per-cell list of near-fibre ids (see :func:`cell_candidate_fibers`); when given,
    the contact search skips the full-fibre sweep.  Returns ``(nodal, per_cell)``.
    """

    nodal = np.zeros_like(network.r)
    per_cell: list[list] = []
    for c, center in enumerate(centers):
        if candidates is not None:
            patches = _cell_patches(network, center, candidates[c])
        else:
            patches = contact_patches(network, center, angle_deg=0.0, half_width_deg=180.0)
        per_cell.append(patches)
        if patches:
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
# Cell-fibre molecular clutch with slippage (opt-in; Bell 1978, Chan & Odde 2008)
# =================================================================================
def _clutch_counter_uniforms(cfg: OrganoidConfig, step: int, channel: int, n_sites: int):
    """Deterministic counter-addressed uniforms over ALL sites (g4_v2 hash)."""
    v = _counter_uniform_vector(cfg.clutch_counter_seed, step, channel,
                                n_sites * cfg.n_clutches_per_site)
    return v.reshape(n_sites, cfg.n_clutches_per_site)


def organoid_clutch_patches(network: Network, centers: np.ndarray, cfg: OrganoidConfig,
                            candidates: list):
    """One grip site per (cell, angular sector): the closest gripped fibre in that sector.

    Generalises g4's per-cell ``n_contact_sectors`` sites to M cells.  Returns a flat
    list of length ``M * n_contact_sectors`` (a ``ContactPatch`` or ``None`` per site,
    ``None`` = that sector has no fibre in reach) plus a matching ``(S, 2)`` array of the
    owning cell centre for each site.
    """
    S = cfg.n_contact_sectors
    patches: list = []
    site_centers: list = []
    for c, center in enumerate(centers):
        per_sector = [None] * S
        best_surf = [np.inf] * S
        for p in _cell_patches(network, center, candidates[c]):
            ang = math.atan2(p.point[1] - center[1], p.point[0] - center[0]) % (2 * math.pi)
            sec = int(ang / (2 * math.pi) * S) % S
            if p.surface_distance < best_surf[sec]:
                best_surf[sec] = p.surface_distance
                per_sector[sec] = p
        patches.extend(per_sector)
        site_centers.extend([center] * S)
    return patches, np.asarray(site_centers, dtype=float)


def _project_site_forces(network: Network, patches: list, site_force: np.ndarray) -> np.ndarray:
    """Project each site's clutch traction onto its fibre beads (first-moment preserving)."""
    nodal = np.zeros_like(network.r)
    for s, patch in enumerate(patches):
        if patch is None or site_force[s] <= 0.0:
            continue
        f = site_force[s] * patch.normal_in            # inward, toward the owning cell
        i, j = network.edges[patch.edge]
        nodal[i] += (1.0 - patch.alpha) * f
        nodal[j] += patch.alpha * f
    nodal[network.fixed] = 0.0
    return nodal


def _clutch_substrate_speeds(network: Network, velocity: np.ndarray, patches: list,
                             site_centers: np.ndarray) -> np.ndarray:
    """Inward speed of each gripped fibre material point toward its owning cell.

    Feeds back into clutch loading: a fast-yielding (soft) fibre loads the clutches
    less (Chan & Odde 2008 motor-clutch).  Empty sites report 0.
    """
    out = np.zeros(len(patches))
    for s, patch in enumerate(patches):
        if patch is None:
            continue
        point = network.material_point(patch.edge, patch.alpha)
        inward = site_centers[s] - point
        inward /= max(float(np.linalg.norm(inward)), 1e-12)
        i, j = network.edges[patch.edge]
        mv = (1.0 - patch.alpha) * velocity[i] + patch.alpha * velocity[j]
        out[s] = float(mv @ inward)
    return out


def _clutch_step(cfg: OrganoidConfig, state: ClutchState, substrate: np.ndarray, step: int):
    """Advance every site's clutch bundle one step (reuses the g4_v2 law verbatim)."""
    S = state.bound.shape[0]
    u_on = _clutch_counter_uniforms(cfg, step, 0, S)
    u_off = _clutch_counter_uniforms(cfg, step, 1, S)
    if cfg.clutch_mode == "shared":
        return _shared_step(cfg, state, substrate, u_on, u_off)
    return _independent_step(cfg, state, substrate, u_on, u_off)


def _new_clutch_state(n_sites: int, cfg: OrganoidConfig) -> ClutchState:
    return ClutchState(
        bound=np.zeros((n_sites, cfg.n_clutches_per_site), dtype=bool),
        extension=np.zeros((n_sites, cfg.n_clutches_per_site)),
        site_extension=np.zeros(n_sites),
    )


def _run_pull_with_clutch(cfg: OrganoidConfig, seed=None, snapshots: bool = False) -> dict:
    """Stage B with a stochastic motor-clutch grip (load -> Bell slip -> reset) at every
    cell sector, instead of a constant pull.  Traction on the collagen EMERGES from the
    clutches, so it is non-constant and shows load-and-fail slip events."""

    network, centers, gap_radius, report = make_organoid(cfg, seed=seed)
    organoid_center = np.zeros(2)
    stepper = OrganoidStepper(network, centers)
    reach = cfg.cell_radius + cfg.contact_width + 2.0
    candidates = cell_candidate_fibers(network, centers, reach)
    patches, site_centers = organoid_clutch_patches(network, centers, cfg, candidates)
    S = len(patches)
    state = _new_clutch_state(S, cfg)
    substrate = np.zeros(S)
    site_force = np.zeros(S)

    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    frames: list[dict] = []
    snaps: list[np.ndarray] = []
    traction_series: list[float] = []

    for step in range(nsteps + 1):
        time = step * cfg.dt
        _, site_force, breaks, binds, site_fail = _clutch_step(cfg, state, substrate, step)
        active = _project_site_forces(network, patches, site_force)
        traction_series.append(float(site_force.sum()))

        if step % every == 0:
            prof = radial_alignment_profile(network, organoid_center)
            gripping = np.asarray([p is not None for p in patches]).reshape(len(centers), -1).any(axis=1)
            frames.append({
                "time": time,
                "global_radial_order": prof["global_radial_order"],
                "shells": prof["shells"],
                "n_gripping_cells": int(gripping.sum()),
                "bound_fraction": float(state.bound.mean()),
                "cumulative_slips": int(state.cumulative_slips),
                "mean_site_force": float(site_force[site_force > 0].mean()) if np.any(site_force > 0) else 0.0,
                "total_traction": float(site_force.sum()),
                "rms_bead_disp": float(np.sqrt(np.mean(np.sum((network.r - network.r0) ** 2, axis=1)))),
            })
            if snapshots:
                snaps.append(network.r.copy())

        if step == nsteps:
            break
        stepper.step(active, cfg.dt)
        substrate = _clutch_substrate_speeds(network, stepper.velocity, patches, site_centers)

    tr = np.asarray(traction_series)
    return {
        "config": asdict(cfg), "centers": centers, "gap_radius": gap_radius,
        "organoid_center": organoid_center, "connectivity": report,
        "n_cells": len(centers), "n_beads": len(network.r), "n_fibers": len(network.fibers),
        "n_crosslinks": len(network.crosslinks), "n_clutch_sites": S,
        "clutch_mode": cfg.clutch_mode,
        "cumulative_slips": int(state.cumulative_slips),
        "cumulative_site_failures": int(state.cumulative_site_failures),
        "traction_cv": float(tr.std() / tr.mean()) if tr.mean() > 0 else 0.0,  # non-constant?
        "frames": frames, "snapshots": np.asarray(snaps) if snapshots else None,
        "final_positions": network.r.copy(), "initial_positions": network.r0.copy(),
        "edges": network.edges.copy(),
        "crosslinks": [(int(x.edge_a), float(x.alpha_a), int(x.edge_b), float(x.alpha_b))
                       for x in network.crosslinks],
    }


# =================================================================================
# Stage B driver: fixed contractile organoid -> radial alignment
# =================================================================================
def run_organoid_pull(cfg: OrganoidConfig = OrganoidConfig(), seed=None,
                      snapshots: bool = False) -> dict:
    """Integrate a fixed (non-translating) contractile organoid pulling collagen.

    Cells stay put; each surface cell reels in nearby fibres.  The gate is whether
    the surrounding collagen's radial-alignment order RISES over the run (Hongbo's
    minimal criterion: pull -> alignment).  Cell-cell forces are computed and
    reported but do not move the fixed cells (they act in Stage D).  With
    ``snapshots=True`` the full bead positions at each sampled frame are stored under
    ``snapshots`` for animation.  With ``cfg.clutch_dynamics`` the constant pull is
    replaced by a stochastic motor-clutch grip (see :func:`_run_pull_with_clutch`).
    """
    if cfg.clutch_dynamics:
        return _run_pull_with_clutch(cfg, seed=seed, snapshots=snapshots)

    network, centers, gap_radius, report = make_organoid(cfg, seed=seed)
    organoid_center = np.zeros(2)
    stepper = OrganoidStepper(network, centers)
    # Fixed cells -> the set of grippable near fibres is stable; cache it once so the
    # contact search skips the O(M*F) full-fibre sweep.
    reach = cfg.cell_radius + cfg.contact_width + 2.0
    candidates = cell_candidate_fibers(network, centers, reach)

    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    contact_every = max(1, int(round(cfg.contact_update_interval / cfg.dt)))
    plastic_every = max(1, int(round(cfg.plasticity_interval / cfg.dt)))
    rng = np.random.default_rng((cfg.seed if seed is None else int(seed)) + 555)
    load_end = cfg.duration - cfg.unload_time
    initial_link_ids = {_link_id(x) for x in network.crosslinks}

    frames: list[dict] = []
    snaps: list[np.ndarray] = []
    plastic_events = {"ruptured": 0, "formed": 0}
    # ``active_full`` is the nodal pull at full magnitude for the current material
    # contacts; forces_from_patches is linear in the force, so each step just scales
    # it by the ramp -- no per-step per-cell loop.  Refreshed only on contact update.
    active_full, per_cell = organoid_active_forces(
        network, centers, cfg.total_pull_force, candidates=candidates
    )

    for step in range(nsteps + 1):
        time = step * cfg.dt
        # load, hold, then (if unload_time>0) drop the pull to zero FAST and hold at
        # zero for the rest of the window so the network can relax -> plasticity index
        # kappa = residual / peak strain (residual is measured after the relaxation).
        if cfg.force_ramp_time and time < cfg.force_ramp_time:
            ramp = time / cfg.force_ramp_time
        elif time <= load_end:
            ramp = 1.0
        else:
            ramp = max(0.0, 1.0 - (time - load_end) / max(cfg.force_ramp_time, cfg.dt))
        if step and step % contact_every == 0:
            active_full, per_cell = organoid_active_forces(
                network, centers, cfg.total_pull_force, candidates=candidates
            )
        active = active_full * ramp

        if cfg.plasticity and step and step % plastic_every == 0:
            rup, form = update_plasticity(network, cfg, cfg.plasticity_interval, rng)
            plastic_events["ruptured"] += rup
            plastic_events["formed"] += form

        if step % every == 0:
            prof = radial_alignment_profile(network, organoid_center)
            frames.append(
                {
                    "time": time,
                    "ramp": ramp,
                    "global_radial_order": prof["global_radial_order"],
                    "shells": prof["shells"],
                    "n_gripping_cells": int(sum(1 for p in per_cell if p)),
                    "n_crosslinks": len(network.crosslinks),
                    "rms_bead_disp": float(
                        np.sqrt(np.mean(np.sum((network.r - network.r0) ** 2, axis=1)))
                    ),
                }
            )
            if snapshots:
                snaps.append(network.r.copy())

        if step == nsteps:
            break

        stepper.step(active, cfg.dt)

    # --- plasticity indices ---------------------------------------------------
    # kappa (rms recovery) is CONFOUNDED by slow large-scale elastic relaxation, so
    # two relaxation-independent measures are reported instead:
    #   kappa_topo  = fraction of the ORIGINAL crosslinks permanently replaced
    #                 (elastic = 0; plastic > 0) -- a pure topology-change signal.
    #   kappa_order = residual of the radial STRUCTURE after unload,
    #                 (S_final - S_0)/(S_peak - S_0): does the radial pattern persist?
    rms = [fr["rms_bead_disp"] for fr in frames]
    peak = max(rms) if rms else 0.0
    kappa = float(rms[-1] / peak) if (cfg.unload_time > 0 and peak > 0) else None
    final_ids = {_link_id(x) for x in network.crosslinks}
    kappa_topo = float(1.0 - len(initial_link_ids & final_ids) / max(len(initial_link_ids), 1))
    orders = [fr["global_radial_order"] for fr in frames]
    o0, opk, of = orders[0], max(orders, key=lambda v: abs(v - orders[0])), orders[-1]
    kappa_order = (float((of - o0) / (opk - o0))
                   if (cfg.unload_time > 0 and abs(opk - o0) > 1e-6) else None)

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
        "plastic_events": plastic_events,
        "kappa": kappa,
        "kappa_topo": kappa_topo,
        "kappa_order": kappa_order,
        "frames": frames,
        "snapshots": np.asarray(snaps) if snapshots else None,
        "final_positions": network.r.copy(),
        "initial_positions": network.r0.copy(),
        "edges": network.edges.copy(),
        "crosslinks": [(int(x.edge_a), float(x.alpha_a), int(x.edge_b), float(x.alpha_b))
                       for x in network.crosslinks],
        "cell_cell_rest_force_max": float(
            np.max(np.linalg.norm(cell_cell_forces(centers, cfg), axis=1))
        ),
    }


# =================================================================================
# Stage D driver: released cells (grip-reel reaction -> invasion; adhesion sets mode)
# =================================================================================
def _run_invasion_with_clutch(cfg: OrganoidConfig, seed=None, snapshots: bool = False) -> dict:
    """Stage D with the molecular clutch (INITIAL coupling).  Each released cell moves
    under the reaction of its ACTUAL clutch traction (not the nominal pull): when a
    cell's clutches slip, its outward reaction drops, so clutch failure feeds directly
    into motility (invasion / escape).  Grip sites are fixed per-cell sectors; their
    gripped fibre is re-selected as cells move (clutch state persists on the sector
    slot -- an approximation, flagged as initial).  Determinism via counter RNG."""

    network, centers, gap_radius, report = make_organoid(cfg, seed=seed)
    centers = centers.copy(); centers0 = centers.copy()
    organoid_center = np.zeros(2)
    stepper = OrganoidStepper(network, centers)
    reach = cfg.cell_radius + cfg.contact_width + 2.0
    n_sec = cfg.n_contact_sectors
    candidates = cell_candidate_fibers(network, centers, reach)
    patches, site_centers = organoid_clutch_patches(network, centers, cfg, candidates)
    S = len(patches)
    state = _new_clutch_state(S, cfg)
    substrate = np.zeros(S)
    site_force = np.zeros(S)

    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    contact_every = max(1, int(round(cfg.contact_update_interval / cfg.dt)))
    frames: list[dict] = []
    bead_snaps: list[np.ndarray] = []
    cell_snaps: list[np.ndarray] = []

    for step in range(nsteps + 1):
        time = step * cfg.dt
        _, site_force, breaks, binds, site_fail = _clutch_step(cfg, state, substrate, step)
        active = _project_site_forces(network, patches, site_force)
        # per-cell reaction = -(actual clutch traction it applies) = outward invasion pull
        reaction = np.zeros((len(centers), 2))
        for s, patch in enumerate(patches):
            if patch is None or site_force[s] <= 0.0:
                continue
            reaction[s // n_sec] -= site_force[s] * patch.normal_in

        if step % every == 0:
            prof = radial_alignment_profile(network, organoid_center)
            frames.append({
                "time": time,
                "global_radial_order": prof["global_radial_order"],
                "mean_cell_radial_disp": float(np.mean(
                    np.linalg.norm(centers, axis=1) - np.linalg.norm(centers0, axis=1))),
                "cell_spread": float(np.mean(np.linalg.norm(centers, axis=1))),
                "max_cell_disp": float(np.max(np.linalg.norm(centers - centers0, axis=1))),
                "bound_fraction": float(state.bound.mean()),
                "cumulative_slips": int(state.cumulative_slips),
                "total_traction": float(site_force.sum()),
            })
            if snapshots:
                bead_snaps.append(network.r.copy()); cell_snaps.append(centers.copy())

        if step == nsteps:
            break
        # move cells under clutch reaction + cell-cell adhesion (overdamped, capped)
        f_cell = reaction + cell_cell_forces(centers, cfg)
        v_cell = f_cell / cfg.cell_drag
        speed = np.linalg.norm(v_cell, axis=1)
        over = speed > cfg.max_cell_speed
        if np.any(over):
            v_cell[over] *= (cfg.max_cell_speed / speed[over])[:, None]
        centers += cfg.dt * v_cell
        stepper.centers[:] = centers
        site_centers = np.repeat(centers, n_sec, axis=0)

        stepper.step(active, cfg.dt)
        substrate = _clutch_substrate_speeds(network, stepper.velocity, patches, site_centers)
        if step and step % contact_every == 0:      # cells moved -> re-select gripped fibres
            candidates = cell_candidate_fibers(network, centers, reach)
            patches, site_centers = organoid_clutch_patches(network, centers, cfg, candidates)

    return {
        "config": asdict(cfg), "centers0": centers0, "centers_final": centers,
        "organoid_center": organoid_center, "connectivity": report,
        "n_cells": len(centers), "n_beads": len(network.r), "n_fibers": len(network.fibers),
        "n_clutch_sites": S, "clutch_mode": cfg.clutch_mode,
        "cumulative_slips": int(state.cumulative_slips),
        "cumulative_site_failures": int(state.cumulative_site_failures),
        "frames": frames,
        "bead_snapshots": np.asarray(bead_snaps) if snapshots else None,
        "cell_snapshots": np.asarray(cell_snaps) if snapshots else None,
        "final_positions": network.r.copy(), "initial_positions": network.r0.copy(),
        "edges": network.edges.copy(),
        "crosslinks": [(int(x.edge_a), float(x.alpha_a), int(x.edge_b), float(x.alpha_b))
                       for x in network.crosslinks],
    }


def per_cell_reaction(per_cell: list, total_force: float) -> np.ndarray:
    """Reaction force on each cell = -(traction it applies to collagen).

    A cell reels its gripped fibres inward, so the reaction points OUTWARD toward the
    ECM it grips -- a grappling pull that drives the cell into the matrix (invasion).
    """

    reac = np.zeros((len(per_cell), 2))
    for c, patches in enumerate(per_cell):
        for p in patches:
            reac[c] -= total_force * p.weight * p.normal_in
    return reac


def run_organoid_invasion(cfg: OrganoidConfig = OrganoidConfig(), seed=None,
                          snapshots: bool = False) -> dict:
    """Stage D: release the cells.  Each cell translates under the reaction of its own
    grip-and-reel traction (pulls it toward the ECM = outward invasion) plus cell-cell
    adhesion.  Strong adhesion -> the front stays cohesive (collective); weak adhesion
    -> cells escape singly.  Returns invasion metrics + optional bead/cell snapshots.
    With ``cfg.clutch_dynamics`` the reaction comes from the ACTUAL molecular-clutch
    traction (slip -> weaker pull -> altered motility); see
    :func:`_run_invasion_with_clutch`.
    """
    if cfg.clutch_dynamics:
        return _run_invasion_with_clutch(cfg, seed=seed, snapshots=snapshots)

    network, centers, gap_radius, report = make_organoid(cfg, seed=seed)
    centers = centers.copy()
    centers0 = centers.copy()
    organoid_center = np.zeros(2)
    stepper = OrganoidStepper(network, centers)
    reach = cfg.cell_radius + cfg.contact_width + 2.0

    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    contact_every = max(1, int(round(cfg.contact_update_interval / cfg.dt)))

    def refresh():
        cand = cell_candidate_fibers(network, centers, reach)
        active_full, per_cell = organoid_active_forces(
            network, centers, cfg.total_pull_force, candidates=cand)
        return active_full, per_cell, per_cell_reaction(per_cell, cfg.total_pull_force)

    active_full, per_cell, reaction = refresh()
    frames: list[dict] = []
    bead_snaps: list[np.ndarray] = []
    cell_snaps: list[np.ndarray] = []

    for step in range(nsteps + 1):
        time = step * cfg.dt
        ramp = min(1.0, time / cfg.force_ramp_time) if cfg.force_ramp_time else 1.0
        if step and step % contact_every == 0:
            active_full, per_cell, reaction = refresh()
        active = active_full * ramp

        if step % every == 0:
            prof = radial_alignment_profile(network, organoid_center)
            radial_disp = float(np.mean(
                np.linalg.norm(centers, axis=1) - np.linalg.norm(centers0, axis=1)))
            spread = float(np.mean(np.linalg.norm(centers - organoid_center, axis=1)))
            frames.append({
                "time": time,
                "global_radial_order": prof["global_radial_order"],
                "shells": prof["shells"],
                "n_gripping_cells": int(sum(1 for p in per_cell if p)),
                "mean_cell_radial_disp": radial_disp,   # >0 = outward invasion
                "cell_spread": spread,
                "max_cell_disp": float(np.max(np.linalg.norm(centers - centers0, axis=1))),
            })
            if snapshots:
                bead_snaps.append(network.r.copy())
                cell_snaps.append(centers.copy())

        if step == nsteps:
            break

        # move cells: reaction (grip-reel) + cell-cell adhesion, overdamped + capped
        f_cell = reaction * ramp + cell_cell_forces(centers, cfg)
        v_cell = f_cell / cfg.cell_drag
        speed = np.linalg.norm(v_cell, axis=1)
        over = speed > cfg.max_cell_speed
        if np.any(over):
            v_cell[over] *= (cfg.max_cell_speed / speed[over])[:, None]
        centers += cfg.dt * v_cell
        stepper.centers[:] = centers          # repulsion follows the moved cells

        stepper.step(active, cfg.dt)

    return {
        "config": asdict(cfg),
        "centers0": centers0,
        "centers_final": centers,
        "organoid_center": organoid_center,
        "connectivity": report,
        "n_cells": len(centers),
        "n_beads": len(network.r),
        "n_fibers": len(network.fibers),
        "frames": frames,
        "bead_snapshots": np.asarray(bead_snaps) if snapshots else None,
        "cell_snapshots": np.asarray(cell_snaps) if snapshots else None,
        "final_positions": network.r.copy(),
        "initial_positions": network.r0.copy(),
        "edges": network.edges.copy(),
        "crosslinks": [(int(x.edge_a), float(x.alpha_a), int(x.edge_b), float(x.alpha_b))
                       for x in network.crosslinks],
    }


def parameter_variant(cfg: OrganoidConfig, **changes) -> OrganoidConfig:
    """Public helper for sweeps / ablations (mirrors the G2 helper)."""

    return replace(cfg, **changes)
