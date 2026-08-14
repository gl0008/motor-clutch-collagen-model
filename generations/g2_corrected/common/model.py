"""Corrected bead--spring collagen mechanics shared by Generation 2.

All lengths are micrometres, forces are nanonewtons, and time is seconds.
The collagen beads obey the professor's overdamped force balance

    zeta * dr_i/dt = F_stretch + F_bend + F_crosslink + F_repulsion + F_active.

There is deliberately no SLS element.  ``zeta`` dissipates motion into the
surrounding fluid/porous matrix; the elastic bonds and freely hinged links
store mechanical energy.  These are separate terms, not two implementations
of the same constitutive law.

The network generator makes finite 20--80 um fibres, discretises each fibre
into visible beads, and anchors *only* beads in the outer boundary band.  A
generated network is accepted only when all cell-contact fibres and at least
85% of all fibres belong to a crosslink-connected boundary component.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from collections import deque
from typing import Iterable, Optional
import math

import numpy as np


def _cross2(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


@dataclass(frozen=True)
class CollagenConfig:
    """Geometry and mechanics for the MDA-MB-231 / collagen-I baseline."""

    domain_size: float = 180.0
    cell_radius: float = 9.0
    n_fibers: int = 99
    min_fiber_length: float = 20.0
    max_fiber_length: float = 80.0
    bead_spacing: float = 0.75
    fiber_diameter: float = 0.30
    collagen_modulus_mpa: float = 32.0
    crosslink_modulus_mpa: float = 0.40
    crosslink_stiffness: float = 75.0
    compression_ratio: float = 0.10
    boundary_width: float = 2.5
    curvature_amplitude: float = 0.9
    contact_width: float = 3.0
    gaussian_sigma: float = 1.5
    contact_update_interval: float = 0.20
    protrusion_half_width_deg: float = 30.0
    total_pull_force: float = 5.0
    force_ramp_time: float = 2.0
    bead_drag: float = 120.0
    repulsion_stiffness: float = 600.0
    cell_clearance: float = 0.20
    dt: float = 0.005
    duration: float = 60.0
    sample_interval: float = 0.60
    required_connected_fraction: float = 0.85
    seed: int = 17
    generation_attempts: int = 120

    @property
    def area(self) -> float:
        return math.pi * (self.fiber_diameter / 2.0) ** 2

    @property
    def second_moment(self) -> float:
        return math.pi * self.fiber_diameter**4 / 64.0

    @property
    def modulus_nN_per_um2(self) -> float:
        # 1 MPa = 1000 nN / um^2.
        return 1000.0 * self.collagen_modulus_mpa

    @property
    def axial_rigidity(self) -> float:
        # EA, in nN.  A bond of rest length l uses k=EA/l.
        return self.modulus_nN_per_um2 * self.area

    @property
    def bending_rigidity(self) -> float:
        # EI, in nN um^2.  The discrete curvature coefficient is EI/l^3.
        return self.modulus_nN_per_um2 * self.second_moment

    def validate(self) -> None:
        if self.domain_size < 6 * self.cell_radius:
            raise ValueError("domain must be at least six cell radii wide")
        if not 0 <= self.compression_ratio <= 1:
            raise ValueError("compression_ratio must be between zero and one")
        if self.n_fibers < 16 or self.bead_spacing <= 0:
            raise ValueError("network resolution is too small")
        if self.dt <= 0 or self.sample_interval < self.dt:
            raise ValueError("invalid integration or output interval")


@dataclass
class Crosslink:
    """Freely hinged link joining two material points on different fibres."""

    edge_a: int
    alpha_a: float
    edge_b: int
    alpha_b: float
    rest_vector: np.ndarray
    stiffness: float
    kind: str = "permanent"
    created_at: float = 0.0


@dataclass
class ContactPatch:
    """Closest eligible material point on one fibre in a protrusion sector."""

    fiber: int
    edge: int
    alpha: float
    point: np.ndarray
    surface_distance: float
    weight: float
    normal_in: np.ndarray


@dataclass
class NetworkSpec:
    positions: np.ndarray
    fibers: list[list[int]]
    fixed: np.ndarray
    contact_fibers: list[int]
    seed_used: int


class Network:
    """A set of bead chains connected by material-point crosslinks."""

    def __init__(
        self,
        spec: NetworkSpec,
        cfg: CollagenConfig,
        *,
        crosslinks: Optional[list[Crosslink]] = None,
    ) -> None:
        self.cfg = cfg
        self.r = np.asarray(spec.positions, dtype=float).copy()
        self.r0 = self.r.copy()
        self.fibers = [np.asarray(ids, dtype=int) for ids in spec.fibers]
        self.fixed = np.asarray(spec.fixed, dtype=bool).copy()
        self.contact_fibers = list(spec.contact_fibers)
        self.seed_used = int(spec.seed_used)

        edges: list[tuple[int, int]] = []
        edge_fiber: list[int] = []
        triplets: list[tuple[int, int, int]] = []
        for fid, ids in enumerate(self.fibers):
            edges.extend(zip(ids[:-1], ids[1:]))
            edge_fiber.extend([fid] * (len(ids) - 1))
            triplets.extend(zip(ids[:-2], ids[1:-1], ids[2:]))
        self.edges = np.asarray(edges, dtype=int)
        self.edge_fiber = np.asarray(edge_fiber, dtype=int)
        self.triplets = np.asarray(triplets, dtype=int)

        i, j = self.edges.T
        self.l0 = np.linalg.norm(self.r[j] - self.r[i], axis=1)
        self.k_tension = cfg.axial_rigidity / np.maximum(self.l0, 1e-12)
        self.k_compression = cfg.compression_ratio * self.k_tension
        a, b, c = self.triplets.T
        self.curvature0 = self.r[a] - 2.0 * self.r[b] + self.r[c]
        local_l = 0.5 * (self.l0[:-1].mean() + self.l0[1:].mean())
        self.bend_coefficient = cfg.bending_rigidity / max(local_l, 1e-12) ** 3
        self.crosslinks = _copy_links(crosslinks) if crosslinks is not None else build_crosslinks(self)
        self.refresh_crosslink_arrays()

    def refresh_crosslink_arrays(self) -> None:
        """Cache link topology so the force calculation stays vectorised."""

        self.link_edge_a = np.asarray([x.edge_a for x in self.crosslinks], dtype=int)
        self.link_edge_b = np.asarray([x.edge_b for x in self.crosslinks], dtype=int)
        self.link_alpha_a = np.asarray([x.alpha_a for x in self.crosslinks], dtype=float)
        self.link_alpha_b = np.asarray([x.alpha_b for x in self.crosslinks], dtype=float)
        self.link_rest = np.asarray([x.rest_vector for x in self.crosslinks], dtype=float).reshape(-1, 2)
        self.link_stiffness = np.asarray([x.stiffness for x in self.crosslinks], dtype=float)

    def clone(self) -> "Network":
        spec = NetworkSpec(
            self.r0.copy(),
            [ids.tolist() for ids in self.fibers],
            self.fixed.copy(),
            list(self.contact_fibers),
            self.seed_used,
        )
        return Network(spec, self.cfg, crosslinks=self.crosslinks)

    def material_point(self, edge: int, alpha: float) -> np.ndarray:
        i, j = self.edges[edge]
        return (1.0 - alpha) * self.r[i] + alpha * self.r[j]

    def crosslink_endpoints(self, link: Crosslink) -> tuple[np.ndarray, np.ndarray]:
        return (
            self.material_point(link.edge_a, link.alpha_a),
            self.material_point(link.edge_b, link.alpha_b),
        )

    def elastic_forces(self, include_crosslinks: bool = True):
        """Evaluate forces as negative gradients of the stored energies.

        Axial bonds use ``k=EA/l0``.  Compression uses 10% of the tensile
        stiffness by default to represent microbuckling.  Bending penalises a
        change from the initially stress-free discrete curvature.  A permanent
        crosslink penalises separation of its two coincident material points.
        """

        fs = np.zeros_like(self.r)
        fb = np.zeros_like(self.r)
        fx = np.zeros_like(self.r)

        i, j = self.edges.T
        d = self.r[j] - self.r[i]
        length = np.linalg.norm(d, axis=1)
        extension = length - self.l0
        stiffness = np.where(extension >= 0.0, self.k_tension, self.k_compression)
        pair = (stiffness * extension / np.maximum(length, 1e-12))[:, None] * d
        np.add.at(fs, i, pair)
        np.add.at(fs, j, -pair)

        a, b, c = self.triplets.T
        curv = self.r[a] - 2.0 * self.r[b] + self.r[c] - self.curvature0
        kb = self.bend_coefficient
        np.add.at(fb, a, -kb * curv)
        np.add.at(fb, b, 2.0 * kb * curv)
        np.add.at(fb, c, -kb * curv)

        crosslink_energy = 0.0
        crosslink_force = np.zeros(len(self.crosslinks), dtype=float)
        if include_crosslinks and len(self.crosslinks):
            ea = self.edges[self.link_edge_a]
            eb = self.edges[self.link_edge_b]
            aa = self.link_alpha_a[:, None]
            ab = self.link_alpha_b[:, None]
            pa = (1.0 - aa) * self.r[ea[:, 0]] + aa * self.r[ea[:, 1]]
            pb = (1.0 - ab) * self.r[eb[:, 0]] + ab * self.r[eb[:, 1]]
            delta = (pb - pa) - self.link_rest
            force = self.link_stiffness[:, None] * delta
            crosslink_force = np.linalg.norm(force, axis=1)
            np.add.at(fx, ea[:, 0], (1.0 - aa) * force)
            np.add.at(fx, ea[:, 1], aa * force)
            np.add.at(fx, eb[:, 0], -(1.0 - ab) * force)
            np.add.at(fx, eb[:, 1], -ab * force)
            crosslink_energy = 0.5 * float(
                np.sum(self.link_stiffness * np.sum(delta * delta, axis=1))
            )

        energy = {
            "stretch": float(0.5 * np.sum(stiffness * extension**2)),
            "bend": float(0.5 * kb * np.sum(curv * curv)),
            "crosslink": float(crosslink_energy),
        }
        return fs, fb, fx, energy, extension / np.maximum(self.l0, 1e-12), crosslink_force

    def repulsion_forces(self, center: np.ndarray) -> np.ndarray:
        """Soft no-penetration constraint at the rigid circular cell surface."""

        radial = self.r - center
        radius = np.linalg.norm(radial, axis=1)
        penetration = np.maximum(
            0.0, self.cfg.cell_radius + self.cfg.cell_clearance - radius
        )
        return (
            self.cfg.repulsion_stiffness
            * penetration[:, None]
            * radial
            / np.maximum(radius, 1e-12)[:, None]
        )

    def advance(
        self,
        active: np.ndarray,
        center: np.ndarray,
        dt: float,
        *,
        include_crosslinks: bool = True,
    ):
        """Explicitly integrate the overdamped bead force balance."""

        fs, fb, fx, energy, strain, link_force = self.elastic_forces(include_crosslinks)
        total = fs + fb + fx + self.repulsion_forces(center) + active
        total[self.fixed] = 0.0
        velocity = total / self.cfg.bead_drag
        self.r[~self.fixed] += dt * velocity[~self.fixed]
        self.r[self.fixed] = self.r0[self.fixed]
        if not np.all(np.isfinite(self.r)):
            raise FloatingPointError("non-finite bead position; reduce dt")
        return velocity, energy, strain, link_force, total


def _copy_links(links: Iterable[Crosslink]) -> list[Crosslink]:
    return [
        Crosslink(
            x.edge_a,
            x.alpha_a,
            x.edge_b,
            x.alpha_b,
            np.asarray(x.rest_vector, dtype=float).copy(),
            x.stiffness,
            x.kind,
            x.created_at,
        )
        for x in links
    ]


def _resample_curve(points: np.ndarray, spacing: float) -> np.ndarray:
    delta = np.diff(points, axis=0)
    arc = np.r_[0.0, np.cumsum(np.linalg.norm(delta, axis=1))]
    n = max(4, int(math.ceil(arc[-1] / spacing)) + 1)
    target = np.linspace(0.0, arc[-1], n)
    return np.column_stack(
        [np.interp(target, arc, points[:, axis]) for axis in range(2)]
    )


def _bezier_contact_fiber(
    phi: float, cfg: CollagenConfig, rng: np.random.Generator
) -> np.ndarray:
    """Create a curved near-cell fibre that terminates in the outer boundary."""

    radial = np.array([math.cos(phi), math.sin(phi)])
    tangent = np.array([-radial[1], radial[0]])
    gap = rng.uniform(0.45, 1.6)
    near = (cfg.cell_radius + gap) * radial
    bend_sign = 1.0 if math.sin(3.0 * phi + 0.4) >= 0 else -1.0
    control = near + 15.0 * bend_sign * tangent + 12.0 * radial
    half = cfg.domain_size / 2.0
    # Do not silently make the fibre longer when a convergence test enlarges
    # the box.  In the 140 um baseline this reaches the boundary band; in a
    # larger box it must reach the boundary through other crosslinked fibres.
    outer_radius = min(
        half - 0.25,
        float(np.linalg.norm(near)) + 0.72 * cfg.max_fiber_length,
    )
    outer = outer_radius * radial + rng.uniform(-10.0, 10.0) * tangent
    outer = np.clip(outer, -half + 0.25, half - 0.25)
    t = np.linspace(0.0, 1.0, 240)[:, None]
    curve = (1 - t) ** 2 * near + 2 * (1 - t) * t * control + t**2 * outer
    return _resample_curve(curve, cfg.bead_spacing)


def _curved_segment(
    start: np.ndarray,
    end: np.ndarray,
    cfg: CollagenConfig,
    rng: np.random.Generator,
) -> np.ndarray:
    vector = end - start
    length = np.linalg.norm(vector)
    tangent = vector / max(length, 1e-12)
    normal = np.array([-tangent[1], tangent[0]])
    phase = rng.uniform(0.0, 2.0 * math.pi)
    t = np.linspace(0.0, 1.0, 180)
    base = start[None, :] + t[:, None] * vector
    wave = cfg.curvature_amplitude * np.sin(math.pi * t) * np.sin(2 * math.pi * t + phase)
    return _resample_curve(base + wave[:, None] * normal, cfg.bead_spacing)


def _random_fiber(
    cfg: CollagenConfig,
    rng: np.random.Generator,
    *,
    boundary_seeded: bool,
) -> Optional[np.ndarray]:
    half = cfg.domain_size / 2.0
    length = rng.uniform(cfg.min_fiber_length, cfg.max_fiber_length)
    if boundary_seeded:
        side = int(rng.integers(4))
        along = rng.uniform(-0.82 * half, 0.82 * half)
        if side == 0:
            start = np.array([-half + 0.25, along]); inward = np.array([1.0, 0.0])
        elif side == 1:
            start = np.array([half - 0.25, along]); inward = np.array([-1.0, 0.0])
        elif side == 2:
            start = np.array([along, -half + 0.25]); inward = np.array([0.0, 1.0])
        else:
            start = np.array([along, half - 0.25]); inward = np.array([0.0, -1.0])
        theta = rng.uniform(-0.75, 0.75)
        rot = np.array(
            [[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]]
        )
        end = start + length * (rot @ inward)
    else:
        center = rng.uniform(-0.55 * half, 0.55 * half, size=2)
        theta = rng.uniform(0.0, math.pi)
        direction = np.array([math.cos(theta), math.sin(theta)])
        start = center - 0.5 * length * direction
        end = center + 0.5 * length * direction
    if np.max(np.abs(end)) > half - 0.15 or np.max(np.abs(start)) > half - 0.15:
        return None
    points = _curved_segment(start, end, cfg, rng)
    if np.min(np.linalg.norm(points, axis=1)) < cfg.cell_radius + cfg.cell_clearance:
        return None
    return points


def _assemble_spec(
    cfg: CollagenConfig, rng: np.random.Generator, seed_used: int
) -> NetworkSpec:
    positions: list[np.ndarray] = []
    fibers: list[list[int]] = []
    fixed: list[bool] = []
    half = cfg.domain_size / 2.0

    def append(points: np.ndarray) -> None:
        ids: list[int] = []
        for p in points:
            ids.append(len(positions))
            positions.append(np.asarray(p, dtype=float))
            fixed.append(bool(np.max(np.abs(p)) >= half - cfg.boundary_width))
        fibers.append(ids)

    # Four contact fibres on each side make V2 and V3 use exactly the same ECM.
    contact_angles = np.deg2rad([-24, -8, 8, 24, 156, 172, 188, 204])
    for phi in contact_angles:
        append(_bezier_contact_fiber(float(phi), cfg, rng))
    contact_fibers = list(range(len(contact_angles)))

    target_boundary = max(18, int(round(0.40 * cfg.n_fibers)))
    attempts = 0
    while len(fibers) < cfg.n_fibers and attempts < 50_000:
        attempts += 1
        points = _random_fiber(
            cfg,
            rng,
            boundary_seeded=len(fibers) < target_boundary,
        )
        if points is not None:
            append(points)
    if len(fibers) != cfg.n_fibers:
        raise RuntimeError("could not construct the requested finite-fibre network")
    return NetworkSpec(
        np.asarray(positions), fibers, np.asarray(fixed), contact_fibers, seed_used
    )


def build_crosslinks(network: Network) -> list[Crosslink]:
    """Place permanent hinged links at intersections of different fibres."""

    links: list[Crosslink] = []
    starts = network.r[network.edges[:, 0]]
    vectors = network.r[network.edges[:, 1]] - starts
    for fa in range(len(network.fibers)):
        ia = np.flatnonzero(network.edge_fiber == fa)
        for fb in range(fa + 1, len(network.fibers)):
            ib = np.flatnonzero(network.edge_fiber == fb)
            p = starts[ia][:, None, :]
            r = vectors[ia][:, None, :]
            q = starts[ib][None, :, :]
            s = vectors[ib][None, :, :]
            den = _cross2(r, s)
            qp = q - p
            safe = np.abs(den) > 1e-10
            denom = np.where(safe, den, 1.0)
            ta = np.where(safe, _cross2(qp, s) / denom, -1.0)
            tb = np.where(safe, _cross2(qp, r) / denom, -1.0)
            hits = np.argwhere(
                safe
                & (ta > 1e-5)
                & (ta < 1.0 - 1e-5)
                & (tb > 1e-5)
                & (tb < 1.0 - 1e-5)
            )
            accepted: list[np.ndarray] = []
            for u, v in hits:
                point = p[u, 0] + ta[u, v] * r[u, 0]
                if any(np.linalg.norm(point - old) < 0.6 for old in accepted):
                    continue
                accepted.append(point)
                links.append(
                    Crosslink(
                        int(ia[u]),
                        float(ta[u, v]),
                        int(ib[v]),
                        float(tb[u, v]),
                        np.zeros(2),
                        network.cfg.crosslink_stiffness,
                    )
                )
    return links


def connectivity_report(network: Network) -> dict:
    """Report fibre connectivity through links to the fixed outer boundary."""

    n = len(network.fibers)
    graph = [set() for _ in range(n)]
    for link in network.crosslinks:
        fa = int(network.edge_fiber[link.edge_a])
        fb = int(network.edge_fiber[link.edge_b])
        graph[fa].add(fb)
        graph[fb].add(fa)
    boundary = {
        fid for fid, ids in enumerate(network.fibers) if np.any(network.fixed[ids])
    }
    reached = set(boundary)
    queue = deque(boundary)
    while queue:
        fid = queue.popleft()
        for other in graph[fid]:
            if other not in reached:
                reached.add(other)
                queue.append(other)
    contact_ok = all(fid in reached for fid in network.contact_fibers)
    return {
        "boundary_fibers": sorted(boundary),
        "boundary_connected_fibers": sorted(reached),
        "connected_fraction": len(reached) / max(n, 1),
        "contact_fibers_connected": contact_ok,
    }


def make_network_spec(
    cfg: CollagenConfig = CollagenConfig(), seed: Optional[int] = None
) -> NetworkSpec:
    """Generate and reject networks until the boundary-percolation gate passes."""

    cfg.validate()
    base = cfg.seed if seed is None else int(seed)
    best: tuple[float, Optional[NetworkSpec]] = (-1.0, None)
    for attempt in range(cfg.generation_attempts):
        seed_used = base + 7919 * attempt
        spec = _assemble_spec(cfg, np.random.default_rng(seed_used), seed_used)
        network = Network(spec, cfg)
        report = connectivity_report(network)
        score = float(report["connected_fraction"])
        if score > best[0]:
            best = (score, spec)
        if (
            report["contact_fibers_connected"]
            and score >= cfg.required_connected_fraction
        ):
            return spec
    raise RuntimeError(
        "network percolation gate failed; best boundary-connected fraction "
        f"was {best[0]:.3f}"
    )


def _angle_in_sector(angle: float, center: float, half_width: float) -> bool:
    return abs((angle - center + math.pi) % (2 * math.pi) - math.pi) <= half_width


def contact_patches(
    network: Network,
    center: np.ndarray,
    angle_deg: float,
    half_width_deg: Optional[float] = None,
) -> list[ContactPatch]:
    """Hard contact selection followed by Gaussian force weighting.

    Only the closest material point on each fibre can enter the cell-surface
    shell.  Eligible points receive ``exp(-d_surface^2/sigma^2)`` weights,
    normalised so that the distributed vectors sum to the requested pull.
    """

    center = np.asarray(center, dtype=float)
    half = math.radians(
        network.cfg.protrusion_half_width_deg
        if half_width_deg is None
        else half_width_deg
    )
    target = math.radians(angle_deg)
    candidates: list[ContactPatch] = []
    for fid in range(len(network.fibers)):
        best: Optional[ContactPatch] = None
        edge_ids = np.flatnonzero(network.edge_fiber == fid)
        pairs = network.edges[edge_ids]
        a = network.r[pairs[:, 0]]
        d = network.r[pairs[:, 1]] - a
        alpha = np.clip(
            np.sum((center - a) * d, axis=1)
            / np.maximum(np.sum(d * d, axis=1), 1e-12),
            0.0,
            1.0,
        )
        point = a + alpha[:, None] * d
        radial = point - center
        radius = np.linalg.norm(radial, axis=1)
        surface = radius - network.cfg.cell_radius
        angle = np.arctan2(radial[:, 1], radial[:, 0])
        angular_delta = np.abs((angle - target + math.pi) % (2 * math.pi) - math.pi)
        eligible = (
            (surface >= 0.0)
            & (surface <= network.cfg.contact_width)
            & (angular_delta <= half)
        )
        if np.any(eligible):
            local = int(np.flatnonzero(eligible)[np.argmin(surface[eligible])])
            best = ContactPatch(
                fid,
                int(edge_ids[local]),
                float(alpha[local]),
                point[local].copy(),
                float(surface[local]),
                0.0,
                -radial[local] / max(float(radius[local]), 1e-12),
            )
        if best is not None:
            candidates.append(best)
    if not candidates:
        return []
    raw = np.exp(
        -np.square([x.surface_distance for x in candidates])
        / network.cfg.gaussian_sigma**2
    )
    raw /= np.sum(raw)
    for patch, weight in zip(candidates, raw):
        patch.weight = float(weight)
    return candidates


def active_forces(
    network: Network,
    center: np.ndarray,
    total_force: float,
    angle_deg: float,
    half_width_deg: Optional[float] = None,
):
    patches = contact_patches(network, center, angle_deg, half_width_deg)
    nodal = np.zeros_like(network.r)
    vectors: list[np.ndarray] = []
    for patch in patches:
        force = total_force * patch.weight * patch.normal_in
        vectors.append(force)
        i, j = network.edges[patch.edge]
        nodal[i] += (1.0 - patch.alpha) * force
        nodal[j] += patch.alpha * force
    nodal[network.fixed] = 0.0
    return nodal, patches, vectors


def forces_from_patches(
    network: Network, patches: list[ContactPatch], total_force: float
):
    """Redistribute a new force magnitude over cached material contacts."""

    nodal = np.zeros_like(network.r)
    vectors: list[np.ndarray] = []
    for patch in patches:
        force = total_force * patch.weight * patch.normal_in
        vectors.append(force)
        i, j = network.edges[patch.edge]
        nodal[i] += (1.0 - patch.alpha) * force
        nodal[j] += patch.alpha * force
    nodal[network.fixed] = 0.0
    return nodal, vectors


def _radial_metrics(network: Network, center: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    i, j = network.edges.T
    seg = network.r[j] - network.r[i]
    tangent = seg / np.maximum(np.linalg.norm(seg, axis=1), 1e-12)[:, None]
    midpoint = 0.5 * (network.r[i] + network.r[j])
    radial = midpoint - center
    radius = np.linalg.norm(radial, axis=1)
    er = radial / np.maximum(radius, 1e-12)[:, None]
    order = 2.0 * np.square(np.sum(tangent * er, axis=1)) - 1.0
    surface = radius - network.cfg.cell_radius
    shells = ((0.0, 10.0), (10.0, 25.0), (25.0, 50.0))
    align, disp = [], []
    bead_surface0 = np.linalg.norm(network.r0 - center, axis=1) - network.cfg.cell_radius
    bead_displacement = np.linalg.norm(network.r - network.r0, axis=1)
    for low, high in shells:
        emask = (surface >= low) & (surface < high)
        bmask = (bead_surface0 >= low) & (bead_surface0 < high)
        align.append(float(np.mean(order[emask])) if np.any(emask) else 0.0)
        disp.append(float(np.mean(bead_displacement[bmask])) if np.any(bmask) else 0.0)
    return np.asarray(align), np.asarray(disp)


def _record_frame(
    out: dict,
    network: Network,
    center: np.ndarray,
    time: float,
    patches: list[ContactPatch],
    vectors: list[np.ndarray],
    velocity: np.ndarray,
    energy: dict,
    strain: np.ndarray,
    link_force: np.ndarray,
    active: np.ndarray,
) -> None:
    align, shell_disp = _radial_metrics(network, center)
    out["time"].append(float(time))
    out["positions"].append(network.r.copy())
    out["bond_strain"].append(strain.copy())
    out["crosslink_force"].append(link_force.copy())
    out["contact_points"].append([x.point.tolist() for x in patches])
    out["contact_weights"].append([x.weight for x in patches])
    out["contact_forces"].append([np.asarray(v).tolist() for v in vectors])
    out["contact_fibers"].append([x.fiber for x in patches])
    out["alignment_by_shell"].append(align)
    out["displacement_by_shell"].append(shell_disp)
    disp = network.r - network.r0
    mobile = ~network.fixed
    out["rms_displacement"].append(
        float(np.sqrt(np.mean(np.sum(disp[mobile] ** 2, axis=1))))
    )
    out["stretch_energy"].append(energy["stretch"])
    out["bend_energy"].append(energy["bend"])
    out["crosslink_energy"].append(energy["crosslink"])
    out["drag_power"].append(float(network.cfg.bead_drag * np.sum(velocity**2)))
    out["active_work_rate"].append(float(np.sum(active * velocity)))
    out["min_cell_gap"].append(
        float(np.min(np.linalg.norm(network.r - center, axis=1)) - network.cfg.cell_radius)
    )


def _empty_output() -> dict:
    names = (
        "time",
        "positions",
        "bond_strain",
        "crosslink_force",
        "contact_points",
        "contact_weights",
        "contact_forces",
        "contact_fibers",
        "alignment_by_shell",
        "displacement_by_shell",
        "rms_displacement",
        "stretch_energy",
        "bend_energy",
        "crosslink_energy",
        "drag_power",
        "active_work_rate",
        "min_cell_gap",
    )
    return {name: [] for name in names}


def _finalise_output(out: dict, network: Network, include_crosslinks: bool) -> dict:
    numeric = (
        "time",
        "positions",
        "bond_strain",
        "alignment_by_shell",
        "displacement_by_shell",
        "rms_displacement",
        "stretch_energy",
        "bend_energy",
        "crosslink_energy",
        "drag_power",
        "active_work_rate",
        "min_cell_gap",
    )
    for name in numeric:
        out[name] = np.asarray(out[name])
    # Every frame has the same permanent-link count in V2/V3.
    out["crosslink_force"] = np.asarray(out["crosslink_force"])
    out.update(
        {
            "initial_positions": network.r0.copy(),
            "fibers": [ids.tolist() for ids in network.fibers],
            "edges": network.edges.copy(),
            "edge_fiber": network.edge_fiber.copy(),
            "fixed": network.fixed.copy(),
            "crosslinks": [asdict(x) for x in network.crosslinks],
            "config": asdict(network.cfg),
            "derived": {
                "axial_rigidity_nN": network.cfg.axial_rigidity,
                "bending_rigidity_nN_um2": network.cfg.bending_rigidity,
                "bead_count": len(network.r),
                "bond_count": len(network.edges),
                "crosslink_count": len(network.crosslinks),
            },
            "connectivity": connectivity_report(network),
            "include_crosslinks": bool(include_crosslinks),
            "seed_used": network.seed_used,
        }
    )
    return out


def run_fixed_pull(
    cfg: CollagenConfig = CollagenConfig(),
    *,
    spec: Optional[NetworkSpec] = None,
    include_crosslinks: bool = True,
    angle_deg: float = 0.0,
) -> dict:
    """Run one fixed-cell pull and return animation plus validation fields."""

    if spec is None:
        spec = make_network_spec(cfg)
    network = Network(spec, cfg)
    out = _empty_output()
    center = np.zeros(2)
    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    velocity = np.zeros_like(network.r)
    energy = {"stretch": 0.0, "bend": 0.0, "crosslink": 0.0}
    strain = np.zeros(len(network.edges))
    link_force = np.zeros(len(network.crosslinks))
    contact_every = max(1, int(round(cfg.contact_update_interval / cfg.dt)))
    patches = contact_patches(network, center, angle_deg)

    for step in range(nsteps + 1):
        time = step * cfg.dt
        ramp = min(1.0, time / cfg.force_ramp_time) if cfg.force_ramp_time else 1.0
        if step and step % contact_every == 0:
            patches = contact_patches(network, center, angle_deg)
        active, vectors = forces_from_patches(
            network, patches, cfg.total_pull_force * ramp
        )
        if step % every == 0:
            _record_frame(
                out,
                network,
                center,
                time,
                patches,
                vectors,
                velocity,
                energy,
                strain,
                link_force,
                active,
            )
        if step == nsteps:
            break
        velocity, energy, strain, link_force, _ = network.advance(
            active, center, cfg.dt, include_crosslinks=include_crosslinks
        )
    return _finalise_output(out, network, include_crosslinks)


def run_fixed_pull_pair(cfg: CollagenConfig = CollagenConfig()) -> dict:
    """Synchronized same-network comparison: isolated vs crosslinked fibres."""

    spec = make_network_spec(cfg)
    return {
        "without_crosslinks": run_fixed_pull(
            cfg, spec=spec, include_crosslinks=False, angle_deg=0.0
        ),
        "with_crosslinks": run_fixed_pull(
            cfg, spec=spec, include_crosslinks=True, angle_deg=0.0
        ),
        "config": asdict(cfg),
        "comparison": "same network, same force, crosslink term only",
    }


def parameter_variant(cfg: CollagenConfig, **changes) -> CollagenConfig:
    """Small public helper used by sensitivity and convergence scripts."""

    return replace(cfg, **changes)
