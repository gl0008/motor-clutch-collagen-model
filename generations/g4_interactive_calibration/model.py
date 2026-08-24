"""Generation 4: calibrated elastic collagen, indirect realignment and motion.

The four stages deliberately add one question at a time:

``G4A`` fixed cell + deterministic pull, with tunable elastic ECM parameters;
``G4B`` the same experiment classified by crosslink-graph distance;
``G4C`` the same fixed-cell ECM with stochastic Bell slip clutches; and
``G4D`` exactly G4C with rigid-cell translation and rotation released.

Units are micrometre (um), nanonewton (nN) and second (s).  Collagen follows

    zeta_b dr_i/dt = F_stretch + F_bend + F_crosslink + F_repulsion + F_active.

Only material points in the cell-surface contact shell receive ``F_active``.
The Gaussian weights forces *inside that eligible shell*; it is never applied
to the whole network.  A non-contact fibre can therefore move only through an
elastic crosslink path (or through excluded-volume contact), which is the G4B
mechanism to be tested rather than assumed from the picture.
"""

from __future__ import annotations

from collections import deque
from dataclasses import asdict, dataclass, replace
from pathlib import Path
import math
import sys
from typing import Iterable

import numpy as np

HERE = Path(__file__).resolve().parent
G2_ROOT = HERE.parent / "g2_corrected"
if str(G2_ROOT) not in sys.path:
    sys.path.insert(0, str(G2_ROOT))

from common.model import (  # noqa: E402
    CollagenConfig,
    ContactPatch,
    Network,
    NetworkSpec,
    _bezier_contact_fiber,
    _random_fiber,
)


@dataclass(frozen=True)
class G4Config(CollagenConfig):
    """Configuration for G4A-D.

    The defaults are a mechanism-calibration starting point, not a parameter
    fit to one tumour.  Website presets vary one coefficient at a time around
    this baseline so a visual change can be assigned to a known cause.
    """

    # G2-scale field, finite collagen fibres and outer-boundary-only anchoring.
    domain_size: float = 180.0
    cell_radius: float = 10.0
    n_fibers: int = 99
    min_fiber_length: float = 20.0
    max_fiber_length: float = 80.0
    bead_spacing: float = 1.0
    boundary_width: float = 2.5
    required_connected_fraction: float = 0.0  # report connectivity; do not hide low-p cases
    generation_attempts: int = 1

    # Axial modulus and geometry follow the softened G3 calibration fixture.
    # Bending is separated from axial stiffness by an explicit multiplier.
    collagen_modulus_mpa: float = 3.0
    bending_multiplier: float = 0.25
    compression_ratio: float = 0.10
    crosslink_probability: float = 0.35
    crosslink_stiffness: float = 10.0

    # Hybrid surface contact: hard eligibility first, Gaussian weights second.
    contact_width: float = 3.0
    gaussian_sigma: float = 1.5
    n_contact_sectors: int = 12
    total_pull_force: float = 24.0
    force_ramp_time: float = 10.0

    # Drag controls the approach time, not the elastic equilibrium.
    bead_drag: float = 180.0
    dt: float = 0.05
    duration: float = 180.0
    sample_interval: float = 3.0

    # Effective, coarse-grained motor-clutch bundle at every surface sector.
    n_clutches_per_site: int = 12
    clutch_stiffness: float = 2.0       # nN / um
    clutch_on_rate: float = 0.055       # 1 / s
    clutch_off_rate0: float = 0.018     # 1 / s
    bell_force: float = 1.5             # nN
    unloaded_actin_speed: float = 0.025 # um / s
    motor_stall_per_site: float = 8.0   # nN
    contact_update_interval: float = 0.50

    # G4D is G4C plus these released rigid-body degrees of freedom.
    cell_drag: float = 600.0            # nN s / um; keeps 12-site path speed in G2 target range
    rotational_drag_factor: float = 1.0
    max_cell_speed: float = 0.012       # um / s safety/biological guard
    max_angular_speed: float = 0.002     # rad / s guard
    seed: int = 41

    def validate(self) -> None:
        super().validate()
        if not 0.0 <= self.crosslink_probability <= 1.0:
            raise ValueError("crosslink_probability must be in [0, 1]")
        if self.bending_multiplier < 0.0 or self.bead_drag <= 0.0:
            raise ValueError("bending multiplier and bead drag must be non-negative")
        if self.n_contact_sectors < 4 or self.n_clutches_per_site < 1:
            raise ValueError("too few contact sectors or clutches")


def _assemble_g4_spec(cfg: G4Config, seed: int) -> NetworkSpec:
    """Make the 99-fibre G2 geometry with direct candidates all around the cell."""

    rng = np.random.default_rng(seed)
    half = cfg.domain_size / 2.0
    positions: list[np.ndarray] = []
    fibers: list[list[int]] = []
    fixed: list[bool] = []

    def append(points: np.ndarray) -> None:
        ids: list[int] = []
        for point in points:
            ids.append(len(positions))
            positions.append(np.asarray(point, dtype=float))
            # Only beads in the outer boundary band are anchored.
            fixed.append(bool(np.max(np.abs(point)) >= half - cfg.boundary_width))
        fibers.append(ids)

    angles = np.linspace(0.0, 2.0 * math.pi, cfg.n_contact_sectors, endpoint=False)
    for phi in angles:
        append(_bezier_contact_fiber(float(phi), cfg, rng))
    contact_fibers = list(range(len(angles)))

    target_boundary = max(20, int(round(0.42 * cfg.n_fibers)))
    attempts = 0
    while len(fibers) < cfg.n_fibers and attempts < 80_000:
        attempts += 1
        points = _random_fiber(cfg, rng, boundary_seeded=len(fibers) < target_boundary)
        if points is not None:
            append(points)
    if len(fibers) != cfg.n_fibers:
        raise RuntimeError("could not construct G4 finite-fibre network")
    return NetworkSpec(np.asarray(positions), fibers, np.asarray(fixed), contact_fibers, seed)


def build_g4_network(
    cfg: G4Config = G4Config(), *, spec: NetworkSpec | None = None
) -> tuple[Network, NetworkSpec, dict]:
    """Build nested probabilistic permanent links on one reusable geometry.

    Every geometric intersection gets a deterministic U(0,1) mark.  A link is
    present iff ``mark < p_x``.  Raising ``p_x`` therefore adds links without
    changing or reshuffling the links already present at lower probability.
    """

    cfg.validate()
    if spec is None:
        spec = _assemble_g4_spec(cfg, cfg.seed)
    full = Network(spec, cfg)
    marks = np.random.default_rng(spec.seed_used + 101).random(len(full.crosslinks))
    full.crosslinks = [x for x, u in zip(full.crosslinks, marks) if u < cfg.crosslink_probability]
    full.refresh_crosslink_arrays()
    # Tune flexural rigidity without silently changing axial rigidity EA.
    full.bend_coefficient *= cfg.bending_multiplier
    graph = fiber_graph(full)
    boundary = {fid for fid, ids in enumerate(full.fibers) if np.any(full.fixed[ids])}
    reached = _reachable(graph, boundary)
    report = {
        "candidate_crosslinks": int(len(marks)),
        "crosslinks": int(len(full.crosslinks)),
        "boundary_connected_fraction": len(reached) / max(len(full.fibers), 1),
        "boundary_fibers": sorted(boundary),
    }
    return full, spec, report


def fiber_graph(network: Network) -> list[set[int]]:
    graph = [set() for _ in network.fibers]
    for link in network.crosslinks:
        fa = int(network.edge_fiber[link.edge_a])
        fb = int(network.edge_fiber[link.edge_b])
        graph[fa].add(fb)
        graph[fb].add(fa)
    return graph


def _reachable(graph: list[set[int]], sources: Iterable[int]) -> set[int]:
    reached = set(int(x) for x in sources)
    queue = deque(reached)
    while queue:
        here = queue.popleft()
        for other in graph[here]:
            if other not in reached:
                reached.add(other)
                queue.append(other)
    return reached


def graph_distance_from_direct(network: Network, direct: Iterable[int]) -> np.ndarray:
    """0=direct, 1=one link away, 2=two or more links, -1=unconnected."""

    graph = fiber_graph(network)
    distance = np.full(len(graph), -1, dtype=int)
    queue = deque()
    for fid in direct:
        distance[int(fid)] = 0
        queue.append(int(fid))
    while queue:
        here = queue.popleft()
        for other in graph[here]:
            if distance[other] < 0:
                distance[other] = distance[here] + 1
                queue.append(other)
    return distance


def direct_contact_patches(network: Network, center: np.ndarray) -> list[ContactPatch]:
    """Hard shell selection plus normalized Gaussian weights on direct fibres.

    The closest continuous point on each designated contact fibre is eligible
    only if it is 0--``contact_width`` um outside the cell surface.  No other
    fibre is assigned active force here.
    """

    center = np.asarray(center, dtype=float)
    candidates: list[ContactPatch] = []
    for fid in network.contact_fibers:
        edge_ids = np.flatnonzero(network.edge_fiber == fid)
        pairs = network.edges[edge_ids]
        a = network.r[pairs[:, 0]]
        d = network.r[pairs[:, 1]] - a
        alpha = np.clip(
            np.sum((center - a) * d, axis=1) / np.maximum(np.sum(d * d, axis=1), 1e-12),
            0.0,
            1.0,
        )
        points = a + alpha[:, None] * d
        radial = points - center
        radius = np.linalg.norm(radial, axis=1)
        surface = radius - network.cfg.cell_radius
        eligible = (surface >= 0.0) & (surface <= network.cfg.contact_width)
        if not np.any(eligible):
            continue
        local = int(np.flatnonzero(eligible)[np.argmin(surface[eligible])])
        candidates.append(ContactPatch(
            fid,
            int(edge_ids[local]),
            float(alpha[local]),
            points[local].copy(),
            float(surface[local]),
            0.0,
            -radial[local] / max(float(radius[local]), 1e-12),
        ))
    if not candidates:
        return []
    raw = np.exp(-np.square([p.surface_distance for p in candidates]) / network.cfg.gaussian_sigma**2)
    raw /= raw.sum()
    for patch, weight in zip(candidates, raw):
        patch.weight = float(weight)
    return candidates


def patch_point(network: Network, patch: ContactPatch) -> np.ndarray:
    i, j = network.edges[patch.edge]
    return (1.0 - patch.alpha) * network.r[i] + patch.alpha * network.r[j]


def active_from_patches(
    network: Network,
    patches: list[ContactPatch],
    scalar_forces: np.ndarray | float,
    center: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map material-point pulls to beads with linear shape functions."""

    active = np.zeros_like(network.r)
    if np.isscalar(scalar_forces):
        magnitudes = np.asarray([float(scalar_forces) * p.weight for p in patches])
    else:
        magnitudes = np.asarray(scalar_forces, dtype=float)
    points, vectors = [], []
    for patch, magnitude in zip(patches, magnitudes):
        point = patch_point(network, patch)
        inward = center - point
        inward /= max(float(np.linalg.norm(inward)), 1e-12)
        force = float(magnitude) * inward
        i, j = network.edges[patch.edge]
        active[i] += (1.0 - patch.alpha) * force
        active[j] += patch.alpha * force
        points.append(point)
        vectors.append(force)
    return active, np.asarray(points), np.asarray(vectors)


def _scatter(n: int, idx: np.ndarray, values: np.ndarray) -> np.ndarray:
    out = np.empty((n, 2))
    out[:, 0] = np.bincount(idx, values[:, 0], minlength=n)
    out[:, 1] = np.bincount(idx, values[:, 1], minlength=n)
    return out


def fast_advance(network: Network, active: np.ndarray, center: np.ndarray, dt: float):
    """Integrate the implemented overdamped force balance for one Euler step."""

    r = network.r
    n = len(r)
    i, j = network.edges.T
    d = r[j] - r[i]
    length = np.linalg.norm(d, axis=1)
    extension = length - network.l0
    stiffness = np.where(extension >= 0.0, network.k_tension, network.k_compression)
    pair = (stiffness * extension / np.maximum(length, 1e-12))[:, None] * d
    stretch = _scatter(n, i, pair) - _scatter(n, j, pair)

    a, b, c = network.triplets.T
    curvature = r[a] - 2.0 * r[b] + r[c] - network.curvature0
    kb = network.bend_coefficient
    bend = (
        _scatter(n, a, -kb * curvature)
        + _scatter(n, b, 2.0 * kb * curvature)
        + _scatter(n, c, -kb * curvature)
    )

    link_force = np.zeros((n, 2))
    link_magnitude = np.zeros(len(network.crosslinks))
    if len(network.crosslinks):
        ea = network.edges[network.link_edge_a]
        eb = network.edges[network.link_edge_b]
        aa = network.link_alpha_a[:, None]
        ab = network.link_alpha_b[:, None]
        pa = (1.0 - aa) * r[ea[:, 0]] + aa * r[ea[:, 1]]
        pb = (1.0 - ab) * r[eb[:, 0]] + ab * r[eb[:, 1]]
        force = network.link_stiffness[:, None] * ((pb - pa) - network.link_rest)
        link_magnitude = np.linalg.norm(force, axis=1)
        link_force = (
            _scatter(n, ea[:, 0], (1.0 - aa) * force)
            + _scatter(n, ea[:, 1], aa * force)
            - _scatter(n, eb[:, 0], (1.0 - ab) * force)
            - _scatter(n, eb[:, 1], ab * force)
        )

    radial = r - center
    radius = np.linalg.norm(radial, axis=1)
    penetration = np.maximum(0.0, network.cfg.cell_radius + network.cfg.cell_clearance - radius)
    repulsion = (
        network.cfg.repulsion_stiffness * penetration[:, None] * radial
        / np.maximum(radius, 1e-12)[:, None]
    )
    total = stretch + bend + link_force + repulsion + active
    total[network.fixed] = 0.0
    velocity = total / network.cfg.bead_drag
    r[~network.fixed] += dt * velocity[~network.fixed]
    r[network.fixed] = network.r0[network.fixed]
    if not np.all(np.isfinite(r)):
        raise FloatingPointError("non-finite bead position; reduce dt")
    return velocity, total, link_magnitude


def radial_order_by_fiber(network: Network, center: np.ndarray) -> np.ndarray:
    values = np.zeros(len(network.fibers))
    for fid in range(len(network.fibers)):
        edge_ids = np.flatnonzero(network.edge_fiber == fid)
        i, j = network.edges[edge_ids].T
        segment = network.r[j] - network.r[i]
        tangent = segment / np.maximum(np.linalg.norm(segment, axis=1), 1e-12)[:, None]
        midpoint = 0.5 * (network.r[i] + network.r[j])
        radial = midpoint - center
        radial /= np.maximum(np.linalg.norm(radial, axis=1), 1e-12)[:, None]
        values[fid] = float(np.mean(2.0 * np.square(np.sum(tangent * radial, axis=1)) - 1.0))
    return values


def displacement_by_graph_class(network: Network, distance: np.ndarray) -> dict[str, float]:
    bead_disp = np.linalg.norm(network.r - network.r0, axis=1)
    fiber_disp = np.asarray([float(np.mean(bead_disp[ids])) for ids in network.fibers])
    classes = {
        "direct": distance == 0,
        "one_hop": distance == 1,
        "two_plus": distance >= 2,
        "unconnected": distance < 0,
    }
    return {name: float(np.mean(fiber_disp[mask])) if np.any(mask) else 0.0 for name, mask in classes.items()}


def _frame(network: Network, center: np.ndarray, patches: list[ContactPatch], vectors: np.ndarray,
           distance: np.ndarray, time: float, *, traction: float = 0.0,
           bound: int = 0, slips: int = 0, theta: float = 0.0,
           traction_drops: int = 0, recoil_events: int = 0) -> dict:
    points = np.asarray([patch_point(network, p) for p in patches]) if patches else np.empty((0, 2))
    return {
        "time": float(time),
        "positions": network.r.copy(),
        "center": np.asarray(center).copy(),
        "theta": float(theta),
        "contact_points": points,
        "contact_vectors": np.asarray(vectors).copy(),
        "traction": float(traction),
        "bound": int(bound),
        "slips": int(slips),
        "traction_drops": int(traction_drops),
        "recoil_events": int(recoil_events),
        "graph_displacement": displacement_by_graph_class(network, distance),
        "radial_order": float(np.mean(radial_order_by_fiber(network, center))),
    }


def run_elastic(
    cfg: G4Config = G4Config(), *, spec: NetworkSpec | None = None
) -> dict:
    """G4A/B: fixed cell, deterministic ramped pull, permanent elastic links."""

    network, spec, report = build_g4_network(cfg, spec=spec)
    center = np.zeros(2)
    patches = direct_contact_patches(network, center)
    direct = [p.fiber for p in patches]
    distance = graph_distance_from_direct(network, direct)
    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    velocity = np.zeros_like(network.r)
    frames: list[dict] = []
    for step in range(nsteps + 1):
        time = step * cfg.dt
        ramp = min(1.0, time / max(cfg.force_ramp_time, cfg.dt))
        active, _, vectors = active_from_patches(
            network, patches, ramp * cfg.total_pull_force, center
        )
        if step % every == 0:
            frames.append(_frame(
                network, center, patches, vectors, distance, time,
                traction=ramp * cfg.total_pull_force,
            ))
        if step == nsteps:
            break
        velocity, _, _ = fast_advance(network, active, center, cfg.dt)
    initial_order = radial_order_by_fiber(Network(spec, cfg), center)
    final_order = radial_order_by_fiber(network, center)
    return {
        "config": asdict(cfg),
        "spec": spec,
        "network": network,
        "frames": frames,
        "report": report,
        "direct_fibers": direct,
        "graph_distance": distance,
        "initial_order": initial_order,
        "final_order": final_order,
        "final_graph_displacement": displacement_by_graph_class(network, distance),
        "max_displacement": float(np.max(np.linalg.norm(network.r - network.r0, axis=1))),
        "mean_speed_final": float(np.mean(np.linalg.norm(velocity, axis=1))),
    }


@dataclass
class ClutchState:
    patches: list[ContactPatch]
    bound: np.ndarray
    extension: np.ndarray
    theta: float = 0.0
    cumulative_slips: int = 0


def bell_off_rate(force: np.ndarray | float, cfg: G4Config) -> np.ndarray:
    """Bell (1978) slip-bond hazard; force raises risk continuously."""

    return cfg.clutch_off_rate0 * np.exp(np.minimum(np.asarray(force) / cfg.bell_force, 12.0))


def bell_summary(cfg: G4Config = G4Config()) -> dict:
    """Human-readable default landmarks; these are not hard thresholds."""

    forces = np.asarray([0.0, cfg.bell_force, 2.0 * cfg.bell_force])
    rates = bell_off_rate(forces, cfg)
    return {
        "forces": forces,
        "extensions": forces / cfg.clutch_stiffness,
        "median_lifetimes": math.log(2.0) / rates,
    }


def _clutch_random_stream(cfg: G4Config, seed: int):
    nsteps = int(round(cfg.duration / cfg.dt))
    shape = (nsteps, cfg.n_contact_sectors, cfg.n_clutches_per_site)
    rng = np.random.default_rng(seed)
    return rng.random(shape), rng.random(shape)


def _patch_substrate_speeds(network: Network, velocity: np.ndarray,
                             patches: list[ContactPatch], center: np.ndarray) -> np.ndarray:
    out = np.zeros(len(patches))
    for index, patch in enumerate(patches):
        point = patch_point(network, patch)
        inward = center - point
        inward /= max(float(np.linalg.norm(inward)), 1e-12)
        i, j = network.edges[patch.edge]
        material_velocity = (1.0 - patch.alpha) * velocity[i] + patch.alpha * velocity[j]
        out[index] = float(material_velocity @ inward)
    return out


def _clutch_step(
    cfg: G4Config,
    state: ClutchState,
    substrate_speed: np.ndarray,
    u_on: np.ndarray,
    u_off: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load, Bell-unbind and bind the effective clutches for one time step."""

    force_before = cfg.clutch_stiffness * state.extension * state.bound
    traction_before = force_before.sum(axis=1)
    actin_speed = cfg.unloaded_actin_speed * np.maximum(
        0.0, 1.0 - traction_before / cfg.motor_stall_per_site
    )
    relative = np.maximum(0.0, actin_speed - substrate_speed)
    state.extension[state.bound] += cfg.dt * np.repeat(
        relative[:, None], cfg.n_clutches_per_site, axis=1
    )[state.bound]

    force = cfg.clutch_stiffness * state.extension
    off_rate = bell_off_rate(force, cfg)
    breaks = state.bound & (u_off < 1.0 - np.exp(-off_rate * cfg.dt))
    state.bound[breaks] = False
    state.extension[breaks] = 0.0

    on_probability = 1.0 - math.exp(-cfg.clutch_on_rate * cfg.dt)
    binds = (~state.bound) & (u_on < on_probability)
    state.bound[binds] = True
    state.extension[binds] = 0.0
    state.cumulative_slips += int(breaks.sum())
    force = cfg.clutch_stiffness * state.extension * state.bound
    return force, force.sum(axis=1), breaks


def run_clutch(
    cfg: G4Config = G4Config(), *, spec: NetworkSpec | None = None,
    moving: bool, random_stream=None,
) -> dict:
    """G4C/D: stochastic symmetric clutches; optionally release rigid cell."""

    network, spec, report = build_g4_network(cfg, spec=spec)
    center = np.zeros(2)
    patches = direct_contact_patches(network, center)
    # The generator creates one direct candidate per sector.  Preserve a fixed
    # rectangular stochastic state even if a pathological geometry loses one.
    if len(patches) != cfg.n_contact_sectors:
        raise RuntimeError(f"expected {cfg.n_contact_sectors} direct patches; found {len(patches)}")
    direct = [p.fiber for p in patches]
    distance = graph_distance_from_direct(network, direct)
    state = ClutchState(
        patches=patches,
        bound=np.zeros((cfg.n_contact_sectors, cfg.n_clutches_per_site), dtype=bool),
        extension=np.zeros((cfg.n_contact_sectors, cfg.n_clutches_per_site)),
    )
    if random_stream is None:
        random_stream = _clutch_random_stream(cfg, cfg.seed + 3001)
    u_on, u_off = random_stream
    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    update_every = max(1, int(round(cfg.contact_update_interval / cfg.dt)))
    bead_velocity = np.zeros_like(network.r)
    substrate_speed = np.zeros(cfg.n_contact_sectors)
    traction = np.zeros(cfg.n_contact_sectors)
    vectors = np.zeros((cfg.n_contact_sectors, 2))
    frames: list[dict] = []
    slip_times: list[float] = []
    path_length = 0.0
    theta = 0.0
    traction_drop_events = 0
    recoil_events = 0

    for step in range(nsteps + 1):
        time = step * cfg.dt
        if step % every == 0:
            frames.append(_frame(
                network, center, patches, vectors, distance, time,
                traction=float(traction.sum()), bound=int(state.bound.sum()),
                slips=state.cumulative_slips, theta=theta,
                traction_drops=traction_drop_events, recoil_events=recoil_events,
            ))
        if step == nsteps:
            break
        traction_before = traction.copy()
        _, traction, breaks = _clutch_step(
            cfg, state, substrate_speed, u_on[step], u_off[step]
        )
        n_slips = int(breaks.sum())
        if n_slips:
            slip_times.extend([time] * n_slips)
            if float(traction.sum()) < float(traction_before.sum()):
                traction_drop_events += 1

        active, points, vectors = active_from_patches(network, patches, traction, center)
        reaction = -active.sum(axis=0)
        if moving:
            cell_velocity = reaction / cfg.cell_drag
            speed = float(np.linalg.norm(cell_velocity))
            if speed > cfg.max_cell_speed:
                cell_velocity *= cfg.max_cell_speed / speed
            # Reaction is applied at cell-surface anchors.  Off-axis forces can
            # therefore rotate the rigid body without prescribing a direction.
            anchors = center + cfg.cell_radius * np.column_stack([
                np.cos(np.linspace(0, 2 * math.pi, cfg.n_contact_sectors, endpoint=False) + theta),
                np.sin(np.linspace(0, 2 * math.pi, cfg.n_contact_sectors, endpoint=False) + theta),
            ])
            site_reactions = -vectors
            arms = anchors - center
            torque = float(np.sum(arms[:, 0] * site_reactions[:, 1] - arms[:, 1] * site_reactions[:, 0]))
            rotational_drag = cfg.rotational_drag_factor * cfg.cell_drag * cfg.cell_radius**2
            angular_velocity = float(np.clip(
                torque / rotational_drag, -cfg.max_angular_speed, cfg.max_angular_speed
            ))
            delta = cfg.dt * cell_velocity
            center += delta
            theta += cfg.dt * angular_velocity
            path_length += float(np.linalg.norm(delta))
        # fixed G4C computes the same reaction but constrains centre and angle.
        bead_velocity, _, _ = fast_advance(network, active, center, cfg.dt)
        substrate_speed = _patch_substrate_speeds(network, bead_velocity, patches, center)
        # Negative inward speed is outward elastic recoil.  Count a recoil only
        # on a site that actually ruptured during this same update.
        if n_slips and np.any((breaks.sum(axis=1) > 0) & (substrate_speed < 0.0)):
            recoil_events += 1

        # Re-select the closest surface material point only after a site's whole
        # clutch bundle has detached.  Bound clutches never jump between beads.
        if step and step % update_every == 0:
            fresh = direct_contact_patches(network, center)
            fresh_by_fiber = {p.fiber: p for p in fresh}
            for site, patch in enumerate(patches):
                if not state.bound[site].any() and patch.fiber in fresh_by_fiber:
                    patches[site] = fresh_by_fiber[patch.fiber]

    return {
        "config": asdict(cfg),
        "spec": spec,
        "network": network,
        "report": report,
        "frames": frames,
        "direct_fibers": direct,
        "graph_distance": distance,
        "slip_times": np.asarray(slip_times),
        "traction_drop_events": traction_drop_events,
        "recoil_events": recoil_events,
        "path_length": path_length,
        "net_displacement": float(np.linalg.norm(center)),
        "moving": moving,
        "bell_summary": bell_summary(cfg),
    }


def run_clutch_pair(cfg: G4Config = G4Config(), *, spec: NetworkSpec | None = None) -> dict:
    """Synchronized G4C fixed / G4D moving pair with identical proposed events."""

    if spec is None:
        spec = _assemble_g4_spec(cfg, cfg.seed)
    stream = _clutch_random_stream(cfg, cfg.seed + 3001)
    return {
        "fixed": run_clutch(cfg, spec=spec, moving=False, random_stream=stream),
        "moving": run_clutch(cfg, spec=spec, moving=True, random_stream=stream),
        "comparison": "same geometry, links and random stream; rigid-body constraint only",
    }


def variant(cfg: G4Config, coefficient: str, value: float) -> G4Config:
    """Build a declared one-coefficient-at-a-time G4A website preset."""

    allowed = {
        "total_pull_force", "bending_multiplier", "bead_drag",
        "crosslink_probability", "crosslink_stiffness",
    }
    if coefficient not in allowed:
        raise ValueError(f"not a G4A tuning coefficient: {coefficient}")
    return replace(cfg, **{coefficient: value})
