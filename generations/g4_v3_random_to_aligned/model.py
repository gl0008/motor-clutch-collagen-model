"""G4 v3 random-to-aligned collagen calibration.

This module deliberately does not call the G4 v1/v2 geometry generator.  That
generator inserts twelve hand-made radial fibres around the cell.  Here the
source-fibre angle is sampled uniformly on ``[0, pi)`` and a circular cell-sized
void is cut out afterwards.  Consequently radial alignment is an output metric,
not an initial condition.

The first benchmark follows the discrete-fibre logic of Abhilash et al. (2014):
the outer boundary is fixed, geometric intersections are permanent freely
hinged links, and the inner circular boundary is prescribed to contract.  It is
not yet a motor-clutch model.  Its single purpose is to test whether the elastic
network can change from bending-dominated deformation to stretching-dominated
radial fibre recruitment before additional biology is added.

Units are micrometres (um), nanonewtons (nN), and seconds (s).
"""

from __future__ import annotations

from collections import deque
from dataclasses import asdict, dataclass, replace
from pathlib import Path
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from generations.g2_corrected.common.model import (  # noqa: E402
    ContactPatch,
    Network,
    NetworkSpec,
)
from generations.g4_interactive_calibration.model import (  # noqa: E402
    active_from_patches,
    fast_advance,
    patch_point,
)
from generations.g4_v2_multiscale.model import G4V2Config  # noqa: E402


@dataclass(frozen=True)
class G4V3Config(G4V2Config):
    """Parameters inherited from G4 v2 plus explicit v3 calibration controls.

    ``n_fibers`` means the number of source line segments sampled before the
    circular void is cut.  A source segment that crosses the void becomes two
    independent physical fragments, so the realised fragment count can be
    slightly larger.

    The defaults are a mechanism-screening baseline, not a molecular fit.
    Values inherited from G4 v2 retain its single-cell length and force scales.
    """

    n_fibers: int = 240
    bead_spacing: float = 1.50
    curvature_amplitude: float = 0.0
    # Keep fragments long enough to possess at least two bending triplets.
    minimum_fragment_length: float = 4.5
    geometry_attempt_factor: int = 40

    # Literature-style prescribed-contraction benchmark.  Full intersection
    # linking is intentional here and is not the later biological baseline.
    benchmark_crosslink_probability: float = 1.0
    benchmark_ramp_time: float = 120.0
    benchmark_relax_time: float = 480.0
    benchmark_sample_interval: float = 12.0
    benchmark_contraction: float = 0.10

    # Deterministic fixed-cell traction stage.  These keep the inherited G4 v2
    # force/contact values but give the random network time to relax.
    traction_duration: float = 480.0
    traction_sample_interval: float = 12.0

    def validate(self) -> None:
        super().validate()
        if self.minimum_fragment_length < 2.0 * self.bead_spacing:
            raise ValueError("minimum_fragment_length is too short for bending")
        if not 0.0 <= self.benchmark_crosslink_probability <= 1.0:
            raise ValueError("benchmark_crosslink_probability must be in [0, 1]")
        if self.benchmark_ramp_time <= 0.0 or self.benchmark_relax_time <= 0.0:
            raise ValueError("benchmark times must be positive")
        if self.benchmark_sample_interval < self.dt:
            raise ValueError("benchmark sample interval must be at least dt")
        if self.traction_duration <= 0.0 or self.traction_sample_interval < self.dt:
            raise ValueError("invalid traction duration or sample interval")


def _clip_segment_to_square(
    start: np.ndarray, end: np.ndarray, half: float
) -> tuple[np.ndarray, np.ndarray] | None:
    """Liang--Barsky clip of a finite segment to the square simulation box."""

    direction = end - start
    t0, t1 = 0.0, 1.0
    for p, q in (
        (-direction[0], start[0] + half),
        (direction[0], half - start[0]),
        (-direction[1], start[1] + half),
        (direction[1], half - start[1]),
    ):
        if abs(float(p)) < 1e-14:
            if q < 0.0:
                return None
            continue
        ratio = float(q / p)
        if p < 0.0:
            t0 = max(t0, ratio)
        else:
            t1 = min(t1, ratio)
        if t0 > t1:
            return None
    return start + t0 * direction, start + t1 * direction


def _outside_disk_fragments(
    start: np.ndarray, end: np.ndarray, radius: float
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Remove the part inside a centred disk without changing line angle."""

    direction = end - start
    aa = float(direction @ direction)
    bb = 2.0 * float(start @ direction)
    cc = float(start @ start) - radius**2
    disc = bb * bb - 4.0 * aa * cc
    cuts = [0.0, 1.0]
    if disc > 0.0 and aa > 1e-14:
        root = math.sqrt(disc)
        for value in ((-bb - root) / (2.0 * aa), (-bb + root) / (2.0 * aa)):
            if 0.0 < value < 1.0:
                cuts.append(float(value))
    cuts = sorted(cuts)
    pieces: list[tuple[np.ndarray, np.ndarray]] = []
    for left, right in zip(cuts[:-1], cuts[1:]):
        middle = start + 0.5 * (left + right) * direction
        if float(middle @ middle) >= radius**2 - 1e-10:
            pieces.append((start + left * direction, start + right * direction))
    return pieces


def _discretise_line(start: np.ndarray, end: np.ndarray, spacing: float) -> np.ndarray:
    length = float(np.linalg.norm(end - start))
    n = max(4, int(math.ceil(length / spacing)) + 1)
    alpha = np.linspace(0.0, 1.0, n)[:, None]
    return start[None, :] + alpha * (end - start)[None, :]


def make_random_void_spec(
    cfg: G4V3Config = G4V3Config(), *, seed: int | None = None
) -> NetworkSpec:
    """Generate isotropically oriented finite fibres around a circular void.

    Sampling rule
    -------------
    1. Draw the source angle uniformly on ``[0, pi)``.
    2. Draw a midpoint uniformly over the square and a length uniformly over
       the G4 v2 20--80 um range.
    3. Clip the segment to the outer square and then cut out the circular void.

    Clipping can shorten a fibre but never rotates it.  We therefore report the
    realised global and near-cell nematic order instead of merely assuming that
    the finite sample is isotropic.
    """

    cfg.validate()
    seed_used = cfg.seed if seed is None else int(seed)
    rng = np.random.default_rng(seed_used)
    half = 0.5 * cfg.domain_size
    void_radius = cfg.cell_radius + cfg.cell_clearance
    positions: list[np.ndarray] = []
    fibers: list[list[int]] = []
    fixed: list[bool] = []
    contact_fibers: list[int] = []
    source_count = 0
    attempts = 0
    maximum_attempts = cfg.geometry_attempt_factor * cfg.n_fibers

    while source_count < cfg.n_fibers and attempts < maximum_attempts:
        attempts += 1
        theta = rng.uniform(0.0, math.pi)
        direction = np.array([math.cos(theta), math.sin(theta)])
        length = rng.uniform(cfg.min_fiber_length, cfg.max_fiber_length)
        midpoint = rng.uniform(-half, half, size=2)
        clipped = _clip_segment_to_square(
            midpoint - 0.5 * length * direction,
            midpoint + 0.5 * length * direction,
            half,
        )
        if clipped is None:
            continue
        fragments = _outside_disk_fragments(*clipped, void_radius)
        kept: list[np.ndarray] = []
        for start, end in fragments:
            if np.linalg.norm(end - start) + 1e-9 < cfg.minimum_fragment_length:
                continue
            kept.append(_discretise_line(start, end, cfg.bead_spacing))
        if not kept:
            continue

        source_count += 1
        for points in kept:
            ids = list(range(len(positions), len(positions) + len(points)))
            fid = len(fibers)
            fibers.append(ids)
            for point in points:
                positions.append(np.asarray(point, dtype=float))
                # Only the outer square is anchored.  Inner-boundary beads are
                # separately constrained only during the contraction benchmark.
                fixed.append(bool(np.max(np.abs(point)) >= half - 1e-8))
            surface = np.linalg.norm(points, axis=1) - void_radius
            if float(np.min(np.abs(surface))) < 1e-6:
                contact_fibers.append(fid)

    if source_count != cfg.n_fibers:
        raise RuntimeError(
            f"could generate only {source_count}/{cfg.n_fibers} source fibres"
        )
    return NetworkSpec(
        positions=np.asarray(positions),
        fibers=fibers,
        fixed=np.asarray(fixed, dtype=bool),
        contact_fibers=contact_fibers,
        seed_used=seed_used,
    )


def _global_nematic_order(network: Network, positions: np.ndarray | None = None) -> float:
    """Return 0 for isotropic and 1 for perfectly parallel 2D segments."""

    r = network.r if positions is None else np.asarray(positions)
    i, j = network.edges.T
    segment = r[j] - r[i]
    angle = np.arctan2(segment[:, 1], segment[:, 0])
    return float(math.hypot(np.mean(np.cos(2.0 * angle)), np.mean(np.sin(2.0 * angle))))


def _radial_order(
    network: Network,
    positions: np.ndarray | None = None,
    shells: tuple[tuple[float, float], ...] = ((0.0, 10.0), (10.0, 25.0), (25.0, 50.0)),
) -> np.ndarray:
    """Compute ``S_r=<2(t.e_r)^2-1>`` in distance-from-cell shells.

    ``S_r=+1`` is radial, ``S_r=-1`` is circumferential and ``S_r~=0`` is
    isotropic.  This signed metric prevents any displacement from being
    mislabeled as radial alignment.
    """

    r = network.r if positions is None else np.asarray(positions)
    i, j = network.edges.T
    segment = r[j] - r[i]
    tangent = segment / np.maximum(np.linalg.norm(segment, axis=1), 1e-12)[:, None]
    midpoint = 0.5 * (r[i] + r[j])
    radius = np.linalg.norm(midpoint, axis=1)
    radial = midpoint / np.maximum(radius, 1e-12)[:, None]
    value = 2.0 * np.square(np.sum(tangent * radial, axis=1)) - 1.0
    surface = radius - network.cfg.cell_radius
    out = []
    for low, high in shells:
        mask = (surface >= low) & (surface < high)
        out.append(float(np.mean(value[mask])) if np.any(mask) else float("nan"))
    return np.asarray(out)


def _inner_boundary_beads(network: Network, tolerance: float = 1e-5) -> np.ndarray:
    target = network.cfg.cell_radius + network.cfg.cell_clearance
    surface = np.abs(np.linalg.norm(network.r0, axis=1) - target)
    endpoint = np.zeros(len(network.r), dtype=bool)
    for ids in network.fibers:
        endpoint[ids[0]] = True
        endpoint[ids[-1]] = True
    return np.flatnonzero(endpoint & (surface <= tolerance))


def build_random_network(
    cfg: G4V3Config = G4V3Config(),
    *,
    spec: NetworkSpec | None = None,
    crosslink_probability: float | None = None,
) -> tuple[Network, NetworkSpec, dict]:
    """Build nested probabilistic links on the new random geometry.

    Every candidate intersection gets a deterministic random mark.  Increasing
    probability adds links without changing geometry or reshuffling links that
    already existed.  This makes cumulative stage comparisons interpretable.
    """

    cfg.validate()
    if spec is None:
        spec = make_random_void_spec(cfg)
    probability = cfg.crosslink_probability if crosslink_probability is None else float(crosslink_probability)
    if not 0.0 <= probability <= 1.0:
        raise ValueError("crosslink_probability must be in [0, 1]")
    network = Network(spec, cfg)
    candidate_count = len(network.crosslinks)
    marks = np.random.default_rng(spec.seed_used + 101).random(candidate_count)
    network.crosslinks = [
        link for link, mark in zip(network.crosslinks, marks) if mark < probability
    ]
    network.refresh_crosslink_arrays()
    network.bend_coefficient *= cfg.bending_multiplier
    inner = _inner_boundary_beads(network)
    graph = [set() for _ in network.fibers]
    for link in network.crosslinks:
        fa = int(network.edge_fiber[link.edge_a])
        fb = int(network.edge_fiber[link.edge_b])
        graph[fa].add(fb)
        graph[fb].add(fa)
    boundary_fibers = {
        fid for fid, ids in enumerate(network.fibers) if np.any(network.fixed[ids])
    }
    contact_fibers = set(int(fid) for fid in spec.contact_fibers)

    def reached_from(sources: set[int]) -> set[int]:
        reached = set(sources)
        queue = deque(sources)
        while queue:
            here = queue.popleft()
            for other in graph[here]:
                if other not in reached:
                    reached.add(other)
                    queue.append(other)
        return reached

    from_boundary = reached_from(boundary_fibers)
    from_contact = reached_from(contact_fibers)
    near_mask = np.linalg.norm(network.r0, axis=1) <= cfg.cell_radius + 15.0
    # A temporary clone is not needed: both metrics read the current initial state.
    report = {
        "source_fibers": int(cfg.n_fibers),
        "realised_fragments": int(len(network.fibers)),
        "beads": int(len(network.r)),
        "candidate_crosslinks": int(candidate_count),
        "crosslinks": int(len(network.crosslinks)),
        "crosslink_probability": probability,
        "outer_fixed_beads": int(network.fixed.sum()),
        "inner_boundary_beads": int(len(inner)),
        "contact_fragments": int(len(spec.contact_fibers)),
        "boundary_connected_fraction": len(from_boundary) / max(len(network.fibers), 1),
        "contact_connected_fraction": len(from_contact) / max(len(network.fibers), 1),
        "contact_fragments_reaching_boundary": int(len(contact_fibers & from_boundary)),
        "initial_global_nematic_order": _global_nematic_order(network),
        "initial_radial_order_by_shell": _radial_order(network).tolist(),
        "near_cell_beads": int(near_mask.sum()),
        "seed": int(spec.seed_used),
    }
    return network, spec, report


def graph_distance_from_contact(
    network: Network, direct_fibers: list[int] | None = None
) -> np.ndarray:
    """Shortest permanent-crosslink distance from a natural cell contact."""

    graph = [set() for _ in network.fibers]
    for link in network.crosslinks:
        fa = int(network.edge_fiber[link.edge_a])
        fb = int(network.edge_fiber[link.edge_b])
        graph[fa].add(fb)
        graph[fb].add(fa)
    distance = np.full(len(network.fibers), -1, dtype=int)
    queue = deque()
    sources = network.contact_fibers if direct_fibers is None else direct_fibers
    for fid in sources:
        distance[int(fid)] = 0
        queue.append(int(fid))
    while queue:
        here = queue.popleft()
        for other in graph[here]:
            if distance[other] < 0:
                distance[other] = distance[here] + 1
                queue.append(other)
    return distance


def surface_contact_patches(
    network: Network, center: np.ndarray | None = None
) -> list[ContactPatch]:
    """Hybrid contact selection followed by Gaussian distance weighting.

    One closest continuous material point is considered per physical fibre
    fragment.  Only points in the 0--``contact_width`` surface band are
    eligible.  The Gaussian therefore distributes force *within* a local
    contact region; it is never evaluated over the whole ECM.
    """

    center = np.zeros(2) if center is None else np.asarray(center, dtype=float)
    candidates: list[ContactPatch] = []
    for fid in range(len(network.fibers)):
        edge_ids = np.flatnonzero(network.edge_fiber == fid)
        pairs = network.edges[edge_ids]
        start = network.r[pairs[:, 0]]
        delta = network.r[pairs[:, 1]] - start
        alpha = np.clip(
            np.sum((center - start) * delta, axis=1)
            / np.maximum(np.sum(delta * delta, axis=1), 1e-12),
            0.0,
            1.0,
        )
        points = start + alpha[:, None] * delta
        radial = points - center
        radius = np.linalg.norm(radial, axis=1)
        surface = radius - network.cfg.cell_radius
        local = int(np.argmin(surface))
        gap = float(surface[local])
        if gap < -1e-7 or gap > network.cfg.contact_width:
            continue
        candidates.append(
            ContactPatch(
                fiber=fid,
                edge=int(edge_ids[local]),
                alpha=float(alpha[local]),
                point=points[local].copy(),
                surface_distance=max(0.0, gap),
                weight=0.0,
                normal_in=-radial[local] / max(float(radius[local]), 1e-12),
            )
        )
    if not candidates:
        return []
    raw = np.exp(
        -np.square([patch.surface_distance for patch in candidates])
        / network.cfg.gaussian_sigma**2
    )
    raw /= raw.sum()
    for patch, weight in zip(candidates, raw):
        patch.weight = float(weight)
    return candidates


def _graph_class_metrics(
    network: Network, distance: np.ndarray, initial_radial_by_fiber: np.ndarray
) -> dict[str, dict[str, float]]:
    """Separate direct and indirectly crosslink-coupled fibre responses."""

    bead_displacement = np.linalg.norm(network.r - network.r0, axis=1)
    fiber_displacement = np.asarray(
        [float(np.mean(bead_displacement[ids])) for ids in network.fibers]
    )
    current = radial_order_by_fiber(network)
    delta = current - initial_radial_by_fiber
    masks = {
        "direct": distance == 0,
        "one_hop": distance == 1,
        "two_plus": distance >= 2,
        "unconnected": distance < 0,
    }
    out = {}
    for label, mask in masks.items():
        out[label] = {
            "n": int(mask.sum()),
            "mean_displacement": float(np.mean(fiber_displacement[mask])) if np.any(mask) else 0.0,
            "mean_delta_radial_order": float(np.mean(delta[mask])) if np.any(mask) else 0.0,
        }
    return out


def radial_order_by_fiber(network: Network) -> np.ndarray:
    """Mean signed radial order of every physical fibre fragment."""

    values = np.empty(len(network.fibers))
    for fid, ids in enumerate(network.fibers):
        positions = network.r[ids]
        segment = np.diff(positions, axis=0)
        tangent = segment / np.maximum(np.linalg.norm(segment, axis=1), 1e-12)[:, None]
        midpoint = 0.5 * (positions[:-1] + positions[1:])
        radial = midpoint / np.maximum(np.linalg.norm(midpoint, axis=1), 1e-12)[:, None]
        values[fid] = float(np.mean(2.0 * np.square(np.sum(tangent * radial, axis=1)) - 1.0))
    return values


def elastic_energy_by_fiber(network: Network) -> tuple[np.ndarray, np.ndarray]:
    """Resolve stretching and bending energy onto physical fibre fragments."""

    i, j = network.edges.T
    delta = network.r[j] - network.r[i]
    length = np.linalg.norm(delta, axis=1)
    extension = length - network.l0
    stiffness = np.where(extension >= 0.0, network.k_tension, network.k_compression)
    edge_energy = 0.5 * stiffness * extension**2
    stretch = np.bincount(
        network.edge_fiber, edge_energy, minlength=len(network.fibers)
    )

    a, b, c = network.triplets.T
    curvature = network.r[a] - 2.0 * network.r[b] + network.r[c] - network.curvature0
    triplet_energy = 0.5 * network.bend_coefficient * np.sum(curvature**2, axis=1)
    # Every triplet lies within one fragment; the middle bead uniquely identifies it.
    bead_to_fiber = np.empty(len(network.r), dtype=int)
    for fid, ids in enumerate(network.fibers):
        bead_to_fiber[ids] = fid
    bend = np.bincount(
        bead_to_fiber[b], triplet_energy, minlength=len(network.fibers)
    )
    return stretch, bend


def _benchmark_frame(
    network: Network,
    time: float,
    contraction: float,
    inner: np.ndarray,
    inner_force: np.ndarray,
) -> dict:
    fs, fb, fx, energy, strain, link_force = network.elastic_forces(True)
    stretch_by_fiber, bend_by_fiber = elastic_energy_by_fiber(network)
    ratio = energy["stretch"] / max(energy["bend"], 1e-18)
    recruited = stretch_by_fiber > bend_by_fiber
    radial = network.r[inner] / np.maximum(
        np.linalg.norm(network.r[inner], axis=1), 1e-12
    )[:, None]
    # Positive means that the deformed network pulls the prescribed boundary
    # outward and the contracting cell would need an equal inward force.
    radial_reaction = float(np.sum(inner_force * radial))
    return {
        "time": float(time),
        "contraction": float(contraction),
        "positions": network.r.copy(),
        "energy": dict(energy),
        "stretch_to_bend_ratio": float(ratio),
        "recruited_fraction": float(np.mean(recruited)),
        "required_inward_boundary_force": radial_reaction,
        "radial_order_by_shell": _radial_order(network),
        "global_nematic_order": _global_nematic_order(network),
        "bond_strain": strain.copy(),
        "crosslink_force": link_force.copy(),
    }


def run_prescribed_contraction(
    cfg: G4V3Config = G4V3Config(),
    *,
    contraction: float | None = None,
    spec: NetworkSpec | None = None,
    crosslink_probability: float | None = None,
) -> dict:
    """Contract the inner void and relax the overdamped elastic network.

    The inner boundary condition is

    ``r_i(t) = [1 - epsilon(t)] r_i(0)``

    and every unconstrained bead follows

    ``zeta dr_i/dt = F_stretch + F_bend + F_crosslink``.

    No motor, clutch, Gaussian active force, Brownian motion, plasticity or cell
    migration is present.  The reaction at the prescribed inner nodes is
    recorded but does not enter a cell equation at this stage.
    """

    cfg.validate()
    target_contraction = cfg.benchmark_contraction if contraction is None else float(contraction)
    if not 0.0 <= target_contraction < 0.5:
        raise ValueError("benchmark contraction must be in [0, 0.5)")
    probability = (
        cfg.benchmark_crosslink_probability
        if crosslink_probability is None
        else float(crosslink_probability)
    )
    network, spec, report = build_random_network(
        cfg, spec=spec, crosslink_probability=probability
    )
    inner = _inner_boundary_beads(network)
    if len(inner) < 4:
        raise RuntimeError("random geometry produced fewer than four inner-boundary contacts")
    outer = network.fixed.copy()
    inner0 = network.r0[inner].copy()
    nsteps = int(round(cfg.benchmark_relax_time / cfg.dt))
    every = max(1, int(round(cfg.benchmark_sample_interval / cfg.dt)))
    frames: list[dict] = []
    cumulative_dissipation = 0.0
    inner_force = np.zeros((len(inner), 2))

    for step in range(nsteps + 1):
        time = step * cfg.dt
        ramp = min(1.0, time / cfg.benchmark_ramp_time)
        current_contraction = target_contraction * ramp
        network.r[inner] = (1.0 - current_contraction) * inner0
        network.r[outer] = network.r0[outer]
        fs, fb, fx, _, _, _ = network.elastic_forces(True)
        total = fs + fb + fx
        inner_force = total[inner].copy()
        constrained = outer.copy()
        constrained[inner] = True
        velocity = total / cfg.bead_drag
        velocity[constrained] = 0.0

        if step % every == 0 or step == nsteps:
            frames.append(
                _benchmark_frame(
                    network, time, current_contraction, inner, inner_force
                )
            )
            frames[-1]["cumulative_drag_dissipation"] = float(cumulative_dissipation)
        if step == nsteps:
            break
        cumulative_dissipation += float(
            cfg.dt * cfg.bead_drag * np.sum(velocity * velocity)
        )
        network.r[~constrained] += cfg.dt * velocity[~constrained]

        if not np.all(np.isfinite(network.r)):
            raise FloatingPointError("non-finite benchmark state; reduce dt")

    return {
        "config": asdict(cfg),
        "spec": spec,
        "network": network,
        "report": report,
        "inner_beads": inner,
        "frames": frames,
        "target_contraction": target_contraction,
        "crosslink_probability": probability,
        "initial_global_nematic_order": float(frames[0]["global_nematic_order"]),
        "initial_radial_order_by_shell": frames[0]["radial_order_by_shell"].copy(),
        "final_radial_order_by_shell": frames[-1]["radial_order_by_shell"].copy(),
        "final_stretch_to_bend_ratio": float(frames[-1]["stretch_to_bend_ratio"]),
        "final_recruited_fraction": float(frames[-1]["recruited_fraction"]),
    }


def run_benchmark_suite(
    cfg: G4V3Config = G4V3Config(),
    *,
    contractions: tuple[float, ...] = (0.05, 0.10, 0.20),
    bending_multipliers: tuple[float, ...] = (0.05, 0.25, 1.0),
) -> list[dict]:
    """One-factor grid on one geometry; no result is selected for appearance."""

    spec = make_random_void_spec(cfg)
    results = []
    for bending in bending_multipliers:
        case_cfg = replace(cfg, bending_multiplier=float(bending))
        for contraction in contractions:
            results.append(
                run_prescribed_contraction(
                    case_cfg, contraction=float(contraction), spec=spec
                )
            )
    return results


def run_fixed_cell_passive_check(
    cfg: G4V3Config = G4V3Config(), *, spec: NetworkSpec | None = None
) -> dict:
    """v3-B: verify that a fixed steric cell alone creates no remodeling."""

    network, spec, report = build_random_network(cfg, spec=spec)
    center = np.zeros(2)
    initial_shell = _radial_order(network)
    initial_repulsion = network.repulsion_forces(center)
    active = np.zeros_like(network.r)
    nsteps = max(1, int(round(1.0 / cfg.dt)))
    for _ in range(nsteps):
        fast_advance(network, active, center, cfg.dt)
    final_shell = _radial_order(network)
    return {
        "config": asdict(cfg),
        "spec": spec,
        "network": network,
        "report": report,
        "initial_max_repulsion": float(np.max(np.linalg.norm(initial_repulsion, axis=1))),
        "max_displacement": float(np.max(np.linalg.norm(network.r - network.r0, axis=1))),
        "delta_radial_order_by_shell": final_shell - initial_shell,
    }


def _traction_frame(
    network: Network,
    center: np.ndarray,
    patches: list[ContactPatch],
    vectors: np.ndarray,
    distance: np.ndarray,
    initial_shell: np.ndarray,
    initial_radial_by_fiber: np.ndarray,
    time: float,
    requested_force: float,
    cell_reaction: np.ndarray,
) -> dict:
    _, _, _, energy, strain, link_force = network.elastic_forces(True)
    stretch_by_fiber, bend_by_fiber = elastic_energy_by_fiber(network)
    shell = _radial_order(network)
    active_vector = -cell_reaction
    balance_error = float(np.linalg.norm(active_vector + cell_reaction))
    return {
        "time": float(time),
        "positions": network.r.copy(),
        "contact_points": np.asarray([patch_point(network, patch) for patch in patches]),
        "contact_vectors": np.asarray(vectors).copy(),
        "requested_force": float(requested_force),
        "distributed_force_magnitude": float(np.sum(np.linalg.norm(vectors, axis=1))),
        "active_force_vector": active_vector.copy(),
        "constrained_cell_reaction": np.asarray(cell_reaction).copy(),
        "action_reaction_error": balance_error,
        "radial_order_by_shell": shell,
        "delta_radial_order_by_shell": shell - initial_shell,
        "energy": dict(energy),
        "stretch_to_bend_ratio": float(
            energy["stretch"] / max(energy["bend"], 1e-18)
        ),
        "recruited_fraction": float(np.mean(stretch_by_fiber > bend_by_fiber)),
        "bond_strain": strain.copy(),
        "crosslink_force": link_force.copy(),
        "graph_class_metrics": _graph_class_metrics(
            network, distance, initial_radial_by_fiber
        ),
    }


def run_fixed_material_traction(
    cfg: G4V3Config = G4V3Config(),
    *,
    total_force: float | None = None,
    spec: NetworkSpec | None = None,
    crosslink_probability: float | None = None,
) -> dict:
    """v3-C: deterministic traction on fixed natural material points.

    Eligibility and weighting are intentionally separate:

    ``0 <= d_surface <= contact_width``

    ``w_a = exp(-d_a^2/sigma_c^2) / sum_b exp(-d_b^2/sigma_c^2)``

    ``F_a = w_a F_total n_in``.

    Non-contact fibres receive no active term.  They can move only through
    permanent crosslinks and collagen stretching/bending.  The equal and
    opposite cell reaction is recorded but the cell remains fixed.
    """

    cfg.validate()
    force = cfg.total_pull_force if total_force is None else float(total_force)
    probability = cfg.crosslink_probability if crosslink_probability is None else float(crosslink_probability)
    network, spec, report = build_random_network(
        cfg, spec=spec, crosslink_probability=probability
    )
    center = np.zeros(2)
    patches = surface_contact_patches(network, center)
    if not patches:
        raise RuntimeError("random geometry produced no contact-band material points")
    direct = [int(patch.fiber) for patch in patches]
    distance = graph_distance_from_contact(network, direct)
    initial_shell = _radial_order(network)
    initial_radial_by_fiber = radial_order_by_fiber(network)
    nsteps = int(round(cfg.traction_duration / cfg.dt))
    every = max(1, int(round(cfg.traction_sample_interval / cfg.dt)))
    frames: list[dict] = []
    vectors = np.zeros((len(patches), 2))

    for step in range(nsteps + 1):
        time = step * cfg.dt
        ramp = min(1.0, time / max(cfg.force_ramp_time, cfg.dt))
        requested = ramp * force
        active, _, vectors = active_from_patches(
            network, patches, requested, center
        )
        reaction = -active.sum(axis=0)
        if step % every == 0 or step == nsteps:
            frames.append(
                _traction_frame(
                    network,
                    center,
                    patches,
                    vectors,
                    distance,
                    initial_shell,
                    initial_radial_by_fiber,
                    time,
                    requested,
                    reaction,
                )
            )
        if step == nsteps:
            break
        fast_advance(network, active, center, cfg.dt)

    return {
        "config": asdict(cfg),
        "spec": spec,
        "network": network,
        "report": report,
        "frames": frames,
        "patches": patches,
        "direct_fibers": direct,
        "graph_distance": distance,
        "total_force": force,
        "crosslink_probability": probability,
        "initial_radial_order_by_shell": initial_shell,
        "final_delta_radial_order_by_shell": frames[-1]["delta_radial_order_by_shell"].copy(),
        "final_graph_class_metrics": frames[-1]["graph_class_metrics"],
        "max_displacement": float(np.max(np.linalg.norm(network.r - network.r0, axis=1))),
    }
