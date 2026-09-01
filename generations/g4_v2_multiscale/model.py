"""G4 v2: long-horizon elastic calibration and shared-load clutch failure.

This module is deliberately separate from :mod:`g4_interactive_calibration`.
The latter remains the frozen, short G4 v1 experiment.  G4 v2 reuses its
99-fibre geometry, hybrid surface contact and permanent probabilistic links,
then adds:

* two-hour elastic experiments with sparse observation frames;
* direct/one-hop/two+-hop orientation and displacement statistics;
* an independent-clutch control and an equal-load-sharing cluster model;
* high-frequency event capture around a complete site failure; and
* fixed/released rigid-cell experiments with no prescribed polarity.

Units are micrometre (um), nanonewton (nN), and second (s).  The ECM equation
is unchanged from G4 v1::

    zeta_b dr_i/dt = F_stretch + F_bend + F_crosslink
                     + F_repulsion + F_active.

The shared-load model follows the cluster state equation

    r_i = i k_off^0 exp(F_site / (i F_b)),
    g_i = (N-i) k_on,

where ``i`` is the number of closed effective clutches at one material-point
site.  A site is visibly detached only when ``i == 0``.
"""

from __future__ import annotations

from collections import deque
from dataclasses import asdict, dataclass, replace
from pathlib import Path
import math
import sys
from typing import Literal

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from generations.g4_interactive_calibration.model import (  # noqa: E402
    G4Config,
    _assemble_g4_spec,
    active_from_patches,
    bell_off_rate,
    build_g4_network,
    direct_contact_patches,
    graph_distance_from_direct,
    patch_point,
    radial_order_by_fiber,
)

try:  # Numba keeps the mandated dt=0.05 s feasible over two hours.
    from numba import njit
except Exception:  # pragma: no cover - requirements install Numba in production.
    njit = None


BoundaryMode = Literal["anchored", "mobile"]
ClutchMode = Literal["independent", "shared"]


@dataclass(frozen=True)
class G4V2Config(G4Config):
    """Configuration shared by G4 v2 A--D.

    Defaults preserve the G4 v1 effective-clutch and ECM values.  Long stage
    durations are set by the builder so tests can still use short runs.
    """

    boundary_mode: BoundaryMode = "anchored"
    metric_sample_interval: float = 4.0
    event_sample_interval: float = 0.10
    event_half_window: float = 10.0
    counter_seed: int = 3042

    def validate(self) -> None:
        super().validate()
        if self.boundary_mode not in ("anchored", "mobile"):
            raise ValueError("boundary_mode must be anchored or mobile")
        if self.metric_sample_interval < self.dt:
            raise ValueError("metric_sample_interval must be at least dt")
        if self.event_sample_interval < self.dt:
            raise ValueError("event_sample_interval must be at least dt")
        if self.event_half_window <= 0.0:
            raise ValueError("event_half_window must be positive")


def build_v2_network(cfg: G4V2Config, *, spec=None):
    """Reuse G4 v1 topology and optionally release only the outer anchors."""

    network, spec, report = build_g4_network(cfg, spec=spec)
    if cfg.boundary_mode == "mobile":
        network.fixed[:] = False
    report = dict(report)
    report["boundary_mode"] = cfg.boundary_mode
    report["fixed_beads"] = int(network.fixed.sum())
    return network, spec, report


if njit is not None:

    @njit(cache=True)
    def _active_numba(r, edges, patch_edges, patch_alpha, magnitudes, center,
                      active, points, vectors):
        active[:, :] = 0.0
        for p in range(patch_edges.shape[0]):
            edge = patch_edges[p]
            i = edges[edge, 0]
            j = edges[edge, 1]
            alpha = patch_alpha[p]
            px = (1.0 - alpha) * r[i, 0] + alpha * r[j, 0]
            py = (1.0 - alpha) * r[i, 1] + alpha * r[j, 1]
            dx = center[0] - px
            dy = center[1] - py
            length = math.sqrt(dx * dx + dy * dy)
            if length < 1e-12:
                length = 1e-12
            fx = magnitudes[p] * dx / length
            fy = magnitudes[p] * dy / length
            active[i, 0] += (1.0 - alpha) * fx
            active[i, 1] += (1.0 - alpha) * fy
            active[j, 0] += alpha * fx
            active[j, 1] += alpha * fy
            points[p, 0] = px
            points[p, 1] = py
            vectors[p, 0] = fx
            vectors[p, 1] = fy

    @njit(cache=True)
    def _advance_numba(
        r, r0, fixed, edges, l0, k_tension, k_compression,
        triplets, curvature0, bend_coefficient,
        link_edge_a, link_edge_b, link_alpha_a, link_alpha_b,
        link_rest, link_stiffness, active, center,
        cell_radius, cell_clearance, repulsion_stiffness,
        bead_drag, dt, velocity, total,
    ):
        """Allocation-free implementation of the exact G4 force balance."""

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
            stiffness = k_tension[e] if extension >= 0.0 else k_compression[e]
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
            dx = r[i, 0] - center[0]
            dy = r[i, 1] - center[1]
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


    @njit(cache=True)
    def _integrate_elastic_numba(
        r, r0, fixed, edges, l0, k_tension, k_compression,
        triplets, curvature0, bend_coefficient,
        link_edge_a, link_edge_b, link_alpha_a, link_alpha_b,
        link_rest, link_stiffness, patch_edges, patch_alpha, patch_weights,
        center, total_pull_force, force_ramp_time,
        cell_radius, cell_clearance, repulsion_stiffness,
        bead_drag, dt, nsteps, output_steps,
    ):
        """Run the fixed-cell elastic loop in compiled code and save sparse states."""

        out = np.empty((output_steps.shape[0], r.shape[0], 2), dtype=np.float64)
        active = np.zeros_like(r)
        points = np.zeros((patch_edges.shape[0], 2))
        vectors = np.zeros((patch_edges.shape[0], 2))
        velocity = np.zeros_like(r)
        total = np.zeros_like(r)
        output_index = 0
        for step in range(nsteps + 1):
            time = step * dt
            ramp = min(1.0, time / max(force_ramp_time, dt))
            magnitudes = patch_weights * (ramp * total_pull_force)
            _active_numba(
                r, edges, patch_edges, patch_alpha, magnitudes, center,
                active, points, vectors,
            )
            if output_index < output_steps.shape[0] and step == output_steps[output_index]:
                out[output_index, :, :] = r
                output_index += 1
            if step < nsteps:
                _advance_numba(
                    r, r0, fixed, edges, l0, k_tension, k_compression,
                    triplets, curvature0, bend_coefficient,
                    link_edge_a, link_edge_b, link_alpha_a, link_alpha_b,
                    link_rest, link_stiffness, active, center,
                    cell_radius, cell_clearance, repulsion_stiffness,
                    bead_drag, dt, velocity, total,
                )
        return out


class FastStepper:
    """Preallocated, compiled overdamped integrator for long G4 v2 runs."""

    def __init__(self, network):
        self.network = network
        self.velocity = np.zeros_like(network.r)
        self.total = np.zeros_like(network.r)

    def step(self, active: np.ndarray, center: np.ndarray, dt: float):
        if njit is None:  # pragma: no cover
            from generations.g4_interactive_calibration.model import fast_advance
            velocity, total, _ = fast_advance(self.network, active, center, dt)
            self.velocity[:] = velocity
            self.total[:] = total
            return self.velocity, self.total
        n = self.network
        _advance_numba(
            n.r, n.r0, n.fixed, n.edges, n.l0, n.k_tension, n.k_compression,
            n.triplets, n.curvature0, n.bend_coefficient,
            n.link_edge_a, n.link_edge_b, n.link_alpha_a, n.link_alpha_b,
            n.link_rest, n.link_stiffness, active, np.asarray(center),
            n.cfg.cell_radius, n.cfg.cell_clearance, n.cfg.repulsion_stiffness,
            n.cfg.bead_drag, dt, self.velocity, self.total,
        )
        if not np.all(np.isfinite(n.r)):
            raise FloatingPointError("non-finite bead position; reduce dt")
        return self.velocity, self.total


class ActiveMapper:
    """Project material-point traction to beads without per-step allocation."""

    def __init__(self, network, patches):
        self.network = network
        self.active = np.zeros_like(network.r)
        self.points = np.zeros((len(patches), 2))
        self.vectors = np.zeros((len(patches), 2))
        self.update_patches(patches)

    def update_patches(self, patches):
        self.patches = patches
        self.patch_edges = np.asarray([p.edge for p in patches], dtype=np.int64)
        self.patch_alpha = np.asarray([p.alpha for p in patches], dtype=float)
        self.weights = np.asarray([p.weight for p in patches], dtype=float)

    def map(self, scalar_forces, center, *, weighted=False):
        if np.isscalar(scalar_forces):
            magnitude = float(scalar_forces)
            magnitudes = magnitude * self.weights if weighted else np.full(len(self.patches), magnitude)
        else:
            magnitudes = np.asarray(scalar_forces, dtype=float)
        if njit is None:  # pragma: no cover
            active, points, vectors = active_from_patches(
                self.network, self.patches,
                float(scalar_forces) if weighted and np.isscalar(scalar_forces) else magnitudes,
                center,
            )
            self.active[:] = active
            self.points[:] = points
            self.vectors[:] = vectors
        else:
            _active_numba(
                self.network.r, self.network.edges, self.patch_edges,
                self.patch_alpha, magnitudes, np.asarray(center),
                self.active, self.points, self.vectors,
            )
        return self.active, self.points, self.vectors


def _fiber_angles(network) -> np.ndarray:
    """Nematic fibre angle from the mean segment orientation."""

    out = np.zeros(len(network.fibers))
    for fid, ids in enumerate(network.fibers):
        seg = np.diff(network.r[ids], axis=0)
        angle = np.arctan2(seg[:, 1], seg[:, 0])
        out[fid] = 0.5 * math.atan2(float(np.mean(np.sin(2 * angle))),
                                    float(np.mean(np.cos(2 * angle))))
    return out


def _nematic_delta(current: np.ndarray, initial: np.ndarray) -> np.ndarray:
    return 0.5 * np.arctan2(np.sin(2.0 * (current - initial)),
                            np.cos(2.0 * (current - initial)))


def graph_class_statistics(network, distance: np.ndarray, initial_angles: np.ndarray) -> dict:
    """Displacement and orientation change by actual crosslink-graph distance."""

    bead_disp = np.linalg.norm(network.r - network.r0, axis=1)
    fiber_disp = np.asarray([float(np.mean(bead_disp[ids])) for ids in network.fibers])
    delta = np.abs(_nematic_delta(_fiber_angles(network), initial_angles))
    masks = {
        "direct": distance == 0,
        "one_hop": distance == 1,
        "two_plus": distance >= 2,
        "unconnected": distance < 0,
    }
    stats: dict[str, dict[str, float]] = {}
    for name, mask in masks.items():
        if np.any(mask):
            d = fiber_disp[mask]
            a = delta[mask]
            stats[name] = {
                "n": int(mask.sum()),
                "mean_displacement": float(np.mean(d)),
                "p90_displacement": float(np.percentile(d, 90)),
                "max_displacement": float(np.max(d)),
                "mean_abs_delta_theta": float(np.mean(a)),
                "p90_abs_delta_theta": float(np.percentile(a, 90)),
            }
        else:
            stats[name] = {
                "n": 0, "mean_displacement": 0.0, "p90_displacement": 0.0,
                "max_displacement": 0.0, "mean_abs_delta_theta": 0.0,
                "p90_abs_delta_theta": 0.0,
            }
    return stats


def _contact_arrays(network, patches, vectors):
    points = np.asarray([patch_point(network, p) for p in patches], dtype=float)
    return points, np.asarray(vectors, dtype=float)


def _frame(
    network, center, theta, time, patches, vectors, distance, initial_angles,
    *, traction=0.0, bound=None, site_force=None, cumulative_slips=0,
    site_failures=0, breaks=None, binds=None,
):
    points, vectors = _contact_arrays(network, patches, vectors)
    angles = _fiber_angles(network)
    return {
        "time": float(time),
        "positions": network.r.copy(),
        "center": np.asarray(center, dtype=float).copy(),
        "theta": float(theta),
        "contact_points": points,
        "contact_vectors": vectors,
        "traction": float(traction),
        "bound": None if bound is None else np.asarray(bound, dtype=np.uint8).copy(),
        "site_force": None if site_force is None else np.asarray(site_force, dtype=float).copy(),
        "cumulative_slips": int(cumulative_slips),
        "site_failures": int(site_failures),
        "breaks": None if breaks is None else np.asarray(breaks, dtype=np.uint8).copy(),
        "binds": None if binds is None else np.asarray(binds, dtype=np.uint8).copy(),
        "radial_order": float(np.mean(radial_order_by_fiber(network, center))),
        "fiber_delta_theta": _nematic_delta(angles, initial_angles),
        "graph_stats": graph_class_statistics(network, distance, initial_angles),
    }


def run_elastic(cfg: G4V2Config, *, spec=None) -> dict:
    """G4 v2 A/B: fixed-cell long-horizon elastic experiment."""

    cfg.validate()
    network, spec, report = build_v2_network(cfg, spec=spec)
    center = np.zeros(2)
    patches = direct_contact_patches(network, center)
    direct = [p.fiber for p in patches]
    distance = graph_distance_from_direct(network, direct)
    initial_angles = _fiber_angles(network)
    mapper = ActiveMapper(network, patches)
    nsteps = int(round(cfg.duration / cfg.dt))
    geometry_every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    metric_every = max(1, int(round(cfg.metric_sample_interval / cfg.dt)))
    frames: list[dict] = []
    metrics: list[dict] = []
    geometry_steps = set(range(0, nsteps + 1, geometry_every)) | {nsteps}
    metric_steps = set(range(0, nsteps + 1, metric_every)) | {nsteps}
    output_steps = np.asarray(sorted(geometry_steps | metric_steps), dtype=np.int64)
    if njit is not None:
        n = network
        positions = _integrate_elastic_numba(
            n.r, n.r0, n.fixed, n.edges, n.l0, n.k_tension, n.k_compression,
            n.triplets, n.curvature0, n.bend_coefficient,
            n.link_edge_a, n.link_edge_b, n.link_alpha_a, n.link_alpha_b,
            n.link_rest, n.link_stiffness, mapper.patch_edges, mapper.patch_alpha,
            mapper.weights, center, cfg.total_pull_force, cfg.force_ramp_time,
            cfg.cell_radius, cfg.cell_clearance, cfg.repulsion_stiffness,
            cfg.bead_drag, cfg.dt, nsteps, output_steps,
        )
    else:  # pragma: no cover
        positions = np.empty((len(output_steps), len(network.r), 2))
        stepper = FastStepper(network)
        output_lookup = {int(s): i for i, s in enumerate(output_steps)}
        for step in range(nsteps + 1):
            time = step * cfg.dt
            ramp = min(1.0, time / max(cfg.force_ramp_time, cfg.dt))
            active, _, _ = mapper.map(ramp * cfg.total_pull_force, center, weighted=True)
            if step in output_lookup:
                positions[output_lookup[step]] = network.r
            if step < nsteps:
                stepper.step(active, center, cfg.dt)
    final_positions = network.r.copy()
    for output_index, step in enumerate(output_steps):
        network.r[:] = positions[output_index]
        time = float(step * cfg.dt)
        ramp = min(1.0, time / max(cfg.force_ramp_time, cfg.dt))
        _, _, vectors = mapper.map(ramp * cfg.total_pull_force, center, weighted=True)
        record = _frame(
            network, center, 0.0, time, patches, vectors, distance, initial_angles,
            traction=ramp * cfg.total_pull_force,
        )
        if int(step) in metric_steps:
            metrics.append(record)
        if int(step) in geometry_steps:
            frames.append(record if int(step) not in metric_steps else dict(record))
    network.r[:] = final_positions
    return {
        "config": asdict(cfg), "spec": spec, "network": network, "report": report,
        "frames": frames, "metrics": metrics, "direct_fibers": direct,
        "graph_distance": distance, "initial_angles": initial_angles,
        "final_graph_stats": graph_class_statistics(network, distance, initial_angles),
        "max_displacement": float(np.max(np.linalg.norm(network.r - network.r0, axis=1))),
    }


@dataclass
class ClutchState:
    bound: np.ndarray
    extension: np.ndarray
    site_extension: np.ndarray
    cumulative_slips: int = 0
    cumulative_site_failures: int = 0


_MASK64 = (1 << 64) - 1
_STEP_CONST = 0x9E3779B97F4A7C15
_CHANNEL_CONST = 0xD1B54A32D192ED03


def _counter_uniform_vector(seed: int, step: int, channel: int, count: int) -> np.ndarray:
    ids = np.arange(count, dtype=np.uint64)
    base = (int(seed) + step * _STEP_CONST + channel * _CHANNEL_CONST) & _MASK64
    with np.errstate(over="ignore"):
        z = ids + np.uint64(base)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        z = z ^ (z >> np.uint64(31))
    return ((z >> np.uint64(11)).astype(np.float64)) * (1.0 / 9007199254740992.0)


def counter_uniforms(cfg: G4V2Config, step: int, channel: int) -> np.ndarray:
    """Counter-addressed uniforms shared across all mechanical comparisons."""

    out = _counter_uniform_vector(
        cfg.counter_seed, step, channel,
        cfg.n_contact_sectors * cfg.n_clutches_per_site,
    )
    return out.reshape(cfg.n_contact_sectors, cfg.n_clutches_per_site)


def _independent_step(cfg, state, substrate_speed, u_on, u_off):
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
    breaks = state.bound & (u_off < 1.0 - np.exp(-bell_off_rate(force, cfg) * cfg.dt))
    before_count = state.bound.sum(axis=1)
    state.bound[breaks] = False
    state.extension[breaks] = 0.0
    on_probability = 1.0 - math.exp(-cfg.clutch_on_rate * cfg.dt)
    binds = (~state.bound) & (u_on < on_probability)
    state.bound[binds] = True
    state.extension[binds] = 0.0
    after_break_count = before_count - breaks.sum(axis=1)
    site_failures = (before_count > 0) & (after_break_count == 0)
    state.cumulative_slips += int(breaks.sum())
    state.cumulative_site_failures += int(site_failures.sum())
    force = cfg.clutch_stiffness * state.extension * state.bound
    return force, force.sum(axis=1), breaks, binds, site_failures


def shared_load_hazard(site_force: np.ndarray, bound_count: np.ndarray, cfg: G4V2Config):
    """Per-clutch and whole-cluster hazards under equal load sharing."""

    count = np.maximum(np.asarray(bound_count, dtype=float), 1.0)
    per_clutch = np.asarray(site_force, dtype=float) / count
    per_rate = bell_off_rate(per_clutch, cfg)
    total_rate = np.asarray(bound_count, dtype=float) * per_rate
    return per_clutch, per_rate, total_rate


def _shared_step(cfg, state, substrate_speed, u_on, u_off):
    before_count = state.bound.sum(axis=1)
    site_stiffness = cfg.n_clutches_per_site * cfg.clutch_stiffness
    traction_before = np.where(before_count > 0, site_stiffness * state.site_extension, 0.0)
    actin_speed = cfg.unloaded_actin_speed * np.maximum(
        0.0, 1.0 - traction_before / cfg.motor_stall_per_site
    )
    relative = np.maximum(0.0, actin_speed - substrate_speed)
    state.site_extension[before_count > 0] += cfg.dt * relative[before_count > 0]
    site_force = np.where(before_count > 0, site_stiffness * state.site_extension, 0.0)
    per_force, per_rate, _ = shared_load_hazard(site_force, before_count, cfg)
    breaks = state.bound & (u_off < 1.0 - np.exp(-per_rate[:, None] * cfg.dt))
    state.bound[breaks] = False
    after_break_count = state.bound.sum(axis=1)
    site_failures = (before_count > 0) & (after_break_count == 0)
    state.site_extension[site_failures] = 0.0
    on_probability = 1.0 - math.exp(-cfg.clutch_on_rate * cfg.dt)
    binds = (~state.bound) & (u_on < on_probability)
    state.bound[binds] = True
    after_count = state.bound.sum(axis=1)
    state.site_extension[after_count == 0] = 0.0
    site_force = np.where(after_count > 0, site_stiffness * state.site_extension, 0.0)
    per_force = np.divide(site_force, np.maximum(after_count, 1), dtype=float)
    force = state.bound * per_force[:, None]
    state.cumulative_slips += int(breaks.sum())
    state.cumulative_site_failures += int(site_failures.sum())
    return force, site_force, breaks, binds, site_failures


def _patch_substrate_speeds(network, velocity, patches, center):
    out = np.zeros(len(patches))
    for index, patch in enumerate(patches):
        point = patch_point(network, patch)
        inward = np.asarray(center) - point
        inward /= max(float(np.linalg.norm(inward)), 1e-12)
        i, j = network.edges[patch.edge]
        material_velocity = (1.0 - patch.alpha) * velocity[i] + patch.alpha * velocity[j]
        out[index] = float(material_velocity @ inward)
    return out


def _local_subset(network, patches=None, radius: float = 14.0):
    """Return only the near-cell beads needed by the event microscope.

    A 14 µm radius includes the 10 µm cell, the full 3 µm contact band and the
    endpoints of every current material-point contact edge.  It intentionally
    does not copy whole 20–80 µm contact fibres into a 20 s high-rate clip.
    """
    mask = np.linalg.norm(network.r0, axis=1) <= radius
    if patches is not None:
        for patch in patches:
            mask[network.edges[patch.edge]] = True
    bead_ids = np.flatnonzero(mask)
    inverse = np.full(len(network.r), -1, dtype=int)
    inverse[bead_ids] = np.arange(len(bead_ids))
    edge_mask = mask[network.edges[:, 0]] & mask[network.edges[:, 1]]
    edge_ids = np.flatnonzero(edge_mask)
    local_edges = inverse[network.edges[edge_ids]]
    return bead_ids, edge_ids, local_edges


def _detail_frame(network, local_beads, center, theta, time, patches, state,
                  site_force, breaks, binds, substrate_speed):
    points = np.asarray([patch_point(network, p) for p in patches])
    return {
        "time": float(time), "positions": network.r[local_beads].copy(),
        "center": np.asarray(center).copy(), "theta": float(theta),
        "contact_points": points, "bound": state.bound.astype(np.uint8).copy(),
        "site_force": np.asarray(site_force).copy(),
        "breaks": np.asarray(breaks, dtype=np.uint8).copy(),
        "binds": np.asarray(binds, dtype=np.uint8).copy(),
        "substrate_speed": np.asarray(substrate_speed).copy(),
    }


def run_clutch(
    cfg: G4V2Config, *, spec=None, mode: ClutchMode, moving: bool,
    capture_event: bool = False, event_center_time: float | None = None,
) -> dict:
    """Run independent or shared-load motor clutches on the same ECM."""

    cfg.validate()
    if mode not in ("independent", "shared"):
        raise ValueError("mode must be independent or shared")
    network, spec, report = build_v2_network(cfg, spec=spec)
    center = np.zeros(2)
    theta = 0.0
    patches = direct_contact_patches(network, center)
    if len(patches) != cfg.n_contact_sectors:
        raise RuntimeError(f"expected {cfg.n_contact_sectors} direct patches; found {len(patches)}")
    direct = [p.fiber for p in patches]
    distance = graph_distance_from_direct(network, direct)
    initial_angles = _fiber_angles(network)
    state = ClutchState(
        bound=np.zeros((cfg.n_contact_sectors, cfg.n_clutches_per_site), dtype=bool),
        extension=np.zeros((cfg.n_contact_sectors, cfg.n_clutches_per_site)),
        site_extension=np.zeros(cfg.n_contact_sectors),
    )
    stepper = FastStepper(network)
    mapper = ActiveMapper(network, patches)
    bead_velocity = np.zeros_like(network.r)
    substrate_speed = np.zeros(cfg.n_contact_sectors)
    site_force = np.zeros(cfg.n_contact_sectors)
    vectors = np.zeros((cfg.n_contact_sectors, 2))
    breaks = np.zeros_like(state.bound)
    binds = np.zeros_like(state.bound)
    nsteps = int(round(cfg.duration / cfg.dt))
    geometry_every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    metric_every = max(1, int(round(cfg.metric_sample_interval / cfg.dt)))
    update_every = max(1, int(round(cfg.contact_update_interval / cfg.dt)))
    event_every = max(1, int(round(cfg.event_sample_interval / cfg.dt)))
    pre_event_frames = max(1, int(round(cfg.event_half_window / cfg.event_sample_interval)))
    frames: list[dict] = []
    metrics: list[dict] = []
    event_frames: list[dict] = []
    ring: deque = deque(maxlen=pre_event_frames + 1)
    event_time = event_center_time
    event_site = -1
    event_finished = False
    local_beads, local_edge_ids, local_edges = _local_subset(network, patches)
    path_length = 0.0
    site_failure_times: list[float] = []

    for step in range(nsteps + 1):
        time = step * cfg.dt
        if step % metric_every == 0 or step == nsteps:
            metrics.append(_frame(
                network, center, theta, time, patches, vectors, distance, initial_angles,
                traction=float(site_force.sum()), bound=state.bound,
                site_force=site_force, cumulative_slips=state.cumulative_slips,
                site_failures=state.cumulative_site_failures, breaks=breaks, binds=binds,
            ))
        if step % geometry_every == 0 or step == nsteps:
            frames.append(_frame(
                network, center, theta, time, patches, vectors, distance, initial_angles,
                traction=float(site_force.sum()), bound=state.bound,
                site_force=site_force, cumulative_slips=state.cumulative_slips,
                site_failures=state.cumulative_site_failures, breaks=breaks, binds=binds,
            ))
        if capture_event and not event_finished and step % event_every == 0:
            detail = _detail_frame(
                network, local_beads, center, theta, time, patches, state,
                site_force, breaks, binds, substrate_speed,
            )
            if event_time is None:
                ring.append(detail)
            elif time >= event_time - cfg.event_half_window - 1e-9:
                event_frames.append(detail)
                if time >= event_time + cfg.event_half_window - 1e-9:
                    event_finished = True
        if step == nsteps:
            break

        u_on = counter_uniforms(cfg, step, 0)
        u_off = counter_uniforms(cfg, step, 1)
        if mode == "shared":
            _, site_force, breaks, binds, site_failures = _shared_step(
                cfg, state, substrate_speed, u_on, u_off
            )
        else:
            _, site_force, breaks, binds, site_failures = _independent_step(
                cfg, state, substrate_speed, u_on, u_off
            )
        if np.any(site_failures):
            for site in np.flatnonzero(site_failures):
                site_failure_times.append(time)
                if capture_event and event_time is None:
                    event_time = time
                    event_site = int(site)
                    event_frames = list(ring)

        active, _, vectors = mapper.map(site_force, center)
        reaction = -active.sum(axis=0)
        if moving:
            cell_velocity = reaction / cfg.cell_drag
            speed = float(np.linalg.norm(cell_velocity))
            if speed > cfg.max_cell_speed:
                cell_velocity *= cfg.max_cell_speed / speed
            angles = np.linspace(0.0, 2.0 * math.pi, cfg.n_contact_sectors, endpoint=False) + theta
            anchors = center + cfg.cell_radius * np.column_stack([np.cos(angles), np.sin(angles)])
            arms = anchors - center
            site_reactions = -vectors
            torque = float(np.sum(arms[:, 0] * site_reactions[:, 1]
                                  - arms[:, 1] * site_reactions[:, 0]))
            rotational_drag = cfg.rotational_drag_factor * cfg.cell_drag * cfg.cell_radius**2
            angular_velocity = float(np.clip(
                torque / rotational_drag, -cfg.max_angular_speed, cfg.max_angular_speed
            ))
            delta = cfg.dt * cell_velocity
            center += delta
            theta += cfg.dt * angular_velocity
            path_length += float(np.linalg.norm(delta))

        bead_velocity, _ = stepper.step(active, center, cfg.dt)
        substrate_speed = _patch_substrate_speeds(network, bead_velocity, patches, center)

        if step and step % update_every == 0:
            fresh = direct_contact_patches(network, center)
            fresh_by_fiber = {p.fiber: p for p in fresh}
            for site, patch in enumerate(patches):
                if not state.bound[site].any() and patch.fiber in fresh_by_fiber:
                    patches[site] = fresh_by_fiber[patch.fiber]
            mapper.update_patches(patches)

    return {
        "config": asdict(cfg), "spec": spec, "network": network, "report": report,
        "mode": mode, "moving": moving, "frames": frames, "metrics": metrics,
        "event_frames": event_frames, "event_time": event_time, "event_site": event_site,
        "local_beads": local_beads, "local_edge_ids": local_edge_ids,
        "local_edges": local_edges, "direct_fibers": direct,
        "graph_distance": distance, "initial_angles": initial_angles,
        "site_failure_times": np.asarray(site_failure_times),
        "path_length": float(path_length), "net_displacement": float(np.linalg.norm(center)),
        "final_graph_stats": graph_class_statistics(network, distance, initial_angles),
        "max_displacement": float(np.max(np.linalg.norm(network.r - network.r0, axis=1))),
    }


def shared_cluster_ensemble(
    cfg: G4V2Config, *, trials: int = 200, duration: float = 1800.0,
    seed0: int = 9000,
) -> dict:
    """Fast isolated-site ensemble used to select a non-extreme demo seed.

    The protocol uses zero substrate speed and the same motor force--velocity,
    Bell, rebinding and shared-load equations as the full ECM.  It is explicitly
    a clutch-cluster calibration, not 200 full collagen-network trajectories.
    """

    nsteps = int(round(duration / cfg.dt))
    lifetimes: list[float] = []
    first_failures = np.full(trials, np.nan)
    for trial in range(trials):
        trial_seed = seed0 + trial
        bound = np.zeros(cfg.n_clutches_per_site, dtype=bool)
        extension = 0.0
        episode_start: float | None = None
        for step in range(nsteps):
            count = int(bound.sum())
            site_stiffness = cfg.n_clutches_per_site * cfg.clutch_stiffness
            force = site_stiffness * extension if count else 0.0
            actin = cfg.unloaded_actin_speed * max(0.0, 1.0 - force / cfg.motor_stall_per_site)
            if count:
                extension += cfg.dt * actin
                force = site_stiffness * extension
                per_rate = float(bell_off_rate(force / count, cfg))
                breaks = bound & (_counter_uniform_vector(
                    trial_seed, step, 1, cfg.n_clutches_per_site)
                                  < 1.0 - math.exp(-per_rate * cfg.dt))
                bound[breaks] = False
            before_bind = int(bound.sum())
            binds = (~bound) & (_counter_uniform_vector(
                trial_seed, step, 0, cfg.n_clutches_per_site)
                                < 1.0 - math.exp(-cfg.clutch_on_rate * cfg.dt))
            bound[binds] = True
            after = int(bound.sum())
            time = (step + 1) * cfg.dt
            if episode_start is None and after > 0:
                episode_start = time
            if count > 0 and before_bind == 0:
                if episode_start is not None:
                    lifetimes.append(time - episode_start)
                if not np.isfinite(first_failures[trial]):
                    first_failures[trial] = time
                episode_start = time if after > 0 else None
                extension = 0.0
            elif after == 0:
                extension = 0.0
    finite = first_failures[np.isfinite(first_failures)]
    median_first = float(np.median(finite)) if len(finite) else math.nan
    if len(finite):
        valid_ids = np.flatnonzero(np.isfinite(first_failures))
        chosen_trial = int(valid_ids[np.argmin(np.abs(first_failures[valid_ids] - median_first))])
    else:
        chosen_trial = 0
    return {
        "trials": int(trials), "duration": float(duration),
        "lifetimes": np.asarray(lifetimes), "first_failures": first_failures,
        "median_lifetime": float(np.median(lifetimes)) if lifetimes else math.nan,
        "median_first_failure": median_first,
        "representative_seed": int(seed0 + chosen_trial),
        "representative_trial": chosen_trial,
    }


def with_counter_seed(cfg: G4V2Config, seed: int) -> G4V2Config:
    return replace(cfg, counter_seed=int(seed))
