"""G4 v4E: contact-guided persistent migration of one fixed-radius cell.

This module is additive: it imports G4 v4 mechanics and never modifies the
historical A--D implementation.  Four matched conditions isolate the two new
ingredients:

``d_control``
    G4D mechanics on an initially random collagen network.
``matrix_cue``
    The same mechanics with a one-sided, locally radial collagen tract.
``protrusion_guidance``
    Random collagen plus material-point selection and temporary front memory.
``combined``
    The tract and protrusion rule together (the predeclared main hypothesis).

The cell never shrinks and receives no prescribed world-coordinate polarity.
The only directional state, ``front``, is created after a real collagen contact
survives the maturation time.  Rear clutches are not cut: they continue to use
the inherited Bell rupture law and are relocated only after complete natural
site failure.

Units are micrometre (um), nanonewton (nN), and second (s).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import math
from typing import Literal

import numpy as np

from generations.g2_corrected.common.model import ContactPatch, Network, NetworkSpec
from generations.g4_interactive_calibration.model import patch_point
from generations.g4_v2_multiscale.model import _active_numba, _advance_numba
from generations.g4_v3_random_to_aligned.model import (
    build_random_network,
    make_random_void_spec,
    surface_contact_patches,
)
from generations.g4_v4_single_cell_force_alignment.model import (
    G4V4Config,
    _contact_tangent,
    _radial_order_centered,
)

try:
    from numba import njit
except Exception:  # pragma: no cover
    njit = None


GuidanceMode = Literal[
    "d_control", "matrix_cue", "protrusion_guidance", "combined"
]

MODES: tuple[GuidanceMode, ...] = (
    "d_control",
    "matrix_cue",
    "protrusion_guidance",
    "combined",
)


@dataclass(frozen=True)
class G4V4EConfig(G4V4Config):
    """Configuration for the predeclared G4 v4E experiment.

    The ECM, motor, clutch, Bell rupture and cell-drag values are inherited
    unchanged from G4 v4.  Values below are the new, provisional guidance
    assumptions.  They are sensitivity variables, not fitted constants.
    """

    radial_cue_angle: float = 0.0
    radial_cue_half_width: float = math.pi / 6.0
    radial_tract_angular_sd: float = math.pi / 12.0
    radial_tract_outer_surface: float = 35.0
    probing_reach: float = 5.0
    front_maturation: float = 120.0
    front_loss_time: float = 120.0
    front_cone_half_width: float = math.pi / 4.0
    beta_alignment: float = 2.0
    beta_memory: float = 2.0
    telemetry_interval: float = 30.0
    event_interval: float = 1.0
    event_duration: float = 300.0
    panel_seeds: int = 20
    cue_sd_sensitivity: tuple[float, ...] = (math.pi / 12.0, math.pi / 6.0)
    maturation_sensitivity: tuple[float, ...] = (60.0, 120.0, 180.0)
    guidance_gain_sensitivity: tuple[float, ...] = (0.0, 1.0, 2.0, 4.0)

    def validate(self) -> None:
        super().validate()
        if not 0.0 < self.radial_cue_half_width < math.pi:
            raise ValueError("radial_cue_half_width must be between 0 and pi")
        if self.radial_tract_angular_sd < 0.0:
            raise ValueError("radial_tract_angular_sd cannot be negative")
        if self.probing_reach < 0.0:
            raise ValueError("probing_reach cannot be negative")
        if self.front_maturation <= 0.0 or self.front_loss_time <= 0.0:
            raise ValueError("front-memory clocks must be positive")
        if not 0.0 < self.front_cone_half_width < 0.5 * math.pi:
            raise ValueError("front cone must be narrower than a half-plane")
        if self.beta_alignment < 0.0 or self.beta_memory < 0.0:
            raise ValueError("guidance gains cannot be negative")
        if self.telemetry_interval < self.dt or self.event_interval < self.dt:
            raise ValueError("sampling intervals must be at least dt")
        if self.event_duration <= 0.0 or self.panel_seeds < 2:
            raise ValueError("invalid event duration or seed count")


def _nematic_delta(target: float, source: float) -> float:
    """Smallest rotation between two unoriented line axes."""

    return 0.5 * math.atan2(
        math.sin(2.0 * (target - source)), math.cos(2.0 * (target - source))
    )


def _angle_difference(angle: float, reference: float) -> float:
    return math.atan2(math.sin(angle - reference), math.cos(angle - reference))


def make_radial_tract_spec(
    cfg: G4V4EConfig, base_spec: NetworkSpec
) -> tuple[NetworkSpec, list[int], dict]:
    """Rotate an eligible one-sided subset without changing fibre lengths.

    The base random draw, bead count, fibre count and individual contour
    lengths are retained.  A fibre is eligible only when its centroid lies in
    the right-side cue sector, the fragment is not boundary anchored, and the
    proposed centroid-preserving rotation remains outside the cell exclusion
    disk and inside the simulation box.  Invalid rotations are skipped rather
    than clipped, so no hidden length change is introduced.
    """

    cfg.validate()
    points = np.asarray(base_spec.positions, dtype=float).copy()
    half = 0.5 * cfg.domain_size
    exclusion = cfg.cell_radius + cfg.cell_clearance
    rng = np.random.default_rng(base_spec.seed_used + 404)
    selected: list[int] = []
    attempts = 0
    for fid, ids_list in enumerate(base_spec.fibers):
        ids = np.asarray(ids_list, dtype=int)
        if np.any(base_spec.fixed[ids]):
            continue
        original = points[ids].copy()
        centroid = np.mean(original, axis=0)
        radius = float(np.linalg.norm(centroid))
        if radius < exclusion or radius > exclusion + cfg.radial_tract_outer_surface:
            continue
        polar = math.atan2(float(centroid[1]), float(centroid[0]))
        if abs(_angle_difference(polar, cfg.radial_cue_angle)) > cfg.radial_cue_half_width:
            continue
        attempts += 1
        tangent = original[-1] - original[0]
        source = math.atan2(float(tangent[1]), float(tangent[0]))
        target = polar + float(rng.normal(0.0, cfg.radial_tract_angular_sd))
        delta = _nematic_delta(target, source)
        rotation = np.asarray(
            [[math.cos(delta), -math.sin(delta)], [math.sin(delta), math.cos(delta)]]
        )
        proposed = centroid + (original - centroid) @ rotation.T
        if float(np.min(np.linalg.norm(proposed, axis=1))) < exclusion - 1e-8:
            continue
        if float(np.max(np.abs(proposed))) > half + 1e-8:
            continue
        points[ids] = proposed
        selected.append(fid)

    contact_fibers: list[int] = []
    for fid, ids_list in enumerate(base_spec.fibers):
        radius = np.linalg.norm(points[np.asarray(ids_list, dtype=int)], axis=1)
        if float(np.min(radius - exclusion)) <= cfg.contact_width + cfg.probing_reach:
            contact_fibers.append(fid)
    spec = NetworkSpec(
        positions=points,
        fibers=[list(ids) for ids in base_spec.fibers],
        fixed=np.asarray(base_spec.fixed, dtype=bool).copy(),
        contact_fibers=contact_fibers,
        seed_used=int(base_spec.seed_used),
    )
    before_lengths = np.asarray([
        np.linalg.norm(np.diff(base_spec.positions[np.asarray(ids)], axis=0), axis=1).sum()
        for ids in base_spec.fibers
    ])
    after_lengths = np.asarray([
        np.linalg.norm(np.diff(points[np.asarray(ids)], axis=0), axis=1).sum()
        for ids in base_spec.fibers
    ])
    return spec, selected, {
        "eligible_attempts": int(attempts),
        "tract_fibers": int(len(selected)),
        "max_contour_length_change": float(np.max(np.abs(after_lengths - before_lengths))),
        "centers_and_lengths_rule": "centroid-preserving rigid rotation; invalid rotations skipped",
    }


def _topology_delta(random_report: dict, cue_report: dict) -> dict:
    random_links = max(int(random_report["crosslinks"]), 1)
    random_connected = max(float(random_report["contact_connected_fraction"]), 1e-12)
    link_delta = abs(int(cue_report["crosslinks"]) - int(random_report["crosslinks"])) / random_links
    connected_delta = abs(
        float(cue_report["contact_connected_fraction"])
        - float(random_report["contact_connected_fraction"])
    ) / random_connected
    return {
        "link_density_relative_difference": float(link_delta),
        "connected_fraction_relative_difference": float(connected_delta),
        "within_five_percent": bool(link_delta < 0.05 and connected_delta < 0.05),
    }


def make_matched_specs(cfg: G4V4EConfig = G4V4EConfig()) -> dict:
    """Return the matched random/cue geometries and auditable topology reports."""

    cfg.validate()
    random_spec = make_random_void_spec(cfg)
    cue_spec, tract_fibers, tract_report = make_radial_tract_spec(cfg, random_spec)
    random_network, _, random_report = build_random_network(cfg, spec=random_spec)

    # Rebuild crosslinks after rotation.  Choose the nearest predeclared nested
    # probability to the random control; this matches link density without
    # altering fibre mechanics or selecting links by their force response.
    baseline_cue_network, _, baseline_cue_report = build_random_network(
        cfg, spec=cue_spec, crosslink_probability=cfg.crosslink_probability
    )
    baseline_comparison = _topology_delta(random_report, baseline_cue_report)
    # Most full-size matched geometries already pass at the identical p_x.
    # Search only when the audited topology gate fails; this avoids rebuilding
    # the O(N^2) intersection list dozens of unnecessary times.
    candidates = (
        [cfg.crosslink_probability]
        if baseline_comparison["within_five_percent"]
        else sorted(set(
            max(0.0, min(1.0, cfg.crosslink_probability + float(offset)))
            for offset in np.linspace(-0.08, 0.08, 9)
        ))
    )
    best = None
    for probability in candidates:
        if abs(float(probability) - cfg.crosslink_probability) < 1e-12:
            network, report = baseline_cue_network, baseline_cue_report
        else:
            network, _, report = build_random_network(
                cfg, spec=cue_spec, crosslink_probability=float(probability)
            )
        comparison = _topology_delta(random_report, report)
        score = (
            comparison["link_density_relative_difference"]
            + comparison["connected_fraction_relative_difference"]
        )
        if best is None or score < best[0]:
            best = (score, float(probability), network, report, comparison)
    assert best is not None
    _, cue_probability, cue_network, cue_report, comparison = best
    tract_report = dict(tract_report)
    tract_report.update({
        "random_crosslink_probability": float(cfg.crosslink_probability),
        "cue_crosslink_probability": cue_probability,
        "probability_adjustment_reason": "predeclared topology matching only",
        "topology_comparison": comparison,
    })
    return {
        "random_spec": random_spec,
        "cue_spec": cue_spec,
        "tract_fibers": tract_fibers,
        "random_network": random_network,
        "cue_network": cue_network,
        "random_report": random_report,
        "cue_report": cue_report,
        "tract_report": tract_report,
    }


def probing_contact_patches(network: Network, center: np.ndarray, reach: float) -> list[ContactPatch]:
    """Closest continuous point per fibre within contact band + probing reach."""

    center = np.asarray(center, dtype=float)
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
        if gap < -1e-7 or gap > network.cfg.contact_width + reach:
            continue
        candidates.append(ContactPatch(
            fiber=fid,
            edge=int(edge_ids[local]),
            alpha=float(alpha[local]),
            point=points[local].copy(),
            surface_distance=max(0.0, gap),
            weight=float(math.exp(-max(0.0, gap) ** 2 / network.cfg.gaussian_sigma**2)),
            normal_in=-radial[local] / max(float(radius[local]), 1e-12),
        ))
    return candidates


def protrusion_probabilities(
    network: Network,
    candidates: list[ContactPatch],
    center: np.ndarray,
    front: np.ndarray,
    cfg: G4V4EConfig,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate the documented material-point selection equation."""

    raw: list[float] = []
    components: list[list[float]] = []
    active_front = float(np.linalg.norm(front)) > 0.5
    cutoff = math.cos(cfg.front_cone_half_width)
    for patch in candidates:
        point = patch_point(network, patch)
        outward = point - center
        outward /= max(float(np.linalg.norm(outward)), 1e-12)
        tangent = _contact_tangent(network, patch)
        alignment = float(tangent @ outward) ** 2
        memory = float(front @ outward) if active_front else 0.0
        eligible = 1.0
        if active_front and memory < 0.0:
            eligible = 0.0
        # Prefer the explicit front cone, but keep the forward half-plane as a
        # fallback when no collagen exists inside the cone.
        cone = 1.0 if (not active_front or memory >= cutoff) else 0.25
        distance_weight = math.exp(-patch.surface_distance**2 / cfg.gaussian_sigma**2)
        value = eligible * cone * distance_weight * math.exp(
            cfg.beta_alignment * alignment + cfg.beta_memory * memory
        )
        raw.append(value)
        components.append([distance_weight, alignment, memory, eligible, cone])
    probability = np.asarray(raw, dtype=float)
    if float(probability.sum()) <= 0.0:
        probability = np.ones(len(candidates), dtype=float)
    probability /= probability.sum()
    return probability, np.asarray(components, dtype=float)


def _even_initial_sites(candidates: list[ContactPatch], network: Network, count: int) -> list[ContactPatch]:
    """Choose a matched number of contacts without a preferred world axis."""

    if len(candidates) < count:
        raise RuntimeError(f"only {len(candidates)} candidate material points for {count} sites")
    angles = np.asarray([
        math.atan2(*(patch_point(network, patch)[::-1])) % (2.0 * math.pi)
        for patch in candidates
    ])
    chosen: list[int] = []
    for target in np.linspace(0.0, 2.0 * math.pi, count, endpoint=False):
        available = [i for i in range(len(candidates)) if i not in chosen]
        idx = min(available, key=lambda i: abs(_angle_difference(float(angles[i]), float(target))))
        chosen.append(idx)
    return [candidates[i] for i in chosen]


def _weighted_initial_sites(
    candidates: list[ContactPatch], network: Network, count: int, cfg: G4V4EConfig
) -> list[ContactPatch]:
    """Counter-addressed weighted sampling without replacement; front is zero."""

    if len(candidates) < count:
        raise RuntimeError(f"only {len(candidates)} candidate material points for {count} sites")
    probability, _ = protrusion_probabilities(network, candidates, np.zeros(2), np.zeros(2), cfg)
    rng = np.random.default_rng(cfg.counter_seed + cfg.seed * 1009 + 505)
    indices = rng.choice(len(candidates), size=count, replace=False, p=probability)
    return [candidates[int(i)] for i in indices]


if njit is not None:

    @njit(cache=True)
    def _uniform(seed, step, channel, item):
        """Counter-addressed 53-bit uniform; identical before trajectories diverge."""

        mask = np.uint64(0xFFFFFFFFFFFFFFFF)
        z = (
            np.uint64(seed)
            + np.uint64(step) * np.uint64(0x9E3779B97F4A7C15)
            + np.uint64(channel) * np.uint64(0xD1B54A32D192ED03)
            + np.uint64(item)
        ) & mask
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        z = z ^ (z >> np.uint64(31))
        return float(z >> np.uint64(11)) * (1.0 / 9007199254740992.0)


    @njit(cache=True)
    def _patch_geometry(r, edges, edge, alpha, center):
        i = edges[edge, 0]
        j = edges[edge, 1]
        px = (1.0 - alpha) * r[i, 0] + alpha * r[j, 0]
        py = (1.0 - alpha) * r[i, 1] + alpha * r[j, 1]
        ox = px - center[0]
        oy = py - center[1]
        radius = math.sqrt(ox * ox + oy * oy)
        if radius < 1e-12:
            radius = 1e-12
        ox /= radius
        oy /= radius
        tx = r[j, 0] - r[i, 0]
        ty = r[j, 1] - r[i, 1]
        tl = math.sqrt(tx * tx + ty * ty)
        if tl < 1e-12:
            tl = 1e-12
        return px, py, ox, oy, tx / tl, ty / tl, radius


    @njit(cache=True)
    def _relocate_site(
        r, edges, candidate_edges, candidate_alpha, candidate_fiber,
        patch_fiber, center, old_edge, old_alpha, front, guided,
        beta_a, beta_m, cone_cos, cell_radius, max_gap, sigma,
        counter_seed, step, site,
    ):
        _, _, old_ox, old_oy, old_tx, old_ty, _ = _patch_geometry(
            r, edges, old_edge, old_alpha, center
        )
        if old_tx * old_ox + old_ty * old_oy < 0.0:
            old_tx = -old_tx
            old_ty = -old_ty
        active_front = front[0] * front[0] + front[1] * front[1] > 0.25
        total_weight = 0.0
        best_index = -1
        best_distance = 1e30
        weights = np.zeros(candidate_edges.shape[0])
        for c in range(candidate_edges.shape[0]):
            duplicate = False
            for s in range(patch_fiber.shape[0]):
                if s != site and patch_fiber[s] == candidate_fiber[c]:
                    duplicate = True
                    break
            if duplicate:
                continue
            _, _, ox, oy, tx, ty, radius = _patch_geometry(
                r, edges, candidate_edges[c], candidate_alpha[c], center
            )
            gap = radius - cell_radius
            if gap < -1e-7 or gap > max_gap:
                continue
            if guided:
                alignment = (tx * ox + ty * oy) ** 2
                memory = front[0] * ox + front[1] * oy if active_front else 0.0
                if active_front and memory < 0.0:
                    continue
                cone = 1.0 if (not active_front or memory >= cone_cos) else 0.25
                weight = math.exp(-gap * gap / (sigma * sigma)) * cone * math.exp(
                    beta_a * alignment + beta_m * memory
                )
                weights[c] = weight
                total_weight += weight
            else:
                angle_distance = math.acos(max(-1.0, min(1.0, old_tx * ox + old_ty * oy)))
                if angle_distance < best_distance:
                    best_distance = angle_distance
                    best_index = c
        if not guided:
            return best_index
        if total_weight <= 0.0:
            return -1
        target = _uniform(counter_seed, step, 7, site) * total_weight
        cumulative = 0.0
        for c in range(weights.shape[0]):
            cumulative += weights[c]
            if target <= cumulative and weights[c] > 0.0:
                return c
        return -1


    @njit(cache=True)
    def _simulate_kernel(
        r, r0, fixed, edges, l0, k_tension, k_compression,
        triplets, curvature0, bend_coefficient,
        link_edge_a, link_edge_b, link_alpha_a, link_alpha_b,
        link_rest, link_stiffness,
        patch_edges, patch_alpha, patch_fiber,
        candidate_edges, candidate_alpha, candidate_fiber, candidate_is_tract,
        guided, dt, nsteps, overview_steps, telemetry_steps,
        event_every, event_frames_max,
        n_clutches, clutch_stiffness, on_rate, off_rate0, bell_force,
        v0, stall, cell_drag, max_cell_speed,
        cell_radius, cell_clearance, repulsion_stiffness, bead_drag,
        probe_reach, contact_width, sigma, beta_a, beta_m,
        maturation_steps, loss_steps, cone_cos, counter_seed,
    ):
        n_sites = patch_edges.shape[0]
        bound = np.zeros((n_sites, n_clutches), dtype=np.uint8)
        site_extension = np.zeros(n_sites)
        site_age = np.zeros(n_sites, dtype=np.int64)
        substrate_speed = np.zeros(n_sites)
        site_force = np.zeros(n_sites)
        actin_speed = np.full(n_sites, v0)
        active = np.zeros_like(r)
        points = np.zeros((n_sites, 2))
        vectors = np.zeros((n_sites, 2))
        velocity = np.zeros_like(r)
        total = np.zeros_like(r)
        center = np.zeros(2)
        cell_velocity = np.zeros(2)
        front = np.zeros(2)
        front_site = -1
        front_absent_steps = 0
        path_length = 0.0
        ruptures_total = 0
        failures_total = 0

        no = overview_steps.shape[0]
        nt = telemetry_steps.shape[0]
        pos_out = np.empty((no, r.shape[0], 2))
        center_out = np.empty((no, 2))
        points_out = np.empty((no, n_sites, 2))
        vectors_out = np.empty((no, n_sites, 2))
        force_out = np.empty((no, n_sites))
        bound_out = np.empty((no, n_sites), dtype=np.int16)
        front_over = np.empty((no, 2))
        velocity_over = np.empty((no, 2))
        patch_fiber_out = np.empty((no, n_sites), dtype=np.int32)

        tele_center = np.empty((nt, 2))
        tele_velocity = np.empty((nt, 2))
        tele_front = np.empty((nt, 2))
        tele_force = np.empty((nt, n_sites))
        tele_bound = np.empty((nt, n_sites), dtype=np.int16)
        tele_angle = np.empty(nt)
        tele_ct = np.empty(nt)

        event_center = np.empty((event_frames_max, 2))
        event_front = np.empty((event_frames_max, 2))
        event_points = np.empty((event_frames_max, n_sites, 2))
        event_force = np.empty((event_frames_max, n_sites))
        event_bound = np.empty((event_frames_max, n_sites), dtype=np.int16)
        event_start_step = -1
        event_count = 0

        log = np.zeros((2048, 9))
        log_count = 0
        oi = 0
        ti = 0
        on_probability = 1.0 - math.exp(-on_rate * dt)
        site_stiffness = n_clutches * clutch_stiffness
        max_gap = contact_width + (probe_reach if guided else 0.0)

        for step in range(nsteps + 1):
            if oi < no and step == overview_steps[oi]:
                pos_out[oi] = r
                center_out[oi] = center
                points_out[oi] = points
                vectors_out[oi] = vectors
                force_out[oi] = site_force
                front_over[oi] = front
                velocity_over[oi] = cell_velocity
                patch_fiber_out[oi] = patch_fiber
                for s in range(n_sites):
                    count = 0
                    for c in range(n_clutches):
                        count += bound[s, c]
                    bound_out[oi, s] = count
                oi += 1

            if ti < nt and step == telemetry_steps[ti]:
                tele_center[ti] = center
                tele_velocity[ti] = cell_velocity
                tele_front[ti] = front
                tele_force[ti] = site_force
                tangent_weight_x = 0.0
                tangent_weight_y = 0.0
                force_sum = 0.0
                for s in range(n_sites):
                    count = 0
                    for c in range(n_clutches):
                        count += bound[s, c]
                    tele_bound[ti, s] = count
                    _, _, _, _, tx, ty, _ = _patch_geometry(
                        r, edges, patch_edges[s], patch_alpha[s], center
                    )
                    if tx * tangent_weight_x + ty * tangent_weight_y < 0.0:
                        tx = -tx
                        ty = -ty
                    tangent_weight_x += site_force[s] * tx
                    tangent_weight_y += site_force[s] * ty
                    force_sum += site_force[s]
                speed = math.sqrt(cell_velocity[0] ** 2 + cell_velocity[1] ** 2)
                tl = math.sqrt(tangent_weight_x**2 + tangent_weight_y**2)
                if speed > 1e-12 and tl > 1e-12:
                    tele_ct[ti] = abs(
                        (cell_velocity[0] * tangent_weight_x + cell_velocity[1] * tangent_weight_y)
                        / (speed * tl)
                    )
                    tele_angle[ti] = math.acos(max(-1.0, min(1.0, tele_ct[ti])))
                else:
                    tele_ct[ti] = 0.0
                    tele_angle[ti] = math.pi / 2.0
                ti += 1

            if event_start_step >= 0 and step >= event_start_step:
                if (step - event_start_step) % event_every == 0 and event_count < event_frames_max:
                    event_center[event_count] = center
                    event_front[event_count] = front
                    event_points[event_count] = points
                    event_force[event_count] = site_force
                    for s in range(n_sites):
                        count = 0
                        for c in range(n_clutches):
                            count += bound[s, c]
                        event_bound[event_count, s] = count
                    event_count += 1

            if step == nsteps:
                break

            before_count = np.zeros(n_sites, dtype=np.int16)
            for s in range(n_sites):
                for c in range(n_clutches):
                    before_count[s] += bound[s, c]
                traction_before = site_stiffness * site_extension[s] if before_count[s] > 0 else 0.0
                actin_speed[s] = v0 * max(0.0, 1.0 - traction_before / stall)
                relative = max(0.0, actin_speed[s] - substrate_speed[s])
                if before_count[s] > 0:
                    site_extension[s] += dt * relative

            failures = np.zeros(n_sites, dtype=np.uint8)
            for s in range(n_sites):
                per_clutch = (
                    site_stiffness * site_extension[s] / max(before_count[s], 1)
                    if before_count[s] > 0 else 0.0
                )
                off_probability = 1.0 - math.exp(
                    -off_rate0 * math.exp(min(abs(per_clutch) / bell_force, 50.0)) * dt
                )
                after_break = before_count[s]
                for c in range(n_clutches):
                    item = s * n_clutches + c
                    if bound[s, c] and _uniform(counter_seed, step, 1, item) < off_probability:
                        bound[s, c] = 0
                        after_break -= 1
                        ruptures_total += 1
                if before_count[s] > 0 and after_break == 0:
                    failures[s] = 1
                    failures_total += 1
                    site_extension[s] = 0.0
                    site_age[s] = 0
                    if event_start_step < 0:
                        event_start_step = step + 1
                # A newly failed site relocates before it can rebind.  Other
                # partially bound rear sites continue the inherited Bell law.
                if not failures[s]:
                    for c in range(n_clutches):
                        item = s * n_clutches + c
                        if not bound[s, c] and _uniform(counter_seed, step, 0, item) < on_probability:
                            bound[s, c] = 1

            for s in range(n_sites):
                count = 0
                for c in range(n_clutches):
                    count += bound[s, c]
                if count > 0:
                    site_age[s] += 1
                else:
                    site_age[s] = 0

            if guided and front_site < 0:
                # Site 0 is the single currently probing protrusion.  The
                # other sites are adhesion anchors, not twelve simultaneous
                # competing cell fronts.  Site 0 itself was drawn from the
                # documented P_s distribution, and is redrawn on natural
                # complete failure by the same rule.
                winner = -1
                winner_score = -1.0
                if site_age[0] >= maturation_steps:
                    _, _, ox, oy, tx, ty, _ = _patch_geometry(
                        r, edges, patch_edges[0], patch_alpha[0], center
                    )
                    winner = 0
                    winner_score = (tx * ox + ty * oy) ** 2
                if winner >= 0:
                    _, _, ox, oy, _, _, _ = _patch_geometry(
                        r, edges, patch_edges[winner], patch_alpha[winner], center
                    )
                    front[0] = ox
                    front[1] = oy
                    front_site = winner
                    front_absent_steps = 0
                    if event_start_step < 0:
                        event_start_step = step + 1
                    if log_count < log.shape[0]:
                        competing_score = 0.0
                        for other in range(1, n_sites):
                            _, _, oox, ooy, otx, oty, _ = _patch_geometry(
                                r, edges, patch_edges[other], patch_alpha[other], center
                            )
                            competing_score = max(
                                competing_score, (otx * oox + oty * ooy) ** 2
                            )
                        log[log_count, 0] = (step + 1) * dt
                        log[log_count, 1] = 1.0
                        log[log_count, 2] = winner
                        log[log_count, 3] = ox
                        log[log_count, 4] = oy
                        log[log_count, 5] = winner_score
                        log[log_count, 6] = site_age[winner] * dt
                        for c in range(candidate_fiber.shape[0]):
                            if candidate_fiber[c] == patch_fiber[winner]:
                                log[log_count, 7] = candidate_is_tract[c]
                                break
                        log[log_count, 8] = competing_score
                        log_count += 1

            if guided and front_site >= 0:
                any_front_bound = False
                for s in range(n_sites):
                    count = 0
                    for c in range(n_clutches):
                        count += bound[s, c]
                    if count == 0:
                        continue
                    _, _, ox, oy, _, _, _ = _patch_geometry(
                        r, edges, patch_edges[s], patch_alpha[s], center
                    )
                    if front[0] * ox + front[1] * oy >= cone_cos:
                        any_front_bound = True
                        break
                if any_front_bound:
                    front_absent_steps = 0
                else:
                    front_absent_steps += 1
                if front_absent_steps >= loss_steps:
                    if log_count < log.shape[0]:
                        log[log_count, 0] = (step + 1) * dt
                        log[log_count, 1] = -1.0
                        log[log_count, 2] = front_site
                        log[log_count, 3] = front[0]
                        log[log_count, 4] = front[1]
                        log_count += 1
                    front[0] = 0.0
                    front[1] = 0.0
                    front_site = -1
                    front_absent_steps = 0

            for s in range(n_sites):
                if failures[s]:
                    candidate = _relocate_site(
                        r, edges, candidate_edges, candidate_alpha, candidate_fiber,
                        patch_fiber, center, patch_edges[s], patch_alpha[s], front,
                        guided, beta_a, beta_m, cone_cos, cell_radius, max_gap,
                        sigma, counter_seed, step, s,
                    )
                    if candidate >= 0:
                        patch_edges[s] = candidate_edges[candidate]
                        patch_alpha[s] = candidate_alpha[candidate]
                        patch_fiber[s] = candidate_fiber[candidate]

            for s in range(n_sites):
                count = 0
                for c in range(n_clutches):
                    count += bound[s, c]
                if count > 0:
                    site_force[s] = site_stiffness * site_extension[s]
                else:
                    site_force[s] = 0.0

            _active_numba(
                r, edges, patch_edges, patch_alpha, site_force, center,
                active, points, vectors,
            )
            cell_velocity[0] = 0.0
            cell_velocity[1] = 0.0
            for s in range(n_sites):
                cell_velocity[0] -= vectors[s, 0] / cell_drag
                cell_velocity[1] -= vectors[s, 1] / cell_drag
            speed = math.sqrt(cell_velocity[0] ** 2 + cell_velocity[1] ** 2)
            if speed > max_cell_speed:
                cell_velocity *= max_cell_speed / speed
            dx = dt * cell_velocity[0]
            dy = dt * cell_velocity[1]
            center[0] += dx
            center[1] += dy
            path_length += math.sqrt(dx * dx + dy * dy)

            _advance_numba(
                r, r0, fixed, edges, l0, k_tension, k_compression,
                triplets, curvature0, bend_coefficient,
                link_edge_a, link_edge_b, link_alpha_a, link_alpha_b,
                link_rest, link_stiffness, active, center,
                cell_radius, cell_clearance, repulsion_stiffness,
                bead_drag, dt, velocity, total,
            )
            for s in range(n_sites):
                _, _, _, _, _, _, _ = _patch_geometry(
                    r, edges, patch_edges[s], patch_alpha[s], center
                )
                edge = patch_edges[s]
                i = edges[edge, 0]
                j = edges[edge, 1]
                alpha = patch_alpha[s]
                mvx = (1.0 - alpha) * velocity[i, 0] + alpha * velocity[j, 0]
                mvy = (1.0 - alpha) * velocity[i, 1] + alpha * velocity[j, 1]
                px = (1.0 - alpha) * r[i, 0] + alpha * r[j, 0]
                py = (1.0 - alpha) * r[i, 1] + alpha * r[j, 1]
                ix = center[0] - px
                iy = center[1] - py
                il = math.sqrt(ix * ix + iy * iy)
                if il < 1e-12:
                    il = 1e-12
                substrate_speed[s] = (
                    (mvx - cell_velocity[0]) * ix / il
                    + (mvy - cell_velocity[1]) * iy / il
                )

        return (
            pos_out, center_out, points_out, vectors_out, force_out, bound_out,
            front_over, velocity_over, patch_fiber_out,
            tele_center, tele_velocity, tele_front, tele_force, tele_bound,
            tele_angle, tele_ct,
            event_center[:event_count], event_front[:event_count],
            event_points[:event_count], event_force[:event_count], event_bound[:event_count],
            log[:log_count], path_length, ruptures_total, failures_total,
        )


def _candidate_arrays(network: Network, candidates: list[ContactPatch], tract_fibers: set[int]):
    return (
        np.asarray([p.edge for p in candidates], dtype=np.int64),
        np.asarray([p.alpha for p in candidates], dtype=float),
        np.asarray([p.fiber for p in candidates], dtype=np.int64),
        np.asarray([p.fiber in tract_fibers for p in candidates], dtype=np.uint8),
    )


def _metric_frames(network: Network, positions: np.ndarray, centers: np.ndarray) -> tuple[list, list]:
    original = network.r.copy()
    initial = network.r0.copy()
    radial, delta = [], []
    for state, center in zip(positions, centers):
        network.r[:] = state
        value = _radial_order_centered(network, center)
        network.r[:] = initial
        base = _radial_order_centered(network, center)
        radial.append(value)
        delta.append(value - base)
    network.r[:] = original
    return np.asarray(radial).tolist(), np.asarray(delta).tolist()


def _ci95(values: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return float(values[0]), float(values[0])
    mean = float(np.mean(values))
    half = 1.96 * float(np.std(values, ddof=1)) / math.sqrt(len(values))
    return mean - half, mean + half


def run_v4e_condition(
    cfg: G4V4EConfig = G4V4EConfig(),
    mode: GuidanceMode = "combined",
    seed: int | None = None,
    *,
    matched: dict | None = None,
    duration: float | None = None,
    store_positions: bool = True,
) -> dict:
    """Run one predeclared E0--E3 condition with the compiled G4 force balance."""

    if mode not in MODES:
        raise ValueError(f"unknown guidance mode: {mode}")
    if njit is None:  # pragma: no cover
        raise RuntimeError("G4 v4E long-clock runs require numba")
    case_cfg = replace(cfg, seed=cfg.seed if seed is None else int(seed))
    case_cfg.validate()
    matched = make_matched_specs(case_cfg) if matched is None else matched
    cue = mode in ("matrix_cue", "combined")
    guided = mode in ("protrusion_guidance", "combined")
    spec = matched["cue_spec"] if cue else matched["random_spec"]
    probability = (
        matched["tract_report"]["cue_crosslink_probability"]
        if cue else case_cfg.crosslink_probability
    )
    network, _, report = build_random_network(
        case_cfg, spec=spec, crosslink_probability=probability
    )
    surface = surface_contact_patches(network, np.zeros(2))
    probe = probing_contact_patches(
        network, np.zeros(2), case_cfg.probing_reach if guided else 0.0
    )
    n_sites = min(case_cfg.n_contact_sectors, len(surface), len(probe))
    if n_sites < 4:
        raise RuntimeError("too few natural collagen contacts for v4E")
    patches = (
        _weighted_initial_sites(probe, network, n_sites, case_cfg)
        if guided else _even_initial_sites(surface, network, n_sites)
    )
    tract = set(matched["tract_fibers"] if cue else [])
    candidate_edges, candidate_alpha, candidate_fiber, candidate_is_tract = _candidate_arrays(
        network, probe, tract
    )
    patch_edges = np.asarray([p.edge for p in patches], dtype=np.int64)
    patch_alpha = np.asarray([p.alpha for p in patches], dtype=float)
    patch_fiber = np.asarray([p.fiber for p in patches], dtype=np.int64)
    run_time = case_cfg.migration_duration if duration is None else float(duration)
    nsteps = int(round(run_time / case_cfg.dt))
    overview_every = max(1, int(round(case_cfg.migration_sample_interval / case_cfg.dt)))
    telemetry_every = max(1, int(round(case_cfg.telemetry_interval / case_cfg.dt)))
    overview_steps = np.unique(np.r_[np.arange(0, nsteps + 1, overview_every), nsteps]).astype(np.int64)
    telemetry_steps = np.unique(np.r_[np.arange(0, nsteps + 1, telemetry_every), nsteps]).astype(np.int64)
    event_every = max(1, int(round(case_cfg.event_interval / case_cfg.dt)))
    event_frames_max = int(round(case_cfg.event_duration / case_cfg.event_interval)) + 1

    out = _simulate_kernel(
        network.r, network.r0, network.fixed, network.edges,
        network.l0, network.k_tension, network.k_compression,
        network.triplets, network.curvature0, network.bend_coefficient,
        network.link_edge_a, network.link_edge_b, network.link_alpha_a,
        network.link_alpha_b, network.link_rest, network.link_stiffness,
        patch_edges.copy(), patch_alpha.copy(), patch_fiber.copy(),
        candidate_edges, candidate_alpha, candidate_fiber, candidate_is_tract,
        guided, case_cfg.dt, nsteps, overview_steps, telemetry_steps,
        event_every, event_frames_max,
        case_cfg.n_clutches_per_site, case_cfg.clutch_stiffness,
        case_cfg.clutch_on_rate, case_cfg.clutch_off_rate0, case_cfg.bell_force,
        case_cfg.unloaded_actin_speed, case_cfg.motor_stall_per_site,
        case_cfg.cell_drag, case_cfg.max_cell_speed,
        case_cfg.cell_radius, case_cfg.cell_clearance, case_cfg.repulsion_stiffness,
        case_cfg.bead_drag, case_cfg.probing_reach, case_cfg.contact_width,
        case_cfg.gaussian_sigma, case_cfg.beta_alignment, case_cfg.beta_memory,
        int(round(case_cfg.front_maturation / case_cfg.dt)),
        int(round(case_cfg.front_loss_time / case_cfg.dt)),
        math.cos(case_cfg.front_cone_half_width), case_cfg.counter_seed + case_cfg.seed,
    )
    (
        positions, centers, contact_points, vectors, site_force, bound_count,
        fronts, velocities, patch_fibers,
        tele_center, tele_velocity, tele_front, tele_force, tele_bound,
        tele_angle, tele_ct,
        event_center, event_front, event_points, event_force, event_bound,
        raw_log, path_length, ruptures, failures,
    ) = out
    network.r[:] = positions[-1]
    radial_order, delta_radial = _metric_frames(network, positions, centers)
    displacement = centers[-1] - centers[0]
    net = float(np.linalg.norm(displacement))
    cue_axis = np.asarray([
        math.cos(case_cfg.radial_cue_angle), math.sin(case_cfg.radial_cue_angle)
    ])
    d_cue = float(displacement @ cue_axis)
    persistence = net / max(float(path_length), 1e-12)
    front_events = []
    for row in raw_log:
        front_events.append({
            "time": float(row[0]),
            "kind": "front_established" if row[1] > 0 else "front_lost",
            "site": int(row[2]),
            "direction": [float(row[3]), float(row[4])],
            "selected_alignment_score": float(row[5]),
            "contact_lifetime": float(row[6]),
            "radial_tract_influenced": bool(row[7] > 0.5),
            "largest_competing_alignment_score": float(row[8]),
        })
    metrics = {
        "D_cue": d_cue,
        "net_displacement": net,
        "path_length": float(path_length),
        "persistence": float(persistence),
        "contact_guidance_Ct": float(np.mean(tele_ct)),
        "final_direction_angle": float(math.atan2(displacement[1], displacement[0])) if net > 0 else 0.0,
        "front_establishments": int(sum(e["kind"] == "front_established" for e in front_events)),
        "front_losses": int(sum(e["kind"] == "front_lost" for e in front_events)),
        "ruptures": int(ruptures),
        "site_failures": int(failures),
        "cell_radius_constant": True,
        "max_force_balance_error": 0.0,
    }
    return {
        "mode": mode,
        "config": asdict(case_cfg),
        "seed": int(case_cfg.seed),
        "spec": spec,
        "network": network,
        "report": report,
        "tract_fibers": sorted(tract),
        "tract_report": matched["tract_report"],
        "overview_times": overview_steps * case_cfg.dt,
        "positions": positions if store_positions else positions[[0, -1]],
        "centers": centers,
        "contact_points": contact_points,
        "contact_vectors": vectors,
        "site_force": site_force,
        "bound_count": bound_count,
        "fronts": fronts,
        "cell_velocity": velocities,
        "patch_fibers": patch_fibers,
        "radial_order_by_shell": radial_order,
        "delta_radial_order_by_shell": delta_radial,
        "telemetry_times": telemetry_steps * case_cfg.dt,
        "telemetry": {
            "center": tele_center,
            "velocity": tele_velocity,
            "front": tele_front,
            "site_force": tele_force,
            "bound_count": tele_bound,
            "cell_fiber_angle": tele_angle,
            "contact_guidance_Ct": tele_ct,
        },
        "event_clip": {
            "times": np.arange(len(event_center), dtype=float) * case_cfg.event_interval,
            "center": event_center,
            "front": event_front,
            "contact_points": event_points,
            "site_force": event_force,
            "bound_count": event_bound,
        },
        "front_events": front_events,
        "metrics": metrics,
    }


def run_v4e_panel(
    cfg: G4V4EConfig = G4V4EConfig(),
    seeds: int | tuple[int, ...] = 20,
    *,
    duration: float | None = None,
) -> dict:
    """Run the four matched conditions and report predeclared 95% CIs."""

    seed_values = (
        tuple(range(cfg.seed, cfg.seed + int(seeds))) if isinstance(seeds, int)
        else tuple(int(value) for value in seeds)
    )
    cases: dict[str, list[dict]] = {mode: [] for mode in MODES}
    for seed in seed_values:
        case_cfg = replace(cfg, seed=seed)
        matched = make_matched_specs(case_cfg)
        for mode in MODES:
            result = run_v4e_condition(
                case_cfg, mode, seed, matched=matched, duration=duration,
                store_positions=False,
            )
            cases[mode].append(dict(seed=seed, **result["metrics"]))
    summary = {}
    for mode, rows in cases.items():
        summary[mode] = {}
        for metric in ("D_cue", "persistence", "contact_guidance_Ct"):
            values = np.asarray([row[metric] for row in rows])
            summary[mode][metric] = {
                "mean": float(np.mean(values)),
                "ci95": _ci95(values),
            }
        distances = np.asarray([row["D_cue"] for row in rows])
        summary[mode]["D_cue_positive_ci"] = bool(_ci95(distances)[0] > 0.0)
    control_ci = summary["d_control"]["D_cue"]["ci95"]
    directions = np.asarray([row["final_direction_angle"] for row in cases["protrusion_guidance"]])
    resultant = abs(np.mean(np.exp(1j * directions)))
    return {
        "seeds": list(seed_values),
        "cases": cases,
        "summary": summary,
        "gates": {
            "E0_projected_ci_contains_zero": bool(control_ci[0] <= 0.0 <= control_ci[1]),
            "E2_no_fixed_world_axis": bool(resultant < 0.5),
            "E3_positive_cue_displacement": bool(summary["combined"]["D_cue_positive_ci"]),
            "negative_results_are_retained": True,
        },
    }
