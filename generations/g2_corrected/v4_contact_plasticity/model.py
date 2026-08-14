"""Generation 2 / V4: load/unload test with stress-free weak-link formation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from common.model import (  # noqa: E402
    CollagenConfig,
    Crosslink,
    Network,
    NetworkSpec,
    _empty_output,
    _finalise_output,
    _record_frame,
    contact_patches,
    forces_from_patches,
    make_network_spec,
)


@dataclass(frozen=True)
class PlasticityConfig(CollagenConfig):
    bead_drag: float = 360.0
    dt: float = 0.05
    duration: float = 900.0
    sample_interval: float = 6.0
    contact_update_interval: float = 0.50
    force_ramp_time: float = 30.0
    load_duration: float = 300.0
    unload_ramp_time: float = 30.0
    weak_link_stiffness: float = 15.0
    weak_link_capture_radius: float = 0.45
    weak_link_min_approach: float = 0.00025
    weak_link_alignment_deg: float = 30.0
    weak_link_formation_radius: float = 30.0
    weak_link_check_interval: float = 5.0
    max_new_links: int = 36
    max_new_links_per_check: int = 3


def force_protocol(cfg: PlasticityConfig, time: float) -> float:
    """Ramp, hold, ramp down, then fully unload."""

    if time < cfg.force_ramp_time:
        return cfg.total_pull_force * time / cfg.force_ramp_time
    if time < cfg.load_duration:
        return cfg.total_pull_force
    if time < cfg.load_duration + cfg.unload_ramp_time:
        return cfg.total_pull_force * (
            1.0 - (time - cfg.load_duration) / cfg.unload_ramp_time
        )
    return 0.0


def _bead_fibre_map(network: Network) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    fibre = np.full(len(network.r), -1, dtype=int)
    edge = np.full(len(network.r), -1, dtype=int)
    alpha = np.zeros(len(network.r), dtype=float)
    for fid, ids in enumerate(network.fibers):
        edge_ids = np.flatnonzero(network.edge_fiber == fid)
        for local, bead in enumerate(ids):
            fibre[bead] = fid
            if local < len(ids) - 1:
                edge[bead] = edge_ids[local]
                alpha[bead] = 0.0
            else:
                edge[bead] = edge_ids[-1]
                alpha[bead] = 1.0
    return fibre, edge, alpha


def _bead_tangents(network: Network) -> np.ndarray:
    tangent = np.zeros_like(network.r)
    for ids in network.fibers:
        tangent[ids[0]] = network.r[ids[1]] - network.r[ids[0]]
        tangent[ids[-1]] = network.r[ids[-1]] - network.r[ids[-2]]
        tangent[ids[1:-1]] = network.r[ids[2:]] - network.r[ids[:-2]]
    return tangent / np.maximum(np.linalg.norm(tangent, axis=1), 1e-12)[:, None]


def form_new_weak_links(
    network: Network,
    cfg: PlasticityConfig,
    time: float,
    *,
    limit: int | None = None,
    require_new_contact: bool = True,
) -> list[Crosslink]:
    """Join newly adjacent, aligned fibres with stress-free weak links.

    This implements the Ban-style candidate mechanism: loading can bring two
    different fibres into contact, and a new weak crosslink is stress-free at
    the deformed geometry where it formed.  Existing links never reset their
    rest vector and no link breaks in this stage.
    """

    remaining = cfg.max_new_links - sum(x.kind == "new-weak" for x in network.crosslinks)
    if remaining <= 0:
        return []
    limit = min(remaining, cfg.max_new_links_per_check if limit is None else limit)
    fibre, edge, alpha = _bead_fibre_map(network)
    tangent = _bead_tangents(network)
    radial = np.linalg.norm(network.r, axis=1) - cfg.cell_radius
    eligible = np.flatnonzero(
        (radial >= 0.0)
        & (radial <= cfg.weak_link_formation_radius)
        & (~network.fixed)
    )
    cell = cfg.weak_link_capture_radius
    buckets: dict[tuple[int, int], list[int]] = {}
    for bead in eligible:
        key = tuple(np.floor(network.r[bead] / cell).astype(int))
        buckets.setdefault(key, []).append(int(bead))

    existing_locations: list[tuple[int, int, np.ndarray]] = []
    for link in network.crosslinks:
        fa = int(network.edge_fiber[link.edge_a])
        fb = int(network.edge_fiber[link.edge_b])
        pa, pb = network.crosslink_endpoints(link)
        existing_locations.append((min(fa, fb), max(fa, fb), 0.5 * (pa + pb)))

    cosine = math.cos(math.radians(cfg.weak_link_alignment_deg))
    candidates: list[tuple[float, int, int]] = []
    seen: set[tuple[int, int]] = set()
    for key, beads in buckets.items():
        neighbourhood: list[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                neighbourhood.extend(buckets.get((key[0] + dx, key[1] + dy), []))
        for a in beads:
            for b in neighbourhood:
                if a >= b or fibre[a] == fibre[b]:
                    continue
                pair = (a, b)
                if pair in seen:
                    continue
                seen.add(pair)
                distance = float(np.linalg.norm(network.r[b] - network.r[a]))
                if distance > cfg.weak_link_capture_radius:
                    continue
                initial_distance = float(np.linalg.norm(network.r0[b] - network.r0[a]))
                if (
                    require_new_contact
                    and initial_distance - distance < cfg.weak_link_min_approach
                ):
                    continue
                if abs(float(tangent[a] @ tangent[b])) < cosine:
                    continue
                fa, fb = sorted((int(fibre[a]), int(fibre[b])))
                midpoint = 0.5 * (network.r[a] + network.r[b])
                if any(
                    xfa == fa and xfb == fb and np.linalg.norm(midpoint - loc) < 1.0
                    for xfa, xfb, loc in existing_locations
                ):
                    continue
                candidates.append((distance, a, b))
    candidates.sort(key=lambda x: x[0])

    created: list[Crosslink] = []
    used_fibre_pairs: set[tuple[int, int]] = set()
    for _, a, b in candidates:
        fpair = tuple(sorted((int(fibre[a]), int(fibre[b]))))
        if fpair in used_fibre_pairs:
            continue
        pa = network.r[a].copy()
        pb = network.r[b].copy()
        link = Crosslink(
            int(edge[a]),
            float(alpha[a]),
            int(edge[b]),
            float(alpha[b]),
            pb - pa,
            cfg.weak_link_stiffness,
            kind="new-weak",
            created_at=float(time),
        )
        network.crosslinks.append(link)
        created.append(link)
        used_fibre_pairs.add(fpair)
        if len(created) >= limit:
            break
    if created:
        network.refresh_crosslink_arrays()
    return created


def run_load_unload(
    cfg: PlasticityConfig,
    *,
    spec: NetworkSpec,
    allow_new_links: bool,
) -> dict:
    network = Network(spec, cfg)
    out = _empty_output()
    out["new_link_count"] = []
    out["protocol_force"] = []
    center = np.zeros(2)
    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    contact_every = max(1, int(round(cfg.contact_update_interval / cfg.dt)))
    link_every = max(1, int(round(cfg.weak_link_check_interval / cfg.dt)))
    patches = contact_patches(network, center, 0.0)
    velocity = np.zeros_like(network.r)
    energy = {"stretch": 0.0, "bend": 0.0, "crosslink": 0.0}
    strain = np.zeros(len(network.edges))
    link_force = np.zeros(len(network.crosslinks))

    for step in range(nsteps + 1):
        time = step * cfg.dt
        protocol = force_protocol(cfg, time)
        if step and step % contact_every == 0 and protocol > 0:
            patches = contact_patches(network, center, 0.0)
        active, vectors = forces_from_patches(network, patches, protocol)
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
            out["new_link_count"].append(
                sum(x.kind == "new-weak" for x in network.crosslinks)
            )
            out["protocol_force"].append(protocol)
        if step == nsteps:
            break
        if (
            allow_new_links
            and protocol > 0.0
            and time >= cfg.force_ramp_time
            and step % link_every == 0
        ):
            form_new_weak_links(network, cfg, time)
            if len(link_force) != len(network.crosslinks):
                link_force = np.pad(
                    link_force, (0, len(network.crosslinks) - len(link_force))
                )
        velocity, energy, strain, link_force, _ = network.advance(
            active, center, cfg.dt, include_crosslinks=True
        )

    final_count = len(network.crosslinks)
    out["crosslink_force"] = [
        np.pad(np.asarray(values), (0, final_count - len(values)))
        for values in out["crosslink_force"]
    ]
    out = _finalise_output(out, network, True)
    out["new_link_count"] = np.asarray(out["new_link_count"])
    out["protocol_force"] = np.asarray(out["protocol_force"])
    out["allow_new_links"] = allow_new_links
    out["load_duration"] = cfg.load_duration
    out["unload_complete_time"] = cfg.load_duration + cfg.unload_ramp_time
    return out


def run_v4_pair(
    cfg: PlasticityConfig = PlasticityConfig(),
    *,
    gates_passed: bool,
) -> dict:
    if not gates_passed:
        raise RuntimeError("V4 is gated until corrected V2 and V3 validations pass")
    spec = make_network_spec(cfg)
    return {
        "elastic": run_load_unload(
            cfg, spec=spec, allow_new_links=False
        ),
        "new_links": run_load_unload(
            cfg, spec=spec, allow_new_links=True
        ),
        "config": asdict(cfg),
        "gates_passed": True,
        "comparison": "same load/unload; stress-free weak-link formation only",
    }
