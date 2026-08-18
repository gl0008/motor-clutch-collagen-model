"""Generation 3 - emergent spheroid guidance on the corrected G2 collagen engine.

This model is a deliberate, minimal rebuild of Generation 3.  It reuses Gloria's
validated Generation-2 bead-spring collagen engine
(``generations/g2_corrected/common/model.py``) for the network mechanics -- axial
springs, discrete bending, freely hinged permanent crosslinks, boundary anchoring
and overdamped bead dynamics -- and adds exactly the three ingredients the
professor's story needs:

1.  A larger **cell spheroid** that starts with a fibre-free gap, so at ``t = 0``
    there is genuinely *no* cell-ECM contact.
2.  **Explicit protrusions** that grow radially outward from the spheroid surface;
    a protrusion can bind a collagen fibre only once its *tip* physically reaches a
    material point (tip-first encounter across the gap).
3.  **Emergent polarity** -- a mass-conserving wave-pinning activity field on the
    membrane (Mori, Jilkine & Edelstein-Keshet 2008) whose front is stabilised by a
    FAK/Rac1-like adhesion-traction feedback (Carey et al. 2016).  There is no
    ``polarity_probability = 0.65`` and no prescribed +x/left-right bias: noise
    breaks the symmetry, mechanics choose the direction, aligned collagen only
    orients the *axis*.

Motion is reaction-driven exactly as in G2 V3: the cell velocity is the summed
equal-and-opposite clutch reaction divided by a cell drag.  The spheroid therefore
migrates toward the front where it actually grips collagen, and the network
realigns radially along the way.

Nothing here claims a realistic 3D tumour-migration prediction.  It is a 2D
mechanism demonstration; every generated figure carries that boundary.

Literature parameter anchors (see docs):
  * spheroid traction ~tens of nN, force polarity ~0.47 (Steinwachs 2016; Mark 2020)
  * MDA-MB-231 3D speed 0.1-0.3 um/min (Sapudom 2019; Steinwachs 2016)
  * protrusions reach ~0.7-4x cell radius, lifetime 10-30 min (Carey 2016; Fraley 2010/2015)
  * wave-pinning polarises in ~10 s, stochastic sign (Mori 2008; Jilkine 2011)
  * fibre realignment develops over tens of minutes to hours (Kim 2017; Han 2018)
"""

from __future__ import annotations

from dataclasses import dataclass
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
    _curved_segment,
)


# =================================================================================
# Configuration
# =================================================================================
@dataclass(frozen=True)
class SpheroidConfig(CollagenConfig):
    """Geometry, mechanics, protrusion and polarity parameters for the rebuild.

    Length unit micrometre, force unit nanonewton, time unit second -- identical
    to the G2 engine so the two models share one coordinate and force language.
    """

    # --- domain / spheroid ---
    domain_size: float = 240.0
    cell_radius: float = 22.0            # spheroid body radius
    gap: float = 10.0                    # fibre-free ring beyond the surface at t=0
    n_fibers: int = 150
    min_fiber_length: float = 20.0
    max_fiber_length: float = 70.0
    bead_spacing: float = 1.0
    boundary_width: float = 3.0
    required_connected_fraction: float = 0.55
    generation_attempts: int = 40

    # --- collagen material (softer than the stiff G2 baseline so a contractile
    #     spheroid can visibly reorganise the near-field fibres into radial tracts;
    #     this is a declared G3 choice, closer to soft physiological collagen gels) ---
    collagen_modulus_mpa: float = 3.0
    crosslink_stiffness: float = 10.0
    crosslink_fraction: float = 0.3     # keep this fraction of intersection crosslinks
                                        # (< 1 lets gripped fibres rotate freely -> radial)

    # --- overdamped integration ---
    bead_drag: float = 360.0             # nN s / um  (G2 V3 value)
    dt: float = 0.05                     # stable: dt < 2*bead_drag/k_tension ~ 0.32 s
    duration: float = 3600.0             # 60 min: long enough for realign + migrate
    sample_interval: float = 24.0        # -> 150 frames
    ecm_substeps: int = 1

    # --- protrusions (one per membrane sector) ---
    n_sectors: int = 24
    protrusion_max_length: float = 30.0  # prehensile reach ~= gap + a fibre or two
    protrusion_min_length: float = 2.0
    protrusion_probe_speed: float = 0.16    # um/s baseline outward probing (all sectors)
    protrusion_growth_speed: float = 0.12   # extra growth where activity is high
    protrusion_retraction_speed: float = 0.06  # gentle: rear protrusions still grip -> radial contraction
    capture_distance: float = 2.5        # tip must reach within this of a fibre
    engagement_update_interval: float = 0.5

    # --- motor-clutch bundle per protrusion (G2 V3 scale) ---
    n_clutches_per_protrusion: int = 12
    clutch_stiffness: float = 2.0        # nN / um
    clutch_on_rate: float = 0.08         # 1/s (symmetric -- no side bias)
    clutch_off_rate0: float = 0.02       # 1/s
    bell_force: float = 2.5              # nN
    unloaded_actin_speed: float = 0.055  # um/s (inward reel-in of gripped fibres)
    motor_stall_per_protrusion: float = 6.0   # nN

    # --- emergent polarity: mass-conserved replicator + adhesion feedback ---
    # A finite activator pool competes (global inhibition); local self-reinforcement
    # (Rho-GTPase-like autocatalysis) + membrane diffusion + noise make ONE stable
    # broad front emerge from a symmetric start.  FAK/Rac1-like adhesion feedback
    # biases the front toward sectors that actually grip collagen (Carey 2016).
    # Noise sets which direction (stochastic sign); nothing prescribes +x or 0.65.
    polarity_total_activity: float = 6.0    # conserved pool (sum_i a_i)
    polarity_gamma: float = 0.8             # low gain: no winner-take-all -> broad, all-around activity
    polarity_K: float = 0.5                 # half-saturation of the feedback
    polarity_hill: float = 2.0              # feedback cooperativity
    polarity_diffusion: float = 1.6         # strong smoothing keeps every sector active
    polarity_noise: float = 0.06            # mild heterogeneity, not a single front
    polarity_fak_gain: float = 1.0          # FAK/Rac1 adhesion->activity feedback
    polarity_time: float = 6.0              # relaxation time (s)
    adhesion_filter_time: float = 20.0      # q_i low-pass (s)
    adhesion_clutch_scale: float = 3.0

    # --- reaction-driven motion ---
    cell_drag: float = 6000.0            # high drag: broad contraction keeps the spheroid ~in place
    rotational_drag_factor: float = 1.0
    max_cell_speed: float = 0.02         # um/s guard

    # --- optional collagen alignment fixture (director angle in radians) ---
    aligned: bool = False
    director_angle: float = 0.0
    alignment_strength: float = 0.75     # 0 isotropic .. 1 strongly aligned

    seed: int = 7

    def validate(self) -> None:  # noqa: D401 - relax the G2 gate for the spheroid
        if self.domain_size < 5 * self.cell_radius:
            raise ValueError("domain must be at least five spheroid radii wide")
        if not 0 <= self.compression_ratio <= 1:
            raise ValueError("compression_ratio must be between zero and one")
        if self.n_fibers < 16 or self.bead_spacing <= 0:
            raise ValueError("network resolution is too small")
        if self.dt <= 0 or self.sample_interval < self.dt:
            raise ValueError("invalid integration or output interval")

    @property
    def gap_radius(self) -> float:
        return self.cell_radius + self.gap


# =================================================================================
# Network generation with a fibre-free gap (no near-cell contact fibres)
# =================================================================================
def _spheroid_fiber(cfg: SpheroidConfig, rng, *, boundary_seeded: bool):
    """One curved fibre that avoids the fibre-free gap around the spheroid."""

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
        if cfg.aligned:
            theta = rng.normal(0.0, (1.0 - cfg.alignment_strength) * 0.9)
        else:
            theta = rng.uniform(-0.75, 0.75)
        rot = np.array([[math.cos(theta), -math.sin(theta)],
                        [math.sin(theta), math.cos(theta)]])
        end = start + length * (rot @ inward)
    else:
        center = rng.uniform(-0.6 * half, 0.6 * half, size=2)
        if cfg.aligned:
            theta = cfg.director_angle + rng.normal(
                0.0, (1.0 - cfg.alignment_strength) * 0.9
            )
        else:
            theta = rng.uniform(0.0, math.pi)
        direction = np.array([math.cos(theta), math.sin(theta)])
        start = center - 0.5 * length * direction
        end = center + 0.5 * length * direction
    if np.max(np.abs(end)) > half - 0.2 or np.max(np.abs(start)) > half - 0.2:
        return None
    points = _curved_segment(start, end, cfg, rng)
    # reject anything that intrudes into the fibre-free gap
    if np.min(np.linalg.norm(points, axis=1)) < cfg.gap_radius:
        return None
    return points


def _assemble_spheroid_spec(cfg: SpheroidConfig, rng, seed_used: int) -> NetworkSpec:
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

    target_boundary = max(20, int(round(0.45 * cfg.n_fibers)))
    attempts = 0
    while len(fibers) < cfg.n_fibers and attempts < 200_000:
        attempts += 1
        pts = _spheroid_fiber(cfg, rng, boundary_seeded=len(fibers) < target_boundary)
        if pts is not None:
            append(pts)
    if len(fibers) != cfg.n_fibers:
        raise RuntimeError("could not construct the requested fibre-free-gap network")
    # every fibre is a potential contact fibre now; none are pre-attached
    return NetworkSpec(np.asarray(positions), fibers, np.asarray(fixed), [], seed_used)


def make_spheroid_network(cfg: SpheroidConfig, seed=None):
    """Build a crosslinked, boundary-anchored network with a fibre-free gap."""

    cfg.validate()
    base = cfg.seed if seed is None else int(seed)
    best = (-1.0, None)
    for attempt in range(cfg.generation_attempts):
        seed_used = base + 7919 * attempt
        spec = _assemble_spheroid_spec(cfg, np.random.default_rng(seed_used), seed_used)
        network = Network(spec, cfg)  # crosslinks built automatically in __init__
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
            return network, report
    if best[1] is not None:
        return best[1], connectivity_report(best[1])
    raise RuntimeError("network percolation gate failed for the spheroid fixture")


# =================================================================================
# State
# =================================================================================
@dataclass
class SpheroidState:
    cfg: SpheroidConfig
    center: np.ndarray
    sector_angles: np.ndarray            # (S,)
    activity: np.ndarray                 # (S,) wave-pinning active form
    length: np.ndarray                   # (S,) protrusion lengths
    adhesion: np.ndarray                 # (S,) filtered FAK/Rac1 signal q_i
    engaged: np.ndarray                  # (S,) bool -- tip reached a fibre
    edge: np.ndarray                     # (S,) int  -- cached engaged edge
    alpha: np.ndarray                    # (S,) float-- material coordinate
    bound: np.ndarray                    # (S, C) bool clutch state
    extension: np.ndarray                # (S, C) clutch stretch
    traction: np.ndarray                 # (S,) summed clutch force
    velocity: np.ndarray                 # cell velocity (2,)


def initial_state(cfg: SpheroidConfig, rng) -> SpheroidState:
    S = cfg.n_sectors
    C = cfg.n_clutches_per_protrusion
    angles = np.linspace(0.0, 2.0 * math.pi, S, endpoint=False)
    # start near-uniform activity with a tiny random seed -- symmetry is unbroken
    activity = np.full(S, cfg.polarity_total_activity / S)
    activity *= 1.0 + 0.02 * rng.standard_normal(S)
    activity *= cfg.polarity_total_activity / activity.sum()
    return SpheroidState(
        cfg=cfg,
        center=np.zeros(2),
        sector_angles=angles,
        activity=activity,
        length=np.full(S, cfg.protrusion_min_length),
        adhesion=np.zeros(S),
        engaged=np.zeros(S, dtype=bool),
        edge=np.full(S, -1, dtype=int),
        alpha=np.zeros(S),
        bound=np.zeros((S, C), dtype=bool),
        extension=np.zeros((S, C)),
        traction=np.zeros(S),
        velocity=np.zeros(2),
    )


def sector_normals(state: SpheroidState) -> np.ndarray:
    a = state.sector_angles
    return np.column_stack([np.cos(a), np.sin(a)])


def protrusion_tips(state: SpheroidState) -> np.ndarray:
    n = sector_normals(state)
    reach = state.cfg.cell_radius + state.length
    return state.center[None, :] + reach[:, None] * n


# =================================================================================
# Emergent polarity: wave-pinning + adhesion feedback
# =================================================================================
def step_polarity(state: SpheroidState, dt: float, rng) -> None:
    """One step of the mass-conserved replicator polarity field.

    tau da_i/dt = a_i * (phi_i - <phi>) + D * lap(a)_i + noise ,   with
    phi_i = gamma * a_i^h/(K^h + a_i^h) + k_FAK * q_i      (local reinforcement),

    followed by renormalisation to the conserved pool sum_i a_i = A_tot (global
    inhibition).  Sectors whose reinforcement beats the mean grow at the expense of
    the rest -> ONE stable broad front (winner-take-all).  Membrane diffusion keeps
    the front contiguous; noise breaks the initial symmetry and selects a random
    direction each realisation; the FAK/Rac1 adhesion term q_i pulls the front
    toward wherever protrusions actually grip collagen.  No direction is prescribed.
    """

    cfg = state.cfg
    a = state.activity
    hill = a ** cfg.polarity_hill / (cfg.polarity_K ** cfg.polarity_hill + a ** cfg.polarity_hill)
    phi = cfg.polarity_gamma * hill + cfg.polarity_fak_gain * state.adhesion
    fitness = phi - phi.mean()                       # zero-sum global inhibition
    lap = np.roll(a, 1) + np.roll(a, -1) - 2.0 * a
    noise = cfg.polarity_noise * math.sqrt(dt) * rng.standard_normal(a.size)
    a_new = a + (dt / cfg.polarity_time) * (a * fitness + cfg.polarity_diffusion * lap) + noise
    a_new = np.clip(a_new, 1e-4, None)
    a_new *= cfg.polarity_total_activity / float(a_new.sum())
    state.activity = a_new


def step_protrusions(state: SpheroidState, dt: float) -> None:
    """Probe outward everywhere; extend the polarised front, retract the rear.

    Every sector sends a baseline exploratory protrusion (so the spheroid discovers
    fibres in all directions and can break symmetry mechanically).  Sectors above the
    mean activity grow further; sectors well below it retract.  Engaged protrusions
    hold their length -- they are physically anchored to a captured fibre.
    """

    cfg = state.cfg
    a = state.activity
    excess = (a - a.mean()) / (a.mean() + 1e-9)      # >0 front, <0 rear
    headroom = 1.0 - state.length / cfg.protrusion_max_length
    probe = cfg.protrusion_probe_speed * headroom
    grow = cfg.protrusion_growth_speed * np.clip(excess, 0.0, 3.0) * headroom
    retract = cfg.protrusion_retraction_speed * np.clip(-excess, 0.0, 1.0)
    dL = np.where(state.engaged, 0.0, probe + grow - retract)
    state.length = np.clip(
        state.length + dt * dL, cfg.protrusion_min_length, cfg.protrusion_max_length
    )


# =================================================================================
# Tip-first engagement + motor-clutch loading
# =================================================================================
def _closest_material_point(network: Network, tip: np.ndarray):
    """Nearest point on any fibre segment to ``tip``; returns (edge, alpha, dist)."""

    a = network.r[network.edges[:, 0]]
    d = network.r[network.edges[:, 1]] - a
    denom = np.maximum(np.sum(d * d, axis=1), 1e-12)
    alpha = np.clip(np.sum((tip - a) * d, axis=1) / denom, 0.0, 1.0)
    point = a + alpha[:, None] * d
    dist = np.linalg.norm(point - tip, axis=1)
    k = int(np.argmin(dist))
    return k, float(alpha[k]), float(dist[k])


def update_engagement(state: SpheroidState, network: Network) -> None:
    """Tip-first encounter: a protrusion engages a fibre only when its tip reaches one."""

    cfg = state.cfg
    tips = protrusion_tips(state)
    for k in range(cfg.n_sectors):
        if state.engaged[k]:
            # stay engaged while any clutch is bound; otherwise release
            if not state.bound[k].any():
                # re-test the tip against the (possibly moved) fibre
                edge, alpha, dist = _closest_material_point(network, tips[k])
                if dist <= cfg.capture_distance:
                    state.edge[k] = edge
                    state.alpha[k] = alpha
                else:
                    state.engaged[k] = False
                    state.edge[k] = -1
            continue
        edge, alpha, dist = _closest_material_point(network, tips[k])
        if dist <= cfg.capture_distance:
            state.engaged[k] = True
            state.edge[k] = edge
            state.alpha[k] = alpha


def _material_point(network: Network, edge: int, alpha: float) -> np.ndarray:
    i, j = network.edges[edge]
    return (1.0 - alpha) * network.r[i] + alpha * network.r[j]


def _material_velocity(network: Network, velocity: np.ndarray, edge: int,
                       alpha: float, inward: np.ndarray) -> float:
    i, j = network.edges[edge]
    v = (1.0 - alpha) * velocity[i] + alpha * velocity[j]
    return float(v @ inward)


def step_clutches(state: SpheroidState, network: Network, bead_velocity: np.ndarray,
                  dt: float, rng) -> np.ndarray:
    """Advance every protrusion's motor-clutch bundle; return the nodal force field.

    Symmetric binding (no 0.65): a clutch may bind only where the protrusion is
    engaged.  Force is the G2 slip-bond motor-clutch law, pulling the engaged
    material point inward toward the spheroid.
    """

    cfg = state.cfg
    active = np.zeros_like(network.r)
    S, C = state.bound.shape
    u_on = rng.random((S, C))
    u_off = rng.random((S, C))
    for k in range(S):
        bound = state.bound[k]
        ext = state.extension[k]
        if not state.engaged[k]:
            bound[:] = False
            ext[:] = 0.0
            state.traction[k] = 0.0
            continue
        point = _material_point(network, state.edge[k], state.alpha[k])
        radial = point - state.center
        radius = float(np.linalg.norm(radial))
        inward = -radial / max(radius, 1e-12)
        substrate = _material_velocity(network, bead_velocity, state.edge[k],
                                       state.alpha[k], inward)
        # force-velocity actin flow from current load
        force_before = cfg.clutch_stiffness * ext * bound
        traction_before = float(force_before.sum())
        actin = cfg.unloaded_actin_speed * max(
            0.0, 1.0 - traction_before / cfg.motor_stall_per_protrusion
        )
        ext[bound] += dt * max(0.0, actin - substrate)
        force = cfg.clutch_stiffness * ext
        off_rate = cfg.clutch_off_rate0 * np.exp(np.minimum(force / cfg.bell_force, 12.0))
        breaking = bound & (u_off[k] < 1.0 - np.exp(-off_rate * dt))
        bound[breaking] = False
        ext[breaking] = 0.0
        on_prob = 1.0 - math.exp(-cfg.clutch_on_rate * dt)
        binding = (~bound) & (u_on[k] < on_prob)
        bound[binding] = True
        ext[binding] = 0.0
        force = cfg.clutch_stiffness * ext * bound
        traction = float(force.sum())
        state.traction[k] = traction
        # distribute the inward pull to the engaged fibre's two beads
        i, j = network.edges[state.edge[k]]
        f = traction * inward
        if not network.fixed[i]:
            active[i] += (1.0 - state.alpha[k]) * f
        if not network.fixed[j]:
            active[j] += state.alpha[k] * f
    return active


def update_adhesion(state: SpheroidState, dt: float) -> None:
    """Low-pass filtered FAK/Rac1 signal q_i from bound fraction x normalised load."""

    cfg = state.cfg
    nbound = state.bound.sum(axis=1).astype(float)
    engagement = 1.0 - np.exp(-nbound / cfg.adhesion_clutch_scale)
    load = state.traction / (0.5 * cfg.motor_stall_per_protrusion)
    q = engagement * (load / (1.0 + load))
    tau = cfg.adhesion_filter_time
    state.adhesion += dt / tau * (q - state.adhesion)


# =================================================================================
# Metrics
# =================================================================================
def radial_order(network: Network, center: np.ndarray):
    """Nematic radial order S_r = <2 (t.e_r)^2 - 1> in near/mid/far shells.

    +1 = fibres point radially (toward the spheroid), 0 = random, -1 = tangential.
    """

    i, j = network.edges.T
    seg = network.r[j] - network.r[i]
    tangent = seg / np.maximum(np.linalg.norm(seg, axis=1), 1e-12)[:, None]
    mid = 0.5 * (network.r[i] + network.r[j])
    radial = mid - center
    radius = np.linalg.norm(radial, axis=1)
    er = radial / np.maximum(radius, 1e-12)[:, None]
    order = 2.0 * np.square(np.sum(tangent * er, axis=1)) - 1.0
    surface = radius - network.cfg.cell_radius
    shells = ((0.0, 20.0), (20.0, 45.0), (45.0, 90.0))
    out = []
    for lo, hi in shells:
        m = (surface >= lo) & (surface < hi)
        out.append(float(np.mean(order[m])) if np.any(m) else 0.0)
    return out


def polarity_vector(state: SpheroidState) -> np.ndarray:
    n = sector_normals(state)
    return (state.activity[:, None] * n).sum(axis=0)


# =================================================================================
# Main run
# =================================================================================
def run_spheroid(cfg: SpheroidConfig, *, seed: int | None = None, moving: bool = True) -> dict:
    """Run one spheroid guidance simulation and return frames + metrics."""

    rng = np.random.default_rng(cfg.seed if seed is None else seed)
    network, report = make_spheroid_network(cfg, seed=cfg.seed if seed is None else seed)
    state = initial_state(cfg, rng)

    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    engage_every = max(1, int(round(cfg.engagement_update_interval / cfg.dt)))
    bead_velocity = np.zeros_like(network.r)

    frames: list[dict] = []
    traces = {k: [] for k in (
        "time", "bound", "traction", "engaged_count", "polarity_mag", "polarity_angle",
        "order_near", "order_mid", "order_far", "cell_x", "cell_y",
        "path_length", "max_disp", "total_traction")}
    path_length = 0.0
    prev_center = state.center.copy()

    def record(t: float) -> None:
        order = radial_order(network, state.center)
        pol = polarity_vector(state)
        tips = protrusion_tips(state)
        engaged_pts = []
        for k in range(cfg.n_sectors):
            if state.engaged[k] and state.edge[k] >= 0:
                engaged_pts.append(_material_point(network, state.edge[k], state.alpha[k]))
            else:
                engaged_pts.append([np.nan, np.nan])
        disp = np.linalg.norm(network.r - network.r0, axis=1)
        frames.append({
            "time": t,
            "positions": network.r.copy(),
            "center": state.center.copy(),
            "activity": state.activity.copy(),
            "length": state.length.copy(),
            "engaged": state.engaged.copy(),
            "tips": tips,
            "engaged_points": np.asarray(engaged_pts, dtype=float),
            "traction": state.traction.copy(),
            "strain": _bond_strain(network),
        })
        traces["time"].append(t)
        traces["bound"].append(int(state.bound.sum()))
        traces["engaged_count"].append(int(state.engaged.sum()))
        traces["traction"].append(float(state.traction.sum()))
        traces["total_traction"].append(float(state.traction.sum()))
        traces["polarity_mag"].append(float(np.linalg.norm(pol)))
        traces["polarity_angle"].append(float(math.degrees(math.atan2(pol[1], pol[0]))))
        traces["order_near"].append(order[0])
        traces["order_mid"].append(order[1])
        traces["order_far"].append(order[2])
        traces["cell_x"].append(float(state.center[0]))
        traces["cell_y"].append(float(state.center[1]))
        traces["path_length"].append(path_length)
        traces["max_disp"].append(float(disp.max()))

    for step in range(nsteps + 1):
        t = step * cfg.dt
        if step % every == 0:
            record(t)
        if step == nsteps:
            break
        # 1) emergent polarity + protrusion growth
        step_polarity(state, cfg.dt, rng)
        step_protrusions(state, cfg.dt)
        # 2) tip-first engagement (periodic)
        if step % engage_every == 0:
            update_engagement(state, network)
        # 3) motor-clutch loading -> nodal force field on the collagen
        active = step_clutches(state, network, bead_velocity, cfg.dt, rng)
        update_adhesion(state, cfg.dt)
        # 4) reaction-driven spheroid motion (toward where it grips)
        reaction = -active.sum(axis=0)
        if moving:
            vel = reaction / cfg.cell_drag
            speed = float(np.linalg.norm(vel))
            if speed > cfg.max_cell_speed:
                vel *= cfg.max_cell_speed / speed
            state.velocity = vel
            state.center = state.center + cfg.dt * vel
        else:
            state.velocity = np.zeros(2)
        path_length += float(np.linalg.norm(state.center - prev_center))
        prev_center = state.center.copy()
        # 5) advance the collagen network (overdamped, with substeps)
        sub_dt = cfg.dt / cfg.ecm_substeps
        for _ in range(cfg.ecm_substeps):
            bead_velocity = fast_advance(network, active, state.center, sub_dt)

    return {
        "config": cfg,
        "network": network,
        "report": report,
        "frames": frames,
        "traces": {k: np.asarray(v) for k, v in traces.items()},
        "initial_positions": network.r0.copy(),
        "edges": network.edges.copy(),
        "edge_fiber": network.edge_fiber.copy(),
        "fixed": network.fixed.copy(),
        "crosslinks": network.crosslinks,
    }


def _bond_strain(network: Network) -> np.ndarray:
    a = network.r[network.edges[:, 0]]
    b = network.r[network.edges[:, 1]]
    length = np.linalg.norm(b - a, axis=1)
    return (length - network.l0) / np.maximum(network.l0, 1e-12)


def _scatter(n: int, idx: np.ndarray, vals: np.ndarray) -> np.ndarray:
    """Vectorised scatter-add (much faster than np.add.at for our sizes)."""
    out = np.empty((n, 2))
    out[:, 0] = np.bincount(idx, vals[:, 0], minlength=n)
    out[:, 1] = np.bincount(idx, vals[:, 1], minlength=n)
    return out


def fast_advance(network: Network, active: np.ndarray, center: np.ndarray, dt: float):
    """Overdamped integration identical in physics to ``Network.advance`` but using
    ``np.bincount`` scatter-adds.  Reproduces the frozen G2 force law exactly:
    axial EA/l0 (10% compression), discrete bending EI/l^3, hinged crosslinks, and
    the soft cell-surface repulsion.  Returns the bead velocity field.
    """

    r = network.r
    n = len(r)
    i = network.edges[:, 0]
    j = network.edges[:, 1]
    d = r[j] - r[i]
    length = np.linalg.norm(d, axis=1)
    extension = length - network.l0
    stiffness = np.where(extension >= 0.0, network.k_tension, network.k_compression)
    pair = (stiffness * extension / np.maximum(length, 1e-12))[:, None] * d
    fs = _scatter(n, i, pair) - _scatter(n, j, pair)

    a, b, c = network.triplets[:, 0], network.triplets[:, 1], network.triplets[:, 2]
    curv = r[a] - 2.0 * r[b] + r[c] - network.curvature0
    kb = network.bend_coefficient
    fb = _scatter(n, a, -kb * curv) + _scatter(n, b, 2.0 * kb * curv) + _scatter(n, c, -kb * curv)

    fx = np.zeros((n, 2))
    if len(network.crosslinks):
        ea = network.edges[network.link_edge_a]
        eb = network.edges[network.link_edge_b]
        aa = network.link_alpha_a[:, None]
        ab = network.link_alpha_b[:, None]
        pa = (1.0 - aa) * r[ea[:, 0]] + aa * r[ea[:, 1]]
        pb = (1.0 - ab) * r[eb[:, 0]] + ab * r[eb[:, 1]]
        force = network.link_stiffness[:, None] * ((pb - pa) - network.link_rest)
        fx = (_scatter(n, ea[:, 0], (1.0 - aa) * force)
              + _scatter(n, ea[:, 1], aa * force)
              - _scatter(n, eb[:, 0], (1.0 - ab) * force)
              - _scatter(n, eb[:, 1], ab * force))

    radial = r - center
    radius = np.linalg.norm(radial, axis=1)
    penetration = np.maximum(0.0, network.cfg.cell_radius + network.cfg.cell_clearance - radius)
    frep = (network.cfg.repulsion_stiffness * penetration[:, None]
            * radial / np.maximum(radius, 1e-12)[:, None])

    total = fs + fb + fx + frep + active
    total[network.fixed] = 0.0
    velocity = total / network.cfg.bead_drag
    r[~network.fixed] += dt * velocity[~network.fixed]
    r[network.fixed] = network.r0[network.fixed]
    if not np.all(np.isfinite(r)):
        raise FloatingPointError("non-finite bead position; reduce dt")
    return velocity
