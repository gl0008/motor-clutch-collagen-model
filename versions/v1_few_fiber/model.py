"""Version 1: fixed circular cell pulling five independent elastic fibers.

This is the archived, mechanism-first model that produced ``demo/index.html``.
It intentionally has no SLS state and no crosslinks.  Connected beads form a
polyline; the displayed fiber is therefore a line, not a set of popping dots.
"""

from dataclasses import asdict, dataclass
import numpy as np


@dataclass(frozen=True)
class Config:
    n_fibers: int = 5
    beads_per_fiber: int = 25
    bead_spacing: float = 0.75
    cell_radius: float = 2.5
    surface_gap: float = 0.65
    fiber_angles_deg: tuple = (0.0, 34.0, 72.0, 116.0, 154.0)
    offset_signs: tuple = (1, -1, 1, -1, 1)
    stretch_stiffness: float = 18.0
    bend_stiffness: float = 0.32
    bead_drag: float = 6.0
    total_pull_force: float = 12.0
    gaussian_sigma: float = 1.1
    force_ramp_time: float = 0.8
    dt: float = 0.004
    duration: float = 8.0
    sample_interval: float = 0.04


class Network:
    def __init__(self, positions, edges, triplets, fibers, fixed, cfg):
        self.r = np.asarray(positions, float)
        self.r0 = self.r.copy()
        self.edges = np.asarray(edges, int)
        self.triplets = np.asarray(triplets, int)
        self.fibers = np.asarray(fibers, int)
        self.fixed = np.asarray(fixed, bool)
        self.ks, self.kb, self.zeta = cfg.stretch_stiffness, cfg.bend_stiffness, cfg.bead_drag
        self.l0 = np.linalg.norm(self.r[self.edges[:, 1]] - self.r[self.edges[:, 0]], axis=1)

    def forces(self):
        """Elastic force: ``F = -d(U_stretch + U_bend)/dr``."""
        f = np.zeros_like(self.r)
        i, j = self.edges.T
        d = self.r[j] - self.r[i]
        length = np.linalg.norm(d, axis=1)
        pair = (self.ks * (length - self.l0) / np.maximum(length, 1e-12))[:, None] * d
        np.add.at(f, i, pair)
        np.add.at(f, j, -pair)
        a, b, c = self.triplets.T
        curvature = self.r[a] - 2 * self.r[b] + self.r[c]
        np.add.at(f, a, -self.kb * curvature)
        np.add.at(f, b, 2 * self.kb * curvature)
        np.add.at(f, c, -self.kb * curvature)
        f[self.fixed] = 0
        return f

    def advance(self, active, dt):
        """Overdamped update ``zeta * dr/dt = F_elastic + F_active``."""
        velocity = (self.forces() + active) / self.zeta
        velocity[self.fixed] = 0
        self.r += dt * velocity
        self.r[self.fixed] = self.r0[self.fixed]
        return velocity


def make_network(cfg=Config()):
    r, edges, triplets, fibers, fixed = [], [], [], [], []
    s = (np.arange(cfg.beads_per_fiber) - (cfg.beads_per_fiber - 1) / 2) * cfg.bead_spacing
    for fid in range(cfg.n_fibers):
        theta = np.deg2rad(cfg.fiber_angles_deg[fid])
        tangent = np.array([np.cos(theta), np.sin(theta)])
        normal = cfg.offset_signs[fid] * np.array([-tangent[1], tangent[0]])
        closest = (cfg.cell_radius + cfg.surface_gap) * normal
        ids = []
        for local, arc in enumerate(s):
            idx = len(r); ids.append(idx); r.append(closest + arc * tangent)
            fixed.append(local in (0, cfg.beads_per_fiber - 1))
            if local: edges.append((idx - 1, idx))
            if local >= 2: triplets.append((idx - 2, idx - 1, idx))
        fibers.append(ids)
    return Network(r, edges, triplets, fibers, fixed, cfg)


def gaussian_forces(network, cfg, traction):
    """V1 all-bead Gaussian: ``w_i ∝ exp[-d_surface²/sigma²]``.

    V2 replaces this with a hard contact shell followed by a Gaussian only
    inside that shell, so this archived behavior must not be mistaken for V2.
    """
    radius = np.linalg.norm(network.r, axis=1)
    raw = np.exp(-((radius - cfg.cell_radius) / cfg.gaussian_sigma) ** 2)
    raw[network.fixed] = 0
    weights = raw / raw.sum()
    inward = -network.r / np.maximum(radius, 1e-12)[:, None]
    return traction * weights[:, None] * inward, weights


def run(cfg=Config()):
    net = make_network(cfg)
    nsteps = round(cfg.duration / cfg.dt)
    every = round(cfg.sample_interval / cfg.dt)
    out = {"time": [], "positions": [], "weights": [], "config": asdict(cfg)}
    weights = np.zeros(len(net.r))
    for step in range(nsteps + 1):
        t = step * cfg.dt
        if step % every == 0:
            out["time"].append(t); out["positions"].append(net.r.copy()); out["weights"].append(weights.copy())
        if step == nsteps: break
        ramp = min(1.0, t / cfg.force_ramp_time)
        active, weights = gaussian_forces(net, cfg, cfg.total_pull_force * ramp)
        net.advance(active, cfg.dt)
    for key in ("positions", "weights"): out[key] = np.asarray(out[key])
    out["time"] = np.asarray(out["time"])
    out["fibers"] = net.fibers
    return out

