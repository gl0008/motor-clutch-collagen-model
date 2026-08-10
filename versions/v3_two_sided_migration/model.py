"""V3: two-sided stochastic motor-clutch imbalance and cell translation.

The code mirrors the equations in the README.  Units are micrometres (length),
an arbitrary force unit (force), and seconds (time).  The force unit is not yet
calibrated to a specific tumor cell; all dimensional geometry is explicit.
"""

from dataclasses import asdict, dataclass
from typing import Optional
import math
import numpy as np


def cross2(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


@dataclass(frozen=True)
class Config:
    domain_size: float = 100.0
    cell_radius: float = 9.0
    n_fibers: int = 30
    min_fiber_length: float = 24.0
    max_fiber_length: float = 78.0
    bead_spacing: float = 0.75
    effective_fiber_diameter: float = 0.30
    curvature_amplitude: float = 0.55
    contact_width: float = 3.0
    gaussian_sigma: float = 1.5
    protrusion_angle_deg: float = 0.0
    protrusion_half_width_deg: float = 30.0
    stretch_stiffness: float = 1.0
    bend_stiffness: float = 5.6e-4  # kappa/(ks*ds^2) = 1e-3
    crosslink_stiffness: float = 1.0
    bead_drag: float = 1.0
    repulsion_stiffness: float = 8.0
    total_pull_force: float = 5.0
    force_ramp_time: float = 1.0
    dt: float = 0.006
    duration: float = 6.0
    sample_interval: float = 0.06
    seed: int = 7
    n_clutches_per_side: int = 12
    clutch_stiffness: float = 0.8
    clutch_on_rate: float = 1.2
    clutch_off_rate0: float = 0.12
    bell_force: float = 1.0
    unloaded_actin_speed: float = 0.32
    motor_stall_force: float = 12.0
    cell_drag: float = 35.0

    def validate(self):
        if self.cell_radius <= 0 or self.domain_size <= 4 * self.cell_radius:
            raise ValueError("domain must be much larger than the cell")
        if self.contact_width <= 0 or self.gaussian_sigma <= 0:
            raise ValueError("contact scales must be positive")
        if self.dt <= 0 or self.duration <= 0:
            raise ValueError("time scales must be positive")


@dataclass
class Crosslink:
    """A freely hinged material-point link between two fiber segments."""

    edge_a: int
    alpha_a: float
    edge_b: int
    alpha_b: float
    rest_vector: np.ndarray


@dataclass
class ContactPatch:
    fiber: int
    edge: int
    alpha: float
    point: np.ndarray
    surface_distance: float
    weight: float
    normal_in: np.ndarray


class Network:
    def __init__(self, positions, fibers, fixed, cfg):
        self.r = np.asarray(positions, float)
        self.r0 = self.r.copy()
        self.fibers = [np.asarray(ids, int) for ids in fibers]
        self.fixed = np.asarray(fixed, bool)
        self.cfg = cfg
        edges, edge_fiber = [], []
        self.triplets = []
        for fid, ids in enumerate(self.fibers):
            for a, b in zip(ids[:-1], ids[1:]):
                edges.append((a, b)); edge_fiber.append(fid)
            for a, b, c in zip(ids[:-2], ids[1:-1], ids[2:]):
                self.triplets.append((a, b, c))
        self.edges = np.asarray(edges, int)
        self.edge_fiber = np.asarray(edge_fiber, int)
        self.triplets = np.asarray(self.triplets, int)
        self.l0 = np.linalg.norm(self.r[self.edges[:, 1]] - self.r[self.edges[:, 0]], axis=1)
        a, b, c = self.triplets.T
        # The generated gentle waviness is stress-free.  Bending resists a
        # *change* from this reference curvature rather than straightening the
        # network before the cell has applied any force.
        self.curvature0 = self.r[a] - 2 * self.r[b] + self.r[c]
        self.crosslinks = build_crosslinks(self)

    def material_point(self, edge, alpha):
        a, b = self.edges[edge]
        return (1 - alpha) * self.r[a] + alpha * self.r[b]

    def elastic_forces(self, include_crosslinks=True):
        """Return stretching, bending, crosslink forces and their energies.

        ``F = -∂U/∂r`` with Hookean axial bonds, a discrete-curvature bending
        penalty, and freely hinged crosslinks between material points.
        """
        fs = np.zeros_like(self.r); fb = np.zeros_like(self.r); fx = np.zeros_like(self.r)
        i, j = self.edges.T
        d = self.r[j] - self.r[i]
        length = np.linalg.norm(d, axis=1)
        extension = length - self.l0
        pair = (self.cfg.stretch_stiffness * extension / np.maximum(length, 1e-12))[:, None] * d
        np.add.at(fs, i, pair); np.add.at(fs, j, -pair)
        a, b, c = self.triplets.T
        curv = self.r[a] - 2 * self.r[b] + self.r[c] - self.curvature0
        np.add.at(fb, a, -self.cfg.bend_stiffness * curv)
        np.add.at(fb, b, 2 * self.cfg.bend_stiffness * curv)
        np.add.at(fb, c, -self.cfg.bend_stiffness * curv)
        ex = 0.0
        if include_crosslinks:
            for link in self.crosslinks:
                pa = self.material_point(link.edge_a, link.alpha_a)
                pb = self.material_point(link.edge_b, link.alpha_b)
                delta = (pb - pa) - link.rest_vector
                force = self.cfg.crosslink_stiffness * delta
                ea = self.edges[link.edge_a]; eb = self.edges[link.edge_b]
                fx[ea[0]] += (1 - link.alpha_a) * force
                fx[ea[1]] += link.alpha_a * force
                fx[eb[0]] -= (1 - link.alpha_b) * force
                fx[eb[1]] -= link.alpha_b * force
                ex += 0.5 * self.cfg.crosslink_stiffness * float(delta @ delta)
        energy = {
            "stretch": 0.5 * self.cfg.stretch_stiffness * float(extension @ extension),
            "bend": 0.5 * self.cfg.bend_stiffness * float(np.sum(curv * curv)),
            "crosslink": ex,
        }
        return fs, fb, fx, energy

    def repulsion_forces(self, center):
        """Prevent beads crossing the rigid cell: soft radial surface penalty."""
        radial = self.r - center
        radius = np.linalg.norm(radial, axis=1)
        penetration = np.maximum(0.0, self.cfg.cell_radius + 0.15 - radius)
        return self.cfg.repulsion_stiffness * penetration[:, None] * radial / np.maximum(radius, 1e-12)[:, None]

    def advance(self, active, center, dt, include_crosslinks=True):
        """Integrate ``zeta*r_dot = Fs + Fb + Fx + Frep + Fact``."""
        fs, fb, fx, energy = self.elastic_forces(include_crosslinks)
        total = fs + fb + fx + self.repulsion_forces(center) + active
        total[self.fixed] = 0
        velocity = total / self.cfg.bead_drag
        self.r[~self.fixed] += dt * velocity[~self.fixed]
        self.r[self.fixed] = self.r0[self.fixed]
        return velocity, energy


def _fiber_points(center, tangent, length, spacing, amplitude, phase):
    n = max(4, int(round(length / spacing)) + 1)
    s = np.linspace(-length / 2, length / 2, n)
    normal = np.array([-tangent[1], tangent[0]])
    curve = amplitude * (np.cos(2 * np.pi * s / length + phase) - np.cos(phase))
    return center + s[:, None] * tangent + curve[:, None] * normal


def make_network(cfg=Config(), seed: Optional[int] = None):
    """Generate finite, slightly curved fibers around a correctly scaled cell.

    Twelve fibers are deliberately placed in the contact shell (four in the
    right protrusion sector); remaining fibers are random but rejected if they
    cross the cell.  This is an explicit initialization assumption, recorded in
    metadata, so a missing contact is never hidden by pulling a remote fiber.
    """
    cfg.validate(); rng = np.random.default_rng(cfg.seed if seed is None else seed)
    positions, fibers, fixed = [], [], []
    contact_angles = np.deg2rad([-24, -8, 8, 24, 62, 104, 156, 172, 188, 204, 266, 316])
    half_box = cfg.domain_size / 2

    def append(points):
        ids = []
        for local, p in enumerate(points):
            ids.append(len(positions)); positions.append(p)
            # Numerical far-field anchoring: only the two distant fiber ends.
            fixed.append(local in (0, len(points) - 1))
        fibers.append(ids)

    for phi in contact_angles:
        radial = np.array([np.cos(phi), np.sin(phi)])
        tangent_angle = phi + np.pi / 2 + rng.normal(0, 0.18)
        tangent = np.array([np.cos(tangent_angle), np.sin(tangent_angle)])
        gap = rng.uniform(0.45, 2.5)
        length = rng.uniform(42, 66)
        closest = (cfg.cell_radius + gap + cfg.curvature_amplitude) * radial
        pts = _fiber_points(closest, tangent, length, cfg.bead_spacing, cfg.curvature_amplitude, rng.uniform(0, 2*np.pi))
        append(np.clip(pts, -half_box + 1, half_box - 1))

    attempts = 0
    while len(fibers) < cfg.n_fibers and attempts < 10000:
        attempts += 1
        angle = rng.uniform(0, np.pi); tangent = np.array([np.cos(angle), np.sin(angle)])
        length = rng.uniform(cfg.min_fiber_length, cfg.max_fiber_length)
        center = rng.uniform(-0.36 * cfg.domain_size, 0.36 * cfg.domain_size, 2)
        phase = rng.uniform(0, 2*np.pi)
        pts = _fiber_points(center, tangent, length, cfg.bead_spacing, cfg.curvature_amplitude, phase)
        if np.max(np.abs(pts)) > half_box - 1: continue
        if np.min(np.linalg.norm(pts, axis=1)) < cfg.cell_radius + 0.25: continue
        append(pts)
    if len(fibers) != cfg.n_fibers:
        raise RuntimeError("could not generate requested fibers")
    return Network(positions, fibers, fixed, cfg)


def build_crosslinks(network):
    """Create permanent links at geometric intersections of different fibers.

    Segment coordinates, rather than bead IDs, make the topology much less
    sensitive to bead resolution.  Crosslinks preserve coincidence but impose
    no preferred crossing angle, so they are freely hinged.
    """
    links = []
    starts = network.r[network.edges[:, 0]]; vecs = network.r[network.edges[:, 1]] - starts
    for fa in range(len(network.fibers)):
        ia = np.flatnonzero(network.edge_fiber == fa)
        for fb in range(fa + 1, len(network.fibers)):
            ib = np.flatnonzero(network.edge_fiber == fb)
            p = starts[ia][:, None, :]; r = vecs[ia][:, None, :]
            q = starts[ib][None, :, :]; s = vecs[ib][None, :, :]
            den = cross2(r, s); qp = q - p
            safe = np.abs(den) > 1e-10
            ta = np.where(safe, cross2(qp, s) / np.where(safe, den, 1), -1)
            tb = np.where(safe, cross2(qp, r) / np.where(safe, den, 1), -1)
            hits = np.argwhere(safe & (ta > 1e-5) & (ta < 1-1e-5) & (tb > 1e-5) & (tb < 1-1e-5))
            last_points = []
            for u, v in hits:
                point = p[u, 0] + ta[u, v] * r[u, 0]
                if any(np.linalg.norm(point - old) < 1.5 for old in last_points): continue
                last_points.append(point)
                links.append(Crosslink(int(ia[u]), float(ta[u, v]), int(ib[v]), float(tb[u, v]), np.zeros(2)))
    return links


def angle_in_sector(angle, center, half_width):
    return abs((angle - center + np.pi) % (2*np.pi) - np.pi) <= half_width


def contact_patches(network, center, angle_deg=0.0, half_width_deg=30.0):
    """Hard contact selection followed by Gaussian weighting within contact.

    One closest segment is retained per fiber.  For segment ``a-b``,
    ``alpha=clip((rc-a)·(b-a)/|b-a|²,0,1)`` and the material contact point is
    ``p=(1-alpha)a+alpha*b``.  Only ``0 <= |p-rc|-R <= d_contact`` is eligible.
    """
    center = np.asarray(center, float); candidates = []
    target = np.deg2rad(angle_deg); half = np.deg2rad(half_width_deg)
    for fid in range(len(network.fibers)):
        edge_ids = np.flatnonzero(network.edge_fiber == fid)
        best = None
        for edge in edge_ids:
            ia, ib = network.edges[edge]; a, b = network.r[ia], network.r[ib]
            d = b-a; alpha = float(np.clip((center-a) @ d / max(d@d, 1e-12), 0, 1))
            p = a + alpha*d; radial = p-center; radius = np.linalg.norm(radial)
            surface = radius-network.cfg.cell_radius; angle = math.atan2(radial[1], radial[0])
            if 0 <= surface <= network.cfg.contact_width and angle_in_sector(angle, target, half):
                if best is None or surface < best.surface_distance:
                    best = ContactPatch(fid, int(edge), alpha, p, float(surface), 0.0, -radial/max(radius,1e-12))
        if best is not None: candidates.append(best)
    if not candidates: return []
    raw = np.exp(-np.square([c.surface_distance for c in candidates]) / network.cfg.gaussian_sigma**2)
    raw /= raw.sum()
    for patch, weight in zip(candidates, raw): patch.weight = float(weight)
    return candidates


def active_forces(network, center, total_force, angle_deg=0.0, half_width_deg=30.0):
    patches = contact_patches(network, center, angle_deg, half_width_deg)
    forces = np.zeros_like(network.r)
    for p in patches:
        f = total_force * p.weight * p.normal_in
        a, b = network.edges[p.edge]
        forces[a] += (1-p.alpha)*f; forces[b] += p.alpha*f
    forces[network.fixed] = 0
    return forces, patches


def radial_alignment(network, center, near=18.0):
    """Return ``S_r=<2(t·e_r)^2-1>`` near the cell (1 radial, -1 tangential)."""
    a, b = network.edges.T; segment = network.r[b]-network.r[a]
    tangent = segment/np.maximum(np.linalg.norm(segment, axis=1),1e-12)[:,None]
    midpoint = .5*(network.r[a]+network.r[b]); radial = midpoint-center
    radius = np.linalg.norm(radial,axis=1); er = radial/np.maximum(radius,1e-12)[:,None]
    mask = (radius-network.cfg.cell_radius >= 0) & (radius-network.cfg.cell_radius <= near)
    return float(np.mean(2*np.square(np.sum(tangent[mask]*er[mask],axis=1))-1))


def run(cfg=Config(), include_crosslinks=True):
    net = make_network(cfg); center = np.zeros(2)
    nsteps = round(cfg.duration/cfg.dt); every=max(1,round(cfg.sample_interval/cfg.dt))
    out = {k: [] for k in ("time","positions","contact_points","contact_weights","alignment","rms_displacement","stretch_energy","bend_energy","crosslink_energy","drag_power","active_work_rate")}
    patches=[]; velocity=np.zeros_like(net.r); active=np.zeros_like(net.r); energy={"stretch":0,"bend":0,"crosslink":0}
    for step in range(nsteps+1):
        t=step*cfg.dt
        if step%every==0:
            out["time"].append(t); out["positions"].append(net.r.copy())
            out["contact_points"].append([p.point.tolist() for p in patches]); out["contact_weights"].append([p.weight for p in patches])
            out["alignment"].append(radial_alignment(net,center))
            disp=net.r-net.r0; out["rms_displacement"].append(float(np.sqrt(np.mean(np.sum(disp[~net.fixed]**2,axis=1)))))
            out["stretch_energy"].append(energy["stretch"]); out["bend_energy"].append(energy["bend"]); out["crosslink_energy"].append(energy["crosslink"])
            out["drag_power"].append(float(cfg.bead_drag*np.sum(velocity*velocity))); out["active_work_rate"].append(float(np.sum(active*velocity)))
        if step==nsteps: break
        ramp=min(1.0,t/cfg.force_ramp_time) if cfg.force_ramp_time else 1
        active,patches=active_forces(net,center,cfg.total_pull_force*ramp,cfg.protrusion_angle_deg,cfg.protrusion_half_width_deg)
        velocity,energy=net.advance(active,center,cfg.dt,include_crosslinks)
        if not np.all(np.isfinite(net.r)): raise FloatingPointError("unstable integration")
    for key in ("time","positions","alignment","rms_displacement","stretch_energy","bend_energy","crosslink_energy","drag_power","active_work_rate"): out[key]=np.asarray(out[key])
    out.update({"fibers":[x.tolist() for x in net.fibers],"edges":net.edges,"fixed":net.fixed,"crosslinks":[asdict(x) for x in net.crosslinks],"config":asdict(cfg),"include_crosslinks":include_crosslinks})
    return out


# ---------------------------------------------------------------------------
# V3 extension: stochastic two-sided motor clutches and a translating cell.


@dataclass
class Clutch:
    bound: bool = False
    edge: int = -1
    alpha: float = 0.0
    extension: float = 0.0
    last_point: Optional[np.ndarray] = None
    last_center: Optional[np.ndarray] = None


def counter_uniform(seed, step, side, clutch, channel):
    """Addressed random number: identical address → identical uniform draw.

    Fixed and moving conditions deliberately omit condition from the address.
    Thus they receive the same proposed bind/unbind random numbers; different
    mechanics can still change whether a proposal crosses its hazard threshold.
    """
    x = np.uint64(seed + 0x9E3779B97F4A7C15)
    for value in (step, side, clutch, channel):
        x ^= np.uint64(value + 0x9E3779B9) + (x << np.uint64(6)) + (x >> np.uint64(2))
    # splitmix64 finalizer with intentional uint overflow
    with np.errstate(over="ignore"):
        x ^= x >> np.uint64(30); x *= np.uint64(0xBF58476D1CE4E5B9)
        x ^= x >> np.uint64(27); x *= np.uint64(0x94D049BB133111EB)
        x ^= x >> np.uint64(31)
    return float(x >> np.uint64(11)) / float(1 << 53)


def _choose_patch(patches, uniform):
    if not patches: return None
    cdf=np.cumsum([p.weight for p in patches])
    return patches[min(int(np.searchsorted(cdf,uniform,side="right")),len(patches)-1)]


def _clutch_forces(network, center, clutches, cfg, step, last_dt):
    """Update Bell slip bonds and spread their spring forces to material points."""
    active=np.zeros_like(network.r); summaries=[]; traction=[]; velocities=[]
    side_specs=((0,0.0),(1,180.0))
    for side,angle in side_specs:
        patches=contact_patches(network,center,angle,30.0)
        group=clutches[side]
        total_bound_force=sum(cfg.clutch_stiffness*c.extension for c in group if c.bound)
        actin_speed=cfg.unloaded_actin_speed*max(0.0,1-total_bound_force/cfg.motor_stall_force)
        velocities.append(actin_speed)
        for cid,c in enumerate(group):
            if c.bound:
                point=network.material_point(c.edge,c.alpha)
                radial=point-center; inward=-radial/max(np.linalg.norm(radial),1e-12)
                if c.last_point is not None:
                    relative_move=(point-c.last_point)-(center-c.last_center)
                    substrate_speed=float(relative_move@inward/last_dt)
                    c.extension=max(0.0,c.extension+last_dt*(actin_speed-substrate_speed))
                force=cfg.clutch_stiffness*c.extension
                koff=cfg.clutch_off_rate0*math.exp(min(30.0,abs(force)/cfg.bell_force))
                p_off=1-math.exp(-koff*last_dt)
                if counter_uniform(cfg.seed,step,side,cid,1)<p_off:
                    c.bound=False;c.edge=-1;c.extension=0;c.last_point=None;c.last_center=None
                    continue
                f=force*inward; a,b=network.edges[c.edge]
                active[a]+=(1-c.alpha)*f;active[b]+=c.alpha*f
                c.last_point=point.copy();c.last_center=center.copy()
            else:
                p_on=1-math.exp(-cfg.clutch_on_rate*last_dt)
                if patches and counter_uniform(cfg.seed,step,side,cid,0)<p_on:
                    patch=_choose_patch(patches,counter_uniform(cfg.seed,step,side,cid,2))
                    c.bound=True;c.edge=patch.edge;c.alpha=patch.alpha;c.extension=0
                    c.last_point=network.material_point(c.edge,c.alpha).copy();c.last_center=center.copy()
        forces=[cfg.clutch_stiffness*c.extension if c.bound else 0.0 for c in group]
        traction.append(float(sum(forces)))
        summaries.append({"bound":[bool(c.bound) for c in group],"force":forces,"contact_count":len(patches)})
    active[network.fixed]=0
    return active,summaries,np.asarray(traction),np.asarray(velocities)


def run_clutch_condition(cfg=Config(), moving=False):
    """Run fixed or translating cell with the same initial network/random stream."""
    net=make_network(cfg);center=np.zeros(2);cell_velocity=np.zeros(2)
    clutches=[[Clutch() for _ in range(cfg.n_clutches_per_side)] for _ in range(2)]
    nsteps=round(cfg.duration/cfg.dt);every=max(1,round(cfg.sample_interval/cfg.dt))
    keys=("time","positions","cell_center","cell_velocity","bound","clutch_force","traction","actin_velocity","reaction","alignment","rms_displacement")
    out={k:[] for k in keys};summary=[{"bound":[False]*cfg.n_clutches_per_side,"force":[0]*cfg.n_clutches_per_side,"contact_count":0} for _ in range(2)]
    traction=np.zeros(2);actin=np.full(2,cfg.unloaded_actin_speed);reaction=np.zeros(2)
    for step in range(nsteps+1):
        t=step*cfg.dt
        if step%every==0:
            out["time"].append(t);out["positions"].append(net.r.copy());out["cell_center"].append(center.copy());out["cell_velocity"].append(cell_velocity.copy())
            out["bound"].append([x["bound"] for x in summary]);out["clutch_force"].append([x["force"] for x in summary]);out["traction"].append(traction.copy());out["actin_velocity"].append(actin.copy());out["reaction"].append(reaction.copy())
            out["alignment"].append(radial_alignment(net,center));disp=net.r-net.r0;out["rms_displacement"].append(float(np.sqrt(np.mean(np.sum(disp[~net.fixed]**2,axis=1)))))
        if step==nsteps:break
        active,summary,traction,actin=_clutch_forces(net,center,clutches,cfg,step,cfg.dt)
        reaction=-np.sum(active,axis=0)
        if moving:
            cell_velocity=reaction/cfg.cell_drag;center=center+cfg.dt*cell_velocity
        else:
            cell_velocity=np.zeros(2)
        net.advance(active,center,cfg.dt,True)
    for key in keys:out[key]=np.asarray(out[key])
    out.update({"fibers":[x.tolist() for x in net.fibers],"edges":net.edges,"fixed":net.fixed,"crosslinks":[asdict(x) for x in net.crosslinks],"config":asdict(cfg),"moving":moving})
    return out


def run_fixed_moving_pair(cfg=Config()):
    """Synchronized comparison; only the cell force-balance switch differs."""
    return {"fixed":run_clutch_condition(cfg,False),"moving":run_clutch_condition(cfg,True),"config":asdict(cfg)}
