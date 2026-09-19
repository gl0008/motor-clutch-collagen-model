"""Generation 5 -- R0 mechanical-consistency gate (additive; imports G5, edits nothing).

This module is the **G5-R0 gate**: a force-consistent multicellular Stage-D-with-clutch
driver that ports Gloria's accepted single-cell **G4D** loop
(``generations/g4_v4_single_cell_force_alignment/model.py``: ``run_motor_clutch``,
``cell_velocity_from_ecm_forces`` @303, ``_relative_substrate_speeds`` @274,
``_relocate_failed_sites`` @596, ``_advance_clutches`` reset @558-574,
``decompose_contact_forces`` @185) onto the M-cell organoid.  It REUSES the validated
G5 primitives from :mod:`generations.g5_organoid.model` unchanged (imports only; the
old G5 stays byte-for-byte and remains the static/perimeter-wide/binary-leader control).

Sites are flattened ``M cells x n_contact_sectors``; site ``s`` belongs to cell
``s // n_contact_sectors``.

What R0 fixes vs the current ``_run_invasion_with_clutch`` (model.py:1205), per the
advisor "cell-ECM force consistency + spaghetti collagen" directive:

1. **Separated force channels** -- each frame reports per-cell ``F_clutch`` (the applied
   grip-reel reaction), ``F_steric`` (diagnostic; the push this cell exerts on nearby
   beads -- NOT applied to the cell, so steric can move cell and ECM the same way), and
   ``F_cell_cell`` (:func:`model.cell_cell_forces`).  A parallel/perp-to-fibre clutch
   decomposition (port of G4D ``decompose_contact_forces``) is also reported.
2. **Force-pair (Newton 3rd law)** -- per-cell reaction ``= -sum_{s in c} site_force[s] *
   normal_in[s]`` uses the SAME per-site vectors that :func:`model._project_site_forces`
   pushes onto the ECM.  ``force_pair_residual = max_c |reaction_c + sum ...|`` is ~0 by
   construction (mirrors G4D ``force_balance_error``); the gate asserts it.
3. **Relative substrate speed** -- :func:`_relative_substrate_speeds_mc` subtracts each
   site's owning-cell velocity, so rigid cell translation is not misread as clutch
   loading (port of G4D ``_relative_substrate_speeds``; actin speed is a cell-frame term).
4. **Event-driven relocation + state reset** -- a site's gripped fibre changes ONLY after
   that site fully fails (``site_failures`` from the clutch step); on relocation a fresh
   eligible patch is chosen by G4D's tangent-guidance rule and the site's clutch state is
   RESET (``bound=False``, ``extension=0``, ``site_extension=0``) so no loading history is
   carried to the new fibre.  A still-loaded site is never teleported (contrast the old
   periodic ``organoid_clutch_patches`` re-selection at model.py:1279-1282).

Baseline (see :func:`r0_config`): ``clutch_mode="shared"`` (Erdmann-Schwarz 2004 shared
load; G4D V4-A012), ``cell_drag=600`` nN*s/um and ``max_cell_speed=0.012`` um/s (G4D
``g4_interactive_calibration`` defaults), ``strain_stiffening=False`` + ``plasticity=False``
(like-for-like elastic gate; G4D V4-A006/A017 exclude them).  Units um/nN/s.

All own-simulation output here is **personal testing, NOT confirmed findings**
(CLAUDE.md 7.5).  No swirling is claimed.  EMT/leader/switching are OUT of R0 scope
(R1-R3); ``leader_*`` are left at their defaults (no-op).

Citations reused from model.py: Bell 1978; Chan & Odde 2008; Adebowale 2021 SI Table 4
(clutch bundle); G4D / g4_v4_single_cell_force_alignment (force-consistent single cell).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math

import numpy as np

from generations.g5_organoid.model import (  # noqa: E402
    OrganoidConfig,
    parameter_variant,
    make_organoid,
    hex_centers,
    build_crosslinks_grid,
    _inside_any_cell,
    OrganoidStepper,
    organoid_clutch_patches,
    _cell_patches,
    cell_candidate_fibers,
    _project_site_forces,
    _clutch_step,
    _clutch_counter_uniforms,
    _new_clutch_state,
    _clutch_active_mask,
    cell_cell_forces,
    leader_adhesion_scale,
    radial_alignment_profile,
)
from common.model import Network, NetworkSpec, connectivity_report  # noqa: E402
from generations.g4_v2_multiscale.model import (  # noqa: E402  (R2: per-site-stall clutch step)
    bell_off_rate, shared_load_hazard)


# =================================================================================
# R0 baseline config (adopts G4D moving-cell params where G5 diverged)
# =================================================================================
R0_BASELINE = dict(
    clutch_mode="shared",     # Erdmann-Schwarz 2004 shared load; G4D V4-A012
    cell_drag=600.0,          # nN*s/um, G4D g4_interactive_calibration default
    max_cell_speed=0.012,     # um/s, G4D g4_interactive_calibration default
    strain_stiffening=False,  # like-for-like elastic gate (G4D V4-A006 excludes stiffening)
    plasticity=False,         # like-for-like elastic gate (G4D V4-A017 excludes plasticity)
)


def r0_config(base: OrganoidConfig | None = None, **overrides) -> OrganoidConfig:
    """OrganoidConfig with the R0 baseline (G4D-provenanced) applied; overrides win."""
    base = base if base is not None else OrganoidConfig()
    return parameter_variant(base, **{**R0_BASELINE, **overrides})


# =================================================================================
# G4D-style ISOTROPIC-RANDOM network (no corona; radial alignment is an OUTPUT)
# =================================================================================
# model.make_organoid seeds a near-field CORONA (_corona_fiber) with a near-radial
# heading, pre-biasing the initial near field toward radial.  G4D instead fills the
# domain with ISOTROPICALLY oriented finite fibres (angle ~ U[0, pi), midpoint ~ U(square),
# length ~ U[min,max]) and cuts a cell-sized VOID, so radial alignment is a genuinely
# emergent OUTPUT (g4_v3 make_random_void_spec, inherited by G4D).  These reproduce that
# for the M-cell organoid using the UNION-of-cell-disks void.  model.py stays untouched
# (its corona builder remains the legacy control).
def _random_isotropic_spec(cfg, centers: np.ndarray, seed_used: int) -> NetworkSpec:
    """Isotropic finite fibres over the domain, union-of-cells void removed (G4D rule).

    Angle ~ U[0, pi) (never rotated by clipping), midpoint ~ U(square), length ~ U[min,max].
    Each fibre is discretised at ``bead_spacing``; beads inside ANY cell disk
    (``cell_radius + cell_clearance``) are dropped and the contiguous OUTSIDE fragments
    kept, so the finite sample is isotropic (t=0 radial order ~ 0).  Outer-square beads are
    anchored; fibres touching a cell surface (within ``contact_width``) are the grippable
    contact set.
    """
    rng = np.random.default_rng(seed_used)
    half = 0.5 * cfg.domain_size
    void = cfg.cell_radius + cfg.cell_clearance
    positions: list = []
    fibers: list = []
    fixed: list = []
    contact: set = set()
    n = 0
    attempts = 0
    min_frag = max(3, int(round(cfg.min_fiber_length / (3.0 * cfg.bead_spacing))))
    max_attempts = 500 * cfg.n_fibers
    while n < cfg.n_fibers and attempts < max_attempts:
        attempts += 1
        theta = rng.uniform(0.0, math.pi)
        d = np.array([math.cos(theta), math.sin(theta)])
        length = rng.uniform(cfg.min_fiber_length, cfg.max_fiber_length)
        mid = rng.uniform(-0.96 * half, 0.96 * half, size=2)
        a = np.clip(mid - 0.5 * length * d, -half + 0.2, half - 0.2)
        b = np.clip(mid + 0.5 * length * d, -half + 0.2, half - 0.2)  # clip keeps direction
        seg = b - a
        L = float(np.linalg.norm(seg))
        if L < cfg.min_fiber_length * 0.4:
            continue
        nb = max(4, int(math.ceil(L / cfg.bead_spacing)) + 1)
        pts = a[None, :] + np.linspace(0.0, 1.0, nb)[:, None] * seg[None, :]
        dmin = np.sqrt(np.min(np.sum((pts[:, None, :] - centers[None, :, :]) ** 2, axis=2), axis=1))
        outside = dmin >= void
        if not outside.any():
            continue
        idx = np.flatnonzero(outside)
        accepted_any = False
        for run in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
            if len(run) < min_frag:
                continue
            frag = pts[run]
            fid = len(fibers)
            fibers.append(list(range(len(positions), len(positions) + len(frag))))
            touches = False
            for p in frag:
                positions.append(np.asarray(p, dtype=float))
                fixed.append(bool(np.max(np.abs(p)) >= half - cfg.boundary_width))
                if abs(float(np.min(np.linalg.norm(p - centers, axis=1))) - void) < cfg.contact_width:
                    touches = True
            if touches:
                contact.add(fid)
            accepted_any = True
        if accepted_any:
            n += 1
    if n < cfg.n_fibers:
        raise RuntimeError("could not build the requested isotropic organoid network")
    return NetworkSpec(np.asarray(positions), fibers, np.asarray(fixed),
                       sorted(contact), seed_used)


def make_random_organoid(cfg, seed=None):
    """Pack cells + build a G4D-style ISOTROPIC random network (no corona).

    Drop-in for :func:`model.make_organoid` used by the R0/R1 driver so the initial fibre
    alignment is genuinely random (t=0 radial order ~ 0) and radial reorganisation is an
    emergent OUTPUT.  Returns ``(network, centers, gap_radius, report)``.
    """
    cfg.validate()
    centers = hex_centers(cfg.organoid_radius, cfg.cell_spacing)
    organoid_outer = float(np.max(np.linalg.norm(centers, axis=1))) + cfg.cell_radius
    gap_radius = organoid_outer + cfg.gap
    base = cfg.seed if seed is None else int(seed)
    best = (-1, -1.0, None)
    for attempt in range(cfg.generation_attempts):
        seed_used = base + 7919 * attempt
        spec = _random_isotropic_spec(cfg, centers, seed_used)
        network = Network(spec, cfg, crosslinks=[])
        network.crosslinks = build_crosslinks_grid(network)
        network.refresh_crosslink_arrays()
        if cfg.crosslink_fraction < 1.0 and network.crosslinks:
            rng = np.random.default_rng(seed_used + 101)
            keep = rng.random(len(network.crosslinks)) < cfg.crosslink_fraction
            network.crosslinks = [x for x, k in zip(network.crosslinks, keep) if k]
            network.refresh_crosslink_arrays()
        report = connectivity_report(network)
        score = float(report["connected_fraction"])
        cc = 1 if report["contact_fibers_connected"] else 0
        if (cc, score) > (best[0], best[1]):
            best = (cc, score, network)
        if report["contact_fibers_connected"] and score >= cfg.required_connected_fraction:
            return network, centers, gap_radius, report
    if best[2] is not None:
        return best[2], centers, gap_radius, connectivity_report(best[2])
    raise RuntimeError("isotropic organoid network percolation gate failed")


# =================================================================================
# Force-consistency primitives (ports of G4D helpers, generalised to M cells)
# =================================================================================
def per_cell_clutch_reaction(patches: list, site_force: np.ndarray, n_sec: int,
                             n_cells: int) -> np.ndarray:
    """Per-cell grip-reel reaction ``= -sum_{s in cell} site_force[s]*normal_in[s]``.

    Uses the SAME per-site inward vectors that :func:`model._project_site_forces`
    applies to the ECM, so the cell receives the exact equal-and-opposite reaction
    (Newton 3rd law).  Points OUTWARD (invasion).  (G4D ``cell_velocity_from_ecm_forces``.)
    """
    reaction = np.zeros((n_cells, 2))
    for s, patch in enumerate(patches):
        if patch is None or site_force[s] <= 0.0:
            continue
        reaction[s // n_sec] -= site_force[s] * patch.normal_in
    return reaction


def _contact_tangent(network, patch) -> np.ndarray:
    """Unit tangent of the gripped fibre segment (G4D ``_contact_tangent``)."""
    i, j = network.edges[patch.edge]
    t = network.r[j] - network.r[i]
    return t / max(float(np.linalg.norm(t)), 1e-12)


def _relative_substrate_speeds_mc(network, bead_velocity: np.ndarray, patches: list,
                                  site_centers: np.ndarray,
                                  cell_velocity_per_site: np.ndarray) -> np.ndarray:
    """Inward ECM material-point speed RELATIVE to the owning cell, per site.

    ``value = (v_material - v_cell) . inward``.  Subtracting the owning cell's velocity
    (G4D ``_relative_substrate_speeds`` @274) is what prevents rigid cell translation
    from being misread as clutch extension.  Empty sites report 0.
    """
    out = np.zeros(len(patches))
    for s, patch in enumerate(patches):
        if patch is None:
            continue
        point = network.material_point(patch.edge, patch.alpha)
        inward = site_centers[s] - point
        inward /= max(float(np.linalg.norm(inward)), 1e-12)
        i, j = network.edges[patch.edge]
        mv = (1.0 - patch.alpha) * bead_velocity[i] + patch.alpha * bead_velocity[j]
        out[s] = float((mv - cell_velocity_per_site[s]) @ inward)
    return out


def _reset_site_state(state, s: int) -> None:
    """Clear a site's clutch history so no loading is carried to a new fibre.

    (G4D resets extension on rupture @560 and site_extension on full failure @563/574;
    on relocation the whole site is cleared.)
    """
    state.bound[s, :] = False
    state.extension[s, :] = 0.0
    state.site_extension[s] = 0.0


def _relocate_site(network, center: np.ndarray, patches: list, s: int, cell: int,
                   n_sec: int, candidate_fibers, cfg: OrganoidConfig):
    """Choose a fresh eligible patch for a fully-failed site by tangent guidance.

    Port of G4D ``_relocate_failed_sites`` @596 (one site): outgoing direction is the
    old fibre's tangent oriented away from the cell; the next patch is the eligible
    surface material point whose outward direction is nearest that guidance.  Uses NO
    global polarization vector.  Returns a ``ContactPatch`` or ``None`` (keep old).
    """
    fresh = _cell_patches(network, center, candidate_fibers)
    if not fresh:
        return None
    old = patches[s]
    if old is None:
        # site had no fibre; just acquire the closest-surface eligible patch
        return min(fresh, key=lambda p: p.surface_distance)
    old_point = network.material_point(old.edge, old.alpha)
    old_outward = old_point - center
    norm = float(np.linalg.norm(old_outward))
    if norm < 1e-12:
        return min(fresh, key=lambda p: p.surface_distance)
    old_outward /= norm
    guidance = _contact_tangent(network, old)
    if float(guidance @ old_outward) < 0.0:
        guidance = -guidance
    # avoid doubling up on a fibre another live site of this cell already grips
    lo, hi = cell * n_sec, (cell + 1) * n_sec
    used = {patches[k].fiber for k in range(lo, hi) if k != s and patches[k] is not None}
    available = [p for p in fresh if p.fiber not in used] or fresh

    def outward_of(p):
        d = network.material_point(p.edge, p.alpha) - center
        return d / max(float(np.linalg.norm(d)), 1e-12)

    return min(available, key=lambda p: math.acos(
        np.clip(float(guidance @ outward_of(p)), -1.0, 1.0)))


# =================================================================================
# Diagnostics (separated force channels + spaghetti readout)
# =================================================================================
def _per_cell_steric(network, centers: np.ndarray, cfg: OrganoidConfig) -> np.ndarray:
    """DIAGNOSTIC per-cell steric push on nearby beads (net), NOT applied to the cell.

    Mirrors :func:`model.multi_cell_repulsion` but resolved per cell so the steric
    channel can be reported separately from the clutch reaction.  Reported to show
    steric is one-directional (ECM only) whereas clutch is a strict reaction pair.
    """
    reach = cfg.cell_radius + cfg.cell_clearance
    out = np.zeros((len(centers), 2))
    for c, center in enumerate(centers):
        radial = network.r - center
        radius = np.linalg.norm(radial, axis=1)
        pen = np.maximum(0.0, reach - radius)
        act = pen > 0.0
        if np.any(act):
            f = (cfg.repulsion_stiffness * pen[act, None] * radial[act]
                 / np.maximum(radius[act], 1e-12)[:, None])
            out[c] = f.sum(axis=0)
    return out


def _clutch_tangent_decomposition(network, patches: list, site_force: np.ndarray):
    """Aggregate clutch force resolved parallel/perp to the local fibre (G4D @185)."""
    par = perp = tot = 0.0
    for s, patch in enumerate(patches):
        if patch is None or site_force[s] <= 0.0:
            continue
        vec = site_force[s] * patch.normal_in
        t = _contact_tangent(network, patch)
        p = float(vec @ t) * t
        par += float(np.linalg.norm(p))
        perp += float(np.linalg.norm(vec - p))
        tot += float(np.linalg.norm(vec))
    return par, perp, tot


def spaghetti_index(network) -> float:
    """Bulk-displacement readout: ``rms(|r - r0|) / bead_spacing`` (modelling choice).

    NOTE: this measures how far beads MOVED, which conflates *coherent* translation
    (cells reeling collagen inward -> large but ordered displacement) with genuine
    floppiness.  Over a long run it can exceed 1 while fibres stay barely strained, so
    it is a POOR spaghetti test on its own -- read it together with
    :func:`floppiness_index`, which is the honest fibre-deformation measure.
    """
    rms = float(np.sqrt(np.mean(np.sum((network.r - network.r0) ** 2, axis=1))))
    return rms / max(float(network.cfg.bead_spacing), 1e-12)


def floppiness_index(network) -> float:
    """Honest "spaghetti" test: fraction of fibre SEGMENTS strained beyond 50%.

    A coherent (strain-transmitting) network keeps every segment near its rest length;
    a "spaghetti" network has segments grossly stretched/buckled.  ``0`` = fully
    coherent.  Uses the network's rest lengths ``l0`` directly, so pure rigid
    translation of collagen (no deformation) contributes nothing (unlike
    :func:`spaghetti_index`).
    """
    i, j = network.edges.T
    length = np.linalg.norm(network.r[j] - network.r[i], axis=1)
    strain = (length - network.l0) / np.maximum(network.l0, 1e-9)
    return float(np.mean(np.abs(strain) > 0.5))


# =================================================================================
# R0 driver
# =================================================================================
def run_r0_invasion(cfg: OrganoidConfig = None, seed=None, snapshots: bool = False,
                    adhesion_scale=None, pair_rule: str = "min",
                    network_mode: str = "random", site_stall=None) -> dict:
    """Force-consistent multicellular Stage-D-with-clutch invasion (the G5-R0 gate).

    Ports G4D's ``run_motor_clutch`` loop ORDER to M cells: sample frame -> clutch step
    -> project active ECM force -> per-cell force-pair reaction + capped overdamped cell
    motion -> ECM step -> relative substrate speed -> event-driven relocation (+reset) of
    only fully-failed sites.  See module docstring.  ``leader_*`` are out of scope (R0).

    ``adhesion_scale`` (R1 hook): ``None`` -> R0 default (:func:`model.leader_adhesion_scale`,
    all-ones at defaults) so the R0 path is byte-for-byte unchanged; otherwise a per-cell
    adhesion multiplier array OR a ``callable(centers0, cfg) -> array`` (R1 passes
    :func:`emt_phenotype`).  ``pair_rule`` sets how a pair's adhesion combines the two
    cells' scales: ``"min"`` (weakest-member; :func:`model.cell_cell_forces`, R0 default) or
    ``"geomean"`` (sqrt product; sensitivity control, :func:`_cell_cell_forces_geomean`).
    NEITHER touches the clutch/traction -- EMT changes ONLY cell-cell adhesion.
    """
    if cfg is None:
        cfg = r0_config()

    # network_mode "random" = G4D-style ISOTROPIC network (t=0 radial order ~0, the honest
    # default; radial alignment is an OUTPUT); "cued" = isotropic + a one-sided imposed radial
    # tract (R2, opt-in, explicitly NOT swirling); "corona" = legacy model.make_organoid (biased).
    _builder = {"random": make_random_organoid, "cued": make_cued_organoid}.get(
        network_mode, make_organoid)
    network, centers, gap_radius, report = _builder(cfg, seed=seed)
    centers = centers.copy()
    centers0 = centers.copy()
    M = len(centers)
    n_sec = cfg.n_contact_sectors
    organoid_center = np.zeros(2)
    stepper = OrganoidStepper(network, centers)
    reach = cfg.cell_radius + cfg.contact_width + 2.0
    if adhesion_scale is None:
        adh_scale = leader_adhesion_scale(centers0, cfg)     # R0 default -> all ones
    elif callable(adhesion_scale):
        adh_scale = np.asarray(adhesion_scale(centers0, cfg), dtype=float)
    else:
        adh_scale = np.asarray(adhesion_scale, dtype=float)
    # pair rule -> cell-cell force (default "min" reuses model.cell_cell_forces exactly)
    cc_force = ((lambda ctr: _cell_cell_forces_geomean(ctr, cfg, adh_scale))
                if pair_rule == "geomean"
                else (lambda ctr: cell_cell_forces(ctr, cfg, adh_scale)))
    # site_stall (R2): None -> global F_stall (R0/R1); array or callable(centers0,cfg) ->
    # per-site F_stall (leaders' higher motor capacity).  Resolve the callable now.
    if callable(site_stall):
        site_stall = np.asarray(site_stall(centers0, cfg), dtype=float)
    candidates = cell_candidate_fibers(network, centers, reach)
    patches, _site_centers = organoid_clutch_patches(network, centers, cfg, candidates)  # INITIAL only
    S = len(patches)
    active_mask = _clutch_active_mask(patches)
    state = _new_clutch_state(S, cfg)
    substrate = np.zeros(S)
    site_force = np.zeros(S)
    v_cell = np.zeros((M, 2))
    reaction = np.zeros((M, 2))
    n_relocations = 0

    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    contact_every = max(1, int(round(cfg.contact_update_interval / cfg.dt)))
    frames: list[dict] = []
    bead_snaps: list[np.ndarray] = []
    cell_snaps: list[np.ndarray] = []
    max_residual = 0.0

    def sample_frame(time: float) -> dict:
        prof = radial_alignment_profile(network, organoid_center)
        f_steric = _per_cell_steric(network, centers, cfg)
        f_cc = cc_force(centers)
        par, perp, tot = _clutch_tangent_decomposition(network, patches, site_force)
        clutch_ecm = np.zeros((M, 2))              # independent recompute of the ECM side
        for s, patch in enumerate(patches):
            if patch is None or site_force[s] <= 0.0:
                continue
            clutch_ecm[s // n_sec] += site_force[s] * patch.normal_in
        residual = float(np.max(np.linalg.norm(reaction + clutch_ecm, axis=1))) if M else 0.0
        return {
            "time": time,
            "global_radial_order": prof["global_radial_order"],
            "shells": prof["shells"],
            "mean_cell_radial_disp": float(np.mean(
                np.linalg.norm(centers, axis=1) - np.linalg.norm(centers0, axis=1))),
            "cell_spread": float(np.mean(np.linalg.norm(centers, axis=1))),
            "max_cell_disp": float(np.max(np.linalg.norm(centers - centers0, axis=1))),
            # cohesion / invasion-mode metrics (R1; Ilina & Friedl 2020):
            "radius_of_gyration": radius_of_gyration(centers),
            "detached_fraction": detached_fraction(centers, cfg),
            "lcc_fraction": largest_connected_component_fraction(centers, cfg),
            # separated force channels (nN):
            "F_clutch_mean": float(np.mean(np.linalg.norm(reaction, axis=1))),
            "F_clutch_max": float(np.max(np.linalg.norm(reaction, axis=1))) if M else 0.0,
            "F_steric_mean": float(np.mean(np.linalg.norm(f_steric, axis=1))),
            "F_cell_cell_mean": float(np.mean(np.linalg.norm(f_cc, axis=1))),
            "clutch_parallel": par, "clutch_perp": perp, "clutch_total": tot,
            "force_pair_residual": residual,
            "spaghetti_index": spaghetti_index(network),          # bulk displacement (see caveat)
            "floppiness_index": floppiness_index(network),        # honest: frac segments >50% strain
            "n_active_sites": int(active_mask.sum()),
            "bound_fraction": float(state.bound[active_mask].mean()) if active_mask.any() else 0.0,
            "total_traction": float(site_force.sum()),
            "cumulative_site_failures": int(state.cumulative_site_failures),
            "n_relocations": int(n_relocations),
        }

    for step in range(nsteps + 1):
        time = step * cfg.dt
        if step % every == 0:
            fr = sample_frame(time)
            max_residual = max(max_residual, fr["force_pair_residual"])
            frames.append(fr)
            if snapshots:
                bead_snaps.append(network.r.copy())
                cell_snaps.append(centers.copy())
        if step == nsteps:
            break

        # 1) clutch bundles load / Bell-slip / rebind (reused g4_v2 law via model._clutch_step).
        # site_stall (R2) = per-site F_stall so leaders' higher traction EMERGES from the
        # force-velocity law v=v0(1-F/F_stall); None -> global cfg.motor_stall_per_site (R0/R1).
        if site_stall is None:
            _, site_force, breaks, binds, site_fail = _clutch_step(
                cfg, state, substrate, step, active_mask)
        else:
            _, site_force, breaks, binds, site_fail = _clutch_step_stall(
                cfg, state, substrate, step, active_mask, site_stall)

        # 2) project the emergent clutch traction onto the ECM
        active = _project_site_forces(network, patches, site_force)

        # 3) force-pair: cell receives -(the traction it applied); move overdamped + capped
        reaction = per_cell_clutch_reaction(patches, site_force, n_sec, M)
        f_cell = reaction + cc_force(centers)
        v_cell = f_cell / cfg.cell_drag
        speed = np.linalg.norm(v_cell, axis=1)
        over = speed > cfg.max_cell_speed
        if np.any(over):
            v_cell[over] *= (cfg.max_cell_speed / speed[over])[:, None]
        centers += cfg.dt * v_cell
        stepper.centers[:] = centers

        # 4) ECM step (steric pushes beads only; clutch active force), then relative substrate
        stepper.step(active, cfg.dt)
        site_centers = np.repeat(centers, n_sec, axis=0)
        v_cell_per_site = np.repeat(v_cell, n_sec, axis=0)
        substrate = _relative_substrate_speeds_mc(
            network, stepper.velocity, patches, site_centers, v_cell_per_site)

        # cheap eligibility refresh only (does NOT touch patches/state)
        if step and step % contact_every == 0:
            candidates = cell_candidate_fibers(network, centers, reach)

        # 5) event-driven relocation + reset of ONLY fully-failed sites
        if np.any(site_fail):
            for s in np.flatnonzero(site_fail):
                c = int(s // n_sec)
                newp = _relocate_site(network, centers[c], patches, int(s), c,
                                      n_sec, candidates[c], cfg)
                if newp is not None:
                    patches[int(s)] = newp
                _reset_site_state(state, int(s))
                n_relocations += 1
            active_mask = _clutch_active_mask(patches)

    return {
        "config": asdict(cfg),
        "clutch_per_cell_final": np.linalg.norm(reaction, axis=1),  # per-cell |reaction| (R2 leader vs follower)
        "centers0": centers0,
        "centers_final": centers,
        "organoid_center": organoid_center,
        "connectivity": report,
        "n_cells": M,
        "n_beads": len(network.r),
        "n_fibers": len(network.fibers),
        "n_clutch_sites": S,
        "clutch_mode": cfg.clutch_mode,
        "max_force_pair_residual": float(max_residual),
        "cumulative_slips": int(state.cumulative_slips),
        "cumulative_site_failures": int(state.cumulative_site_failures),
        "n_relocations": int(n_relocations),
        "frames": frames,
        "bead_snapshots": np.asarray(bead_snaps) if snapshots else None,
        "cell_snapshots": np.asarray(cell_snaps) if snapshots else None,
        "final_positions": network.r.copy(),
        "initial_positions": network.r0.copy(),
        "edges": network.edges.copy(),
    }


# =================================================================================
# R1 -- discrete partial-EMT phenotype (changes ONLY cell-cell adhesion) + cohesion
#       metrics.  Decouples EMT-PHENOTYPE (this) from LEADER-ROLE (R2).  Built on the
#       force-consistent R0 driver; clutch / drag / motor budget are untouched.
#       Cell-cell adhesion sets collective <-> single-cell (Ilina & Friedl 2020).
# =================================================================================
@dataclass(frozen=True)
class G5RevisionConfig(OrganoidConfig):
    """OrganoidConfig + G5-revision (R1..R3) fields.  Defaults reduce to homogeneous /
    R0-equivalent behaviour (``emt_fraction=0`` -> every cell epithelial -> no-op)."""

    # --- R1: discrete partial-EMT composition (phenotype, NOT leader role) ---
    emt_fraction: float = 0.0          # f_pEMT: fraction of cells that are partial-EMT
    emt_adhesion_factor: float = 0.5   # a_pEMT: partial-EMT adhesion multiplier (1 = epithelial)
    emt_assignment: str = "random"     # "random" (matched, seeded) | "boundary" (outer-enriched control)
    emt_seed: int = 12345
    emt_pair_rule: str = "min"         # "min" (weakest-member) | "geomean" (sensitivity control)
    emt_traction_factor: float = 1.0   # partial-EMT cells' per-site motor-stall multiplier
                                       #   (biology audit B: EMT raises traction/protrusion too, not
                                       #   adhesion alone). 1.0 = R1 adhesion-only (decoupled); >1 =
                                       #   COUPLED EMT-leader (same cells loosen junctions AND pull harder)

    # --- R2: static leader NUMBER + LOCATION (functional role, decoupled from EMT) ---
    n_leaders: int = 0                 # number of ACTIVE leader cells (0 -> R1 baseline, no boost)
    leader_location: str = "cue_front" # "cue_front" (localized adjacent front on the cue side) |
                                       #   "perimeter" (outermost n_leaders spread around; control)
    leader_stall_factor: float = 1.0   # leaders' per-site F_stall multiple (>1 = higher motor
                                       #   capacity -> higher EMERGENT traction; Reffay 2014)
    budget_mode: str = "matched_total" # "matched_total" (total leader stall budget fixed vs N_L:
                                       #   isolates number/geometry) | "fixed_per_leader" (each leader
                                       #   keeps leader_stall_factor; total grows with N_L)
    # --- R2: opt-in one-sided imposed radial ECM cue (EXPLICIT assumption, NOT swirling) ---
    radial_cue: bool = False           # off -> isotropic random (default); on -> one-sided tract
    cue_angle: float = 0.0             # rad, azimuth of the cue front
    cue_half_width: float = 0.6        # rad, angular half-width of the cue sector
    cue_band: float = 40.0             # um, radial thickness of the cued near-field band
    cue_angular_sd: float = 0.35       # rad, jitter of the imposed radial tract (v4E)
    cue_seed: int = 404

    # --- R3: energy-based DYNAMIC leader switching (persistence layer; Zhang 2019) ---
    # off -> static R2 leaders.  On: a front leader that keeps doing reel work DRAINS a
    # per-cell energy; below energy_off it DEMOTES and a fresher front-follower (E>energy_on)
    # is PROMOTED (relay handoff).  n_leaders = MAX simultaneous active leaders.
    leader_switching: bool = False
    energy_E0: float = 1.0             # resting per-cell energy (dimensionless)
    energy_tau_rec: float = 300.0      # s, recovery time toward E0 when idle
    energy_cap: float = 130.0          # nN*um of motor work per unit energy drained (drain rate).
                                       #   Calibrated so at the typical R3 regime (cued, low adhesion,
                                       #   ~x6 leader; measured P_leader~0.34, P_follower~0.17 nN*um/s)
                                       #   the steady energy E*=E0-P*tau_rec/cap puts an ACTIVE leader
                                       #   below energy_off (~0.25) but followers above energy_on (~0.63)
                                       #   -> leader lifetime ~minutes (Zhang 2019 ORDER, compressed for
                                       #   demonstration; calibrated to LIFETIME, not measured ATP).
    energy_on: float = 0.5             # promote a candidate only if E > energy_on ...
    energy_off: float = 0.3            # ... demote an active leader when E < energy_off (hysteresis)


def r1_config(base: OrganoidConfig | None = None, **overrides) -> "G5RevisionConfig":
    """G5RevisionConfig with the R0 baseline (shared clutch, G4D drag/speed, elastic)
    applied, so R1 experiments inherit the validated R0 mechanics; overrides win."""
    base = base if base is not None else G5RevisionConfig()
    return parameter_variant(base, **{**R0_BASELINE, **overrides})


def emt_phenotype(centers: np.ndarray, cfg) -> np.ndarray:
    """Per-cell adhesion multiplier from the EMT composition (R1).

    ``1.0`` for epithelial cells, ``emt_adhesion_factor`` for the ``emt_fraction`` subset.
    ``emt_assignment="random"`` picks ``round(f*M)`` cells deterministically (``emt_seed``);
    ``"boundary"`` picks the outermost ``round(f*M)`` (an enrichment control).  Returns
    all-ones when ``emt_fraction==0`` (homogeneous epithelial baseline).  This is the ONLY
    channel EMT touches -- it feeds ``run_r0_invasion(adhesion_scale=...)``, i.e. cell-cell
    adhesion; the clutch/traction is phenotype-independent.
    """
    centers = np.asarray(centers, dtype=float)
    M = len(centers)
    scale = np.ones(M)
    n = int(round(float(getattr(cfg, "emt_fraction", 0.0)) * M))
    if n <= 0:
        return scale
    if getattr(cfg, "emt_assignment", "random") == "boundary":
        idx = np.argsort(-np.linalg.norm(centers, axis=1))[:n]
    else:
        rng = np.random.default_rng(int(getattr(cfg, "emt_seed", 12345)))
        idx = rng.choice(M, size=n, replace=False)
    scale[idx] = float(getattr(cfg, "emt_adhesion_factor", 0.5))
    return scale


def emt_cell_ids(centers: np.ndarray, cfg) -> np.ndarray:
    """Indices of the partial-EMT cells -- the SAME selection as :func:`emt_phenotype`
    (so the low-adhesion cells and the coupled high-traction cells coincide)."""
    centers = np.asarray(centers, dtype=float)
    M = len(centers)
    n = int(round(float(getattr(cfg, "emt_fraction", 0.0)) * M))
    if n <= 0:
        return np.zeros(0, dtype=int)
    if getattr(cfg, "emt_assignment", "random") == "boundary":
        return np.argsort(-np.linalg.norm(centers, axis=1))[:n]
    rng = np.random.default_rng(int(getattr(cfg, "emt_seed", 12345)))
    return np.asarray(sorted(rng.choice(M, size=n, replace=False)), dtype=int)


def emt_site_stall(centers: np.ndarray, cfg) -> np.ndarray:
    """Per-site F_stall (length M*n_sec): partial-EMT cells get ``emt_traction_factor`` x base.

    Biology audit (B): a real partial-EMT cell raises traction / protrusive matrix-pulling, not
    cell-cell adhesion alone; and leaders often ARE partial-EMT cells (biorxiv 2025).  So the
    COUPLED variant boosts the SAME cells that :func:`emt_phenotype` weakens.  ``emt_traction_factor
    == 1.0`` -> all base -> decoupled (identical to the R1 adhesion-only path)."""
    n_sec = cfg.n_contact_sectors
    base = float(cfg.motor_stall_per_site)
    M = len(centers)
    stall = np.full(M * n_sec, base)
    f = float(getattr(cfg, "emt_traction_factor", 1.0))
    if f != 1.0:
        for c in emt_cell_ids(centers, cfg):
            stall[int(c) * n_sec:(int(c) + 1) * n_sec] = base * f
    return stall


def _cell_cell_forces_geomean(centers: np.ndarray, cfg: OrganoidConfig,
                              adhesion_scale: np.ndarray) -> np.ndarray:
    """Cell-cell force with a GEOMETRIC-MEAN adhesion pair rule (sensitivity control).

    Identical to :func:`model.cell_cell_forces` except a pair's adhesion is scaled by
    ``sqrt(s_i * s_j)`` instead of ``min(s_i, s_j)`` -- the weakest-member rule is an
    uncalibrated assumption, so geomean is the alternative to test robustness (plan R1).
    Repulsion is unscaled (steric).  Returns ``(M, 2)``.
    """
    centers = np.asarray(centers, dtype=float)
    d_eq = cfg.cell_spacing
    cutoff = d_eq + cfg.cc_adhesion_range
    d = centers[None, :, :] - centers[:, None, :]
    dist = np.linalg.norm(d, axis=2)
    safe = np.maximum(dist, 1e-12)
    unit = d / safe[:, :, None]
    rep = np.where(dist < d_eq, cfg.cc_repulsion * (d_eq - dist), 0.0)
    adh = np.where((dist >= d_eq) & (dist < cutoff), cfg.cc_adhesion * (dist - d_eq), 0.0)
    s = np.asarray(adhesion_scale, dtype=float)
    adh = adh * np.sqrt(np.maximum(s[:, None] * s[None, :], 0.0))
    mag = rep - adh
    mag[dist < 1e-9] = 0.0
    return -np.sum(mag[:, :, None] * unit, axis=1)


# --- cohesion / invasion-mode metrics (documented modelling choices) --------------
def radius_of_gyration(centers: np.ndarray) -> float:
    """RMS distance of cell centres from their centre of mass (cluster size)."""
    c = np.asarray(centers, dtype=float)
    if len(c) == 0:
        return 0.0
    return float(np.sqrt(np.mean(np.sum((c - c.mean(0)) ** 2, axis=1))))


def _adhesion_cutoff(cfg) -> float:
    return cfg.cell_spacing + cfg.cc_adhesion_range


def detached_fraction(centers: np.ndarray, cfg) -> float:
    """Fraction of cells whose NEAREST neighbour is beyond the adhesion cutoff.

    Past ``cell_spacing + cc_adhesion_range`` the piecewise adhesion is exactly zero, so
    such a cell feels no cohesion from any neighbour -- it has left the pack (single-cell
    escape; Ilina & Friedl 2020).  Returns 0 for a single-cell organoid.
    """
    c = np.asarray(centers, dtype=float)
    M = len(c)
    if M < 2:
        return 0.0
    d = np.linalg.norm(c[None, :, :] - c[:, None, :], axis=2)
    np.fill_diagonal(d, np.inf)
    return float(np.mean(d.min(axis=1) > _adhesion_cutoff(cfg)))


def largest_connected_component_fraction(centers: np.ndarray, cfg) -> float:
    """Size of the largest cohesive cluster / M (cells within the adhesion cutoff are
    linked).  ``1.0`` = fully cohesive front; small = fragmented / dispersed."""
    c = np.asarray(centers, dtype=float)
    M = len(c)
    if M == 0:
        return 0.0
    d = np.linalg.norm(c[None, :, :] - c[:, None, :], axis=2)
    adj = (d > 0.0) & (d <= _adhesion_cutoff(cfg))
    seen = np.zeros(M, dtype=bool)
    best = 0
    for s0 in range(M):
        if seen[s0]:
            continue
        stack = [s0]
        seen[s0] = True
        size = 0
        while stack:
            u = stack.pop()
            size += 1
            for v in np.flatnonzero(adj[u]):
                if not seen[v]:
                    seen[v] = True
                    stack.append(int(v))
        best = max(best, size)
    return best / M


def run_r1_invasion(cfg: "G5RevisionConfig" = None, seed=None, snapshots: bool = False) -> dict:
    """R1 invasion: the R0 force-consistent driver with a discrete partial-EMT phenotype
    that scales ONLY cell-cell adhesion (``emt_phenotype``), via the ``adhesion_scale``
    hook.  Clutch / drag / motor budget are identical to R0, so any change in invasion
    mode comes from adhesion alone.  Returns the R0 dict (frames already carry
    ``radius_of_gyration`` / ``detached_fraction`` / ``lcc_fraction``) plus the phenotype.
    """
    if cfg is None:
        cfg = r1_config()
    out = run_r0_invasion(cfg, seed=seed, snapshots=snapshots,
                          adhesion_scale=emt_phenotype, pair_rule=cfg.emt_pair_rule)
    out["emt_phenotype"] = emt_phenotype(out["centers0"], cfg)
    out["emt_fraction"] = float(cfg.emt_fraction)
    out["emt_adhesion_factor"] = float(cfg.emt_adhesion_factor)
    return out


def run_coupled_emt_invasion(cfg: "G5RevisionConfig" = None, seed=None,
                             snapshots: bool = False) -> dict:
    """COUPLED EMT-leader variant (biology-audit motivated): the SAME partial-EMT cells get
    BOTH reduced cell-cell adhesion (:func:`emt_phenotype`) AND elevated per-site motor stall
    (:func:`emt_site_stall`, via ``emt_traction_factor``) -- i.e. the partial-EMT cells ARE the
    high-traction leaders, matching the biology that leaders are often partial-EMT cells that
    loosen junctions AND pull harder / protrude more (biorxiv 2025; BIOLOGY_AUDIT.md Q4/B).

    ``emt_traction_factor == 1.0`` reduces EXACTLY to :func:`run_r1_invasion` (adhesion-only),
    so the two are a clean decoupled-vs-coupled comparison at matched emt_fraction/adhesion.
    Force-pair stays 0 (traction still emerges from the clutch force-velocity law, not a
    post-hoc multiply).  Runs on the isotropic-random network (or cued if ``radial_cue``).
    """
    if cfg is None:
        cfg = r1_config()
    net_mode = "cued" if getattr(cfg, "radial_cue", False) else "random"
    ss = emt_site_stall if float(getattr(cfg, "emt_traction_factor", 1.0)) != 1.0 else None
    out = run_r0_invasion(cfg, seed=seed, snapshots=snapshots,
                          adhesion_scale=emt_phenotype,
                          pair_rule=getattr(cfg, "emt_pair_rule", "min"),
                          network_mode=net_mode, site_stall=ss)
    out["emt_phenotype"] = emt_phenotype(out["centers0"], cfg)
    out["emt_cell_ids"] = emt_cell_ids(out["centers0"], cfg).tolist()
    out["emt_fraction"] = float(cfg.emt_fraction)
    out["emt_adhesion_factor"] = float(cfg.emt_adhesion_factor)
    out["emt_traction_factor"] = float(getattr(cfg, "emt_traction_factor", 1.0))
    return out


# =================================================================================
# R2 -- static leader NUMBER + LOCATION (functional role; per-site motor capacity)
# =================================================================================
# Leaders are a FUNCTIONAL role (front-cell traction), separate from EMT phenotype (R1).
# Their larger traction EMERGES from a higher PER-SITE F_stall via the force-velocity law
# v = v0 (1 - F/F_stall) (Chan & Odde 2008), NOT a post-hoc force multiply.  These are
# stall-aware copies of g4_v2 _independent_step / _shared_step: byte-identical to the
# originals except ``cfg.motor_stall_per_site`` -> the per-site ``stall`` array.  Reffay 2014
# (front cells exert larger traction) motivates leader_stall_factor>1 (a modelling choice).
def _independent_step_stall(cfg, state, substrate_speed, u_on, u_off, stall):
    force_before = cfg.clutch_stiffness * state.extension * state.bound
    traction_before = force_before.sum(axis=1)
    actin_speed = cfg.unloaded_actin_speed * np.maximum(0.0, 1.0 - traction_before / stall)
    relative = np.maximum(0.0, actin_speed - substrate_speed)
    state.extension[state.bound] += cfg.dt * np.repeat(
        relative[:, None], cfg.n_clutches_per_site, axis=1)[state.bound]
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


def _shared_step_stall(cfg, state, substrate_speed, u_on, u_off, stall):
    before_count = state.bound.sum(axis=1)
    site_stiffness = cfg.n_clutches_per_site * cfg.clutch_stiffness
    traction_before = np.where(before_count > 0, site_stiffness * state.site_extension, 0.0)
    actin_speed = cfg.unloaded_actin_speed * np.maximum(0.0, 1.0 - traction_before / stall)
    relative = np.maximum(0.0, actin_speed - substrate_speed)
    state.site_extension[before_count > 0] += cfg.dt * relative[before_count > 0]
    site_force = np.where(before_count > 0, site_stiffness * state.site_extension, 0.0)
    _pf, per_rate, _ = shared_load_hazard(site_force, before_count, cfg)
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


def _clutch_step_stall(cfg, state, substrate, step, active_mask, site_stall):
    """model._clutch_step with a per-site F_stall array (R2).  Same empty-sector masking."""
    S = state.bound.shape[0]
    u_on = _clutch_counter_uniforms(cfg, step, 0, S)
    u_off = _clutch_counter_uniforms(cfg, step, 1, S)
    if active_mask is not None:
        empty = ~active_mask
        u_on[empty] = 2.0
        state.bound[empty] = False
        state.extension[empty] = 0.0
        state.site_extension[empty] = 0.0
    stall = np.maximum(np.asarray(site_stall, dtype=float), 1e-9)
    if cfg.clutch_mode == "shared":
        return _shared_step_stall(cfg, state, substrate, u_on, u_off, stall)
    return _independent_step_stall(cfg, state, substrate, u_on, u_off, stall)


def leader_ids(centers: np.ndarray, cfg) -> np.ndarray:
    """Cell indices of the n_leaders leaders (deterministic).

    ``leader_location="cue_front"``: a LOCALIZED adjacent front -- among the outer (boundary)
    cells, the n_leaders whose outward azimuth is nearest ``cue_angle``.  ``"perimeter"``: the
    outermost n_leaders spread around the whole rim (control).  Empty when n_leaders<=0.
    """
    centers = np.asarray(centers, dtype=float)
    M = len(centers)
    n = int(getattr(cfg, "n_leaders", 0))
    if n <= 0 or M == 0:
        return np.zeros(0, dtype=int)
    n = min(n, M)
    radii = np.linalg.norm(centers, axis=1)
    if getattr(cfg, "leader_location", "cue_front") == "perimeter":
        return np.argsort(-radii)[:n]
    outer = np.argsort(-radii)[:max(n, int(round(0.5 * M)))]
    ang = np.arctan2(centers[outer, 1], centers[outer, 0])
    dang = np.abs(np.arctan2(np.sin(ang - cfg.cue_angle), np.cos(ang - cfg.cue_angle)))
    return outer[np.argsort(dang)[:n]]


def leader_site_stall(centers: np.ndarray, cfg) -> np.ndarray:
    """Per-site F_stall (length M*n_sec): base for followers, elevated for leader sites.

    ``budget_mode="fixed_per_leader"``: leader sites = base*leader_stall_factor (total front
    capacity GROWS with N_L).  ``"matched_total"``: the TOTAL leader stall budget is held
    CONSTANT vs N_L -- leader sites = base*leader_stall_factor/N_L -- so more leaders each get
    proportionally less, isolating leader NUMBER/geometry from total front traction (at large
    N_L a matched-total leader can drop toward/below base -- the intended budget trade-off).
    """
    centers = np.asarray(centers, dtype=float)
    M = len(centers)
    n_sec = cfg.n_contact_sectors
    base = float(cfg.motor_stall_per_site)
    stall = np.full(M * n_sec, base)
    lead = leader_ids(centers, cfg)
    if len(lead) == 0:
        return stall
    f = float(getattr(cfg, "leader_stall_factor", 1.0))
    per = base * f / len(lead) if getattr(cfg, "budget_mode", "matched_total") == "matched_total" else base * f
    for c in lead:
        stall[int(c) * n_sec:(int(c) + 1) * n_sec] = per
    return stall


def _nematic_rotation(source: float, target: float) -> float:
    """Signed rotation source->target modulo pi (fibre orientation is headless)."""
    delta = math.atan2(math.sin(target - source), math.cos(target - source))
    if delta > math.pi / 2:
        delta -= math.pi
    elif delta < -math.pi / 2:
        delta += math.pi
    return delta


def _cue_positions(positions, fibers, fixed, centers, cfg, organoid_outer):
    """Impose a ONE-SIDED radial tract by centroid-preserving, LENGTH-PRESERVING rigid
    rotation of eligible near-field fibres (port of v4E ``make_radial_tract_spec``,
    model.py:124).  Eligible = not boundary-anchored, centroid in the cue sector
    (``cue_angle`` +/- ``cue_half_width``) and radial band [organoid_outer, +cue_band].
    Rotations that would enter a cell disk or leave the box are SKIPPED, never clipped, so no
    hidden length change is introduced.  Returns ``(new_positions, tract_fids, report)``.  The
    cue is an EXPLICIT temporary assumption, NOT swirling-derived (Kolade sees no swirling).
    """
    pts = np.asarray(positions, dtype=float).copy()
    orig_pts = np.asarray(positions, dtype=float)
    fixed = np.asarray(fixed, dtype=bool)
    half = 0.5 * cfg.domain_size
    void = cfg.cell_radius + cfg.cell_clearance
    rng = np.random.default_rng(int(getattr(cfg, "cue_seed", 404)))
    tract: list = []
    attempts = 0
    lo, hi = organoid_outer, organoid_outer + cfg.cue_band
    for fid, ids_list in enumerate(fibers):
        ids = np.asarray(ids_list, dtype=int)
        if bool(np.any(fixed[ids])):
            continue
        seg = pts[ids]
        centroid = seg.mean(axis=0)
        r = float(np.linalg.norm(centroid))
        if r < lo or r > hi:
            continue
        polar = math.atan2(float(centroid[1]), float(centroid[0]))
        if abs(math.atan2(math.sin(polar - cfg.cue_angle),
                          math.cos(polar - cfg.cue_angle))) > cfg.cue_half_width:
            continue
        attempts += 1
        tangent = seg[-1] - seg[0]
        source = math.atan2(float(tangent[1]), float(tangent[0]))
        target = polar + float(rng.normal(0.0, cfg.cue_angular_sd))
        delta = _nematic_rotation(source, target)
        rot = np.array([[math.cos(delta), -math.sin(delta)], [math.sin(delta), math.cos(delta)]])
        proposed = centroid + (seg - centroid) @ rot.T
        dmin = np.sqrt(np.min(np.sum((proposed[:, None, :] - centers[None, :, :]) ** 2, axis=2), axis=1))
        if float(dmin.min()) < void - 1e-8 or float(np.max(np.abs(proposed))) > half + 1e-8:
            continue
        pts[ids] = proposed
        tract.append(fid)

    def contour(P):
        return np.array([float(np.linalg.norm(np.diff(P[np.asarray(i, dtype=int)], axis=0), axis=1).sum())
                         for i in fibers])
    dL = float(np.max(np.abs(contour(pts) - contour(orig_pts)))) if fibers else 0.0
    return pts, tract, {"eligible_attempts": attempts, "tract_fibers": len(tract),
                        "max_contour_length_change": dL}


def make_cued_organoid(cfg, seed=None):
    """:func:`make_random_organoid` + a one-sided imposed radial tract (R2, opt-in).

    Identical isotropic build, then :func:`_cue_positions` rotates the eligible cue-sector
    near-field fibres toward radial (length-preserving) BEFORE the Network is built, so the
    cued geometry is the rest state.  The cue is an EXPLICIT assumption, NOT swirling-derived.
    Returns ``(network, centers, gap_radius, report)`` with ``report["cue"]``.
    """
    cfg.validate()
    centers = hex_centers(cfg.organoid_radius, cfg.cell_spacing)
    organoid_outer = float(np.max(np.linalg.norm(centers, axis=1))) + cfg.cell_radius
    gap_radius = organoid_outer + cfg.gap
    base = cfg.seed if seed is None else int(seed)
    best = (-1, -1.0, None, None)
    for attempt in range(cfg.generation_attempts):
        seed_used = base + 7919 * attempt
        spec = _random_isotropic_spec(cfg, centers, seed_used)
        pts, _tract, cue_report = _cue_positions(spec.positions, spec.fibers, spec.fixed,
                                                 centers, cfg, organoid_outer)
        spec = NetworkSpec(pts, spec.fibers, spec.fixed, spec.contact_fibers, seed_used)
        network = Network(spec, cfg, crosslinks=[])
        network.crosslinks = build_crosslinks_grid(network)
        network.refresh_crosslink_arrays()
        if cfg.crosslink_fraction < 1.0 and network.crosslinks:
            rng = np.random.default_rng(seed_used + 101)
            keep = rng.random(len(network.crosslinks)) < cfg.crosslink_fraction
            network.crosslinks = [x for x, k in zip(network.crosslinks, keep) if k]
            network.refresh_crosslink_arrays()
        report = connectivity_report(network)
        report["cue"] = cue_report
        score = float(report["connected_fraction"])
        cc = 1 if report["contact_fibers_connected"] else 0
        if (cc, score) > (best[0], best[1]):
            best = (cc, score, network, report)
        if report["contact_fibers_connected"] and score >= cfg.required_connected_fraction:
            return network, centers, gap_radius, report
    if best[2] is not None:
        return best[2], centers, gap_radius, best[3]
    raise RuntimeError("cued organoid network percolation gate failed")


# --- R2 morphology metrics (documented modelling choices) -------------------------
def aspect_ratio(centers: np.ndarray) -> float:
    """Organoid elongation: sqrt(max/min principal variance of the cell cloud) (>=1)."""
    c = np.asarray(centers, dtype=float)
    if len(c) < 2:
        return 1.0
    w = np.linalg.eigvalsh(np.cov((c - c.mean(0)).T))
    w = np.maximum(w, 1e-12)
    return float(np.sqrt(w.max() / w.min()))


def leader_follower_separation(centers: np.ndarray, centers0: np.ndarray, cfg) -> float:
    """Mean OUTWARD advance along the cue axis of leaders minus followers (um).  >0 = leaders
    lead the front; the follower_lag is the same signed quantity (leaders ahead of followers)."""
    c = np.asarray(centers, dtype=float)
    c0 = np.asarray(centers0, dtype=float)
    lead = leader_ids(c0, cfg)
    if len(lead) == 0 or len(c) == 0:
        return 0.0
    axis = np.array([math.cos(cfg.cue_angle), math.sin(cfg.cue_angle)])
    adv = (c - c0) @ axis
    mask = np.zeros(len(c), dtype=bool)
    mask[lead] = True
    lead_adv = float(adv[mask].mean())
    foll_adv = float(adv[~mask].mean()) if (~mask).any() else 0.0
    return lead_adv - foll_adv


def follower_lag(centers: np.ndarray, centers0: np.ndarray, cfg) -> float:
    """Alias of :func:`leader_follower_separation` (how far followers lag the leaders, um)."""
    return leader_follower_separation(centers, centers0, cfg)


def strand_metrics(centers: np.ndarray, cfg):
    """Strand readout: connected chains of PROTRUDING cells (radius > median + cell_spacing)
    linked within the adhesion cutoff.  Returns (n_strands, mean_len_cells, max_len_cells)."""
    c = np.asarray(centers, dtype=float)
    M = len(c)
    if M == 0:
        return 0, 0.0, 0
    radii = np.linalg.norm(c, axis=1)
    idx = np.flatnonzero(radii > (np.median(radii) + cfg.cell_spacing))
    if len(idx) == 0:
        return 0, 0.0, 0
    cutoff = cfg.cell_spacing + cfg.cc_adhesion_range
    d = np.linalg.norm(c[idx][None, :, :] - c[idx][:, None, :], axis=2)
    adj = (d > 0.0) & (d <= cutoff)
    seen = np.zeros(len(idx), dtype=bool)
    sizes: list = []
    for s0 in range(len(idx)):
        if seen[s0]:
            continue
        stack = [s0]
        seen[s0] = True
        sz = 0
        while stack:
            u = stack.pop()
            sz += 1
            for v in np.flatnonzero(adj[u]):
                if not seen[v]:
                    seen[v] = True
                    stack.append(int(v))
        sizes.append(sz)
    return len(sizes), float(np.mean(sizes)), int(max(sizes))


def run_r2_invasion(cfg: "G5RevisionConfig" = None, seed=None, snapshots: bool = False) -> dict:
    """R2: force-consistent invasion with LOCALIZED leaders (number + location) whose higher
    traction EMERGES from per-site motor capacity (:func:`leader_site_stall`), optionally on a
    one-sided imposed radial cue (``radial_cue``).  EMT (R1 adhesion) is a SEPARATE knob: by
    default leaders are NOT low-adhesion (set ``emt_fraction`` for the coupled control).
    Force-pair stays 0.  Returns the R0/R1 dict + leader ids, per-site stall, cue report, and
    strand / leader-follower / aspect metrics.
    """
    if cfg is None:
        cfg = r1_config()
    net_mode = "cued" if getattr(cfg, "radial_cue", False) else "random"
    adh = emt_phenotype if float(getattr(cfg, "emt_fraction", 0.0)) > 0.0 else None
    boost = (int(getattr(cfg, "n_leaders", 0)) > 0
             and float(getattr(cfg, "leader_stall_factor", 1.0)) != 1.0)
    ss = leader_site_stall if boost else None
    out = run_r0_invasion(cfg, seed=seed, snapshots=snapshots, adhesion_scale=adh,
                          pair_rule=getattr(cfg, "emt_pair_rule", "min"),
                          network_mode=net_mode, site_stall=ss)
    c0, cf = out["centers0"], out["centers_final"]
    out["leader_ids"] = leader_ids(c0, cfg).tolist()
    out["site_stall"] = leader_site_stall(c0, cfg) if boost else None
    out["cue_report"] = out.get("connectivity", {}).get("cue")
    n_strands, mean_len, max_len = strand_metrics(cf, cfg)
    out["strand_count"] = int(n_strands)
    out["strand_mean_len"] = float(mean_len)
    out["strand_max_len"] = int(max_len)
    out["leader_follower_separation"] = leader_follower_separation(cf, c0, cfg)
    out["aspect_ratio"] = aspect_ratio(cf)
    return out


# =================================================================================
# R3 -- energy-based DYNAMIC leader switching (the PERSISTENCE layer; Zhang 2019)
# =================================================================================
# A front leader that keeps doing mechanical (reel) work DRAINS a per-cell energy; when it
# drops below a threshold the leader DEMOTES and a fresher front-follower (higher energy,
# gripping) is PROMOTED -- a relay-like handoff (Zhang et al. 2019: leader replacement in
# MDA-MB-231 spheroid/organoid collagen invasion, energy-depletion driven, leader lifetime
# ~120-480 min).  MODELLING CHOICE + explicit HYPOTHESIS: the motor mechanical power
# P = sum_s F_s * v_motor,s is used as an ATP-consumption PROXY -- it is NOT a measured or
# calibrated equation.  The energy params are calibrated to leader LIFETIME (so a handoff is
# visible within a run; here compressed to ~minutes vs Zhang's 120-480 min ORDER for
# demonstration), NOT to measured ATP.  EMT phenotype (adhesion) is UNCHANGED by switching.
# Built on the R2 localized-leader / per-site-stall machinery; force-pair stays 0.
def motor_power_per_cell(cfg, patches, site_force, site_stall, n_sec, M):
    """Per-cell motor mechanical power (ATP-consumption PROXY -- a HYPOTHESIS):
    ``P_i = sum_{s in i} F_s * v_motor,s`` with ``v_motor,s = v0*max(0, 1 - F_s/F_stall,s)``
    (the SAME force-velocity law the clutch step uses; Chan & Odde 2008).  Gripping loaded
    sites only; a non-gripping cell has ``P = 0``.  Units nN*um/s."""
    v0 = float(cfg.unloaded_actin_speed)
    stall = np.maximum(np.asarray(site_stall, dtype=float), 1e-9)
    P = np.zeros(M)
    for s, patch in enumerate(patches):
        if patch is None or site_force[s] <= 0.0:
            continue
        v = v0 * max(0.0, 1.0 - float(site_force[s]) / float(stall[s]))
        P[s // n_sec] += float(site_force[s]) * v
    return P


def _front_pool_mask(centers, cfg):
    """Boolean mask of FRONT-eligible cells: the outer (boundary) cells, restricted to the
    cue sector (|azimuth - cue_angle| <= cue_half_width) when ``leader_location='cue_front'``
    (all outer cells for ``'perimeter'``)."""
    centers = np.asarray(centers, dtype=float)
    M = len(centers)
    mask = np.zeros(M, dtype=bool)
    if M == 0:
        return mask
    radii = np.linalg.norm(centers, axis=1)
    outer = np.argsort(-radii)[:max(1, int(round(0.5 * M)))]
    if getattr(cfg, "leader_location", "cue_front") == "perimeter":
        mask[outer] = True
        return mask
    ang = np.arctan2(centers[outer, 1], centers[outer, 0])
    d = np.abs(np.arctan2(np.sin(ang - cfg.cue_angle), np.cos(ang - cfg.cue_angle)))
    mask[outer[d <= cfg.cue_half_width]] = True
    return mask


def update_active_leaders(centers, energy, current_leaders, gripping_mask, cfg):
    """Dynamic active-leader set with HYSTERESIS (Zhang 2019 relay).  A candidate must be a
    FRONT cell (:func:`_front_pool_mask`) AND currently gripping collagen.  A current leader
    STAYS active until its energy < ``energy_off`` (or it loses grip / leaves the front) then
    DEMOTES; free slots (< ``n_leaders``) are filled by the highest-energy eligible candidate
    with energy > ``energy_on`` (PROMOTE).  ``energy_on > energy_off`` prevents flip-flop.
    Returns ``(active_ids sorted, events)`` where events is a list of ``('demote'|'promote', cell)``."""
    M = len(centers)
    n_max = int(getattr(cfg, "n_leaders", 0))
    if n_max <= 0 or M == 0:
        return np.zeros(0, dtype=int), []
    eligible = _front_pool_mask(centers, cfg) & np.asarray(gripping_mask, dtype=bool)
    E = np.asarray(energy, dtype=float)
    e_on = float(cfg.energy_on)
    e_off = float(cfg.energy_off)
    events = []
    active = []
    for c in [int(x) for x in current_leaders]:
        if eligible[c] and E[c] >= e_off:
            active.append(c)
        else:
            events.append(("demote", c))
    if len(active) < n_max:
        cand = [int(c) for c in np.flatnonzero(eligible)
                if int(c) not in active and E[c] > e_on]
        cand.sort(key=lambda c: -E[c])
        for c in cand[: n_max - len(active)]:
            active.append(int(c))
            events.append(("promote", int(c)))
    return np.asarray(sorted(active), dtype=int), events


def _site_stall_from_ids(active_ids, cfg, M):
    """Per-site F_stall (length ``M*n_sec``) elevating the CURRENT active leaders -- dynamic
    version of :func:`leader_site_stall` (same matched_total / fixed_per_leader budget)."""
    n_sec = cfg.n_contact_sectors
    base = float(cfg.motor_stall_per_site)
    stall = np.full(M * n_sec, base)
    ids = [int(c) for c in active_ids]
    if not ids:
        return stall
    f = float(getattr(cfg, "leader_stall_factor", 1.0))
    per = (base * f / len(ids)
           if getattr(cfg, "budget_mode", "matched_total") == "matched_total" else base * f)
    for c in ids:
        stall[c * n_sec:(c + 1) * n_sec] = per
    return stall


def run_r3_invasion(cfg: "G5RevisionConfig" = None, seed=None, snapshots: bool = False) -> dict:
    """R3: force-consistent invasion with ENERGY-BASED DYNAMIC leader switching (persistence).

    Each step the per-cell motor power ``P`` (ATP proxy; a HYPOTHESIS) drains a per-cell energy
    ``E`` that recovers toward ``energy_E0`` with ``energy_tau_rec``; the ACTIVE leader set is
    recomputed with hysteresis (:func:`update_active_leaders`) so a drained front leader hands
    off to a fresher front-follower, and the per-site ``F_stall`` is rebuilt from the CURRENT
    leaders.  ``leader_switching=False`` -> identical to :func:`run_r2_invasion`.  Force-pair
    stays 0.  Imposed cue NOT swirling; EMT adhesion unchanged by switching.  Personal testing."""
    if cfg is None:
        cfg = r1_config()
    if not getattr(cfg, "leader_switching", False):
        return run_r2_invasion(cfg, seed=seed, snapshots=snapshots)

    net_mode = "cued" if getattr(cfg, "radial_cue", False) else "random"
    _builder = {"random": make_random_organoid, "cued": make_cued_organoid}.get(
        net_mode, make_random_organoid)
    network, centers, gap_radius, report = _builder(cfg, seed=seed)
    centers = centers.copy()
    centers0 = centers.copy()
    M = len(centers)
    n_sec = cfg.n_contact_sectors
    organoid_center = np.zeros(2)
    stepper = OrganoidStepper(network, centers)
    reach = cfg.cell_radius + cfg.contact_width + 2.0
    # R1 EMT adhesion hook; leaders are NOT auto-low-adhesion (phenotype separate from role)
    adh = emt_phenotype if float(getattr(cfg, "emt_fraction", 0.0)) > 0.0 else None
    adh_scale = (leader_adhesion_scale(centers0, cfg) if adh is None
                 else np.asarray(adh(centers0, cfg), dtype=float))
    cc_force = ((lambda ctr: _cell_cell_forces_geomean(ctr, cfg, adh_scale))
                if getattr(cfg, "emt_pair_rule", "min") == "geomean"
                else (lambda ctr: cell_cell_forces(ctr, cfg, adh_scale)))
    candidates = cell_candidate_fibers(network, centers, reach)
    patches, _sc = organoid_clutch_patches(network, centers, cfg, candidates)
    S = len(patches)
    active_mask = _clutch_active_mask(patches)
    state = _new_clutch_state(S, cfg)
    substrate = np.zeros(S)
    site_force = np.zeros(S)
    v_cell = np.zeros((M, 2))
    reaction = np.zeros((M, 2))
    n_relocations = 0

    # R3 energy + dynamic leaders (seed with the R2 static pick)
    energy = np.full(M, float(cfg.energy_E0))
    active_leaders = [int(c) for c in leader_ids(centers0, cfg)]
    n_switches = 0
    switch_log: list = []
    leader_since = {c: 0.0 for c in active_leaders}
    cue_axis = np.array([math.cos(cfg.cue_angle), math.sin(cfg.cue_angle)])

    nsteps = int(round(cfg.duration / cfg.dt))
    every = max(1, int(round(cfg.sample_interval / cfg.dt)))
    contact_every = max(1, int(round(cfg.contact_update_interval / cfg.dt)))
    frames: list = []
    bead_snaps: list = []
    cell_snaps: list = []
    max_residual = 0.0

    def front_advance():
        fp = _front_pool_mask(centers0, cfg)
        return float(np.mean((centers[fp] - centers0[fp]) @ cue_axis)) if fp.any() else 0.0

    def sample_frame(time):
        prof = radial_alignment_profile(network, organoid_center)
        clutch_ecm = np.zeros((M, 2))
        for s, patch in enumerate(patches):
            if patch is None or site_force[s] <= 0.0:
                continue
            clutch_ecm[s // n_sec] += site_force[s] * patch.normal_in
        residual = float(np.max(np.linalg.norm(reaction + clutch_ecm, axis=1))) if M else 0.0
        lead_E = [float(energy[c]) for c in active_leaders]
        return {
            "time": time,
            "global_radial_order": prof["global_radial_order"],
            "mean_cell_radial_disp": float(np.mean(
                np.linalg.norm(centers, axis=1) - np.linalg.norm(centers0, axis=1))),
            "max_cell_disp": float(np.max(np.linalg.norm(centers - centers0, axis=1))),
            "radius_of_gyration": radius_of_gyration(centers),
            "detached_fraction": detached_fraction(centers, cfg),
            "lcc_fraction": largest_connected_component_fraction(centers, cfg),
            "force_pair_residual": residual,
            "floppiness_index": floppiness_index(network),
            "total_traction": float(site_force.sum()),
            # R3:
            "active_leaders": list(active_leaders),
            "n_active_leaders": len(active_leaders),
            "energy_min": float(energy.min()),
            "energy_mean": float(energy.mean()),
            "leader_energy_mean": float(np.mean(lead_E)) if lead_E else 0.0,
            "n_switches": int(n_switches),
            "front_advance": front_advance(),
            "leader_follower_separation": leader_follower_separation(centers, centers0, cfg),
        }

    for step in range(nsteps + 1):
        time = step * cfg.dt
        if step % every == 0:
            fr = sample_frame(time)
            max_residual = max(max_residual, fr["force_pair_residual"])
            frames.append(fr)
            if snapshots:
                bead_snaps.append(network.r.copy())
                cell_snaps.append(centers.copy())
        if step == nsteps:
            break

        # dynamic per-site F_stall from the CURRENT active leaders (rebuilt every step)
        site_stall = _site_stall_from_ids(active_leaders, cfg, M)

        # 1) clutch step (per-site stall) -> emergent traction
        _, site_force, breaks, binds, site_fail = _clutch_step_stall(
            cfg, state, substrate, step, active_mask, site_stall)
        # 2) project onto ECM
        active = _project_site_forces(network, patches, site_force)
        # 3) force-pair reaction + capped overdamped motion
        reaction = per_cell_clutch_reaction(patches, site_force, n_sec, M)
        f_cell = reaction + cc_force(centers)
        v_cell = f_cell / cfg.cell_drag
        speed = np.linalg.norm(v_cell, axis=1)
        over = speed > cfg.max_cell_speed
        if np.any(over):
            v_cell[over] *= (cfg.max_cell_speed / speed[over])[:, None]
        centers += cfg.dt * v_cell
        stepper.centers[:] = centers
        # 4) ECM step + relative substrate speed
        stepper.step(active, cfg.dt)
        site_centers = np.repeat(centers, n_sec, axis=0)
        v_cell_per_site = np.repeat(v_cell, n_sec, axis=0)
        substrate = _relative_substrate_speeds_mc(
            network, stepper.velocity, patches, site_centers, v_cell_per_site)
        if step and step % contact_every == 0:
            candidates = cell_candidate_fibers(network, centers, reach)
        # 5) event-driven relocation + reset of fully-failed sites
        if np.any(site_fail):
            for s in np.flatnonzero(site_fail):
                c = int(s // n_sec)
                newp = _relocate_site(network, centers[c], patches, int(s), c,
                                      n_sec, candidates[c], cfg)
                if newp is not None:
                    patches[int(s)] = newp
                _reset_site_state(state, int(s))
                n_relocations += 1
            active_mask = _clutch_active_mask(patches)

        # --- R3: energy update (ATP proxy) + dynamic leader switching ---
        P = motor_power_per_cell(cfg, patches, site_force, site_stall, n_sec, M)
        energy += cfg.dt * ((cfg.energy_E0 - energy) / cfg.energy_tau_rec - P / cfg.energy_cap)
        np.clip(energy, 0.0, None, out=energy)
        gripping = active_mask.reshape(M, n_sec).any(axis=1)
        new_active, events = update_active_leaders(centers, energy, active_leaders, gripping, cfg)
        for kind, c in events:
            switch_log.append((float(time), kind, int(c)))
            if kind == "promote":
                n_switches += 1
                leader_since[c] = float(time)
        active_leaders = [int(c) for c in new_active]

    lifetimes = [float(cfg.duration - leader_since.get(c, cfg.duration)) for c in active_leaders]
    ns, ml, mx = strand_metrics(centers, cfg)
    return {
        "config": asdict(cfg),
        "centers0": centers0, "centers_final": centers,
        "organoid_center": organoid_center, "connectivity": report,
        "n_cells": M, "n_beads": len(network.r), "n_fibers": len(network.fibers),
        "n_clutch_sites": S, "clutch_mode": cfg.clutch_mode,
        "max_force_pair_residual": float(max_residual),
        "cumulative_slips": int(state.cumulative_slips),
        "cumulative_site_failures": int(state.cumulative_site_failures),
        "n_relocations": int(n_relocations),
        "frames": frames,
        "bead_snapshots": np.asarray(bead_snaps) if snapshots else None,
        "cell_snapshots": np.asarray(cell_snaps) if snapshots else None,
        "final_positions": network.r.copy(), "initial_positions": network.r0.copy(),
        "edges": network.edges.copy(),
        # R3 summary
        "n_switches": int(n_switches),
        "switch_log": switch_log,
        "final_active_leaders": list(active_leaders),
        "initial_leaders": [int(c) for c in leader_ids(centers0, cfg)],
        "final_leader_lifetimes": lifetimes,
        "energy_final": energy.copy(),
        "strand_count": int(ns), "strand_mean_len": float(ml), "strand_max_len": int(mx),
        "aspect_ratio": aspect_ratio(centers),
        "leader_follower_separation": leader_follower_separation(centers, centers0, cfg),
    }
