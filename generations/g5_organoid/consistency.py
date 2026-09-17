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
    _new_clutch_state,
    _clutch_active_mask,
    cell_cell_forces,
    leader_adhesion_scale,
    radial_alignment_profile,
)
from common.model import Network, NetworkSpec, connectivity_report  # noqa: E402


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
                    network_mode: str = "random") -> dict:
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
    # default; radial alignment is an OUTPUT); "corona" = legacy model.make_organoid (biased).
    _builder = make_random_organoid if network_mode == "random" else make_organoid
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

        # 1) clutch bundles load / Bell-slip / rebind (reused g4_v2 law via model._clutch_step)
        _, site_force, breaks, binds, site_fail = _clutch_step(
            cfg, state, substrate, step, active_mask)

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
