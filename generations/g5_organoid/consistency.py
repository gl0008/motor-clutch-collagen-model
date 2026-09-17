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

from dataclasses import asdict
import math

import numpy as np

from generations.g5_organoid.model import (  # noqa: E402
    OrganoidConfig,
    parameter_variant,
    make_organoid,
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
def run_r0_invasion(cfg: OrganoidConfig = None, seed=None, snapshots: bool = False) -> dict:
    """Force-consistent multicellular Stage-D-with-clutch invasion (the G5-R0 gate).

    Ports G4D's ``run_motor_clutch`` loop ORDER to M cells: sample frame -> clutch step
    -> project active ECM force -> per-cell force-pair reaction + capped overdamped cell
    motion -> ECM step -> relative substrate speed -> event-driven relocation (+reset) of
    only fully-failed sites.  See module docstring.  ``leader_*`` are out of scope (R0).
    """
    if cfg is None:
        cfg = r0_config()

    network, centers, gap_radius, report = make_organoid(cfg, seed=seed)
    centers = centers.copy()
    centers0 = centers.copy()
    M = len(centers)
    n_sec = cfg.n_contact_sectors
    organoid_center = np.zeros(2)
    stepper = OrganoidStepper(network, centers)
    reach = cfg.cell_radius + cfg.contact_width + 2.0
    adh_scale = leader_adhesion_scale(centers0, cfg)         # defaults -> all ones (R0)
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
        f_cc = cell_cell_forces(centers, cfg, adh_scale)
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
        f_cell = reaction + cell_cell_forces(centers, cfg, adh_scale)
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
