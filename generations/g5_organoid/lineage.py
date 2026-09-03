"""G4 -> G5 parameter lineage: frozen-parameter preset and transfer controls.

Motivation (advisor discussion, 2026-09): the exploratory G5 default changed several
ECM / contact / mobility parameters *at the same time* as going from 1 cell to M cells,
so a raw G4->G5 difference cannot be attributed to multicellularity alone.  This module
does NOT touch the exploratory G5 defaults (version-preservation); it adds:

  * ``g4_frozen_config`` -- an OrganoidConfig whose every non-multicellular parameter is
    pinned to its G4 value, so the ONLY things G5 adds are cell number, packing geometry,
    cell-cell repulsion/adhesion, and per-cell ECM reaction.
  * ``fiber_density`` -- fibers per um^2, so networks are compared by DENSITY not by count.
  * ``matched_fiber_count`` -- fiber count that reproduces a target density in a given box.

Verified G4 values (generations/g4_interactive_calibration/model.py) vs current G5 default:

    parameter          G4      G5-default   category
    bead_drag          180.0   120.0        collagen timescale   -> FREEZE
    bead_spacing       1.0     0.75         fibre discretisation -> FREEZE
    cell_drag          600.0   400.0        mobility             -> FREEZE
    max_cell_speed     0.012   0.02         mobility             -> FREEZE
    total_pull_force   24.0    12.0         per-cell traction    -> FREEZE (constant-pull control)
    cell_radius        10.0    9.0          geometry             -> FREEZE (or justify 9 = 18um dia)
    contact_width      3.0     8.0          contact geometry     -> FREEZE (fix corona, not width)
    crosslink retain   0.35(*) 0.85         DIFFERENT MECHANISM  -> match by density, not value
    gaussian_sigma     1.5     1.5          consistent
    collagen_modulus   3.0     3.0          consistent
    crosslink_stiffness10.0    10.0         consistent
    compression_ratio  0.10    0.10         consistent
    clutch params      ==      ==           consistent (Adebowale 2021 SI T4)

    (*) G4 uses ``crosslink_probability`` (keep fraction of auto links); G5 uses
        ``crosslink_fraction`` -- the generators differ, so the raw number is not
        directly comparable.  Match measured density instead.
"""

from __future__ import annotations

from dataclasses import replace

from generations.g5_organoid.model import OrganoidConfig

# --- G4 reference values (source: g4_interactive_calibration/model.py) --------------
G4_FROZEN = dict(
    bead_drag=180.0,          # G4 line 83 (vs G5 parent 120)
    bead_spacing=1.0,         # G4 line 62 (vs G5 parent 0.75)
    cell_drag=600.0,          # G4 line 99
    max_cell_speed=0.012,     # G4 line 101
    total_pull_force=24.0,    # G4 line 79 (per cell; constant-pull control)
    cell_radius=10.0,         # G4 line 58
    contact_width=3.0,        # G4 line 76
    gaussian_sigma=1.5,       # already consistent
    collagen_modulus_mpa=3.0, # already consistent
    crosslink_stiffness=10.0, # already consistent
    compression_ratio=0.10,   # already consistent
)

# G4 network scale, for density matching
G4_N_FIBERS = 99
G4_DOMAIN = 180.0


def fiber_density(n_fibers: int, domain_size: float) -> float:
    """Fibres per um^2 -- the density-based way to compare networks across box sizes."""
    return n_fibers / (domain_size ** 2)


def matched_fiber_count(target_density: float, domain_size: float) -> int:
    """Fibre count that reproduces ``target_density`` in a ``domain_size`` box."""
    return int(round(target_density * domain_size ** 2))


def g4_frozen_config(**overrides) -> OrganoidConfig:
    """OrganoidConfig with every non-multicellular parameter pinned to its G4 value.

    Multicellular knobs (leader_fraction, cc_adhesion, cc_repulsion, organoid_radius,
    cell_spacing, n_corona_fibers, ...) and scale knobs (domain_size, n_fibers) stay at
    whatever you pass in.  Fibre count should be set via :func:`matched_fiber_count` so the
    DENSITY -- not the count -- matches G4.  Pass ``**overrides`` to layer an experiment on
    top (e.g. clutch_dynamics=True, leader_fraction=0.1)."""
    params = dict(G4_FROZEN)
    params.update(overrides)
    return OrganoidConfig(**params)


G4_FIBER_DENSITY = fiber_density(G4_N_FIBERS, G4_DOMAIN)  # ~3.06e-3 /um^2
