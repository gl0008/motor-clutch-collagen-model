# G5 Code Assessment — engine audit for organoid scale-up

> Audit of the existing generations to decide the G5 foundation (2026-09-02). Feeds `docs/G5_organoid_plan.md` §1, §6.

## Shared core (the key fact)
All generations sit on ONE physics engine via subclassing — not forks:
```
CollagenConfig / Network  (generations/g2_corrected/common/model.py, 857 lines, 92 tests)
   ├── SpheroidConfig  → g3   (reuses g2 Network unchanged; adds spheroid+protrusions; fast_advance = bit-identical physics)
   └── G4Config → G4V2Config → g4_v2 (imports g2 Network; adds Numba integrator + Bell/shared-load clutch + rigid-body motion)
```

## g2 core (reuse as-is)
- Overdamped: `zeta·dr/dt = F_stretch + F_bend + F_crosslink + F_repulsion + F_active` (`model.py:6`).
- Data: `r/r0` (N,2) beads; `fibers` bead chains; `edges`/`triplets`; crosslinks pin material points
  `(edge, alpha)` on two fibers; `fixed` mask anchors outer boundary. **2D only.**
- Forces vectorized (`elastic_forces`, `np.add.at`); stretch has 10% compression asymmetry (microbuckling);
  bending harmonic about rest curvature; Gaussian kernel is a **projection tool, not a spring**.
- Forward-Euler, dt=0.005 s. Units µm/nN/s. Pure NumPy, no Numba. Baseline 99 fibers / 180 µm.

## g4_v2 additions (lift these parts)
Numba `_advance_numba` (allocation-free, same 5 forces) · Bell `k_off=k_off0·exp(F/F_b)` · shared-load
vs independent clutch cluster · crosslink-graph distance classification · counter-based deterministic RNG ·
released rigid cell (translation+rotation). **But** `run_clutch` is a 140-line monolith and there are
**two integrator copies** to keep synced — messy; mine parts, don't extend in place.

## The blocker: single-cell everywhere
`center`/`cell_radius`/`theta` are scalars; `repulsion_forces` assumes one center; clutch state arrays
are `(sectors, clutches)` = one cell. **No multi-cell support anywhere** (grep n_cells/centers/organoid
→ nothing). g3's "spheroid" = one big rigid ball, not an assembly.

## Verdict
Reuse g2 `Network` (validated ECM core) → lift g4_v2 Numba integrator + clutch equations → write a NEW
multi-cell driver in `generations/g5_organoid/`. Port g3's fibre-free-gap generator.

## Performance walls (setup + contact detection, NOT the force law)
- `build_crosslinks` O(F²) Python double loop; `make_network_spec` whole-network retry; `contact_patches`
  loops every fiber per update (×M cells = M×F).
- Fix: spatial **neighbor grid / cell list** for crosslinking + per-cell contact; generalize Numba
  repulsion to M centers; keep Numba integrator + g3 `np.bincount` scatter. Ship a perf smoke test as the
  G5-A gate.

## Parameter provenance
Full tables in `docs/g4_parameter_provenance.md` + CLAUDE.md §8 (Saraswathibhatla 2025 SI Table 2;
Adebowale 2021 SI Table 4; Bell 1978). g4 baseline: R=10 µm, σ=1.5 µm, contact shell 3 µm, bead drag
180 nN·s/µm, softened collagen modulus 3 MPa, crosslink prob 0.35. Clutch: N=12/site, k_c=2 nN/µm,
k_on=0.055, k_off0=0.018 s⁻¹, F_b=1.5 nN, v0=0.025 µm/s, F_stall=8 nN/site, cell drag 600 nN·s/µm.
