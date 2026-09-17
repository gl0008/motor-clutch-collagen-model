# G5-R0 — mechanical-consistency gate (findings)

> **Personal testing, NOT confirmed findings** (CLAUDE.md §7.5). No swirling is claimed.
> R0 is the advisor-requested consistency gate that must pass before any EMT/leader work
> (R1–R3). Code: `generations/g5_organoid/consistency.py` (new; imports G5, edits nothing —
> `model.py`/`visualize.py`/`test_g5.py` stay byte-for-byte, so old G5 is the control).

## What R0 does (ports Gloria's accepted **G4D** single-cell force-consistent loop to M cells)

`run_r0_invasion` replaces the buggy moving-cell clutch path (`model._run_invasion_with_clutch`,
which carried loading history across periodic re-selection) with a faithful multicellular port of
G4D (`generations/g4_v4_single_cell_force_alignment/model.py`):

1. **Force-pair (Newton 3rd law):** per-cell reaction = −Σ(site_force·normal_in), the SAME vectors
   projected onto the ECM → the cell feels the exact equal-and-opposite reaction.
2. **Relative substrate speed:** subtracts each site's owning-cell velocity, so rigid cell
   translation is not misread as clutch loading (G4D `_relative_substrate_speeds`).
3. **Event-driven relocation + reset:** a grip changes fibre ONLY after that site fully fails, then
   picks a fresh patch by G4D's tangent-guidance rule and RESETS clutch state (no history carried).
4. **Separated force channels:** reports per-cell `F_clutch / F_steric / F_cell_cell` (steric is
   ECM-only — cell + collagen may move the same way — while clutch is a strict pair) + a
   parallel/perp-to-fibre clutch decomposition.

**Locked baseline** (`r0_config()`): `clutch_mode="shared"` (Erdmann-Schwarz 2004; G4D V4-A012),
`cell_drag=600` nN·s/µm, `max_cell_speed=0.012` µm/s (G4D defaults), softened collagen (3 MPa,
links 10 nN/µm), `strain_stiffening=False`, `plasticity=False` (like-for-like elastic gate).

## Gate results — ALL PASS (independently re-verified)

`python -m pytest generations/g5_organoid/tests/test_r0_consistency.py generations/g5_organoid/tests/test_g5.py`
→ **22 passed** (7 R0 gate + 15 unchanged G5).

| Invariant | Result |
|---|---|
| Force-pair residual `max_c \|reaction_c + Σ clutch vectors\|` | **0.000e+00** (exact, by construction) — held every frame through the full 2 h |
| Symmetric reaction ⇒ no net translation | PASS (unit) · symmetric-organoid centroid drift < 0.5 µm |
| Passive (`clutch_on_rate=0`) ⇒ ECM inert, cells frozen | PASS (bead disp < 1e-6, 0 failures) |
| Event-driven relocation resets state; `n_relocations == cumulative_site_failures` | PASS (1018 == 1018 at 2 h) |
| Rigid co-moving translation ⇒ ~0 relative substrate speed | PASS (< 1e-12) |

## 2 h diagnostic (headline) — coherent strain, NOT spaghetti

`python generations/g5_organoid/experiments/r0_diagnostic.py headline`
(19 cells, 19 122 beads, 228 sites, shared clutch, 7200 s, seed 23; 348 s wall)

- **max_force_pair_residual = 0.000e+00** across all 121 frames.
- **floppiness_index = 0.0000** (fraction of fibre segments strained > 50%); segment strain
  mean 0.1 %, std 0.8 %, **max 9 %** → fibres barely deform: **coherent strain, confirmed
  visually** (near-field collagen reeled into straight radial tracts, organoid stays cohesive).
- **Honest metric caveat:** the raw `spaghetti_index` (= rms bead displacement / bead_spacing)
  rose 0 → 1.48, but that measures *coherent* bulk displacement (cells reel collagen inward,
  reach 125 µm), **not** floppiness. Fixed by adding `floppiness_index` (fibre-deformation based);
  use it as the spaghetti test.
- Traction is now realistic: **F_clutch mean 11.8 / max 26.6 nN/cell** (closer to Mark 2020 than
  the old ~10 nN), F_steric 3.5 nN, F_cell_cell 11.7 nN; clutch ∥/⊥ = 256/63 nN (mostly along fibre).
- Mean invasion **+2.4 µm**, max 7.3 µm over 2 h (elastic, speed-capped).
- Artifacts: `output/r0_diagnostic_2h.gif` (121 frames = 2 h), `output/r0_diagnostic_2h.npz`.

## Baseline sweep — coherent across the ECM window

`python generations/g5_organoid/experiments/r0_diagnostic.py sweep` (150 s, seed 23): over
collagen_modulus {1.5,3,6} × crosslink_fraction {0.6,0.85} × pull {12,24}, **force-pair residual
= 0 and spaghetti_index < 0.13 in every cell**; `collagen_modulus` barely matters (stretch stiffness
is not the lever here — consistent with prior G5 findings). Relocations == failures throughout.

## Honest limitations / scope

- Moderate scale (~19 cells) for the 2 h run (machine is memory-tight; a 12-run sweep got
  OOM-killed). Full 43-cell scale is a rerun away.
- Elastic, speed-capped (no plasticity/stiffening by design for a like-for-like gate) → invasion
  saturates and stays small; the point of R0 is *consistency*, not large motility.
- R0 reuses G5's `organoid_clutch_patches` for INITIAL site selection only; the periodic buggy
  re-selection is replaced by the event-driven rule.
- EMT/leader/switching are OUT of R0 scope (R1–R3); `leader_*` left at defaults (no-op).

## Reproduce

```bash
python -m pytest generations/g5_organoid/tests/test_r0_consistency.py   # 7 gate invariants
python -m pytest generations/g5_organoid/tests/test_g5.py               # 15 (unchanged control)
python generations/g5_organoid/experiments/r0_diagnostic.py sweep       # lock baseline (fast)
python generations/g5_organoid/experiments/r0_diagnostic.py headline    # 2 h gif + npz + numbers
```
