# Saved default-run results

These are deterministic outputs for the checked-in seeds and parameters.  They
verify implementation behavior; they are not experimental predictions.

| Version | Default saved result | Interpretation boundary |
|---|---|---|
| V1 | local alignment 0.310 → 0.501; RMS displacement 0.237 model length | dimensionless five-fiber demonstration; Gaussian acts on all beads |
| V2 | `S_r` −0.11854 → −0.10328; RMS displacement 0.10482 µm; 128 permanent links | weak radial reorientation under a single right-sector pull |
| V3 fixed | cell displacement exactly zero; final reaction `(-2.198, 0.289)` model force | reaction is supplied by the fixed-position constraint |
| V3 moving | final cell center `(-0.18043, 0.03093)` µm | one stochastic trajectory, not a migration-speed estimate |
| V4 elastic | residual `ΔS_r = 0.008921`; RMS displacement 0.07762 µm | the 4 s unload is not a full equilibrium recovery |
| V4 plastic | 7 breaks, 2 re-forms; residual `ΔS_r = 0.008910`; RMS displacement 0.07763 µm | indistinguishable from elastic control: current turnover law does not establish permanent remodeling |

## What should be reported next

- V2: radial profiles of displacement/alignment, no-crosslink control, domain
  and bead-resolution convergence, and energy/work balance.
- V3: distributions across many seeds of displacement, persistence, occupancy,
  traction imbalance, and reaction force.
- V4: a longer equilibrium-unload protocol plus preregistered comparisons of
  same-neighborhood rest reset, new-link formation, sliding/merging, and fiber
  plastic lengthening.
