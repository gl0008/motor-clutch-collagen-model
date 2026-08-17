# G3 validation targets

These are causal-mechanism and numerical gates. Passing them does not validate realistic 3D
tumor migration. Seeds 0–19 are reserved for calibration/debugging; seeds 1000–1099 are the
untouched final-validation set.

## G3A — spatial material-point clutches

| Control | Gate |
|---|---|
| Material-point persistence | `(fiber_id, segment_id, alpha)` remains fixed until Bell unbinding |
| Material interpolation | The attachment position always equals `(1-alpha) r_a + alpha r_b` |
| No nearby fibre | Clutch remains unbound and traction is zero |
| Force conservation | Relative error below `1e-10` |
| First-moment/torque conservation | Relative error below `1e-8` |
| Rotation covariance | Rigidly rotating the fixture rotates point and projected forces identically |
| Pull response | Pull-induced FOI change exceeds no-pull numerical drift by at least 10× |

## G3B — protrusion selection

| Control | Gate |
|---|---|
| Balanced fixture | No fixed `+x` peak |
| Isotropic ensemble | Ensemble polar resultant / mean individual resultant below 0.1 |
| Aligned ECM | Protrusion axis follows the collagen director |
| Nematic symmetry | Positive/negative director choices each occupy 40–60% of trials |
| Rotated aligned ECM | Estimated director rotates `30° ± 5°` with the matrix |
| Feedback ablation | `beta_geometry=beta_traction=0` reduces guidance by at least 50% |
| Empty ECM | Protrusions turn over, but no clutch traction develops |

## G3C — rigid-body motion

| Control | Gate |
|---|---|
| Empty ECM | Numerical-zero displacement and rotation |
| Symmetric attachments | Net torque approaches zero |
| Asymmetric attachments | Nonzero rotation |
| Mirror fixture | Rotation sign and lateral trajectory reverse |
| Rotated ECM | Trajectory distribution rotates with the matrix |
| Isotropic ensemble | Mean displacement approaches zero while individual runs may persist |
| Drag sweep | Speed changes, but axis and symmetry conclusions do not |
| Hidden-drive audit | No `v0`, fixed `+x` force, or `polarity_probability=0.65` |

## Elastic load–unload diagnostic

FOI and κ follow Nam et al. 2016. For a finite synthetic fixture, κ uses the measured initial
FOI rather than the ideal random-network value `2/pi`. When the pull signal is resolved, the
engineering gate is `kappa < 0.1`. A trajectory that has not equilibrated by 600 s is reported
as `unresolved_recovery`, not plasticity.

The complete 100-seed, 600-s G3B/G3C ensemble has not yet been run. Current saved outputs are
mechanism smoke tests only.
