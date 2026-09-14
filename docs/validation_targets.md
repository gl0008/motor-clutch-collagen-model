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
| Contact exclusion | Maximum bead penetration remains below `0.1 um` in every run |
| Contact action--reaction | Bead contact force plus cell reaction is numerical zero |
| Contact ablation | Disabling contact reproduces the preserved overlap failure |
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

The 20-seed, 600-s G3B calibration has been run. Its aligned, rotated, feedback-off and
no-fibre controls support the intended collagen-guidance mechanism, but 7/20 isotropic runs
failed the rigid-cell overlap guard. The complete 100-seed G3B/G3C final ensemble has not been
opened, and G3C formal validation is halted until the contact-enabled G3B campaign is rerun from
clean checkpoints and all calibration gates pass.

# G4 v4 validation status

G4 v4 is cumulative and does not overwrite G4 v3. The full generated build uses
a 7,200-s fixed-cell clock, a 21,600-s released-cell clock, 0.05-s solver steps,
four-minute fixed-cell output and ten-minute released-cell output.

| Gate | Full-build status |
|---|---|
| Passive fixed cell has zero traction and negligible alignment drift | **Pass** |
| Biological cell radius is constant | **Pass** |
| Force decomposition and cell/ECM action–reaction identities | **Pass** |
| Timestamped rupture and complete shared-load site failure | **Pass** |
| Near-cell alignment differs from passive zero across seeds 41–45 | **Pass**: mean \(\Delta S_r=0.0227\), 95% CI [0.0084, 0.0370] |
| One-hop displacement exceeds unconnected displacement | **Not yet**: mean difference \(1.40\times10^{-2}\,\mu\mathrm m\), 95% CI [\(-8.22\times10^{-4}\), \(2.88\times10^{-2}\)] |
| 8/10/12-µm radius sensitivity | **Implemented, not yet run as a full ensemble** |
| 180/270/360-µm constant-density domain convergence | **Implemented, not yet run** |
| Full 0.05/0.025-s timestep sensitivity | **Not yet run** |

Because the indirect-transmission interval still includes zero, the G4 v4-D
released-cell trajectory is a diagnostic output, not an accepted migration
prediction. The website preserves this failed gate rather than visually tuning
it away.
