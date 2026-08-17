# G3 model status and evidence boundary

## Scientific claim

G3 is a **2D minimal mechanism model of emergent cell–collagen guidance**. It tests whether
spatial motor-clutch attachments, local collagen geometry, and traction-dependent protrusion
persistence can generate a migration axis without a prescribed global polarity.

It is not a realistic three-dimensional tumor-migration model. It does not contain calibrated
collagen concentration/pore topology, a deformable cell or nucleus, proteolysis, EMT, cell-cell
interactions, SLS, transient crosslinks, or irreversible plasticity.

## Active stages

1. **G3A — material-point clutches.** A fixed rigid cell pulls a few elastic fibres. Each
   clutch stores `(fiber_id, segment_id, alpha)` and remains attached to that material point
   until Bell-law unbinding. A local Gaussian projects the point force to nearby beads on the
   same fibre while preserving net force and first moment.
   The single-fibre gate uses a tangential fibre with both endpoints fixed; the default run is
   15 s because a longer prescribed pull reaches the explicitly monitored no-steric-overlap
   boundary of this minimal model.
2. **G3B — emergent protrusions.** Twenty-four candidate surface sectors probe the matrix.
   Formation depends on local collagen availability/alignment; persistence is reinforced by
   bound-clutch fraction and traction success. Initial sector selection is uniform and contains
   no `+x` preference.
3. **G3C — rigid-body motion.** The cell translates and rotates only from equal-and-opposite
   clutch reactions. No prescribed self-propulsion velocity is present.

## Current validation state (2026-08-18)

G3B calibration has been completed for all six preregistered conditions using seeds 0--19
and 600 s per run (120 runs total). The configuration is **not frozen** and validation seeds
1000--1099 have **not** been opened.

The calibration supports the intended aligned-ECM mechanism:

- `aligned_8`: mean nematic guidance 0.726 (bootstrap 95% CI 0.626--0.819), with a 55/45
  split between the two directions of the collagen axis.
- `aligned_8_rotated_30`: the inferred protrusion axis rotated with the fixture; covariance
  error was 3.03 degrees.
- `aligned_feedback_off`: mean nematic guidance was -0.033 (95% CI -0.195--0.127), so
  removing geometry/traction feedback eliminated the aligned preference.
- `no_fibre`: zero bound clutches, force, and torque in all 20 runs.

G3B nevertheless **fails its calibration gate** because 7 of 20 `isotropic_random_8` runs
entered `invalid_geometry_overlap`: collagen beads were actively pulled more than 0.1 um
inside the rigid cell. The valid-run subset must not be treated as an unbiased isotropic
ensemble. At 20 seeds, the balanced and valid-isotropic residual polar ratios (0.182 and
0.276) also cannot establish the preregistered <0.1 final-ensemble criterion.

A concurrently launched `--stage both` job produced 169 partial G3C calibration records
before this G3B failure was summarized. Those records are retained for audit only and are
excluded from formal evidence. The job was stopped; G3C calibration/final validation must
not resume until the G3B geometry-overlap mechanism is resolved and G3B calibration passes.

See `docs/results/g3b_calibration_2026-08-18.md` and
`results/g3_validation/calibration/` for the complete summary and per-seed records.

## Relationship to Generation 2

The existing `generations/g2_corrected/` implementation remains frozen and unchanged.

- V2 supports force transmission through permanent crosslinks in its specified 2D network.
- V3 supports reaction-driven motion under a prescribed persistent clutch imbalance.
- Neither result establishes emergent direction selection or realistic 3D migration.

G3 is a new branch-level extension in the same repository. It does not rewrite G2 files or
reinterpret G2 V2/V3 as realistic tumor migration.

## Active/superseded map

| Path | Status |
|---|---|
| `src/g3/` | Active G3 implementation |
| `src/config/params_g3.yaml` | Active G3 configuration |
| `generations/g3_emergent_guidance/` | G3 scientific-version entry and run instructions |
| `generations/g2_corrected/` | Frozen G2 mechanism demonstrations |
| `versions/` | Frozen Generation 1 archive |

## Interpretation rule

All generated figures must include:

> Personal simulation / mechanism validation; not yet a realistic 3D tumor-migration prediction.
