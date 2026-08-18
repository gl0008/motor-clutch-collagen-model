# G3 model status and evidence boundary

## Scientific claim

G3 is a **2D minimal mechanism model testing cell–collagen guidance**. It tests whether
spatial motor-clutch attachments, local collagen geometry, and traction-dependent protrusion
persistence can generate a migration axis without a prescribed global polarity.

It is not a realistic three-dimensional tumor-migration model. It does not contain calibrated
collagen concentration/pore topology, a deformable cell or nucleus, proteolysis, EMT, cell-cell
interactions, SLS, transient crosslinks, or irreversible plasticity.

## Active stages

1. **G3A-R — protrusion-tip material-point clutches.** A fixed rigid cell grows a visible
   protrusion into a crosslinked 99-fibre network. Each
   clutch stores `(fiber_id, segment_id, alpha)` and remains attached to that material point
   until Bell-law unbinding. A local Gaussian projects the point force to nearby beads on the
   same fibre while preserving net force and first moment.
   The displayed and freshly validated condition uses the shared crosslinked 99-fibre fixture,
   with boundary-only anchors and eight controlled near-cell contact fibres. The diagnostic
   protrusion is prescribed only in G3A so growth, capture and network loading can be isolated.
2. **G3B-R — cell-intrinsic polarity.** A conserved noisy activity field on 24 membrane
   sectors breaks symmetry without collagen-direction input. Physical attachment and traction
   can stabilize activity after contact; there is no `+x` preference.
3. **G3C — rigid-body motion.** The cell translates and rotates only from equal-and-opposite
   clutch reactions. No prescribed self-propulsion velocity is present.

## Current validation state (2026-08-18)

### G3-R 99-fibre revision

A fresh short mechanism suite has now been run for the explicit-protrusion, cell-intrinsic
polarity revision. All displayed G3A/B/C conditions use a crosslinked 99-fibre network rather
than the old few-fibre fixtures. All nine mechanism/execution gates passed, including matched
motor-off, fixed-cell and empty-network controls. The full result table and interpretation are
in `docs/results/g3_revision_validation_2026-08-18.md`.

This revision does not erase the preserved result below and does not yet clear a guidance
claim. Its four-seed aligned ensemble has nematic order 0.0142, whereas the four-seed isotropic
ensemble has 0.0160. A matched 20 s feedback OFF/ON pair also follows the same sector sequence
and finishes at 30.02 degrees. The new model therefore needs a larger preregistered aligned-versus-
isotropic calibration, rotated-axis covariance and parameter sensitivity before final
validation seeds are opened.

### Superseded original contact-free calibration — audit only

The original contact-free G3B calibration was completed for all six preregistered conditions
using seeds 0--19 and 600 s per run (120 runs total). It is preserved for audit but is not
current evidence because one comparison condition failed its geometry-validity gate. The
configuration is **not frozen** and validation seeds 1000--1099 have **not** been opened.

The valid subsets produced the following historical signals, which motivated the revision but
must not be reported as a successful guidance result:

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

A concurrently launched `--stage both` job produced 173 partial G3C calibration records
before this G3B failure was summarized. Those records are retained for audit only and are
excluded from formal evidence. The old campaign was stopped. The fresh G3C-R short mechanism
run is a separate control-based execution test; predictive migration validation must not
resume until the revised G3B guidance gate passes.
`results/g3_validation/HALT.json` enforces this boundary in the campaign runner so stale
resume commands cannot silently continue the campaign.

The active configuration now contains an explicit conservative cell--collagen
contact law. A bead that penetrates the rigid circle receives a radial force derived from
`0.5 * contact_stiffness * penetration^2`, and the cell receives the exact opposite reaction.
The previously failing isotropic seed 2 reproduces its 48.1 s overlap when contact is disabled;
with contact enabled it completes a 60 s diagnostic with 0.00375 um maximum penetration.
This was an implementation diagnostic only. The new 99-fibre short suite has since been run,
but its aligned/isotropic results remain negative. A larger preregistered rerun in a new
version-preserving output directory is required before the guidance stage can be cleared.

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
