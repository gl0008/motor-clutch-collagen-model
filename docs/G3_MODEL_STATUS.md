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
