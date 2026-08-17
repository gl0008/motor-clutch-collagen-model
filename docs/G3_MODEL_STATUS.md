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

## External G2 reference

Gloria's G2 code is frozen and external to this implementation.

- V2 supports force transmission through permanent crosslinks in its specified 2D network.
- V3 supports reaction-driven motion under a prescribed persistent clutch imbalance.
- Neither result establishes emergent direction selection or realistic 3D migration.

The nested checkout must not be modified or committed into this parent repository.

## Active/superseded map

| Path | Status |
|---|---|
| `src/g3/` | Active G3 implementation |
| `src/ecm/forces.py` | Active reusable elastic force laws |
| `src/simulation/integrator.py` | Active reusable overdamped integrator |
| `src/cells/traction.py` | Superseded prescribed dipole |
| `src/cells/motility.py` | Superseded direct polarity torque and self-propulsion |
| `src/simulation/cell_ecm.py` | Superseded Phase-B coupled loop |
| `versions/v1/` | Frozen archive |
| `test/v1_5/` | Archived SLS experiments |

## Interpretation rule

All generated figures must include:

> Personal simulation / mechanism validation; not yet a realistic 3D tumor-migration prediction.
