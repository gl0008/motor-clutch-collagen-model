# G3 — Emergent cell–collagen guidance

**Status:** active 2D mechanism-validation model. Not a realistic 3D tumor-migration model.

G3 replaces the superseded prescribed dipole/direct-polarity cell path with spatial stochastic
clutches. Each clutch binds a collagen segment material coordinate, pulls through an evolving
vector, and projects its force locally to beads on the same fibre. Twenty-four candidate
surface sectors provide unbiased probes; local collagen availability/alignment affects
formation, while bound fraction and traction affect persistence. G3C releases the rigid cell
and applies only equal-and-opposite clutch reactions.

A one-sided harmonic contact potential now prevents collagen material nodes from crossing the
rigid circular boundary. Contact is conservative, frictionless, and non-adhesive; its exact
opposite reaction is included in G3C cell motion. The contact-enabled ensemble has not yet
passed calibration.

## Stage boundaries

- **G3A:** prescribed test protrusion, fixed cell, material-point and conservation gates.
- **G3B:** two active protrusions, fixed cell, emergent nematic-axis and ablation gates.
- **G3C:** mobile rigid cell, reaction-only translation/rotation and drag controls.

## Interpretation

The model can test whether directionality emerges from the specified coarse-grained feedback.
Absolute speed remains a calibrated output because cell drag is provisional. Fibre traversal,
pore confinement, nuclear deformation, proteolysis, 3D topology, rheology-calibrated collagen
concentration, stress relaxation, and plastic remodeling are deferred to later generations.
