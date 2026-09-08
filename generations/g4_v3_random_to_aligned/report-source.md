# Random-to-aligned collagen: evidence and G4 v3 execution report

**Audience:** Gloria Liu and her advisor

**Date:** 2026-09-08

**Scope:** Whether initially random collagen can become radially aligned, and
how to test that mechanism cumulatively without altering frozen G4 versions.

## Direct answer

There is solid experimental and computational precedent that an initially
disordered fibrous collagen network can become aligned under cell-scale strain
or contractility. What is not established is that the exact proposed
combination—two-dimensional bead fibres, the present effective force scale,
stochastic clutch clusters and repeated nearest-material-point relocation—will
reproduce that transition. G4 v3 therefore validates the elastic network first
and stops before clutch kinetics when deterministic traction remains too weak
and local.

## Consequential evidence

- Abhilash et al. model contractile cells in random linked fibre networks and
  report local alignment, stretching and shape-dependent long-range force
  transmission. This supports the prescribed-contraction mechanics benchmark,
  not the present motor-force calibration.
- Vader et al. experimentally observe strain-induced alignment in collagen.
  Crosslinked samples align reversibly, showing that plasticity is not required
  to test alignment itself.
- Provenzano et al. observe tumor explants realigning initially random collagen
  and invading along radial fibres. This supports radial alignment as a relevant
  output, not a particular causal equation.
- Kopanska et al. report spheroid-driven collagen contraction, surface-parallel
  and radial organization, and collagen retraction after ablation. This supports
  a tension-bearing network interpretation over hour-to-day experiments.
- Böhringer et al. compare initially disordered network simulations with 3D
  collagen experiments, measuring orientation in distance shells and finding a
  relationship between local contractility and fibre alignment/density.
- Ley et al. (2026) provide a closely related 3D bead-spring, tractor-based,
  cell–ECM-feedback model. Its radial reorientation result includes repeated
  cycles and plasticity, so it is recorded as a later-stage comparison rather
  than imported into the present elastic baseline.

## Implemented comparison

1. Generate straight finite fibres with uniform random angles and cut the same
   cell-sized void for every cumulative block.
2. Verify initial global and shell-wise orientation rather than assuming the
   realised finite network is exactly isotropic.
3. Run a prescribed-contraction benchmark with full permanent linking.
4. Replace full linking with nested probabilistic permanent links while keeping
   geometry and random marks fixed.
5. Insert a passive fixed cell as a zero-force control.
6. Apply deterministic, force-conserving traction to only the nearest material
   point on each fibre inside the hard contact band, with Gaussian weighting
   within that eligible set.
7. Stop before stochastic clutch failure or contact relocation unless the
   deterministic force-controlled alignment gate passes.

## Result and limitation

The fully linked 20% contraction benchmark produced reproducible but localized
near-cell radial alignment in five seeds. The inherited 12–48 nN deterministic
traction range at effective crosslink probability 0.35 produced only weak
near-cell alignment and negligible middle/far response. Halving the numerical
time step did not change that conclusion. Doubling duration increased the local
response modestly but did not establish a long-range effect or a two-hour
steady state.

The next experiment must be source-matched: measured boundary displacement,
documented nonspherical/contractile-dipole geometry, or calibration to one
collagen preparation. These alternatives must not be mixed into a single
uninformative tuning exercise.

## Claim-to-source record

Full titles, authors, DOI links, open-access links and the exact imported claim
are recorded in [REFERENCES.md](REFERENCES.md). All seven foundation papers were
located with open full text. The 2026 Ley et al. article metadata and abstract
were accessible; full-text availability was not required because no mechanism
from it was implemented in this block.

## Research stopping reason

The key phenomenon, metrics and alternative explanations have primary support;
the remaining uncertainty is parameter/experiment matching rather than a lack
of additional papers. Further broad searching is unlikely to decide among the
three next calibration experiments without first choosing the target biological
system and measurement.
