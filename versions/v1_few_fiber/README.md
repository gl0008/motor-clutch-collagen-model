# Version 1 — Few-fiber overdamped baseline

## Research question

Can a fixed circular cell reorganize a few connected collagen-like fibers under
elastic stretching, bending, active pulling, and viscous drag alone?

## What changed from V0

- Removed SLS from every collagen bond.
- Replaced a triangular bulk substrate with five explicit bead chains.
- Kept the cell fixed and prescribed a total inward pull.
- Displayed each bead chain as a connected polyline.

## Equations

Each free bead obeys

\[
\zeta\dot{\mathbf r}_i=
\mathbf F_i^{stretch}+\mathbf F_i^{bend}+\mathbf F_i^{active}.
\]

Stretching uses Hookean bonds; bending penalizes the discrete curvature
`r[i-1] - 2 r[i] + r[i+1]`.  Drag dissipates energy at
`P_drag = sum_i zeta |v_i|^2`.  The elastic fibers store energy; drag is the
component that dissipates it.

## Important limitation

This archived version gives every bead a nonzero Gaussian weight.  It has no
physical micrometre calibration and no crosslinks, so it cannot represent
force transmission from a directly contacted fiber into a remote fiber.  V2
fixes both problems.

## Verification

- Fixed endpoints remain fixed.
- The normalized weights sum to one.
- The saved run increased local radial alignment from 0.310 to 0.501.

[Open the V1 animation](demo/index.html)

