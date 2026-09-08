# G4 v3 decision log

This file prevents a visual impression from silently becoming a model claim.
Every cumulative block records its new assumption, predicted consequence,
measurement and allowed interpretation.

## Block 0 — random geometry check

- **New assumption:** 240 straight source fibres have independent angles drawn
  uniformly from \([0,\pi)\); a circular cell-sized void is cut afterward.
- **Why needed:** radial alignment must be an output rather than a property of
  hand-placed contact fibres.
- **Prediction:** the finite realization should have low global nematic order,
  although individual shells can have sampling bias.
- **Result:** for seed 41, global nematic order was 0.043 and initial radial
  orders were -0.058, 0.074 and -0.008 in the three shells.
- **Decision:** pass for the demonstration geometry. Always use
  \(\Delta S_r\), not final \(S_r\), and report each seed's initial values.

## Benchmark — prescribed contraction with full linking

- **New assumption:** every geometric intersection is a permanent elastic
  hinge and the cell-sized inner boundary contracts by 5%, 10% or 20%.
- **Why needed:** this isolates whether stretching, bending and crosslink force
  transmission can reproduce the established bending-to-stretching recruitment
  mechanism before active-force calibration.
- **Prediction:** sufficiently strong contraction should raise near-cell radial
  order and stretching contribution more than far-field order.
- **Result:** at 20% contraction, five seeds gave positive near-cell
  \(\Delta S_r\) with mean 0.101; middle and far means were 0.012 and 0.002.
- **Decision:** pass as a mechanics benchmark. It does not establish a
  biological crosslink density or a cell-generated force.

## v3-A — probabilistic permanent links

- **New assumption:** every candidate intersection has a fixed random mark and
  is linked when that mark is below \(p_x\).
- **Why needed:** not every 2D intersection should transmit force, and the
  nested construction changes connectivity without changing fibre geometry.
- **Prediction:** increasing \(p_x\) should increase network connectivity,
  recruitment and indirect deformation.
- **Result:** contact-connected fraction rose from 0.037 at \(p_x=0\) to 0.983
  at \(p_x=0.60\). At 20% prescribed contraction, near-cell \(\Delta S_r\)
  rose from 0.014 to 0.054 over the same range.
- **Decision:** pass as an effective two-dimensional connectivity experiment.
  \(p_x\) is not a biochemical crosslinker concentration.

## v3-B — passive fixed cell

- **New assumption:** the cell is a fixed steric circle and applies no active
  traction.
- **Why needed:** it is the negative control for false alignment caused by cell
  insertion or numerical overlap correction.
- **Prediction:** displacement, traction and \(\Delta S_r\) should be zero.
- **Result:** maximum displacement was zero and all shell-wise
  \(\Delta S_r\) values were zero; residual repulsion was numerical roundoff.
- **Decision:** pass.

## v3-C — deterministic fixed-material-point traction

- **New assumption:** only the nearest continuous material point on each fibre
  inside a 3 um contact band is eligible; normalized Gaussian weights with
  \(\sigma_c=1.5\) um distribute a fixed total pull. The same material points
  remain attached throughout the run.
- **Why needed:** this tests active force transmission before attachment
  stochasticity or a moving contact rule is introduced.
- **Prediction:** force should be conserved exactly, contact fibres should move
  first, and crosslinked indirect fibres may respond through elastic paths.
- **Result:** 12, 24 and 48 nN at \(p_x=0.35\) gave near-cell
  \(\Delta S_r=0.0078,0.0102,0.0141\), with negligible middle/far alignment.
  At 24 nN, one-hop mean displacement was 0.0014 um and two-plus-hop mean
  displacement was 0.000041 um. Action-reaction error was zero.
- **Numerical check:** halving \(\Delta t\) changed the near-cell result by only
  \(7.1\times10^{-6}\) relative. Doubling duration to 960 s raised the
  near-cell result to 0.0139 but did not create middle/far alignment.
- **Decision:** mechanics implementation passes conservation and time-step
  checks, but the biological alignment gate does not pass. Do not add clutch
  failure or contact relocation yet.

## Next decision

The next cumulative block must use one selected experimental comparison to
decide among displacement-controlled contraction, nonspherical/contractile-
dipole cell geometry, or material/contact calibration. Its new assumption,
expected bias and acceptance metric must be written before it is implemented.
