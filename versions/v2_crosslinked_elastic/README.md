# Version 2 — Crosslinked elastic network, fixed cell

## Research question

After correcting the cell-to-collagen scale, are stretching, bending, permanent
crosslink force transmission and overdamped drag sufficient to produce local
fiber displacement or alignment around a fixed tumor cell?

## What changed from V1

- Geometry is in micrometres: a 9 µm cell radius in a 100 × 100 µm field;
  fibers are 24–78 µm long and have a displayed effective diameter of 0.30 µm.
- Fibers are finite, gently curved, connected bead chains (0.75 µm bead spacing).
- Geometric fiber intersections become permanent, freely hinged crosslinks.
- Direct cell force uses a **contact-shell + Gaussian hybrid**, not a Gaussian
  over the whole matrix.
- A soft rigid-cell repulsion prevents fibers from passing through the cell.

## Governing equation

For each free bead,

\[
\zeta\dot{\mathbf r}_i=\mathbf F_i^{stretch}
+\mathbf F_i^{bend}+\mathbf F_i^{crosslink}
+\mathbf F_i^{repulsion}+\mathbf F_i^{active}.
\]

There is no SLS state.  Stretching, bending, and crosslinks store elastic
energy.  The local drag term dissipates mechanical power:

\[
P_{drag}=\sum_i\zeta|\dot{\mathbf r}_i|^2.
\]

## Contact-shell + Gaussian implementation

For every segment `a–b`, the closest material coordinate to the cell center is

\[
\lambda=\mathrm{clip}\left(
\frac{(\mathbf r_c-\mathbf a)\cdot(\mathbf b-\mathbf a)}
{|\mathbf b-\mathbf a|^2},0,1\right),\quad
\mathbf p=(1-\lambda)\mathbf a+\lambda\mathbf b.
\]

Only a closest segment satisfying `0 ≤ |p-r_c|-R ≤ 3 µm` is kept per fiber and
only if it lies inside the ±30° protrusion sector.  Within that selected set,

\[
w_j=\frac{\exp[-(d_j/1.5\,\mu m)^2]}
{\sum_k\exp[-(d_k/1.5\,\mu m)^2]}.
\]

Outside the shell, direct active force is exactly zero.  Remote fibers move
only through crosslinks.  A material-point force is distributed to its two
segment endpoints by `(1-λ, λ)`, reducing bead-resolution bias.

## Assumptions to discuss with the professor

- Force, stiffness and drag are mechanism-first and not yet calibrated to a
  particular cell line or collagen concentration.
- The two distant ends of every finite fiber are fixed numerical far-field
  anchors.  This implements the requested fading-to-fixed far field, but a
  larger percolated domain should later test boundary sensitivity.
- Twelve fibers are deliberately initialized near the cell, including four in
  the right protrusion sector.  Remaining fibers are random and rejected if
  they cross the cell.
- Crosslinks never break, slide, or prefer an angle in V2.

## Verification and interpretation

- V2's saved run has 30 fibers and 128 permanent crosslinks.
- The hard-shell test finds at least three contact fibers and normalized weights.
- With zero active force the initial network is stress-free and does not move.
- Fixed endpoints remain fixed.
- Alignment uses `S_r = <2(t·e_r)^2-1>`: +1 radial, −1 tangential, 0 mixed.

The default run changes near-cell `S_r` from −0.119 to −0.103 and reaches
0.105 µm RMS displacement.  This small but positive change is intentionally
reported without retuning.  Use the no-crosslink control to test whether remote
deformation actually depends on network coupling.

[Open the V2 interactive lab](demo/index.html)

## Repository record

- **Lineage:** G1 V1 → **G1 V2** → G1 V3; corrected later by G2 V2.
- **Stable tag:** `g1-v2`; the corrected network lives in a different folder.
- **Notebook:** [purpose, complete equations, evidence boundaries and HOWEVER](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g1-v2).
- **References used:** R1 Lee; R2 Abhilash; R3 Wang; R4 Notbohm. See the
  [version-to-reference registry](../../references/README.md).
