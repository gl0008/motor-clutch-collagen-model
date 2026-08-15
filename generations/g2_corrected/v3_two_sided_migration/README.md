# Generation 2 · V3 — two-sided clutch imbalance and cell motion

## What changes from corrected V2?

The same connected bead--spring collagen network is used.  The fixed cell now
has a left and a right clutch ensemble (12 effective clutches on each side).
The synchronized comparison changes one constraint only:

- **fixed:** the rigid cell centre is held at its initial position;
- **moving:** the centre follows overdamped translation along the declared
  left/right protrusion axis under the reaction force from collagen.

Crosslinks remain permanent.  V3 does not use crosslink rupture to manufacture
motion.

## Motor--clutch law

For bound clutch \(j\),

\[
F_j=k_c x_j, \qquad
k_{off,j}=k_{off}^0\exp(F_j/F_b),
\]

and the actin loading speed on each side is

\[
v_a=v_0\max(0,1-F_{side}/F_{stall}).
\]

A mild, declared polarity \(\psi_{pol}=0.65\) increases the right-side binding
rate and decreases the left-side rate.  This follows the need for persistent
polarity in the 1-D cell-migration framework; random binding alone gives only a
wandering, zero-mean imbalance.

The collagen receives the two Gaussian-distributed side tractions.  Its
opposite reaction drives the released rigid cell:

\[
\zeta_c\dot{x}_c=-\sum_i F_{i,x}^{active},\qquad \dot y_c=0.
\]

This one-axis constraint is intentional.  Full 2-D translation and rotation
would add a new degree of freedom before the left/right causal experiment has
been validated.

Clutch attachment sites are recomputed from the nearest eligible material
points every 0.5 s as both cell and fibres move.  Fixed and moving runs use the
same generated ECM and the same indexed random stream.

## Calibration boundary

The animation lasts 30 min and stores frames every 6 s.  Cell mobility is
checked over 20 fast clutch-ensemble seeds against an MDA-MB-231 target range of
0.2--0.4 µm/min; the full animation then includes collagen deformation.  This
is a mechanism-first calibration, not a fit to one trajectory.

The visualization draws the initial cell outline and an accumulated trajectory,
so a several-micrometre translation cannot be hidden by a 180 µm field of view.
Actual bound clutches appear as spokes between the cell surface and the current
left/right collagen contact points.

## Repository record

- **Lineage:** G2 V2 → **G2 V3** → G2 V4.
- **Stable tag:** `g2-v3`; the G1 V3 prototype remains separately accessible.
- **Notebook:** [purpose, complete equations, evidence boundaries and current limitation](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g2-v3).
- **References used:** R5 Bell; R6 Bangasser & Odde; R7 Prahl; R8
  Steinwachs; R9 Aguilar-Rojas. See the
  [version-to-reference registry](../../../references/README.md).
