# Version 3 — Two-sided clutch imbalance and cell translation

## Research question

Can stochastic left/right motor–clutch dynamics create a persistent force
imbalance, and does that imbalance move a rigid cell when the cell is released?

## What changed from V2

- Added 12 effective clutches to each of two opposing ±30° sectors.
- Each clutch binds a V2 contact patch and retains that segment's material
  coordinate until it unbinds; contact candidates are recomputed as fibers move.
- Added Bell slip-bond rupture and a linear motor force–velocity relation.
- Runs two synchronized conditions from the same collagen network and
  counter-addressed random stream: fixed cell and translating cell.
- Crosslinks remain permanent.  Cell motion does **not** require collagen
  crosslink rupture in this baseline.

## Clutch and motor equations

For a bound clutch,

\[
F_c=k_c x_c,\qquad
k_{off}=k_{off}^0\exp(|F_c|/F_b).
\]

For each side, actin retrograde speed slows with total clutch force:

\[
v_a=v_0\max\left(0,1-\frac{\sum F_c}{F_{stall}}\right).
\]

Fiber material-point motion relieves clutch loading.  When the cell is released,

\[
\zeta_c\dot{\mathbf r}_c=\mathbf F_{reaction},\qquad
\mathbf F_{reaction}=-\sum_i\mathbf F_i^{active}.
\]

In the fixed condition `r_c` is held exactly constant and the same reaction is
recorded as the force the constraint must supply.

## Assumptions to discuss

- The 12 clutches per side are effective adhesion units, not individual
  integrin molecules.
- Binding/off rates, clutch stiffness, motor speed, stall force, and cell drag
  are mechanism-first.  They require calibration before a predicted migration
  speed can be interpreted biologically.
- Only translation is allowed; the cell cannot rotate or deform.
- Bell slip bonds are the simplest starting law.  Catch-bond or glassy-clutch
  kinetics should be added only if data require them.

## Controlled comparison

Both simulations use the same seed and addressed random number
`u(seed, time step, side, clutch, event channel)`.  Before mechanics causes the
hazards to diverge, bind/unbind proposals are identical.  Therefore the only
intentional difference is whether the cell force balance is integrated.

[Open the synchronized V3 lab](demo/index.html)

