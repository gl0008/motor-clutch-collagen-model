# D011 -- Conservative cell--collagen contact

## Status

Implemented on 2026-08-18. Calibration is open; final-validation seeds remain sealed.

## Problem exposed by the contact-free calibration

Seven of twenty `isotropic_random_8` G3B calibration runs pulled at least one collagen bead
more than 0.1 um inside the nominally rigid cell. The code correctly labeled those runs
`invalid_geometry_overlap`. The rigid-cell assumption was therefore incomplete: clutch
traction existed, but excluded-volume mechanics did not.

## Decision

Use a one-sided harmonic penalty potential for every collagen bead:

\[
d_i=\|\mathbf r_i-\mathbf r_c\|,
\qquad
\delta_i=\max(0,R-d_i),
\]

\[
U_i^{\rm contact}=\frac12 k_{\rm contact}\delta_i^2,
\qquad
\mathbf F_i^{\rm contact}=k_{\rm contact}\delta_i
\frac{\mathbf r_i-\mathbf r_c}{d_i}.
\]

The cell reaction is `-sum(F_contact)`. The contact force is central and therefore contributes
zero ideal torque about a circular cell center. It is included in the same global force and
moment conservation diagnostics as clutch traction.

The functional form is adapted from the repulsive signed-distance branch of the linear
elastic contact model in Runser, Vetter & Iber, *Nature Computational Science* 4, 299--309
(2024), Methods Eq. 2: https://doi.org/10.1038/s43588-024-00620-9.

## Active assumptions

- `contact_stiffness = 4.0e-3 N/m`, provisionally equal to the active collagen extensional
  spring scale rather than fitted to a measured cortex modulus.
- Collagen beads are zero-radius material nodes; contact begins at the rigid-cell radius.
- Contact is repulsive only. There is no contact adhesion, tangential friction, cortex
  deformation, fibre thickness, or nucleus.
- The explicit ECM stability check includes the added contact stiffness.
- The original `invalid_geometry_overlap` guard remains active as a safety gate. Contact does
  not convert a deeply overlapping trajectory into an accepted one.

## Verification completed before ensemble rerun

- force points outward and equals the negative gradient of contact energy;
- bead force plus cell reaction is numerical zero;
- global contact torque is numerical zero for a circular cell;
- Numba and NumPy implementations agree;
- disabling contact reproduces isotropic seed 2 failure at 48.1 s;
- enabling contact lets the same seed complete 60 s with maximum penetration 0.00375 um.

## Required next gate

Run all G3B conditions with calibration seeds 0--19 for 600 s in a new output root. Do not
reuse the contact-free checkpoints. Report contact count, maximum penetration, contact energy,
the original guidance/symmetry metrics, and a contact-stiffness sensitivity control. Freeze
the configuration only if all geometry and preregistered G3B gates pass.
