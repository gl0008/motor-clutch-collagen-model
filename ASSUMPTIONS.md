# Assumption and evidence register

This file separates facts transferred from papers from choices made because the original simulation
code and complete calibration data are unavailable.

| ID | Current choice | Status and evidence | How it will be tested or replaced |
|---|---|---|---|
| A1 | Use a 2D network and cell. | Modeling choice. Adebowale's CMS is a 2D migration simulator, while Zhao demonstrates that the same filament framework can be run in 2D and 3D. | Reproduce the 2D trends first; add a 3D network only after calibration. |
| A2 | Use a diluted triangular lattice with freely hinged crosslinks. | Zhao uses a bond-diluted hexagonal/triangular lattice in 2D; Licup and related collagen-network models use diluted subisostatic fiber networks. | Compare against image-derived or random off-lattice collagen geometry later. |
| A3 | Apply bending only to surviving collinear segment pairs. | Directly follows Zhao's definition that remaining collinear lattice bonds belong to the same filament and that crosslinks are freely hinged. | Unit-test the undeformed network and bending response. |
| A4 | Replace every intrafiber axial spring by one SLS element. | Requested hypothesis. Single collagen fibrils are intrinsically viscoelastic, but experiments often require multiple relaxation modes (Shen 2011; Chasiotis group 2022/2023). | Begin with one relaxation time so the causal effect is identifiable; generalize to a Prony series if one SLS cannot fit data. |
| A5 | Keep bending elasticity time-independent. | Simplifying choice. The user specifically requested replacement of bead-to-bead springs; no evidence currently identifies the correct bending-relaxation law. | Add viscoelastic bending only if bulk relaxation cannot be fitted by axial SLS bonds. |
| A6 | Do not allow fiber or crosslink breaking/rebinding. | Zhao neglects plasticity for reversible chromatin deformation. Collagen can be plastic through crosslink and entanglement dynamics, but that is a distinct mechanism. | Add dynamic crosslinks only in a later plastic-remodeling experiment. |
| A7 | Use overdamped bead dynamics without inertia, hydrodynamic coupling, or thermal noise. | Direct transfer from Zhao's minimal bead dynamics. It is appropriate for an athermal mesoscopic network prototype. | Check timestep convergence and, if needed, add Brownian noise or hydrodynamic mobility. |
| A8 | Clutches bind beads rather than a continuum substrate. | Required coupling choice. It is the discrete counterpart of the clutch–substrate link in Adebowale. | Later allow binding to any point on a segment and distribute force to its endpoint beads to remove mesh bias. |
| A9 | A simulated clutch is a coarse adhesion unit. | Necessary coarse-graining choice: Adebowale's CMS uses many molecular clutches, while this prototype uses tens of effective clutches. | Fit effective `kon`, `koff0`, `Fb`, and `kc` to bond-lifetime and traction data. |
| A10 | Default fast and slow SLS times are 5 s and 50 s. | Literature-informed starting values, not a fit. Isolated fibrils show fast relaxation around 2–7 s and slower modes around 60–160 s. | Replace with the exact Adebowale substrate values for reproduction, then with collagen relaxation data for the biological model. |
| A11 | Use linear axial SLS forces in tension and compression. | Deliberate reproduction baseline. Real collagen networks buckle easily in compression and stiffen through alignment/stretching. | After baseline verification, introduce a small compression/tension stiffness ratio and test whether nonlinear collagen behavior is required. |
| A12 | Default network stiffness is numerical, not yet mapped to 2 kPa. | A collagen fibril modulus cannot be equated with bulk gel modulus; network geometry and connectivity strongly renormalize the response. | Run virtual shear/step-strain tests and optimize segment parameters until network-level `E0`, `Einf`, and relaxation curve match target data. |
| A13 | The 1D fibre has 10 beads and nine identical serial SLS bonds, with per-bond stiffness multiplied by nine. | Numerical discretization choice. Axial segment stiffness follows `k = EA/L`; shortening a fixed fibre segment by a factor of nine increases its stiffness by nine (Lee 2014). | Confirm identical end-to-end step-relaxation responses with 5, 9, and 19 bonds. |
| A14 | The 1D internal beads are massless and remain at uniform axial strain. | Mechanism-isolation choice. It removes an additional bead-drag timescale while retaining SLS relaxation, force loading, and rupture. | Introduce consistently scaled bead drag only in the later 2D network stage and perform resolution convergence. |
| A15 | Use `Fb`, `kc`, and `v0` as force, stiffness, and velocity scales. | Motor–clutch nondimensionalization: `L*=Fb/kc` and `t_load=Fb/(kc*v0)`. Defaults compare `tau/t_load` and `Kchain/kc`, not fitted physical units. | Map the scales to pN, nm, and s only after clutch and collagen calibration. |
| A16 | A clutch transmits tension but not compression. | Adhesion-link modeling choice; a detached or slack molecular link should not push the collagen fibre. | Compare with a signed linear clutch as a numerical control. |
| A17 | Cluster episodes still active at the end of a run are right-censored. | Statistical requirement; treating the 30-`t_load` cutoff as a failure would bias lifetime downward. | Report Kaplan–Meier medians and observed-failure counts; extend duration if the median is not reached. |

## Parameter calibration order

1. Match single-bond analytical relaxation (implementation check).
2. Match the network's instantaneous bulk modulus by adjusting axial and bending stiffnesses.
3. Match the long-time/instantaneous modulus ratio by adjusting `kinf/k0`.
4. Match the network relaxation curve by adjusting `tau` (or adding Prony modes).
5. Fit clutch kinetics to bond lifetime and traction-force distributions.
6. Only then compare migration speed, MSD, and persistence between fast and slow conditions.

## Excluded from version 0.1

- condensate phase fields and all Cahn–Hilliard terms from Zhao;
- collagen degradation or synthesis;
- crosslink rupture, rebinding, or plastic rest-length remodeling;
- ligand-density heterogeneity;
- catch-bond integrin kinetics;
- steric contact between fibers;
- fluid poroelasticity;
- cell shape, nucleus, lamellipodia, and explicit filopodium length dynamics.
