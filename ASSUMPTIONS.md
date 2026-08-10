# Assumption and evidence register

| ID | Version | Choice | Evidence status | Required follow-up |
|---|---|---|---|---|
| A1 | V1–V4 | Overdamped bead dynamics; no inertia | Professor's requested baseline. Drag dissipates energy but is not SLS material memory. | timestep/drag convergence and experimental timescale fit |
| A2 | V2–V4 | 2D 100 µm domain; cell radius 9 µm | representative cell-scale geometry, not cell-line-specific | replace using microscopy for the selected tumor/immune cell |
| A3 | V2–V4 | finite 24–78 µm fibers, 0.30 µm effective diameter | inside Lee et al.'s reported 20–200 µm fiber and 0.20–0.35 µm thickness ranges | image-derived length/orientation/pore distributions |
| A4 | V2–V4 | 0.75 µm bead spacing | numerical resolution choice | repeat at 0.50 and 1.0 µm with segment-level observables |
| A5 | V2–V4 | initialize 12 near-cell fibers, four per required side sector | avoids silently coupling to remote collagen; deliberate initialization | sample image-derived networks and report failed contact realizations |
| A6 | V2–V4 | 3 µm hard contact shell then 1.5 µm Gaussian | professor-approved hybrid; numerical contact width not fitted | sweep contact width and sigma separately |
| A7 | V2–V4 | one closest material segment per fiber contact patch | reduces bead-density bias | compare multi-patch clustering and contact-area data |
| A8 | V2–V4 | two distant ends of each finite fiber fixed | requested far-field anchoring; stronger than a true percolating bulk boundary | enlarge domain and test boundary sensitivity |
| A9 | V2–V4 | intersection crosslinks are permanent, freely hinged | elastic baseline inspired by explicit crosslink bead models | compare crosslink density/stiffness to shear data |
| A10 | V2–V4 | linear stretching and bending | mechanism-isolation baseline | add buckling and strain stiffening after baseline validation |
| A11 | V3 | 12 effective clutches per side | coarse adhesion units, not molecules | calibrate to adhesion/traction measurements |
| A12 | V3 | Bell slip-bond off-rate and linear motor force–velocity law | standard minimal motor–clutch starting laws | test catch bonds/glassy kinetics only if data require |
| A13 | V3 | rigid translating cell; no rotation/deformation | isolates reaction-force-driven motion | add torque, rotation, then a deformable boundary |
| A14 | V3 | same counter-addressed random stream in both conditions | variance-reduction design | ensemble over seeds and report uncertainty |
| A15 | V4 | 50% dynamic links; force-accelerated rupture and same-neighborhood re-formation | new mechanism-first hypothesis, not literature-calibrated | sweep kinetics and fit load–unload/residual-alignment data |
| A16 | V4 | new rest vector equals geometry at re-formation | minimal plastic rest-state reset | replace with evidence-based sliding, merging, or new-partner search |

The default bead drag is `ζ = 1` model force·s/µm.  It was selected so that
crosslink-mediated propagation is resolvable during the 6–8 s demonstration
window; it is not a measured collagen drag.  A diagnostic run at `ζ = 8`
confined nearly all visible motion to direct contacts.  Drag and observation
time must therefore be calibrated together before physical times are claimed.

## Required validation order

1. Verify elastic forces, fixed boundaries, zero-force equilibrium, and timestep convergence.
2. Fit network mechanics to collagen shear/tensile or microrheology data.
3. Verify V2 spatial decay: direct force is zero outside contact; remote motion
   vanishes or decreases strongly in the no-crosslink control.
4. Run a V3 seed ensemble; report displacement, persistence, reaction force,
   clutch occupancy, and left/right imbalance distributions.
5. Only then fit V4 using load–unload data and quantify residual alignment and
   densification after force removal.
