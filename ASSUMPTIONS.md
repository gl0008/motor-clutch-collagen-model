# Assumption and evidence register

Generation 1 assumptions are preserved inside the legacy version folders and
the Git tag `legacy-generation-1-2026-08-14`.  This table governs the corrected
Generation 2 baseline.

| ID | Stage | Choice | Evidence / reason | Remaining follow-up |
|---|---|---|---|---|
| G2-A1 | V2–V4 | Overdamped bead force balance; no inertia and no SLS | professor's requested baseline; drag dissipates bead motion, while springs store energy | fit effective drag to time-resolved collagen displacement |
| G2-A2 | V2–V4 | MDA-MB-231 cell radius 9 µm | representative 18 µm suspended-cell diameter | replace with microscopy from the intended experiment |
| G2-A3 | V2–V4 | 180 × 180 µm field, 99 fibres | 140/180 µm short-run response differed by more than the declared 5% boundary tolerance, so the larger field was adopted and fibre count scaled by area | compare 180 with a still larger field before quantitative prediction |
| G2-A4 | V2–V4 | finite 20–80 µm fibres; 0.30 µm effective diameter | within Lee et al.'s reported fibre-level 20–200 µm and 0.20–0.35 µm ranges | use image-derived length, diameter, curvature and pore distributions |
| G2-A5 | V2–V4 | 0.75 µm bead spacing | computational resolution; 0.50/0.75/1.0 sensitivity is saved and is not fully converged for pointwise displacement | refine the discrete bending/contact implementation before parameter fitting |
| G2-A6 | V2–V4 | only beads inside the outer boundary band are fixed | approximates a finite observation window embedded in bulk; no finite fibre end is fixed merely because it is an end | compare fixed band, periodic boundary and larger-domain results |
| G2-A7 | V2–V4 | reject networks unless all contact fibres and ≥85% of all fibres reach a fixed boundary through links | prevents a visually connected but mechanically floating contact cluster | ensemble over image-derived network realizations |
| G2-A8 | V2–V4 | wet fibre modulus 32 MPa; `k=EA/l0`, bending `EI/l0^3` | Lee et al. scale and beam discretization | calibrate network shear/tensile response, not just single-fibre modulus |
| G2-A9 | V2–V4 | compression stiffness = 10% of tension | minimum microbuckling baseline; 0/10/100% sensitivity saved | use fibre-level compression/buckling data |
| G2-A10 | V2–V4 | permanent freely hinged links at intersections; 75 nN/µm penalty | explicit inter-fibre coupling; 400 kPa literature modulus does not uniquely map to a 2-D point spring | calibrate link density and stiffness to shear data |
| G2-A11 | V2–V4 | 3 µm surface shell followed by Gaussian `sigma=1.5 µm` | professor-approved hybrid; Gaussian is normalized only over eligible contacts | fit contact width to protrusion/adhesion microscopy |
| G2-A12 | V2/V4 | total pull 5 nN | mechanism-first value within a 2.5/5/10 nN sensitivity and cell-traction scale | fit to 3-D traction reconstruction for the selected condition |
| G2-A13 | V3 | 12 effective clutches per side; Bell slip bond; linear motor force–velocity | standard minimal motor–clutch laws | test catch bonds or glassy kinetics only if data require |
| G2-A14 | V3 | persistent polarity probability 0.65 | a symmetric stochastic pair has no persistent migration direction; Prahl et al. explicitly requires polarity for directional migration | derive polarity dynamics from cell data rather than prescribe it |
| G2-A15 | V3 | cell moves only on the left/right protrusion axis | isolates the requested two-sided imbalance; an early diagnostic 2-D run produced unsupported y drift and excessive strain | add full 2-D translation and torque as a later stage |
| G2-A16 | V3 | cell drag 300 nN·s/µm; 30 min run; frames every 6 s | 20-seed clutch ensemble gives median path speed 0.280 µm/min, inside the 0.2–0.4 MDA-MB-231 range | calibrate jointly with collagen concentration and full trajectory statistics |
| G2-A17 | V4 | new weak link requires different fibres, distance ≤0.45 µm, alignment ≤30°, near-cell radius ≤30 µm, and ≥0.25 nm approach from initial geometry | Ban et al. motivates stretch/approach-dependent weak-link formation; the approach test prevents pre-existing neighbours from being mislabeled | fit formation kinetics/capture to creep-recovery or imaging data |
| G2-A18 | V4 | new link is stress-free at formation and has 15 nN/µm stiffness | minimal Ban-style coarse-grained hypothesis; weaker than permanent links | sweep stiffness/rate only after a load regime creates measurable new contacts |
| G2-A19 | future | Kolade/Adebowale U937 outward-pushing model is separate | different cell type, force direction and adhesion dependence from the professor's MDA-MB-231 pulling baseline | implement as a separate generation/biology track |

## Generation 3 assumptions

| ID | Stage | Choice | Evidence / reason | Remaining follow-up |
|---|---|---|---|---|
| G3-A1 | G3A–C | Eight 40 µm fibres, 1 µm bead spacing, fixed endpoints | deliberately inspectable few-fibre fixtures | replace with 3D image/concentration-calibrated topology only in G4 |
| G3-A2 | G3A | Clutch binds one segment material coordinate until Bell unbinding | prevents unphysical nearest-bead jumping | validate against spatial adhesion tracking when available |
| G3-A3 | G3A–C | Engagement begins surface-normal; subsequent force follows the evolving clutch vector | preserves a physical moment arm and permits rigid-cell torque | compare with an explicit cortex/adhesion geometry later |
| G3-A4 | G3A–C | Same-fibre Gaussian force projection with first-moment correction | represent a point clutch on discrete beads while conserving force and torque | test sensitivity to projection width/resolution |
| G3-A5 | G3B–C | 24 unbiased sectors, at most two active protrusions | minimal 2D probing model with no prescribed `+x` polarity | validate protrusion counts and lifetimes from microscopy |
| G3-A6 | G3B–C | Geometry and traction reinforce protrusion persistence | Carey et al. 2016 supports the feedback direction | sweep all provisional rates/gains; do not treat them as measured constants |
| G3-A7 | G3C | Rigid circular cell, reaction-only translation/rotation | isolates mechanical guidance without hidden self-propulsion | deformable cell, nucleus and pore confinement are deferred |
| G3-A8 | G3C | 300 nN·s/µm translational drag and dimensional rotational drag | calibrated timescale closure, not an independent prediction | report direction/symmetry before speed; sweep drag |
| G3-A9 | all | Permanent elastic ECM; no SLS or irreversible remodeling | isolates guidance before material memory | report unresolved recovery as unresolved, not plasticity |
| G3-A10 | G3A–C | One-sided conservative bead–cell contact with `k_contact=4×10⁻³ N/m` | Runser et al. 2024 supports the signed-distance repulsive form; required after the preserved overlap failure | complete 2/4/8×10⁻³ N/m sensitivity and all-condition recalibration; add cortex/fibre geometry later |

Full values and provenance are in [`docs/parameter_provenance.md`](docs/parameter_provenance.md).

## Generation 4 assumptions

| ID | Stage | Choice | Evidence / reason | Remaining follow-up |
|---|---|---|---|---|
| G4-A1 | A–D | 99 finite fibres, 180 µm field, outer-boundary-only anchoring | inherit the accepted G2-scale geometry so G4 changes topology/contact rather than the whole field at once | use image-derived 3D geometry |
| G4-A2 | A–D | 3 MPa axial modulus and independently tunable bending multiplier | inherit the active G3 softened realignment fixture; professor requested lower bending | fit shear/tension and time-resolved fibre rotation |
| G4-A3 | A–D | retain an intersection when a fixed random mark is below `p_x` | avoids G2's every-intersection link while keeping nested, comparable topologies | infer crosslink density/topology from imaging or rheology |
| G4-A4 | A–D | retained links are permanent, linear and freely hinged | isolates elastic transmission before crosslink breaking/plasticity | add breaking only in a later generation |
| G4-A5 | A–D | hard 3 µm contact eligibility then Gaussian `sigma=1.5 µm` only on direct fibres | professor-approved hybrid; prevents whole-network Gaussian leakage | fit contact width to cell-surface imaging |
| G4-A6 | B | indirect fibre has zero active force and is classified by retained-link graph distance | makes indirect realignment falsifiable against `p_x=0` | repeat over network seeds and Euclidean shells |
| G4-A7 | C–D | 12 symmetric sites × 12 effective clutches; Bell slip law with G2 V3 starting values | Bell/Bangasser/Prahl support the mechanism; existing G2 provides a reproducible starting calibration | fit dwell, force and recoil distributions to videos |
| G4-A8 | C–D | update a site only when its entire bundle is detached | bound clutches must remain at one material coordinate; professor requested updated sites as fibres/cell move | compare alternative capture/search rules |
| G4-A9 | D | released rigid cell translates and rotates only from equal-and-opposite ECM reaction | isolates movement without `+x`, prescribed velocity or `polarity_probability` | add persistence feedback only after unbiased ensembles |
| G4-A10 | A–D | one-coefficient-at-a-time website cases are exact precomputed trajectories | preserves cause/effect and avoids a duplicate browser solver | add multi-parameter response surfaces after one-axis calibration |

Full G4 values, derived slip landmarks and evidence boundaries are in
[`docs/g4_parameter_provenance.md`](docs/g4_parameter_provenance.md).

## Important negative result

At the accepted 5 nN baseline, V4 creates only one pair that satisfies the
new-approach criterion, and its excess residual over the elastic control is
below the saved numerical resolution.  Therefore this model currently supports
**no claim of permanent collagen remodeling**.  Larger or repeated loads,
longer maturation, fibre sliding/merging, or a different network calibration
must be tested as explicit hypotheses rather than introduced only to make the
animation dramatic.
