# G3 revision plan: cell-intrinsic polarity, explicit protrusions, and collagen feedback

**Status:** first implementation and short 99-fibre mechanism suite completed on 2026-08-18;
statistical guidance validation remains pending. The old outputs and failed/contact-enabled
calibration history remain preserved. See `docs/results/g3_revision_validation_2026-08-18.md`.

## Decision in one sentence

Revise G3 so that the **cell creates a polarity activity patch and grows a physical protrusion**,
the protrusion can attach only when its tip reaches a collagen material point, motor--clutch
traction then deforms the fibre network, and successful adhesion feeds back to stabilize the
cell's polarity. No direction or 0.65 side probability is prescribed.

The causal chain is:

```text
intracellular stochastic symmetry breaking
    -> protrusion nucleation and growth
    -> geometric encounter with a collagen segment
    -> spatial clutch binding at (fibre, segment, alpha)
    -> actin loading and equal/opposite traction
    -> local and crosslink-mediated fibre reorientation
    -> adhesion/traction feedback stabilizes the intracellular front
    -> cell translation/rotation only after the fixed-cell gates pass
```

## What the present G3 already has, and what must change

| Component | Present implementation | Revision needed |
|---|---|---|
| G3A attachment | A clutch stores `(fibre_id, segment_id, alpha)` and conservatively projects force to same-fibre beads | Keep this. Add a protrusion shaft, tip, length and explicit tip--fibre encounter before binding |
| G3B polarity | Two active sectors are initialized randomly; collagen geometry directly biases replacement and persistence | Replace the active-sector switch with a cell-intrinsic membrane activity field. Collagen may stabilize a contacted front through adhesion feedback, but must not directly assign a preferred global direction |
| G3B fibres | Eight isolated straight fibres, both endpoints fixed | Retain only as controlled fixtures. Add a crosslinked local network, then reuse the validated G2 network topology for scale-up |
| G3C motion | Equal-and-opposite clutch/contact reactions move and rotate a rigid circle | Keep as the first motion gate. Recalibrate mobility only after forces and network response pass; absolute speed remains calibrated, not predicted |
| Cell--collagen contact | Conservative, frictionless repulsion | Keep as an excluded-volume safeguard; do not treat it as adhesion or cortex mechanics |
| Visualization | Precomputed static GIF inside the notebook | Rebuild with the G2 visual grammar and interactive playback: pause, step, zoom, overlay toggles, synchronized controls and true-scale geometry |

## Revised stage definitions

### G3A-R -- one protrusion reaches, attaches, and pulls one fibre

**Question:** can one visible protrusion grow from a fixed cell, encounter a fibre, bind a
specific material point, load clutches, and reorient the fibre without violating force or moment
conservation?

This stage intentionally triggers one intracellular activity patch at a declared membrane
location. That trigger is a diagnostic fixture, not a directional migration parameter.

Add the following state for protrusion `p`:

- base angle `theta_p` and base point on the cell surface;
- length `L_p`, tip position `x_tip,p`, and state
  `growing | searching | attached | retracting`;
- the contacted `(fibre_id, segment_id, alpha)` and tip-to-material-point distance;
- a pool of clutches assigned to this protrusion;
- bound fraction, traction, actin speed, and attachment age.

Minimal kinematics:

\[
\mathbf x_{base,p}=\mathbf x_c+R_c\mathbf n_p,\qquad
\mathbf x_{tip,p}=\mathbf x_{base,p}+L_p\mathbf n_p,
\]

\[
\dot L_p=v_g a_p\left(1-\frac{L_p}{L_{max}}\right)-v_r(1-a_p).
\]

The tip may bind only when the closest collagen segment point is within `capture_distance_tip`.
Before contact there is no clutch force. After contact, preserve the current material-coordinate
rule and Bell/motor--clutch mechanics. The protrusion is a geometric and mechanical connector;
it must not apply an additional hidden pulling force.

Required controls:

- no fibre within reach -> protrusion grows/retracts, zero clutches and zero traction;
- binding disabled -> tip reaches fibre but cannot pull it;
- motor force set to zero -> attachment occurs without fibre deformation;
- rotated fixture -> all points and vectors rotate covariantly;
- force and first-moment residuals remain below the current G3A gates.

### G3B-R -- cell-intrinsic polarity with adhesion feedback

**Question:** can an initially symmetric cell spontaneously select a front, probe collagen, and
have successful fibre attachment stabilize that front without any `+x`, left/right probability,
or collagen-assigned sign?

Represent the cell perimeter with the existing sectors, but give every sector a continuous
intracellular polarity activity `a_i` rather than a Boolean active flag. Use a finite activator
pool (or equivalent global inhibition), local positive feedback, nearest-neighbour membrane
diffusion, and noise:

\[
\dot a_i = b\left[k_0+k_{fb}\frac{a_i^h}{K_a^h+a_i^h}
+k_{FAK}q_i\right]-k_d a_i+D_\theta\Delta_\theta a_i+\sigma_a\xi_i(t),
\]

\[
b=\max\!\left(0,\frac{A_{tot}-\sum_i a_i\Delta s}{V_{cyto}}\right).
\]

Here `q_i` is not collagen alignment itself. It is a filtered **successful adhesion signal**:

\[
\tau_q\dot q_i=
\left(\frac{N_{bound,i}}{N_{clutch,i}}\right)
\frac{T_i}{T_i+T_*}-q_i.
\]

This separation is important:

1. basal noise plus intracellular feedback can break symmetry in an isotropic/no-fibre setting;
2. collagen geometry affects whether a growing tip can physically encounter a fibre;
3. bound clutches and traction create the FAK-like stabilizing feedback;
4. aligned collagen can orient the polarity **axis**, while stochastic intracellular dynamics
   choose one of the two signs along that axis.

The wave-pinning literature supports this model class, and Carey et al. support the sign of the
alignment/FAK/Rac1 feedback. Neither paper supplies the coefficients above. Calibrate activity
turnover, protrusion lifetime and length distributions to imaging; do not select coefficients to
force a desired migration direction.

Required ensemble controls:

- identical parameters and zero initial directional bias in every sector;
- isotropic ECM: ensemble mean polarity vector approaches zero, although individual cells polarize;
- aligned ECM: polarity is nematically aligned with the collagen director, with 40--60% sign split;
- rotate the ECM by 30 degrees: the polarity-axis distribution rotates by 30 degrees;
- no fibre: cell-intrinsic activity can polarize, but adhesion signal, clutch force and ECM response are zero;
- `k_FAK = 0`: front formation may remain, but collagen-dependent stabilization is strongly reduced;
- `k_fb = 0`: no stable cell-intrinsic front should survive from noise alone.

### G3C-R -- crosslinked multi-fibre network and reaction-driven cell motion

**Question:** after G3A-R and G3B-R pass, can the same protrusion--adhesion loop reorganize a
connected collagen network and generate rigid-cell translation/rotation from reaction forces?

Use two network tiers:

1. **Local crosslinked fixture:** 8--20 fibres around the cell for debugging, with explicit
   crosslinks and at least one indirect force-transmission path.
2. **G2-scale network:** reuse/adapt the validated 99-fibre, boundary-connected G2 generator,
   crosslink topology, compression-softened fibres, near/middle/far metrics and rendering.

Do not create a crosslink at every 2D line intersection in the eventual 3D model. In the 2D
mechanism stage, label intersection links as a declared planar-network approximation and sweep
link density independently.

Cell motion remains overdamped and reaction-driven:

\[
\zeta_c\dot{\mathbf x}_c=-\sum_c\mathbf f_c+\mathbf F_{contact},\qquad
\zeta_r\dot\phi_c=\sum_c(\mathbf m_c-\mathbf x_c)\times(-\mathbf f_c)+\tau_{contact}.
\]

This stage should first use the rigid circle so that motion can be attributed to the validated
protrusion--clutch--network loop. Deformable cortex, nucleus and pore passage are a subsequent
generation because they add a different scientific question and many unidentifiable parameters.

Required controls:

- fixed versus released cell using the same random stream;
- crosslinks off versus on using the same network and cell history;
- no fibre and no-attachment controls produce zero motion;
- mirror geometry reverses transverse displacement and rotation;
- rotated geometry rotates the trajectory distribution;
- changing cell drag changes speed, but not symmetry or director-relative conclusions;
- report direct-fibre and remote-fibre alignment/displacement distributions, not only their means.

## Size, mass, and mobility: what to model

At cell and collagen length scales, inertia relaxes far faster than protrusion, clutch and matrix
timescales. The governing balance should therefore remain overdamped. A cell/fibre **mass ratio**
should not be added to the active model unless an inertial term is explicitly introduced and shown
to matter.

The useful ratios are instead:

| Ratio | Meaning | Use |
|---|---|---|
| `R_cell / bead_spacing` | numerical/geometric resolution of the cell--network interface | require convergence when bead spacing changes |
| `L_fibre / R_cell` | whether fibres span only the local contact zone or transmit to the boundary | choose domain/network scale |
| `L_max_protrusion / R_cell` | physical probing reach | calibrate to protrusion imaging |
| `pore_size / D_nucleus` | confinement and arrest criterion | later deformable-cell/nucleus generation |
| `zeta_cell / (N_moving zeta_bead)` | relative mobility of cell and the locally recruited fibre cluster | sweep and calibrate with cell and bead velocities |
| `N_bound k_clutch / K_ECM,local` | clutch-to-matrix stiffness competition | controls load rate, bond turnover and deformation |
| `T_contractile / T_ECM` | contraction cycle versus matrix relaxation time | needed only when reversible stress relaxation is added |

For every mobility result, report the assumed `zeta_cell`, the number of appreciably moving beads,
and a drag sensitivity. This makes clear that speed is calibrated while direction/symmetry may be
a model prediction.

## Priority order

| Priority | Work item | Why it comes here | Completion gate |
|---:|---|---|---|
| P0 | Preserve current G3 outputs; mark them as superseded mechanism baselines when G3-R begins | Prevents the new polarity law from silently changing old evidence | old results remain reproducible and version-labeled |
| P1 | G3A-R explicit protrusion geometry and tip-first attachment | Directly answers the requested visible mechanism and isolates attachment physics | reach -> attach -> load -> fibre deformation is visible and conservative |
| P2 | G3B-R cell-intrinsic activity/polarity plus adhesion feedback | Removes the hidden directional input while keeping collagen as a cue, not the source of the cell's polarity machinery | isotropic zero-mean, aligned nematic guidance, sign symmetry and ablations pass |
| P3 | Local crosslinked network, then G2-scale 99-fibre reuse | A single/isolated fibre cannot establish network realignment or remote force transfer | direct and remote response distributions plus link-off control pass |
| P4 | Nonlinear fibre response and reversible stress-relaxation module | Needed once strains leave the small-strain regime; distinct from permanent plasticity | rheology/load-hold-unload tests identify parameters before cell runs |
| P5 | Rigid-cell G3C-R ensemble and drag calibration | Motion is interpretable only after polarity and matrix gates pass | mirror/rotation/no-fibre/drag gates pass |
| P6 | 3D topology with concentration, fibre and pore distributions | Necessary before claiming realistic collagen migration | network-only rheology and image-statistics targets pass |
| P7 | Deformable cortex, nucleus and pore passage, then MMP if required | Necessary for true pore traversal and arrest, but it is a new scientific layer | pore-size/nuclear-deformation controls reproduce qualitative arrest/rescue regimes |

P4 can be developed in parallel as a network-only branch, but it should not be coupled to cell
migration until its material parameters are independently identifiable.

## Visualization specification: match the G1/G2 notebook language

The notebook structure, type, colours and explanatory cards remain consistent with Gloria's
G1/G2 pages. This is now implemented with the same architecture: Python-solved positions are
quantized into `g3-web-data.js`, then a coordinate-faithful SVG renderer supplies play, step,
scrub, full/near view and vector-overlay controls. The GIFs are secondary downloads generated
from those same saved solver frames.

### Shared encoding

- collagen: grey/blue beads joined by visible springs;
- crosslinks: gold diamonds/links;
- fixed boundary beads: black squares;
- cell: peach body with true-scale radius;
- intracellular polarity activity: a membrane heat map, separate from collagen colouring;
- protrusion shaft/tip: green line and filled tip marker;
- contacted material point: magenta marker labelled by fibre/segment/alpha in the detail panel;
- bound clutch springs: red spokes; cell reaction: blue arrow; fibre force: orange arrow;
- bond tension/compression: the existing G2 red/blue convention;
- initial geometry: faint background; optional displacement vectors: green with a declared scale.

### Stage panels

**G3A-R:** synchronized `binding disabled` versus `spatial clutches enabled`, with an event banner
showing `nucleate -> grow -> contact -> bind -> load -> detach/retract`.

**G3B-R:** synchronized `FAK feedback off` versus `feedback on` using the same intracellular noise
stream and collagen fixture. Show the membrane activity field, every protrusion length, contacted
tips and the evolving polarity vector.

**G3C-R:** offer `fixed/released` and `crosslinks off/on` as separate declared comparisons rather
than a four-way visual overload. Show cell path, body angle, direct versus remote fibres and the
local collagen director.

Every viewer needs pause, step, timeline scrub, full-field/near-cell zoom, true-scale geometry,
force/displacement overlay toggles, and a visible note when arrows or vectors are display-scaled.

### Metrics that should remain on screen

- current intracellular polarity magnitude and angle;
- protrusion state, length and lifetime;
- tip--fibre gap, bound clutches and adhesion signal;
- clutch traction and equal/opposite cell reaction;
- local director/FOI before and after pull;
- direct, intermediate and remote bead displacement distributions;
- cell position, path length, speed and angle for G3C-R;
- maximum strain, contact penetration, force error and torque error.

## Literature-to-model map

| Source | What it supports | Reasonable to include now | What it does not justify |
|---|---|---|---|
| Carey et al. 2016 | Aligned 3D collagen changes protrusion frequency, length and persistence; FAK/Rac1 are required for anisotropic response | adhesion/traction feedback stabilizing contacted protrusions; aligned-ECM and FAK-ablation validation targets | direct use of collagen alignment as a hardcoded direction probability or the coefficients in the proposed activity equation |
| Mori, Jilkine & Edelstein-Keshet 2008/2011 | Mass-conserved reaction--diffusion can spontaneously create a stable front from a symmetric state | a finite intracellular activator pool and local positive feedback as the source of cell polarity | MDA-MB-231-specific parameter values or proof that wave-pinning is the only polarity mechanism |
| Bell 1978; Bangasser & Odde 2013 | Force-dependent clutch rupture and motor--clutch load sharing | retain the current spatial Bell/motor--clutch core | a complete 3D integrin adhesion model or current coarse-grained clutch count |
| Fraley et al. 2010; 2015 | Protrusion/adhesion regulation is central in 3D; fibre alignment predicts protrusion orientation and motility across collagen conditions | protrusion-angle/lifetime metrics and alignment controls | treating detectable 2D focal-adhesion plaques as the required 3D adhesion structure |
| Lee et al. 2014 | 3D bead--spring networks, independent crosslink density/strength, fibre-scale geometry and concentration-dependent network construction | scale-up architecture and independent link-density parameter | making every 2D projected intersection a physical crosslink or transferring all fitted values without matching the experimental gel |
| Abhilash et al. 2014; Wang et al. 2014 | Cell contractility produces heterogeneous alignment and tension-dominated long-range force transmission | network-scale alignment/displacement distributions, compression softening and cell-shape sensitivity | claiming broad realignment from the current highly localized G2 mean alone |
| Licup et al. 2015 | Collagen networks stiffen nonlinearly under stress and network geometry controls onset | add a network-only nonlinear mechanics calibration before high-strain cell runs | one universal 10% threshold or a fitted law for this network without rheology |
| Steinwachs et al. 2016 | MDA-MB-231 forces in nonlinear 3D collagen and alternating contractility/elongation/speed phases | traction-scale and future migration-cycle validation | using one drag value to turn speed into an independent prediction |
| Nam et al. 2016 (strain-enhanced stress relaxation) | Relaxation accelerates with strain; force-dependent interfibre bond unbinding/rebinding is a plausible mechanism | a separately calibrated reversible junction kinetics module after the elastic/nonlinear baseline | copying the old per-bond SLS model or calling reversible relaxation permanent remodeling |
| Wolf et al. 2013; Krause et al. 2019 | Pore size, nuclear deformation, adhesion and traction constrain 3D migration phenotype and arrest | later pore/nucleus stage and its validation gates | adding a nucleus to the present attachment gate before network topology is calibrated |
| Aguilar-Rojas et al. 2024 | MDA-MB-231 motile fraction and speed vary strongly with collagen concentration | concentration-specific speed comparison after 3D network calibration | one universal target speed or a claim that current 2D fibre count equals mg/mL |

## Parameter and inference policy

Classify every parameter as one of:

1. **measured/imported:** comes from a specific experiment with compatible cell, matrix and units;
2. **calibrated:** fitted to a declared training target such as protrusion lifetime or traction scale;
3. **numerical:** resolution, timestep or visualization choice with convergence evidence;
4. **hypothesis/swept:** not known, varied across a preregistered range.

Do not use the same target both to calibrate and validate. A practical split is:

- calibrate intracellular turnover to protrusion lifetime/length distributions;
- calibrate clutch/force scale to traction or local fibre displacement;
- calibrate cell drag to path speed;
- validate direction sign symmetry, rotation covariance, FAK ablation, remote network response and
  concentration-dependent trends on held-out conditions.

## Implementation sequence in the repository

1. Create a version-preserving `src/g3r/` package or a new generation directory; do not mutate the
   saved G3 calibration records.
2. Port the current `ClutchState` material-coordinate and conservation kernels with tests unchanged.
3. Add `PolarityFieldState` and `ProtrusionState(length, tip, phase, contact)` independently of ECM.
4. Add tip-to-segment encounter and attachment state-machine tests.
5. Generate G3A-R frames and the synchronized interactive control before scaling the network.
6. Add cell-intrinsic polarity, adhesion feedback and G3B-R ensemble controls.
7. Port the G2 network/crosslink interface behind a common network protocol; add local and 99-fibre fixtures.
8. Add G3C-R motion only after G3A-R/G3B-R gates pass.
9. Run calibration seeds, freeze parameters, and only then open held-out validation seeds.

## Evidence boundary after this revision

If G3A-R through G3C-R pass, the strongest justified claim is:

> In a two-dimensional crosslinked collagen mechanism model, cell-intrinsic stochastic polarity,
> explicit protrusion--fibre encounter, spatial motor--clutch adhesion and traction feedback can
> jointly generate director-relative guidance, collagen reorientation and reaction-driven rigid-cell
> motion without a prescribed global direction.

It would still not be a realistic 3D tumor-migration prediction until 3D topology, concentration,
pore/nucleus mechanics and migration-cycle validation are added.
