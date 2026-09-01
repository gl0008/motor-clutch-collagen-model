# Generation 4 v2 — long-horizon mechanics and visible clutch failure

G4 v2 does **not** replace
[`g4_interactive_calibration`](../g4_interactive_calibration/).  The earlier
folder and [`g4-lab.html`](../../docs/g4-lab.html) remain G4 v1: the original
2–5 minute experiment that revealed why real mechanical differences were hard
to see.  G4 v2 changes the observation design while retaining a 10 µm-radius
single cell and the declared nN–µm–s parameter scale.

The interactive page is
[`docs/g4-v2.html`](../../docs/g4-v2.html).  It lazy-loads exact Python results;
the browser contains no second mechanics implementation.

## The user request this version records

1. Extend time far enough to see cell-induced collagen changes without making
   forces or stiffness arbitrarily large.
2. Explain why fixed/moving and linked/unlinked cases previously looked alike.
3. Make individual clutch slip, complete attachment loss, traction drop and
   fibre recoil visible.
4. Preserve every earlier model and its reasoning on the same website.
5. Keep the explanation, equations and citations readable while moving large
   generated arrays into on-demand data chunks.

## Why the old picture looked static

- G4 v1 stored only 2–5 minutes.  Riching et al. observed MDA-MB-231-driven
  collagen displacement for 2 h at 4 min intervals.
- In the G4 v1 baseline, mean one-hop displacement was commonly
  (10^{-3})–(10^{-2}) µm.  At a 180 µm full-field scale this is subpixel.
- G4 v1 stored geometry every 4–5 s.  Many individual stochastic rupture and
  rebinding events happened between website frames.
- The original 12 clutches were independent.  One rupture was a real stochastic
  slip, but eleven parallel clutches could keep the material-point site bound.

G4 v2 therefore uses time, local zoom, clearly labelled displacement-vector
magnification, event-aligned output and an explicit shared-load comparison.
The true bead positions are never magnified.

## Stage sequence

### G4A v2 — mechanics calibration

The cell is fixed.  The primary ECM has outer-boundary anchoring and permanent
probabilistic links.  Each exact 2 h trajectory changes one coefficient from the
baseline: pull, bending multiplier, collagen modulus, bead drag, crosslink
probability, crosslink stiffness, or Gaussian contact width.  A mobile-boundary
control detects whole-network drift.

### G4B v2 — indirect realignment

Only the 12 direct contact fibres receive the active cell force.  Every other
fibre has zero active force and is classified by its shortest path through the
actual retained-link graph:

- distance 0: direct contact;
- distance 1: one retained crosslink away;
- distance 2+: longer connected path;
- distance −1: unconnected.

For each class, the simulation saves mean, 90th-percentile and maximum bead
displacement plus mean and 90th-percentile absolute nematic orientation change.
The no-link control uses the identical initial geometry.

### G4C v2 — independent versus shared-load clutch sites

The left condition retains the G4 v1 independent Bell clutches.  The right
condition implements an equal-load-sharing cluster.  If (i) of (N) clutches
remain bound at one material-point site,

\[
f_i=\frac{F_{\mathrm{site}}}{i},\qquad
r_i=i k_{\mathrm{off}}^0
\exp\!\left(\frac{F_{\mathrm{site}}}{iF_b}\right),\qquad
g_i=(N-i)k_{\mathrm{on}}.
\]

One rupture lowers (i), raises the force and off-rate of each survivor, and
may start a cascade.  Rebinding can still rescue the cluster.  Complete site
detachment is the first passage to (i=0), not every individual rupture.

The shared site uses the fully bound independent-bundle stiffness
(K_{site}=Nk_c) while at least one effective clutch remains.  This is a
coarse-grained equal-load-sharing hypothesis: it says that compliance outside
the individual bonds maintains the site load long enough to redistribute it.
It is displayed beside—not substituted silently for—the independent model.

Both conditions retain the G4 effective defaults:

| parameter | value | interpretation |
|---|---:|---|
| (N) | 12/site | effective clutches, not 12 counted integrin molecules |
| (k_c) | 2 nN/µm | effective clutch stiffness |
| (k_{on}) | 0.055 s⁻¹ | rebinding hazard per open clutch |
| (k_{off}^0) | 0.018 s⁻¹ | zero-force rupture hazard |
| (F_b) | 1.5 nN | effective Bell force scale |
| (v_0) | 0.025 µm/s | unloaded actin speed |
| (F_{stall}) | 8 nN/site | motor stall force per site |

The page provides a 30 min overview and a 20 s event microscope sampled every
0.1 s.  An isolated-site 200-trial ensemble selects the counter seed whose
first failure is nearest the median; no rare dramatic trajectory is chosen.

### G4D v2 — released rigid cell

G4D keeps the shared-load ECM and clutch equations and releases only rigid-cell
translation and rotation:

\[
\gamma_c\dot{\mathbf r}_c=-\sum_i\mathbf F_i^{cell\to ECM},
\qquad
\gamma_\theta\dot\theta=
\sum_i(\mathbf r_i-\mathbf r_c)\times(-\mathbf F_i).
\]

Fixed and moving runs share the initial network, retained links, parameter
values and counter-addressed random proposals.  There is no (+x) force,
preferred velocity or `polarity_probability`.  A mobile-ECM control distinguishes
cell movement from rigid drift of the entire matrix.

## Complete implemented mechanics

### Overdamped beads — key equation

\[
\zeta_b\dot{\mathbf r}_i=
\mathbf F_i^{stretch}+\mathbf F_i^{bend}+\mathbf F_i^{xlink}
+\mathbf F_i^{repulsion}+\mathbf F_i^{active}.
\]

There is no SLS in G4.  Drag dissipates motion; bonds and crosslinks store
elastic energy.

### Axial bond

\[
\mathbf F_{ij}=\frac{EA}{\ell^0_{ij}}
(\ell_{ij}-\ell^0_{ij})\hat{\mathbf e}_{ij}.
\]

Compression uses 10% of tensile stiffness to represent microbuckling.

### Discrete bending

\[
E_b=\frac12\frac{EI}{\ell^3}
\left\|\mathbf r_{i-1}-2\mathbf r_i+\mathbf r_{i+1}
-\boldsymbol\kappa_i^0\right\|^2.
\]

### Permanent material-point crosslink

\[
\mathbf F_i^x=k_x[(\mathbf p_b-\mathbf p_a)-\boldsymbol\ell_x^0],
\qquad link\iff u_{intersection}<p_x.
\]

### Hybrid direct coupling — key equation

\[
0\le d_{surface}\le3\,\mu m,
\qquad w_m\propto\exp(-d_m^2/\sigma_c^2).
\]

The Gaussian is evaluated only among points that first pass the contact-shell
test; it is not a whole-network force field.

### Motor and independent Bell slip

\[
v_a=v_0\max(0,1-F_{site}/F_{stall}),
\qquad
P_{off}=1-\exp[-k_{off}^0e^{F/F_b}\Delta t].
\]

### Alignment metrics

\[
S_r=\langle2(\hat{\mathbf t}\cdot\hat{\mathbf r})^2-1\rangle,
\qquad
\Delta\theta=\tfrac12\operatorname{atan2}
[\sin2(\theta-\theta_0),\cos2(\theta-\theta_0)].
\]

## Data and token-efficient website organization

- `g4-v2-manifest.js` contains geometry metadata, case filenames and the compact
  ensemble summary.  It is the only simulation metadata loaded initially.
- `docs/g4-v2-data/*.json` contains one selected trajectory per file.  Position
  differences are quantized into `Int16` typed-array payloads and decoded only
  after the user selects that case.
- Long overviews store sparse geometry frames.  The event microscope stores
  only a near-cell subset at high temporal resolution.
- [`generated_summary.json`](generated_summary.json) records the numerical
  results without raw frames.  Future Codex work should read this summary by
  default and open raw data only for a specific visualization bug.
- G4 v1's large inline `data.js` remains untouched as an archive.

## Run and verify

```bash
python3 -m unittest discover -s generations/g4_v2_multiscale/tests -v
python3 generations/g4_v2_multiscale/build_demo.py --workers 4
```

The unit tests check the accelerated force balance against the frozen G4 v1
step, boundary release, counter-stream reproducibility, shared-load hazard
increase and the (i=0) detachment rule.  The builder also creates the 200-trial
cluster ensemble and compact per-case summaries.

## Evidence and limitations

- Riching et al., *3D Collagen Alignment Limits Protrusions to Enhance Breast
  Cancer Cell Persistence* (2014), PMCID: PMC4255204 — supports 2 h / 4 min
  observation scale.
- Erdmann & Schwarz, *Adhesion clusters under shared linear loading* (2004),
  arXiv:cond-mat/0403552 — supports equal-load-sharing cluster kinetics.
- Gao et al., *Probing mechanical principles of focal contacts...* (2011),
  PMCID: PMC3140725 — supports stochastic–elastic load redistribution and
  possible catastrophic cluster failure.
- Bell (1978), DOI: 10.1126/science.347575 — supports exponential slip-bond
  force sensitivity.

G4 v2 remains 2D, elastic and single-cell.  Crosslinks are permanent; there is
no plasticity, proteolysis, nonlinear strain stiffening, cell deformation,
nucleus or 3D pore constraint.  The effective nN clutch values are a
mechanism-calibration scale and must not be presented as single-molecule fits.

## HOWEVER — what must come after G4 v2

If the long-horizon model and event microscope pass their stated controls, the
next generation can calibrate against experimental displacement/slip/recoil
videos and then add one missing biology at a time.  Plastic crosslink turnover
or collagen damage remains a separate later hypothesis; it must not be used to
rescue a failed elastic baseline.
