# Motor–clutch collagen model: versioned research path

This repository records the model as a sequence of controlled scientific
questions.  Each version changes one conceptual layer, keeps its own source,
tests, assumptions, and animation, and states what cannot yet be concluded.

> **Current recommendation:** use V2 to explain collagen force transmission,
> V3 to explain fixed-versus-moving cell mechanics, and V4 only as an explicit
> hypothesis for permanent remodeling.  V0 and V1 are archived baselines.

## Version map

| Version | Question isolated | Added mechanism | Cell | Crosslinks | Visualization |
|---|---|---|---|---|---|
| [V0](versions/v0_sls_prototypes/) | Can SLS relaxation alter clutch lifetime? | 1D SLS chain and old SLS lattice | fixed / early point cell | lattice | archived templates |
| [V1](versions/v1_few_fiber/) | Can elastic fibers reorganize under overdamped pull? | five connected bead chains; no SLS | fixed | none | [five-fiber lab](versions/v1_few_fiber/demo/index.html) |
| [V2](versions/v2_crosslinked_elastic/) | How does local force spread through collagen? | µm scale, finite fibers, permanent hinged links, contact-shell + Gaussian | fixed | permanent | [crosslinked lab](versions/v2_crosslinked_elastic/demo/index.html) |
| [V3](versions/v3_two_sided_migration/) | Can left/right clutch imbalance move the cell? | Bell clutches, force–velocity motor, cell force balance | fixed vs moving | permanent | [synchronized lab](versions/v3_two_sided_migration/demo/index.html) |
| [V4](versions/v4_plastic_remodeling/) | What remains after force removal? | stress-activated link turnover and rest-state reset | fixed load–unload test | permanent vs dynamic | [remodeling lab](versions/v4_plastic_remodeling/demo/index.html) |

The visualizations draw each fiber as a continuous polyline.  Beads are its
numerical discretization, not independent particles that appear and disappear.

## Causal progression for a presentation

```text
V0  SLS hypothesis and first motor–clutch implementation
 ↓  remove bond-memory to isolate professor's overdamped equation
V1  stretching + bending + active pull + drag, but no network coupling
 ↓  correct physical scale and add permanent links
V2  contacted fiber → crosslink force transmission → remote elastic response
 ↓  add opposing stochastic adhesion systems and release cell constraint
V3  left/right force imbalance → rigid-cell translation
 ↓  only after elastic and migration baselines are understood, test unloading
V4  stress-activated link turnover → candidate residual remodeling
```

## Run and test

Python and NumPy are the only model dependencies.

```bash
python3 -m unittest discover -s versions/v2_crosslinked_elastic/tests -v
python3 -m unittest discover -s versions/v3_two_sided_migration/tests -v
python3 -m unittest discover -s versions/v4_plastic_remodeling/tests -v

python3 versions/v2_crosslinked_elastic/run.py
python3 versions/v3_two_sided_migration/run.py
python3 versions/v4_plastic_remodeling/run.py
```

Each `run.py` precomputes the physics into `demo/data.js`.  The HTML pages only
play those data; there is no hidden second physics implementation in JavaScript.

## Default V2 geometry

- cell radius: 9 µm (18 µm diameter);
- domain: 100 × 100 µm;
- fiber contour length: 24–78 µm;
- effective displayed fiber diameter: 0.30 µm;
- bead spacing: 0.75 µm;
- direct contact shell: 3 µm;
- within-shell Gaussian width: 1.5 µm;
- far ends fixed; crosslinks permanent and freely hinged.

The collagen length and thickness ranges originate from Lee et al.'s
cell-scale bead-network paper.  The 18 µm cell is a transparent representative
tumor-cell geometry, not a universal value.  Mechanical parameters remain
uncalibrated model units.  See [saved default results](RESULTS.md),
[ASSUMPTIONS.md](ASSUMPTIONS.md), and [CITATIONS.md](CITATIONS.md).

## GitHub Pages

After this branch is merged and GitHub Pages is configured to use **GitHub
Actions**, the included workflow publishes a presentation landing page and all
four interactive labs.  The expected URL is:

`https://gl0008.github.io/motor-clutch-collagen-model/`
