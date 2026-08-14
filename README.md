# Motor–clutch collagen model

This repository now keeps **two complete model generations on one website**.
The first generation is preserved exactly as the historical development path;
the second generation corrects the boundary conditions, network coupling,
physical units, migration comparison, and visualization without overwriting the
older work.

- Website: <https://gl0008.github.io/motor-clutch-collagen-model/>
- Frozen pre-correction state: Git tag `legacy-generation-1-2026-08-14`

## Generation 1 — original V0–V4 archive

The original folders remain under [`versions/`](versions/).  They document how
the project progressed from the SLS hypothesis to the first few-fibre,
crosslink, moving-cell, and rest-state-reset prototypes.  Their old URLs remain
valid.  They are useful for explaining the research history, but their V2–V4
animations should not be used as the corrected physical baseline.

| Stage | Historical question | Animation |
|---|---|---|
| V0 | Can SLS relaxation alter one protrusion's clutch lifetime? | archived code/templates |
| V1 | Can five elastic bead chains respond to active pull? | [legacy V1](versions/v1_few_fiber/demo/index.html) |
| V2 | First fixed-cell crosslinked network | [legacy V2](versions/v2_crosslinked_elastic/demo/index.html) |
| V3 | First fixed-versus-moving comparison | [legacy V3](versions/v3_two_sided_migration/demo/index.html) |
| V4 | First same-partner rest-reset plasticity test | [legacy V4](versions/v4_plastic_remodeling/demo/index.html) |

## Generation 2 — corrected baseline

Generation 2 uses one shared bead–spring collagen engine under
[`generations/g2_corrected/`](generations/g2_corrected/).  There is no SLS in
these models.

| Stage | Question isolated | Difference between synchronized panels | Animation |
|---|---|---|---|
| [V2](generations/g2_corrected/v2_crosslink_transmission/) | Do permanent links transmit a local cell pull to non-contact fibres? | crosslink force off vs on | [corrected V2](generations/g2_corrected/v2_crosslink_transmission/demo/index.html) |
| [V3](generations/g2_corrected/v3_two_sided_migration/) | Does a persistent left/right clutch imbalance move a released rigid cell? | fixed centre vs x-axis translation | [corrected V3](generations/g2_corrected/v3_two_sided_migration/demo/index.html) |
| [V4](generations/g2_corrected/v4_contact_plasticity/) | Does the validated baseline bring new fibre pairs close enough to form weak links and leave excess post-unload change? | elastic vs approach-gated new links | [corrected V4](generations/g2_corrected/v4_contact_plasticity/demo/index.html) |

### Shared corrected mechanics

- 180 × 180 µm domain, adopted after the 140/180 µm boundary check exceeded
  the preregistered 5% tolerance;
- 99 finite 20–80 µm collagen fibres, 0.30 µm effective diameter;
- 0.75 µm bead spacing (with 0.50/0.75/1.0 sensitivity saved);
- only beads in the outer boundary band are fixed;
- all cell-contact fibres and at least 85% of fibres must connect to the outer
  boundary through the crosslink graph; the saved seed reaches 100%;
- axial stiffness `EA/l0`, bending coefficient `EI/l0^3`, and 10% compression
  stiffness as the microbuckling baseline;
- hard 3 µm contact shell plus Gaussian weighting with `sigma = 1.5 µm`;
- 5 nN total active pull (2.5/5/10 nN sensitivity saved);
- overdamped bead dynamics; no inertial term and no SLS material memory.

## What is visible in every corrected animation

- every circle is a collagen bead and every connecting segment is a spring;
- gold diamonds are permanent crosslinks; purple diamonds are newly formed weak
  links in V4;
- black squares are the **only** fixed beads and occur at the outside boundary;
- orange halos/arrows are the current direct cell contacts and their
  Gaussian-weighted forces;
- red/blue bonds indicate tension/compression;
- the faint network is the initial geometry;
- geometry always stays at true 1× scale.  The optional 5× control magnifies
  displacement vectors only;
- near-cell zoom is available without changing the mechanics.

The orange contacts are present in frame zero, so they do not disappear once
per loop as they did in the legacy V2 animation.

## Validation status

[`validation_summary.json`](generations/g2_corrected/validation_summary.json)
records the accepted gates and sensitivity runs.

- V2 network gate: **pass** — 7,173 beads, 7,074 bonds, 383 permanent links,
  and 100% boundary-connected fibres.
- V2 mechanics gate: **pass** — conserved Gaussian force, no cell penetration,
  99th-percentile strain 0.029%, maximum strain 0.104%, and overdamped energy
  balance error below 0.01%.
- Crosslinks increase the intermediate-shell response by 21.6% relative to the
  same no-crosslink network.
- V3 mobility gate: **pass** — 20-seed median path speed 0.280 µm/min for the
  MDA-MB-231 calibration target of 0.2–0.4 µm/min.
- V3 saved moving trajectory: 8.18 µm path and 4.76 µm net x displacement in
  30 min; the fixed centre remains exactly zero.
- V4 currently forms only one genuinely newly approached weak link and produces
  no resolved excess residual over the elastic control.  This is a negative
  baseline result, not evidence of permanent remodeling.

## Run and test

```bash
python3 -m unittest discover -s generations/g2_corrected/v2_crosslink_transmission/tests -v
python3 -m unittest discover -s generations/g2_corrected/v3_two_sided_migration/tests -v
python3 -m unittest discover -s generations/g2_corrected/v4_contact_plasticity/tests -v

python3 generations/g2_corrected/validate.py
python3 generations/g2_corrected/v2_crosslink_transmission/run.py
python3 generations/g2_corrected/v3_two_sided_migration/run.py
python3 generations/g2_corrected/v4_contact_plasticity/run.py
```

Python precomputes all physics into each `demo/data.js`.  The web pages only
play saved frames; JavaScript does not implement a second hidden model.

See [`ASSUMPTIONS.md`](ASSUMPTIONS.md), [`RESULTS.md`](RESULTS.md), and
[`CITATIONS.md`](CITATIONS.md) for interpretation boundaries and sources.

## Separate future biology track

The Kolade/Adebowale U937 monocyte mechanism is kept conceptually separate.
That work emphasizes adhesion-dispensable outward protrusive path opening,
whereas Generation 2 V2/V3 implements the professor's inward, attachment-based
collagen pulling experiment.  A later U937 model should not silently combine
the two force directions or cell types.
