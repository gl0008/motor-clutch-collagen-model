# Saved results and interpretation boundaries

## Generation 1 — historical outputs

These values remain useful for documenting why the model was rebuilt.

| Legacy stage | Saved output | Why it is not the corrected baseline |
|---|---|---|
| V1 | local alignment 0.310 → 0.501; RMS displacement 0.237 model length | dimensionless five-fibre demonstration |
| V2 | radial order −0.11854 → −0.10328; RMS displacement 0.10482 µm | every finite fibre end was fixed; contacts vanished in frame zero; most motion localized into a few bonds |
| V3 fixed | centre exactly zero | valid constraint, but paired moving trajectory was only 0.18 µm and visually unresolved |
| V3 moving | final centre (−0.18043, 0.03093) µm | no cell-line speed calibration or persistent polarity |
| V4 | 7 breaks, 2 reforms; plastic and elastic residuals indistinguishable | same-partner rest reset did not establish remodeling |

## Generation 2 — accepted corrected outputs

### V2 network and force transmission

- domain: 180 × 180 µm;
- 99 fibres, 7,173 beads, 7,074 axial bonds;
- 383 permanent intersection crosslinks;
- 135 fixed beads, all in the outer boundary band;
- 100% of fibres in a boundary-connected crosslink component;
- six frame-zero contact fibres; Gaussian weights sum to one and distributed
  force error is numerically zero;
- minimum cell gap: 0.520 µm;
- final 99th-percentile absolute bond strain: 0.0287%;
- final maximum absolute bond strain: 0.1038%;
- overdamped energy-balance relative error: 0.0092%;
- crosslinked/no-crosslink intermediate-shell displacement ratio: 1.216;
- crosslinked near/intermediate/far mean displacement:
  0.00621 / 0.000567 / 0.0000948 µm.

The true-geometry displacement is small compared with the 180 µm field.  The
near-cell zoom and optional 5× **vector** overlay expose this spatial decay
without deforming the fibre geometry for display.

### V3 fixed versus released cell

- 20-seed mobility calibration: median path speed 0.280 µm/min, range
  0.251–0.312 µm/min;
- fixed cell: centre exactly `(0, 0)` for all 301 frames;
- moving cell: path 8.18 µm and net x displacement 4.76 µm over 30 min;
- saved moving-trajectory mean sampled speed: 0.289 µm/min;
- minimum cell–collagen gap: 0.171 µm;
- final 99th-percentile absolute bond strain: 2.51%;
- final maximum absolute strain: 12.99%; global saved maximum 14.33%;
- cell y coordinate is exactly zero by the declared V3 one-axis constraint.

The moving difference is now visible in true geometry via the initial dashed
cell outline and trajectory.  It is not a display magnification.

### V4 load–unload plasticity candidate

- permanent links: 383;
- genuinely newly approached weak links: 1;
- elastic current post-force RMS: 0.00512871 µm;
- new-link current post-force RMS: 0.00512872 µm;
- excess residual between conditions: below saved positional resolution;
- no penetration; global 99th-percentile strain 0.019%, maximum 0.233%.

The active force is exactly zero during the final interval, but the elastic
network is still slowly relaxing.  Permanent change must therefore be judged
as an excess over the synchronized elastic control.  No such excess is resolved
in the default run, so V4 is a documented negative baseline.

## Sensitivity and convergence notes

The complete machine-readable values are in
[`generations/g2_corrected/validation_summary.json`](generations/g2_corrected/validation_summary.json).

- 2.5/5/10 nN responses change monotonically without approaching the 5% strain
  gate in the short runs.
- 0/10/100% compression-stiffness cases have similar small-strain short-time
  outputs; this does not validate the microbuckling law at large deformation.
- 0.50/0.75/1.0 µm bead spacing changes pointwise short-time displacement, so
  quantitative parameter fitting must wait for stronger resolution convergence.
- the 140/180 µm check exceeded the chosen 5% tolerance; 180 µm was adopted.
