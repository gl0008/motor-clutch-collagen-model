# G3B calibration report -- 2026-08-18

> Personal simulation / mechanism validation; not yet a realistic 3D tumor-migration
> prediction.

## Decision

G3B calibration is complete, but the stage is **not cleared for final validation**. Do not
freeze the configuration, inspect seeds 1000--1099, or resume G3C until the isotropic-fixture
geometry-overlap failure is resolved and the full G3B calibration gate is rerun.

## Campaign definition

- Git commit at launch: `1ea4f31a271363d92dd42e6f4f219258ad623ba4`
- Phase: calibration only
- Seeds: 0--19
- Simulated duration: 600 s per condition/seed
- Conditions: balanced, isotropic, aligned, aligned rotated by 30 degrees, aligned with
  geometry/traction feedback disabled, and no-fibre control
- Total G3B records: 120
- Final-validation seeds 1000--1099 used: none

The resolved SI-unit configuration and campaign manifest are stored alongside the raw
records under `results/g3_validation/calibration/`. Because the checkpointed campaign was
resumed several times, the manifest's `workers` field describes the last resume invocation;
worker counts of 4, 6, and 8 were used across launches. Parallel worker count does not alter
the per-seed random stream or resolved physics configuration.

## Results

| Condition | Valid / total | Mean nematic guidance | Bootstrap 95% CI | Axis / symmetry result |
|---|---:|---:|---:|---|
| balanced | 20 / 20 | 0.012 | [-0.144, 0.171] | no statistically resolved axis; residual polar ratio 0.182 |
| isotropic | 13 / 20 | 0.254 | [0.063, 0.445] | **invalid subset: 7 overlap failures; do not interpret as an ensemble** |
| aligned | 20 / 20 | 0.726 | [0.626, 0.819] | inferred axis 179.23 degrees, equivalent to -0.77 degrees; 55/45 sign split |
| aligned rotated 30 degrees | 20 / 20 | 0.608 | [0.407, 0.756] | inferred axis 32.26 degrees; covariance error 3.03 degrees |
| aligned, feedback off | 20 / 20 | -0.033 | [-0.195, 0.127] | no resolved director preference |
| no fibre | 20 / 20 | 0.146 | [-0.012, 0.302] | zero bound clutches, force, and torque; angle turnover is stochastic only |

The aligned-guidance signal therefore depends on the implemented local geometry/traction
feedback and transforms correctly under a 30-degree rigid rotation. The no-fibre control
also confirms that protrusion turnover alone does not create traction.

These positive mechanism checks do not override the failed isotropic geometry gate. Seven
isotropic runs (seeds 0, 2, 3, 5, 7, 8, and 15) crossed the explicit rigid-cell overlap
threshold after 48.1--366.7 s. Their maximum simultaneous bound-clutch counts were 65--99,
compared with 27 in the aligned condition. The current model has attachment traction but no
cell--fibre excluded-volume force, so some random local geometries can recruit fibres through
the rigid cell boundary.

## Gate interpretation

Passed calibration checks:

- aligned guidance bootstrap CI is above zero;
- aligned positive/negative direction split is within 40--60%;
- 30-degree rotational covariance error is below 5 degrees;
- feedback ablation reduces the guidance metric by more than 50%;
- no-fibre traction is exactly zero;
- all non-isotropic G3B runs completed without worker errors or geometry overlap.

Failed or unresolved checks:

- isotropic all-runs-valid gate: failed (13/20 valid);
- balanced residual polar ratio <0.1: not met at n=20 (0.182);
- isotropic residual polar ratio <0.1: not testable from the biased 13-run valid subset.

The residual polar-ratio threshold was preregistered for the 100-seed final ensemble. Its
20-seed calibration value is diagnostic, not by itself evidence of hidden polarity. The
overlap failure is different: it is a structural invalidation and must be resolved before
opening the final seeds.

## Premature G3C records

A separate campaign process launched with `--stage both` before the G3B gate was available.
It produced 169 G3C calibration records, including 8/20 overlap failures in the isotropic
condition, and was then stopped. These files are retained as exploratory audit records only.
They are incomplete, were generated before G3B clearance, and must not be cited as G3C
validation.

## Software checks

- Active G3 unit suite: 41 passed (`PYTHONPATH=src`, `tests/unit`).
- Frozen G2 V2 suite: 5 passed.
- Frozen G2 V3 suite: 4 passed.
- Frozen G2 V4 suite: 2 passed.

A repository-root unscoped pytest collection is not a supported test command because archived
generations contain repeated top-level module and test names. Each frozen generation was
therefore checked in its own directory, while G3 was checked with its documented source path.

## Required next design decision

The next change must be explicit and documented. Two scientifically different options are:

1. Add a conservative cell--collagen excluded-volume/contact law and calibrate its range and
   stiffness. This treats the rigid cell boundary as mechanically impenetrable, but adds a
   new force term that was not part of the preregistered G3 plan.
2. Redesign the minimal isotropic fixture so its fibres remain a non-penetrating contact
   geometry over the tested load range. This preserves the current force model, but narrows
   the question to a controlled mechanism fixture rather than unrestricted random geometry.

Neither correction should be hidden as a numerical patch. After choosing one, rerun all
calibration seeds 0--19 for every G3B condition, summarize the gates, freeze the configuration
only if they pass, and then run the untouched validation seeds 1000--1099. G3C follows only
after that sequence.
