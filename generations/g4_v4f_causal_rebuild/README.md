# G4 F_causal — one declared change at a time

This generation does not replace v4E, historical v4F, or the F_new lineage
audit. It rebuilds the transition from v4E so that the effect of each change
can be identified.

## Why this exists

Historical v4F changed several numerical and physical rules in one jump.
Therefore, the raw difference `v4F - v4E` cannot be interpreted as the effect
of collagen steric reaction, dynamic contact updating, or memory removal.

F_causal uses two experiments:

1. a cumulative **numerical microstage ladder**, with exactly one new switch
   per row;
2. a paired **2 × 2 × 2 mechanics factorial**, after the numerical baseline is
   locked.

## Gate 0 — exact v4E replay

`E_REPLAY` has every new switch off. It is compared directly against the
frozen v4E implementation using the same seed and duration. The gate compares
cell positions, collagen bead positions, site forces, bound-clutch counts, and
front state. No causal comparison is accepted unless this replay matches to
floating-point tolerance.

## Gate 1 — numerical microstage ladder

| Stage | Only newly enabled change | Classification |
|---|---|---|
| `E_REPLAY` | none | exact historical v4E |
| `N1_RNG` | independently mixed stochastic address | numerical repair |
| `N2_INTEGRATOR` | cell and ECM use the same time-level geometry | numerical repair |
| `N3_FRONT_SELECTION` | every mature site may compete for front formation | model-definition repair |
| `N4_INITIAL_CONTACTS` | guided and unguided cases share initial contacts | experimental-design repair |
| `N5_SHARED_REACH` | every condition has the same probing reach | experimental-design repair |

The ladder runs both a random control and a guided sentinel. Some changes,
especially front selection, are dormant in an unguided control and would be
missed by testing E0 alone.

## Gate 2 — mechanics factorial

After N1–N5 are held constant, three binary physical/model assumptions are
varied:

| Factor | Off | On |
|---|---|---|
| **S** | cell does not receive collagen steric reaction | equal-and-opposite steric reaction acts on cell |
| **C** | failed sites use the initial frozen candidate list | eligible material points are recomputed around the current cell position |
| **M** | unguided rebinding retains the inherited directional rule | unguided rebinding is Gaussian distance-only and memoryless |

This gives eight arms: `S0C0M0` through `S1C1M1`. Every arm uses the same ECM,
clutch, motor, force, drag, observation clock, seed list, and N1–N5 numerical
baseline. Main effects and pairwise interactions are calculated from paired
seed differences.

Operationally, **C = 1** also prevents binding to an invalid old point. If no
eligible current point exists, the site remains unavailable and retries every
5 s. The clock is an uncalibrated implementation assumption. Therefore C is
reported as the **dynamic-contact module**; its effect must not be described as
the effect of coordinate updating alone. Search-clock sensitivity is required
before interpreting C biologically.

The cell equation for this comparison is

\[
\zeta_c\dot{\mathbf r}_c
=
-\sum_i\mathbf F_i^{\mathrm{clutch}}
-S\sum_i\mathbf F_i^{\mathrm{steric}}.
\]

When **C = 1**, a failed site searches continuous fiber material points using
the current cell center and deformed collagen geometry. When **M = 1**, the
unguided selection probability is

\[
P_s^{0}
=
\frac{I_s\exp(-d_s^2/\sigma_c^2)}
{\sum_q I_q\exp(-d_q^2/\sigma_c^2)}.
\]

## Interpretation boundary

- N1–N5 are not biological treatments.
- S is a mechanical consistency question.
- C changes what collagen is physically reachable after the cell moves.
- M changes the rebinding rule and is explicitly a model assumption.
- A single seed is a lineage diagnostic, not biological evidence.
- Scientific inference uses paired seeds and confidence intervals.
- The four guidance conditions E0–E3 are rerun only after this shared baseline
  is selected; guidance is not mixed into the mechanics-calibration decision.

## Run

```bash
PYTHONPATH=. .venv/bin/python -m generations.g4_v4f_causal_rebuild.run_analysis --quick
PYTHONPATH=. .venv/bin/pytest -q generations/g4_v4f_causal_rebuild/tests
```

The quick analysis is a short causal smoke test. It is not a substitute for
the predeclared long-clock, 20-seed inference run.

## However

This reconstruction makes attribution possible but does not validate the
underlying 2D rigid-cell biology. The steric penalty, search clock, probing
reach, and memoryless probability are still assumptions that require
sensitivity checks or experimental calibration. A factor can also interact
with another factor; this is why the factorial reports interactions rather
than treating the three effects as automatically additive.
