# G5 coupled EMT-leader variant — does EMT also raising traction change invasion?

> **Personal testing, NOT confirmed findings** (CLAUDE.md §7.5). Motivated by `BIOLOGY_AUDIT.md`
> finding B: a real partial-EMT cell does NOT change cell-cell adhesion alone — it also raises
> traction / protrusive matrix-pulling, and leaders are often the partial-EMT cells (biorxiv 2025).
> Uses only the committed R1 (adhesion) + R2 (per-site stall) machinery via
> `consistency.run_coupled_emt_invasion`; `model.py` untouched.

## The variant
`emt_traction_factor` (new `G5RevisionConfig` field): the SAME partial-EMT cells that
`emt_phenotype` weakens (↓ cell-cell adhesion) now ALSO get `emt_traction_factor × F_stall`
per site (↑ motor traction). `emt_traction_factor = 1.0` reduces EXACTLY to R1 (adhesion-only) —
verified identical `centers_final`. So the comparison is a clean decoupled-vs-coupled test at
matched `emt_fraction`/`emt_adhesion_factor`. Force-pair stays 0 (traction still EMERGES from the
clutch force-velocity law, not a post-hoc multiply).

## Result — coupling traction to EMT does NOT increase invasion (slightly DECREASES it)
2 h, isotropic-random network, ~19 cells, 10 partial-EMT cells boundary-enriched (front),
`emt_adhesion_factor=0.15`, `cc_adhesion=6`, seed 23:

| | mean invasion | max invasion | detached | Rg | LCC | fpair |
|---|--:|--:|--:|--:|--:|--:|
| **adhesion-only (R1, ×1)** | **+25.5 µm** | 68.2 | **0.26** | **61.4** | 0.47 | 0 |
| **coupled (traction ×6)** | +21.6 µm | 59.6 | 0.16 | 55.7 | 0.58 | 0 |

Coupled is lower on **every** invasion/escape metric. The 1200 s traction-factor sweep (×1,2,4,6)
was flat (mean ~4 µm, detached 0) — the divergence only appears at 2 h.

## Why — emergent motor-clutch LOAD-AND-FAIL (not hand-coded)
Short (600 s) diagnostic, same cells, ×1 vs ×6:

| traction × | instantaneous total traction | clutch **site-failures** | cumulative slips | net invasion |
|--:|--:|--:|--:|--:|
| 1 | 259.5 nN | **67** | 6711 | +3.17 µm |
| 6 | 313.6 nN | **486 (7×)** | 8764 | +2.71 µm |

The stronger motors **do** generate higher *instantaneous* traction (314 vs 260 nN), **but** they
drive the clutch bundles past what they can hold → the Bell-law off-rate `k_off0·exp(F/F_b)`
spikes → **7× more complete site failures (load-and-fail)** → the grip repeatedly ruptures and the
cell **slips** instead of advancing. This is the biphasic motor-clutch behaviour (Chan & Odde 2008;
Bangasser et al. 2017): beyond an optimal motor strength, more motor = less productive traction.
It is EMERGENT from the validated clutch physics, not imposed.

## Interpretation & implications
- **Adhesion is the DOMINANT invasion-mode knob; adding traction to the EMT cells is second-order
  (even mildly counterproductive) here.** This partially **VINDICATES R1's adhesion-only
  simplification**: although the biology co-varies traction with adhesion, in this model the
  *invasion-mode* outcome is set by the adhesion change, so R1 captures the dominant effect.
- **Testable prediction:** an EMT/leader cell that merely maximises motor traction is NOT optimal
  for invasion — there is a clutch-limited optimum. Consistent with the biology that leader
  traction is regulated (RhoA/ROCK), not maximised.

## Honest limitations
- **Single seed (23); modest effect.** The direction (coupled ≤ adhesion-only) is consistent across
  all metrics and mechanistically explained (site-failures), but robustness across seeds and the
  location of the clutch optimum (depends on `clutch_stiffness`, `n_clutches_per_site`, `bell_force`)
  are not swept here.
- Biology audit B is only PARTIALLY addressed: EMT also changes integrins/protrusions/MMP-proteolysis
  — this variant adds the traction axis only. Proteolysis remains deferred.

## Reproduce
```bash
python generations/g5_organoid/experiments/r1_coupled_emt.py sweep      # traction-factor sweep
python generations/g5_organoid/experiments/r1_coupled_emt.py headline   # 2 h adhesion-only vs coupled + gifs
```
GIFs: `output/coupled_adhesion_only_2h.gif`, `output/coupled_coupled_2h.gif` (red = partial-EMT cells).
