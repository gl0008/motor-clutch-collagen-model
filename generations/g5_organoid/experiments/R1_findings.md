# G5-R1 — discrete partial-EMT composition (adhesion only) + isotropic-random network

> **Personal testing, NOT confirmed findings** (CLAUDE.md §7.5). No swirling. Built on the R0
> force-consistent driver; `model.py`/`visualize.py`/`test_g5.py` byte-for-byte untouched (old G5 =
> control). Code: `generations/g5_organoid/consistency.py` (extended).

## Two changes in R1

**1. Decoupled EMT phenotype from leader role.** `emt_phenotype(centers, cfg)` assigns a per-cell
adhesion multiplier — `1.0` epithelial, `emt_adhesion_factor` (a_pEMT) for the `emt_fraction`
(f_pEMT) subset (random-seeded or boundary-enriched) — that scales **ONLY cell-cell adhesion**
(`run_r0_invasion(adhesion_scale=...)`), never clutch/drag/motor budget. Verified: an isolated
cell's clutch traction is **invariant** to emt_fraction (gate test). Pair rule `min` (weakest-member)
+ `geomean` sensitivity control. This is the plan's core: phenotype ≠ leader-function.

**2. Isotropic-RANDOM initial network (fixes the "too corona" bias).** `make_random_organoid`
builds a G4D-faithful network — isotropic fibres (angle ~ U[0,π)) + union-of-cell-disks void
(`g4_v3 make_random_void_spec`, which G4D inherits) — so **radial alignment is an OUTPUT, not seeded**.
`run_r0_invasion(network_mode="random")` is now the default; `model.make_organoid`'s corona is the
legacy control. Measured t=0 global radial order: **random ≈ 0** (−0.06 at radius 32, +0.005 at radius
40) vs **corona +0.335** (strongly radial-biased) — and random grips MORE cells (12/13 vs 6/13). Needs
~300 fibres (vs 200). All R1 numbers below are on the random network.

## Gate — `test_r1_emt.py` (7) + `test_r0_consistency.py` (7) + `test_g5.py` (15) = **29 pass**
Force-pair residual 0 and floppiness 0 in every run (mechanics stay consistent on the random network).

## Sweeps (900 s, cc=6 base, seed 23, ~13 cells)

**(a) Homogeneous adhesion window** (f_pEMT=0): low adhesion → faster spread, cluster stays cohesive.
| k_adh | Rg | detached | LCC | mean_inv µm | rate µm/h |
|--:|--:|--:|--:|--:|--:|
| 1 | 31.2 | 0.00 | 1.00 | +6.5 | 26.1 |
| 6 | 26.5 | 0.00 | 1.00 | +2.0 | 8.1 |
| 60 | 24.7 | 0.00 | 1.00 | +0.5 | 2.2 |
→ **cohesive-but-deformable across the whole range at 900 s** (LCC=1, detached=0); adhesion sets
invasion RATE (12× from k=60→1). (The corona had spurious detachment at k=1 — a near-field artifact.)

**(b) Partial-EMT composition** (a_pEMT=0.25): f_pEMT 0→1 raises rate **8→23 µm/h (~3×)**, Rg 26→30.
**(c) EMT adhesion factor** (f_pEMT=0.5): a_pEMT 0.75→0.15 raises mean invasion **+2.5→+6.0 µm**.
Both monotonic (reduced cohesion → more invasion; Ilina & Friedl 2020); connected at 900 s.

## Headline — 2 h, radius-40 (~19 cells), SAME base adhesion cc=6, ONLY EMT differs

| Case | Rg | detached | LCC | mean_inv µm | max_inv µm | fpair | flop |
|---|--:|--:|--:|--:|--:|--:|--:|
| **cohesive** (f_pEMT=0) | 30.6 | 0.00 | **1.00** | +1.95 | 4.45 | 0 | 0 |
| **partial-EMT escape** (f_pEMT=0.5, a_pEMT=0.15) | **58.7** | **0.053** | **0.63** | **+25.2** | 58.1 | 0 | 0 |

→ Flipping half the cells to partial-EMT (weak adhesion) at fixed base adhesion drives a genuine
**collective → single-cell escape**: Rg nearly doubles, the cluster **fragments** (LCC 1.0→0.63),
invasion is **13×**, and cells detach — on a genuinely isotropic matrix, with the force-consistent
mechanics intact (force-pair 0, no spaghetti). GIFs: `output/r1_cohesive_2h.gif`,
`output/r1_emt_escape_2h.gif` (121 frames each = 2 h).

## Honest limitations
- Elastic, speed-capped, no plasticity/stiffening (R0 baseline); moderate scale (~13–19 cells;
  memory-tight machine). Absolute detached-fraction is small — the strong mode signal is Rg + LCC +
  invasion, not hard detachment count.
- EMT assignment is discrete (epithelial / partial-EMT) with random placement (+ boundary control);
  a continuous EMT score is a later refinement.
- Partial-EMT here changes ONLY adhesion (first R1 layer, by design). Clutch/protrusion coupling of
  EMT is deferred to keep attribution clean.

## Reproduce
```bash
python -m pytest generations/g5_organoid/tests/test_r1_emt.py generations/g5_organoid/tests/test_r0_consistency.py generations/g5_organoid/tests/test_g5.py
python generations/g5_organoid/experiments/r1_emt_adhesion.py adhesion      # (a)
python generations/g5_organoid/experiments/r1_emt_adhesion.py composition   # (b)
python generations/g5_organoid/experiments/r1_emt_adhesion.py factor        # (c)
python generations/g5_organoid/experiments/r1_emt_adhesion.py headline      # 2 h gifs
```
