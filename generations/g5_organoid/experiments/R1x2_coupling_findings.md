# G5 R1×R2 coupling — does low adhesion + localized leaders pull a strand (without R3)?

> **Personal testing, NOT confirmed findings** (CLAUDE.md §7.5). Uses ONLY the committed R1
> (EMT/adhesion) + R2 (localized per-site-stall leaders) code via `consistency.run_r2_invasion` on the
> cued network — no new model code. Imposed cue, NOT swirling.

## Question
R2 showed strong STATIC leaders at *cohesive* adhesion don't pull a strand. Does **low cell–cell
adhesion + localized cue-front leaders** — combining the two knobs — make a directional strand
**without** dynamic switching (R3)?

## Result — YES for DIRECTIONALITY (2 h, cc_adhesion=1.5, cued, ~19 cells)
| Case | leader–follower sep | aspect | strands | mean_inv | max_inv | LCC | detached |
|---|--:|--:|--:|--:|--:|--:|--:|
| **low-adh + 3 cue-front leaders (×6)** | **+34.8 µm** | **1.28** | 6 | +29.1 | 69.4 | 0.47 | 0.21 |
| low-adh, no leaders (control) | 0.0 | 1.10 | 1 | +29.1 | 68.6 | 0.47 | 0.11 |

Both invade the same *amount* (+29 µm mean, ~69 µm max) and fragment equally (LCC 0.47), but **only with
leaders is the invasion directional**: the 3 leaders pull **+34.8 µm ahead** of the followers and the
organoid **elongates toward the cue (aspect 1.28 vs 1.10)**, with more detachment (0.21 vs 0.11). Low
adhesion **alone** disperses isotropically (no direction). Force-pair = 0 throughout both.

## Interpretation
- **R1 (low adhesion)** enables motility + detachment (cells *can* move/leave), but undirected.
- **R2 (localized leaders)** supplies the DIRECTION (a cue-facing front that leads by +35 µm).
- **Together** → a directed, elongated, leader-led protrusion. **Neither knob alone does it**: R2 leaders
  at cohesive cc=6 gave no strand (R2_findings); R1 low-adhesion alone gave isotropic dispersal.

## Honest caveat (metric)
`strand_metric` counts cells protruding > median + `cell_spacing` (18 µm). At the 1200 s sweep (invasions
~10–13 µm) it read 0 even where a mild leader front existed — too strict for short runs. **Use
`leader_follower_separation` + `aspect_ratio` as the sensitive directionality metrics**; `strand_count`
only becomes meaningful once invasions exceed ~18 µm (here, at 2 h).

## Implication for R3
The DIRECTED strand does NOT require dynamic switching — low adhesion + static localized leaders already
produce it. **R3 (energy-based leader switching, Zhang 2019) is therefore the PERSISTENCE layer**: keep
the front advancing as the front leaders fatigue (drain motor energy) and fresher followers take over —
not a prerequisite for strand *formation*. This sharpens R3's purpose.

## Reproduce
```bash
python generations/g5_organoid/experiments/r1x2_coupling.py sweep      # cc x leaders (1200 s)
python generations/g5_organoid/experiments/r1x2_coupling.py headline   # 2 h gifs (leaders red)
```
GIFs: `output/r1x2_lowadh_leaders_2h.gif`, `output/r1x2_lowadh_control_2h.gif`.
