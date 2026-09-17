# G5-R2 — static leader NUMBER + LOCATION (functional role, decoupled from EMT)

> **Personal testing, NOT confirmed findings** (CLAUDE.md §7.5). The imposed one-sided radial cue is
> an EXPLICIT temporary assumption, **NOT swirling-derived** (Kolade sees no swirling). Built on the
> R0/R1 force-consistent driver; `model.py`/`visualize.py`/`test_g5.py` byte-for-byte untouched.

## What R2 adds (all in `consistency.py`, additive; R0/R1 defaults unchanged)
1. **Leader = functional role, separate from EMT phenotype.** `n_leaders`, `leader_location`
   (`"cue_front"` = a LOCALIZED adjacent front of the boundary cells nearest `cue_angle`; `"perimeter"`
   = spread control). Leaders are NOT auto-low-adhesion (EMT stays the R1 knob).
2. **Leader traction EMERGES from per-site motor capacity** (`leader_site_stall` → `_clutch_step_stall`):
   leader sites get a higher `F_stall` in the force-velocity law `v=v0·(1−F/F_stall)` (Chan & Odde 2008),
   so traction is genuinely higher from the motor-clutch — **NOT** the post-hoc `×leader_pull_factor`
   multiply my exp2 used. `budget_mode="matched_total"` holds the total front stall constant vs N_L
   (isolates leader number/geometry); `"fixed_per_leader"` keeps each leader's capacity fixed.
3. **Imposed one-sided radial cue** (`make_cued_organoid` / `_cue_positions`) — faithful port of Gloria's
   v4E `make_radial_tract_spec`: rigid **centroid- and length-preserving** rotation of eligible cue-sector
   near-field fibres toward radial (invalid rotations skipped, not clipped). `radial_cue=False` default
   → pure isotropic random. Verified `max_contour_length_change ≈ 4e-14`.
4. Morphology metrics: `strand_metrics`, `leader_follower_separation`, `aspect_ratio`.

## Gate — `test_r2_leader.py` (7) + R1 (7) + R0 (7) + `test_g5.py` (15) = **36 pass**
Cue length-preserving; cue raises t=0 order on the cue side only; n_leaders=0 == R1; leaders localized
for `cue_front`; matched-total stall sum constant across N_L; **leader traction comes from stall**
(factor 1→3 raises a leader 10.06→14.01 nN; factor=1 removes it); force-pair 0 with leaders on.

## Sweeps (900 s, cued network, seed 23, ~13 cells, factor 4)

**(a) Leader count** (matched-total budget) — fewer/stronger leaders lead the front more:
| N_L | leader–follower sep µm | max_inv µm | LCC | aspect | fpair |
|--:|--:|--:|--:|--:|--:|
| 0 | +0.00 | 4.10 | 1.00 | 1.04 | 0 |
| 1 | **+3.52** | 4.08 | 1.00 | 1.04 | 0 |
| 3 | +2.24 | 5.03 | 1.00 | 1.08 | 0 |
| 5 | −0.29 | 5.03 | 1.00 | 1.05 | 0 |
→ at fixed total budget, a single concentrated leader advances most (+3.5 µm); 5 weak leaders (each
base·4/5) barely exceed followers.

**(b) Budget** @ N_L=3: matched-total sep **+2.24** vs fixed-per-leader **+1.33** (stronger fixed leaders
drag the whole front along → less *separation*).
**(c) Location** @ N_L=3 (fixed): **cue-front +1.33 sep, aspect 1.06** vs **perimeter +0.75, aspect 1.03**
→ localized leaders lead more and elongate the organoid more than spread-out ones.

## Headline — 2 h (cued, ~19 cells, cc_adhesion=6)
| Case | leader–follower sep | max_inv µm | mean_inv µm | strands | aspect | fpair |
|---|--:|--:|--:|--:|--:|--:|
| **3 cue-front leaders, stall ×6** | +1.70 | **6.89** | +2.19 | 0 | 1.03 | 0 |
| no-leader control | 0.00 | 4.72 | +1.72 | 0 | 1.03 | 0 |

## Honest key finding (reproduces the meeting observation)
Localized cue-front leaders **do** form the leading edge (red cells ahead on the cue side, higher max
reach), and per-site stall gives them genuinely higher traction — **but strong STATIC leaders alone do
NOT pull a stable finger/strand** at cohesive adhesion (strands 0, aspect ~1.03 even at stall ×6 over
2 h). This independently reproduces the meeting note *"高 traction leader 自己跑得快，但沒有穩定拖出
follower strand；不能只繼續提高 leader force."* Cranking `leader_stall_factor` raises leader reach but not
strand formation. GIFs: `output/r2_leader_front_2h.gif` (leaders red), `output/r2_no_leader_2h.gif`.

**Implication for R3 / next:** a real strand likely needs (i) lower cell–cell adhesion so leaders can
pull ahead (R1×R2 coupling — a partial-EMT leader front), and/or (ii) **dynamic leader switching (R3)** so
a fatigued leader hands off and the front advances persistently. R2 gives the clean static baseline +
the honest negative that motivates R3.

## Reproduce
```bash
python -m pytest generations/g5_organoid/tests/test_r2_leader.py generations/g5_organoid/tests/test_r1_emt.py generations/g5_organoid/tests/test_r0_consistency.py generations/g5_organoid/tests/test_g5.py
python generations/g5_organoid/experiments/r2_leader.py count      # (a)
python generations/g5_organoid/experiments/r2_leader.py budget     # (b)
python generations/g5_organoid/experiments/r2_leader.py location   # (c)
python generations/g5_organoid/experiments/r2_leader.py headline   # 2 h gifs (leaders red)
```
