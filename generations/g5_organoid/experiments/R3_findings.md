# G5-R3 — energy-based dynamic leader switching (the persistence layer)

> **Personal testing, NOT confirmed findings** (CLAUDE.md §7.5). The "motor mechanical work = ATP
> consumption" proxy is an explicit **HYPOTHESIS**, not a measured/calibrated law. Imposed cue = NOT
> swirling. Built on R0/R1/R2; `model.py`/`visualize.py`/`test_g5.py` byte-for-byte untouched.

> ## ⚠️ CORRECTION (2026-09-19) — no follower strand here either; timescale off
> Adhesion-graph re-audit of `r3_switching_2h`: the switching relay fires (leaders on/off; a promoted
> leader briefly advances ~+23 µm) but **followers never advance along the cue axis (~0) and scatter**;
> the cluster fragments to 6 components at 2 h — **no persistent follower strand**. The premise "R3 is a
> persistence layer on top of an R1×R2 directional strand" is **wrong**: R1×R2 produced no strand (see
> `R1x2_coupling_findings.md` correction). Also **Zhang 2019 leader lifetime is 120–480 min** (2–8 h),
> so at a 2 h scale switching is minor — R3's compressed-to-minutes calibration over-weights it. The
> missing "followers follow" mechanism is **R4 follower contact-guidance** (`R4_findings.md`), which is
> a prerequisite, not R3. Re-run R3 ON TOP OF R4 once R4 is settled.

## What R3 adds (`consistency.py`, additive; R0/R1/R2 unchanged)
A per-cell energy that drains as the cell does motor work and recovers when idle, driving a **relay**:
a tired front leader demotes and a fresher front-follower takes over (Zhang et al. 2019 relay-like
leader replacement). Built on the R1×R2 directional-strand front (weak glue + localized cue-front
leaders); leader **identity is now time-varying**, EMT phenotype (adhesion) is unchanged by switching.

- **Motor power (ATP proxy — HYPOTHESIS):** `P_i = Σ_{s∈i} F_s·v_motor,s`, `v_motor,s = v0·max(0,1−F_s/F_stall,s)`
  (the SAME force-velocity law the clutch step uses; Chan & Odde 2008). Gripping loaded sites only.
- **Energy:** `E_i += dt·((E0−E_i)/τ_rec − P_i/energy_cap)`, clipped ≥0.
- **Switching (hysteresis):** a candidate must be a FRONT cell (outer + within `cue_half_width` of
  `cue_angle`) AND gripping collagen. A current leader stays until `E<energy_off` (or loses grip/front)
  → demote; free slots filled by the highest-energy eligible candidate with `E>energy_on`.
  `energy_on(0.5) > energy_off(0.3)` prevents flip-flop. Per-site `F_stall` is rebuilt each step from
  the CURRENT active leaders. `leader_switching=False` → delegates to `run_r2_invasion` exactly.
- **Calibration — to LIFETIME, not ATP:** measured `P_leader≈0.34`, `P_follower≈0.17` nN·µm/s. Defaults
  `τ_rec=300 s, energy_cap=130` give an active leader steady energy below `energy_off` (drains → hands
  off in ~minutes) while followers stay above `energy_on` — a leader lifetime of Zhang 2019's 120–480 min
  ORDER, compressed for a 2 h demo. A real relay needs ≥3 front-gripping cells (organoid_radius ≥ ~36);
  fewer → the recovered original re-wins (self-cycle).

## Gate — `test_r3_switching.py` (6) + R2 (7) + R1 (7) + R0 (7) + `test_g5.py` (15) = **42 pass**
switching off == R2; motor power >0 only when gripping+loaded; energy drains under load & recovers idle;
hysteresis prevents flip-flop; ≥1 real hand-off (demote + promote of a different cell); force-pair ≈0.

## Result — relay sustains the front (2 h, cued, cc_adhesion=2, n_leaders=2, stall ×6, r=40)
| | front advance (µm) | hand-offs | distinct leaders | max_inv | fpair |
|---|--:|--:|--:|--:|--:|
| **relay ON** | **+17.6** | 47 | 3 (cells 10, 11, 15) | 59.3 | 0 |
| static OFF (energy frozen) | +13.2 | 0 | 1 fixed set | 60.3 | 0 |

- **The relay pushes the front ~34% further** (+17.6 vs +13.2 µm). The two curves **track together for
  the first ~60 min** (leaders not yet tired), then the relay **pulls ahead and keeps climbing** while
  the static leader plateaus — persistence appears exactly when the static leader would fatigue.
- **Leader energy sawtooths**: drains to the demote line (0.3) → hand-off → a fresh leader jumps in →
  drains again — 47 times among 3 front cells taking turns. `output/r3_persistence.png` (curves +
  energy trace with switch markers); `output/r3_switching_2h.gif` (red = the CURRENT leader — watch the
  baton pass). Force-pair = 0 throughout.

## Honest limitations
- **Modest effect (+34%)** at this scale/regime; overall spread (max_inv ~59–60 µm) is similar — the
  relay's gain is in the *directed front advance*, not total dispersal.
- The energy law + parameters are **calibrated to leader-lifetime ORDER (Zhang 2019), not to measured
  ATP**; "motor work = ATP" is a hypothesis. Different `energy_cap`/`τ_rec` shift hand-off frequency.
- Elastic, speed-capped, no plasticity; ~19 cells (memory-tight machine). Needs ≥3 front-gripping cells
  for a genuine relay (not a self-cycle).

## Reproduce
```bash
python -m pytest generations/g5_organoid/tests/test_r3_switching.py generations/g5_organoid/tests/test_r2_leader.py generations/g5_organoid/tests/test_r1_emt.py generations/g5_organoid/tests/test_r0_consistency.py generations/g5_organoid/tests/test_g5.py
python generations/g5_organoid/experiments/r3_switching.py headline   # 2 h relay vs static + gif + plot
```
