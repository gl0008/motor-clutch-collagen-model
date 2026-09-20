# G5 — the leader-led strand (and its R4 guidance rescue) is CROSSLINK-DENSITY GATED

> **Personal testing, NOT confirmed findings** (CLAUDE.md §7.5). Found 2026-09-20 while auditing why
> the R0–R3 viewer felt "off" (Minnie: *leader cells 好像沒什麼用, 也沒看到 follower strand*). No swirling;
> imposed radial cue only. `model.py` untouched. Builds on the committed R4 follower-guidance layer.

## TL;DR
Minnie's skepticism was right — and deeper than "leaders do nothing": **whether localized leaders
detach and whether follower contact-guidance (R4) rescues a strand is entirely controlled by ECM
crosslink density.** It works in DENSE collagen and fails in the SPARSE G4D-parity collagen we had
switched to. So the leader-strand phenomenon is a *regime*, not a robust result.

## The experiment
Exact R4 headline config (2 h, seed 23, cc_adhesion=1.5 = the low-adhesion "failure case", 3 cue-front
leaders at per-site stall ×6), run guided vs unguided at TWO crosslink densities. `leader_comp_followers`
/ `leaders_detached` come from the honest adhesion-graph `strand_report` (not the retracted
separation/aspect metrics).

| crosslink_fraction | crosslinks/fibre | guidance | leader-comp followers | leaders detached | LCC | detached | mean inv |
|---|---|---|--:|:--:|--:|--:|--:|
| **0.85 (dense)** | 2.66 | none | **0** | **yes** | 0.47 | 0.21 | +29.1 |
| **0.85 (dense)** | 2.66 | R4 | **10** | no | 0.63 | 0.16 | +28.2 |
| **0.62 (G4D parity)** | 1.66 | none | **11** | no | 0.68 | 0.11 | +21.7 |
| **0.62 (G4D parity)** | 1.66 | R4 | **1** | no | 0.58 | **0.26** | +26.5 |

## What it means
- **Dense ECM (0.85):** high-traction leaders, transmitted through a well-crosslinked network, advance
  far and **rip away from the cohesive bulk** (unguided → 0 followers, LCC 0.47). Follower
  contact-guidance then pulls followers along the leader-remodelled aligned track → **10 followers kept**
  (LCC 0.63). **R4 works** — reproduces the committed R4 headline.
- **Sparse ECM (0.62, G4D parity):** the same leader traction is dissipated locally; leaders do **not**
  detach (unguided already keeps **11 followers**, LCC 0.68). There is no failure to rescue, so the extra
  guidance force just over-drives followers outward → **more scatter** (detached 0.11 → 0.26, LCC 0.68 →
  0.58). **R4 hurts.**
- **Robust across density:** force-pair = 0 everywhere; and R1 (partial-EMT lowers adhesion → single-cell
  escape) holds regardless. Those are the trustworthy results.
- **Interpretation:** ECM crosslink density is a *control parameter* for leader-follower cohesion —
  denser matrix → leaders detach → guidance matters; sparser matrix → collective holds → guidance is
  counterproductive. This is a genuine (small) mechanical phase-dependence, directly on-theme for the
  project's mechanical phase-diagram goal, NOT a bug.

## Consequence for the earlier headlines
Switching to G4D-parity crosslinks (0.62) for visual consistency **moved the system out of the regime
where leaders detach**, which is why the 0.62 R2/R3/R4 headlines looked like "leaders do nothing / no
strand". The leader physics is alive at 0.85, dormant at 0.62. Both are shown in the viewer.

## Interactive viewer narrative (docs/g5-organoid.html)
① cohesive baseline → ② partial-EMT escape (robust) → ③ leaders on cohesive matrix (no strand) →
④ DENSE + leaders, no guidance (leaders detach) → ⑤ DENSE + guidance (followers follow — R4 rescue) →
⑥ SPARSE G4D-parity + leaders (organoid holds) → ⑦ SPARSE + guidance (over-scatters — R4 fails).

## Honest limitations
- Single seed (23); one adhesion (cc=1.5) and one leader set (3 × stall 6). A full cc × crosslink phase
  map is the natural next step to draw the boundary of the "leaders detach" regime.
- Guidance params (strength 20 nN, range 25 µm, align_min 0.30) are tuned, not calibrated.
- The R3 relay's 16 switches / 2 h over-weights leader switching (leader lifetime is 120–480 min) — kept
  out of the phase-story viewer; treat as a separate timescale caveat.

## Reproduce
```bash
python generations/g5_organoid/experiments/gloria_headlines.py all r4        # sparse 0.62
python generations/g5_organoid/experiments/gloria_headlines.py all r4dense    # dense 0.85
python generations/g5_organoid/experiments/export_web_viewer.py               # -> docs viewer
```
