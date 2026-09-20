# G5-R4 — FOLLOWER contact-guidance: fixing "leaders detach, followers scatter"

> **Personal testing, NOT confirmed findings** (CLAUDE.md §7.5). No swirling; imposed cue only.
> Built on the R0/R1/R2 force-consistent driver; `model.py`/`visualize.py`/`test_g5.py` untouched.
> Code: `generations/g5_organoid/consistency.py` (extended, additive).

## The problem R4 fixes (honest re-audit of R1×R2 and R3)

A visual + connectivity re-audit of `r1x2_lowadh_leaders_2h` and `r3_switching_2h` showed the
"directional strand" claim was a **metric artifact**: `leader_follower_separation` = (leader outward
advance − follower outward advance) is large **even when the leaders simply DETACH and escape while
followers scatter/retreat** — it contains no connectivity. Adhesion-graph tracing of the committed
R1×R2 run (cc=1.5, 2 h) confirms the leaders break into a **2-cell island with 0 followers by frame
90** (8 components at 2 h; `leaders_detached=True`). No follower strand ever formed. The model had
**no mechanism making followers follow** — cell velocity was only `clutch_reaction + cc_force`, and
once R1 lowered adhesion the tow-line was cut.

## Biology (deep-research synthesis; see session notes for full citations)

Followers follow via a **combination**, not one knob:
1. **Retained (partial-EMT) cell–cell adhesion** — full EMT → single-cell dissemination; hybrid E/M
   retains E-cadherin and stays collective (Pastushenko/Blanpain; Jolly). *Too-weak adhesion = single
   cells* is the exact failure mode of R1×R2 (Clusters-fingers-singles phase diagram, PLOS Comp Biol
   2025: weak adhesion→singles, intermediate→fingers, strong→no invasion).
2. **Contact guidance along the leader-remodelled aligned ECM path** — leaders specialise in *creating*
   an aligned track (force-bundled radial fibres / TACS-3 + MMP microtracks); followers specialise in
   *responding* to it (Friedl & Gilmour 2009; Ray 2017; PLOS One 2024; Nat Rev Cancer 2021).
3. **Followers' own cryptic-lamellipodia migration** along the path — not purely passive towing (MBoC 2021).
- Leader switching (Zhang 2019) is real but slow: **leader lifetime 120–480 min**, so at a 2 h scale
  switching is minor — R3's compressed-to-minutes timescale over-weights it.

## What R4 adds (all in `consistency.py`, additive; `follower_guidance=False` ⇒ R0..R3 byte-for-byte)

`contact_guidance_forces(centers, network, centers0, cfg)`: a per-cell FOLLOWER force along the local
collagen **nematic director**, oriented outward, scaled by the local 2D nematic order `S∈[0,1]`:
`F = guidance_strength · max(0, S − guidance_align_min) · n_hat`. `S≈0` on isotropic matrix ⇒ **no
force**, so guidance is LOCALIZED to genuine tracks (does not isotropically inflate the organoid).
Applied to followers only (leaders build the track). It is an **external motility force on the CELL,
never projected onto the ECM**, so the clutch reaction pair is untouched → **force-pair stays 0**.
`run_r4_invasion` = R2 localized per-site-stall leaders + this guidance.

Also added an **honest strand metric** replacing `leader_follower_separation`:
`connected_components` + `strand_report` (n_components, `leader_comp_followers`, `leaders_detached`,
`strand_reach_um`) — reads the adhesion graph directly, so "leaders escaping alone" cannot masquerade
as a strand.

## Headline result — 2 h, cc=1.5 (the EXACT R1×R2 failure case), guidance OFF vs ON

| Case | comps | leader-comp followers | leaders_detached | LCC | reach µm | sep µm (old metric) |
|---|--:|--:|:--:|--:|--:|--:|
| **unguided (control)** | 8 | **0** | **True** | 0.47 | 19.3 | +34.8 |
| **guided (R4 fix)** | 6 | **10** | **False** | 0.63 | 28.6 | +19.9 |

Adhesion-graph trace (leader-containing component, L=leaders/F=followers over 2 h):
- **unguided:** 3L+16F → 2L+12F → **2L+0F (frame 90) → 2L+0F** — leaders break off alone.
- **guided:** 3L+16F → 2L+13F → 2L+11F → **2L+10F** — leaders stay connected to 10 followers.

→ Guidance converts **"leaders detach, followers scatter"** into a **connected leader-led group**:
leaders keep 10 followers, cluster stays more cohesive (LCC 0.47→0.63), fewer fragments (8→6).
Force-pair = 0 throughout both. GIFs: `output/r4_guided_2h.gif`, `output/r4_unguided_2h.gif`
(leaders red; in unguided the reds are isolated with a gap to the blue bulk, in guided a blue
follower chain trails the reds).

**Why the old metrics got smaller (and that's correct):** guidance LOWERS `sep` (34.8→19.9) and
`aspect` (1.28→1.14) — because the unguided "large" values were the artifact of leaders escaping while
followers retreated. A real connected strand has a *smaller* separation. This is the clean
demonstration that `leader_follower_separation`/`aspect_ratio` are the wrong strand metrics; use
`strand_report`.

## Honest limitations
- Not a picture-perfect single finger: still 6 components / LCC 0.63 at 2 h (a few singles peel off).
  The fix is a genuine, MODERATE improvement in follower-following/cohesion, not a dramatic strand.
- **Scale/time-dependent** — the benefit needs the control to have fragmented first. At a SMALL/short
  fixture (250 fibres, ~13 cells, 1500 s) NEITHER run has fragmented yet, so guidance shows no benefit
  and can even be marginally negative (leader-comp followers 2 guided vs 3 unguided) from isotropic
  outward drift. The 10-vs-0 benefit is a 2 h / radius-40 effect. The R4 gate test therefore asserts
  only the robust invariants (force-pair 0; guidance does not create detachment the control lacked),
  NOT the emergent benefit.
- Guidance direction = local nematic director oriented radially outward (the invasion axis). It is a
  contact-guidance *modelling choice*; it is not a measured force. `guidance_strength=20 nN`,
  `guidance_range=25 µm`, `guidance_align_min=0.30` are tuned, not calibrated (log in provenance).
- Elastic, speed-capped R0 baseline; ~19 cells (memory-tight). cc=1.5 is deliberately low (the failure
  case); the intermediate-adhesion window (Clusters-fingers-singles) is a separate sweep to run next.
- R4 builds on the STATIC-leader R2 path; combining guidance with R3 dynamic switching is a follow-up.

## Reproduce
```bash
python generations/g5_organoid/experiments/r4_follower_guidance.py compare    # off vs on x adhesion (1500 s)
python generations/g5_organoid/experiments/r4_follower_guidance.py headline   # 2 h gifs (guided vs unguided)
```
