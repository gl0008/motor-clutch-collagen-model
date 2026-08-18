# G3 model status and evidence boundary

## Scientific claim

G3 is a **2D minimal mechanism model of emergent cell–collagen guidance**. It tests whether
a cell **spheroid** can, *without any prescribed direction or polarity probability*, break
symmetry, send protrusions out to grip collagen across an initially fibre-free gap, remodel
the surrounding fibre network by motor–clutch traction, and migrate toward the side where it
actually grips.

It is **not** a realistic three-dimensional tumour-migration model. It does not contain a
calibrated collagen concentration/pore topology, a deformable cell or nucleus, proteolysis,
EMT, cell–cell interactions, stress relaxation, or irreversible plasticity.

> Personal simulation / mechanism validation; not yet a realistic 3D tumour-migration prediction.

## Current implementation (2026-08-19 rebuild)

The active implementation is **`generations/g3_spheroid_guidance/`**, a clean rebuild on top
of Gloria's frozen, validated Generation-2 collagen engine
(`generations/g2_corrected/common/model.py`). The engine — crosslinked, boundary-anchored
bead–spring network with overdamped dynamics — is reused unchanged; G3 adds only:

1. **A spheroid with a fibre-free gap**, so `t = 0` has no cell–ECM contact.
2. **Explicit protrusions** that probe outward and can bind a fibre only when their *tip*
   reaches a material point (tip-first encounter), then load it with the G2 slip-bond
   motor–clutch law.
3. **Emergent polarity** — a mass-conserved replicator membrane-activity field with local
   autocatalysis, diffusion, noise, and a FAK/Rac1-like adhesion-traction feedback. One
   stable broad front emerges from a symmetric start; noise sets the direction; adhesion
   feedback pulls the front toward wherever protrusions grip. **The old G2
   `polarity_probability = 0.65` side-bias is removed** — the binding bias is now an emergent
   consequence of the mechanics, not a parameter.
4. **Reaction-driven migration** (as in G2 V3): spheroid velocity is the summed
   equal-and-opposite clutch reaction over a cell drag.

Design rationale, literature grounding, and run instructions are in
`generations/g3_spheroid_guidance/README.md`.

## Result of the rebuilt reference run

150-fibre / 240 µm field / 22 µm spheroid / 12 µm gap / 60 min, seed 7
(`results/g3_spheroid_guidance/main/`, GIFs in `figures/g3_spheroid_guidance/`):

| Quantity | Value | Literature target |
|---|---|---|
| Peak collective traction | ~36 nN | tens of nN, ≤100 (Steinwachs 2016) |
| Net migration in 60 min | ~19 µm (net ≈ path — directional) | — |
| Migration speed | ~0.32 µm/min | 0.1–0.3 µm/min (Sapudom 2019) |
| Engaged protrusions at end | ~10 | — |
| Near-shell radial order | −0.62 → −0.31 (fibres pulled radially) | realignment over tens of min–h |
| Max bead displacement | ~5 µm | — |

**Emergence check** (5 seeds, short runs): the migration direction is 84°, −84°, −8°, −99°,
64° — spread ≈ 75°. The direction is genuinely stochastic and set per realisation, confirming
there is no built-in directional bias.

Two GIF view modes are rendered, matching the G2-v3 notebook toggle: `*_full.gif` (full
180 µm-style field) and `*_follow.gif` (follow-cell zoom).

## Controls that make the mechanism falsifiable

- `--fixed`: fibres still get pulled/realigned but the spheroid does not migrate.
- Zero motor stall force: protrusions attach without deforming fibres.
- No reachable fibre in the gap: zero engagement, zero traction, zero motion.

## Relationship to Generation 2

`generations/g2_corrected/` remains frozen and unchanged. G3 reuses its engine but does not
rewrite G2 files or reinterpret G2 V2/V3 as realistic tumour migration.

## Superseded G3 v1

The earlier `src/g3/` package and all its outputs (mis-scaled, sub-minute runs, chaotic
animations, a failed guidance test with aligned nematic order 0.0142 < isotropic 0.0160) have
been archived under **`legacy/g3_v1_superseded/`**. They are kept for audit only and are not
current evidence.

## Interpretation rule

All generated figures carry:

> Mechanism demo — not a 3D tumour-migration prediction.
