# G3 model status and evidence boundary

## Scientific claim

G3 is a **2D minimal mechanism model of collagen remodelling by a contractile cell spheroid**.
It tests whether a cell **spheroid** can, *without any prescribed direction or polarity
probability*, send protrusions out **in every direction** to grip collagen across an initially
fibre-free gap, and reorganise the surrounding fibre network **from a disordered/tangential
tangle into a radial pattern (an aster)** by motor–clutch traction — while the cell itself
stays essentially in place. The headline result is the matrix remodelling, not cell motion.

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
3. **Broad intrinsic activity (no 0.65)** — a mass-conserved membrane-activity field kept
   deliberately broad (low autocatalytic gain, strong diffusion), so protrusions probe and
   grip **all around** the spheroid rather than at a single front. **The old G2
   `polarity_probability = 0.65` side-bias is removed**; gripping is symmetric with no
   prescribed direction.
4. **Softened, lightly crosslinked collagen** — the modulus is lowered and only a fraction of
   the intersection crosslinks are kept (declared G3 material choices), so a contractile cell
   can visibly rotate the near-field fibres into radial tracts. Excluded-volume repulsion keeps
   beads off the cell body, and the symmetric reaction leaves the spheroid essentially fixed.

Design rationale, literature grounding, and run instructions are in
`generations/g3_spheroid_guidance/README.md`.

## Result of the rebuilt reference run

120-fibre / 220 µm field / 20 µm spheroid / 9 µm gap / 90 min, seed 5
(`results/g3_spheroid_guidance/main/`, GIFs in `figures/g3_spheroid_guidance/`):

| Quantity | Value | Literature target |
|---|---|---|
| Engaged protrusions at end | 24 of 24 — all around | — |
| Peak collective traction | ~123 nN | tens of nN, ≤~100+ (Steinwachs 2016; Mark 2020) |
| **Near-shell radial order** | **−0.49 → −0.04** (tangential → approaching radial) | realignment over tens of min–h |
| Max fibre displacement | ~12 µm | — |
| Net spheroid displacement | ~0.6 µm in 90 min (essentially fixed) | the point is remodelling, not motion |

The collagen visibly reorganises from a disordered/tangential tangle into a radial aster, and
the near-shell radial order climbs monotonically (it is still rising at 90 min). Two GIF view
modes are rendered, matching the G2-v3 notebook toggle: `*_full.gif` (full 180 µm-style field)
and `*_follow.gif` (follow-cell zoom). The interactive page shows a live radial-order readout.

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
