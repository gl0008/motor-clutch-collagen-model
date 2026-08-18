# G3 — contractile spheroid remodels collagen into a radial pattern (rebuild on the G2 engine)

This is a clean rebuild of Generation 3. It is a **2D mechanism demonstration** of how a
cancer-cell **spheroid** can, without any prescribed direction, send protrusions out **in every
direction** to grip collagen across an initially fibre-free gap, and reorganise the surrounding
fibres **from a disordered/tangential tangle into a radial aster** by motor–clutch traction —
while the cell itself stays essentially in place. The headline is the matrix remodelling, not
cell motion. It replaces the earlier `src/g3/` package (now archived under
`legacy/g3_v1_superseded/`), which was mis-scaled, ran for far too short a time to show
anything, and produced chaotic animations that did not match Gloria's G2.

> The model can also be run in a **polarised/migrating** regime (higher autocatalytic gain,
> stiffer matrix); the default here is the contractile radial-remodelling regime that the
> figures show.

> Mechanism / personal simulation. **Not** a realistic 3D tumour-migration prediction.

## Why the rebuild

The old G3 built a brand-new spatial-clutch engine from scratch and could not reproduce the
two things that actually matter: **visible fibre realignment** and **cell motion**. Its
runs were 15–60 s — but the biophysics literature says fibre realignment develops over
**tens of minutes to hours** and migration takes **hours**, so nothing had time to happen.
Its polarity field flickered between 2 of 24 sectors every few seconds ("twitchy"), and its
guidance test actually failed (aligned nematic order 0.0142 < isotropic 0.0160).

Meanwhile Gloria's **corrected Generation 2** engine
(`generations/g2_corrected/common/model.py`) already does the hard part — a crosslinked,
boundary-anchored bead–spring collagen network with overdamped dynamics, validated at
180 µm / 99 fibres / 5 nN and a realistic 0.28 µm/min migration in V3. So this rebuild
**reuses that engine unchanged** and adds only what the story needs.

## What was added on top of the G2 engine

1. **A cell spheroid with a fibre-free gap.** The network generator leaves an empty ring
   around the spheroid, so at `t = 0` there is genuinely no cell–ECM contact.
2. **Explicit protrusions.** Every membrane sector grows an exploratory protrusion; a
   protrusion can bind a fibre only once its **tip** physically reaches a material point
   across the gap (tip-first encounter). Binding then uses the G2 slip-bond motor–clutch law.
3. **Emergent polarity (no 0.65).** A mass-conserved replicator activity field on the
   membrane — local autocatalytic reinforcement + membrane diffusion + noise + a
   FAK/Rac1-like adhesion-traction feedback — makes **one stable broad front** emerge from a
   symmetric start. Noise chooses the direction (stochastic sign); the adhesion feedback
   pulls the front toward wherever protrusions actually grip collagen; aligned collagen only
   orients the *axis*. Nothing prescribes `+x`, a left/right probability, or `0.65`.
4. **Reaction-driven migration.** Exactly as in G2 V3, the spheroid velocity is the summed
   equal-and-opposite clutch reaction over a cell drag, so it migrates toward its front.

The old G2 hardcoded `polarity_probability = 0.65` is **gone**: here the binding bias is not
a parameter, it is the emergent consequence of which protrusions find and load fibres.

## The narrative the GIF shows

`t = 0` fibre-free gap → protrusions probe outward in all directions → some tips reach fibres
and bind → motor–clutches load and pull those fibres inward (network realigns radially at the
front) → a single front is selected and reinforced by adhesion feedback → the spheroid
migrates toward that front, remodelling collagen along its path.

Two view modes are rendered (matching the G2-v3 notebook toggle):
`*_full.gif` (full 180 µm-style field) and `*_follow.gif` (follow-cell zoom).

## Literature grounding (targets, not fits)

| Quantity | Model target | Source |
|---|---|---|
| Spheroid traction | tens of nN (≤100) | Steinwachs 2016; Mark 2020 |
| Force polarity | ~0.47 | Steinwachs 2016 |
| MDA-MB-231 3D speed | 0.1–0.3 µm/min | Sapudom 2019; Steinwachs 2016 |
| Protrusion reach | ~0.7–4× cell radius, lifetime 10–30 min | Carey 2016; Fraley 2010/2015 |
| Emergent polarity | wave-pinning, polarises ~10 s, stochastic sign | Mori 2008; Jilkine 2011 |
| Alignment axis, not sign | nematic FAK/Rac1 bias | Carey 2016 |
| Fibre realignment time | tens of minutes to hours | Kim 2017; Han 2018 |

## Files

- `model.py` — `SpheroidConfig`, the fibre-free-gap network generator, the emergent-polarity
  and protrusion dynamics, tip-first engagement + motor–clutch loading, a `bincount`-based
  `fast_advance` that is bit-identical to G2's `Network.advance`, and `run_spheroid`.
- `run.py` — run one simulation, save `frames.pkl`, `meta.pkl`, `traces.npz`, `summary.json`.
- `render.py` — render `*_full.gif` and `*_follow.gif` in the G2 visual grammar.

## Run

```bash
cd generations/g3_spheroid_guidance
python run.py    --out ../../results/g3_spheroid_guidance/main --fibers 150 --domain 240 \
                 --radius 22 --gap 12 --duration 3600 --dt 0.06 --seed 7
python render.py --run ../../results/g3_spheroid_guidance/main --out ../../figures/g3_spheroid_guidance
```

## Controls (for verifying the mechanism is real)

- `--fixed` clamps the cell: fibres still get pulled/realigned, but there is no migration.
- Setting the motor stall force to zero (in `SpheroidConfig`) gives attachment without fibre
  deformation.
- An empty gap with no reachable fibres yields zero engagement, zero traction, zero motion.
