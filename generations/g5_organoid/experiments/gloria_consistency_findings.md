# G5 ↔ Gloria g4-v4e consistency: why our fibres "looked different", and the fix

> **Personal testing, NOT confirmed findings** (CLAUDE.md §7.5). Motivated by Minnie's question
> (2026-09-19): *"為甚麼我們的 fiber 看起來差這麼多? 我想要一致的 setting 跟 the same visualization"*,
> referencing <https://gl0008.github.io/motor-clutch-collagen-model/g4-v4e.html>. `model.py` /
> `visualize.py` / `test_g5.py` untouched. No swirling claim implied ([[project-swirling-caveat]]).

## TL;DR
The fibre **geometry already matched G4D**. What made our movies look different was (a) our old
animation **colour-mapped** the fibres (coolwarm by radial order) instead of drawing neutral grey
lines, and (b) our **crosslink density** was higher. Both are now fixed: a Gloria-style renderer +
an event microscope, and a G4D-parity crosslink density that we verified keeps the mechanics clean.

## 1. Diagnosis — the fibre geometry was already consistent
Drawn in the **same style at true scale** (`output/_compare_g4d_vs_g5_network.png`), G4D's
single-cell network and our G5 organoid network are the same kind of mesh, because our
`consistency._random_isotropic_spec` is a port of g4_v3 `make_random_void_spec` (which G4D inherits):

| property | Gloria G4D | our G5 organoid | consistent? |
|---|---|---|---|
| fibre angle | isotropic U[0, π) | isotropic U[0, π) | ✅ |
| fibre shape | straight (`curvature_amplitude=0`) | straight (random build adds no wave) | ✅ |
| fibre length | 20–80 µm (mean **45.5**) | 20–80 µm (mean **47.4**) | ✅ |
| bead areal density | 0.237 /µm² | 0.248 /µm² | ✅ |
| void | one circular void (single cell) | union-of-cell-disks void (organoid) | ⚠️ by design (organoid) |
| **crosslinks / fibre** | **1.67** | 2.66 (old default) → **1.66 (new)** | ✅ after fix |
| domain | 180 µm | 280–340 µm (holds the organoid) | ⚠️ by design |

So the "太不一樣" impression came from **rendering**, not the fibres.

## 2. Fix A — clean Gloria-style renderer
New module `experiments/gloria_style_viz.py` (does not touch `model.py` / `visualize.py`).

**Final readable style** (Minnie 2026-09-19: the first cut had ~228 per-grip arrows + gold crosslink
dots + red × everywhere = unreadable; the event microscope's force trace failed+rebound every frame
= a forest of red lines). The shipped look is deliberately minimal — `render_gloria_style` /
`render_frame_png` with `arrows="none"`:

- neutral **light grey** collagen fibres at true scale (faint background);
- cells as filled circles — **blue** = cohesive, **red** = partial-EMT / leader (`red_per_frame`
  moves the red circle for R3's baton);
- faint **displacement trails** (t0 → now) so "who moved / who escaped" is obvious at a glance;
- white background, 50 µm scale bar, one-line title (stage, red meaning, mean invasion).

The renderer keeps optional knobs for the dense diagnostic look (`arrows="net"` = one resultant
arrow per cell, `arrows="site"` = per grip, `show_crosslinks`, `show_failures`) but the headlines
use none of them. The **event microscope** (`render_event_microscope`) was **dropped** from the
headlines as unhelpful (it stays in the module, unused).

These consume **new rich per-site snapshots** added to `consistency.run_r0_invasion` **and**
`run_r3_invasion` (only when `snapshots=True`; backward-compatible, dynamics byte-identical, and
the `snapshots=False` path is unchanged — verified by tests): `site_force_snapshots (F,S)`,
`site_point_snapshots (F,S,2)` (grip points interpolated on the current beads via
`_current_site_geometry`), `site_normal_snapshots`, `site_failed_snapshots` (per-interval failure
flags), `crosslink_edge/alpha`, `n_contact_sectors`.

## 3. Fix B — crosslink density = G4D parity, verified mechanically safe
G4D uses ~**1.67 crosslinks/fibre**. Naively lowering ours trips the **binary
`contact_fibers_connected` gate** (our *larger organoid void* means the fibres gripping the cells
must cross a wider annulus to reach the anchored outer boundary than G4D's small single-cell void).
But that flag is **over-conservative**: an empirical test of the *actual* R0 mechanics (30-min runs)
showed the network stays clean well below it —

| crosslink_fraction | crosslinks/fibre | contact_conn flag | force-pair | floppiness (t0→tT) | spaghetti | mean invasion |
|---|---|---|---|---|---|---|
| 0.55 (300/32/280) | 1.68 | False | 0 | 0.000 → 0.000 | 0.773 | +1.89 |
| **0.62 (400/40/340 headline)** | **1.66** | False | **0** | **0.000 → 0.000** | 0.525 | +2.02 |
| 0.85 (old default) | 2.66 | True | 0 | 0.000 → 0.000 | 0.781 | +1.78 |

So **`XL_FRAC = 0.62` → 1.66 crosslinks/fibre = G4D parity**, with force-pair 0, floppiness 0, and
invasion/spaghetti indistinguishable from the old 0.85 network. Locked for the Gloria re-runs.

## 4. Re-run results (2 h, seed 23, G4D-parity crosslinks, Gloria style) — TO FILL
Driver `experiments/gloria_headlines.py` (two-phase: `sim` saves `gv_<tag>_full.npz`, `render`
draws from it, so a render tweak never re-runs a 2 h sim). New `gv_*` filenames — the old 0.85 gifs
are **not** overwritten.

| stage | condition | mean inv (µm) | max inv | detached | Rg | LCC | force-pair | grip-fails |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| R0 | force-consistent | +2.08 | 7.47 | 0.00 | 30.4 | 1.00 | 0 | 890 |
| R1 | cohesive (no EMT) | +2.08 | 7.47 | 0.00 | 30.4 | 1.00 | 0 | 890 |
| R1 | partial-EMT escape | +17.93 | 46.0 | 0.158 | 49.7 | 0.74 | 0 | 955 |
| R2 | 3 cue-front leaders | +2.21 | 5.67 | 0.00 | 30.9 | 1.00 | 0 | 1841 |
| R2 | no leaders | +2.31 | 7.85 | 0.00 | 30.8 | 1.00 | 0 | 1131 |
| R3 | energy relay (16 switches) | +18.24 | 48.9 | 0.105 | 50.8 | 0.79 | 0 | 1318 |

*(R1 at G4D-parity crosslinks reproduces the committed story: cohesive stays collective
(LCC 1.0, detached 0), partial-EMT breaks up (LCC 0.74, detached 0.16, ~8.6× invasion). force-pair 0.)*

## Reproduce
```bash
python generations/g5_organoid/experiments/gloria_headlines.py all r1   # r2 / r3 / r0
# GV_SMOKE=1 ... = fast full-pipeline check (not a headline)
```
GIFs: `output/gv_<tag>_2h.gif` (clean, no-arrow) + `output/gv_<tag>_last.png` (hero).
Network compare: `output/_compare_g4d_vs_g5_network.png`. Re-render only (no re-sim) from the saved
`gv_<tag>_full.npz`: `python .../gloria_headlines.py render <stage>`.

## Honest limitations
- Single seed (23) per condition, as with the committed headlines.
- The event microscope is **site-level** (our clutch bundles are modelled per grip site, not per
  individual molecular clutch), so it shows site-level bind/load/slip/rebind, honestly labelled.
- Crosslink parity is matched at the **headline scale (400/40/340)**; other scales need their own
  `crosslink_fraction` for 1.7/fibre (see the calibration in the driver docstring).
