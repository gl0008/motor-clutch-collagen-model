# G3-R fresh mechanism validation — 2026-08-18

> Personal simulation / mechanism validation; not yet a realistic 3D tumor-migration
> prediction.

## Decision

The revised G3A/B/C implementation passes its predeclared **mechanism and execution gates**
on a crosslinked 99-fibre network. It is cleared for visualization and the next calibration
round. It is not cleared for a realistic migration claim or a statistically supported claim
that aligned collagen guides the new cell-intrinsic polarity rule.

## What changed

- Every displayed G3A/B/C condition uses 99 fibres, 20–80 µm fibre lengths, 0.75 µm bead
  spacing and Generation-2 crosslink/network mechanics.
- Protrusions now have explicit continuous lengths and tips. Clutches attach at a tip-near
  collagen material point, rather than at an abstract cell-surface sector.
- A 4 s minimum lifetime plus relative/absolute hysteresis prevents instantaneous active-sector
  relabeling. Replaced protrusions retract instead of teleporting to the new direction.
- Polarity is a conserved noisy activity field on 24 membrane sectors. Its symmetry breaking
  is cell intrinsic; collagen geometry is not read before contact. Attachment/traction can
  reinforce a sector only after physical coupling.
- G3C translation and rotation remain equal-and-opposite reaction driven. No prescribed
  migration velocity is present.

The revised coarse-grained runs use a 10× motor/Bell/clutch-force scale relative to the old
minimal fixture. This is an explicit calibration assumption, not a measured molecular value.

## Network audit

| Quantity | Fresh network |
|---|---:|
| Fibres | 99 |
| Beads | 7,089 |
| Crosslinks | 334 |
| Boundary-connected fraction | 98 / 99 = 0.9899 |
| Cell-contact fibres connected to boundary | yes |

## Fresh results

| Stage/control | Result | Interpretation |
|---|---:|---|
| G3A first attachment | 5.0 s | explicit protrusion growth precedes 0.5 µm-range capture |
| G3A maximum protrusion length | 2.418 µm | tip visibly reaches the contact fibre |
| G3A maximum bound clutches | 9 / 200 | attachment is sparse, not whole-surface locking |
| G3A bead displacement, motor on | 0.001171 µm | transmitted network deformation |
| G3A bead displacement, motor off | 1.38×10⁻¹⁴ µm | matched no-motor control |
| G3A maximum contact-fibre angle change | 0.199° | real sub-degree reorientation; GIF vectors are display-only ×1000 |
| G3B feedback-ON final activity axis | 30.02° | matched 20 s single-run diagnostic |
| G3B feedback-OFF final activity axis | 30.02° | no OFF/ON divergence in this run |
| G3B aligned, 4 seeds | nematic order 0.0142 | execution/symmetry-breaking check only |
| G3B isotropic, 4 seeds | nematic order 0.0160 | effectively indistinguishable at n=4 |
| G3C fixed cell | 0 µm | matched fixed control |
| G3C released cell | 0.003122 µm, 0.01675° | reaction-driven motion is nonzero but tiny |
| G3C released empty network | 0 µm | no hidden self-propulsion |

## Gate result

All nine automated gates passed: completion, 99-fibre scale, boundary connectivity,
growth-before-attachment, motor-on deformation above motor-off, positive (but very small)
aligned nematic order, valid isotropic execution, zero hidden motion and nonzero released
reaction motion.

Passing these gates does **not** mean the revised guidance hypothesis is validated. The matched
OFF/ON pair follows the same sector sequence and finishes at the same axis; with only four
seeds per ensemble, aligned and isotropic nematic orders are nearly identical. The
next scientific gate must therefore compare aligned versus isotropic confidence intervals on
a preregistered larger ensemble, repeat a rotated-axis covariance condition, and sweep the
polarity/adhesion gains before any aligned-guidance claim.

## Reproduction

```powershell
$env:PYTHONPATH='src'
py -3 run_g3_revision_experiments.py --ensemble-seeds 4 --ensemble-duration 15
```

Machine-readable results are in `results/g3_revision/validation_summary.json`. The notebook uses
Gloria's G2 precomputed-data/SVG architecture (`g3-web-data.js` + `g3-gloria.js`): there are no
numerical axes or grids, and the renderer interpolates material-point crosslinks and solved bead,
cell, protrusion and clutch states. Synchronized downloadable GIFs from the same saved frames are
in `figures/g3_revision/`.
