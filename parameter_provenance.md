# ECM parameter provenance — current Gloria G2 baseline

This document distinguishes the **parameters currently used by the corrected
Generation 2 (G2) Gloria model** from the Saraswathibhatla acIPN experimental
condition that motivates a future calibration.  Values labelled *not modelled*
must not be presented as inputs to the current simulation.

All G2 lengths are in micrometres (µm), forces in nanonewtons (nN), and time
in seconds.  The baseline is a planar, discrete collagen network pulled by an
MDA-MB-231 cell; it is not yet a 3-D acIPN tumour-spheroid model.

## 1. Collagen concentration and bulk stiffness

| Quantity | Current G2 setting | Saraswathibhatla acIPN reference / intended future target | Status |
|---|---:|---:|---|
| Type-I collagen concentration | Not modelled | 1.5–1.6 mg/mL | Not an input to G2 |
| Alginate concentration | Not modelled | 4.8 mg/mL | Not an input to G2 |
| Bulk storage modulus, G′ | Not computed or calibrated | 0.6 / 1.2 / 3 kPa | Requires virtual rheology calibration |
| Pure-collagen comparison range | Not modelled | 1–4 mg/mL | Literature context only |

Consequently, G2 cannot currently be described as a 1.5 mg/mL acIPN model or
as representing any of the three acIPN stiffness conditions.  Its fibre and
link mechanics must first be calibrated against a bulk shear/tensile test.

## 2. Fibre geometry and discretisation

| Parameter | Current G2 setting | Saraswathibhatla / literature reference | Status |
|---|---:|---:|---|
| Effective fibre diameter | 0.30 µm (300 nm) | 100 nm in Saraswathibhatla; reported measurements span 62–325 nm | Used |
| Fibre contour length | 20–80 µm | 8–100+ µm reported; Lee used 8 µm in its model | Used |
| Nominal bead spacing | 0.75 µm | 1.0 µm | Used |
| Beads per fibre | Variable, resampled by contour length; saved seed: 7,173 beads / 99 fibres = 72.5 on average | Lee: 5; previous project note: 3 | Used; neither 3 nor 5 |
| Fibre centreline | Straight-to-curved random segments, curvature amplitude 0.9 µm | Not a Saraswathibhatla SI input | Used |

The 0.75 µm value is a numerical segment spacing, not the physical collagen
fibre diameter.  G2 uses a much finer chain representation than a five-bead
fibre.

## 3. Fibre elasticity

| Quantity | Current G2 setting | Saraswathibhatla / literature reference | Status |
|---|---:|---:|---|
| Single-fibre Young's modulus, E | 32 MPa | 32 MPa simulation value; reported range roughly 30–800 MPa | Used |
| Fibre axial rule | `k = EA / l0` | Direct spring `κ_s,f = 4.0×10⁻³ N/m` | Used, but not the same parametrisation |
| Nominal axial spring, at `l0 = 0.75 µm` | 3,016 nN/µm = 3.02 N/m | 0.004 N/m | Derived current value |
| Compression/tension ratio | 0.10 | Not specified in supplied acIPN table | Used microbuckling approximation |
| Beam bending rigidity, EI | 12.72 nN·µm² = `1.27×10⁻²⁰ N·m²` | `κ_b,f = 8.27×10⁻²⁰ N·m` | Used; do not compare directly without matching force laws |
| Bulk G′ | Not computed | 1 mg/mL: ~3 Pa; 2 mg/mL: ~44–200 Pa; 3 mg/mL: ~97–550 Pa | Not modelled |

The nominal G2 axial spring follows from its 32 MPa modulus, 300 nm effective
diameter, and 0.75 µm segment length.  It is approximately 754 times the
literal 0.004 N/m Saraswathibhatla spring value.  This is a material-model
difference that requires calibration, not a value that can be exchanged alone.

## 4. Bead–spring and crosslink parameters

| Parameter | Current G2 setting | Saraswathibhatla SI Table 2 reference | Status |
|---|---:|---:|---|
| Fibre segment rest length | Each resampled segment, nominally 0.75 µm | 1.0 µm | Used |
| Fibre diameter | 300 nm | 100 nm | Used |
| Fibre stretch constant | Derived as `EA/l0`, nominally 3.02 N/m | `4.0×10⁻³ N/m` | Used, different formulation/value |
| Fibre bending law | Discrete beam `EI/l0³` | `8.27×10⁻²⁰ N·m` | Used, different convention |
| Permanent crosslink location | Exact intersection of two fibre segments | 20 nm arm | Used; arm length is effectively zero |
| Crosslink diameter | Not modelled | 10 nm | Not modelled |
| Permanent crosslink stiffness | 75 nN/µm = 0.075 N/m | `2.0×10⁻³ N/m` | Used |
| Crosslink binding-site spacing | Not modelled | One site per 100 nm | Not modelled |
| Crosslink type | Permanent, freely hinged, stress-free at the initial intersection | Explicit crosslinker arm | Different model |
| V4 newly formed weak-link stiffness | 15 nN/µm = 0.015 N/m | Not a supplied Saraswathibhatla baseline parameter | V4 only |

`crosslink_modulus_mpa = 0.40` remains in the G2 configuration for provenance,
but it is not used by the force calculation; permanent links use
`crosslink_stiffness = 75 nN/µm`.

## 5. Network density and geometry

| Parameter | Current G2 setting | Saraswathibhatla / Lee reference | Status |
|---|---:|---:|---|
| Simulation domain | 180 × 180 µm planar square | Annulus: inner diameter 50 µm, outer diameter 200 µm, height 1 µm | Different geometry |
| Number of fibres | 99 | 3-D density ~2.85 fibres/µm³ | Used |
| 2-D areal density | 0.00306 fibres/µm² (99 / 180²) | No direct equivalent | Derived; no thickness is defined |
| Crosslinks in saved seed 17 | 383 permanent intersection links | Crosslinker/fibre ratio ~5.72 | Used snapshot |
| Links per fibre in saved seed 17 | 3.87 (383 / 99) | ~5.72 | Derived snapshot |
| Crosslink count rule | Generated from intersections; changes with network geometry/seed | Ratio-controlled crosslinker model | Different rule |
| Boundary condition | Only beads within a 2.5 µm outer band are fixed | Not the acIPN annular boundary condition | Used |

The current model is 2-D.  A 3-D fibre density and a crosslinker/fibre ratio
from the acIPN network therefore cannot be transferred directly without
choosing an effective thickness and a network-generation rule.

## 6. Cell or tumour-spheroid size

| Quantity | Current G2 setting | acIPN tumour-spheroid reference | Status |
|---|---:|---:|---|
| Biological object | One MDA-MB-231 cell | 250-cell aggregated tumour spheroid | Different biology |
| Radius | 9 µm | ~50 µm (day-0 diameter ~100 µm) | Used, not spheroid-equivalent |
| Initial inner boundary | None; a circular cell exclusion/pull surface | Annular inner diameter 50 µm | Not modelled |
| Active pull | 5 nN total in V2/V4 (2.5/5/10 nN sensitivity); motor–clutch rules in V3 | Not the supplied spheroid force protocol | Used |

## 7. Displacement-field validation

| Validation quantity | Current G2 state | Saraswathibhatla target | Status |
|---|---|---|---|
| Radial displacement fit | Not implemented | `u_radial ∼ r⁻¹·²` | Pending |
| Tangential displacement fit | Not implemented | `u_tangential ∼ r⁻²·²` | Pending |
| Current spatial outputs | Mean displacement and alignment in near/intermediate/far shells | Full radial/tangential displacement profiles | Partial only |

The power-law exponents are appropriate future validation targets, but the
current G2 validation does not calculate or fit them.  They must not be marked
as passed.

## Code provenance

- Base geometry, material parameters, loading, and integration defaults:
  `generations/g2_corrected/common/model.py`, `CollagenConfig`.
- Axial and bending conversions: `CollagenConfig.axial_rigidity` and
  `CollagenConfig.bending_rigidity` in the same file.
- Intersection-only permanent crosslink construction:
  `build_crosslinks` in the same file.
- Saved network counts for seed 17:
  `generations/g2_corrected/validation_summary.json`.
- V3 motor–clutch overrides:
  `generations/g2_corrected/v3_two_sided_migration/model.py`,
  `MigrationConfig`.
- V4 weak-link overrides:
  `generations/g2_corrected/v4_contact_plasticity/model.py`,
  `PlasticityConfig`.

The Saraswathibhatla, Lee, Picu, Stein, and Licup values in this document are
research references supplied for model comparison.  They are not automatically
active G2 inputs unless explicitly listed as a current G2 setting above.
