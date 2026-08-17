# Validation Targets

Quantitative results from the literature that each version must reproduce. These are the
**acceptance criteria** per version (CLAUDE.md §9.2). Each entry names the target value,
the tolerance (to be set when the test is written), and the paper figure it corresponds
to. Tolerances marked *TBD* are pinned down when the corresponding test file is authored.

---

## G3 — active mechanism gates

These are implementation and causal-mechanism gates, not validation of realistic 3D migration.
Calibration seeds are 0–19; final ensemble seeds are 1000–1099 and are not used for tuning.

| Stage | Control | Gate |
|---|---|---|
| G3A | Material-point persistence | `(fiber_id, segment_id, alpha)` remains fixed until Bell unbinding |
| G3A | Internal force/torque | Relative force error < 1e-10; relative first-moment/torque error < 1e-8 |
| G3A | Rotation covariance | Rigidly rotating geometry rotates all point/projected forces identically |
| G3B | Isotropic ensemble | Ensemble polar resultant / mean individual resultant < 0.1 |
| G3B | Aligned ECM | Positive nematic order about the fibre director; +/− fraction in 0.4–0.6 |
| G3B | Rotated aligned ECM | Director follows a 30° matrix rotation within 5° |
| G3B | Feedback ablation | `beta_geometry=beta_traction=0` reduces guidance by at least 50% |
| G3C | Empty ECM | Exactly zero displacement and rotation (no hidden `v0`) |
| G3C | Asymmetric attachments | Nonzero torque; mirror control reverses its sign |
| G3C | Drag sweep | Speed changes with drag while directional symmetry conclusions remain unchanged |
| G3A unload | Elastic recovery | If pull-induced FOI signal is resolved, κ < 0.1; otherwise report `insufficient_foi_signal` |

FOI and κ follow Nam et al. 2016. For finite synthetic fixtures, κ uses the measured initial
FOI as its reference; using the ideal random value 2/π for a deliberately aligned fixture would
incorrectly label perfect elastic recovery as plasticity.

---

## v1 — fiber only (no cells)

| Test file | Target | Tolerance | Source |
|-----------|--------|-----------|--------|
| `test_negative_normal_stress.py` | Applied shear → normal stress is **negative**, comparable in magnitude to the shear stress | sign + order of magnitude, *TBD* | Saraswathibhatla 2025 Fig 4h–i |
| `test_deformation_exponents.py` | Under dipole-shear boundary forcing: radial deformation $u_r \sim r^{-1.2}$ | exponent ± *TBD* | Saraswathibhatla 2025 Fig 4p |
| `test_deformation_exponents.py` | Tangential deformation $u_t \sim r^{-2.2}$ | exponent ± *TBD* | Saraswathibhatla 2025 Fig 4p |
| `test_alignment_index_evolution.py` | Radial alignment index rises from ~0.5 to ~0.75 over ~2 simulated days | *TBD* | Saraswathibhatla 2025 Fig 4q–r |

## v1 — with cells

| Test | Target | Source |
|------|--------|--------|
| Interface velocity | Cells at the spheroid–ECM interface show tangential velocity ≫ radial | Saraswathibhatla 2025 Fig 3c, f |
| High cell–cell adhesion | Collective invasion strand forms (stiff/3 kPa HP condition) | Saraswathibhatla 2025 Fig 1m |
| Low cell–cell adhesion | Single-cell dispersion (soft/0.6 kPa HP condition) | Saraswathibhatla 2025 Fig 1m |

## v1.5 — SLS viscoelasticity

| Test | Target | Source |
|------|--------|--------|
| Single-SLS stress relaxation | $\sigma(t) = \varepsilon_0[E_\infty + E_1 e^{-t/\tau_R}]$ (matches analytical form) | SLS/Zener; Adebowale 2021 SI |
| Traction on SLS ECM | Shorter $\tau_R$ → larger deformation per unit time | Adebowale 2021 (conceptual) |

## v2 — plasticity (Bell's-law crosslinker unbinding)

| Test | Target | Source |
|------|--------|--------|
| LP condition (low $k_{-,xl}^0$) | No persistent alignment despite swirling | Saraswathibhatla 2025 (LP acIPN) |
| HP condition (high $k_{-,xl}^0$) | Persistent alignment develops | Saraswathibhatla 2025 (HP acIPN) |
| Cell-lysis / recoil | After removing cells, ECM partially recoils but retains some deformation | Saraswathibhatla 2025 Extended Fig 5e |

---

*Target numbers **verified against Saraswathibhatla 2025 Fig. 4 and its caption (2026-07-06)**,
figure image at `papers/figures/saraswathibhatla_2025_fig4.jpg`:*
- *Deformation exponents (Fig. 4p / caption): tangential n = **2.2**, radial n = **1.2**.*
- *Radial alignment index (Fig. 4r violin): median rises from **~0.5** (time=0) to **~0.75**
  (time=end), p<0.0001. Note the simulation time is in model units — the "~2 days" is the
  experimental timescale (Fig. 4l, days 1 & 3), not a figure-verified simulation duration.*
- *Negative normal stress (Fig. 4h–i): under 20% oscillatory shear, normal stress is negative
  (contractile), reaching ~−200 Pa and normalizing to ~−1 at peak shear; mechanism (Fig. 4g):
  ~50% of fibers stretch, ~50% compress.*
- *Drag law (Eq. 5) verified: ζ = 3πμ·r_c(3 + 2r₀/r_c)/5.*
