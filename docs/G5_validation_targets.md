# G5 Validation Targets — radial collagen alignment around invading tumours

> Compiled 2026-09-03 from two deep-research passes (2nd pass heavily server-rate-limited) + our own
> structure-tensor analysis of Kolade/Gloria's imaging. Marks each number **verified** (3-0 adversarial)
> vs **unverified** (extracted but the verifier step was rate-limited — a lead, not disproven).
> Own-simulation output is *personal testing*, not a confirmed finding (CLAUDE.md §7.5).

## 0. The metric matches the field (good)

The published per-fibre metric is **cos²(angle between the radial vector and the fibre)** → 0 tangential,
1 radial, **0.5 = random/isotropic** (PMC4773852, verified). Our model/analysis uses
**S = 2·cos²θ − 1** → −1 tangential, +1 radial, **0 = isotropic**. They are the same quantity rescaled
(their 0.5 ↔ our 0), so our radial-order numbers are **directly comparable** to the literature and to the
Kolade imaging.

## 1. Radial alignment vs distance — model vs imaging vs published

| Quantity | Published | Kolade imaging (ours, approx.) | G5 model |
|---|---|---|---|
| Right at the surface | **tangential** (TACS-2 ~0°; spheroid parallel ~45°) ✅ PLoS ONE 0156442; Provenzano/Conklin | slightly tangential at centre (−0.16) | config-dependent |
| Peak radial alignment | **radial** (TACS-3 ~90°) ✅; unverified spheroid peak ~0.8 (4T1) | **+0.28 @ ~50 µm** | near-field **+0.2 … +0.4** |
| Reach (persists above baseline) | **100 µm → 2.65 mm** (~5× radius) ✅ PMC4773852 | positive to **~130 µm** | **~100 µm** (mid-field) |
| Isotropic baseline | 0.5 (cos²) / 0 (our S); unverified 4T1 baseline ~0.71 | ~0 | 0 |

**Reading:** near-field peak magnitude and the tangential-surface→radial-farther shape **agree** across
all three. The model's **reach is short** (~100 µm) vs published (100 µm–mm). That is the main gap.

## 2. Matrix displacement / strain (model is weak here)

| Quantity | Published (verified, Mark 2020 eLife 51912) | G5 model |
|---|---|---|
| Deformation reach | **~200 µm** from surface at 24 h | µm-scale |
| Near-surface strain | **>50 %** | small |
| Displacement power-law exponent | **α = −2 (linear, ~1 Pa) → −0.2 (>1000 Pa, strain-stiffening)** | not yet measured |

**Reading:** the model's collagen barely displaces vs the real ~200 µm / >50 % strain. The published
α = −2 → −0.2 shift is the signature of **strain-stiffening** carrying force far — the mechanism our
Stage C should supply.

## 3. Tissue-matched numbers — PDAC & oral-squamous (VERIFIED, 2nd pass)

**PDAC (Ray et al. 2022, JCI Insight art. 150330 — verified 3-0):**
- TACS-3 (fibers **~90°±30°** to the ductal boundary) present in **~40 % of stage-I ductal structures**.
- Early PanIN: **~28 % of partially-extruded cells perpendicular (TACS-3) vs ~53 % parallel (TACS-2)**.
- Quantified by **CurveAlign / CT-FIRE on SHG**; aligned collagen guides dissemination by **FAK-dependent
  contact guidance** (FAK inhibition cuts single-cell extrusion + liver mets).

**PDAC anisotropy (Hamilton et al. 2022, Front. Oncol. PMC9623060 — verified 3-0):**
- 40 patients / 544 SHG images: cancer collagen has higher anisotropy factor Fa than normal at small
  scales (**~2.5–3 µm**); all tissue types converge at **~21 µm**. *(Fa is a WTMM anisotropy factor, not a
  0–1 alignment coefficient, and gives no radial-vs-distance curve.)*

**Oral-squamous / HNSCC spheroids (Kim/Ahn et al. 2019, Acta Biomater. 84:280, PMID 30500449 —
verified 3-0) — the closest experimental analog to G5, same tissue class as Kolade's 4MOSC1:**
- Highly-invasive HNSCC spheroids (OECM-1, SAS), EGF-stimulated, **generate strong contraction
  deformation of the surrounding collagen in the very early stage and align the fibres RADIALLY about
  the spheroid centre.**
- **The magnitude of that initial contraction positively correlates with the extent of subsequent
  invasion.** Radial alignment happens **BEFORE cells invade**, via direct contraction/pulling.
  *(Qualitative — no 0–1 coefficient, no µm/h speed reported.)*

*Dropped:* the breast-spheroid values (4T1 peak ~0.8, baseline ~0.71, reach ~400/250 µm) **failed
verification** in both passes — do not use them.

## 4. Mechanism: swirling vs direct pulling (RESOLVED — supports our model)

**Prof. Kolade reportedly does not observe the "swirling" of Saraswathibhatla et al. 2025** (Minnie's
recollection — still worth confirming with Kolade). The published evidence now **directly supports
treating swirling as system-specific, not required**:

- **Kim/Ahn 2019 (HNSCC — the tissue class matching Kolade's):** radial alignment is produced by
  **direct spheroid contraction/pulling, BEFORE cells invade**, with no swirling invoked. Direct radial
  traction is thus an **independently sufficient** mechanism (verified).
- Classic breast work (Riching 2014, Provenzano/Conklin) attributes alignment to **Rho/ROCK actomyosin
  contractility** and defines/prognosticates TACS-3 **without swirling**.
- Swirling (Saraswathibhatla 2025, breast) is one recent, system-specific route.

**Decision for the model:** do **not** add a swirling/tangential-motility term — G5's mechanism (surface
cells directly reel collagen inward = contractile radial pull) is exactly the **HNSCC-verified**
mechanism for our tissue class. The real gaps — short reach and small displacement — are closed with
**strain-stiffening (Stage C)** and **plasticity (Stage E)** (better force transmission + persistence),
not a new motion mode. A testable prediction we already reproduce qualitatively: **contraction magnitude
↔ invasion extent** (Kim/Ahn 2019). See `memory/project_swirling_caveat.md`.

## 5. Concrete targets to aim the model at

1. Near-field radial order **≈ +0.25–0.3 peaking ~50 µm** outside the organoid (matches Kolade + TACS-3). ✅ already ~met.
2. Radial alignment **persisting to ≥100 µm**, ideally further — **needs strain-stiffening** to extend reach.
3. Collagen displacement reaching **~200 µm with high near-surface strain** — currently far too small.
4. Displacement-vs-distance exponent moving from **≈ −2 toward ≈ −0.2** as stiffening engages — a clean
   Stage-C ablation target.

## Caveats

Second research pass was server-rate-limited (synthesis + most verification failed), so §3 and the
spheroid peak/baseline/reach numbers are **leads to re-verify**, not settled. Kolade imaging is a single
2D image, structure-tensor estimate, auto-calibrated scale — approximate. Re-run the focused literature
search when the API is not rate-limiting to confirm the PDAC/HNSCC-specific values.

## Sources
Provenzano 2006 / Conklin 2011 (TACS-3) · PLoS ONE 0156442 (spheroid orientation, OrientationJ) ·
PMC4773852 (cos² metric, 2.65 mm reach) · Mark 2020 eLife 51912 (200 µm / >50 % strain, α = −2→−0.2) ·
Riching 2014 PMC4255204 (contractility mechanism, prestrain alignment) · Saraswathibhatla 2025 (swirling —
treat as system-specific).
