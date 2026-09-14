# Generation 5 (G5) — Tumor-Organoid Invasion into Fibrous Collagen

> **Status:** proposed plan (not yet built). Drafted 2026-09-02.
> **Goal:** scale the validated single-cell motor-clutch + bead-spring-collagen core up to a
> multicellular tumor **organoid** (~150–250 µm) to study how cell–cell + cell–ECM interaction
> sets the invasion mode (confined / collective / single-cell escape).
> All numbers below carry a source; own-sim outputs are *personal testing*, not confirmed findings
> (CLAUDE.md §7.5). Literature synthesis lives in `docs/G5_research_findings.md`.
>
> **Decisions locked 2026-09-02 (Minnie):** (1) cell model = **Option A**, N motor-clutch disks +
> cell–cell adhesion (§2); (2) first milestone = **Stage B** — contractile pull → radial alignment,
> the validate-first target (§7). Stages C–F follow only after B validates. Advisor sign-off on the
> cell–cell adhesion term (§2, §10 Q1) still pending before building.

---

## 0. Target phenomenon (from Kolade/Gloria imaging + literature)

The experimental movies (4MOSC1 OSCC & PDAC organoids in reconstituted collagen, confocal
reflection + collagen probe; `memory/reference_microscopy_scale_pptx.md`) show, over Day2→Day5:

1. **Organoid** ~200–250 µm diameter, individual cells ~15–20 µm, embedded in a disordered
   collagen network.
2. **Radial leader/protrusion strands** reach out, grip collagen, and pull it inward.
3. Surrounding fibers reorganize from disordered → **radial alignment (aster / TACS-3-like)**,
   with **clearance zones** near the organoid.
4. Mostly **collective** advance with occasional **single-cell escape**.

G5 must reproduce (2)+(3) first (matches Hongbo's minimal criterion: pull → alignment), then (4).
This maps onto the project hypotheses H1 (alignment ↔ invasion mode) and H3 (alignment × EMT →
collective vs escape); CLAUDE.md §1.

---

## 1. Foundation decision — reuse the g2 core, lift g4_v2 parts, write a NEW multi-cell driver

Code audit (full report `docs/G5_code_assessment.md`) established that **all generations already
share one physics core** through subclassing:

```
CollagenConfig / Network  (generations/g2_corrected/common/model.py, validated, 92 tests)
   ├── SpheroidConfig      → g3   (single large rigid ball)
   └── G4Config → G4V2Config → g4_v2 (Numba integrator + Bell/shared-load clutch + rigid-body)
```

**The hard blocker:** the cell is a *single rigid circle* baked in everywhere (`center`,
`cell_radius`, `theta` scalars; `repulsion_forces` assumes one center; clutch state arrays are one
cell's worth). **There is no multi-cell support anywhere in the repo.** g3's "spheroid" is still
one big ball, not an assembly of cells.

**Decision:**

| Choice | Verdict |
|--------|---------|
| ECM core | **Reuse g2 `Network` / `elastic_forces` unchanged** — the single validated collagen engine. |
| Numba integrator + Bell / shared-load clutch kinetics | **Lift from `g4_v2_multiscale/model.py` as vetted building blocks** (`_advance_numba`, `bell_off_rate`, shared-load hazard). |
| Multi-cell organoid loop | **Write fresh** in `generations/g5_organoid/` — do NOT extend g4_v2's 140-line `run_clutch` monolith. |
| Fibre-free-gap network generator + protrusion idea | **Port from g3** (closest precedent for a large body with a t=0 gap). |

Rationale: keeps the validated collagen physics, avoids inheriting the g4_v2 mess the user flagged,
and isolates the genuinely new work (M cells, cell–cell contact, neighbor lists).

---

## 2. Cell representation — organoid = cluster of N motor-clutch disks + cell–cell adhesion

**Recommended (Option A):** represent the organoid as **N discrete rigid disks** (radius ~8–10 µm),
each a motor-clutch cell that grips nearby fibers via the existing Gaussian kernel, plus a
**cell–cell adhesion/repulsion** interaction between disks (a soft attractive-repulsive potential,
i.e. tissue surface tension). Cells move under overdamped rigid-body dynamics (reuse g4_v2's
released-cell branch, generalized to M centers).

Why this is the right first choice:
- **Minimal delta** from current code: the disk + Gaussian-kernel + clutch machinery already exists
  and is validated; we generalize scalars → arrays over M cells.
- **Directly produces the target science:** collective advance vs single-cell escape is an *emergent*
  outcome of the balance between cell–cell adhesion (holds cluster together) and cell–ECM traction
  (pulls cells out) — exactly the H3 axis. Ilina & Friedl 2020 (Nat Cell Biol) rank cadherin-based
  cell–cell adhesion + matrix confinement as the two 1st-order determinants of the jamming/collective
  ↔ single-cell transition — both are naturally captured by (adhesion potential) × (fiber network).
- **Extensible:** proliferation (add disks), EMT (per-cell lower adhesion / higher motility), and
  leader/follower heterogeneity all become per-cell parameters later.

Alternatives considered and deferred:
- *Single large deformable body / active-gel* — loses the individual-cell escape phenomenon; can't
  represent single-cell dissemination. Deferred.
- *Vertex / Cellular-Potts cells* (Tsingos 2023 CPM+bead-spring; Zhang/Schwarz 2024 vertex+fiber
  linkers) — richer cell shape & true confluent tissue mechanics, but a large rewrite and overkill
  for the first question. Keep as the Stage-F upgrade if disk-cluster proves too coarse.

> **Advisor check (Hongbo/Kolade):** cell–cell interaction is currently *disregarded* per CLAUDE.md
> §5 ("open question for advisors"). Adding a cell–cell adhesion potential is the central new physics
> in G5 and should be confirmed with advisors before building.

---

## 3. Scale, size, and domain

| Quantity | G4 v2 (single cell) | **G5 (organoid)** | Basis |
|----------|---------------------|-------------------|-------|
| Organoid diameter | — (one 10 µm cell) | **150–250 µm** (start ~180 µm) | Kolade/Gloria movies; ~200 µm |
| Cell diameter | 20 µm | **16–20 µm** | movies (nuclei ~15–20 µm) |
| N cells (2D disk) | 1 | **~19–37** (hex packing, ⌀ ~9 cells) | 200 µm / 20 µm ≈ 10 across → ~⌊π·5²/…⌋; start small (19) |
| Domain | 180 µm | **≥ 600 µm** (≈ organoid + 200 µm remodeling halo ×2) | force propagates ~8–10 cell diam ≈ 200 µm (Shenoy 2014; Tsingos 2023) |
| Fiber count | 99 | **~1000–3000** (density-matched to larger domain) | Lee 2014 density; scale by area |
| Beads | few ×10³ | **~10⁴–few ×10⁴** | bead spacing 1 µm × total fiber length |

**Consequence:** ~10–30× more beads/fibers and M cells. This is the performance wall — see §6.
Start 2D (matches the confocal-plane imaging and all prior generations); 3D is Stage-F.

---

## 4. Physics to add — ranked 1st-order vs deferred (from research synthesis)

Ranked by the deep-research evidence (`docs/G5_research_findings.md`; each item cites its source):

### 1st-order (build in Stages A–D)
1. **Cell–cell adhesion + tissue surface tension** — sets collective vs single-cell (Ilina & Friedl
   2020). *New.* The defining G5 physics.
2. **Force-induced radial collagen alignment** — the headline output. Emerges from many cells pulling
   + nonlinear fiber mechanics. Already partly present (10% compression → microbuckling in g2
   `elastic_forces`); **add tension-driven fiber reorientation / strain-stiffening** so the pull
   propagates ~200 µm and aligns fibers radially (Shenoy 2014; Saraswathibhatla 2025; Su/Kim 2021 —
   radial vs circumferential alignment around the *same* spheroid changes invasion).
3. **Nonlinear fiber-network mechanics (strain-stiffening + buckling)** — required for long-range
   (~8–10 cell diameter) transmission that a linear-elastic matrix cannot produce (Shenoy 2014;
   Mark 2020 eLife: radial displacement exponent shifts α=−2 → −0.2 as matrix stiffens; >20× local
   stiffening near a contractile spheroid). Partly a *parameter/force-law* change to the existing
   stretch term, not a new module.

### 2nd-order (Stage E, after A–D validate)
4. **Matrix plasticity / irreversible remodeling** — makes the radial aster *persist* instead of
   elastically recovering. Bell's-law crosslinker rupture/reformation (CLAUDE.md §5; Nam 2016
   κ≈0.82; Wisdom 2018: high vs low plasticity → 5× migration). This is the "presentable results"
   requirement from lab meeting. **Note:** SLS viscoelasticity stays *deferred* (CLAUDE.md §3a) —
   plasticity here is crosslink topology change, not per-spring SLS.
5. **Spatially heterogeneous stiffening** — local (Loxl3-type) stiffening promotes collective
   invasion; global stiffening suppresses it (Ray 2022). Model stiffness as a field, not one knob.

### Deferred (Stage F / later versions)
- MMP proteolytic matrix degradation / pore opening (Wisdom 2018) — needed for true confined
  single-cell channels; couple to plasticity later.
- Proliferation-driven internal organoid pressure (open question: does it change alignment vs
  contractile pull alone?).
- Durotaxis / explicit contact-guidance steering; EMT continuum; 3D; immune/IL-6 (v3).

---

## 5. Parameters and provenance

**Inherited (keep from g2/g4, already sourced — CLAUDE.md §8, `docs/g4_parameter_provenance.md`):**
fiber stretch κ_s=4.0×10⁻³ N/m, bending κ_b=8.27×10⁻²⁰ N·m, crosslink stiffness, bead spacing 1 µm,
Gaussian σ≈1.5–2 µm, bead drag ζ, compression ratio 0.10, and the full motor-clutch set
(N=12/site, k_c=2 nN/µm, k_on=0.055 s⁻¹, k_off0=0.018 s⁻¹, F_b=1.5 nN, v0=0.025 µm/s,
F_stall=8 nN/site; Adebowale 2021 SI Table 4 scale, Bell 1978).

**New parameters G5 introduces (each needs a source before use):**

| Symbol | Meaning | Proposed value | Source / status |
|--------|---------|----------------|-----------------|
| N_cells | cells in organoid | 19 → 37 (2D) | movies; start small for perf |
| R_cell | cell radius | 8–10 µm | movies |
| ε_cc, σ_cc | cell–cell adhesion depth / range | **TBD — fit** | tune to keep cluster cohesive at rest; Ilina 2020 qualitative |
| γ_tissue | effective tissue surface tension | **TBD — derive** from ε_cc | — |
| strain-stiffening onset / exponent | fiber tension nonlinearity | **TBD** | Steinwachs 2016 / Mark 2020 collagen fits |
| k_off,xl, F_b,xl | crosslink Bell rupture (plasticity, Stage E) | **TBD** | Nam 2016 κ≈0.82 target; check bond strains first (CLAUDE.md §5) |

> Discipline: do **not** invent ε_cc/σ_cc silently — either fit them to a resting-organoid target
> (cluster neither disperses nor collapses) or ask Kolade for a cadherin-scale estimate. Log every
> value in `parameter_provenance.md`.

---

## 6. Performance plan (the real risk)

The engine cost walls are **not** the per-step force law (it Numba-compiles fine and scales
O(E+T+links)); they are **network setup and contact detection**, which are Python O(F²) loops today:
- `build_crosslinks` nests over all fiber pairs → O(F²) with full segment×segment tests.
- `make_network_spec` retries whole networks against a percolation gate.
- `contact_patches` loops over every fiber per contact update → ×M cells for an organoid = M×F.

**Required before scale-up:**
1. **Spatial neighbor grid (cell list)** for (a) crosslink construction and (b) per-cell contact
   detection — replaces both O(F²) and M×F Python loops with O(F) / O(M·local).
2. Generalize g4_v2's `_advance_numba` **repulsion loop to M cell centers**.
3. Keep the Numba integrator; reuse g3's `np.bincount` scatter where Numba isn't used.

Ship a **perf smoke test** (Stage A gate): 2000 fibers + 30 cells must build + step at an
acceptable rate before adding kinetics.

---

## 7. Staged implementation plan (mirrors CLAUDE.md §3b discipline — validate in order)

| Stage | Deliverable | Gate criterion |
|-------|-------------|----------------|
| **G5-A** | Multi-cell scaffolding: N disks, cell–cell adhesion potential, neighbor grid, large network generator, Numba repulsion over M centers | Resting organoid is cohesive & force-balanced; perf smoke test passes (2–3k fibers, ~30 cells) |
| **G5-B** | Contractile organoid pulls collagen (all cells grip via Gaussian kernel; no cell motion yet) | **Radial alignment emerges** — radial-order metric rises; displacement halo extends ~200 µm |
| **G5-C** | Add strain-stiffening / tension-driven fiber reorientation | Displacement scaling & alignment persistence match literature targets (§8) within order-of-magnitude |
| **G5-D** | Release cells (translation+rotation) with cell–cell adhesion → **collective vs single-cell** | Sweep adhesion strength → mode transition (cohesive front ↔ escaping cells) |
| **G5-E** | Matrix plasticity (crosslink Bell rupture/reformation) so aster **persists** | κ metric > 0; aster remains after relaxing pull (vs elastic recovery control) |
| **G5-F** | Deferred richness: proteolysis, proliferation-pressure, 3D, vertex cells | only after A–E validate |

Each stage: preserve prior version (CLAUDE.md §7.2), add tests, regenerate LaTeX equation summary,
update `parameter_provenance.md`.

---

## 8. Validation targets (experimental numbers → model metrics)

These are the **quantitative anchors** for "is the model getting closer?" and the basis for ablation.
Breast/TNBC/PyMT numbers are cross-tissue **order-of-magnitude anchors**, not PDAC-exact (caveat from
research). Prefer the experimentally-measured metrics (alignment index, displacement magnitude,
contractility, plasticity %) over model-output power-law exponents.

| Metric (model) | Experimental target | Source |
|----------------|--------------------|--------|
| Radial alignment index near organoid | rises Day1→3, **plateaus** Day3→5; highest at 0–50 µm border; **persists >100 µm**; >0.2 = anisotropic; ~0.5→0.75 rise | Saraswathibhatla 2025; Ray 2022; Lee 2017 (method) |
| Radial vs tangential displacement scaling | exponents **n_r≈1.2** (radial, farther) vs **n_t≈2.2** (tangential) | Saraswathibhatla 2025 |
| Radial displacement exponent vs matrix nonlinearity | **α: −2 (linear, ~1 Pa) → −0.2 (>1000 Pa, stiffened)** | Mark 2020 eLife 51912 |
| Aggregate organoid contractility / pressure | **344±35 µN / 677±68 Pa** (TNBC ~4000-cell spheroid, 24 h, 1.2 mg/mL) | Mark 2020 |
| Near-surface strain / stiffening | ~**200 µm surface deformation (>50% strain)**, **>20× local stiffening** | Mark 2020 |
| Matrix plasticity index | **κ≈0.82** (collagen); 10–30% permanent strain at 100 Pa | Nam 2016; Wisdom 2018 |
| Plasticity → migration | high vs low plasticity → **~5× migration** | Wisdom 2018 |
| Radial vs circumferential fiber orientation → invasion | radial side → more disseminated clusters (same spheroid) | Su/Kim 2021 |

> Missing PDAC/OSCC-specific numbers (traction, alignment-rise rate, strand outgrowth µm/h,
> single-cell speed) are an **open question** — ask Kolade/Gloria whether their movies can be
> quantified (they can: they have the timelapse + collagen probe). That would give *tissue-matched*
> targets instead of breast proxies.

---

## 9. Ablation study design (mechanism on/off → metric → compare to experiment)

The strongest experimental ablations map directly to in-silico switches:

| Ablation (turn X off/on) | Model switch | Metric | Expected (experiment) |
|--------------------------|-------------|--------|-----------------------|
| **Radial vs circumferential** initial fiber orientation around same organoid | seed network anisotropy | dissemination / escaped-cell count | radial → more invasion (Su/Kim 2021) |
| **Cell–cell adhesion** high vs low | ε_cc | collective front vs single-cell escape fraction | low adhesion → single-cell (Ilina 2020) |
| **Strain-stiffening** on/off | linear vs nonlinear stretch law | displacement exponent α, alignment reach | nonlinear → far-reaching alignment (Mark 2020) |
| **Plasticity** high vs low | crosslink k_off,xl | κ; aster persistence after unload; migration | high → persistent aster + ~5× migration (Wisdom 2018) |
| **Global vs local stiffening** | stiffness field | collective invasion extent | local promotes, global suppresses (Ray 2022) |

Each ablation changes **one** switch and compares **one** metric to a literature target — clean
Sobol-style attribution of what matters (H1–H4). This is the deliverable for Kolade's lab talk.

---

## 10. Open questions for advisors

1. **Cell–cell interaction** (Hongbo/Kolade): confirm adding a cell–cell adhesion potential now — it's
   currently *disregarded* (CLAUDE.md §5). Disk-cluster vs vertex/CPM cells?
2. **PDAC/OSCC-specific targets** (Kolade/Gloria): can the existing movies be quantified for
   alignment-index rise, strand outgrowth speed (µm/h), and single-cell speed, to replace breast
   proxies?
3. **Plasticity timing** (Kolade): activate crosslink-rupture plasticity in G5-E, or keep the aster
   purely elastic for the first organoid milestone? Check bond-strain magnitudes first (Gloria's V4
   negative result: negligible weak-bond formation at 5 nN).
4. **Proliferation pressure** (Hongbo): does internal organoid growth pressure materially change
   collagen alignment vs contractile pull alone, or is it deferrable?

---

*Companion docs to generate on approval: `docs/G5_research_findings.md` (full cited literature
synthesis), `docs/G5_code_assessment.md` (engine audit), and a per-stage LaTeX equation summary.*
