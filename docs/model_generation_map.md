# Model generation & version map

> Companion to [`model_generation_map.drawio`](model_generation_map.drawio) (open in the VS Code
> *Draw.io Integration* extension or app.diagrams.net). Same content; this text file is for quick
> reading and for **appending meeting feedback**.
>
> The diagram is a **tree** with two lanes so each advisor can read their own column:
> 🟩 **left = BIOLOGY** (what the cell/organoid does — for **Kolade**) ·
> 🟦 **right = PHYSICS** (the equations — for **Hongbo**).
> The central **spine** is the lineage; each arrow is the limitation that motivates the next version.
> Nothing is overwritten — every version is preserved ([`VERSION_MAP.md`](../VERSION_MAP.md)).

Status: 🟦 frozen · 🟨 current · 🟥 planned/prototype

```
G1 ──▶ G2 ──▶ G3 ──▶ G4 v1 ──▶ G4 v2 ──▶ G5
                     └── (both G4 versions preserved) ──┘
```

---

## G1 · V0–V4 — development archive 🟦 `frozen`

- **🟩 Biology (Kolade):** the question forming. How does a cell mechanically talk to collagen?
  Prototyped in order — SLS relaxation → a few single fibres → first fibre–fibre crosslinks →
  first moving cell → first "does the matrix stay deformed after unloading?" test. *Kept only to
  show how the idea evolved — not results (CLAUDE.md §7.5).*
- **🟦 Physics (Hongbo):** no single active equation set (archived prototypes). SLS was inserted
  *before* elastic collagen + crosslink transmission were validated → corrected in G2.
- **⚠ Limitation → next:** correct BCs, units, coupling; drop premature SLS.

**Feedback log:** _(archive — usually none)_

---

## G2 · V2–V4 — corrected elastic baseline (the shared engine, 92 tests) 🟦 `frozen`

- **🟩 Biology (Kolade):** one cell grips collagen and the pull spreads — a local 5 nN pull reaches
  distant fibres *only* through fibre–fibre crosslinks (21.6% larger response with links). A released
  cell then migrates (0.28 µm/min). Plasticity tested and honestly **not** seen at 5 nN.
- **🟦 Physics (Hongbo):** overdamped bead dynamics

  $$\zeta\,\dot{\mathbf r}_i = \mathbf F_i^{\text{stretch}} + \mathbf F_i^{\text{bend}} + \mathbf F_i^{\text{clutch}}$$

  Elastic stretch (Hooke) + standard elastic bend. **No SLS, no plasticity.**
- **⚠ Limitation → next:** migration direction is still prescribed (imposed imbalance).

**Feedback log:**

---

## G3 · spheroid — contractile spheroid → radial collagen aster 🟦 `frozen`

- **🟩 Biology (Kolade):** reel collagen into a radial aster. A round spheroid, not touching collagen
  at t=0, sends protrusions out in *all* directions, grips fibres and reels them inward → disordered
  tangle becomes a radial pattern (TACS-3-like; radial order −0.49 → −0.04). Symmetric, no prescribed
  side; the point is matrix remodelling, not cell motion.
- **🟦 Physics (Hongbo):**
  - material point: $\mathbf x_j=(1-\alpha_j)\mathbf r_a+\alpha_j\mathbf r_b$
  - Bell off-rate (clutch): $r_{\text{off}}=r_{\text{off}}^{0}\exp(|\mathbf F_j|/F_b)$
  - Gaussian projection: $w_{ij}\propto\exp(-d^2/2\sigma^2),\ \sum_i w_{ij}=1$
  - rigid body: $\gamma_c\,\dot{\mathbf r}_c=-\sum_j\mathbf F_j$
  - contact penalty: $\mathbf F^{\text{contact}}=k\,\delta\,(\mathbf r_i-\mathbf r_c)/d$
- **⚠ Limitation → next:** calibrate stiffness; prove the indirect path; release the cell.

**Feedback log:**

---

## G4 v1 · A–D — interactive calibration (short-time) 🟦 `frozen`

- **🟩 Biology (Kolade):** calibrate the single cell in four steps — A tune matrix stiffness ·
  B show a non-contact fibre moves *only* via crosslinks · C clutches stochastically slip & let go ·
  D release the cell so it can crawl with no prescribed direction.
- **🟦 Physics (Hongbo):** same overdamped core + Bell clutch kinetics; **12 independent** clutches
  per material point. (2–5 min horizon → changes came out sub-pixel.)
- **⚠ Limitation → next:** observe 2 h; resolve individual clutch failure.

**Feedback log:**

---

## G4 v2 · A–D — long-time + visible clutch failure 🟨 **CURRENT**

- **🟩 Biology (Kolade):** watch the same cell for 2 h (like the real imaging). Individual clutches
  visibly load, slip and fully detach; fibres recoil. Compare 12 *independent* clutches vs a
  *shared-load* team that fails as a cascade. Fixed / moving / mobile-ECM controls.
- **🟦 Physics (Hongbo):** shared load on $i$ of $N$ bound clutches at one site

  $$f_i=\frac{F_{\text{site}}}{i},\qquad r_i=i\,k_{\text{off}}^{0}\exp\!\Big(\frac{F_{\text{site}}}{i\,F_b}\Big),\qquad g_i=(N-i)\,k_{\text{on}}$$
- **⚠ Limitation → next:** still one rigid cell → scale to a multicellular organoid.

**Feedback log:** _(date / who / point)_

---

## G5 · organoid — multicellular tumour organoid invasion 🟥 **PLANNED / PROTOTYPE**

- **🟩 Biology (Kolade):** a whole organoid, not one cell — ~200 µm, ~19–43 motor–clutch cells held
  together by **cell–cell adhesion**. Together they pull the collagen radial, then **invade** outward.
  The balance of cell–cell adhesion vs cell–matrix traction decides **collective advance vs
  single-cell escape** (grounded in Kolade/Gloria organoid movies).
- **🟦 Physics (Hongbo):**
  - force balance per bead: $\zeta\dot{\mathbf r}_i=\mathbf F^{\text{stretch}}+\mathbf F^{\text{bend}}+\mathbf F^{\text{crosslink}}+\mathbf F^{\text{repulsion}}+\mathbf F^{\text{active}}$
  - cell $m$ reaction (invasion): $\gamma_c\,\dot{\mathbf r}_c^{(m)}=-\sum_j\mathbf F_j^{(m)}+\sum \mathbf F_{cc}$
  - **new — cell–cell:** $\mathbf F_{cc}=$ soft attract–repel $U(|\mathbf r_a-\mathbf r_b|)$
  - **Stage E plasticity:** crosslink Bell rupture + re-weld
- **⚠ Open questions:** advisor sign-off on the cell–cell term; performance wall (10–30× beads);
  2D only; SLS/plasticity still deferred.

**Feedback log:**
- **Advisor sign-off PENDING — Q1:** exact cell–cell adhesion form? (CLAUDE.md §5)
- **Q2:** swirling not seen in Kolade lab — don't build G5 on the Saraswathibhatla swirling mechanism.
- _(date / who / point)_

---

## Table 1 · Two force couplings — does either break?

Clears up the "which spring ruptures?" confusion (CLAUDE.md §5: **Bell-on-clutch ≠ Bell-on-crosslinker**).

| Coupling | Breaks? | Governed by |
|----------|---------|-------------|
| **Cell–fibre clutch** (cell grips collagen) | **No** (default): constant 12 nN | Prescribed traction (`total_pull_force`). Opt-in slip mode: Bell $k_{\text{off}}^{0}e^{F/F_b}$, $k_{\text{off}}^{0}=0.018$/s, $F_b=1.5$ nN |
| **Crosslinker** (fibre–fibre weld) | **Yes** — Stage E (`plasticity=True`) | Bell rupture $k_{\text{off}}^{0}e^{F/F_b}$, $k_{\text{off}}^{0}=2\times10^{-4}$/s, $F_b=3.0$ nN; re-welds at current crossings → irreversible strain |
| **Fibre segment** (bead–bead spring) | **No** | Elastic Hooke $F=\kappa_s\cdot\text{ext}$, $\kappa_s=4\times10^{-3}$ N/m |

## Table 2 · Generation status at a glance

| Gen | What it is | Status |
|-----|-----------|--------|
| G1 | V0–V4 history archive | frozen |
| G2 | corrected elastic baseline (engine) | frozen |
| G3 | spheroid → radial aster | frozen |
| G4 v1 | single-cell calibration | frozen |
| G4 v2 | 2 h + visible clutch failure | **CURRENT** |
| G5 | organoid remodel + invade | planned / prototype |

## Table 3 · Stage A–F roadmap (orthogonal to G1–G5)

The mechanics-first staging *inside* the current focus (CLAUDE.md §3b) — not separate models.

| Stage | Adds | Gate |
|-------|------|------|
| A | minimal fibre network | forces balance at rest |
| B | overdamped dynamics | alignment from pull alone |
| C | fixed rigid cell + Gaussian | material-point + FOI response |
| D | moving attachment sites | attachment tracks deforming fibre |
| E | mobile rigid cell | cell migrates; reactions balance |
| F | deformable cell / richer biology | only after A–E validate |

---

## How to add a professor's comment

1. **Diagram:** type into the yellow strip under the relevant generation (`date / who / point`).
   Kolade reads the green (left) lane; Hongbo reads the blue (right) lane.
2. **This file:** add a bullet under that generation's **Feedback log**.
3. **New generation?** In the diagram copy a band, recolour the trunk, extend the spine arrow;
   here, copy a section block. Then follow the "rule for adding a future model" in
   [`VERSION_MAP.md`](../VERSION_MAP.md).
