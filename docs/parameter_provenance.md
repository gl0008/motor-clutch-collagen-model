# Parameter Provenance

**Single source of truth for every numerical parameter used in simulation.**
(CLAUDE.md §6 and §12.1.)

Rules:
- Every parameter used anywhere in `src/` gets an entry here. **No exceptions.**
- If a value is "just a guess," that must be stated explicitly, with the range considered.
- If the value in code differs from the value here, one of them is wrong — reconcile immediately.
- Default to Saraswathibhatla et al. 2025 parameters when a value is ambiguous, until an advisor overrides (CLAUDE.md §14).

## Entry format

```markdown
## <symbol> — <short name>
- **Value**: <value> <unit>
- **Source**: <paper, table/figure/equation> OR "guess — see note"
- **Physical meaning**: <one line>
- **Considered range**: <lo> to <hi> <unit>
- **Sensitivity**: <tested? planned sweep? not yet tested>
- **Related parameters**: <symbols>
```

## ⚠️ Provenance-integrity note (read before using these values)

The v1 ECM values below were transcribed from **Saraswathibhatla et al. 2025,
Supplementary Table 2**, via automated text-extraction of the PDF on 2026-07-06 (the table
layout was scrambled by extraction; the symbol↔value mapping was reconstructed and 7 values
were cross-checked against the running text — high confidence, but **not** yet checked
against the rendered table image). **Action before any production run:** open the
Supplementary Table 2 PDF and verify each value, especially the drag-law prefactors (Eq. 5)
and the Bell's-law exponent form (Eq. 9), which were partly garbled in extraction. Tracked
in `logs/open_questions.md`.

Status by version:
- **v1 (elastic)** — populated below from Saraswathibhatla Suppl. Table 2 + Methods.
- **v1.5 (SLS)** — Adebowale 2021 *reference* values recorded, flagged "to be mapped": their
  values are for a motor-clutch cell-adhesion SLS in pN/nm, not yet mapped onto our ECM
  fiber-spring SLS (a v1.5 modeling decision).
- **v2 (plasticity)** — the Bell's-law crosslinker parameters (`k0_off_xl`, `lambda_xl`)
  already appear under ECM crosslinkers; v2 sweeps `k0_off_xl` (hypothesis H4).
- **Cell agents (v1)** — project extension **not** present in Saraswathibhatla's fiber-only
  model; listed as pending, no direct source.

---

# v1 — Geometry

## r_inner — annular inner radius (spheroid–matrix interface)
- **Value**: 5.0e-5 m (50 μm)
- **Source**: Saraswathibhatla et al. 2025, Methods "Theoretical model"; Suppl. Table 2. See `docs/design_decisions/D001_geometry_annular.md`.
- **Physical meaning**: Radius of the inner circular boundary representing the tumor spheroid surface, where traction enters the ECM.
- **Considered range**: 50 μm fixed for apples-to-apples reproduction; spheroid radius could later vary 30–100 μm.
- **Sensitivity**: Not yet tested. Held fixed in v1.
- **Related parameters**: `r_outer`, boundary parameters (`kappa_s_bnd`, `kappa_r_bnd`, `k_on_bnd`).

## r_outer — annular outer radius (far field)
- **Value**: 2.0e-4 m (200 μm)
- **Source**: Saraswathibhatla et al. 2025, Methods; Suppl. Table 2. D001.
- **Physical meaning**: Radius of the outer circular boundary (far field), held fixed.
- **Considered range**: 200 μm fixed (matches reference). Must exceed the deformation decay length.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `r_inner`.

## h_slab — slab height (quasi-2D interpretation)
- **Value**: 1.0e-6 m (1 μm) in the reference; **our v1 is genuinely 2D**.
- **Source**: Saraswathibhatla et al. 2025, Methods (1 μm-thick annular slab = spheroid equatorial plane).
- **Physical meaning**: Out-of-plane thickness of the reference's quasi-2D slab. In our 2D model there is no z-extent; this value only matters for converting the reference's per-μm³ densities to 2D (see `C_f`).
- **Considered range**: n/a (2D).
- **Sensitivity**: n/a. Revisit if we go to 3D (D002).
- **Related parameters**: `C_f`, `R_xl`.

---

# v1 — ECM fibers

## kappa_s_f — fiber extensional (Hookean) stiffness
- **Value**: 4.0e-3 N/m
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2. Force from harmonic potential U_s = ½ κ_s (r − r0)² (their Eq. 7).
- **Physical meaning**: Spring constant of one ~1 μm fiber segment resisting stretch.
- **Considered range**: 1e-3 to 1e-2 N/m (order of magnitude for reconstituted col1).
- **Sensitivity**: Not yet tested. Candidate for a mechanics sweep.
- **Related parameters**: `kappa_b_f`, `kappa_s_xl`, `r0_f`.

## kappa_b_f — fiber bending stiffness
- **Value**: 8.27e-20 N·m  (≈ 20 k_BT)
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2. Harmonic bending potential U_b = ½ κ_b (θ − θ0)² (their Eq. 8).
- **Physical meaning**: Resistance of a three-bead fiber segment to bending away from the equilibrium angle θ0_f = 0.
- **Considered range**: as published; sets the fiber persistence length.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `kappa_s_f`, `theta0_f`.

## r0_f — fiber segment (bead spacing) rest length
- **Value**: 1.0e-6 m (1 μm)
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2.
- **Physical meaning**: Equilibrium length of one cylindrical fiber element (bead–bead spacing).
- **Considered range**: as published (also the coarse-graining length scale).
- **Sensitivity**: Not yet tested; sets spatial resolution.
- **Related parameters**: `rc_f`, `kappa_s_f`.

## rc_f — fiber diameter
- **Value**: 1.0e-7 m (100 nm)
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2.
- **Physical meaning**: Cylindrical fiber-element diameter; enters the drag coefficient (their Eq. 5) and steric/crosslink spacing.
- **Considered range**: as published.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `zeta` (drag), `r0_f`.

## theta0_f — fiber equilibrium bending angle
- **Value**: 0 rad
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2 (adjacent fiber elements prefer to be collinear).
- **Physical meaning**: Rest angle in the bending potential; 0 → straight fibers are unstressed.
- **Considered range**: fixed at 0.
- **Sensitivity**: n/a.
- **Related parameters**: `kappa_b_f`.

---

# v1 — ECM crosslinkers

## kappa_s_xl — crosslinker extensional stiffness
- **Value**: 2.0e-3 N/m
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2.
- **Physical meaning**: Spring constant of a crosslinker arm (harmonic, their Eq. 7).
- **Considered range**: as published (~half of `kappa_s_f`).
- **Sensitivity**: Not yet tested.
- **Related parameters**: `kappa_s_f`, `r0_xl`, `k0_off_xl`.

## r0_xl — crosslinker arm length
- **Value**: 2.0e-8 m (20 nm)
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2.
- **Physical meaning**: Equilibrium length of each of the two crosslinker arms joined at their midpoint.
- **Considered range**: as published.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `rc_xl`, `kappa_s_xl`.

## rc_xl — crosslinker arm diameter
- **Value**: 1.0e-8 m (10 nm)
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2.
- **Physical meaning**: Crosslinker arm diameter (enters drag).
- **Considered range**: as published.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `r0_xl`.

## k_on_xl — crosslinker binding rate
- **Value**: 1.0e2 M⁻¹·s⁻¹
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2. Binding to sites spaced every 100 nm along fibers.
- **Physical meaning**: On-rate for a crosslinker arm to bind an available fiber site.
- **Considered range**: as published.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `k0_off_xl`, `R_xl`.

## k0_off_xl — crosslinker zero-force unbinding rate  (k⁰₋,xl)
- **Value**: 1.0e-6 s⁻¹
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2. Bell's-law base off-rate (their Eq. 9).
- **Physical meaning**: Unloaded crosslinker dissociation rate; force accelerates it (Bell's law). Low value → near-permanent crosslinks (elastic-like) in v1.
- **Considered range**: **This is the v2 plasticity sweep parameter (hypothesis H4)** — vary to move between "no persistent alignment" (low) and "persistent alignment" (high). Range TBD from Saraswathibhatla LP/HP discussion.
- **Sensitivity**: The key plasticity knob. Held at 1e-6 s⁻¹ (near-elastic) in v1; swept in v2.
- **Related parameters**: `lambda_xl`, `kappa_s_xl`, `k_on_xl`.

## lambda_xl — Bell's-law force sensitivity (bond length)  (λ₋,xl)
- **Value**: 4.0e-10 m
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2; Bell 1978. Appears in the Eq. 9 exponent.
- **Physical meaning**: Characteristic bond length setting how strongly tensile force accelerates unbinding.
- **Considered range**: as published. ⚠️ The exact exponent form in Eq. 9 was garbled in extraction — verify the dimensionless argument (F·λ/k_BT vs the extracted F·λ·κ_s/k_BT) against the paper before coding v2.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `k0_off_xl`, `kappa_s_xl`, `kBT`.

---

# v1 — Network composition

## C_f — fiber density
- **Value**: ~2.85 fibers/μm³  (reference's 1 μm-thick 3D slab)
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2 (approximate; derived).
- **Physical meaning**: Number density of fibers seeded during self-assembly.
- **Considered range**: n/a yet. ⚠️ **2D conversion needed:** for our genuinely-2D model the per-μm³ value must be converted using the 1 μm slab height (`h_slab`) to a per-μm² areal density. Pin this down in `network_init.py`.
- **Sensitivity**: Not yet tested; controls network connectivity/percolation.
- **Related parameters**: `R_xl`, `mean_L_f`, `h_slab`.

## R_xl — crosslinker-to-fiber ratio
- **Value**: ~5.72 crosslinkers per fiber
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2 (approximate).
- **Physical meaning**: Sets network crosslink density → elastic modulus and rigidity-percolation state.
- **Considered range**: as published.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `C_f`, `k_on_xl`.

## mean_L_f — average fiber length
- **Value**: ~1 μm
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2 (approximate; fibers are short, ~1 segment).
- **Physical meaning**: Mean length of assembled fibers (nucleation + polymerization, no depolymerization).
- **Considered range**: as published.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `r0_f`, `C_f`.

---

# v1 — Network initialization & alignment seeding (D008)

## S_alignment — target 2D nematic order (H1 control)
- **Value**: 0.0 default (swept)
- **Source**: Modeling choice (D008). Not a Saraswathibhatla parameter — it is the H1 sweep axis.
- **Physical meaning**: Degree of initial fiber alignment about the director; S = ⟨cos 2Δφ⟩, 0 = isotropic, →1 = fully aligned.
- **Considered range**: [0, 1) — the primary H1 sweep.
- **Sensitivity**: This IS the H1 factor; swept by design.
- **Related parameters**: `director`.

## director — alignment reference field
- **Value**: `radial` default
- **Source**: Modeling choice (D008); radial is the H1-relevant case (radial alignment → collective invasion, Saraswathibhatla).
- **Physical meaning**: Field the nematic order is measured/seeded against (`radial` | `circumferential` | `horizontal`).
- **Considered range**: the three named options.
- **Sensitivity**: Not swept in the primary H1 study (fixed radial).
- **Related parameters**: `S_alignment`.

## n_fibers — number of fibers (dev)
- **Value**: 3000 (explicit dev value)
- **Source**: Modeling choice for tractable development. Full reference density via `n_fibers_from_density(C_f, h_slab, r_inner, r_outer)` ≈ 3.4×10⁵.
- **Physical meaning**: Fiber count placed in the annulus.
- **Considered range**: 10²–10³ (dev) up to ~3.4×10⁵ (reference density; production).
- **Sensitivity**: Affects network connectivity; production runs should approach reference density (optimize first, §12.8).
- **Related parameters**: `C_f`, `h_slab`, `R_xl`.

## n_beads_per_fiber — beads per fiber
- **Value**: 3  (**modeling choice / to revisit**)
- **Source**: Guess for v1. Saraswathibhatla's *average* fiber is ~1 segment (2 beads); we use 3 so bending is exercised. Not matched to their fiber-length distribution yet.
- **Physical meaning**: Chain length of each fiber (≥3 gives an internal bending angle).
- **Considered range**: 2 (min) to ~10; revisit against the paper's length distribution.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `r0_f`, `mean_L_f`.

## xl_cutoff — crosslinker formation distance
- **Value**: 5.0e-6 m  (**modeling choice / to revisit**)
- **Source**: Guess for the v1 single-spring crosslinker model (D008). Saraswathibhatla's two-arm crosslinker binds sites ~100 nm apart; our cutoff is larger because it acts on a sparse dev network.
- **Physical meaning**: Max bead–bead distance (different fibers) at which a crosslinker spring forms.
- **Considered range**: ~1e-7 to ~5e-6 m depending on network density.
- **Sensitivity**: Controls crosslink count with `R_xl`; not yet tested.
- **Related parameters**: `R_xl`, `kappa_s_xl`, `r0_xl`.

## max_place_attempts — fiber placement retries
- **Value**: 20
- **Source**: Numerical (rejection sampling to fit a fiber inside the annulus). Not physical.
- **Physical meaning**: Retry budget before a fiber is placed as-is.
- **Considered range**: n/a.
- **Sensitivity**: n/a (large enough that fibers fit given fiber length ≪ annulus width).
- **Related parameters**: `r_inner`, `r_outer`, `n_beads_per_fiber`, `r0_f`.

---

# v1 — Langevin (overdamped) dynamics

## eta_medium — effective medium viscosity
- **Value**: 8.6 kg/(m·s)  (= 8.6 Pa·s)
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2.
- **Physical meaning**: Effective viscosity of the surrounding medium; enters the per-element drag coefficient via **Eq. 5 (verified)**: ζ = 3π·μ·r_c·(3 + 2·r₀/r_c)/5 ≈ 3.73×10⁻⁵ kg/s for fiber elements (see `integrator.bead_drag_coefficient`).
- **Considered range**: ⚠️ This is **~8600× water** — an *effective* coarse-grained value, **not** literal solvent viscosity. Do not "correct" it to water's value.
- **Sensitivity**: Sets the overall timescale; not a physics sweep target.
- **Related parameters**: `zeta` (derived), `dt`, `rc_f`, `rc_xl`.

## dt — integration time step
- **Value**: 4.0e-4 s
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2 and Eq. 6 (forward Euler).
- **Physical meaning**: Time step of the overdamped Euler integrator.
- **Considered range**: must satisfy overdamped stability (ζ/κ_max ≫ dt). Verify a stability estimate when the integrator design-decision is written.
- **Sensitivity**: Numerical, not physical. Check convergence by halving dt.
- **Related parameters**: `eta_medium`, `kappa_s_f`, integrator choice (open question).

## kBT — thermal energy
- **Value**: 4.142e-21 J  (≈ 300 K)
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2; fluctuation–dissipation (their Eq. 4).
- **Physical meaning**: Thermal noise amplitude in the Langevin equation.
- **Considered range**: fixed (room/physiological temperature).
- **Sensitivity**: n/a.
- **Related parameters**: `eta_medium`, `lambda_xl`, `kappa_b_f`.

---

# v1 — Boundary conditions

## kappa_s_bnd — boundary-binding spring stiffness
- **Value**: 4.0e-3 N/m
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2.
- **Physical meaning**: Spring tethering a fiber endpoint that has bound to a wall (within 1 μm).
- **Considered range**: as published (= `kappa_s_f`).
- **Sensitivity**: Not yet tested.
- **Related parameters**: `k_on_bnd`, `kappa_r_bnd`.

## kappa_r_bnd — boundary repulsion strength
- **Value**: 1.69e-3 N/m
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2.
- **Physical meaning**: Strength of the repulsive normal force walls exert on overlapping matrix elements (∝ overlap distance).
- **Considered range**: as published.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `kappa_s_bnd`.

## k_on_bnd — fiber→boundary binding rate
- **Value**: 1.0e3 M⁻¹·s⁻¹
- **Source**: Saraswathibhatla et al. 2025, Suppl. Table 2.
- **Physical meaning**: Rate at which a fiber endpoint within 1 μm of a wall binds irreversibly to it.
- **Considered range**: as published.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `kappa_s_bnd`, `binding_distance`.

## binding_distance — wall-binding capture distance
- **Value**: 1.0e-6 m (1 μm)
- **Source**: Saraswathibhatla et al. 2025, Methods ("within 1 μm of either boundary bind").
- **Physical meaning**: A fiber endpoint within this distance of the inner or outer wall binds and is tethered by a `kappa_s_bnd` spring (`boundary.py`).
- **Considered range**: as published.
- **Sensitivity**: Not yet tested; sets how many endpoints anchor.
- **Related parameters**: `kappa_s_bnd`, `k_on_bnd`, `r_inner`, `r_outer`.

---

# v1 — Dipole traction boundary condition (fiber-only validation)

*Used to reproduce Saraswathibhatla's fiber-only results before cell agents are coupled.
In the full model this is replaced by explicit cell-agent traction (see D001 divergence note).*

## bc_rotation_rate — inner-anchor rotation increment
- **Value**: 0.001 degree per time step
- **Source**: Saraswathibhatla et al. 2025, Methods (~"dipole" boundary scheme).
- **Physical meaning**: Angular step by which each inner-boundary anchor rotates (CW/CCW) to shear the matrix, until its fiber reaches the tension threshold.
- **Considered range**: as published.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `bc_tension_threshold`, `dt`.

## bc_tension_threshold — anchor stop-rotating tension
- **Value**: 1.0e-10 N (100 pN)
- **Source**: Saraswathibhatla et al. 2025, Methods.
- **Physical meaning**: Tensile force on the attached fiber at which an inner anchor stops rotating.
- **Considered range**: as published.
- **Sensitivity**: Not yet tested.
- **Related parameters**: `bc_rotation_rate`.

## bc_dipole_rule — rotation-direction pattern
- **Value**: For an anchor at initial angle φ (CCW from +x): take φ mod 60°; rotate **clockwise if remainder < 30°, counterclockwise if ≥ 30°**.
- **Source**: Saraswathibhatla et al. 2025, Methods. Only this "dipole" rule (vs stochastic/uniform controls) reproduced the experimental deformation profile (their Extended Fig. 9b–e).
- **Physical meaning**: Creates alternating shear "dipoles" around the interface, mimicking tangentially swirling cells.
- **Considered range**: dipole (chosen) vs stochastic vs uniform (controls).
- **Sensitivity**: Qualitatively decisive (see source).
- **Related parameters**: `bc_rotation_rate`, `bc_tension_threshold`.

---

# v1 — Cell agents  (PENDING — project extension, no direct Saraswathibhatla source)

Saraswathibhatla's published model has **no explicit cell agents** (traction is the boundary
condition above). The following are needed for our cell-agent architecture and have **no
direct source yet** — values to be set from other references (motor-clutch traction from
Adebowale 2021; adhesion/jamming from Ilina 2020) or as explicit, documented guesses.

- [ ] `traction_dipole_mag` — magnitude of the contractile dipole a cell exerts on the ECM. *Candidate source: Adebowale 2021 traction ~100 Pa upper bound; Saraswathibhatla 100 pN anchor threshold.*
- [ ] `adhesion_strength` — cell–cell adhesion (EMT-modulated); sets collective vs single-cell mode. *Candidate source: Ilina 2020 jamming framework.*
- [ ] `polarity_relax` / `guidance_gain` — polarity update + contact-guidance coupling to local fiber orientation.
- [ ] EMT dynamics rate constants for de/dt.
- [ ] `dbscan_eps` — cluster-detection radius for analysis = 20 μm (CLAUDE.md §11; Saraswathibhatla Fig. 1o–p). *(analysis parameter, not a model force)*

---

# v1.5 — Viscoelasticity (SLS)  (REFERENCE VALUES — to be mapped onto ECM springs)

From Adebowale et al. 2021 (motor-clutch cell-adhesion SLS, in pN/nm). **Not yet adopted:**
their SLS sits under a cell adhesion, not our fiber springs. Mapping onto the ECM
bead-spring SLS is a v1.5 modeling decision (design-decision file D004 forthcoming).

Adebowale SLS relations: **E₀ = kₐ + k_l** (instantaneous), **E_∞ = k_l** (residual),
**E₁ = kₐ** (relaxing), **τ_R = η/kₐ** (relaxation time). Relaxation modulus
E(t) = E_∞ + (E₀ − E_∞) e^(−t/τ_R).

- Reference `k_a` (additional/relaxing stiffness) ≈ 0.9–1.0 pN/nm (Adebowale Suppl. Tables 3–4).
- Reference `k_l` (long-term stiffness) ≈ 0.1 pN/nm.
- Reference `eta_dashpot` ≈ 1–100 pN·s/nm.
- **`tau_R` (relaxation time)** — Adebowale sweep **1–1,000 s**; representative 1 s (fast), 10 s (medium), 100 s (slow). **This is the v1.5 sweep parameter (hypothesis H2).**
- Experimental context: initial substrate modulus held ≈ 2 kPa; stress-relaxation half-times ~100 s (fast) / ~240 s (medium) / ~2,200 s (slow). Tumor relevance (one example): pancreatic transformation lowers ECM relaxation time (~93 s → 66 s; Rubiano 2018 via Adebowale) — a case where faster relaxation is expected to enhance migration; flag for Kolade.

---

# v2 — Plasticity

Plasticity in the ECM is realized through the **existing** Bell's-law crosslinker
parameters above (`k0_off_xl`, `lambda_xl`): v2 sweeps `k0_off_xl` (hypothesis H4) to move
between transient and persistent radial alignment. See Nam 2016 and Wisdom 2018 for the
viscoplasticity rationale (design-decision file forthcoming). No new numerical parameters
introduced beyond those already listed, pending the v2 design.

---

# G3 — Active emergent-guidance model

G3 uses SI internally. Published pN/nm motor-clutch values are converted once at config load.
Provisional feedback and drag choices are not biological measurements and must be swept rather
than presented as predictions.

| Parameter | Active value | Provenance | Status |
|---|---:|---|---|
| `cell_radius` | 10 µm | Minimal rigid-cell geometry agreed in the Stage C/E design | Coarse-grained assumption |
| `fibre_length` | 40 µm | Few-fibre numerical fixture; within Lee 2014's broad 20–200 µm range but not concentration/image calibrated | Numerical fixture only |
| `bead_spacing` | 1 µm | Saraswathibhatla et al. 2025 SI Table 2 | Published |
| `kappa_s_f` | 4.0×10⁻³ N/m | Saraswathibhatla et al. 2025 SI Table 2 | Published |
| `kappa_b_f` | 8.27×10⁻²⁰ N·m | Saraswathibhatla et al. 2025 SI Table 2 | Published |
| `gaussian_sigma` | 2 µm | Advisor-approved G3 projection scale; inherited from the v1.0b coupling test | Modeling choice |
| `capture_distance` | 2 µm | Set equal to the Gaussian width for the minimal attachment test | Provisional; sweep 1–3 µm |
| `n_clutches` | 200 | Adebowale et al. 2021 SI Table 4 | Published one-module count; shared across active protrusions |
| `n_motors` | 200 | Adebowale et al. 2021 SI Table 4 | Published one-module count; shared across active protrusions |
| `motor_force` | 2 pN/motor | Adebowale et al. 2021 SI Table 4 | Published |
| `bell_force` | 2 pN | Adebowale et al. 2021 SI Table 4 | Published |
| `bind_rate` | 0.2 s⁻¹ | Adebowale et al. 2021 SI Table 4 | Published |
| `unbind_rate` | 0.02 s⁻¹ | Adebowale et al. 2021 SI Table 4 | Published |
| `clutch_stiffness` | 5 pN/nm = 5×10⁻³ N/m | Adebowale et al. 2021 SI Table 4 | Published |
| `unloaded_actin_speed` | 24 nm/s | Adebowale et al. 2021 SI Table 4 | Published |
| `dt` | 0.005 s | Adebowale et al. 2021 SI Table 4 motor-clutch update | Published |
| `ecm_substeps` | 2 (0.0025 s each) | Explicit bead-chain stability requires Δt_ECM < ζ/(2κ_s) | Numerical stability choice |
| `n_sectors` | 24 | 15° angular discretization so the 30° covariance control maps exactly | Numerical resolution choice |
| `n_active_protrusions` | 2 | Minimal low-protrusion-count hypothesis informed qualitatively by Fraley et al. 2010 | Coarse-grained assumption |
| `protrusion_lifetime` | 120 s | No transferable MDA-MB-231 kinetic constant in Carey 2016 | Provisional; sweep 60/120/300 s |
| `beta_geometry` | 2 | Carey 2016 motivates the feedback sign, not this coefficient | Provisional; sweep 0/1/2/4 |
| `beta_traction` | 2 | Coarse-grained adhesion/traction reinforcement | Provisional; sweep 0/1/2/4 |
| `feedback_time` | 30 s | Numerical smoothing timescale | Provisional sensitivity parameter |
| `cell_drag` | 300 nN·s/µm | Calibrated G2-scale reference, not an independent prediction | Sweep 150/300/600 nN·s/µm |
| `rotational_drag` | `cell_drag × radius²` | Dimensional rigid-body closure | Provisional; sweep factor 0.5/1/2 |
| `duration_g3a` | 15 s | Long enough to resolve FOI in the tangential fixture; ~23 s reaches the no-steric-overlap guard | Numerical scope limit |

**G3 exclusions:** no SLS, transient crosslinks, permanent plasticity, WLC, 3D concentration
mapping, nucleus, MMP, EMT, or prescribed self-propulsion velocity.
