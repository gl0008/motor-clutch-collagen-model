# G4 v4E — contact-guided persistent single-cell migration

## Purpose

G4D showed a fixed-radius single cell moving under unbalanced collagen traction, but its path wandered because every newly encountered direction was equivalent. G4 v4E asks one narrower question:

> If one side of a single MDA-MB-231 cell encounters locally radial collagen, and a collagen-aligned probing protrusion is allowed to mature into a temporary front, does wandering become persistent migration?

G4 v4A–D remain frozen. This directory adds four matched conditions on top of the accepted G4D collagen, permanent elastic crosslinks, shared-load motor–clutch sites, Bell rupture, and overdamped cell translation.

| Condition | Matrix | Protrusion rule | Question isolated |
|---|---|---|---|
| E0 `d_control` | Random | G4D relocation | What does the inherited model do without a directional cue? |
| E1 `matrix_cue` | One-sided radial tract | G4D relocation | Is matrix geometry alone sufficient? |
| E2 `protrusion_guidance` | Random | Contact selection + temporary front memory | Does a cell-internal persistence rule work without a world-axis bias? |
| E3 `combined` | One-sided radial tract | Contact selection + temporary front memory | Do the environmental cue and cell response reinforce each other? |

## Cumulative mechanics retained from G4D

Every collagen bead remains overdamped:

$$
\zeta_b\dot{\mathbf r}_i=
\mathbf F_i^{\mathrm{stretch}}+
\mathbf F_i^{\mathrm{bend}}+
\mathbf F_i^{\mathrm{crosslink}}+
\mathbf F_i^{\mathrm{contact}}+
\mathbf F_i^{\mathrm{steric}}.
$$

The fixed-radius cell receives the equal-and-opposite reaction:

$$
\zeta_c\dot{\mathbf r}_c=-\sum_s\mathbf F_s^{\mathrm{cell\to ECM}}.
$$

Shared-load clutch clusters keep the G4C/D loading and Bell rupture equations. No force constant is increased in E1–E3.

## What v4E adds

Eligible protrusions attach only to real continuous collagen material points within the contact band plus a provisional 5 µm probing reach. Their normalized selection probability is

$$
P_s=
\frac{I_sw_s\exp\!\left[
\beta_a(\hat{\mathbf t}_s\!\cdot\!\hat{\mathbf e}_s)^2+
\beta_m\hat{\mathbf p}\!\cdot\!\hat{\mathbf e}_s
\right]}
{\sum_q I_qw_q\exp\!\left[
\beta_a(\hat{\mathbf t}_q\!\cdot\!\hat{\mathbf e}_q)^2+
\beta_m\hat{\mathbf p}\!\cdot\!\hat{\mathbf e}_q
\right]}.
$$

The first probing contact that remains continuously bound for 120 s establishes the temporary front $\hat{\mathbf p}$. This is not a prescribed $+x$ polarization: E2 has no matrix cue, and its front direction must be distributed across world coordinates over seeds. Other motor–clutch sites remain adhesion anchors. Rear sites are never cut. After a rear site reaches zero bound clutches through the original Bell law, its next candidate is restricted to the forward half-plane. If no front-cone site remains bound for 120 s, memory clears.

## Primary measurements

$$
D_{\mathrm{cue}}=[\mathbf r_c(T)-\mathbf r_c(0)]\cdot\hat{\mathbf e}_{\mathrm{cue}},
\qquad
P=\frac{D_{\mathrm{net}}}{L_{\mathrm{path}}},
$$

$$
C_t=\left\langle
|\hat{\mathbf v}_c\cdot\hat{\mathbf t}_{\mathrm{contact}}|
\right\rangle.
$$

The predeclared endpoint is the 20-seed E3 confidence interval, compared with E0–E2. Negative results are retained.

## Full-clock result (20 seeds, 6 h)

| Condition | $D_{cue}$, µm (mean [95% CI]) | Persistence (mean [95% CI]) | $C_t$ (mean [95% CI]) |
|---|---:|---:|---:|
| E0 | −9.642 [−19.709, 0.426] | 0.104 [0.047, 0.161] | 0.635 [0.613, 0.657] |
| E1 | −2.764 [−8.422, 2.895] | 0.085 [0.057, 0.113] | 0.627 [0.606, 0.647] |
| E2 | −2.206 [−15.109, 10.696] | 0.247 [0.174, 0.320] | 0.635 [0.621, 0.648] |
| E3 | −3.904 [−13.060, 5.252] | 0.197 [0.156, 0.238] | 0.610 [0.589, 0.632] |

The control cue-axis CI contains zero, and E2 front directions are not concentrated on a fixed world axis. The primary E3 gate is **not supported**: its cue-axis CI includes zero and its mean is negative. Protrusion guidance increases mean path persistence in E2 and E3 relative to E0, but adding the present one-sided tract does not turn that persistence into reliable motion toward the tract. This is retained as a mechanism-level negative result; force, seed, and endpoints were not changed after inspection.

## Declared short-clock audit

The 10 min numerical/coordinate-bias audit passes timestep direction and early front-event ordering at $\Delta t=0.05$ versus $0.025$ s. It does **not** pass cue rotation: the 90° cue gives −5.684 µm projected displacement. It also does **not** pass the domain gate: near-cell $\Delta S_r$ changes by 90.6% between the density-matched 180 and 240 µm fields. Radius sensitivity changes direction at $R=12\,\mu m$. The declared 15°/30° tract-width and 60/120/180 s maturation checks are similar for this seed, whereas the $\beta_a=\beta_m$ sweep is non-monotonic. Therefore v4E remains exploratory even apart from its negative primary endpoint. This short audit identifies where the mechanism is not robust; it is not a substitute for full-clock sensitivity.

The website separates three observation clocks. The ECM overview stores one frame every 10 min for 6 h; protrusion telemetry stores one frame every 30 s; and a five-minute event microscope stores one frame every second beginning at the first complete site failure or front-establishment event. The event view shows every bound effective clutch, complete site failures, returns from a zero-bound state, and site force. These are sampled views of the same Python trajectory, not separate simulations.

## Run

```bash
python -m generations.g4_v4e_contact_guided_single_cell.build_demo --quick
python -m generations.g4_v4e_contact_guided_single_cell.build_demo
python -m generations.g4_v4e_contact_guided_single_cell.validate_model --duration 600
pytest generations/g4_v4e_contact_guided_single_cell/tests
```

The full build uses 6 h at $\Delta t=0.05$ s. The validation command is a declared short-clock numerical and coordinate-bias audit; it does not replace the 6 h biological comparison and cannot by itself validate migration. The browser only plays Python-generated states; it contains no second physics implementation.
