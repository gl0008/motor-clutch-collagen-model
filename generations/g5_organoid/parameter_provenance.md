# G5 parameter provenance

Units: µm, nN, s. G5 reuses the G2 collagen parameters (see repo-root
`parameter_provenance.md` and CLAUDE.md §8) and the G4 motor-clutch parameters
(`docs/g4_parameter_provenance.md`). This file records the **cell–fibre molecular clutch**
parameters added to `OrganoidConfig` for `clutch_dynamics=True`. Own-simulation output is
personal testing, not a confirmed finding (CLAUDE.md §7.5).

## Cell–fibre molecular clutch (opt-in) — inherited from G4, do not invent

| Symbol (field) | Value | Unit | Source |
|---|---|---|---|
| N clutches/site (`n_clutches_per_site`) | 12 | — | Effective bundle, not 12 integrins. Adebowale 2021 Nat Mater SI Table 4 scale; inherited from G4 (`docs/g4_parameter_provenance.md`). |
| grip sites/cell (`n_contact_sectors`) | 12 | — | Angular resolution of the cell surface; G4 `n_contact_sectors`. |
| clutch stiffness k_c (`clutch_stiffness`) | 2.0 | nN/µm | Effective clutch stiffness (not single-integrin); G2 V3 / G4 provenance. |
| bind rate k_on (`clutch_on_rate`) | 0.055 | 1/s | Adebowale 2021 SI Table 4; G4 value. |
| zero-force off-rate k_off0 (`clutch_off_rate0`) | 0.018 | 1/s | Adebowale 2021 SI Table 4; G4 value (zero-force median lifetime ≈ 38.5 s). |
| Bell force F_b (`bell_force`) | 1.5 | nN | Bell 1978 slip-bond e-fold scale; Adebowale 2021; G4 value. |
| unloaded actin speed v0 (`unloaded_actin_speed`) | 0.025 | µm/s | ≈ 24 nm/s, Adebowale 2021; Chan & Odde 2008 motor. |
| motor stall F_stall/site (`motor_stall_per_site`) | 8.0 | nN/site | Effective-site stall; Chan & Odde 2008 force–velocity scale; G4 value. |
| counter-RNG seed (`clutch_counter_seed`) | 3042 | — | Deterministic counter-addressed RNG (g4_v2). |
| mode (`clutch_mode`) | independent / shared | — | Independent vs equal-load-sharing cluster (Erdmann & Schwarz 2004; g4_v2). |

**Governing laws (all inherited verbatim from g4):**
- Motor force–velocity: `v = v0·max(0, 1 − F_site/F_stall)` — Chan & Odde 2008 Science.
- Clutch force: `F = k_c·extension`.
- Bell slip-bond off-rate: `k_off = k_off0·exp(F/F_b)` — Bell 1978.
- Shared-load hazard: `r_i = i·k_off0·exp(F_site/(i·F_b))` — Erdmann & Schwarz 2004 (g4_v2 `shared_load_hazard`).

**New G5 wiring (not a new parameter, a coupling choice):** M cells × sectors are flattened into
one "sites" axis so g4_v2's `_independent_step` / `_shared_step` apply unchanged; traction per site is
projected to collagen beads by the existing first-moment-preserving Gaussian kernel; in Stage D the
per-cell reaction uses the actual summed clutch traction. See `generations/g5_organoid/model.py`.

## G5D v2 leader / EMT-like heterogeneity — modelling assumptions

| Field | Value used | Status / provenance |
|---|---:|---|
| `leader_fraction` | default 0; 0.10–0.20 in exploratory sweeps; 0.12 in the three-mode comparison | A deterministic modelling knob: the outermost cells are selected. The range is not calibrated to a measured leader fraction. |
| `leader_adhesion_factor` | 0.15 | Coarse proxy for partial-EMT adhesion loss. Qualitatively motivated by Kalluri & Weinberg 2009 and Aiello et al. 2018; the numerical multiplier is a personal-test assumption, not a fitted E-cadherin measurement. |
| ordinary `cc_adhesion` in the three-mode comparison | 60 / 8 / 1 nN/µm | Stuck / collective / single-cell demonstration settings. They map a model phase space and are not direct cadherin-force measurements. |
| total pull in the three-mode comparison | 6 / 45 / 45 nN | Stuck / collective / single-cell demonstration settings. Own-simulation inputs, not yet calibrated to the Kolade-lab organoid movies. |

For a pair of cells, adhesion is scaled by the weaker member:
`k_adh,ij = cc_adhesion × min(s_i, s_j)`, where `s_i = 0.15` for a
selected leader and `1` otherwise. Repulsion is unchanged. A current leader has
lower adhesion only; increased leader traction, MMP secretion, and a continuous
EMT phenotype remain outside this version.
