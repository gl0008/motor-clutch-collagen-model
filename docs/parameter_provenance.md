# G3 parameter provenance

G3 uses SI units internally. Published pN/nm motor-clutch values are converted once when the
configuration is loaded. Provisional feedback and drag values are assumptions to sweep, not
measured biological constants or independent predictions.

| Parameter | Active value | Provenance | Status |
|---|---:|---|---|
| `cell_radius` | 10 µm | Minimal rigid-cell fixture | Coarse-grained assumption |
| `fibre_count` | 8 | Few-fibre debugging fixture | Numerical choice |
| `fibre_length` | 40 µm | Within the broad Lee et al. 2014 fibre scale; not image calibrated | Numerical fixture |
| `bead_spacing` | 1 µm | Saraswathibhatla et al. 2025 SI Table 2 | Published discretization |
| `kappa_s_f` | 4.0×10⁻³ N/m | Saraswathibhatla et al. 2025 SI Table 2 | Published |
| `kappa_b_f` | 8.27×10⁻²⁰ N·m | Saraswathibhatla et al. 2025 SI Table 2 | Published |
| `gaussian_sigma` | 2 µm | Approved local projection scale | Modeling choice |
| `capture_distance` | 2 µm | Equal to the projection width in the minimal attachment test | Provisional; sweep 1–3 µm |
| `n_clutches` | 200 | Adebowale et al. 2021 SI Table 4 | Published one-module count |
| `n_motors` | 200 | Adebowale et al. 2021 SI Table 4 | Published one-module count |
| `motor_force` | 2 pN/motor | Adebowale et al. 2021 SI Table 4 | Published |
| `bell_force` | 2 pN | Adebowale et al. 2021 SI Table 4 | Published |
| `bind_rate` | 0.2 s⁻¹ | Adebowale et al. 2021 SI Table 4 | Published |
| `unbind_rate` | 0.02 s⁻¹ | Adebowale et al. 2021 SI Table 4 | Published |
| `clutch_stiffness` | 5 pN/nm | Adebowale et al. 2021 SI Table 4 | Published |
| `unloaded_actin_speed` | 24 nm/s | Adebowale et al. 2021 SI Table 4 | Published |
| `dt` | 0.005 s | Adebowale et al. 2021 SI Table 4 | Published clutch update step |
| `ecm_substeps` | 2 | Required by the explicit bead-chain stability estimate | Numerical safeguard |
| `n_sectors` | 24 | 15° sectors resolve the preregistered 30° covariance test exactly | Numerical choice |
| `n_active_protrusions` | 2 | Minimal competing-protrusion representation | Coarse-grained assumption |
| `protrusion_lifetime` | 120 s | Carey et al. 2016 motivates feedback, not this rate | Sweep 60/120/300 s |
| `beta_geometry` | 2 | Coarse-grained geometry reinforcement | Sweep 0/1/2/4 |
| `beta_traction` | 2 | Coarse-grained traction reinforcement | Sweep 0/1/2/4 |
| `feedback_time` | 30 s | EMA smoothing timescale | Provisional |
| `cell_drag` | 300 nN·s/µm | G2-scale calibration reference | Sweep 150/300/600 nN·s/µm |
| `rotational_drag` | `cell_drag × radius²` | Dimensional rigid-body closure | Sweep factor 0.5/1/2 |
| `duration_g3a` | 15 s | Resolves FOI before the no-steric-overlap guard is reached | Numerical scope limit |

## Explicit exclusions

G3 contains no SLS, transient crosslinks, permanent plasticity, WLC, 3D collagen,
concentration/pore calibration, nucleus, deformable cortex, MMP, EMT, cell-cell interaction,
CAF, interstitial flow, or prescribed self-propulsion velocity.
