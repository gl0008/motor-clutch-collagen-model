# G4 parameter provenance and tuning boundaries

G4 uses µm, nN and s. Values labelled **starting choice** are not literature
measurements; they are exposed as controlled sensitivities or explicitly left
for experimental fitting. A paper can support an equation or an order of
magnitude without validating the exact G4 number.

## Geometry and elastic ECM

| Parameter | Baseline / exact website values | Provenance | Status |
|---|---:|---|---|
| domain | 180 × 180 µm | accepted G2 boundary-size check | inherited numerical fixture |
| cell radius | 10 µm | G2/G3 single-cell-scale fixture | replace with experiment-specific segmentation |
| fibres | 99 finite fibres | accepted G2 density fixture | image-derived topology remains future work |
| fibre length | 20–80 µm | inside Lee et al. (2014) 20–200 µm effective fibre range | distribution not fitted |
| effective diameter | 0.30 µm | Lee et al. fibre-scale 0.20–0.35 µm | coarse-grained fibre, not a 1.5 nm fibril |
| bead spacing | 1.0 µm | computation/visual resolution | convergence still required |
| collagen modulus | 3 MPa | active G3 softened reorientation fixture | declared starting choice |
| bending multiplier | 0.10, **0.25**, 0.50, 1.0 | professor requested lower bending stiffness; separated from axial EA | one-axis G4A sensitivity |
| compression ratio | 0.10 | G2 microbuckling baseline; Notbohm et al. motivates soft compression | exact ratio not fitted |
| crosslink probability | 0, 0.20, **0.35**, 0.60 | professor requested not every intersection be linked | topology sensitivity, not measured density |
| crosslink stiffness | 5, **10**, 25, 75 nN/µm | spans softened G3 to stiff G2 link penalty | compliance sensitivity, not molecular stiffness |
| bead drag | 90, **180**, 360 nN·s/µm | brackets G2/G3 overdamped response times | fit to time-resolved displacement |
| pull | 12, **24**, 48 nN total scalar force | tens-of-nN cell/spheroid traction context (Steinwachs et al.) | mechanism range, not force fit |
| direct contact shell | 0–3 µm | professor-approved hybrid coupling | fit to protrusion/contact microscopy |
| Gaussian width | 1.5 µm | professor-approved local weighting choice | applies only after hard contact eligibility |

G4B uses 48 nN, `p_x=0.60`, `k_x=25 nN/µm` and 240 s as a declared
stronger-visibility validation case. Its job is to ask whether graph-class
motion is detectable, not to claim that this combination is already a fitted
tumour microenvironment.

## Effective motor-clutch and rigid cell

| Parameter | G4C/D default | Provenance | Interpretation boundary |
|---|---:|---|---|
| clutches per site | 12 | G2 V3 effective-clutch calibration | coarse-grained bundle |
| sites | 12 symmetric sectors | matches the 12 all-around direct candidates | numerical angular resolution |
| clutch stiffness `k_c` | 2 nN/µm | inherited G2 V3 | not single-integrin stiffness |
| on-rate `k_on` | 0.055 s⁻¹ | inherited G2 V3 | sweep against bound-fraction/dwell data |
| zero-force off-rate `k_off0` | 0.018 s⁻¹ | inherited G2 V3; near Adebowale et al. effective rate scale | zero-force median = 38.5 s |
| Bell force `F_b` | 1.5 nN | inherited G2 V3 effective-clutch scale | e-fold hazard landmark, not hard threshold |
| unloaded actin speed `v0` | 0.025 µm/s | G2 V3; close to Adebowale et al. 24 nm/s scale | site-level loading speed |
| stall force per site | 8 nN | G2 V3 effective-site calibration | not single-motor stall force |
| cell drag | 600 nN·s/µm | upper end of the preregistered 150/300/600 sweep; 12-site reference gives ~0.37 µm/min path speed | ensemble calibration still required |
| rotational drag | `cell_drag × R²` | dimensional rigid-body closure | prefactor not fitted |

Bell (1978) supports `k_off(F)=k_off0 exp(F/F_b)`. Bangasser & Odde
(2013) supports stochastic motor-clutch load-and-fail coupling. Prahl et al.
(2020) supports reaction-driven motion from spatially opposed modules. None of
those papers makes the exact G4 coarse-grained values a molecular fit.

## Derived slip landmarks, not thresholds

For `k_c=2 nN/µm`, `k_off0=0.018 s⁻¹`, and `F_b=1.5 nN`:

| constant force | extension | median lifetime `ln(2)/k_off(F)` |
|---:|---:|---:|
| 0 nN | 0 µm | 38.5 s |
| 1.5 nN | 0.75 µm | 14.2 s |
| 3.0 nN | 1.50 µm | 5.2 s |

A clutch can detach below these forces or survive above them. The table only
shows how the probability distribution shifts as load increases.
