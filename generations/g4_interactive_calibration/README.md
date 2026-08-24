# Generation 4 · interactive elastic calibration → slippage → moving cell

G4 is a cumulative experiment, not four unrelated animations:

1. **G4A — tune the basic mechanics.** A rigid cell is fixed. A prescribed,
   symmetric pull acts only on collagen material points in a 3 µm contact shell;
   Gaussian weights distribute the finite total force among those eligible
   points. The website varies one coefficient at a time: pull, bending rigidity,
   bead drag, crosslink probability, and crosslink stiffness.
2. **G4B — verify indirect realignment.** The cell remains fixed. Direct fibres
   receive the active force; indirect fibres receive **zero direct Gaussian
   force**. Their only mechanical path is
   `direct fibre → permanent elastic crosslink → indirect fibre`. The output
   labels each fibre by crosslink-graph distance (0, 1, 2+, or unconnected) and
   reports displacement and radial order by class.
3. **G4C — add motor–clutch slippage.** The accepted elastic network is unchanged
   and the cell remains fixed. Each surface site has 12 stochastic effective
   clutches. Force builds through clutch extension and rupture follows Bell's
   continuous force-dependent hazard; fibre recoil is measured after traction
   drops.
4. **G4D — release the cell.** This is **everything in G4C plus** rigid-body
   translation and rotation under equal-and-opposite ECM reactions. No +x force,
   prescribed velocity, or `polarity_probability=0.65` is added. Fixed and
   moving conditions share the same network, crosslinks, and proposed stochastic
   events.

## Why this generation exists

The planning diagnosis was:

> G3A/B/C 只有 3–8 條、沒有 crosslink、而且每條 fiber 兩端固定；G2 雖然有
> 99-fiber outer-boundary network，卻把每個 2D intersection 都設為 crosslink。
> 所以下一版最合理的做法不是直接繼續調 G3，而是保留 G3，建立一個新的
> elastic calibration generation：沿用 G2 的網路幾何與 boundary anchoring、
> G3 的 material-point clutch/contact，再把 crosslink 改成有機率的永久
> elastic links。先固定 cell 把 ECM 調對，再恢復 protrusion 與 cell movement。

That paragraph records the gap found in the **superseded early G3 fixtures**.
The active G3 spheroid was subsequently rebuilt with 120–150 fibres,
outer-boundary anchoring and a 0.3 link fraction. G4 does not erase that history:
it turns the same scientific concern into an explicit calibration and validation
sequence with controlled coefficients and graph-distance readouts.

## G4A coefficients and what they mean

| Website control | Implemented quantity | What changing it tests |
|---|---|---|
| Pull intensity | total scalar force `F_pull` shared by eligible direct contacts | active loading amplitude; not a longer display arrow |
| Bending multiplier | `EI / EI₀`, applied only to the discrete curvature force | how easily a fibre rotates/straightens without changing axial `EA` |
| Bead drag | `ζ_b` in `ζ_b ṙ = ΣF` | response time and dissipation; it does not change elastic equilibrium |
| Crosslink probability | keep an intersection if a fixed random mark is `< p_x` | how much of the geometric intersection network transmits force |
| Crosslink stiffness | `k_x` in `F_x = k_x δ_x` | compliance of each retained link; it is not link number |

All slider cases use the same geometry. Crosslink probability uses nested sets:
increasing `p_x` adds links without reshuffling the links retained at lower
probability. The browser selects exact Python-precomputed cases; it does not run
a second or interpolated mechanics engine.

## G4C: when does a clutch slip?

There is deliberately no hard “too long” or “too heavy” threshold. For a bound
effective clutch,

\[
F=k_c x,\qquad k_\mathrm{off}(F)=k_\mathrm{off}^0e^{F/F_b},\qquad
P_\mathrm{off}=1-e^{-k_\mathrm{off}\Delta t}.
\]

The defaults (`k_c=2 nN/µm`, `k_off^0=0.018 s⁻¹`, `F_b=1.5 nN`) are inherited
from the G2 V3 effective-clutch mechanism calibration. They give intuitive
landmarks, **not deterministic thresholds**:

| clutch force | extension `F/k_c` | median lifetime at constant force |
|---:|---:|---:|
| 0 nN | 0 µm | 38.5 s |
| 1.5 nN = `F_b` | 0.75 µm | 14.2 s |
| 3.0 nN = `2F_b` | 1.50 µm | 5.2 s |

These are coarse-grained effective sites, not literal single integrin molecules.
Bell (1978) supports the exponential slip-bond law; Bangasser & Odde (2013)
supports coupling it to motor force–velocity; Prahl et al. (2020) supports
opposing motor–clutch modules. The exact G4 values remain parameters to sweep
against experimental pulling–slip–recoil videos.

## Core equations

The code comments identify each formula where it is evaluated. The full equation
set and evidence boundaries are also rendered in the web notebook.

\[
\zeta_b\dot{\mathbf r}_i=\mathbf F_i^s+\mathbf F_i^b+
\mathbf F_i^x+\mathbf F_i^{rep}+\mathbf F_i^{active}
\]

\[
\mathbf F_x=k_x[(\mathbf p_b-\mathbf p_a)-\boldsymbol\ell_x^0]
\]

\[
\gamma_c\dot{\mathbf r}_c=-\sum_i\mathbf F_i^{active},\qquad
\gamma_\theta\dot\theta=\sum_m(\mathbf a_m-\mathbf r_c)\times
\mathbf F_m^{reaction}.
\]

## Run and test

```bash
python3 -m unittest discover -s generations/g4_interactive_calibration/tests -v
python3 generations/g4_interactive_calibration/build_demo.py
```

The generated GitHub Pages lab is `docs/g4-lab.html`.
