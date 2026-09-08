# G4 v3 — random collagen to measured alignment

## Why this branch exists

G4 v1/v2 answered useful questions about ECM stiffness, drag, crosslink force
transmission, clutch failure and a released cell. However, their geometry
generator inserted twelve curved fibres around the cell at evenly spaced
angles. That made contact reliable, but it also meant the near-cell geometry
was partly radial before the simulation started.

G4 v3 is a separate branch from the frozen G4 v2 commit. It removes that
geometric assumption. Every source-fibre angle is sampled uniformly on
\([0,\pi)\), then a cell-sized circular void is cut from the network. Radial
alignment is measured from the result; it is not drawn into the initial state.

No G6 is created. No existing G4 file is overwritten.

## Question and cumulative construction

The scientific question is deliberately narrow:

> Can an initially near-isotropic elastic collagen network develop local radial
> alignment under cell-scale contraction or traction, and which added mechanism
> first produces a measurable change?

The development sequence is cumulative rather than ablation-based:

| Block | New factor | Question answered | Status |
|---|---|---|---|
| Benchmark | Prescribed circular contraction; all geometric intersections linked | Can the base stretching/bending mechanics reproduce the published bending-to-stretching recruitment mechanism? | Implemented first |
| v3-0 | Random finite bead-fibres around an inert circular void | Is the realised finite network actually near-isotropic? | Implemented |
| v3-A | Probabilistic permanent elastic crosslinks | Is there a connected path for indirect force transmission? | Implemented; connectivity/alignment dose measured |
| v3-B | Fixed rigid cell and steric exclusion | Does merely adding the cell create no false alignment? | Passed: no force, displacement or alignment |
| v3-C | Fixed material-point traction | Does force-driven deformation increase radial order? | Implemented; 12–48 nN response is weak and local |
| v3-D | Existing G4 Bell/shared-load clutch kinetics | Does load-dependent attachment lifetime produce pull–slip–recoil cycles? | Paused: v3-C mechanics gate not yet passed |
| v3-E | Dynamic contact relocation | Does following the deformed nearest collagen change sustained remodeling? | Hypothesis stage; not claimed as literature-validated |

Each stage must reuse the same geometry and random seed as the previous stage.
The later stage is not accepted because its animation looks more dramatic; it
must change a predeclared metric in the predicted direction.

## Benchmark equations

For every unconstrained bead \(i\),

\[
\zeta_b\dot{\mathbf r}_i =
\mathbf F_i^{\mathrm{stretch}}+
\mathbf F_i^{\mathrm{bend}}+
\mathbf F_i^{\mathrm{crosslink}}.
\]

There is no SLS element, Brownian motion, active force, clutch, plasticity or
cell translation in this first benchmark.

The cell-shaped inner boundary is prescribed as

\[
\mathbf r_i(t)=[1-\varepsilon(t)]\mathbf r_i(0),
\qquad
\varepsilon(t)=\varepsilon_{\max}
\min\left(1,{t\over t_{\mathrm{ramp}}}\right).
\]

The primary alignment metric is

\[
S_r=\left\langle 2(\hat{\mathbf t}\cdot\hat{\mathbf e}_r)^2-1\right\rangle.
\]

Thus \(S_r=1\) is radial, \(S_r=-1\) is circumferential and \(S_r\simeq0\)
is isotropic. The relevant response is \(\Delta S_r=S_r(t)-S_r(0)\), not the
final value alone.

The deformation-mode diagnostic is

\[
E_R={E_{\mathrm{stretch}}\over E_{\mathrm{bend}}}.
\]

\(E_R>1\) indicates a stretching-dominated network-level state. This is a
mechanical classification, not by itself proof of collagen remodeling.

## Acceptance gates

1. Initial whole-field nematic order must be below 0.10 for the demonstration
   seed and reported for every seed.
2. Initial radial order in every shell must be reported; the simulation must
   use \(\Delta S_r\) to remove finite-sample bias.
3. No bead may begin inside the cell-sized void.
4. Higher crosslink probability must add links without changing geometry or
   removing links already present at lower probability.
5. Prescribed contraction must produce finite elastic energy and a positive
   required inward boundary force.
6. A single circular-cell result is not generalized until a seed ensemble is
   run. Weak or absent alignment is a valid result.

## Current result and stopping gate

At 20% prescribed contraction with every intersection linked, five independent
geometries produced positive near-cell \(\Delta S_r\) (mean 0.101), while the
25--50 um shell remained approximately unchanged. Reducing the nested link
probability reduced both connectivity and alignment.

The inherited force-controlled contact model did **not** reproduce that
response: 12, 24 and 48 nN total pull at \(p_x=0.35\) produced near-cell
\(\Delta S_r=0.0078, 0.0102, 0.0141\), respectively, with negligible middle-
and far-field change. Therefore v3-D clutch kinetics is intentionally not yet
added. Intermittent attachment would make this weak deterministic baseline
harder, not easier, to interpret.

This stopping result is not a time-step artefact. At 24 nN, repeating the
480 s calculation with \(\Delta t=0.025\) s instead of 0.05 s changed the
near-cell \(\Delta S_r\) by only \(7.1\times10^{-6}\) relative. Extending the
same case from 480 to 960 s increased near-cell \(\Delta S_r\) from 0.0102 to
0.0139, while the middle and far shells remained approximately unchanged. The
960 s case was still changing slowly, so this is a convergence diagnostic—not
a claim about the two-hour steady state.

## Decision before the next cumulative block

The next step is **not** to add clutch failure or contact relocation. It is to
choose and document one source-matched way to close the gap between the
prescribed-contraction benchmark and the weak force-controlled result:

1. use measured cell/spheroid boundary displacement as a displacement-controlled
   input and compare its reaction force with the present traction range;
2. replace the circularly symmetric cell by a documented elongated or
   contractile-dipole geometry, because cell shape changes the alignment field;
3. calibrate fibre modulus, bending rigidity, drag and contact geometry to one
   selected collagen preparation and imaging experiment.

These are alternative calibration questions, not interchangeable tuning knobs.
Whichever is selected must introduce one declared assumption at a time and
retain the present random geometry and metrics for comparison.

## Run

```bash
python -m generations.g4_v3_random_to_aligned.build_demo
```

The command writes a small human-readable summary beside the model and stores
the larger animation arrays separately under `docs/g4-v3-data/`.

See [ASSUMPTIONS.md](ASSUMPTIONS.md) for the explicit assumption ledger and
[REFERENCES.md](REFERENCES.md) for the claim-to-source map, and
[DECISION_LOG.md](DECISION_LOG.md) for the pass/fail interpretation of every
cumulative block.
