# G4 v4 — literature-grounded single-cell force–collagen alignment

## Why this generation exists

G4 v3 established an initially random, finite bead–fiber network and showed a
critical gap: a prescribed contracting inner boundary produced measurable local
radial alignment, whereas inherited 12–48 nN fixed traction was weak and local.
Its website drew the contracting boundary as a shrinking cell, which made a
mechanics benchmark look like a biological prediction.

G4 v4 is a new branch. G4 v3 is not overwritten. This generation separates:

1. a **contractile-cavity benchmark**, retained only to test ECM mechanics;
2. a **fixed-size biological cell**, where traction comes from surface contacts;
3. a **motor–clutch cell**, where actin loading continues after physical contact;
4. a **released cell**, where motion follows reaction-force imbalance.

The primary comparison is the MDA-MB-231 single-cell timescale used by Riching
et al.: collagen displacement over two hours and migration over six hours.
These times are observation clocks, not a claim that all molecular processes
share one timestep.

## Cumulative blocks

| Block | Added mechanism | Allowed interpretation |
|---|---|---|
| v4-A passive | Fixed 10 µm-radius cell, steric exclusion, no force | Negative control for geometry-induced alignment |
| v4-A cavity | Prescribed contraction of an inner cavity | ECM solver benchmark; **not cell shrinkage** |
| v4-B normal | Fixed cell, contact band, Gaussian distance weights | Deterministic force transmission and alignment |
| v4-B dipole | Two opposite traction sectors | Literature comparison for contraction anisotropy |
| v4-C independent control | Independent Bell rupture with clutch redundancy | Shows why individual ruptures need not detach a site |
| v4-C shared-load primary | Equal sharing of site load, Bell rupture, rebinding, relocation after complete site failure | Makes collective pull–slip–recoil testable without changing cell radius |
| v4-D | Release cell translation | Mechanical consequence of attachment imbalance; no prescribed global polarity |

Every later block reuses the same random source geometry and crosslinks unless
its result explicitly states otherwise.

## Key equations

The collagen bead dynamics are

\[
\zeta_b\dot{\mathbf r}_i=
\mathbf F_i^{\mathrm{stretch}}+
\mathbf F_i^{\mathrm{bend}}+
\mathbf F_i^{\mathrm{crosslink}}+
\mathbf F_i^{\mathrm{contact}}+
\mathbf F_i^{\mathrm{steric}}.
\]

Deterministic contact eligibility and Gaussian weighting are

\[
0\le d_s\le d_c,\qquad
w_s=\frac{\exp(-d_s^2/\sigma_c^2)}
{\sum_q\exp(-d_q^2/\sigma_c^2)}.
\]

For a local fiber tangent \(\mathbf t_s\), the same contact force is resolved as

\[
\mathbf F_{\parallel,s}=(\mathbf F_s\!\cdot\!\mathbf t_s)\mathbf t_s,
\qquad
\mathbf F_{\perp,s}=\mathbf F_s-\mathbf F_{\parallel,s}.
\]

In the clutch block, continued loading comes from relative actin/ECM motion:

\[
\dot x_s=v_0\max\!\left(0,1-\frac{F_s}{F_{\mathrm{stall}}}\right)
-(\mathbf v_{\mathrm{ECM},s}-\mathbf v_c)\!\cdot\!\hat{\boldsymbol\ell}_s,
\qquad F_c=k_cx_c,
\]

with slip-bond rupture

\[
k_{\mathrm{off}}(F_c)=k_{\mathrm{off}}^0
\exp\!\left(\frac{|F_c|}{F_b}\right).
\]

The shared-load primary case uses

\[
f_c=\frac{F_{\mathrm{site}}}{i},\qquad
r_i=i k_{\mathrm{off}}^0
\exp\!\left(\frac{F_{\mathrm{site}}}{iF_b}\right),
\]

where \(i\) is the number of bound clutches. The independent-clutch case is
retained as a control because its redundant bonds can hide complete site
failure even when many individual ruptures occur.

The released circular cell translates according to

\[
\zeta_c\dot{\mathbf r}_c=-\sum_s\mathbf F_s^{\mathrm{cell\to ECM}}.
\]

## Important scientific boundary

The implementation exposes the local collagen tangent at every contact and
measures radial versus boundary-parallel orientation. It does not claim that
the parallel/perpendicular decomposition alone proves migration. The fixed
cell first tests the proposed mechanical sequence:

\[
F_\perp\;\longrightarrow\;\text{bending/rotation},\qquad
F_\parallel\;\longrightarrow\;\text{axial tension/longer-range transmission}.
\]

Those are predeclared hypotheses. Failure to observe them is a valid result.

## Run

```bash
python -m generations.g4_v4_single_cell_force_alignment.build_demo --quick
python -m generations.g4_v4_single_cell_force_alignment.build_demo
```

The first command checks the pipeline with explicitly shortened clocks. The
second generates the two-hour fixed-cell and six-hour moving-cell datasets.
Generated files always record their actual duration and carry a `quick` flag.

See [LITERATURE_REVIEW.md](LITERATURE_REVIEW.md),
[PAPER_LEDGER.md](PAPER_LEDGER.md),
[ASSUMPTIONS.md](ASSUMPTIONS.md), [REFERENCES.md](REFERENCES.md), and
[DECISION_LOG.md](DECISION_LOG.md).
