# G3 G2-scale equations: provenance and boundary

## What this suite is

G3A/B/C uses the existing G2 99-fibre, boundary-anchored, permanently crosslinked collagen network. A protrusion must first grow to, and bind, a material point; only then does a **controlled total 5 nN traction ramp** begin. This is a G2-mechanics bridge for visual/mechanistic comparison, **not** a claim that the 200-motor module itself produces 5 nN.

G3S is an explicitly labelled sensitivity model: eight coarse ECM-facing surface patches on a 35-um radius spheroid, no interior cells and no cell-cell mechanics. It uses 20 nN total (nominally 2.5 nN per patch). It is not a fitted spheroid-force prediction.

## Equation-to-source map

| Term implemented | Equation / algorithm | Origin | What is assumed here |
|---|---|---|---|
| Overdamped collagen mechanics | `zeta_i rdot_i = Fstretch + Fbend + Fcrosslink + Fclutch + Fcontact` | G2 implementation; bead-spring collagen architecture follows Lee et al., PLOS ONE 2014, doi:10.1371/journal.pone.0111896 | 2D, athermal, no viscoelastic SLS or remodelling.
| Axial fibre elasticity | `F = (EA/l0)(l-l0) e`; compression is softened | G2's frozen network model; collagen-network modelling convention in Lee et al. | Fibre diameter/modulus are G2 values, not image-fitted to this sample.
| Bending and permanent material-point crosslinks | discrete curvature energy plus `0.5 k_xl |p_b-p_a-rest|^2` | G2 model; material-point crosslink mechanics is a numerical collagen-network representation | Links are permanent; no weak/transient link rupture or plasticity.
| Point clutch to beads | Gaussian interpolation plus a first-moment correction | Numerical force-conservation discretisation introduced in this repository | It preserves net force/torque locally; it is not a molecular adhesion law.
| Motor-clutch binding, force-velocity and rupture | stochastic binding; Bell slip-bond `koff = k0 exp(F/Fb)`; linear force-velocity relation | Adebowale et al. 2021 SI Table 4 / SI Eq. 2, as already used in G3 config | 200 motors/clutches describe one coarse module.
| Explicit protrusion | `L(t+dt)=clip(L +/- v dt,0,Lmax)` and a tip can capture only inside the configured range | Coarse modelling choice; biology motivation from Carey et al., Nat Cell Biol 2016, doi:10.1038/ncb3407 | Not a measured protrusion-growth law for this cell line.
| Cell-intrinsic polarity (G3B) | normalized discrete activity update with local self/adhesion terms and noise, `sum_s a_s=1` | Inspired by mass-conserved wave-pinning: Mori, Jilkine & Edelstein-Keshet 2008, Biophys J, doi:10.1529/biophysj.107.120824 | A discrete coarse approximation, **not** that paper's calibrated PDE.
| Cell motion (G3C) | `rdot_c=Fcell/zeta_c`, `phidot_c=taucell/zeta_r`, with exact opposite clutch/contact reaction | Overdamped rigid-body force balance | Drag is a G2-scale numerical closure, not measured cytoplasmic drag.
| G2-scale load controller | after first attachment, `sum_bound |f_c| = Ftarget min((t-tattach)/tramp,1)` | Deliberate protocol bridge to G2's `total_pull_force=5 nN` | It overrides the magnitude applied to ECM while retaining attachment identity/kinetics. It must never be called a molecular motor prediction.
| Spheroid shell | 8 independent surface traction patches, no interior interaction | Conceptually motivated by spheroid-induced radial collagen alignment (Kopanska et al. 2016, doi:10.1016/j.mattod.2016.06.006) | 20 nN is a conservative sensitivity input. Collective measurements can be much larger and are system-dependent; no direct calibration is claimed.

## Why the deformation is still small at 5 nN

The frozen G2 baseline itself gives only nanometre-to-tens-of-nanometre displacement over short runs because a boundary-connected collagen network is stiff. In the regenerated movies, orange clutch force is physical and green displacement arrows are labelled **x25 visual aid**; the metric cards always show the unscaled displacement and angle. The 20 nN spheroid sensitivity is included precisely to test larger, collective loading without relabelling it as a single-cell result.

## Missing physics

Permanent plastic remodelling and weak crosslink unbinding are intentionally excluded. Nam et al. show why weak-link mechanics can create irreversible collagen plasticity (Biophys J 2016, doi:10.1016/j.bpj.2016.03.031; 2018, doi:10.1016/j.bpj.2017.10.048). Adding that would be a new model version, not a cosmetic GIF improvement.
