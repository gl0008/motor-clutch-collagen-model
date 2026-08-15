# Version 4 — Crosslink turnover and permanent remodeling

## Research question

After the elastic V2/V3 baselines work, can stress-activated crosslink turnover
leave a residual fiber configuration after the cell stops pulling?

## Why load–unload is essential

Fiber displacement while force is still applied is not evidence of permanent
remodeling.  V4 therefore applies the V2 right-sector pull for 4 s, sets it to
zero for 4 s, and compares:

1. **Elastic control:** all crosslinks remain permanent.
2. **Plastic hypothesis:** half the links are eligible for force-accelerated
   dissociation and stochastic re-formation.

The cell is fixed in this comparison so that the only changed law is crosslink
turnover.  V3 remains the cleaner model for studying cell movement.

## Added crosslink kinetics

For an eligible bound link of force magnitude `F_x`,

\[
k_{off}^{x}=k_{off,0}^{x}\exp\left[
\frac{\max(0,F_x-F_x^*)}{F_x^{scale}}\right].
\]

An unbound link reforms with constant rate `k_on^x`.  At re-formation its rest
vector becomes the current separation of the two material neighborhoods.  This
rest-state reset erases stored crosslink strain, supplying a minimal plastic
mechanism that can preserve deformation after unloading.

## Important scientific boundary

This is a **mechanism-first hypothesis**, not a calibrated molecular collagen
crosslink law.  The present link reforms between the same two material
neighborhoods; it does not yet search for a new molecular partner, slide along
a fibril, merge fibers, or represent enzymatic degradation.  Those are separate
candidate mechanisms that require experimental evidence and parameters.

## Default assumptions

- Dynamic fraction: 0.50 of geometric crosslinks.
- Force threshold: 0.015 model force units.
- Force scale: 0.035 model force units.
- Unloaded off rate: 0.010 s⁻¹; re-form rate: 0.35 s⁻¹.
- All geometry, elastic mechanics, active loading and random seed match the
  permanent control.

The relevant endpoints are not the number of animated breaks alone.  Report
residual `ΔS_r`, residual RMS displacement, crosslink energy after unloading,
and the cumulative broken/re-formed counts.  A parameter sweep and experimental
calibration are required before biological conclusions.

## Default result: a useful negative control

The saved run produced 7 break events and 2 re-form events, but its end-of-run
alignment change (`ΔS_r = 0.008910`) is essentially the same as the permanent
control (`0.008921`).  RMS displacement is also indistinguishable.  Therefore
the current same-neighborhood rest reset is **not evidence of permanent
remodeling**.  It suggests that turnover alone is insufficient at this loading,
and that a next hypothesis should explicitly test new-partner crosslinks,
stress-driven sliding/merging, or fiber plastic lengthening.  This negative
result is retained rather than tuning parameters until the animation looks
different.

[Open the V4 load–unload lab](demo/index.html)

## Repository record

- **Lineage:** G1 V3 → **G1 V4**; its flaw motivates the separate G2 V4 rule.
- **Stable tag:** `g1-v4`; the negative result remains preserved.
- **Notebook:** [purpose, complete equations, evidence boundaries and HOWEVER](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g1-v4).
- **References used:** R2 Abhilash; R11 Ban; R12 Kim. See the
  [version-to-reference registry](../../references/README.md).
