# G3 mechanism-demo figures

These files are smoke demonstrations of the implemented mechanisms. They are not the
preregistered 100-seed biological validation and must not be presented as realistic 3D tumor
migration.

## G3A — material-point clutches

- Fixture: one tangential 40 µm fibre, both endpoints fixed; cell fixed.
- Seed/duration: 4 / 15 s.
- FOI: 0.6046 → 0.6239.
- Maximum bead displacement: 0.395 µm.
- Maximum bound effective clutches: 173/200.
- Force/torque residuals: below 4.0×10⁻¹⁶ / 9.0×10⁻¹⁸ relative.
- A five-seed, 5 s smoke gate gave ΔFOI = 0 for no pull and +0.00641 at 1× motor capacity.
- Extending the prescribed pull to ~23 s reaches the cell-overlap guard; the model does not
  contain steric pore traversal.

Files: `g3a_material_point_clutches.gif`, `g3a_summary.png`.

## G3B — emergent protrusion turnover

- Fixture: eight fibres along a horizontal nematic director; cell fixed.
- Seed/duration: 23 / 120 s with the unaccelerated 120 s protrusion lifetime.
- Active sectors changed at ~59 s and ~70 s; no global `+x` term is present.
- Maximum bound effective clutches: 11/200.
- This single run demonstrates turnover and spatial attachment only. It does not establish the
  ensemble directionality/ablation gates.

Files: `g3b_emergent_protrusions.gif`, `g3b_summary.png`.

## G3C — reaction-driven translation and rotation

- Asymmetric fixture seed/duration: 3 / 30 s.
- Net displacement: 0.00187 µm.
- Final body rotation: 6.68×10⁻⁵ rad.
- The nonzero values demonstrate that full-vector attachments can generate reaction force and
  torque. Their tiny magnitude is not realistic migration and speed is not an independent
  prediction.
- `g3c_direction_controls.gif` compares 20 s isotropic, aligned, and 30°-rotated aligned
  fixtures; it is a visualization control, not the 100-seed covariance result.

Files: `g3c_translation_rotation.gif`, `g3c_summary.png`,
`g3c_direction_controls.gif`.

## Elastic load–unload diagnostic

A 15 s G3A pull followed by the full 600 s recovery window retained κ ≈ 0.489 and elastic
energy ≈1.26% of peak, so the result is correctly marked `unresolved_recovery`. This is not
evidence of plasticity: the active model contains no irreversible mechanism. The κ < 0.1 gate
has therefore not yet passed under the configured recovery window.

> Personal simulation / mechanism validation; not yet a realistic 3D tumor-migration prediction.
