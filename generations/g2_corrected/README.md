# Generation 2 — corrected collagen-model series

This folder is the corrected successor to the original `versions/` series.
The legacy folders and website pages are preserved; Generation 2 does not
overwrite them.

## Shared engine

[`common/model.py`](common/model.py) contains the only collagen mechanics used
by corrected V2–V4:

\[
\zeta\dot{\mathbf r}_i=
\mathbf F_i^{stretch}+\mathbf F_i^{bend}+\mathbf F_i^{xlink}
+\mathbf F_i^{repulsion}+\mathbf F_i^{active}.
\]

The source comments identify which lines implement each term.  There is no SLS
in this engine.  `common/serialization.py` only compresses precomputed bead
positions for the website and does not calculate physics.

## Stage map

1. [`v2_crosslink_transmission/`](v2_crosslink_transmission/) — fixed cell,
   same-network no-crosslink/crosslinked comparison.
2. [`v3_two_sided_migration/`](v3_two_sided_migration/) — two clutch ensembles,
   fixed versus one-axis released rigid cell.
3. [`v4_contact_plasticity/`](v4_contact_plasticity/) — gated load–unload test
   of newly approached, stress-free weak links.

## Gate order

Run `python3 generations/g2_corrected/validate.py`.  It saves
[`validation_summary.json`](validation_summary.json).  Corrected V4 refuses to
run unless `v2_network`, `v2_mechanics`, and `v3_mobility` are all true.

The current V4 result is intentionally negative: only one genuinely new contact
forms at the accepted load, and it produces no resolved excess residual over
the elastic control.
