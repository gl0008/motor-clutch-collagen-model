# Generation 2 · V4 — gated contact-induced plasticity

V4 is intentionally downstream of the corrected elastic-network and migration
gates.  It will not run unless `validation_summary.json` reports passing V2
network/mechanics checks and V3 mobility calibration.

## What changes from V2/V3?

Existing collagen bonds and permanent intersection crosslinks retain their
original rest states and never break.  During the loaded interval only, two
different fibres may acquire a **new weak crosslink** when they

1. are within 0.45 µm,
2. are locally aligned within 30°, and
3. lie within 30 µm of the cell surface, and
4. have approached one another by at least 0.00025 µm (0.25 nm) relative to
   their initial separation.

The last condition prevents the code from relabelling two fibres that were
already close before loading as newly remodelled contact.  New-link checks
begin only after the force ramp is complete.

The new link is stress-free at its formation geometry,

\[
\mathbf d_{new,0}=\mathbf p_b(t_{form})-\mathbf p_a(t_{form}),
\]

then contributes

\[
U_{new}=\tfrac12k_{new}\lVert
(\mathbf p_b-\mathbf p_a)-\mathbf d_{new,0}\rVert^2.
\]

This follows the constructive part of Ban et al.'s plastic-remodelling idea:
loading brings fibres into new contact and weak crosslinks form in that loaded
configuration.  It replaces the legacy V4 rule that simply reset the rest state
of the same broken/reformed partner.

## Required interpretation

The animation compares elastic-only and new-link runs under the same ramp,
hold, ramp-down and 570 s fully unloaded relaxation.  Permanent links are gold;
new weak links are purple.  Residual displacement or alignment is reported only
after the active force has returned to exactly zero.

The molecular identity and kinetics of these weak links remain a hypothesis.
This stage is not evidence that all collagen plasticity occurs by this route.
