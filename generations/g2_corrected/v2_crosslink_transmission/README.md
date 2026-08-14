# Generation 2 · V2 — elastic crosslink transmission

## Question isolated

When a fixed MDA-MB-231-sized cell pulls only collagen in a local contact
sector, do permanent crosslinks transmit some of that deformation to fibres
which the cell never touches directly?

The animation is a synchronized comparison of exactly the same initial network
and exactly the same 5 nN pulling history:

- **left:** the crosslink force term is switched off;
- **right:** all intersection crosslinks are permanent, elastic and freely
  hinged.

No bond breaks, no plastic rest length changes, and no SLS element is present.

## Equations implemented

Every collagen bead satisfies

\[
\zeta\dot{\mathbf r}_i = \mathbf F_i^{stretch}
+\mathbf F_i^{bend}+\mathbf F_i^{xlink}
+\mathbf F_i^{repulsion}+\mathbf F_i^{active}.
\]

For an axial bond of rest length \(\ell_0\),

\[
k_t=EA/\ell_0, \qquad
k_c=0.1k_t,
\]

where the smaller compressive stiffness is the baseline microbuckling
assumption.  The discrete bending coefficient is \(EI/\ell_0^3\).  At a
crosslink, two interpolated material points are joined by

\[
U_x=\tfrac12 k_x\lVert
(\mathbf p_b-\mathbf p_a)-\mathbf d_{x,0}\rVert^2.
\]

Cell coupling is deliberately hybrid.  First, the nearest material point on
each fibre must lie in the 3 µm surface shell and ±30° protrusion sector.  Only
those contact candidates receive Gaussian weights

\[
w_i=\frac{\exp[-d_i^2/\sigma_c^2]}
{\sum_j\exp[-d_j^2/\sigma_c^2]},\qquad \sigma_c=1.5\ \mu m.
\]

The total force is therefore conserved; Gaussian weighting does not apply a
small force to every fibre in the box.

## What the visualization must show

- circles are the numerical collagen beads;
- bead-to-bead segments are axial bonds, coloured by tension/compression;
- gold diamonds are crosslinks; a connector appears if their two material
  points separate under load;
- orange halos and arrows are the current Gaussian-weighted direct contacts;
- black squares are the only fixed beads and occur only at the outer boundary;
- the faint initial network is a reference, not another collagen layer;
- the 5× option scales displacement **arrows only**.  Fibre geometry always
  remains at true 1× scale.

## Acceptance gates

The generated network is rejected unless every cell-contact fibre reaches the
outer boundary through the crosslink graph and at least 85% of all fibres are
in a boundary-connected component.  Tests also check force normalization,
frame-zero contacts, no cell penetration, and non-contact-fibre response.

The default 32 MPa wet collagen-fibre modulus and 20–80 µm fibre lengths follow
the Lee et al. bead-network scale; the 400 kPa crosslink modulus is used only as
a literature-scale reference because mapping it to one 2-D point-link spring is
not unique.  The effective 75 nN/µm link penalty is therefore still an explicit
coarse-graining assumption.
