# G3 equations and talking points

## One-sentence model

G3 asks whether a rigid cell can read a few elastic collagen fibres through spatial stochastic
clutches, stabilize successful protrusions, and move from clutch reaction forces without a
prescribed global polarity.

## Elastic ECM

\[
\zeta_i\dot{\mathbf r}_i
=\mathbf F_i^{\rm stretch}+\mathbf F_i^{\rm bend}+\mathbf F_i^{\rm clutch}.
\]

Stretching and bending are purely elastic. There is no SLS or plasticity in G3.

## Conservative cell--collagen contact

For a zero-radius collagen bead at distance

\[
d_i=\|\mathbf r_i-\mathbf r_c\|,
\qquad
\delta_i=\max(0,R-d_i),
\]

the one-sided penalty potential and bead force are

\[
U_i^{\rm contact}=\frac12 k_{\rm contact}\delta_i^2,
\qquad
\mathbf F_i^{\rm contact}
=k_{\rm contact}\delta_i
\frac{\mathbf r_i-\mathbf r_c}{d_i}.
\]

The cell receives the exact opposite reaction. The force is central, so ideal circular contact
adds no torque about the cell center. This is the repulsive signed-distance branch of the
linear elastic contact approach in Runser, Vetter & Iber 2024, Methods Eq. 2, specialized to
a rigid 2D circle and point-like ECM beads. G3 includes no contact adhesion or friction.

## Persistent material-point attachment

For clutch \(j\) bound to segment \((a,b)\),

\[
\mathbf x_j=(1-\alpha_j)\mathbf r_a+\alpha_j\mathbf r_b,
\qquad 0\leq\alpha_j\leq1.
\]

`fiber_id`, `segment_id`, and \(\alpha_j\) remain unchanged until unbinding. This is the key
difference from repeatedly selecting the nearest bead.

## Motor-clutch kinetics

\[
r_{{\rm off},j}=r_{\rm off}^0
\exp\!\left(\frac{|\mathbf F_j|}{F_b}\right),
\qquad
v_{r,k}=v_u\max\!\left(0,1-\frac{T_k}{N_{m,k}f_m}\right).
\]

Source: Adebowale et al. 2021 SI motor-clutch equations and Table 4.

## Local Gaussian projection

\[
w_{ij}\propto\exp\!\left(-\frac{d_{ij}^2}{2\sigma^2}\right),
\qquad \sum_i w_{ij}=1.
\]

Only beads on the attached fibre within \(3\sigma\) receive the point force. A zero-net-force
correction pair preserves the point force's first moment. The Gaussian is a numerical
representation, not a crosslinker energy.

## Emergent protrusion feedback

\[
G_k=\frac{C_k}{\max_l C_l+\epsilon}A_k,
\qquad
h_k=\tau_p^{-1}\exp(-\beta_GG_k-\beta_QQ_k).
\]

\(C_k\) is local collagen availability, \(A_k\) is nematic alignment, and \(Q_k\) is smoothed
traction success. Carey et al. 2016 supports the feedback direction but does not provide the
chosen coarse-grained rate constants.

## Rigid-body motion

\[
\gamma_c\dot{\mathbf r}_c=-\sum_j\mathbf F_j,
\qquad
\gamma_\theta\dot\theta_c
=\sum_j(\mathbf y_j-\mathbf r_c)\times(-\mathbf F_j).
\]

A clutch engages along the surface normal, then follows its evolving full vector. If every
force were constrained to remain radial on a circular cell, the torque would be identically
zero. No independent self-propulsion velocity is present.

With contact active, the translational force also includes
\(-\sum_i\mathbf F_i^{\rm contact}\). This preserves action--reaction balance and prevents
the rigid cell and collagen nodes from occupying the same space.

## FOI and κ

\[
{\rm FOI}_{\rm random}=2/\pi,
\qquad
\kappa=\frac{{\rm RMS}({\rm FOI}_{\rm post}-{\rm FOI}_{0})}
{{\rm RMS}({\rm FOI}_{\rm pre-unload}-{\rm FOI}_{0})}.
\]

Source: Nam et al. 2016 Eq. 7. For a finite synthetic fixture, \({\rm FOI}_0\) is its measured
initial value. G3 expects reversible alignment and κ near zero only when the loading signal is
large enough to resolve.
