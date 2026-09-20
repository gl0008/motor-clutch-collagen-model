# G4 v4F — Corrected OFAT contact-guidance experiment

G4 v4F preserves G4 v4E and corrects three omissions that belonged in every
condition rather than in an experimental treatment: cell reaction to collagen
steric force, contact search around the current cell position, and a genuinely
memoryless control.  It then repeats E0–E3 as a one-factor-at-a-time (OFAT)
comparison.

## Experimental matrix

| Condition | ECM geometry | Protrusion memory | Difference from its baseline |
|---|---|---:|---|
| E0 corrected control | Random | No | Shared corrected mechanics only |
| E1 matrix cue | One-sided radial tract | No | E1 − E0 isolates matrix geometry |
| E2 protrusion memory | Random | Yes | E2 − E0 isolates the cell response |
| E3 combined | One-sided radial tract | Yes | Combination; E3 − E1 and E3 − E2 isolate each factor conditionally |

E0 and E2 use the same random network. E1 and E3 use the same cue network.
All four conditions use the same initial contacts, clutch constants, cell drag,
motor budget, contact reach, counter-addressed random stream and observation
clock. Before the first front-memory event, E0/E2 and E1/E3 are numerically
identical. Seed, timestep, stochastic event channel and clutch index are mixed
with separate 64-bit constants; this prevents a seed change from acting like a
shift in clutch/site index.

## Shared corrected baseline

### Collagen and cell force balance

The collagen remains an overdamped elastic bead–spring network:

\[
\zeta_b\dot{\mathbf r}_i=
\mathbf F_i^{\mathrm{stretch}}+
\mathbf F_i^{\mathrm{bend}}+
\mathbf F_i^{\mathrm{crosslink}}+
\mathbf F_i^{\mathrm{clutch}}+
\mathbf F_i^{\mathrm{steric}}.
\]

The cell now receives both clutch and steric reactions:

\[
\zeta_c\dot{\mathbf r}_c=
-\sum_i\mathbf F_i^{\mathrm{clutch}}
-\sum_i\mathbf F_i^{\mathrm{steric}}.
\]

The second term is a mechanical consistency correction. It is not a new
biological guidance mechanism. Experimental work shows that collagen pore size
and steric hindrance materially constrain MDA-MB-231 invasion, but the precise
linear penalty used here remains a model assumption.

### Dynamic material-point contact

A bound clutch remains attached to the same collagen material point until the
site fails. A completely failed site searches the *current*, deformed network
around the *current* cell center. One closest continuous point per fiber is
eligible within the common finite reach. E0/E1 choose among these points using
Gaussian distance weight only:

\[
P_s^{0}=\frac{I_s\exp(-d_s^2/\sigma_c^2)}
{\sum_q I_q\exp(-d_q^2/\sigma_c^2)}.
\]

No previous tangent, attachment direction or world-coordinate direction enters
this baseline probability. The default failed-site search interval is 5 s,
which is shorter than the effective \(1/k_{on}\approx18\,s\) binding timescale.
Because no paper calibrates this search clock directly, it is predeclared as a
model assumption and checked at 0.5, 2 and 5 s.

## Added factors

E1 rotates a local one-sided subset of collagen fibers into a radial tract
without changing their contour lengths. This is an environmental cue, not a
cell polarity value.

E2 leaves the ECM random. A front vector appears only after a real bound contact
survives the maturation clock. Once that front exists, failed sites use:

\[
P_s^{\mathrm{memory}}\propto
I_sw_s\exp\!\left[
\beta_a(\hat{\mathbf t}_s\cdot\hat{\mathbf e}_s)^2+
\beta_m\hat{\mathbf p}\cdot\hat{\mathbf e}_s
\right].
\]

E3 contains both E1 and E2. It is intentionally the only non-OFAT condition.

## Interpretation rule

The control is not required to remain exactly at the origin in every random
realization. Its null expectation is no reproducible world-axis bias across
seeds. Claims are based on paired E1−E0, E2−E0, E3−E1 and E3−E2 effects with
confidence intervals, not on one animation.

The short first-front diagnostic is kept separately from the six-hour outcome.
It caught and corrected an address alias in the inherited counter stream:
``seed + 1`` had previously been equivalent to ``clutch item - 1``. With
independent tuple mixing, the 20-seed first-front directional resultant falls
from 0.675 to 0.221. This is a numerical baseline repair, not a fitted guidance
effect.

## Run

```bash
PYTHONPATH=. .venv/bin/python -m generations.g4_v4f_corrected_ofat_guidance.build_demo
PYTHONPATH=. .venv/bin/python -m generations.g4_v4f_corrected_ofat_guidance.validate_model
PYTHONPATH=. .venv/bin/python -m generations.g4_v4f_corrected_ofat_guidance.diagnose_front_bias
PYTHONPATH=. .venv/bin/pytest -q generations/g4_v4f_corrected_ofat_guidance/tests
```

## However

The cell is still a rigid 2D circle. The linear steric penalty is not a
calibrated 3D pore-crossing law, the 120 s front clock and guidance gains remain
provisional, and there is no cell deformation, nuclear mechanics, proteolysis,
plasticity or crosslink rupture. Negative E1–E3 results must be retained.

The completed 20-seed, six-hour experiment supports one limited conclusion:
E2−E0 increases net displacement by 6.20 micrometres (95% CI 4.12–8.27) and
persistence by 0.493 (0.379–0.607), but its direction remains seed dependent.
Neither E1−E0 nor E3 establishes positive cue-directed displacement.

The short numerical audit also leaves three predeclared gates **not passed**:
halving the timestep changes the displacement direction for the audited seed,
the 180-degree cue rotation does not reverse the trajectory, and the near-cell
result changes by more than the five-percent domain criterion. Front-event
ordering, pair-force balance, dynamic-contact reach, the memoryless E0 event
gate and pre-memory E0/E2 identity pass. Therefore v4F is an exploratory
mechanism result, not a validated directional-migration model.
