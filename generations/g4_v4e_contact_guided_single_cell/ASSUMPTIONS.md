# G4 v4E assumption ledger

| ID | Assumption | Evidence | Expected impact | Required check |
|---|---|---|---|---|
| V4E-A001 | A locally aligned collagen tract can provide a migration direction. | **Direct mechanism support.** Riching et al. (2014) and Carey et al. (2016) report more persistent migration/protrusion stabilization along aligned collagen. | E1/E3 may gain positive cue-axis displacement. | Rotate the cue 90° and 180°; the response must rotate. |
| V4E-A002 | Protrusions near a locally aligned fiber are more likely to persist. | **Direct qualitative MDA-MB-231 support.** Carey et al. observed short probing structures stabilize and grow near radially aligned fibers. | Raises the probability of an aligned contact becoming the front. | Sweep $\beta_a=0,1,2,4$. |
| V4E-A003 | One probing contact lasting 120 s establishes a temporary front. | **Indirect support only.** Riching et al. observed fewer, longer protrusions in aligned collagen; no 120 s calibration was reported. | Converts momentary contact into directional persistence. | Repeat at 60 and 180 s. |
| V4E-A004 | Rear sites are not forcibly cut, but a naturally failed rear site cannot rebind behind the current front. | **Model assumption.** General migration cycles require rear release, but the exact rule is not calibrated by the two MDA-MB-231 collagen papers. | Prevents restored rear symmetry after natural Bell failure. | Report the rule explicitly; compare E2 with E0. |
| V4E-A005 | A 5 µm reach represents a probing protrusion. | **Scale support, not a force calibration.** Carey et al. describe initial probing structures shorter than 5 µm. | Enlarges the set of physically encountered material points. | Keep contact eligibility visible; never apply a Gaussian to the full ECM. |
| V4E-A006 | $\beta_a=\beta_m=2$, a 45° front cone, 30° cue half-width, and 15° tract SD are provisional. | **No direct quantitative calibration.** | Determines how sharply contacts concentrate around alignment/front memory. | Predeclared gain and width sweeps; do not tune by animation appearance. |
| V4E-A007 | One cell has one currently probing protrusion plus multiple adhesion sites. | **Coarse-grained model assumption.** It prevents every clutch site from being mislabeled as a simultaneous cell front. | Gives a unique maturation clock while retaining the G4D motor budget. | Website distinguishes protrusion from ordinary adhesion spokes. |
| V4E-A008 | Tract fibers are made by centroid-preserving rigid rotations of eligible random fibers; invalid rotations are skipped. | **Geometric experimental control, not biology.** | Keeps bead count, fiber count, centers and lengths matched while changing local orientation. | Crosslink density/connectivity within 5%; report tract count and length error. |

## Explicit exclusions

No organoid, leader cell, cell growth, deformable membrane, plasticity, proteolysis, crosslink breaking, Brownian motion, SLS element, or global polarization is introduced. The model is 2D and is not claimed as a quantitative 3D migration fit.
