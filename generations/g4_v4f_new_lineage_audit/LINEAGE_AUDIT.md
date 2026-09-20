# G4 v4F_new lineage audit

## Audit verdict

| Question | Verdict |
|---|---|
| Does v4F descend from the preserved v4E implementation? | **Yes.** `d49858b` is the parent implementation checkpoint and precedes `921e8fc`. |
| Were the published scalar motor, clutch, cell-drag and guidance defaults preserved? | **Yes**, except for the newly introduced 5 s contact-search interval. |
| Did v4F change only steric reaction and current-position contact search? | **No.** Nine scientifically relevant changes are listed below. |
| Is v4E → v4F a valid one-factor comparison? | **No.** Mechanics, stochastic stream, initial contacts, contact reach, front selection and representative seed changed together. |
| Are E0–E3 comparisons inside v4F matched? | **Mostly yes**, under the common v4F baseline and its existing timestep, rotation and domain limitations. |

## Exact Git lineage

```text
d49858b  G4 v4E implementation
   ↓
c1248c2  G4 v4E evidence record
   ↓
921e8fc  G4 v4F implementation
   ↓
75d7dd9  G4 v4F evidence record
   ↓
f26751c  merge of v4F into main
```

The history is preserved correctly. The problem is not ancestry; it is that
the child commit contains more changed mechanisms than the short version label
made visible.

## Complete change register

| ID | Category | Change from v4E to v4F | Requested? | Changes dynamics? | Interpretation |
|---|---|---|---|---|---|
| FNEW-REQ-001 | Requested mechanics | Add equal-and-opposite cell reaction to the inherited collagen steric penalty. | Yes | Yes | Prevents the cell from receiving clutch reaction while ignoring contact resistance. Exact rigid-circle linear penalty remains provisional. |
| FNEW-REQ-002 | Requested mechanics | After complete site failure, search continuous material points around the current cell in the currently deformed ECM; retry unavailable sites every 5 s. | Yes | Yes | Removes frozen-world contact candidates. The 5 s retry clock is an added assumption. |
| FNEW-REQ-003 | Requested experimental design | Remove directional/tangent memory from E0/E1 relocation; reserve front memory for E2/E3. | Yes | Yes | Makes the within-v4F E0/E1 baseline conceptually memoryless. |
| FNEW-HARM-001 | Baseline harmonisation | Give E0/E1 the same 5 µm probing reach as E2/E3. v4E gave probing reach only to guided cases. | No—implementation-added | Yes | Changes the accessible contact set in the controls and is not a steric-only correction. |
| FNEW-NUM-001 | Numerical correction | Replace additive counter RNG addressing with independently mixed seed/step/channel/item addresses. | No—implementation-added | Yes, stochastic trajectories | Removes a seed–clutch-index alias, but all rupture/rebinding histories change. |
| FNEW-NUM-002 | Bias correction | Replace fixed site 0 as the only possible front with weighted competition among all mature bound sites. | No—implementation-added | Yes | Removes an index/angular-order bias but changes the front process itself. |
| FNEW-INIT-001 | Initial-condition correction | Replace guided weighted initial contacts with the same evenly distributed surface contacts in E0–E3. | No—implementation-added | Yes | Improves within-v4F matching, but parent and child no longer share the same initial contacts. |
| FNEW-NUM-003 | Integrator correction | Advance ECM and cell from the same \(t_n\) geometry, moving the cell only after the ECM update. | No—implementation-added | Yes | Removes operator-splitting imbalance but changes every subsequent state. |
| FNEW-PRES-001 | Presentation/statistics | Change representative animation seed 45 → 51 and selection from E3-typical to jointly typical E0–E3; add paired contrasts and telemetry. | No—implementation-added | Animation/interpretation | The displayed v4E and v4F trajectories use different networks and cannot be compared visually as a matched run. |

## Parameters demonstrably preserved

The published manifests agree on the following values:

| Parameter | v4E | v4F |
|---|---:|---:|
| timestep | 0.05 s | 0.05 s |
| cell radius | 10 µm | 10 µm |
| cell drag | 600 nN·s/µm | 600 nN·s/µm |
| clutches per site | 12 | 12 |
| clutch stiffness | 2 nN/µm | 2 nN/µm |
| on-rate | 0.055 s⁻¹ | 0.055 s⁻¹ |
| zero-force off-rate | 0.018 s⁻¹ | 0.018 s⁻¹ |
| Bell force | 1.5 nN | 1.5 nN |
| stall force per site | 8 nN | 8 nN |
| probing reach | 5 µm | 5 µm |
| front maturation/loss | 120/120 s | 120/120 s |
| guidance gains | 2/2 | 2/2 |
| newly defined search interval | — | 5 s |

Preserving these scalar values does not make the trajectories matched when the
RNG algorithm, initial contacts, reachable contacts, front rule and force
balance also change.

## Why the animations look dramatically different

The published representative runs are:

| Condition | v4E seed 45 net displacement | v4F seed 51 net displacement |
|---|---:|---:|
| E0 control | 74.762 µm | 0.405 µm |
| E1 matrix cue | 25.568 µm | 1.165 µm |
| E2 protrusion memory | 16.499 µm | 4.276 µm |
| E3 combined | 19.899 µm | 4.200 µm |

These numbers answer “what did the two published demonstrations show?” They do
**not** answer “what is the effect of steric reaction?” because:

1. the representative seeds and therefore displayed networks differ;
2. the stochastic stream differs;
3. initial attachment rules differ;
4. contact availability differs;
5. front selection differs;
6. steric reaction and integration order differ.

The mechanistic observation that v4F clutch and steric reactions nearly cancel
is useful for diagnosing v4F, but it is not a clean estimate of a steric effect.

## Allowed and disallowed claims

### Allowed

- v4F repairs several identifiable omissions in the v4E shared baseline.
- v4F uses two-way contact reaction and current-position material-point search.
- Within v4F, E1 adds matrix cue to E0, E2 adds memory to E0, and E3 combines
  them, subject to the documented numerical/domain limitations.
- The v4F ensemble supports a persistence/displacement effect of its memory
  module but not reproducible cue-directed migration.

### Not allowed

- “Steric reaction alone reduced E0 displacement from 74.8 to 0.4 µm.”
- “v4F is the same v4E run with one added term.”
- “The side-by-side v4E and v4F animations use the same seed and network.”
- “Every v4F correction was an explicitly requested biological change.”

## Contract for the next model revision

1. Name the exact parent commit and create a declared change allowlist.
2. Use seed 45 for the primary parent–child lineage view. Do not reselect that
   seed based on the child outcome.
3. Use seeds 41–60 for inference, separately from the lineage view.
4. Preserve geometry, initial contacts, RNG version and integration order unless
   the current microstage explicitly changes one of them.
5. Put RNG, integrator, bias and initial-condition corrections into independent
   microstages with their own result and commit.
6. Snapshot configuration, geometry hash, initial-contact hash and RNG version.
7. Fail validation if a non-allowlisted field or hash changes.

