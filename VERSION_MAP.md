# Version lineage and preservation policy

This file is the map for reading the repository without confusing a scientific
model version with a temporary Git working branch.

## The lineage

```text
G1 V0  →  G1 V1  →  G1 V2  →  G1 V3  →  G1 V4
                                      │
                                      └── G2 corrected baseline
                                           G2 V2 → G2 V3 → G2 V4
                                                       │
                                                       └── G3 emergent guidance
                                                            G3 spheroid remodelling
                                                                 │
                                                                 └── G4 calibration
                                                                      G4A → G4B → G4C → G4D
```

Generation 1 is the conceptual mainline: it records how the question developed
from the SLS prototype to elastic fibres, crosslinks, cell motion and the first
plasticity hypothesis. Generation 2 branches from the frozen G1 archive and
repeats V2–V4 with corrected mechanics, boundary conditions, calibration and
visual encoding. Generation 3 preserves G2 as evidence and replaces its prescribed
directional imbalance with spatial clutches and unbiased protrusion feedback.
Generation 4 preserves G3 and formalizes the next calibration sequence: tune
elastic ECM, verify indirect graph transmission, add slippage, then release the
rigid cell.

The repository's default Git branch remains `main` because it is the published
catalogue and GitHub Pages source. It contains both generations so that the
website can compare them. This does **not** mean that G2 overwrote G1.

## Git branches and immutable tags

| Git name | Purpose | Mutation policy |
|---|---|---|
| `main` | Published catalogue containing every preserved model | receives reviewed merges |
| `generation/g1` | Frozen G1 lineage at the pre-correction archive | do not modify |
| `generation/g2` | Corrected branch descended from G1 | may receive later G2-only work |
| `codex/g3-emergent-guidance` | Historical G3 construction branch | preserved implementation history |
| `agent/g4-interactive-calibration` | G4A–D construction/review branch | merge after tests and site validation |
| `g1-v0` … `g1-v4` | Stable pointers to the archived G1 release | immutable tags |
| `g2-v2` … `g2-v4` | Stable pointers to the documented corrected release | immutable tags |

The five G1 folders were introduced together when the earlier work was
reorganized, so their stable tags point to the same frozen G1 repository
snapshot and identify different folders inside that snapshot. Likewise, the
three G2 tags identify different model folders in one validated corrected
snapshot. This is more honest than inventing separate historical commits that
never existed.

Feature branches such as `agent/...` are temporary construction history. They
are not scientific model versions.

## Model directory map

| Model | Permanent directory | Website notebook | Stable tag |
|---|---|---|---|
| G1 V0 | [`versions/v0_sls_prototypes/`](versions/v0_sls_prototypes/) | [`?model=g1-v0`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g1-v0) | `g1-v0` |
| G1 V1 | [`versions/v1_few_fiber/`](versions/v1_few_fiber/) | [`?model=g1-v1`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g1-v1) | `g1-v1` |
| G1 V2 | [`versions/v2_crosslinked_elastic/`](versions/v2_crosslinked_elastic/) | [`?model=g1-v2`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g1-v2) | `g1-v2` |
| G1 V3 | [`versions/v3_two_sided_migration/`](versions/v3_two_sided_migration/) | [`?model=g1-v3`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g1-v3) | `g1-v3` |
| G1 V4 | [`versions/v4_plastic_remodeling/`](versions/v4_plastic_remodeling/) | [`?model=g1-v4`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g1-v4) | `g1-v4` |
| G2 V2 | [`generations/g2_corrected/v2_crosslink_transmission/`](generations/g2_corrected/v2_crosslink_transmission/) | [`?model=g2-v2`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g2-v2) | `g2-v2` |
| G2 V3 | [`generations/g2_corrected/v3_two_sided_migration/`](generations/g2_corrected/v3_two_sided_migration/) | [`?model=g2-v3`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g2-v3) | `g2-v3` |
| G2 V4 | [`generations/g2_corrected/v4_contact_plasticity/`](generations/g2_corrected/v4_contact_plasticity/) | [`?model=g2-v4`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g2-v4) | `g2-v4` |
| G3 active spheroid | [`generations/g3_spheroid_guidance/`](generations/g3_spheroid_guidance/) | [`?model=g3`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g3) | preserved current G3 |
| G3 early fixtures | [`legacy/g3_v1_superseded/`](legacy/g3_v1_superseded/) | historical notebook data retained | superseded archive |
| G4A | [`generations/g4_interactive_calibration/`](generations/g4_interactive_calibration/) | [`?model=g4-a`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g4-a) | pending reviewed G4 tag |
| G4B | [`generations/g4_interactive_calibration/`](generations/g4_interactive_calibration/) | [`?model=g4-b`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g4-b) | pending reviewed G4 tag |
| G4C | [`generations/g4_interactive_calibration/`](generations/g4_interactive_calibration/) | [`?model=g4-c`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g4-c) | pending reviewed G4 tag |
| G4D | [`generations/g4_interactive_calibration/`](generations/g4_interactive_calibration/) | [`?model=g4-d`](https://gl0008.github.io/motor-clutch-collagen-model/model-notebook.html?model=g4-d) | pending reviewed G4 tag |

## Rule for adding a future model

1. Create a new directory; never replace an earlier version directory.
2. Put its purpose, changed assumptions, complete equations, evidence boundary,
   run instructions, tests and result interpretation in that directory's
   `README.md`.
3. Reuse shared code only through an explicitly named generation-level module.
4. Precompute visualization data from Python; the website only plays those
   frames.
5. Add the version to the website and to
   [`references/README.md`](references/README.md).
6. After tests pass, merge to `main` and create an immutable model tag.

These rules preserve every result, including negative results, and prevent a
later hypothesis from silently changing the meaning of an earlier model.
