# Version lineage and preservation policy

This file is the map for reading the repository without confusing a scientific
model version with a temporary Git working branch.

For a human-readable visual map with the limitation → next-question reasoning,
equations, parameter provenance, collapsible professor feedback and Git evidence
under every version, open the [research notebook](https://gl0008.github.io/motor-clutch-collagen-model/#evolution).

## The lineage

```text
G1 cumulative prototype path
V0 → V1 → V2 → V3 → V4

G2 shared corrected engine (parallel controls, not a V2→V3→V4 chain)
├── V2 transmission
├── V3 migration
└── V4 plasticity
       │
       └── G3 spheroid-guided radial remodelling

G4: read large revisions across; read experiment stages down

v1 short-time       ⇢ v2 long-time      → v3/v4 rebuild       ⇢ guidance revisions
A mechanics           A mechanics         v3 random ECM          v4E guidance test
↓                     ↓                   ↓                      ↓
B transmission        B transmission      v4A controls           v4F corrected OFAT
↓                     ↓                   ↓                      ↓
C clutch              C clutch            v4B force direction    v4F_new lineage audit
↓                     ↓                   ↓
D motion              D motion            v4C clutch cycling
                      └─→ G5 branch        ↓
                                           v4D motion

G5: read evidence blocks across; read their stages down

A–C core build      → 0-series audit     → D revisions          → parallel tests
A scaffold             0A scale             D v1 baseline          ├ D v3 leader
↓                      ↓                    ↓                       ├ E plasticity
B pull                 0B cell number       D v2 modes             └ ablation
↓                      ↓
C stiffness            0C adhesion
```

In this diagram, horizontal `→` means a major model change and horizontal `⇢`
means a revision of the same broad question. Vertical `↓` means the experimental
stages inside that release. Forks are parallel questions and must not be rendered
as a single linear upgrade chain.

Generation 1 is the conceptual mainline: it records how the question developed
from the SLS prototype to elastic fibres, crosslinks, cell motion and the first
plasticity hypothesis. Generation 2 branches from the frozen G1 archive and
repeats V2–V4 with corrected mechanics, boundary conditions, calibration and
visual encoding. Generation 3 preserves G2 as evidence and replaces its prescribed
directional imbalance with spatial clutches and unbiased protrusion feedback.
Generation 4 preserves G3 and formalizes the next calibration sequence: tune
elastic ECM, verify indirect graph transmission, add slippage, then release the
rigid cell. G4 v1 is the short-time implementation; G4 v2 branches from it and
adds two-hour observation, multiscale rendering, explicit individual-clutch
state, a shared-load alternative and fixed/moving/mobile-ECM controls.
G4 v3 then rebuilt the ECM from uniform random finite fibers. G4 v4 branches
from that preserved state and adds a literature ledger, separates a contractile
cavity benchmark from the biological cell, keeps the biological cell radius
constant, decomposes surface force, adds clutch loading/failure, and only then
releases single-cell translation.
G4 v4E then preserves A–D and tests contact-guided persistence as a four-arm
additive experiment: G4D control, matrix cue only, protrusion memory only, and
their combination. It adds neither global polarity nor extra force. Its
predeclared E3 cue-direction endpoint is not supported, while E2/E3 show higher
mean persistence; failed cue-rotation and domain gates keep the result
exploratory.
G4 v4F preserves v4E and repairs its shared baseline: two-way collagen steric
reaction, current-position material-point search and memoryless E0/E1
rebinding. E1 and E2 then each add exactly one factor; E3 is their declared
combination. Its 20-seed experiment finds that memory increases displacement
and persistence without selecting a reproducible direction. The matrix-cue
effect is not supported, and timestep-direction, 180-degree cue-rotation and
domain gates remain not passed; v4F is therefore preserved as an exploratory
mechanism result rather than accepted directional migration.
G4 v4F_new adds no mechanics and reruns no simulation. It audits the v4E → v4F
transition and records that v4F bundled requested steric/contact corrections
with contact-reach harmonisation, an RNG correction, a front-index correction,
matched initial contacts, an integrator-order correction and a different
representative seed. The parent-child result is therefore not a one-factor
comparison. F_new also establishes a fixed lineage seed, a separate matched
ensemble and mandatory microstages for future numerical corrections.
Generation 5 first runs G5-0A scale-transfer, G5-0B cell-number-only and G5-0C
adhesion controls with the G4 molecular clutch frozen. Its A–E sequence then
adds the organoid scaffold, collective contraction, strain stiffening,
released-cell invasion and crosslink plasticity. G5D v2 compares adhesion modes;
G5D v3 preserves the earlier high-leader-traction test. The newest matched
ablation changes one mechanism at a time around the molecular-clutch collective
baseline. None of these versions overwrites the earlier G5 results.

The repository's default Git branch remains `main` because it is the reviewed
synthesis and GitHub Pages source. It can contain every generation so the website
can compare them, but it never replaces their permanent checkpoints.

## Permanent branch rule

- Every scientific generation has its own `generation/gN` branch, even after its
  reviewed files are merged into `main`.
- `main` is the complete readable synthesis; it is not a sixth generation and is
  not the only surviving copy of a model.
- Important `agent/...` and `codex/...` branches are retained when they contain a
  meaningful experiment or documentation trail. Merging them does not imply
  deletion.
- A version may be placed under a generation on the website only when a preserved
  branch and an evidence commit support that placement. The experiment log must
  also name its parent, question, change, result, limitation and next question.
- A future generation receives a new directory and new `generation/gN` branch;
  it must not be folded into the previous generation merely to simplify Git.

## Git branches and immutable tags

| Git name | Purpose | Mutation policy |
|---|---|---|
| `main` | Published catalogue containing every preserved model | receives reviewed merges |
| `generation/g1` | Frozen G1 lineage at the pre-correction archive | do not modify |
| `generation/g2` | Corrected branch descended from G1 | may receive later G2-only work |
| `generation/g3` | Accepted spheroid-guidance evidence checkpoint | preserve after integration |
| `generation/g4` | Accepted G4 v2 long-time/multiscale checkpoint | preserve after integration |
| `generation/g5` | Current molecular-clutch, lineage-control and ablation checkpoint | preserve after integration |
| `codex/g3-emergent-guidance` | Historical G3 construction branch | preserved implementation history |
| `agent/g4-interactive-calibration` | G4A–D construction/review branch | merge after tests and site validation |
| `agent/g4-v2-multiscale` | G4 v2 long-time/multiscale construction branch | merge after tests and site validation |
| `revision/g4-v3-random-to-aligned` | Preserved random-to-aligned mechanics branch | do not overwrite |
| `revision/g4-v4-single-cell-force-alignment` | G4 v4 literature-grounded fixed-radius single-cell branch | current review branch |
| `revision/g4-v4e-contact-guided-single-cell` | G4 v4E one-sided matrix cue + contact-created front-memory branch | current review branch |
| `revision/g4-v4f-corrected-ofat-guidance` | G4 v4F corrected shared baseline + OFAT guidance branch | current review branch |
| `revision/g4-v4f-new-lineage-audit` | Audit-only record of every v4E → v4F change; no new physics | preserve after review |
| `agent/g5-organoid-plan` | G5-0A–0C lineage controls, molecular-clutch organoid A–E, mode comparison and ablation | preserve staged results; merge after review |
| `codex/model-evolution-notebook` | Human-readable website and experiment-record development | preserve documentation history after merge |
| `g1-v0` … `g1-v4` | Stable pointers to the archived G1 release | immutable tags |
| `g2-v2` … `g2-v4` | Stable pointers to the documented corrected release | immutable tags |

The five G1 folders were introduced together when the earlier work was
reorganized, so their stable tags point to the same frozen G1 repository
snapshot and identify different folders inside that snapshot. Likewise, the
three G2 tags identify different model folders in one validated corrected
snapshot. This is more honest than inventing separate historical commits that
never existed.

Working branches such as `agent/...` are implementation history, not scientific
model versions themselves. Important ones remain available after merge because
their commit sequence is useful evidence; the permanent `generation/gN` branch
is still the stable scientific checkpoint.

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
| G4 v1 A–D | [`generations/g4_interactive_calibration/`](generations/g4_interactive_calibration/) | [preserved G4 v1 lab](https://gl0008.github.io/motor-clutch-collagen-model/g4-lab.html) | preserved short-time generation |
| G4 v2 A | [`generations/g4_v2_multiscale/`](generations/g4_v2_multiscale/) | [mechanics calibration](https://gl0008.github.io/motor-clutch-collagen-model/g4-v2.html#a) | current long-time generation |
| G4 v2 B | [`generations/g4_v2_multiscale/`](generations/g4_v2_multiscale/) | [realignment paths](https://gl0008.github.io/motor-clutch-collagen-model/g4-v2.html#b) | current long-time generation |
| G4 v2 C | [`generations/g4_v2_multiscale/`](generations/g4_v2_multiscale/) | [clutch failure](https://gl0008.github.io/motor-clutch-collagen-model/g4-v2.html#c) | current long-time generation |
| G4 v2 D | [`generations/g4_v2_multiscale/`](generations/g4_v2_multiscale/) | [fixed/moving cell](https://gl0008.github.io/motor-clutch-collagen-model/g4-v2.html#d) | current long-time generation |
| G4 v3 | [`generations/g4_v3_random_to_aligned/`](generations/g4_v3_random_to_aligned/) | [random-to-aligned benchmark](https://gl0008.github.io/motor-clutch-collagen-model/g4-v3.html) | preserved parent of v4 |
| G4 v4 A–D | [`generations/g4_v4_single_cell_force_alignment/`](generations/g4_v4_single_cell_force_alignment/) | [single-cell force–alignment lab](https://gl0008.github.io/motor-clutch-collagen-model/g4-v4.html) | current review branch |
| G4 v4E E0–E3 | [`generations/g4_v4e_contact_guided_single_cell/`](generations/g4_v4e_contact_guided_single_cell/) | [contact-guided persistence lab](https://gl0008.github.io/motor-clutch-collagen-model/g4-v4e.html) | additive review branch |
| G4 v4F E0–E3 | [`generations/g4_v4f_corrected_ofat_guidance/`](generations/g4_v4f_corrected_ofat_guidance/) | [corrected OFAT guidance lab](https://gl0008.github.io/motor-clutch-collagen-model/g4-v4f.html) | corrected review branch |
| G4 v4F_new audit | [`generations/g4_v4f_new_lineage_audit/`](generations/g4_v4f_new_lineage_audit/) | [v4E → v4F lineage audit](https://gl0008.github.io/motor-clutch-collagen-model/g4-v4f-new.html) | audit-only branch |
| G5-0A / 0B / 0C | [`generations/g5_organoid/lineage.py`](generations/g5_organoid/lineage.py) | [scale-transfer, cell-number and adhesion controls](https://gl0008.github.io/motor-clutch-collagen-model/g5.html#lineage-controls) | latest lineage audit |
| G5 A–E | [`generations/g5_organoid/`](generations/g5_organoid/) | [organoid stages](https://gl0008.github.io/motor-clutch-collagen-model/g5.html) | implemented organoid sequence |
| G5D v2 / v3 | [`generations/g5_organoid/`](generations/g5_organoid/) | [current adhesion modes and preserved leader-traction test](https://gl0008.github.io/motor-clutch-collagen-model/g5.html#leader-update) | current mode comparison + earlier ablation |
| G5 ablation | [`generations/g5_organoid/`](generations/g5_organoid/) | [eight-way mechanism comparison](https://gl0008.github.io/motor-clutch-collagen-model/g5.html#clutch-baseline) | newest G5 result |

## Rule for adding a future model

1. Name the exact parent commit, create a new directory, and create a separate
   generation/revision branch; never replace an earlier version directory.
2. Put its purpose, changed assumptions, complete equations, evidence boundary,
   run instructions, tests and result interpretation in that directory's
   `README.md`.
3. Reuse shared code only through an explicitly named generation-level module.
4. Precompute visualization data from Python; the website only plays those
   frames.
5. Record the version's preserved branch and evidence commit, plus the human
   question → change → result → limitation → next-question chain.
6. Add the version to the website and to
   [`references/README.md`](references/README.md) only after that evidence exists.
   Declare whether its map edge adds an experimental block, revises the same
   question, or creates a parallel branch; never infer lineage from array or
   commit order.
7. Before running a parent-child comparison, declare an allowlist of changed
   components and snapshot the full configuration, geometry hash, initial-contact
   hash and RNG version. Use a fixed lineage seed separately from the ensemble.
8. Treat an RNG, integrator, initial-condition or index-bias correction as its
   own microstage; do not bundle it with a biological mechanism.
9. After tests pass, merge to `main`, retain the generation and important working
   branches, and create an immutable model tag when a release is declared.

These rules preserve every result, including negative results, and prevent a
later hypothesis from silently changing the meaning of an earlier model.
