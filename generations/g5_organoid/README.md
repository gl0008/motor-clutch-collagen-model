# Generation 5 — tumor organoid invading a collagen network

G5 scales the validated single-cell motor-clutch + G2 bead-spring collagen engine up to a
**multicellular organoid**: N motor-clutch **disks** held together by a **simplified cell–cell
adhesion** potential, all coupled to one shared G2 fibre network. It asks the organoid-scale version
of Hongbo's minimal question — *does the collective pull of many cells reorganise the surrounding
collagen into a radial (TACS-3-like) aster?* — grounded in Kolade/Gloria's 4MOSC1 OSCC/PDAC organoid
movies (~200 µm organoid, radial collagen alignment, clearance zones).

Full architecture, scale, parameters, validation targets and ablation design:
[`docs/G5_organoid_plan.md`](../../docs/G5_organoid_plan.md).

## What is built (Stages A–B)

- **Stage A — multi-cell scaffold.** `hex_centers` packs N disks; `cell_cell_forces` is a
  piecewise-linear soft adhesive-disk potential (equilibrium at `cell_spacing`, repulsion below,
  short-range adhesion above). A hex packing at that pitch is force-balanced at rest and cohesive
  under displacement. `make_organoid` builds a boundary-anchored G2 network with a **fibre-free
  organoid gap** plus a **near-field corona** of grippable collagen hugging the organoid.
  `multi_cell_repulsion` keeps beads out of every disk.
- **Stage B — contractile organoid → radial alignment.** `organoid_active_forces` has every surface
  cell grip nearby fibres (G2 Gaussian kernel, all-around) and reel them inward. `run_organoid_pull`
  integrates the fixed (non-translating) organoid and records the **radial-alignment order** in
  distance shells. Gate: the near-field radial order rises under pull.

## Design choices (decisions locked 2026-09-02)

- **Cell model = N disks + cell–cell adhesion** (not one deformable blob, not vertex/CPM yet). Chosen
  so collective vs single-cell escape can emerge in Stage D from the adhesion ↔ traction balance.
- **Softened, lightly crosslinked collagen** (3 MPa modulus, 10 nN/µm links) — same prof-requested
  softening as G3/G4, so the pull can visibly reorganise the near field. Not the stiff 32 MPa G2
  default.
- **Foundation:** reuses `g2_corrected/common/model.py` unchanged (`Network`, `contact_patches`,
  bead-spring energies). Does **not** extend g4_v2's `run_clutch` monolith.

## Performance (done)

All three walls from plan §6 are addressed, so a full-scale organoid (43 cells, ~24k beads, 420
fibres, 2400 steps) runs in **~8 s** (was ~300 s):

- **Numba integrator** (`_advance_organoid_numba` / `OrganoidStepper`) — the g4_v2 kernel generalised
  so the repulsion loop sums over M cell centres. Bit-identical to the NumPy force law (0.0 diff).
- **O(E) spatial-grid crosslinker** (`build_crosslinks_grid`) replaces g2's O(F²) fibre-pair build
  (nf=520 builds in <1 s, conn 0.92).
- **Cached candidate contacts** (`cell_candidate_fibers`) — fixed cells have a stable set of grippable
  near fibres, so contact detection skips the O(M·F) sweep.

`python generations/g5_organoid/visualize.py output/<run>.npz` renders a before/after radial-order
figure.

## Not yet built (later stages)

Cell translation/rotation (**Stage D**, collective vs single-cell), strain-stiffening (**C**, to
propagate the pull further and sharpen the aster), matrix plasticity / persistent aster (**E**). The
current near-field alignment is real but modest (sparse corona → only surface cells grip); Stage C and
a denser corona are the levers to a dramatic aster.

## Status

Own simulation output is **personal testing, not confirmed findings** (CLAUDE.md §7.5). Quantitative
validation targets (alignment index, displacement scaling) are in the plan; most are breast-cancer
proxies pending PDAC/OSCC numbers from Kolade/Gloria's movies.

## Run

```
python -m pytest generations/g5_organoid/tests/test_g5.py      # Stage A + B gates
python -c "from generations.g5_organoid.model import run_organoid_pull, OrganoidConfig; \
           run_organoid_pull(OrganoidConfig())"                # full-scale run
```
