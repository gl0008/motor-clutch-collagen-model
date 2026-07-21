# Motor–clutch cell on an SLS bead–collagen network

This repository is a standalone research prototype built from three pieces:

1. the stochastic motor–clutch concept used by Adebowale et al.;
2. the discrete stretching/bending/crosslinked filament framework used by Zhao et al.;
3. collagen-specific bead-network and single-fibril viscoelastic measurements.

It is intentionally a minimal model. It is not the unpublished Adebowale CMS code and it is not a
validated reproduction yet. Its purpose is to make the proposed coupling explicit, executable, and
testable before calibration against the original experiments.

## Single-protrusion 1D experiment

The first causal experiment is implemented separately from the exploratory 2D model. A fixed cell
loads twelve stochastic clutches against the terminal bead of a 10-bead collagen fibre. The nine
serial axial bonds are identical SLS elements. Each bond uses `9 * K0_chain` and `9 * Kinf_chain`,
so the end-to-end fibre stiffness is unchanged if the same fibre is redrawn with a different number
of beads. Internal beads are massless, quasi-static points in this 1D stage.

The public entry points are:

```python
from collagen_model import (
    SingleProtrusionConfig,
    run_single_protrusion_pair,
    run_lifetime_ensemble,
    run_mechanism_sweep,
)

config = SingleProtrusionConfig()
pair = run_single_protrusion_pair(config, seed=17)
ensemble = run_lifetime_ensemble(config, trials=200, seed=1000)
sweep = run_mechanism_sweep(config)
```

Fast and slow conditions use a counter-addressed common random stream. A cluster episode starts at
the first binding event and ends only when the bound-clutch count returns to zero. An episode still
active at 30 loading times is retained as right-censored rather than incorrectly counted as a
completed lifetime.

Generate the paired trajectory, 200-trial ensemble tables, and playback data with:

```bash
PY=/Users/glorialiu/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3
$PY run_single_protrusion.py --trials 200 --workers 8
```

## Implemented model

- A two-dimensional diluted triangular lattice represents a crosslinked collagen network.
- Lattice vertices are freely hinged crosslinks.
- Collinear bead triplets resist bending through a discrete-curvature energy.
- Every surviving bead-to-bead collagen segment is a standard linear solid (SLS), with instantaneous
  stiffness `k0`, long-time stiffness `kinf`, and relaxation time `tau`.
- Network beads follow overdamped dynamics.
- A point-like cell has radial protrusion modules. Each module contains stochastic clutches that bind
  nearby collagen beads, are loaded by inward actin flow, and detach with a Bell-type force-dependent
  rate.
- Equal and opposite clutch forces act on the collagen bead and cell.

The axial SLS bond uses

```text
F = kinf * extension + q
dq/dt = (k0 - kinf) * extension_rate - q/tau
```

where `q` is the force carried by the Maxwell branch. The implementation uses an exponential update,
not forward Euler, for the internal SLS state.

## Run

The project only needs Python and NumPy. In the Codex desktop workspace the bundled Python is:

```bash
PY=/Users/glorialiu/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3
$PY run_experiments.py --output results
$PY -m unittest discover -s tests -v
```

The experiment command writes:

- `sls_relaxation.csv`: numerical and analytic single-bond stress relaxation;
- `network_relaxation.csv` and `.svg`: virtual step-strain relaxation of the complete bead network;
- `trajectories.csv`: cell positions and clutch observables for every replicate;
- `condition_summary.csv`: ensemble mean squared displacement (MSD), clutch occupancy, and force;
- `summary.svg`: dependency-free visual summary;
- `representative_network.svg` and node/edge CSV files: one coupled network before/after loading;
- `run_metadata.json`: exact parameters and limitations.

The single-protrusion command writes its outputs under `results/single_protrusion/`:

- `paired_timeseries.csv`: synchronized fast/slow trajectory observables;
- `ensemble_trials.csv`: paired loading-rate and completed-cluster counts;
- `cluster_episodes.csv`: observed and right-censored episode durations;
- `mechanism_sweep.csv`: the prespecified `De` and `Kchain/kc` diagnostic sweep;
- `ensemble_summary.json`: acceptance checks and censoring-aware median lifetimes;
- `run_metadata.json`: exact dimensionless parameters and random seeds.

For a faster smoke run:

```bash
$PY run_experiments.py --output results_smoke --duration 10 --replicates 2
```

## Interpretation boundary

The default parameters are a stable, literature-informed starting point, not a fitted biological
parameter set. In particular, a coarse clutch represents an adhesion unit rather than one integrin,
and the network spring constants have not yet been homogenized to the 2 kPa substrate used by
Adebowale et al. See [ASSUMPTIONS.md](ASSUMPTIONS.md) for the assumption register and calibration plan.

## Primary references

- K. Adebowale et al., *Enhanced substrate stress relaxation promotes filopodia-mediated cell
  migration*, Nature Materials 20, 1290–1299 (2021),
  https://doi.org/10.1038/s41563-021-00981-w
- H. Zhao et al., *Condensate-driven chromatin organization via elastocapillary interactions*, bioRxiv
  (2025), https://doi.org/10.1101/2025.06.12.659369
- B. Lee et al., *A Three-Dimensional Computational Model of Collagen Network Mechanics*, PLOS ONE
  9:e111896 (2014), https://doi.org/10.1371/journal.pone.0111896
- A. J. Licup et al., *Stress controls the mechanics of collagen networks*, PNAS 112, 9573–9578
  (2015), https://doi.org/10.1073/pnas.1504258112
- Z. L. Shen et al., *Viscoelastic Properties of Isolated Collagen Fibrils*, Biophysical Journal 100,
  3008–3015 (2011), https://doi.org/10.1016/j.bpj.2011.04.052
