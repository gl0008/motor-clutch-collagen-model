"""G5-R0 diagnostic driver: lock a coherent (non-spaghetti) baseline + 2 h gif.

R0 (the mechanical-consistency gate) is implemented in
``generations/g5_organoid/consistency.py`` (force-consistent multicellular driver
ported from Gloria's G4D).  This script does NOT change any model code; it only
runs ``run_r0_invasion`` to (a) sweep a few ECM/clutch knobs and confirm force
propagates as COHERENT STRAIN rather than "spaghetti", and (b) render the headline
**2 h (7200 s)** diagnostic gif + numbers.

All output is personal testing, NOT confirmed findings (CLAUDE.md 7.5).  No swirling
is claimed.  New output filenames only (prefix ``r0_``); nothing is overwritten.

    python generations/g5_organoid/experiments/r0_diagnostic.py sweep      # fast, lock baseline
    python generations/g5_organoid/experiments/r0_diagnostic.py headline   # 2 h gif + npz + numbers
    python generations/g5_organoid/experiments/r0_diagnostic.py all
"""

from __future__ import annotations

from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.consistency import run_r0_invasion, r0_config  # noqa: E402
from generations.g5_organoid import visualize as viz  # noqa: E402


def _reach_um(frame, thresh=0.05):
    """Farthest shell (outer radius, um) whose mean bead displacement > thresh um."""
    reach = 0.0
    for sh in frame["shells"]:
        if sh["mean_displacement"] > thresh:
            reach = sh["shell"][1]
    return reach


def sweep(outdir="output"):
    """Short runs over collagen_modulus x crosslink_fraction x total_pull.

    Reports the R0 consistency diagnostics per cell so we can lock a window where
    force propagates coherently (spaghetti_index << 1) with the force-pair intact.
    """
    import gc
    base = dict(n_fibers=90, n_corona_fibers=32, organoid_radius=28.0, domain_size=260.0,
                boundary_width=5.0, generation_attempts=15, dt=0.05, duration=150.0,
                sample_interval=50.0)
    # one-knob-at-a-time about a centre (mod 3.0, xlink 0.85, pull 24) to keep it light
    combos = [(1.5, 0.85, 24.0), (3.0, 0.85, 24.0), (6.0, 0.85, 24.0),
              (3.0, 0.60, 24.0), (3.0, 0.85, 12.0)]
    grid = []
    for modulus, xlink, pull in combos:
        cfg = r0_config(collagen_modulus_mpa=modulus, crosslink_fraction=xlink,
                        total_pull_force=pull, **base)
        out = run_r0_invasion(cfg, seed=23)
        f = out["frames"][-1]
        grid.append((modulus, xlink, pull, out["max_force_pair_residual"],
                     f["spaghetti_index"], _reach_um(f), f["mean_cell_radial_disp"],
                     out["cumulative_site_failures"], out["n_relocations"]))
        del out
        gc.collect()
    print("\n=== R0 baseline sweep (shared clutch, cell_drag=600, 300 s, seed 23) ===")
    print(" mod  xlink  pull | fpair_resid  spaghetti  reach_um  mean_inv_um  fails  relocs")
    for r in grid:
        print(" %.1f  %.2f  %4.0f | %10.2e  %8.4f  %7.0f  %+10.3f  %5d  %5d" % r)
    ok = all(r[3] < 1e-9 and r[4] < 1.0 for r in grid)
    print("\nforce-pair residual < 1e-9 AND spaghetti_index < 1.0 in EVERY cell: %s" % ok)
    print("(spaghetti_index << 1 everywhere -> the R0 fixes give coherent strain, not spaghetti)")
    return grid


def headline(outdir="output"):
    """Full 2 h (7200 s) run at the locked baseline -> diagnostic gif + npz + numbers."""
    Path(outdir).mkdir(exist_ok=True)
    # locked R0 baseline: default organoid scale (~43 cells), softened collagen, shared
    # clutch, G4D drag/speed.  2 h simulated so cells actually translate a visible amount.
    # ISOTROPIC-RANDOM network (run_r0_invasion default network_mode="random"): ~19 cells,
    # ~400 fibres for grip/percolation at this scale (t=0 radial order ~0, not corona-biased),
    # so the radial reorganisation is a genuine OUTPUT.  New filenames (does NOT overwrite the
    # earlier corona-network r0_diagnostic_2h.*).
    cfg = r0_config(organoid_radius=40.0, n_fibers=400, domain_size=340.0, boundary_width=6.0,
                    generation_attempts=20, duration=7200.0, sample_interval=60.0,
                    total_pull_force=24.0)
    t0 = time.time()
    out = run_r0_invasion(cfg, seed=23, snapshots=True)
    frames = out["frames"]
    f = frames[-1]
    resid = out["max_force_pair_residual"]
    spa = [fr["spaghetti_index"] for fr in frames]
    print("run %.0fs | cells %d beads %d sites %d mode %s" % (
        time.time() - t0, out["n_cells"], out["n_beads"], out["n_clutch_sites"], out["clutch_mode"]))
    print("max_force_pair_residual = %.3e  (gate <1e-9)" % resid)
    print("spaghetti_index: start %.4f -> end %.4f  (max %.4f; <<1 = coherent)" % (spa[0], spa[-1], max(spa)))
    print("site_failures %d == n_relocations %d ; slips %d" % (
        out["cumulative_site_failures"], out["n_relocations"], out["cumulative_slips"]))
    print("F_clutch mean/max %.2f/%.2f nN | F_steric %.2f | F_cell_cell %.2f nN" % (
        f["F_clutch_mean"], f["F_clutch_max"], f["F_steric_mean"], f["F_cell_cell_mean"]))
    print("clutch parallel/perp/total %.1f/%.1f/%.1f nN" % (
        f["clutch_parallel"], f["clutch_perp"], f["clutch_total"]))
    print("mean invasion %+.3f um | max %.3f um | reach %.0f um" % (
        f["mean_cell_radial_disp"], f["max_cell_disp"], _reach_um(f)))

    edges = out["edges"]
    np.savez_compressed(f"{outdir}/r0_diagnostic_random_2h.npz",
                        bead_snapshots=out["bead_snapshots"], cell_snapshots=out["cell_snapshots"],
                        edges=edges, centers0=out["centers0"], centers_final=out["centers_final"],
                        spaghetti=np.asarray(spa),
                        max_force_pair_residual=resid)
    span = float(np.max(np.abs(out["bead_snapshots"][0]))) * 1.02
    gif = viz.animate_invasion(
        out["bead_snapshots"], out["cell_snapshots"], edges, f"{outdir}/r0_diagnostic_random_2h.gif",
        cell_radius=cfg.cell_radius, span=span, fps=8,
        title="G5-R0 force-consistent invasion (2 h, ISOTROPIC-random ECM): coherent strain, force-pair=0 (personal testing)")
    print("wrote", gif)
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("sweep", "all"):
        sweep()
    if mode in ("headline", "all"):
        headline()
