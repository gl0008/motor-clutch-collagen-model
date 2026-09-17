"""G5-R1 experiments: adhesion window + discrete partial-EMT composition (adhesion only).

Built on the R0 force-consistent driver (``consistency.run_r1_invasion``): EMT changes
ONLY cell-cell adhesion (Ilina & Friedl 2020), never clutch/drag/motor budget, so any
change in invasion mode is attributable to adhesion alone.  Run each phase as its OWN
process (the machine is memory-tight):

    python generations/g5_organoid/experiments/r1_emt_adhesion.py adhesion     # k_adh window
    python generations/g5_organoid/experiments/r1_emt_adhesion.py composition  # f_pEMT sweep
    python generations/g5_organoid/experiments/r1_emt_adhesion.py factor        # a_pEMT sweep
    python generations/g5_organoid/experiments/r1_emt_adhesion.py headline      # 2 h gifs

Personal testing, NOT confirmed findings (CLAUDE.md 7.5).  No swirling.  New filenames.
"""

from __future__ import annotations

from pathlib import Path
import gc
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.consistency import run_r1_invasion, r1_config  # noqa: E402
from generations.g5_organoid import visualize as viz  # noqa: E402

# ISOTROPIC-RANDOM network (network_mode="random" is the run_r1_invasion default): needs
# ~300 fibres at this scale for grip + percolation (t=0 radial order ~0, not corona-biased).
SWEEP = dict(n_fibers=300, organoid_radius=32.0, domain_size=280.0, boundary_width=6.0,
             generation_attempts=20, dt=0.05, duration=900.0, sample_interval=100.0)


def _row(out):
    f = out["frames"][-1]
    rate = f["mean_cell_radial_disp"] / max(out["frames"][-1]["time"] / 3600.0, 1e-9)  # um/h
    return (f["radius_of_gyration"], f["detached_fraction"], f["lcc_fraction"],
            f["mean_cell_radial_disp"], f["max_cell_disp"], rate,
            f["force_pair_residual"], f["floppiness_index"])


def adhesion_sweep():
    print("\n=== R1 (a) homogeneous adhesion window (f_pEMT=0, 900 s, seed 23) ===")
    print(" k_adh | Rg     detached  LCC    mean_inv  max_inv  rate_um/h  fpair   flop")
    for cc in (1.0, 2.0, 4.0, 6.0, 8.0, 15.0, 30.0, 60.0):
        cfg = r1_config(cc_adhesion=cc, emt_fraction=0.0, **SWEEP)
        out = run_r1_invasion(cfg, seed=23)
        rg, det, lcc, mi, mx, rate, fp, flop = _row(out)
        print(" %5.1f | %5.2f  %6.3f   %5.3f  %+7.3f  %6.3f  %8.2f  %.0e  %.3f" % (
            cc, rg, det, lcc, mi, mx, rate, fp, flop))
        del out; gc.collect()
    print("cohesive-but-deformable = the k_adh where LCC=1/detached=0 yet the cluster still"
          " deforms (Rg grows, cells invade). Force-pair stays 0 throughout.")


def composition_sweep():
    cc = 6.0
    print("\n=== R1 (b) partial-EMT composition (cc_adhesion=%.0f, a_pEMT=0.25, 900 s) ===" % cc)
    print(" f_pEMT | Rg     detached  LCC    mean_inv  max_inv  rate_um/h")
    for f in (0.0, 0.1, 0.25, 0.5, 0.75, 1.0):
        cfg = r1_config(cc_adhesion=cc, emt_fraction=f, emt_adhesion_factor=0.25, **SWEEP)
        out = run_r1_invasion(cfg, seed=23)
        rg, det, lcc, mi, mx, rate, fp, flop = _row(out)
        print(" %6.2f | %5.2f  %6.3f   %5.3f  %+7.3f  %6.3f  %8.2f" % (f, rg, det, lcc, mi, mx, rate))
        del out; gc.collect()


def factor_sweep():
    cc, f = 6.0, 0.5
    print("\n=== R1 (c) partial-EMT adhesion factor a_pEMT (cc=%.0f, f_pEMT=%.2f, 900 s) ===" % (cc, f))
    print(" a_pEMT | Rg     detached  LCC    mean_inv  max_inv  (0.15 = legacy extreme)")
    for a in (0.75, 0.50, 0.25, 0.15):
        cfg = r1_config(cc_adhesion=cc, emt_fraction=f, emt_adhesion_factor=a, **SWEEP)
        out = run_r1_invasion(cfg, seed=23)
        rg, det, lcc, mi, mx, rate, fp, flop = _row(out)
        print(" %6.2f | %5.2f  %6.3f   %5.3f  %+7.3f  %6.3f" % (a, rg, det, lcc, mi, mx))
        del out; gc.collect()


def _headline_one(tag, title, outdir, **cfg_over):
    base = dict(n_fibers=400, organoid_radius=40.0, domain_size=340.0,
                boundary_width=6.0, generation_attempts=20, dt=0.05, duration=7200.0,
                sample_interval=60.0)
    cfg = r1_config(**{**base, **cfg_over})
    t0 = time.time()
    out = run_r1_invasion(cfg, seed=23, snapshots=True)
    f = out["frames"][-1]
    print("[%s] run %.0fs | cells %d | fpair %.0e flop %.3f | Rg %.2f detached %.3f LCC %.3f | mean_inv %+.2f max %.2f um | pEMT %d/%d" % (
        tag, time.time() - t0, out["n_cells"], f["force_pair_residual"], f["floppiness_index"],
        f["radius_of_gyration"], f["detached_fraction"], f["lcc_fraction"],
        f["mean_cell_radial_disp"], f["max_cell_disp"],
        int(round(out["emt_fraction"] * out["n_cells"])), out["n_cells"]))
    edges = out["edges"]
    np.savez_compressed(f"{outdir}/r1_{tag}_2h.npz",
                        cell_snapshots=out["cell_snapshots"], edges=edges,
                        centers0=out["centers0"], centers_final=out["centers_final"],
                        emt_phenotype=out["emt_phenotype"])
    span = float(np.max(np.abs(out["bead_snapshots"][0]))) * 1.02
    gif = viz.animate_invasion(out["bead_snapshots"], out["cell_snapshots"], edges,
                               f"{outdir}/r1_{tag}_2h.gif", cell_radius=cfg.cell_radius,
                               span=span, fps=8, title=title)
    print("wrote", gif)
    del out; gc.collect()


def headline(outdir="output"):
    Path(outdir).mkdir(exist_ok=True)
    # SAME base adhesion (cc=6) -> the ONLY difference is EMT composition, isolating its effect.
    _headline_one("cohesive", "G5-R1 cohesive (cc=6, f_pEMT=0) -- epithelial front (personal testing)",
                  outdir, cc_adhesion=6.0, emt_fraction=0.0)
    _headline_one("emt_escape", "G5-R1 partial-EMT (cc=6, f_pEMT=0.5, a_pEMT=0.15) -- escape (personal testing)",
                  outdir, cc_adhesion=6.0, emt_fraction=0.5, emt_adhesion_factor=0.15)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "adhesion"
    {"adhesion": adhesion_sweep, "composition": composition_sweep,
     "factor": factor_sweep, "headline": headline}.get(mode, adhesion_sweep)()
