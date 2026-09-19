"""G5 coupled EMT-leader variant: partial-EMT cells that ALSO pull harder (biology-audit B).

BIOLOGY_AUDIT.md finding B: a real partial-EMT cell does NOT change cell-cell adhesion alone
-- it also raises traction / protrusive matrix-pulling, and leaders are often the partial-EMT
cells (biorxiv 2025). This script compares, at matched emt_fraction / adhesion:
  * ADHESION-ONLY (R1): emt_traction_factor = 1  (partial-EMT lowers cell-cell adhesion only)
  * COUPLED:            emt_traction_factor > 1  (same cells ALSO get higher per-site motor stall)
Uses consistency.run_coupled_emt_invasion (reuses R1 adhesion + R2 per-site-stall machinery) on
the isotropic-random network; partial-EMT cells are boundary-enriched (at the invasive front).

    python .../r1_coupled_emt.py sweep      # emt_traction_factor sweep (1200 s)
    python .../r1_coupled_emt.py headline    # 2 h: adhesion-only vs coupled + gifs (EMT cells red)

Personal testing, NOT confirmed findings (CLAUDE.md 7.5). No swirling.
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

from generations.g5_organoid.consistency import run_coupled_emt_invasion, r1_config  # noqa: E402
from generations.g5_organoid.experiments.r2_leader import animate_leaders  # noqa: E402

# partial-EMT at the front (boundary-enriched), half the cells, strong EMT (adhesion x0.15).
BASE = dict(emt_fraction=0.5, emt_adhesion_factor=0.15, emt_assignment="boundary",
            cc_adhesion=6.0, radial_cue=False)
SWEEP = dict(n_fibers=300, organoid_radius=32.0, domain_size=280.0, boundary_width=6.0,
             generation_attempts=20, dt=0.05, duration=1200.0, sample_interval=120.0)


def _row(out):
    f = out["frames"][-1]
    return (f["mean_cell_radial_disp"], f["max_cell_disp"], f["detached_fraction"],
            f["radius_of_gyration"], f["lcc_fraction"], f["force_pair_residual"])


def sweep():
    print("\n=== Coupled EMT-leader: emt_traction_factor sweep (boundary EMT, f=0.5, a=0.15, 1200 s) ===")
    print(" traction_x | mean_inv max_inv  detached   Rg    LCC   fpair    (1 = R1 adhesion-only)")
    for tf in (1.0, 2.0, 4.0, 6.0):
        cfg = r1_config(emt_traction_factor=tf, **BASE, **SWEEP)
        out = run_coupled_emt_invasion(cfg, seed=23)
        mi, mx, det, rg, lcc, fp = _row(out)
        print(" %9.1f | %+7.3f  %6.3f  %7.3f  %5.2f  %4.2f  %.0e" % (tf, mi, mx, det, rg, lcc, fp))
        del out; gc.collect()
    print("Q: does coupling traction to EMT (same cells loosen AND pull) increase invasion/escape")
    print("   vs the adhesion-only R1 (traction_x=1)?")


def _one(tag, title, outdir, **over):
    base = dict(n_fibers=400, organoid_radius=40.0, domain_size=340.0, boundary_width=6.0,
                generation_attempts=20, dt=0.05, duration=7200.0, sample_interval=60.0)
    cfg = r1_config(**{**BASE, **base, **over})
    t0 = time.time()
    out = run_coupled_emt_invasion(cfg, seed=23, snapshots=True)
    mi, mx, det, rg, lcc, fp = _row(out)
    print("[%s] run %.0fs | traction_x %.1f | fpair %.0e | mean_inv %+.2f max %.2f | detached %.3f | Rg %.2f | LCC %.2f | pEMT %d/%d" % (
        tag, time.time() - t0, out["emt_traction_factor"], fp, mi, mx, det, rg, lcc,
        len(out["emt_cell_ids"]), out["n_cells"]))
    edges = out["edges"]
    np.savez_compressed(f"{outdir}/coupled_{tag}_2h.npz", cell_snapshots=out["cell_snapshots"],
                        edges=edges, centers0=out["centers0"], centers_final=out["centers_final"],
                        emt_cell_ids=np.asarray(out["emt_cell_ids"]))
    span = float(np.max(np.abs(out["bead_snapshots"][0]))) * 1.02
    gif = animate_leaders(out["bead_snapshots"], out["cell_snapshots"], edges, out["emt_cell_ids"],
                          f"{outdir}/coupled_{tag}_2h.gif", cell_radius=cfg.cell_radius, span=span, title=title)
    print("wrote", gif)
    del out; gc.collect()


def headline(outdir="output"):
    Path(outdir).mkdir(exist_ok=True)
    # SAME emt_fraction/adhesion; ONLY difference = whether EMT also raises traction. EMT cells red.
    _one("adhesion_only", "Coupled-EMT control: partial-EMT lowers ADHESION ONLY (R1) -- red=pEMT (personal testing)",
         outdir, emt_traction_factor=1.0)
    _one("coupled", "COUPLED EMT-leader: partial-EMT lowers adhesion AND raises traction x6 -- red=pEMT (personal testing)",
         outdir, emt_traction_factor=6.0)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "sweep"
    {"sweep": sweep, "headline": headline}.get(mode, sweep)()
