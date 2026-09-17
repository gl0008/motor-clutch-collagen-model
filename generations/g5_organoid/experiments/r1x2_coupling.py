"""G5 R1xR2 coupling check: does LOW cell-cell adhesion + LOCALIZED leaders pull a strand
WITHOUT dynamic switching (R3)?  Uses the existing R1 (emt/adhesion) + R2 (localized
per-site-stall leaders) code via consistency.run_r2_invasion on the cued network.

    python .../r1x2_coupling.py sweep      # cc_adhesion x leaders(on/off), 1200 s
    python .../r1x2_coupling.py headline    # 2 h gifs: low-adhesion + leaders vs control

Personal testing, NOT confirmed findings (CLAUDE.md 7.5).  Imposed cue, NOT swirling.
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

from generations.g5_organoid.consistency import run_r2_invasion, r1_config  # noqa: E402
from generations.g5_organoid.experiments.r2_leader import animate_leaders  # noqa: E402

CUED = dict(n_fibers=300, organoid_radius=32.0, domain_size=280.0, boundary_width=6.0,
            generation_attempts=20, dt=0.05, radial_cue=True, cue_angle=0.0)
LEAD = dict(n_leaders=3, leader_location="cue_front", leader_stall_factor=6.0,
            budget_mode="fixed_per_leader")


def sweep():
    print("\n=== R1xR2: cc_adhesion x leaders (cued, cue_front stall x6, 1200 s, seed 23) ===")
    print(" cc  lead | strands max_len  lead-foll_sep  aspect  LCC   detached  mean_inv max_inv  fpair")
    for cc in (1.0, 2.0, 4.0, 6.0):
        for nl in (0, 3):
            over = dict(cc_adhesion=cc, duration=1200.0, sample_interval=120.0, **CUED)
            if nl:
                over.update(LEAD)
            else:
                over.update(n_leaders=0)
            out = run_r2_invasion(r1_config(**over), seed=23)
            f = out["frames"][-1]
            print(" %3.0f  %4d | %6d %6d  %+12.3f  %5.2f  %5.3f  %6.3f  %+7.3f  %6.3f  %.0e" % (
                cc, nl, out["strand_count"], out["strand_max_len"], out["leader_follower_separation"],
                out["aspect_ratio"], f["lcc_fraction"], f["detached_fraction"],
                f["mean_cell_radial_disp"], f["max_cell_disp"], f["force_pair_residual"]))
            del out; gc.collect()
    print("Q: does LOW cc + leaders make strands/elongation that cohesive cc + leaders (R2) did not?")


def _one(tag, title, outdir, **over):
    base = dict(n_fibers=400, organoid_radius=40.0, domain_size=340.0, boundary_width=6.0,
                generation_attempts=20, dt=0.05, duration=7200.0, sample_interval=60.0,
                radial_cue=True, cue_angle=0.0)
    cfg = r1_config(**{**base, **over})
    t0 = time.time()
    out = run_r2_invasion(cfg, seed=23, snapshots=True)
    f = out["frames"][-1]
    print("[%s] run %.0fs | leaders %s cc %.1f | fpair %.0e | strands %d max_len %d | sep %+.2f | aspect %.2f | mean_inv %+.2f max %.2f | LCC %.2f det %.2f" % (
        tag, time.time() - t0, out["leader_ids"], cfg.cc_adhesion, f["force_pair_residual"],
        out["strand_count"], out["strand_max_len"], out["leader_follower_separation"], out["aspect_ratio"],
        f["mean_cell_radial_disp"], f["max_cell_disp"], f["lcc_fraction"], f["detached_fraction"]))
    edges = out["edges"]
    np.savez_compressed(f"{outdir}/r1x2_{tag}_2h.npz", cell_snapshots=out["cell_snapshots"],
                        edges=edges, centers0=out["centers0"], centers_final=out["centers_final"],
                        leader_ids=np.asarray(out["leader_ids"]))
    span = float(np.max(np.abs(out["bead_snapshots"][0]))) * 1.02
    gif = animate_leaders(out["bead_snapshots"], out["cell_snapshots"], edges, out["leader_ids"],
                          f"{outdir}/r1x2_{tag}_2h.gif", cell_radius=cfg.cell_radius, span=span, title=title)
    print("wrote", gif)
    del out; gc.collect()


def headline(outdir="output"):
    Path(outdir).mkdir(exist_ok=True)
    _one("lowadh_leaders", "G5 R1xR2: low adhesion (cc=1.5) + 3 cue-front leaders (personal testing)",
         outdir, cc_adhesion=1.5, **LEAD)
    _one("lowadh_control", "G5 R1xR2: low adhesion (cc=1.5), NO leaders (control) (personal testing)",
         outdir, cc_adhesion=1.5, n_leaders=0)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "sweep"
    {"sweep": sweep, "headline": headline}.get(mode, sweep)()
