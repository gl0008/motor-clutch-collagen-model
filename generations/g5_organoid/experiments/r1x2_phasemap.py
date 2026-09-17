"""G5 R1xR2 phase map: cc_adhesion x n_leaders -> directionality of invasion.

Quantifies WHEN low adhesion + localized leaders produce a DIRECTIONAL, elongated strand,
using the committed R1+R2 driver (consistency.run_r2_invasion) on the cued network -- no new
model code.  Metrics: leader-follower separation (um, along the cue axis), aspect_ratio
(elongation), detached_fraction.  Renders annotated heatmaps + saves the grid.

    python generations/g5_organoid/experiments/r1x2_phasemap.py

Personal testing, NOT confirmed findings (CLAUDE.md 7.5).  Imposed cue, NOT swirling.
"""

from __future__ import annotations

from pathlib import Path
import gc
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.consistency import run_r2_invasion, r1_config  # noqa: E402

CC = [1.0, 2.0, 4.0, 6.0, 8.0]          # rows (cell-cell adhesion, nN/um)
NL = [0, 1, 2, 3, 5]                     # cols (leader number, cue-front)
BASE = dict(n_fibers=300, organoid_radius=32.0, domain_size=280.0, boundary_width=6.0,
            generation_attempts=20, dt=0.05, duration=1500.0, sample_interval=150.0,
            radial_cue=True, cue_angle=0.0, leader_location="cue_front",
            leader_stall_factor=6.0, budget_mode="fixed_per_leader")


def run_grid(outdir="output"):
    Path(outdir).mkdir(exist_ok=True)
    sep = np.zeros((len(CC), len(NL)))
    asp = np.zeros((len(CC), len(NL)))
    det = np.zeros((len(CC), len(NL)))
    mxi = np.zeros((len(CC), len(NL)))
    fpair = 0.0
    print("=== R1xR2 phase map (cued, cue_front stall x6, 1500 s, seed 23) ===")
    print(" cc  N_L | lead-foll_sep  aspect  detached  max_inv")
    for i, cc in enumerate(CC):
        for j, nl in enumerate(NL):
            cfg = r1_config(cc_adhesion=cc, n_leaders=nl, **BASE)
            out = run_r2_invasion(cfg, seed=23)
            f = out["frames"][-1]
            sep[i, j] = out["leader_follower_separation"]
            asp[i, j] = out["aspect_ratio"]
            det[i, j] = f["detached_fraction"]
            mxi[i, j] = f["max_cell_disp"]
            fpair = max(fpair, f["force_pair_residual"])
            print(" %3.0f  %3d | %+12.3f  %5.2f  %7.3f  %6.2f" % (cc, nl, sep[i, j], asp[i, j], det[i, j], mxi[i, j]))
            del out; gc.collect()
    print("max force_pair_residual over the whole grid: %.1e" % fpair)
    np.savez(f"{outdir}/r1x2_phasemap.npz", cc=np.array(CC), nl=np.array(NL),
             leader_follower_sep=sep, aspect_ratio=asp, detached_fraction=det, max_inv=mxi)
    _render(sep, asp, det, f"{outdir}/r1x2_phasemap.png")
    return sep, asp, det


def _panel(ax, Z, title, cmap):
    im = ax.imshow(Z, origin="lower", aspect="auto", cmap=cmap)
    ax.set_xticks(range(len(NL))); ax.set_xticklabels(NL)
    ax.set_yticks(range(len(CC))); ax.set_yticklabels([("%g" % c) for c in CC])
    ax.set_xlabel("n_leaders (cue-front)"); ax.set_ylabel("cc_adhesion (nN/um)")
    ax.set_title(title, fontsize=10)
    for i in range(len(CC)):
        for j in range(len(NL)):
            ax.text(j, i, "%.2f" % Z[i, j], ha="center", va="center", fontsize=7,
                    color="k")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def _render(sep, asp, det, path):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    _panel(axes[0], sep, "leader-follower separation (um)  [directionality]", "viridis")
    _panel(axes[1], asp, "aspect_ratio  [elongation, 1=round]", "magma")
    _panel(axes[2], det, "detached_fraction  [escape]", "cividis")
    fig.suptitle("G5 R1xR2 phase map: low adhesion + localized leaders -> directional strand "
                 "(personal testing, NOT confirmed; imposed cue, not swirling)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


if __name__ == "__main__":
    run_grid()
