"""G5-R2 experiments: static leader NUMBER + LOCATION (functional role, not EMT).

Built on the R0/R1 force-consistent driver (`consistency.run_r2_invasion`): leaders are a
LOCALIZED front on an explicitly-imposed one-sided radial cue (v4E `make_radial_tract_spec`
port; NOT swirling), and their higher traction EMERGES from per-site motor stall (Chan &
Odde 2008), NOT a post-hoc force multiply.  EMT (R1 adhesion) is a separate knob.  Run each
phase as its own process (memory-tight machine):

    python .../r2_leader.py count      # N_L sweep (matched-total budget)
    python .../r2_leader.py budget     # matched_total vs fixed_per_leader
    python .../r2_leader.py location   # cue_front vs perimeter
    python .../r2_leader.py headline   # 2 h gifs (leader front vs no-leader control)

Personal testing, NOT confirmed findings (CLAUDE.md 7.5).  No swirling.  New filenames.
"""

from __future__ import annotations

from pathlib import Path
import gc
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.animation import FuncAnimation, PillowWriter

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.consistency import run_r2_invasion, r1_config  # noqa: E402

# cued isotropic-random network; ~300 fibres for grip/percolation; one-sided cue toward +x.
SWEEP = dict(n_fibers=300, organoid_radius=32.0, domain_size=280.0, boundary_width=6.0,
             generation_attempts=20, dt=0.05, duration=900.0, sample_interval=100.0,
             radial_cue=True, cue_angle=0.0, cc_adhesion=6.0)


def _row(out):
    f = out["frames"][-1]
    return (out["strand_count"], out["strand_max_len"], out["leader_follower_separation"],
            out["aspect_ratio"], f["lcc_fraction"], f["detached_fraction"],
            f["mean_cell_radial_disp"], f["max_cell_disp"], f["force_pair_residual"])


def count_sweep():
    print("\n=== R2 (a) leader-count sweep (cued, cue_front, MATCHED-total budget, factor 4, 900 s) ===")
    print(" N_L | strands max_len  lead-foll_sep  aspect  LCC   detached  mean_inv max_inv  fpair")
    for nl in (0, 1, 2, 3, 5):
        cfg = r1_config(n_leaders=nl, leader_location="cue_front", leader_stall_factor=4.0,
                        budget_mode="matched_total", **SWEEP)
        out = run_r2_invasion(cfg, seed=23)
        st, ml, lfs, ar, lcc, det, mi, mx, fp = _row(out)
        print(" %3d | %6d %6d  %+12.3f  %5.2f  %5.3f  %6.3f  %+7.3f  %6.3f  %.0e" % (
            nl, st, ml, lfs, ar, lcc, det, mi, mx, fp))
        del out; gc.collect()
    print("matched-total: total front stall held constant vs N_L -> isolates leader NUMBER/geometry.")


def budget_compare():
    print("\n=== R2 (b) budget mode @ N_L=3 (cued, cue_front, factor 4, 900 s) ===")
    print(" mode             | strands max_len  lead-foll_sep  mean_inv max_inv")
    for mode in ("matched_total", "fixed_per_leader"):
        cfg = r1_config(n_leaders=3, leader_location="cue_front", leader_stall_factor=4.0,
                        budget_mode=mode, **SWEEP)
        out = run_r2_invasion(cfg, seed=23)
        st, ml, lfs, ar, lcc, det, mi, mx, fp = _row(out)
        print(" %-16s | %6d %6d  %+12.3f  %+7.3f  %6.3f" % (mode, st, ml, lfs, mi, mx))
        del out; gc.collect()


def location_compare():
    print("\n=== R2 (c) leader LOCATION @ N_L=3 (cued, fixed_per_leader, factor 4, 900 s) ===")
    print(" location   | strands max_len  lead-foll_sep  aspect  mean_inv max_inv")
    for loc in ("cue_front", "perimeter"):
        cfg = r1_config(n_leaders=3, leader_location=loc, leader_stall_factor=4.0,
                        budget_mode="fixed_per_leader", **SWEEP)
        out = run_r2_invasion(cfg, seed=23)
        st, ml, lfs, ar, lcc, det, mi, mx, fp = _row(out)
        print(" %-10s | %6d %6d  %+12.3f  %5.2f  %+7.3f  %6.3f" % (loc, st, ml, lfs, ar, mi, mx))
        del out; gc.collect()


def _radial_order(pos, edges, center):
    seg = pos[edges[:, 1]] - pos[edges[:, 0]]
    tan = seg / np.maximum(np.linalg.norm(seg, axis=1), 1e-12)[:, None]
    mid = 0.5 * (pos[edges[:, 0]] + pos[edges[:, 1]])
    rad = mid - center
    er = rad / np.maximum(np.linalg.norm(rad, axis=1), 1e-12)[:, None]
    return 2.0 * np.square(np.sum(tan * er, axis=1)) - 1.0


def animate_leaders(bead_snaps, cell_snaps, edges, leaders, out_path, *, cell_radius=9.0,
                    span=None, fps=8, title=""):
    """Stage-D gif with LEADER cells red, followers blue (so the localized front is visible)."""
    bead_snaps = np.asarray(bead_snaps); cell_snaps = np.asarray(cell_snaps)
    c0 = cell_snaps[0]; center = np.zeros(2)
    lead = set(int(i) for i in leaders)
    if span is None:
        span = float(np.max(np.abs(bead_snaps[0]))) * 1.02
    fig, ax = plt.subplots(figsize=(6.6, 6.6))

    def draw(k):
        ax.clear()
        pos = bead_snaps[k]
        order = _radial_order(pos, edges, center)
        segs = np.stack([pos[edges[:, 0]], pos[edges[:, 1]]], axis=1)
        lc = LineCollection(segs, cmap="coolwarm", norm=plt.Normalize(-1, 1), linewidths=0.5, alpha=0.8)
        lc.set_array(order); ax.add_collection(lc)
        cells = cell_snaps[k]
        ax.add_collection(LineCollection(np.stack([c0, cells], axis=1), colors="#111", linewidths=0.6, alpha=0.35))
        for i, c in enumerate(cells):
            ax.add_patch(plt.Circle(c, cell_radius, color=("#d62728" if i in lead else "#2c7fb8"),
                                    alpha=0.9 if i in lead else 0.8, lw=0))
        invaded = float(np.mean(np.linalg.norm(cells, axis=1) - np.linalg.norm(c0, axis=1)))
        ax.set_xlim(-span, span); ax.set_ylim(-span, span); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("%s\nframe %d/%d  red=leaders  invaded %+.1f um" % (title, k + 1, len(bead_snaps), invaded), fontsize=9)

    FuncAnimation(fig, draw, frames=len(bead_snaps), interval=1000 / fps).save(out_path, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return out_path


def _headline_one(tag, title, outdir, **over):
    base = dict(n_fibers=400, organoid_radius=40.0, domain_size=340.0, boundary_width=6.0,
                generation_attempts=20, dt=0.05, duration=7200.0, sample_interval=60.0,
                radial_cue=True, cue_angle=0.0, cc_adhesion=6.0)
    cfg = r1_config(**{**base, **over})
    t0 = time.time()
    out = run_r2_invasion(cfg, seed=23, snapshots=True)
    f = out["frames"][-1]
    print("[%s] run %.0fs | leaders %s | fpair %.0e | strands %d max_len %d | lead-foll_sep %+.2f | aspect %.2f | mean_inv %+.2f max %.2f um" % (
        tag, time.time() - t0, out["leader_ids"], f["force_pair_residual"], out["strand_count"],
        out["strand_max_len"], out["leader_follower_separation"], out["aspect_ratio"],
        f["mean_cell_radial_disp"], f["max_cell_disp"]))
    edges = out["edges"]
    np.savez_compressed(f"{outdir}/r2_{tag}_2h.npz", cell_snapshots=out["cell_snapshots"],
                        edges=edges, centers0=out["centers0"], centers_final=out["centers_final"],
                        leader_ids=np.asarray(out["leader_ids"]))
    span = float(np.max(np.abs(out["bead_snapshots"][0]))) * 1.02
    gif = animate_leaders(out["bead_snapshots"], out["cell_snapshots"], edges, out["leader_ids"],
                          f"{outdir}/r2_{tag}_2h.gif", cell_radius=cfg.cell_radius, span=span, title=title)
    print("wrote", gif)
    del out; gc.collect()


def headline(outdir="output"):
    Path(outdir).mkdir(exist_ok=True)
    _headline_one("leader_front", "G5-R2 3 localized cue-front leaders (per-site stall x6) -- finger (personal testing)",
                  outdir, n_leaders=3, leader_location="cue_front", leader_stall_factor=6.0,
                  budget_mode="fixed_per_leader")
    _headline_one("no_leader", "G5-R2 no leaders (control, cued ECM) -- uniform front (personal testing)",
                  outdir, n_leaders=0)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "count"
    {"count": count_sweep, "budget": budget_compare,
     "location": location_compare, "headline": headline}.get(mode, count_sweep)()
