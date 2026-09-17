"""G5-R3 experiment: energy-based DYNAMIC leader switching = PERSISTENCE.

Does a relay of leaders (tired front leader -> demote -> fresh follower takes over; Zhang
2019) keep the invasion front ADVANCING vs a STATIC leader (R2)?  Uses the committed R3
driver (consistency.run_r3_invasion) on the R1xR2 sweet spot (cued, low adhesion, few
localized leaders).  The gif colours the CURRENT active leaders each frame, so the baton
hand-off is visible.

    python generations/g5_organoid/experiments/r3_switching.py headline   # 2 h ON vs OFF + gif

Personal testing, NOT confirmed findings (CLAUDE.md 7.5).  Motor-work = ATP proxy is a
HYPOTHESIS, not a measured law.  Imposed cue, NOT swirling.
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

from generations.g5_organoid.consistency import run_r3_invasion, r1_config  # noqa: E402
from generations.g5_organoid.experiments.r2_leader import _radial_order  # noqa: E402

# R1xR2 sweet spot: cued, low adhesion, few localized cue-front leaders; radius>=40 so the
# front pool has >=3 gripping cells (needed for a real relay, not a self-cycle).
BASE = dict(n_fibers=400, organoid_radius=40.0, domain_size=340.0, boundary_width=6.0,
            generation_attempts=20, dt=0.05, duration=7200.0, sample_interval=60.0,
            radial_cue=True, cue_angle=0.0, cc_adhesion=2.0, n_leaders=2,
            leader_location="cue_front", leader_stall_factor=6.0, budget_mode="fixed_per_leader")


def animate_switching(bead_snaps, cell_snaps, edges, leaders_per_frame, out_path, *,
                      cell_radius=9.0, span=None, fps=8, title=""):
    """Stage-D gif; the CURRENT active leaders (per frame) are red -> the relay is visible."""
    bead_snaps = np.asarray(bead_snaps); cell_snaps = np.asarray(cell_snaps)
    c0 = cell_snaps[0]; center = np.zeros(2)
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
        lead = set(int(i) for i in leaders_per_frame[k])
        for i, c in enumerate(cells):
            ax.add_patch(plt.Circle(c, cell_radius, color=("#d62728" if i in lead else "#2c7fb8"),
                                    alpha=0.9 if i in lead else 0.8, lw=0))
        ax.set_xlim(-span, span); ax.set_ylim(-span, span); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("%s\nframe %d/%d  red=CURRENT leaders (baton passes)  leaders now: %s"
                     % (title, k + 1, len(bead_snaps), sorted(lead)), fontsize=8)

    FuncAnimation(fig, draw, frames=len(bead_snaps), interval=1000 / fps).save(out_path, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return out_path


def headline(outdir="output"):
    Path(outdir).mkdir(exist_ok=True)
    # ON = dynamic switching; OFF = static leader (R2) via leader_switching=False.
    t0 = time.time()
    # ON = relay (energy drains -> hand-offs). OFF = STATIC control through the SAME R3 loop
    # with energy_cap huge (energy never drains -> the initial leaders never demote), so it is
    # apples-to-apples (identical code path; only the switching differs) and records front_advance.
    on = run_r3_invasion(r1_config(leader_switching=True, **BASE), seed=23, snapshots=True)
    off = run_r3_invasion(r1_config(leader_switching=True, energy_cap=1e12, **BASE), seed=23, snapshots=False)
    fon, foff = on["frames"], off["frames"]
    t = np.array([fr["time"] for fr in fon])
    adv_on = np.array([fr.get("front_advance", 0.0) for fr in fon])
    adv_off = np.array([fr.get("front_advance", 0.0) for fr in foff])
    nsw = fon[-1]["n_switches"]
    distinct = sorted({int(c) for fr in fon for c in fr["active_leaders"]})
    print("run %.0fs" % (time.time() - t0))
    print("SWITCHING ON : front_advance %+.2f um | n_switches %d | distinct leaders %s | fpair %.0e | max_inv %.2f" % (
        adv_on[-1], nsw, distinct, on["max_force_pair_residual"], fon[-1]["max_cell_disp"]))
    print("STATIC   OFF : front_advance %+.2f um | (1 fixed leader set) | fpair %.0e | max_inv %.2f" % (
        adv_off[-1], off["max_force_pair_residual"], foff[-1]["max_cell_disp"]))

    # persistence plot: front advance (ON vs OFF) + leader energy & switch markers
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    ax[0].plot(t / 60, adv_on, "-", color="#d62728", label="switching ON (relay)")
    ax[0].plot(t / 60, adv_off, "--", color="#2c7fb8", label="static leader (OFF)")
    ax[0].set_xlabel("time (min)"); ax[0].set_ylabel("front advance along cue (um)")
    ax[0].set_title("Persistence: does the relay keep the front advancing?"); ax[0].legend(); ax[0].grid(alpha=0.3)
    le = np.array([fr["leader_energy_mean"] for fr in fon])
    ns = np.array([fr["n_switches"] for fr in fon])
    ax[1].plot(t / 60, le, "-", color="#d62728", label="active-leader energy")
    for i in range(1, len(ns)):
        if ns[i] > ns[i - 1]:
            ax[1].axvline(t[i] / 60, color="#888", lw=0.8, alpha=0.7)
    ax[1].axhline(BASE_off := 0.3, color="k", ls=":", lw=0.8, label="energy_off (demote)")
    ax[1].set_xlabel("time (min)"); ax[1].set_ylabel("leader energy"); ax[1].set_ylim(0, 1.05)
    ax[1].set_title("Leader energy drains -> hand-off (grey lines = switches)"); ax[1].legend(); ax[1].grid(alpha=0.3)
    fig.suptitle("G5-R3 energy-based dynamic leader switching (personal testing; ATP proxy = hypothesis)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(f"{outdir}/r3_persistence.png", dpi=130, bbox_inches="tight")
    plt.close(fig); print("wrote", f"{outdir}/r3_persistence.png")

    leaders_per_frame = [fr["active_leaders"] for fr in fon]
    edges = on["edges"]
    np.savez_compressed(f"{outdir}/r3_switching_2h.npz", cell_snapshots=on["cell_snapshots"], edges=edges,
                        centers0=on["centers0"], centers_final=on["centers_final"],
                        leaders_per_frame=np.array(leaders_per_frame, dtype=object),
                        front_advance_on=adv_on, front_advance_off=adv_off, times=t)
    span = float(np.max(np.abs(on["bead_snapshots"][0]))) * 1.02
    gif = animate_switching(on["bead_snapshots"], on["cell_snapshots"], edges, leaders_per_frame,
                            f"{outdir}/r3_switching_2h.gif", cell_radius=9.0, span=span,
                            title="G5-R3 relay: tired leader hands off to a fresh one (personal testing)")
    print("wrote", gif)
    del on, off; gc.collect()


if __name__ == "__main__":
    (headline if (len(sys.argv) > 1 and sys.argv[1] == "headline") else headline)()
