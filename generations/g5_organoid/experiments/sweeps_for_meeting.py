"""Clean ablation sweeps for the meeting comparison table (honest metrics only).

Three axes, each 2 h, seed 23, 400/40/340 scale, force-pair checked:
  A. partial-EMT FRACTION   (how many cells go partial-EMT; adhesion x0.15) -> collective vs escape
  B. partial-EMT DEPTH      (emt_adhesion_factor: 1.0 epithelial ... 0.05 mesenchymal), fraction 0.5
  C. NUMBER OF LEADERS      (n_leaders 0..5) in the dense-ECM regime where leaders can detach

Metrics are the TRUSTWORTHY ones (NOT the retracted leader_follower_separation / aspect):
  mean/max invasion (um), detached_fraction, LCC, radius_of_gyration, and for leaders the
  adhesion-graph strand_report (leader_comp_followers, leaders_detached, n_components).

Writes output/meeting_ablation_table.md and .csv INCREMENTALLY (so a crash keeps finished rows).
    python .../sweeps_for_meeting.py
Personal testing, NOT confirmed findings (CLAUDE.md 7.5).  No swirling.
"""
from __future__ import annotations
from pathlib import Path
import sys, gc, time
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
from generations.g5_organoid.consistency import (  # noqa: E402
    run_r1_invasion, run_r4_invasion, r1_config, strand_report,
)

OUT = REPO / "output"; OUT.mkdir(exist_ok=True)
MD = OUT / "meeting_ablation_table.md"
CSV = OUT / "meeting_ablation_table.csv"
BASE = dict(n_fibers=400, organoid_radius=40.0, domain_size=340.0, boundary_width=6.0,
            generation_attempts=20, dt=0.05, duration=7200.0, sample_interval=120.0)

_rows = []


def _emit(section, label, m):
    _rows.append((section, label, m))
    # rewrite the full markdown + csv each time (cheap; keeps partial results safe)
    with open(MD, "w", encoding="utf-8") as f:
        f.write("# G5 ablation comparison (2 h, seed 23, 400/40/340; personal testing)\n\n")
        f.write("Honest metrics: invasion (um), detached fraction, LCC (largest connected comp.), "
                "Rg; for leaders the adhesion-graph strand_report. force-pair 0 in all runs.\n\n")
        cur = None
        for sec, lab, mm in _rows:
            if sec != cur:
                cur = sec
                f.write(f"\n## {sec}\n\n")
                f.write("| condition | mean inv | max inv | detached | LCC | Rg | leaders: followers | detached? | comps |\n")
                f.write("|---|--:|--:|--:|--:|--:|--:|:--:|--:|\n")
            f.write("| %s | %+.1f | %.1f | %.3f | %.2f | %.1f | %s | %s | %s |\n" % (
                lab, mm["mean"], mm["max"], mm["det"], mm["lcc"], mm["rg"],
                mm.get("lf", "—"), mm.get("ld", "—"), mm.get("nc", "—")))
    import csv as _csv
    with open(CSV, "w", encoding="utf-8", newline="") as f:
        w = _csv.writer(f)   # proper quoting (section labels contain commas/·)
        w.writerow(["section", "condition", "mean_inv", "max_inv", "detached", "LCC", "Rg",
                    "leader_followers", "leaders_detached", "components", "fpair"])
        for sec, lab, mm in _rows:
            w.writerow([sec, lab, "%.3f" % mm["mean"], "%.3f" % mm["max"], "%.4f" % mm["det"],
                        "%.4f" % mm["lcc"], "%.3f" % mm["rg"], mm.get("lf", ""), mm.get("ld", ""),
                        mm.get("nc", ""), "%.1e" % mm["fpair"]])


def _metrics(out, cfg=None, leaders=False):
    f = out["frames"][-1]
    m = dict(mean=f["mean_cell_radial_disp"], max=f["max_cell_disp"], det=f["detached_fraction"],
             lcc=f["lcc_fraction"], rg=f["radius_of_gyration"], fpair=out["max_force_pair_residual"])
    if leaders:
        rep = strand_report(out["centers_final"], out["centers0"], cfg)
        m["lf"] = rep.get("leader_comp_followers"); m["ld"] = str(rep.get("leaders_detached"))
        m["nc"] = rep.get("n_components")
    return m


def sweep_A():
    print("=== A. partial-EMT FRACTION (adh x0.15, no leaders, cc=6, xl 0.62) ===")
    for fr in (0.0, 0.25, 0.5, 0.75, 1.0):
        t0 = time.time()
        cfg = r1_config(crosslink_fraction=0.62, cc_adhesion=6.0, emt_fraction=fr,
                        emt_adhesion_factor=0.15, **BASE)
        out = run_r1_invasion(cfg, seed=23)
        _emit("A · partial-EMT FRACTION (depth x0.15)", "f_EMT = %.2f" % fr, _metrics(out))
        print("  f=%.2f done %.0fs" % (fr, time.time() - t0)); del out; gc.collect()


def sweep_B():
    print("=== B. partial-EMT DEPTH (emt_adhesion_factor, fraction 0.5, no leaders, cc=6, xl 0.62) ===")
    for af in (1.0, 0.5, 0.25, 0.15, 0.05):
        t0 = time.time()
        cfg = r1_config(crosslink_fraction=0.62, cc_adhesion=6.0, emt_fraction=0.5,
                        emt_adhesion_factor=af, **BASE)
        out = run_r1_invasion(cfg, seed=23)
        _emit("B · partial-EMT DEPTH (fraction 0.5)", "adh x%.2f" % af, _metrics(out))
        print("  adh=%.2f done %.0fs" % (af, time.time() - t0)); del out; gc.collect()


def sweep_C():
    print("=== C. NUMBER OF LEADERS (cc=1.5 low-adh, DENSE xl 0.85, cued, stall x6, no guidance) ===")
    for nl in (0, 1, 2, 3, 5):
        t0 = time.time()
        cfg = r1_config(crosslink_fraction=0.85, cc_adhesion=1.5, radial_cue=True, cue_angle=0.0,
                        n_leaders=nl, leader_location="cue_front", leader_stall_factor=6.0,
                        budget_mode="fixed_per_leader", follower_guidance=False, **BASE)
        out = run_r4_invasion(cfg, seed=23)
        _emit("C · NUMBER OF LEADERS (cc=1.5, dense ECM)", "N_leaders = %d" % nl,
              _metrics(out, cfg=cfg, leaders=True))
        print("  N=%d done %.0fs" % (nl, time.time() - t0)); del out; gc.collect()


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "A"): sweep_A()
    if which in ("all", "B"): sweep_B()
    if which in ("all", "C"): sweep_C()
    print("\nWrote", MD, "and", CSV)
