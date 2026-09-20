"""Re-run the G5 R0-R3 headlines in Gloria g4-v4e VISUAL STYLE (consistency layer).

Same physics / configs as the committed R0-R3 headline experiments, but at G4D-parity crosslink
density (``XL_FRAC``) and rendered with ``gloria_style_viz`` (neutral grey fibres, gold crosslinks,
open-circle cells with red partial-EMT / leaders, capped traction arrows, red x grip failures) +
a pull-slip-rebind EVENT MICROSCOPE per run.

TWO-PHASE (a 2 h sim must never be lost to a render bug):
    python .../gloria_headlines.py sim    r1     # run sim(s), save gv_<tag>_full.npz (rich, float32)
    python .../gloria_headlines.py render r1     # load npz -> main gif + event gif + last PNG
    python .../gloria_headlines.py all    r1     # sim then render
Stages: r0 | r1 (cohesive+escape) | r2 (leader_front+no_leader) | r3.  One stage per process (memory).

Crosslink note: G4D uses ~1.7 crosslinks/fibre; our larger organoid void trips the conservative
binary contact_conn flag at that density, but the floppiness/force-pair test showed the mechanics
stay clean, so XL_FRAC=0.62 -> 1.66/fibre at the 400/40/340 scale = G4D parity.
Personal testing, NOT confirmed findings (CLAUDE.md 7.5).  No swirling claim implied.
"""

from __future__ import annotations

from pathlib import Path
import gc
import os
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.consistency import (  # noqa: E402
    run_r0_invasion, run_r1_invasion, run_r2_invasion, run_r3_invasion, run_r4_invasion,
    r0_config, r1_config, emt_cell_ids,
)
from generations.g5_organoid.experiments.gloria_style_viz import (  # noqa: E402
    render_gloria_style, render_event_microscope, render_frame_png, run_out_from_npz,
)

XL_FRAC = 0.62  # G4D parity at the headline scale (see module docstring + floppiness decision)
OUT = REPO / "output"
BASE2H = dict(n_fibers=400, organoid_radius=40.0, domain_size=340.0, boundary_width=6.0,
              generation_attempts=20, dt=0.05, duration=7200.0, sample_interval=60.0,
              crosslink_fraction=XL_FRAC)
if os.environ.get("GV_SMOKE"):   # fast full-pipeline check (sim->npz->render), NOT a headline
    BASE2H.update(n_fibers=200, organoid_radius=24.0, domain_size=230.0, generation_attempts=15,
                  duration=90.0, sample_interval=9.0)


def _save_full(tag, out, *, red_ids=(), red_per_frame=None, red_label, title, cell_radius):
    """Persist EVERYTHING the renderer needs (float32) so rendering is decoupled from the sim."""
    OUT.mkdir(exist_ok=True)
    f = out["frames"][-1]
    grip = np.asarray(out["site_failed_snapshots"])
    payload = dict(
        bead_snapshots=np.asarray(out["bead_snapshots"], dtype=np.float32),
        cell_snapshots=np.asarray(out["cell_snapshots"], dtype=np.float32),
        edges=np.asarray(out["edges"], dtype=np.int32),
        site_force_snapshots=np.asarray(out["site_force_snapshots"], dtype=np.float32),
        site_point_snapshots=np.asarray(out["site_point_snapshots"], dtype=np.float32),
        site_normal_snapshots=np.asarray(out["site_normal_snapshots"], dtype=np.float32),
        site_failed_snapshots=grip,
        crosslink_edge=np.asarray(out["crosslink_edge"], dtype=np.int32),
        crosslink_alpha=np.asarray(out["crosslink_alpha"], dtype=np.float32),
        times=np.asarray([fr["time"] for fr in out["frames"]], dtype=np.float32),
        n_contact_sectors=np.int32(out["n_contact_sectors"]),
        red_ids=np.asarray(list(red_ids), dtype=np.int32),
        cell_radius=np.float32(cell_radius), title=np.array(title), red_label=np.array(red_label),
        mean_inv=np.float32(f["mean_cell_radial_disp"]), max_inv=np.float32(f["max_cell_disp"]),
        detached=np.float32(f.get("detached_fraction", np.nan)),
        rg=np.float32(f.get("radius_of_gyration", np.nan)),
        lcc=np.float32(f.get("lcc_fraction", np.nan)),
        fpair=np.float32(f["force_pair_residual"]), grip_fails=np.int32(grip.sum()),
    )
    if red_per_frame is not None:
        payload["red_per_frame"] = np.array([list(r) for r in red_per_frame], dtype=object)
    path = OUT / f"gv_{tag}_full.npz"
    np.savez_compressed(path, **payload)
    print("[sim %s] saved %s | fpair %.0e mean_inv %+.2f max %.2f detached %.3f Rg %.2f LCC %.2f grip-fails %d"
          % (tag, path.name, f["force_pair_residual"], f["mean_cell_radial_disp"], f["max_cell_disp"],
             f.get("detached_fraction", float("nan")), f.get("radius_of_gyration", float("nan")),
             f.get("lcc_fraction", float("nan")), int(grip.sum())))
    return path


# ------------------------------- SIM PHASE -------------------------------
def sim_r0():
    t0 = time.time(); cfg = r0_config(**BASE2H)
    out = run_r0_invasion(cfg, seed=23, snapshots=True)
    _save_full("r0", out, red_label="(none)", cell_radius=cfg.cell_radius,
               title="G5-R0 force-consistent invasion (isotropic-random, G4D-parity crosslinks) -- personal testing")
    print("R0 sim %.0fs" % (time.time() - t0)); del out; gc.collect()


def sim_r1():
    for tag, over, title in [
        ("r1_cohesive", dict(cc_adhesion=6.0, emt_fraction=0.0),
         "G5-R1 cohesive (no EMT) -- collective front (personal testing)"),
        ("r1_escape", dict(cc_adhesion=6.0, emt_fraction=0.5, emt_adhesion_factor=0.15),
         "G5-R1 partial-EMT (f=0.5, adh x0.15) -> single-cell escape (personal testing)"),
    ]:
        t0 = time.time(); cfg = r1_config(**{**BASE2H, **over})
        out = run_r1_invasion(cfg, seed=23, snapshots=True)
        red = list(emt_cell_ids(out["centers0"], cfg)) if over["emt_fraction"] > 0 else []
        _save_full(tag, out, red_ids=red, red_label="partial-EMT", cell_radius=cfg.cell_radius, title=title)
        print("R1[%s] sim %.0fs" % (tag, time.time() - t0)); del out; gc.collect()


def sim_r2():
    for tag, over, title in [
        ("r2_leader_front", dict(radial_cue=True, cue_angle=0.0, cc_adhesion=6.0, n_leaders=3,
                                 leader_location="cue_front", leader_stall_factor=6.0, budget_mode="fixed_per_leader"),
         "G5-R2 3 cue-front leaders (per-site stall x6) -- finger (personal testing)"),
        ("r2_no_leader", dict(radial_cue=True, cue_angle=0.0, cc_adhesion=6.0, n_leaders=0),
         "G5-R2 no leaders (control, cued ECM) -- uniform front (personal testing)"),
    ]:
        t0 = time.time(); cfg = r1_config(**{**BASE2H, **over})
        out = run_r2_invasion(cfg, seed=23, snapshots=True)
        red = list(out.get("leader_ids", []))
        _save_full(tag, out, red_ids=red, red_label="leader", cell_radius=cfg.cell_radius, title=title)
        print("R2[%s] sim %.0fs | leaders %s" % (tag, time.time() - t0, red)); del out; gc.collect()


def sim_r3():
    t0 = time.time()
    cfg = r1_config(radial_cue=True, cue_angle=0.0, cc_adhesion=2.0, n_leaders=2,
                    leader_location="cue_front", leader_stall_factor=6.0,
                    budget_mode="fixed_per_leader", leader_switching=True, **BASE2H)
    out = run_r3_invasion(cfg, seed=23, snapshots=True)
    red_pf = [fr["active_leaders"] for fr in out["frames"]]
    _save_full("r3_switching", out, red_per_frame=red_pf, red_label="CURRENT leader (baton)",
               cell_radius=cfg.cell_radius,
               title="G5-R3 energy-based leader relay (red=current leader; ATP proxy=hypothesis) -- personal testing")
    print("R3 sim %.0fs | switches %d" % (time.time() - t0, out["n_switches"])); del out; gc.collect()


_R4 = dict(radial_cue=True, cue_angle=0.0, cc_adhesion=1.5, n_leaders=3,
           leader_location="cue_front", leader_stall_factor=6.0, budget_mode="fixed_per_leader")
_R4_UNG = dict(follower_guidance=False, guidance_strength=0.0)
_R4_GUID = dict(follower_guidance=True, guidance_strength=20.0, guidance_range=25.0, guidance_align_min=0.30)


def _sim_r4(xl_frac, prefix, dens_label):
    """R4 (cc=1.5 low-adhesion + cue-front leaders) guided vs unguided at crosslink density xl_frac.
    The follower-guidance benefit is CROSSLINK-DEPENDENT: dense (0.85) leaders detach -> guidance
    rescues a follower group; sparse G4D-parity (0.62) the organoid holds -> guidance over-scatters."""
    for suffix, over, story in [
        ("unguided", _R4_UNG, "NO guidance"), ("guided", _R4_GUID, "+ follower contact-guidance")]:
        tag = f"{prefix}{suffix}"
        t0 = time.time()
        cfg = r1_config(**{**BASE2H, "crosslink_fraction": xl_frac, **_R4, **over})
        out = run_r4_invasion(cfg, seed=23, snapshots=True)
        red = list(out.get("leader_ids", []))
        title = f"G5-R4 {dens_label} ECM (xl {xl_frac}) + leaders, {story} (personal testing)"
        _save_full(tag, out, red_ids=red, red_label="leader", cell_radius=cfg.cell_radius, title=title)
        print("R4[%s] sim %.0fs | leaders %s | leader-comp followers %s | detached %s | LCC %.2f" % (
            tag, time.time() - t0, red, out.get("leader_comp_followers"), out.get("leaders_detached"),
            out["frames"][-1].get("lcc_fraction", float("nan"))))
        del out; gc.collect()


def sim_r4():        # sparse = G4D-parity 0.62 (tags r4_unguided / r4_guided)
    _sim_r4(0.62, "r4_", "sparse G4D-parity")


def sim_r4dense():   # dense 0.85 = the regime where leaders detach and R4 guidance rescues them
    _sim_r4(0.85, "r4_dense_", "dense")


# ------------------------------ RENDER PHASE ------------------------------
def _render_tag(tag):
    path = OUT / f"gv_{tag}_full.npz"
    if not path.exists():
        print("[render %s] MISSING %s -- run sim first" % (tag, path.name)); return
    out, meta = run_out_from_npz(path)
    span = float(np.max(np.abs(out["bead_snapshots"][0]))) * 1.02
    R = meta["cell_radius"]; rp = meta["red_per_frame"]; ri = meta["red_ids"]
    # CLEANEST look (per Minnie): light fibres + cells + displacement trails, NO arrows at all,
    # no crosslink dots, no per-site failure marks, no event microscope.
    render_frame_png(out, str(OUT / f"gv_{tag}_last.png"), red_ids=ri, red_per_frame=rp,
                     red_label=meta["red_label"], cell_radius=R, span=span, frame=-1,
                     title=meta["title"], arrows="none")
    render_gloria_style(out, str(OUT / f"gv_{tag}_2h.gif"), red_ids=ri, red_per_frame=rp,
                        red_label=meta["red_label"], cell_radius=R, span=span, fps=10,
                        title=meta["title"], arrows="none")
    s = meta["summary"]
    print("[render %s] clean gif+png | fpair %.0e mean_inv %+.2f max %.2f detached %.2f Rg %.2f LCC %.2f grip-fails %d"
          % (tag, s.get("fpair", float("nan")), s.get("mean_inv", float("nan")),
             s.get("max_inv", float("nan")), s.get("detached", float("nan")), s.get("rg", float("nan")),
             s.get("lcc", float("nan")), s.get("grip_fails", -1)))


STAGES = {"r0": (sim_r0, ["r0"]), "r1": (sim_r1, ["r1_cohesive", "r1_escape"]),
          "r2": (sim_r2, ["r2_leader_front", "r2_no_leader"]), "r3": (sim_r3, ["r3_switching"]),
          "r4": (sim_r4, ["r4_unguided", "r4_guided"]),
          "r4dense": (sim_r4dense, ["r4_dense_unguided", "r4_dense_guided"])}


if __name__ == "__main__":
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    stage = sys.argv[2] if len(sys.argv) > 2 else "r1"
    sim_fn, tags = STAGES[stage]
    if phase in ("sim", "all"):
        sim_fn()
    if phase in ("render", "all"):
        for t in tags:
            _render_tag(t)
