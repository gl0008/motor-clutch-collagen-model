"""G5-R4: does FOLLOWER contact-guidance (followers migrate ALONG the leader-remodelled
aligned collagen track) turn "leaders detach + followers scatter" into a connected,
directional leader-follower STRAND?

Diagnosis it fixes (R1xR2 / R3): with low cell-cell adhesion the strong leaders DETACH and
escape while followers scatter randomly -- the model had NO mechanism making followers follow.
Biology (deep-research synthesis): leaders specialise in CREATING an aligned ECM path
(force-bundled radial fibres / TACS-3 + proteolytic microtracks); followers specialise in
RESPONDING to it by contact guidance / cryptic-lamellipodia migration ALONG the path (Friedl &
Gilmour 2009; Ray 2017 Biophys J; PLOS One 2024 leader-aligned-collagen; Nat Rev Cancer 2021).
R4 = R2 localized per-site-stall leaders + `contact_guidance_forces` on the followers.

Honest metric: `leader_follower_separation` / `aspect_ratio` CANNOT tell a connected strand from
leaders escaping alone (a large value arises even when followers move BACKWARD).  We therefore
read the ADHESION-GRAPH connectivity (`strand_report`): n_components, whether the leaders' cluster
still carries followers (`leader_comp_followers`, `leaders_detached`), and its outward reach.

    python .../r4_follower_guidance.py compare    # guidance off vs on x adhesion (1500 s)
    python .../r4_follower_guidance.py headline    # 2 h gifs: guided vs unguided (leaders red)

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

from generations.g5_organoid.consistency import run_r4_invasion, r1_config  # noqa: E402
from generations.g5_organoid.experiments.r2_leader import animate_leaders  # noqa: E402

CUED_LEAD = dict(n_fibers=300, organoid_radius=32.0, domain_size=280.0, boundary_width=6.0,
                 generation_attempts=20, dt=0.05, radial_cue=True, cue_angle=0.0,
                 n_leaders=3, leader_location="cue_front", leader_stall_factor=6.0,
                 budget_mode="fixed_per_leader")
GUID = dict(follower_guidance=True, guidance_strength=20.0, guidance_range=25.0,
            guidance_align_min=0.30)


def compare():
    print("\n=== R4: follower contact-guidance OFF vs ON (cued, 3 cue-front leaders x6, 1500 s) ===")
    print(" cc  guid | fpair | sep aspect | comps lgfrac Lcomp_foll detached reach | meaninv maxinv LCC det")
    for cc in (1.5, 4.0):
        for guid in (False, True):
            over = dict(CUED_LEAD)
            over.update(cc_adhesion=cc, duration=1500.0, sample_interval=150.0)
            over.update(GUID if guid else dict(follower_guidance=False, guidance_strength=0.0))
            out = run_r4_invasion(r1_config(**over), seed=23)
            f = out["frames"][-1]
            print(" %3.1f %5s | %.0e | %+5.2f %5.2f | %5d %6.2f %6d %8s %5.1f | %+6.2f %6.2f %.2f %.2f" % (
                cc, guid, out["max_force_pair_residual"], out["leader_follower_separation"],
                out["aspect_ratio"], out["n_components"], out["largest_comp_fraction"],
                out["leader_comp_followers"], out["leaders_detached"], out["strand_reach_um"],
                f["mean_cell_radial_disp"], f["max_cell_disp"], f["lcc_fraction"], f["detached_fraction"]))
            del out
            gc.collect()
    print("Q: does guidance keep followers CONNECTED to the leaders (leader_comp_followers>0,")
    print("   leaders_detached=False, fewer components) instead of leaders escaping alone?")


def _one(tag, title, outdir, **over):
    base = dict(n_fibers=400, organoid_radius=40.0, domain_size=340.0, boundary_width=6.0,
                generation_attempts=20, dt=0.05, duration=7200.0, sample_interval=60.0,
                radial_cue=True, cue_angle=0.0, n_leaders=3, leader_location="cue_front",
                leader_stall_factor=6.0, budget_mode="fixed_per_leader")
    cfg = r1_config(**{**base, **over})
    t0 = time.time()
    out = run_r4_invasion(cfg, seed=23, snapshots=True)
    f = out["frames"][-1]
    print("[%s] run %.0fs | leaders %s cc %.1f guid %s | fpair %.0e | comps %d Lcomp_foll %d detached %s reach %.1f | sep %+.2f aspect %.2f | mean_inv %+.2f max %.2f | LCC %.2f det %.2f" % (
        tag, time.time() - t0, out["leader_ids"], cfg.cc_adhesion, cfg.follower_guidance,
        f["force_pair_residual"], out["n_components"], out["leader_comp_followers"],
        out["leaders_detached"], out["strand_reach_um"], out["leader_follower_separation"],
        out["aspect_ratio"], f["mean_cell_radial_disp"], f["max_cell_disp"],
        f["lcc_fraction"], f["detached_fraction"]))
    edges = out["edges"]
    np.savez_compressed(f"{outdir}/r4_{tag}_2h.npz", cell_snapshots=out["cell_snapshots"],
                        edges=edges, centers0=out["centers0"], centers_final=out["centers_final"],
                        leader_ids=np.asarray(out["leader_ids"]))
    span = float(np.max(np.abs(out["bead_snapshots"][0]))) * 1.02
    gif = animate_leaders(out["bead_snapshots"], out["cell_snapshots"], edges, out["leader_ids"],
                          f"{outdir}/r4_{tag}_2h.gif", cell_radius=cfg.cell_radius, span=span, title=title)
    print("wrote", gif)
    del out
    gc.collect()


def headline(outdir="output"):
    Path(outdir).mkdir(exist_ok=True)
    _one("guided", "G5-R4: low adhesion + 3 leaders + FOLLOWER contact-guidance (personal testing)",
         outdir, cc_adhesion=1.5, **GUID)
    _one("unguided", "G5-R4 control: low adhesion + 3 leaders, NO guidance (personal testing)",
         outdir, cc_adhesion=1.5, follower_guidance=False, guidance_strength=0.0)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "compare"
    {"compare": compare, "headline": headline}.get(mode, compare)()
