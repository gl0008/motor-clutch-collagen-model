"""Build a clearer G5 demo: longer time, denser corona, animated GIF + figure.

Runs a fixed contractile organoid, saves per-frame snapshots, and renders both a
before/after PNG and an animated GIF (fibres coloured by radial order, crosslinks
drawn in yellow, near-field zoom).  Personal testing, not a validated finding.

    python generations/g5_organoid/build_demo.py
"""

from __future__ import annotations

from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.model import OrganoidConfig, run_organoid_pull  # noqa: E402
from generations.g5_organoid import visualize as viz  # noqa: E402


# Demo config: denser grippable corona + stronger pull + longer time so the
# near-field aster is clearly visible.  Stage-C stiffening is OFF (it barely moves
# this gentle-contraction regime; keep it as an ablation flag).
DEMO = OrganoidConfig(
    n_fibers=460,
    n_corona_fibers=170,
    corona_band=10.0,
    gap=1.0,
    organoid_radius=63.0,
    domain_size=440.0,
    total_pull_force=20.0,
    duration=300.0,
    sample_interval=10.0,     # 31 frames
    contact_update_interval=4.0,
)


def main(outdir: str = "output") -> None:
    Path(outdir).mkdir(exist_ok=True)
    t0 = time.time()
    out = run_organoid_pull(DEMO, seed=23, snapshots=True)
    f = out["frames"]
    print("run %.1fs | cells %d beads %d fibers %d links %d conn %.2f grip %d" % (
        time.time() - t0, out["n_cells"], out["n_beads"], out["n_fibers"],
        out["n_crosslinks"], out["connectivity"]["connected_fraction"],
        f[-1]["n_gripping_cells"]))
    print("near-field (50-75um) radial order:",
          [round(fr["shells"][2]["radial_order"], 3) for fr in f[:: max(1, len(f) // 8)]])

    snaps = out["snapshots"]
    edges = out["edges"]
    centers = out["centers"]
    links = out["crosslinks"]
    organoid_outer = float(np.max(np.linalg.norm(centers, axis=1))) + DEMO.cell_radius
    zoom = organoid_outer + 70.0

    # save positions for later re-rendering
    np.savez_compressed(f"{outdir}/g5_demo_seed23.npz",
                        snapshots=snaps, edges=edges, centers=centers,
                        crosslinks=np.asarray(links, dtype=float))

    png = viz.render(f"{outdir}/g5_demo_seed23.npz", f"{outdir}/g5_demo_beforeafter.png",
                     cell_radius=DEMO.cell_radius) if False else None  # render() expects initial/final keys

    gif = viz.animate(snaps, edges, centers, f"{outdir}/g5_demo.gif",
                      crosslinks=links, cell_radius=DEMO.cell_radius, span=zoom, fps=6)
    print("wrote", gif)


if __name__ == "__main__":
    main()
