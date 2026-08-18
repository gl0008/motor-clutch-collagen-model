"""Single-question G3 collective spheroid collagen-remodelling experiment."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import shutil

from g3.config import G3Config
from g3.run import save_run
from g3.simulation import run_g3
from g3.visualization import make_spheroid_remodelling_animation


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "g3_spheroid_remodelling"
ASSET = ROOT / "docs" / "assets" / "g3_spheroid_remodelling.gif"


def main():
    base = G3Config.from_yaml(ROOT / "src" / "config" / "params_g3.yaml")
    # 1 uN is a collective whole-spheroid sensitivity point, within the
    # micro-Newton-scale contractility reported for tumour spheroids.  It is
    # deliberately not a calibration of a 12-patch shell to a real spheroid.
    cfg = replace(
        base,
        cell_radius=base.spheroid_radius,
        n_clutches=base.n_clutches * base.spheroid_surface_patches,
        n_motors=base.n_motors * base.spheroid_surface_patches,
        traction_target=1.0e-6,
        traction_ramp_time=5.0,
        traction_sectorwise=True,
        duration_g3s=30.0,
        frame_interval=1.0,
        metrics_interval=0.1,
    )
    result = run_g3("g3s", "scaled_spheroid_99", cfg, seed=18, duration=30.0)
    save_run(result, OUT, make_gif=False)
    gif = OUT / "g3_spheroid_remodelling.gif"
    make_spheroid_remodelling_animation(result, gif, max_frames=22)
    ASSET.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(gif, ASSET)
    print(result.summary)


if __name__ == "__main__":
    main()
