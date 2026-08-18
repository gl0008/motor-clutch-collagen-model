"""Regenerate the G3ABC-R G2-scale attachment-to-alignment demonstrations.

This is a controlled loading study, not a calibrated tumour prediction.  G3A/B/C
use a 5 nN total traction target only after explicit tip attachment.  G3S is a
coarse eight-patch, traction-only spheroid shell at 20 nN (2.5 nN per nominal
patch); it contains neither cell--cell mechanics nor an interior-cell model.
"""

from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import shutil

import numpy as np

from g3.config import G3Config
from g3.run import _json_ready, save_run
from g3.simulation import run_g3
from g3.visualization import make_comparison_animation, make_stage_animation


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "g3_g2scale"
ASSETS = ROOT / "docs" / "assets" / "g3_g2scale"
DURATION = 20.0


def _snapshot_metrics(result):
    initial = result.initial_positions
    maximum_force = 0.0
    maximum_displacement = 0.0
    for frame in result.snapshots:
        maximum_force = max(
            maximum_force,
            float(np.linalg.norm(frame.clutch_forces, axis=1).sum() * 1.0e9),
        )
        maximum_displacement = max(
            maximum_displacement,
            float(np.linalg.norm(frame.positions - initial, axis=1).max(initial=0.0) * 1.0e6),
        )
    first_attachment = next(
        (float(frame.time) for frame in result.snapshots if frame.bound_points.size),
        None,
    )
    return {
        "first_attachment_s_at_frame_resolution": first_attachment,
        "max_applied_clutch_load_nN": maximum_force,
        "max_bead_displacement_um_across_saved_frames": maximum_displacement,
    }


def _save(result, name: str):
    output = RESULTS / name
    save_run(result, output, make_gif=False)
    metrics = _snapshot_metrics(result)
    with (output / "g2scale_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(_json_ready(metrics), handle, indent=2)
    return output, metrics


def main():
    base = G3Config.from_yaml(ROOT / "src" / "config" / "params_g3.yaml")
    # G2's collagen mechanics and its total 5 nN loading magnitude, but gated
    # by a genuine G3 protrusion-tip attachment before the three-second ramp.
    cell_cfg = replace(
        base,
        traction_target=5.0e-9,
        traction_ramp_time=3.0,
        duration_g3a=DURATION,
        duration_g3b=DURATION,
        duration_g3c=DURATION,
        frame_interval=0.5,
        metrics_interval=0.1,
    )

    g3a = run_g3("g3a", "scaled_isotropic_99", cell_cfg, seed=3, duration=DURATION)
    g3b_off = run_g3("g3b", "scaled_aligned_99", cell_cfg, seed=4, duration=DURATION,
                      feedback_enabled=False)
    g3b_on = run_g3("g3b", "scaled_aligned_99", cell_cfg, seed=4, duration=DURATION,
                     feedback_enabled=True)
    g3c_fixed = run_g3("g3b", "scaled_aligned_99", cell_cfg, seed=6, duration=DURATION,
                       feedback_enabled=True)
    g3c_released = run_g3("g3c", "scaled_aligned_99", cell_cfg, seed=6, duration=DURATION,
                          feedback_enabled=True)

    # A deliberately coarse surface-cell shell: 8 ECM-facing patches, each
    # represented by a 200-clutch module.  Twenty nN is a visual/mechanical
    # sensitivity condition, not a fitted spheroid force measurement.
    shell_cfg = replace(
        cell_cfg,
        cell_radius=base.spheroid_radius,
        n_clutches=base.n_clutches * base.spheroid_surface_patches,
        n_motors=base.n_motors * base.spheroid_surface_patches,
        traction_target=20.0e-9,
        traction_ramp_time=4.0,
        duration_g3s=DURATION,
    )
    g3s = run_g3("g3s", "scaled_isotropic_99", shell_cfg, seed=8, duration=DURATION)

    runs = {
        "g3a_5nN": g3a,
        "g3b_feedback_off_5nN": g3b_off,
        "g3b_feedback_on_5nN": g3b_on,
        "g3c_fixed_5nN": g3c_fixed,
        "g3c_released_5nN": g3c_released,
        "g3s_surface_shell_20nN": g3s,
    }
    record = {name: _save(result, name)[1] for name, result in runs.items()}

    RESULTS.mkdir(parents=True, exist_ok=True)
    ASSETS.mkdir(parents=True, exist_ok=True)
    gifs = {
        "g3a_g2scale_loading.gif": lambda path: make_stage_animation(g3a, path, fps=8, max_frames=18),
        "g3b_g2scale_feedback.gif": lambda path: make_comparison_animation(
            [g3b_off, g3b_on], path, fps=8, max_frames=18,
            panel_labels=["Intrinsic polarity only", "+ adhesion feedback"],
            title="G3B diagnostic · no OFF/ON separation in this short matched run",
        ),
        "g3c_g2scale_released.gif": lambda path: make_comparison_animation(
            [g3c_fixed, g3c_released], path, fps=8, max_frames=18,
            panel_labels=["Rigid cell (fixed)", "Released rigid cell"],
            title="G3C question · does releasing the cell change reaction-driven motion?",
        ),
        "g3s_surface_shell_20nN.gif": lambda path: make_stage_animation(g3s, path, fps=8, max_frames=18),
    }
    for filename, render in gifs.items():
        local = RESULTS / filename
        render(local)
        shutil.copy2(local, ASSETS / filename)

    manifest = {
        "purpose": "G2-scale controlled-loading G3ABC-R visual demonstrations",
        "duration_s": DURATION,
        "g3abc_total_traction_target_nN": 5.0,
        "traction_ramp_after_first_attachment_s": 3.0,
        "spheroid": {
            "kind": "traction-only coarse surface shell; no interior or cell-cell mechanics",
            "surface_patches": 8,
            "nominal_clutches_per_patch": 200,
            "total_traction_target_nN": 20.0,
            "nominal_target_per_patch_nN": 2.5,
            "radius_um": 35.0,
        },
        "runs": record,
    }
    with (RESULTS / "experiment_manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
