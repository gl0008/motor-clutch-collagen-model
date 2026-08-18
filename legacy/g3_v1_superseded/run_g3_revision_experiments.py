"""Run the version-preserving G3-R multi-fibre validation/visualization suite."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from g3.config import G3Config
from g3.run import _json_ready, load_saved_run, save_run
from g3.simulation import run_g3
from g3.validation import summarize_ensemble
from g3.visualization import make_comparison_animation, make_stage_animation


def _first_positive_time(run, key):
    index = np.flatnonzero(run.traces[key] > 0)
    return None if index.size == 0 else float(run.traces["time"][index[0]])


def _axis_angle(run):
    x = float(run.traces["protrusion_x"][-1])
    y = float(run.traces["protrusion_y"][-1])
    return float(np.arctan2(y, x))


def _run_representative(job):
    """Process-pool entry point for layout-only GIF regeneration."""
    name, stage, fixture, cfg, seed, duration, feedback_enabled = job
    return name, run_g3(
        stage, fixture, cfg, seed=seed, duration=duration,
        feedback_enabled=feedback_enabled,
    )


def _render_representative_gifs(runs, figures):
    """Render the three public GIFs from already available run objects."""
    make_stage_animation(runs["g3a"], figures / "g3a_multifibre_protrusion.gif", fps=8)
    make_comparison_animation(
        [runs["g3b_off"], runs["g3b_on"]],
        figures / "g3b_intrinsic_polarity_feedback.gif",
        fps=8,
        panel_labels=["Adhesion feedback OFF", "Adhesion feedback ON"],
        title="G3B-R · does adhesion feedback stabilize a different protrusion?",
    )
    make_comparison_animation(
        [runs["fixed"], runs["released"]],
        figures / "g3c_fixed_released.gif",
        fps=8,
        panel_labels=["Fixed cell", "Released cell"],
        title="G3C-R · same polarity/ECM loop, reaction-driven rigid-cell release",
    )


def _rerender_representative_gifs(cfg, figures, workers):
    jobs = [
        ("g3a", "g3a", "scaled_isotropic_99", cfg, 3, cfg.duration_g3a, True),
        ("g3b_on", "g3b", "scaled_aligned_99", cfg, 4, cfg.duration_g3b, True),
        ("g3b_off", "g3b", "scaled_aligned_99", cfg, 4, cfg.duration_g3b, False),
        ("fixed", "g3b", "scaled_aligned_99", cfg, 6, cfg.duration_g3c, True),
        ("released", "g3c", "scaled_aligned_99", cfg, 6, cfg.duration_g3c, True),
    ]
    with ProcessPoolExecutor(max_workers=min(max(workers, 1), len(jobs))) as pool:
        runs = dict(pool.map(_run_representative, jobs))
    _render_representative_gifs(runs, figures)
    return runs


def main(argv=None):
    parser = argparse.ArgumentParser(description="G3-R 99-fibre experiments")
    parser.add_argument("--config", default="src/config/params_g3_scaled.yaml")
    parser.add_argument("--output", default="results/g3_revision")
    parser.add_argument("--figures", default="figures/g3_revision")
    parser.add_argument("--ensemble-seeds", type=int, default=4)
    parser.add_argument("--ensemble-duration", type=float, default=15.0)
    parser.add_argument(
        "--representative-gifs-only", action="store_true",
        help="rerun only the five displayed trajectories and replace the three GIFs",
    )
    parser.add_argument(
        "--render-saved-frames", action="store_true",
        help="replace the three GIFs from saved numerical frames without rerunning physics",
    )
    parser.add_argument("--workers", type=int, default=5)
    args = parser.parse_args(argv)

    cfg = G3Config.from_yaml(args.config)
    output = Path(args.output)
    figures = Path(args.figures)
    output.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    if args.render_saved_frames:
        runs = {
            "g3a": load_saved_run(output / "g3a_multifibre"),
            "g3b_on": load_saved_run(output / "g3b_feedback_on"),
            "g3b_off": load_saved_run(output / "g3b_feedback_off"),
            "fixed": load_saved_run(output / "g3c_fixed_control"),
            "released": load_saved_run(output / "g3c_released"),
        }
        _render_representative_gifs(runs, figures)
        print(json.dumps({
            "rendered_from_saved_frames": True,
            "statuses": {name: run.status for name, run in runs.items()},
        }, indent=2))
        return 0

    if args.representative_gifs_only:
        runs = _rerender_representative_gifs(cfg, figures, args.workers)
        print(json.dumps({
            "rendered": [
                "g3a_multifibre_protrusion.gif",
                "g3b_intrinsic_polarity_feedback.gif",
                "g3c_fixed_released.gif",
            ],
            "statuses": {name: run.status for name, run in runs.items()},
            "layout": "fixed axes; no per-frame tight_layout",
        }, indent=2))
        return 0 if all(run.status == "complete" for run in runs.values()) else 2

    # Run the complete fresh suite in a bounded process pool. Each job still
    # executes an independent deterministic trajectory; parallelism changes
    # only wall time, not its seed stream or physics.
    seeds = list(range(args.ensemble_seeds))
    jobs = [
        ("g3a", "g3a", "scaled_isotropic_99", cfg, 3, cfg.duration_g3a, True),
        ("g3a_motor_off", "g3a", "scaled_isotropic_99",
         replace(cfg, motor_force=0.0), 3, cfg.duration_g3a, True),
        ("g3b_on", "g3b", "scaled_aligned_99", cfg, 4, cfg.duration_g3b, True),
        ("g3b_off", "g3b", "scaled_aligned_99", cfg, 4, cfg.duration_g3b, False),
        ("fixed", "g3b", "scaled_aligned_99", cfg, 6, cfg.duration_g3c, True),
        ("released", "g3c", "scaled_aligned_99", cfg, 6, cfg.duration_g3c, True),
        ("empty", "g3c", "empty", cfg, 6, 1.0, True),
    ]
    jobs.extend(
        (f"aligned_{seed}", "g3b", "scaled_aligned_99", cfg,
         seed, args.ensemble_duration, True)
        for seed in seeds
    )
    jobs.extend(
        (f"isotropic_{seed}", "g3b", "scaled_isotropic_99", cfg,
         seed, args.ensemble_duration, True)
        for seed in seeds
    )
    if args.workers == 1:
        runs = dict(map(_run_representative, jobs))
    else:
        with ProcessPoolExecutor(max_workers=min(max(args.workers, 1), len(jobs))) as pool:
            runs = dict(pool.map(_run_representative, jobs))

    g3a = runs["g3a"]
    g3a_motor_off = runs["g3a_motor_off"]
    g3b_on = runs["g3b_on"]
    g3b_off = runs["g3b_off"]
    fixed = runs["fixed"]
    released = runs["released"]
    empty = runs["empty"]
    aligned = [runs[f"aligned_{seed}"] for seed in seeds]
    isotropic = [runs[f"isotropic_{seed}"] for seed in seeds]

    # G3A-R: same multi-fibre network, attachment mechanics on versus motor off.
    save_run(g3a, output / "g3a_multifibre", make_gif=False)
    save_run(g3a_motor_off, output / "g3a_motor_off", make_gif=False)
    make_stage_animation(g3a, figures / "g3a_multifibre_protrusion.gif", fps=8)

    # G3B-R: same seed/ECM with intracellular adhesion feedback on versus off.
    save_run(g3b_on, output / "g3b_feedback_on", make_gif=False)
    save_run(g3b_off, output / "g3b_feedback_off", make_gif=False)
    make_comparison_animation(
        [g3b_off, g3b_on],
        figures / "g3b_intrinsic_polarity_feedback.gif",
        fps=8,
        panel_labels=["Adhesion feedback OFF", "Adhesion feedback ON"],
        title="G3B-R · does adhesion feedback stabilize a different protrusion?",
    )

    # Small held-out mechanism ensemble. This is a fresh revised-model calibration
    # suite, not the old sealed 100-seed final validation campaign.
    aligned_metrics = summarize_ensemble(aligned, 0.0).to_dict()
    isotropic_metrics = summarize_ensemble(isotropic, 0.0).to_dict()

    # G3C-R: fixed versus released on the same 99-fibre aligned fixture/seed.
    save_run(fixed, output / "g3c_fixed_control", make_gif=False)
    save_run(released, output / "g3c_released", make_gif=False)
    make_comparison_animation(
        [fixed, released],
        figures / "g3c_fixed_released.gif",
        fps=8,
        panel_labels=["Fixed cell", "Released cell"],
        title="G3C-R · same polarity/ECM loop, reaction-driven rigid-cell release",
    )

    report = {
        "scope": "fresh G3-R calibration/experiment suite; not realistic 3D validation",
        "config": cfg.to_dict(),
        "network": {
            "n_fibres": g3a.summary["n_fibres"],
            "n_beads": g3a.summary["n_beads"],
            "n_crosslinks": len(g3a.fixture.network.external_network.crosslinks),
            "connectivity": g3a.fixture.connectivity_report,
        },
        "g3a": {
            "first_attachment_time_s": _first_positive_time(g3a, "bound_count"),
            "max_bound_clutches": g3a.summary["max_bound_clutches"],
            "max_protrusion_length_um": float(
                g3a.traces["max_protrusion_length"].max(initial=0.0) * 1.0e6
            ),
            "max_bead_displacement_um": g3a.summary["bead_displacement"]["max"] * 1.0e6,
            "motor_off_max_bead_displacement_um": (
                g3a_motor_off.summary["bead_displacement"]["max"] * 1.0e6
            ),
        },
        "g3b": {
            "feedback_on_final_axis_deg": float(np.degrees(_axis_angle(g3b_on))),
            "feedback_off_final_axis_deg": float(np.degrees(_axis_angle(g3b_off))),
            "feedback_on_peak_adhesion_signal": float(max(
                np.max(snapshot.traction_scores) for snapshot in g3b_on.snapshots
            )),
            "aligned_ensemble": aligned_metrics,
            "isotropic_ensemble": isotropic_metrics,
        },
        "g3c": {
            "fixed_displacement_um": fixed.summary["cell_net_displacement_m"] * 1.0e6,
            "released_displacement_um": released.summary["cell_net_displacement_m"] * 1.0e6,
            "released_rotation_deg": float(np.degrees(released.cell.body_angle)),
            "empty_displacement_um": empty.summary["cell_net_displacement_m"] * 1.0e6,
        },
    }
    report["gates"] = {
        "all_primary_runs_complete": all(
            run.status == "complete"
            for run in (g3a, g3a_motor_off, g3b_on, g3b_off, fixed, released, empty)
        ),
        "scaled_99_fibres": report["network"]["n_fibres"] == 99,
        "boundary_connectivity": (
            report["network"]["connectivity"]["connected_fraction"]
            >= cfg.scaled_required_connected_fraction
        ),
        "g3a_grow_before_attach": (
            report["g3a"]["first_attachment_time_s"] is not None
            and report["g3a"]["first_attachment_time_s"] > 0.0
        ),
        "g3a_motor_deforms_more_than_motor_off": (
            report["g3a"]["max_bead_displacement_um"]
            > report["g3a"]["motor_off_max_bead_displacement_um"] + 1.0e-12
        ),
        "g3b_aligned_nematic_positive": aligned_metrics["nematic_order"] > 0.0,
        "g3b_isotropic_all_valid": isotropic_metrics["n_valid"] == isotropic_metrics["n_runs"],
        "g3c_no_hidden_motion": report["g3c"]["empty_displacement_um"] == 0.0,
        "g3c_reaction_motion_nonzero": report["g3c"]["released_displacement_um"] > 0.0,
    }
    with open(output / "validation_summary.json", "w", encoding="utf-8") as handle:
        json.dump(_json_ready(report), handle, indent=2, ensure_ascii=False)
    with open(figures / "metrics.json", "w", encoding="utf-8") as handle:
        json.dump(_json_ready(report), handle, indent=2, ensure_ascii=False)
    print(json.dumps(_json_ready(report), indent=2, ensure_ascii=False))
    return 0 if all(report["gates"].values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
