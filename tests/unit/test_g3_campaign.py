from dataclasses import replace

from g3.campaign import (
    ALL_CONDITIONS,
    G3B_CONDITIONS,
    G3C_CONDITIONS,
    campaign_fingerprint,
    evaluate_gates,
    selected_conditions,
    summarize_condition,
)
from g3.config import G3Config


def test_campaign_conditions_are_unique_and_cover_preregistered_controls():
    keys = [(condition.stage, condition.name) for condition in ALL_CONDITIONS]
    assert len(keys) == len(set(keys))
    assert {condition.name for condition in G3B_CONDITIONS} >= {
        "balanced", "isotropic", "aligned", "aligned_rotated_30",
        "aligned_feedback_off", "no_fibre",
    }
    assert {condition.name for condition in G3C_CONDITIONS} >= {
        "no_fibre", "isotropic", "aligned", "aligned_rotated_30",
        "asymmetric_torque", "asymmetric_torque_mirror_x",
        "aligned_drag_150", "aligned_drag_600",
    }


def test_campaign_fingerprint_changes_if_locked_physics_changes():
    cfg = G3Config()
    baseline = campaign_fingerprint(cfg, 600.0)
    assert campaign_fingerprint(replace(cfg, beta_geometry=4.0), 600.0) != baseline
    assert campaign_fingerprint(cfg, 300.0) != baseline


def test_condition_selection_respects_stage_and_name():
    selected = selected_conditions({"g3b"}, {"aligned"})
    assert [(condition.stage, condition.name) for condition in selected] == [("g3b", "aligned")]


def test_summary_uses_independent_run_vectors_and_reports_invalid_overlap():
    records = [
        {
            "status": "complete", "axis_cos2": 1.0, "axis_sin2": 0.0,
            "polar_x": 1.0, "polar_y": 0.0, "direction_rad": 0.0,
            "net_displacement_m": 1.0, "path_length_m": 2.0, "duration_s": 60.0,
            "final_cell_angle_rad": 0.1, "cell_x_m": 1.0, "cell_y_m": 0.0,
            "max_bound_clutches": 2, "max_cell_force_N": 3.0,
            "max_abs_cell_torque_N_m": 4.0, "max_force_error": 1e-15,
            "max_torque_error": 1e-15, "wall_time_s": 1.0,
        },
        {"status": "invalid_geometry_overlap"},
    ]
    summary = summarize_condition(records, director=0.0)
    assert summary["n_valid"] == 1
    assert summary["n_invalid_overlap"] == 1
    assert summary["mean_nematic_alignment"] == 1.0
    assert summary["positive_fraction"] == 1.0


def test_gate_evaluator_handles_partial_campaign_without_inventing_missing_gates():
    summary = {
        "g3b/aligned": {
            "n_records": 20, "n_valid": 20, "n_invalid_overlap": 0,
            "n_worker_error": 0, "nematic_alignment_ci95": [0.2, 0.8],
            "positive_fraction": 0.5, "estimated_axis_rad": 0.0,
        }
    }
    gates = evaluate_gates(summary)
    assert gates["g3b_aligned_guidance_ci_above_zero"]
    assert gates["g3b_aligned_plus_minus_40_60"]
    assert gates["g3b/aligned_all_runs_valid"]
    assert "g3b_feedback_reduces_guidance_50pct" not in gates
