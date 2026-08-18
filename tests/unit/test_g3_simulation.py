from dataclasses import replace

import numpy as np

from g3.config import G3Config
from g3.simulation import run_g3


def _fast_config():
    return replace(
        G3Config(),
        n_clutches=20,
        n_motors=20,
        bind_rate=100.0,
        unbind_rate=1.0e-9,
        metrics_interval=0.01,
        frame_interval=0.01,
        geometry_update_interval=0.01,
    )


def test_g3a_short_run_binds_and_conserves_internal_force_and_torque():
    result = run_g3("g3a", "single_fibre", _fast_config(), seed=7, duration=0.05)
    assert result.status == "complete"
    assert result.summary["max_bound_clutches"] > 0
    assert result.summary["max_force_error"] < 1.0e-10
    assert result.summary["max_torque_error"] < 1.0e-8
    assert np.allclose(result.traces["cell_x"], 0.0)
    assert np.allclose(result.traces["cell_y"], 0.0)


def test_g3_seeded_short_run_is_replayable():
    cfg = _fast_config()
    first = run_g3("g3a", "single_fibre", cfg, seed=11, duration=0.03)
    second = run_g3("g3a", "single_fibre", cfg, seed=11, duration=0.03)
    assert np.array_equal(first.final_positions, second.final_positions)
    assert np.array_equal(first.traces["bound_count"], second.traces["bound_count"])


def test_g3b_keeps_cell_fixed_while_g3c_releases_reaction_motion():
    cfg = replace(_fast_config(), protrusion_lifetime=0.1)
    fixed = run_g3("g3b", "asymmetric_torque", cfg, seed=2, duration=0.5)
    mobile = run_g3("g3c", "asymmetric_torque", cfg, seed=2, duration=0.5)
    assert np.allclose(fixed.traces["cell_x"], 0.0)
    assert np.allclose(fixed.traces["cell_angle"], 0.0)
    assert mobile.summary["cell_net_displacement_m"] > 0.0
    assert abs(mobile.summary["cell_final_angle_rad"]) > 0.0


def test_g3c_empty_fixture_has_no_hidden_self_propulsion_or_rotation():
    result = run_g3("g3c", "empty", _fast_config(), seed=9, duration=0.05)
    assert result.summary["cell_net_displacement_m"] == 0.0
    assert result.summary["cell_path_length_m"] == 0.0
    assert result.summary["cell_final_angle_rad"] == 0.0
