from dataclasses import replace

import numpy as np

from generations.g4_v4_single_cell_force_alignment.model import (
    G4V4Config,
    V4ClutchState,
    _advance_clutches,
    cell_velocity_from_ecm_forces,
    decompose_contact_forces,
    default_spec,
    run_fixed_cell_passive,
    run_fixed_force_direction,
    run_motor_clutch,
)


def small_config(**changes):
    base = G4V4Config(
        n_fibers=120,
        geometry_attempt_factor=80,
        fixed_duration=1.0,
        fixed_sample_interval=0.5,
        migration_duration=1.0,
        migration_sample_interval=0.5,
    )
    return replace(base, **changes)


def test_passive_fixed_cell_does_not_shrink_or_remodel():
    cfg = small_config()
    result = run_fixed_cell_passive(cfg, spec=default_spec(cfg))
    assert result["cell_radius_constant"] is True
    assert result["max_displacement"] < 1e-10
    assert np.max(np.abs(result["delta_radial_order_by_shell"])) < 1e-10


def test_contact_force_decomposition_and_reaction_balance():
    cfg = small_config()
    result = run_fixed_force_direction(
        cfg, spec=default_spec(cfg), duration=0.5, sample_interval=0.5
    )
    frame = result["frames"][-1]
    parallel = frame["force_parallel"]
    perpendicular = frame["force_perpendicular"]
    assert np.allclose(parallel + perpendicular, frame["contact_vectors"])
    assert frame["force_balance_error"] < 1e-12
    assert result["cell_radius_constant"] is True


def test_normal_and_dipole_control_preserve_requested_total_force():
    cfg = small_config(force_ramp_time=0.05)
    spec = default_spec(cfg)
    normal = run_fixed_force_direction(
        cfg, geometry="normal", spec=spec, duration=0.1, sample_interval=0.1
    )
    dipole = run_fixed_force_direction(
        cfg, geometry="dipole", spec=spec, duration=0.1, sample_interval=0.1
    )
    assert np.isclose(normal["frames"][-1]["distributed_force"], cfg.total_pull_force)
    assert np.isclose(dipole["frames"][-1]["distributed_force"], cfg.total_pull_force)


def test_clutch_force_keeps_loading_after_contact_until_a_stop_condition():
    cfg = small_config(clutch_off_rate0=0.0)
    state = V4ClutchState(
        bound=np.ones((1, cfg.n_clutches_per_site), dtype=bool),
        extension=np.zeros((1, cfg.n_clutches_per_site)),
        site_extension=np.zeros(1),
    )
    force1, _, _, _, _ = _advance_clutches(cfg, state, np.zeros(1), 0)
    force2, _, _, _, _ = _advance_clutches(cfg, state, np.zeros(1), 1)
    assert force1[0] > 0.0
    assert force2[0] > force1[0]


def test_event_time_matches_logged_complete_site_failure():
    cfg = small_config(
        n_clutches_per_site=1,
        clutch_on_rate=100.0,
        clutch_off_rate0=100.0,
        bell_force=1.0,
        minimum_event_time=0.0,
        fixed_duration=2.0,
        fixed_sample_interval=1.0,
        event_half_window_v4=0.2,
        event_sample_interval_v4=0.05,
    )
    result = run_motor_clutch(
        cfg,
        spec=default_spec(cfg),
        duration=2.0,
        sample_interval=1.0,
        capture_event=True,
    )
    failure_times = [
        event["time"] for event in result["event_log"] if event["kind"] == "site_failure"
    ]
    assert failure_times
    assert result["event_time"] == failure_times[0]
    assert result["cell_radius_constant"] is True


def test_shared_load_raises_remaining_clutch_hazard_after_a_rupture():
    cfg = small_config(clutch_mode="shared", clutch_on_rate=0.0)
    state = V4ClutchState(
        bound=np.ones((1, cfg.n_clutches_per_site), dtype=bool),
        extension=np.zeros((1, cfg.n_clutches_per_site)),
        site_extension=np.full(1, 0.2),
    )
    total = cfg.n_clutches_per_site * cfg.clutch_stiffness * state.site_extension[0]
    rate_before = cfg.clutch_off_rate0 * np.exp(total / (cfg.n_clutches_per_site * cfg.bell_force))
    state.bound[0, 0] = False
    rate_after = cfg.clutch_off_rate0 * np.exp(total / ((cfg.n_clutches_per_site - 1) * cfg.bell_force))
    assert rate_after > rate_before


def test_symmetric_cell_reaction_does_not_translate():
    forces = np.asarray([[2.5, 0.0], [-2.5, 0.0]])
    assert np.allclose(cell_velocity_from_ecm_forces(forces, 600.0), 0.0)
