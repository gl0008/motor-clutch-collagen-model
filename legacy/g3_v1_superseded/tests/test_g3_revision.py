import numpy as np
import pytest

from g3.config import G3Config
from g3.fixtures import build_fixture
from g3.mechanics import apply_traction_protocol, clutch_geometry
from g3.protrusions import step_intrinsic_polarity
from g3.state import ClutchState, ProtrusionState, RigidCellState


def test_cell_intrinsic_polarity_breaks_symmetry_and_conserves_activity_pool():
    cfg = G3Config()
    rng = np.random.default_rng(14)
    state = ProtrusionState.initialize(cfg.n_sectors, cfg.n_active_protrusions, rng)
    for _ in range(200):
        step_intrinsic_polarity(state, cfg, rng, dt=0.1, feedback_enabled=False)
    assert state.activity.sum() == pytest.approx(cfg.polarity_total_activity)
    assert state.activity.max() > 0.60
    assert state.active.sum() == cfg.n_active_protrusions


def test_expressed_protrusions_use_lifetime_and_hysteresis_instead_of_teleporting():
    cfg = G3Config()
    rng = np.random.default_rng(9)
    state = ProtrusionState.initialize(cfg.n_sectors, 2, rng, prescribed=[0, 12])
    state.activity[:] = 0.0
    state.activity[[0, 12]] = [0.45, 0.40]
    state.activity[5] = 0.15
    # A challenger cannot replace a still-growing protrusion even if its
    # activity briefly becomes the largest value.
    state.activity[5] = 0.80
    state.activity[12] = 0.05
    step_intrinsic_polarity(state, cfg, rng, dt=0.1, feedback_enabled=False)
    assert set(np.flatnonzero(state.active)) == {0, 12}

    # Once mature, replacement still needs a declared activity margin.
    state.active_age[[0, 12]] = cfg.protrusion_min_lifetime
    state.activity[:] = 0.0
    state.activity[0] = 0.50
    state.activity[12] = 0.10
    state.activity[5] = 0.50
    state.traction_score[:] = 0.0
    step_intrinsic_polarity(state, cfg, rng, dt=0.1, feedback_enabled=False)
    assert 5 in np.flatnonzero(state.active)
    assert 0 in np.flatnonzero(state.active)


def test_clutch_motor_endpoint_starts_at_explicit_protrusion_tip():
    cfg = G3Config()
    rng = np.random.default_rng(2)
    protrusions = ProtrusionState.initialize(cfg.n_sectors, 1, rng, prescribed=[0])
    protrusions.length[0] = 2.5e-6
    clutches = ClutchState.empty(1)
    clutches.sector_id[0] = 0
    clutches.body_anchor_angle[0] = 0.0
    cell = RigidCellState.at_origin(cfg.cell_radius)
    fixture = build_fixture("single_fibre", cfg, rng)
    _, tip, motor, _ = clutch_geometry(
        clutches, cell, fixture.initial_positions, fixture.network.fiber_bonds, protrusions
    )
    expected = np.array([cfg.cell_radius + protrusions.length[0], 0.0])
    assert tip[0] == pytest.approx(expected)
    assert motor[0] == pytest.approx(expected)


def test_scaled_fixture_has_99_fibres_crosslinks_and_boundary_percolation():
    cfg = G3Config()
    fixture = build_fixture("scaled_isotropic_99", cfg, np.random.default_rng(3))
    assert fixture.network.n_fibers == 99
    assert fixture.network.n_beads > 6000
    assert len(fixture.network.external_network.crosslinks) > 100
    assert fixture.connectivity_report["contact_fibers_connected"]
    assert (
        fixture.connectivity_report["connected_fraction"]
        >= cfg.scaled_required_connected_fraction
    )


def test_g2_scale_traction_ramp_persists_across_clutch_turnover():
    cfg = G3Config(traction_target=5.0e-9, traction_ramp_time=2.0)
    clutches = ClutchState.empty(2)
    clutches.bound[0] = True
    clutches.binding_time[0] = 1.0
    clutches.force_vector[0] = [2.0e-12, 0.0]
    assert apply_traction_protocol(clutches, cfg, 2.0) == pytest.approx(2.5e-9)
    assert np.linalg.norm(clutches.force_vector[0]) == pytest.approx(2.5e-9)

    # The first clutch ruptures; the protocol must not restart from the new
    # clutch's binding time.
    clutches.detach(np.array([0]))
    clutches.bound[1] = True
    clutches.binding_time[1] = 3.0
    clutches.force_vector[1] = [3.0e-12, 0.0]
    assert apply_traction_protocol(clutches, cfg, 3.5) == pytest.approx(5.0e-9)
    assert np.linalg.norm(clutches.force_vector[1]) == pytest.approx(5.0e-9)
