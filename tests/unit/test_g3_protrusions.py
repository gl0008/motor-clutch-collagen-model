from dataclasses import replace

import numpy as np
import pytest

from g3.config import G3Config
from g3.fixtures import build_fixture
from g3.protrusions import geometry_scores, step_protrusion_turnover, update_traction_scores
from g3.state import ClutchState, ProtrusionState, RigidCellState


def test_aligned_geometry_is_nematic_and_rotates_by_two_sectors():
    cfg = G3Config()
    rng = np.random.default_rng(0)
    cell = RigidCellState.at_origin(cfg.cell_radius)
    state = ProtrusionState.initialize(cfg.n_sectors, 2, rng)
    aligned = build_fixture("aligned_8", cfg, rng)
    _, _, score = geometry_scores(state, cell, aligned.initial_positions, aligned, cfg)
    assert score[0] == pytest.approx(score[12])
    assert np.argmax(score[:12]) == 0
    assert np.argmax(score[12:]) + 12 == 12

    rotated = build_fixture("aligned_8_rotated_30", cfg, rng)
    _, _, rotated_score = geometry_scores(state, cell, rotated.initial_positions, rotated, cfg)
    assert rotated_score[2] == pytest.approx(rotated_score[14])
    assert np.argmax(rotated_score[:12]) == 2
    assert np.argmax(rotated_score[12:]) == 2


def test_feedback_biases_replacement_but_ablation_is_uniform():
    cfg = replace(G3Config(), n_active_protrusions=1, protrusion_lifetime=1.0e-9,
                  beta_geometry=12.0, beta_traction=0.0)

    def count_target(feedback):
        rng = np.random.default_rng(123)
        selected = 0
        for _ in range(400):
            state = ProtrusionState.initialize(cfg.n_sectors, 1, rng, prescribed=[12])
            state.geometry_score[0] = 1.0
            step_protrusion_turnover(state, cfg, rng, feedback_enabled=feedback)
            selected += int(state.active[0])
        return selected

    assert count_target(True) > 380
    assert count_target(False) < 40


def test_traction_success_saturates_with_engaged_clutches_and_load():
    cfg = replace(G3Config(), n_clutches=10, n_motors=10, n_active_protrusions=1,
                  feedback_time=0.005)
    # feedback_time == dt makes the exponential moving average equal the instant score.
    rng = np.random.default_rng(0)
    state = ProtrusionState.initialize(cfg.n_sectors, 1, rng, prescribed=[3])
    clutches = ClutchState.empty(cfg.n_clutches)
    clutches.sector_id[:] = 3
    clutches.bound[:5] = True
    stall = cfg.n_motors * cfg.motor_force
    clutches.force_vector[:5, 0] = stall / 10.0
    instant = update_traction_scores(state, clutches, cfg)
    expected = 1.0 - np.exp(-5.0 / cfg.adhesion_clutch_scale)
    assert instant[3] == pytest.approx(expected)
    assert state.traction_score[3] == pytest.approx(expected)
