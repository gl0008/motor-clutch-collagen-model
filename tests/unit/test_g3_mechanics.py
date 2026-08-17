from dataclasses import replace

import numpy as np
import pytest

from g3.config import G3Config
from g3.fixtures import build_fixture, rotated_fixture
from g3.mechanics import (
    bell_off_rate,
    clutch_geometry,
    closest_material_point,
    conservation_errors,
    force_velocity,
    project_clutch_forces,
    update_spatial_clutches,
)
from g3.state import ClutchState, ProtrusionState, RigidCellState


def _single_bound_clutch(cfg, fixture, alpha=0.37):
    clutches = ClutchState.empty(1)
    clutches.bound[0] = True
    clutches.fiber_id[0] = 0
    clutches.segment_id[0] = 0
    clutches.alpha[0] = alpha
    clutches.body_anchor_angle[0] = 0.0
    clutches.force_vector[0] = np.array([-3.0e-12, 1.2e-12])
    return clutches


def test_closest_point_returns_segment_material_coordinate():
    cfg = G3Config()
    fixture = build_fixture("single_fibre", cfg, np.random.default_rng(0))
    anchor = fixture.initial_positions[0] + np.array([0.25e-6, 0.7e-6])
    segment, fibre, alpha, point, distance = closest_material_point(
        anchor, fixture.initial_positions, fixture.network.fiber_bonds,
        fixture.segment_fiber_id)
    assert segment == 0
    assert fibre == 0
    assert alpha == pytest.approx(0.25)
    assert point[0] == pytest.approx(fixture.initial_positions[0, 0] + 0.25e-6)
    assert distance == pytest.approx(0.7e-6)


def test_material_point_tracks_same_segment_and_alpha_after_deformation():
    cfg = G3Config()
    fixture = build_fixture("single_fibre", cfg, np.random.default_rng(0))
    clutches = _single_bound_clutch(cfg, fixture, alpha=0.37)
    cell = RigidCellState.at_origin(cfg.cell_radius)
    positions = fixture.initial_positions.copy()
    positions[0] += np.array([0.2e-6, -0.3e-6])
    positions[1] += np.array([-0.1e-6, 0.5e-6])
    material, _, _, _ = clutch_geometry(clutches, cell, positions,
                                        fixture.network.fiber_bonds)
    expected = 0.63 * positions[0] + 0.37 * positions[1]
    assert np.allclose(material[0], expected)
    assert clutches.segment_id[0] == 0
    assert clutches.alpha[0] == pytest.approx(0.37)


def test_bell_rate_and_force_velocity_have_expected_limits():
    cfg = G3Config()
    rates = bell_off_rate(np.array([0.0, cfg.bell_force]), cfg)
    assert rates[0] == pytest.approx(cfg.unbind_rate)
    assert rates[1] == pytest.approx(np.e * cfg.unbind_rate)
    assert force_velocity(0.0, cfg.n_motors, cfg) == pytest.approx(cfg.unloaded_actin_speed)
    assert force_velocity(cfg.n_motors * cfg.motor_force, cfg.n_motors, cfg) == 0.0


def test_gaussian_projection_preserves_force_and_first_moment():
    cfg = G3Config()
    fixture = build_fixture("single_fibre", cfg, np.random.default_rng(0))
    clutches = _single_bound_clutch(cfg, fixture)
    cell = RigidCellState.at_origin(cfg.cell_radius)
    material, _, motor, _ = clutch_geometry(clutches, cell, fixture.initial_positions,
                                            fixture.network.fiber_bonds)
    connector = motor[0] - material[0]
    clutches.force_vector[0] = 3.0e-12 * connector / np.linalg.norm(connector)
    bead_forces = project_clutch_forces(
        clutches, material, fixture.initial_positions, fixture, cfg)
    assert np.allclose(bead_forces.sum(axis=0), clutches.force_vector.sum(axis=0), atol=1e-25)
    point_moment = material[0, 0] * clutches.force_vector[0, 1] - material[0, 1] * clutches.force_vector[0, 0]
    projected_moment = np.sum(
        fixture.initial_positions[:, 0] * bead_forces[:, 1]
        - fixture.initial_positions[:, 1] * bead_forces[:, 0])
    assert projected_moment == pytest.approx(point_moment, abs=1e-25)
    force_error, torque_error = conservation_errors(
        clutches, material, motor, bead_forces, fixture.initial_positions, cell)
    assert force_error < 1e-10
    assert torque_error < 1e-8


def test_projection_rotates_covariantly():
    cfg = G3Config()
    angle = np.radians(30.0)
    c, s = np.cos(angle), np.sin(angle)
    rotation = np.array([[c, -s], [s, c]])
    fixture = build_fixture("single_fibre", cfg, np.random.default_rng(0))
    rotated = rotated_fixture(fixture, angle)
    clutch = _single_bound_clutch(cfg, fixture)
    rotated_clutch = _single_bound_clutch(cfg, rotated)
    rotated_clutch.force_vector[:] = clutch.force_vector @ rotation.T
    material, _, _, _ = clutch_geometry(clutch, RigidCellState.at_origin(cfg.cell_radius),
                                        fixture.initial_positions, fixture.network.fiber_bonds)
    rotated_cell = RigidCellState.at_origin(cfg.cell_radius)
    rotated_cell.body_angle = angle
    rotated_clutch.body_anchor_angle[:] = clutch.body_anchor_angle
    rotated_material, _, _, _ = clutch_geometry(
        rotated_clutch, rotated_cell, rotated.initial_positions, rotated.network.fiber_bonds)
    force = project_clutch_forces(clutch, material, fixture.initial_positions, fixture, cfg)
    force_rotated = project_clutch_forces(
        rotated_clutch, rotated_material, rotated.initial_positions, rotated, cfg)
    assert np.allclose(force_rotated, force @ rotation.T, atol=1e-24)


def test_binding_initializes_zero_load_and_empty_fixture_cannot_bind():
    cfg = replace(G3Config(), bind_rate=1.0e9, n_clutches=20, n_motors=20)
    rng = np.random.default_rng(4)
    cell = RigidCellState.at_origin(cfg.cell_radius)
    protrusions = ProtrusionState.initialize(cfg.n_sectors, 1, rng, prescribed=[0])
    fixture = build_fixture("single_fibre", cfg, rng)
    clutches = ClutchState.empty(cfg.n_clutches)
    update_spatial_clutches(clutches, protrusions, cell, fixture.initial_positions,
                             fixture, cfg, 0.0, rng)
    assert clutches.bound.any()
    assert np.allclose(clutches.force_vector[clutches.bound], 0.0)

    empty = build_fixture("empty", cfg, rng)
    none = ClutchState.empty(cfg.n_clutches)
    update_spatial_clutches(none, protrusions, cell, empty.initial_positions,
                             empty, cfg, 0.0, rng)
    assert not none.bound.any()
