from dataclasses import replace

import numpy as np
import pytest

from g3.config import G3Config
from g3.fixtures import build_fixture, rotated_fixture
from g3.mechanics import (
    bell_off_rate,
    cell_contact_forces,
    cell_contact_forces_numpy,
    clutch_geometry,
    clutch_geometry_numpy,
    closest_material_point,
    conservation_errors,
    force_velocity,
    project_clutch_forces,
    project_clutch_forces_numpy,
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
    direction = (fixture.initial_positions[1] - fixture.initial_positions[0]) / cfg.bead_spacing
    normal = np.array([-direction[1], direction[0]])
    anchor = fixture.initial_positions[0] + 0.25e-6 * direction + 0.7e-6 * normal
    segment, fibre, alpha, point, distance = closest_material_point(
        anchor, fixture.initial_positions, fixture.network.fiber_bonds,
        fixture.segment_fiber_id)
    assert segment == 0
    assert fibre == 0
    assert alpha == pytest.approx(0.25)
    assert np.allclose(point, fixture.initial_positions[0] + 0.25e-6 * direction)
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


def test_conservative_cell_contact_is_outward_and_equal_opposite():
    cfg = G3Config()
    cell = RigidCellState.at_origin(cfg.cell_radius)
    penetration = 0.2e-6
    positions = np.array([
        [cfg.cell_radius - penetration, 0.0],
        [0.0, cfg.cell_radius + 0.3e-6],
    ])
    bead, reaction, torque, maximum, count, energy = cell_contact_forces_numpy(
        positions, cell, cfg)
    expected_force = cfg.contact_stiffness * penetration
    assert bead[0] == pytest.approx([expected_force, 0.0])
    assert bead[1] == pytest.approx([0.0, 0.0])
    assert bead.sum(axis=0) + reaction == pytest.approx([0.0, 0.0], abs=1e-25)
    assert torque == pytest.approx(0.0, abs=1e-25)
    assert maximum == pytest.approx(penetration)
    assert count == 1
    assert energy == pytest.approx(0.5 * cfg.contact_stiffness * penetration**2)

    # The force is the negative spatial gradient of the one-sided contact energy.
    step = 1.0e-10
    shifted_plus = positions.copy()
    shifted_minus = positions.copy()
    shifted_plus[0, 0] += step
    shifted_minus[0, 0] -= step
    energy_plus = cell_contact_forces_numpy(shifted_plus, cell, cfg)[-1]
    energy_minus = cell_contact_forces_numpy(shifted_minus, cell, cfg)[-1]
    numerical_force = -(energy_plus - energy_minus) / (2.0 * step)
    assert numerical_force == pytest.approx(expected_force, rel=1e-10)


def test_accelerated_contact_matches_numpy_and_conserves_global_moment():
    cfg = G3Config()
    cell = RigidCellState.at_origin(cfg.cell_radius)
    positions = np.array([
        [0.6 * (cfg.cell_radius - 0.1e-6), 0.8 * (cfg.cell_radius - 0.1e-6)],
        [-cfg.cell_radius - 0.1e-6, 0.0],
    ])
    accelerated = cell_contact_forces(positions, cell, cfg)
    reference = cell_contact_forces_numpy(positions, cell, cfg)
    for actual, expected in zip(accelerated, reference):
        assert actual == pytest.approx(expected, rel=1e-13, abs=1e-25)

    clutches = ClutchState.empty(0)
    zero = np.zeros_like(positions)
    force_error, torque_error = conservation_errors(
        clutches, np.empty((0, 2)), np.empty((0, 2)), zero, positions, cell,
        contact_bead_forces=accelerated[0],
        contact_cell_force=accelerated[1],
        contact_cell_torque=accelerated[2],
    )
    assert force_error < 1e-12
    assert torque_error < 1e-12


def test_contact_can_be_disabled_for_ablation():
    cfg = replace(G3Config(), contact_enabled=False)
    cell = RigidCellState.at_origin(cfg.cell_radius)
    result = cell_contact_forces(np.array([[0.5 * cfg.cell_radius, 0.0]]), cell, cfg)
    assert np.allclose(result[0], 0.0)
    assert np.allclose(result[1], 0.0)
    assert result[3:] == (0.0, 0, 0.0)


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


def test_accelerated_clutch_geometry_and_projection_match_numpy_reference():
    cfg = G3Config()
    fixture = build_fixture("single_fibre", cfg, np.random.default_rng(0))
    clutches = _single_bound_clutch(cfg, fixture, alpha=0.37)
    cell = RigidCellState.at_origin(cfg.cell_radius)
    accelerated = clutch_geometry(
        clutches, cell, fixture.initial_positions, fixture.network.fiber_bonds
    )
    reference = clutch_geometry_numpy(
        clutches, cell, fixture.initial_positions, fixture.network.fiber_bonds
    )
    for actual, expected in zip(accelerated, reference):
        assert np.allclose(actual, expected, equal_nan=True, rtol=1.0e-14, atol=1.0e-25)
    accelerated_force = project_clutch_forces(
        clutches, accelerated[0], fixture.initial_positions, fixture, cfg
    )
    reference_force = project_clutch_forces_numpy(
        clutches, reference[0], fixture.initial_positions, fixture, cfg
    )
    assert np.allclose(accelerated_force, reference_force, rtol=1.0e-13, atol=1.0e-25)


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
