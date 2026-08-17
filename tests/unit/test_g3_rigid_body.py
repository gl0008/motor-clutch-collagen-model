import numpy as np
import pytest

from g3.config import G3Config
from g3.mechanics import advance_cell, reaction_force_and_torque
from g3.state import ClutchState, RigidCellState


def test_rigid_cell_has_only_reaction_driven_translation_and_rotation():
    cfg = G3Config()
    cell = RigidCellState.at_origin(cfg.cell_radius)
    force = np.array([3.0e-12, -2.0e-12])
    torque = 4.0e-18
    advance_cell(cell, force, torque, cfg)
    assert np.allclose(cell.velocity, force / cfg.cell_drag)
    assert np.allclose(cell.center, (force / cfg.cell_drag) * cfg.dt)
    assert cell.angular_velocity == pytest.approx(torque / cfg.rotational_drag)


def test_surface_normal_force_has_zero_torque_but_evolving_vector_can_rotate():
    cfg = G3Config()
    cell = RigidCellState.at_origin(cfg.cell_radius)
    clutches = ClutchState.empty(1)
    clutches.bound[0] = True

    motor = np.array([[cfg.cell_radius, 0.0]])
    clutches.force_vector[0] = np.array([-2.0e-12, 0.0])
    _, radial_torque = reaction_force_and_torque(clutches, motor, cell)
    assert radial_torque == pytest.approx(0.0)

    clutches.force_vector[0] = np.array([-2.0e-12, -1.0e-12])
    _, off_axis_torque = reaction_force_and_torque(clutches, motor, cell)
    assert off_axis_torque > 0.0


def test_mirroring_attachment_and_force_reverses_torque_sign():
    cfg = G3Config()
    cell = RigidCellState.at_origin(cfg.cell_radius)
    clutches = ClutchState.empty(1)
    clutches.bound[0] = True
    motor = np.array([[cfg.cell_radius, 1.0e-6]])
    clutches.force_vector[0] = np.array([-1.0e-12, -2.0e-12])
    _, torque = reaction_force_and_torque(clutches, motor, cell)

    mirrored_motor = motor.copy()
    mirrored_motor[:, 1] *= -1.0
    clutches.force_vector[:, 1] *= -1.0
    _, mirrored_torque = reaction_force_and_torque(clutches, mirrored_motor, cell)
    assert mirrored_torque == pytest.approx(-torque)
