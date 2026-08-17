import numpy as np

from g3.elastic import (
    bending_energy,
    bending_forces,
    bending_forces_numpy,
    extensional_energy,
    extensional_forces,
    harmonic_bond_forces_numpy,
    overdamped_step,
)


def _numerical_force(energy, positions, bead, coordinate, epsilon=1.0e-9):
    plus = positions.copy()
    minus = positions.copy()
    plus[bead, coordinate] += epsilon
    minus[bead, coordinate] -= epsilon
    return -(energy(plus) - energy(minus)) / (2.0 * epsilon)


def test_extensional_force_is_equal_opposite_and_matches_energy_gradient():
    positions = np.array([[0.0, 0.0], [1.2e-6, 0.2e-6]])
    bonds = np.array([[0, 1]], dtype=int)
    stiffness = 4.0e-3
    force = extensional_forces(positions, bonds, 1.0e-6, stiffness)
    energy = lambda x: extensional_energy(x, bonds, 1.0e-6, stiffness)
    assert np.allclose(force.sum(axis=0), 0.0, atol=1.0e-20)
    assert np.isclose(force[1, 1], _numerical_force(energy, positions, 1, 1), rtol=1.0e-6)


def test_bending_force_conserves_internal_force_and_matches_energy_gradient():
    positions = np.array([[0.0, 0.0], [1.0e-6, 0.0], [1.8e-6, 0.5e-6]])
    triples = np.array([[0, 1, 2]], dtype=int)
    stiffness = 8.27e-20
    force = bending_forces(positions, triples, 0.0, stiffness)
    energy = lambda x: bending_energy(x, triples, 0.0, stiffness)
    assert np.allclose(force.sum(axis=0), 0.0, atol=1.0e-24)
    assert np.isclose(force[2, 1], _numerical_force(energy, positions, 2, 1), rtol=1.0e-5)


def test_overdamped_step_uses_force_over_drag_and_has_no_hidden_noise():
    positions = np.array([[1.0, 2.0]])
    forces = np.array([[4.0, -2.0]])
    updated = overdamped_step(positions, forces, drag=2.0, dt=0.5)
    assert np.array_equal(updated, np.array([[2.0, 1.5]]))


def test_accelerated_elastic_kernels_match_numpy_reference():
    rng = np.random.default_rng(18)
    positions = np.column_stack((np.arange(7) * 1.0e-6, rng.normal(0.0, 0.08e-6, 7)))
    bonds = np.column_stack((np.arange(6), np.arange(1, 7))).astype(int)
    triples = np.column_stack((np.arange(5), np.arange(1, 6), np.arange(2, 7))).astype(int)
    assert np.allclose(
        extensional_forces(positions, bonds, 1.0e-6, 4.0e-3),
        harmonic_bond_forces_numpy(positions, bonds, 1.0e-6, 4.0e-3),
        rtol=1.0e-13,
        atol=1.0e-25,
    )
    assert np.allclose(
        bending_forces(positions, triples, 0.0, 8.27e-20),
        bending_forces_numpy(positions, triples, 0.0, 8.27e-20),
        rtol=1.0e-12,
        atol=1.0e-25,
    )
