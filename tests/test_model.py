import math
import unittest
from dataclasses import replace

import numpy as np

from collagen_model.clutch import MotorClutchParameters
from collagen_model.network import make_diluted_triangular_network
from collagen_model.simulation import SimulationConfig, run_condition
from collagen_model.sls import StandardLinearSolid
from collagen_model.single_protrusion import (
    SingleProtrusionConfig,
    counter_uniform,
    run_single_protrusion_pair,
)


class StandardLinearSolidTests(unittest.TestCase):
    def test_step_strain_relaxation_matches_analytic_solution(self):
        bond = StandardLinearSolid(k0=10.0, kinf=2.0, tau=3.0)
        bond.set_extension(0.5, instantaneous=True)
        dt = 0.01
        for step in range(1001):
            time = step * dt
            expected = (2.0 + 8.0 * math.exp(-time / 3.0)) * 0.5
            self.assertAlmostEqual(bond.force, expected, places=10)
            if step < 1000:
                bond.set_extension(0.5, dt=dt)

    def test_invalid_parameters_are_rejected(self):
        with self.assertRaises(ValueError):
            StandardLinearSolid(k0=1.0, kinf=2.0, tau=1.0)


class NetworkTests(unittest.TestCase):
    def test_reference_network_has_zero_internal_force(self):
        network = make_diluted_triangular_network(nx=8, ny=7, seed=4)
        forces = network.internal_forces(0.01)
        self.assertLess(float(np.max(np.abs(forces))), 1.0e-10)

    def test_network_contains_filament_bending_triplets(self):
        network = make_diluted_triangular_network(nx=8, ny=7, seed=4)
        self.assertGreater(network.n_edges, network.n_beads)
        self.assertGreater(len(network.bending_triplets), 0)

    def test_affinely_held_network_relaxes(self):
        network = make_diluted_triangular_network(nx=8, ny=7, seed=4, k0=10.0, kinf=2.0, tau=1.0)
        low, _ = network.bounds
        network.positions[:, 0] = low[0] + 1.05 * (network.positions[:, 0] - low[0])
        network.instantaneous_sync()
        initial = float(np.mean(network.axial_force_magnitudes()))
        for _ in range(500):
            network.internal_forces(0.01, zero_fixed=False)
        relaxed = float(np.mean(network.axial_force_magnitudes()))
        self.assertLess(relaxed, initial)
        self.assertAlmostEqual(relaxed / initial, 0.2 + 0.8 * math.exp(-5.0), places=9)


class CoupledSimulationTests(unittest.TestCase):
    def test_short_simulation_is_finite(self):
        config = SimulationConfig(duration=0.5, dt=0.02, nx=10, ny=8)
        result = run_condition(
            tau=5.0,
            seed=2,
            config=config,
            clutch_parameters=MotorClutchParameters(n_modules=4, clutches_per_module=2),
        )
        self.assertTrue(np.all(np.isfinite(result["x"])))
        self.assertTrue(np.all(np.isfinite(result["final_network_positions"])))
        self.assertEqual(result["time"][0], 0.0)
        self.assertAlmostEqual(result["time"][-1], 0.5)


class SingleProtrusionTests(unittest.TestCase):
    def test_serial_sls_stiffness_is_resolution_independent(self):
        total_extension = 0.4
        time = 0.73
        target = (0.5 + (2.0 - 0.5) * math.exp(-time / 0.1)) * total_extension
        for n_bonds in (5, 9, 19):
            bond = StandardLinearSolid(
                k0=n_bonds * 2.0,
                kinf=n_bonds * 0.5,
                tau=0.1,
            )
            bond.set_extension(total_extension / n_bonds, instantaneous=True)
            bond.set_extension(total_extension / n_bonds, dt=time)
            self.assertAlmostEqual(bond.force, target, places=12)

    def test_counter_random_stream_is_addressable(self):
        draw = counter_uniform(17, 23, 4, 1)
        self.assertEqual(draw, counter_uniform(17, 23, 4, 1))
        self.assertNotEqual(draw, counter_uniform(17, 23, 4, 0))
        self.assertGreaterEqual(draw, 0.0)
        self.assertLess(draw, 1.0)

    def test_no_binding_means_zero_traction(self):
        config = SingleProtrusionConfig(
            on_rate=0.0,
            duration=0.2,
            dt=0.002,
            sample_interval=0.01,
        )
        pair = run_single_protrusion_pair(config, seed=2)
        np.testing.assert_allclose(pair["fast"]["total_traction"], 0.0)
        np.testing.assert_allclose(pair["slow"]["total_traction"], 0.0)

    def test_elastic_limit_removes_tau_difference(self):
        config = SingleProtrusionConfig(
            k0_chain=1.5,
            kinf_chain=1.5,
            duration=1.0,
            dt=0.002,
            sample_interval=0.01,
        )
        pair = run_single_protrusion_pair(config, seed=8)
        for key in (
            "bead_positions",
            "clutch_bound",
            "clutch_force",
            "total_traction",
            "actin_velocity",
        ):
            np.testing.assert_allclose(pair["fast"][key], pair["slow"][key])

    def test_force_independent_off_rate_gives_identical_cluster_states(self):
        config = SingleProtrusionConfig(
            force_dependent_off=False,
            duration=2.0,
            dt=0.002,
            sample_interval=0.01,
        )
        pair = run_single_protrusion_pair(config, seed=11)
        np.testing.assert_array_equal(
            pair["fast"]["clutch_bound"],
            pair["slow"]["clutch_bound"],
        )
        self.assertEqual(pair["fast"]["episodes"], pair["slow"]["episodes"])

    def test_single_protrusion_result_has_requested_observables(self):
        config = replace(
            SingleProtrusionConfig(),
            duration=0.5,
            dt=0.002,
            sample_interval=0.01,
        )
        result = run_single_protrusion_pair(config, seed=5)["fast"]
        n_samples = len(result["time"])
        self.assertEqual(result["bead_positions"].shape, (n_samples, 10))
        self.assertEqual(result["sls_q"].shape, (n_samples, 9))
        self.assertEqual(result["clutch_bound"].shape, (n_samples, 12))
        self.assertEqual(result["clutch_force"].shape, (n_samples, 12))
        self.assertEqual(result["phase"].shape, (n_samples,))
        self.assertTrue(np.all(np.isfinite(result["total_traction"])))


if __name__ == "__main__":
    unittest.main()
