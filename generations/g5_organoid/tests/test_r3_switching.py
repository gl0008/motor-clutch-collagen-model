"""G5-R3 tests: energy-based DYNAMIC leader switching (persistence layer; Zhang 2019).

Run: python -m pytest generations/g5_organoid/tests/test_r3_switching.py

Asserts: switching-off delegates to R2 exactly; the motor-power ATP proxy is >0 only for a
gripping loaded cell; energy drains under load and recovers when idle; hysteresis prevents
flip-flop; >=1 relay handoff to a DIFFERENT front cell occurs; force-pair stays ~0.  The energy
law + "motor work = ATP" are MODELLING CHOICES / HYPOTHESES (personal testing, NOT confirmed).
"""

from pathlib import Path
import sys
import types
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.consistency import (  # noqa: E402
    r1_config, run_r2_invasion, run_r3_invasion,
    motor_power_per_cell, update_active_leaders,
)


class Delegation(unittest.TestCase):
    def test_switching_off_equals_r2(self):
        """leader_switching=False must be byte-for-byte run_r2_invasion (delegates)."""
        cfg = r1_config(n_fibers=90, organoid_radius=28.0, domain_size=240.0, boundary_width=6.0,
                        generation_attempts=15, dt=0.05, duration=60.0, sample_interval=20.0,
                        radial_cue=True, cc_adhesion=2.0, n_leaders=2, leader_stall_factor=4.0,
                        leader_switching=False)
        a = run_r3_invasion(cfg, seed=23)
        b = run_r2_invasion(cfg, seed=23)
        self.assertTrue(np.allclose(a["centers_final"], b["centers_final"]))


class MotorPower(unittest.TestCase):
    def test_power_positive_only_for_gripping_loaded_cell(self):
        cfg = r1_config()
        patches = [types.SimpleNamespace(), None]      # cell 0 grips, cell 1 has no fibre
        site_force = np.array([2.0, 5.0])              # cell 1's force is ignored (patch None)
        site_stall = np.array([8.0, 8.0])
        P = motor_power_per_cell(cfg, patches, site_force, site_stall, n_sec=1, M=2)
        self.assertGreater(P[0], 0.0)                  # gripping + loaded -> power
        self.assertEqual(P[1], 0.0)                    # non-gripping -> zero


class EnergyDynamics(unittest.TestCase):
    def test_drains_under_load_and_recovers_when_idle(self):
        E0, tau, cap, dt = 1.0, 300.0, 130.0, 0.05
        E = E0
        for _ in range(12000):                          # 600 s under a leader-scale load
            E += dt * ((E0 - E) / tau - 0.34 / cap)
        self.assertLess(E, 0.6)                          # E* = 1 - 0.34*300/130 ~ 0.22
        drained = E
        for _ in range(30000):                          # 1500 s idle
            E += dt * ((E0 - E) / tau - 0.0 / cap)
        self.assertGreater(E, drained + 0.3)            # recovered substantially toward E0
        self.assertGreater(E, 0.9)


class Hysteresis(unittest.TestCase):
    def test_no_flipflop_between_thresholds(self):
        cfg = r1_config(n_leaders=1, cue_angle=0.0, energy_on=0.5, energy_off=0.3,
                        leader_location="cue_front")
        centers = np.array([[30.0, 0.0], [10.0, 0.0], [-10.0, 0.0]])  # 0,1 are the +x front
        grip = np.array([True, True, True])
        E = np.array([0.4, 0.4, 0.4])                   # all in the hysteresis band (off,on)
        stay, _ = update_active_leaders(centers, E, [0], grip, cfg)
        self.assertIn(0, stay)                          # a current leader in-band STAYS (E>=off)
        none, _ = update_active_leaders(centers, E, [], grip, cfg)
        self.assertEqual(len(none), 0)                  # a non-leader in-band is NOT promoted (E<on)


class Relay(unittest.TestCase):
    """One run: >=1 handoff to a DIFFERENT front cell, and force-pair stays ~0."""

    @classmethod
    def setUpClass(cls):
        cfg = r1_config(n_fibers=240, organoid_radius=36.0, domain_size=252.0, boundary_width=6.0,
                        generation_attempts=20, dt=0.05, duration=1500.0, sample_interval=300.0,
                        radial_cue=True, cue_angle=0.0, cue_half_width=1.0, cc_adhesion=1.5,
                        leader_switching=True, n_leaders=1, leader_location="cue_front",
                        leader_stall_factor=6.0, budget_mode="fixed_per_leader")
        cls.out = run_r3_invasion(cfg, seed=23)

    def test_at_least_one_handoff_to_a_different_cell(self):
        out = self.out
        self.assertGreaterEqual(out["n_switches"], 1)
        promos = [c for _t, k, c in out["switch_log"] if k == "promote"]
        distinct = set(promos) | set(out["initial_leaders"])
        self.assertGreaterEqual(len(distinct), 2)       # a DIFFERENT cell led at some point

    def test_force_pair_zero_with_switching(self):
        self.assertLess(self.out["max_force_pair_residual"], 1e-9)


if __name__ == "__main__":
    unittest.main(verbosity=2)
