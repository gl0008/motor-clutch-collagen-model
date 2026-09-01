from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO))

from generations.g4_interactive_calibration.model import (  # noqa: E402
    _assemble_g4_spec,
    active_from_patches,
    build_g4_network,
    direct_contact_patches,
    fast_advance,
)
from generations.g4_v2_multiscale.model import (  # noqa: E402
    ClutchState,
    FastStepper,
    G4V2Config,
    _shared_step,
    build_v2_network,
    counter_uniforms,
    run_elastic,
    shared_cluster_ensemble,
    shared_load_hazard,
)


class G4V2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = G4V2Config(
            n_fibers=32, bead_spacing=1.5, duration=0.2,
            sample_interval=0.1, metric_sample_interval=0.1,
        )
        cls.spec = _assemble_g4_spec(cls.cfg, cls.cfg.seed)

    def test_compiled_integrator_matches_frozen_g4_v1_step(self):
        old, _, _ = build_g4_network(self.cfg, spec=self.spec)
        new, _, _ = build_v2_network(self.cfg, spec=self.spec)
        center = np.zeros(2)
        patches = direct_contact_patches(old, center)
        active_old, _, _ = active_from_patches(old, patches, 8.0, center)
        active_new = active_old.copy()
        v_old, total_old, _ = fast_advance(old, active_old, center, self.cfg.dt)
        v_new, total_new = FastStepper(new).step(active_new, center, self.cfg.dt)
        np.testing.assert_allclose(new.r, old.r, rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(v_new, v_old, rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(total_new, total_old, rtol=1e-11, atol=1e-11)

    def test_mobile_boundary_releases_only_anchor_constraint(self):
        anchored, _, anchored_report = build_v2_network(self.cfg, spec=self.spec)
        mobile_cfg = replace(self.cfg, boundary_mode="mobile")
        mobile, _, mobile_report = build_v2_network(mobile_cfg, spec=self.spec)
        self.assertGreater(anchored.fixed.sum(), 0)
        self.assertEqual(mobile.fixed.sum(), 0)
        np.testing.assert_allclose(anchored.r0, mobile.r0)
        np.testing.assert_array_equal(anchored.edges, mobile.edges)
        self.assertEqual(anchored_report["boundary_mode"], "anchored")
        self.assertEqual(mobile_report["boundary_mode"], "mobile")

    def test_shared_hazard_increases_after_one_failure_at_fixed_site_load(self):
        force = np.asarray([8.0, 8.0])
        _, per12, total12 = shared_load_hazard(force, np.asarray([12, 12]), self.cfg)
        _, per11, total11 = shared_load_hazard(force, np.asarray([11, 11]), self.cfg)
        self.assertTrue(np.all(per11 > per12))
        # The per-survivor risk is the cascade mechanism; the total hazard need
        # not be monotone for every possible i and load.
        self.assertTrue(np.all(total12 > 0.0))
        self.assertTrue(np.all(total11 > 0.0))

    def test_shared_site_detaches_only_at_zero_bound_count(self):
        state = ClutchState(
            bound=np.ones((self.cfg.n_contact_sectors, self.cfg.n_clutches_per_site), dtype=bool),
            extension=np.zeros((self.cfg.n_contact_sectors, self.cfg.n_clutches_per_site)),
            site_extension=np.full(self.cfg.n_contact_sectors, 0.5),
        )
        substrate = np.zeros(self.cfg.n_contact_sectors)
        u_on = np.ones_like(state.extension)
        u_off = np.ones_like(state.extension)
        u_off[0, 0] = 0.0
        _, site_force, _, _, failed = _shared_step(
            self.cfg, state, substrate, u_on, u_off
        )
        self.assertEqual(state.bound[0].sum(), self.cfg.n_clutches_per_site - 1)
        self.assertFalse(failed[0])
        self.assertGreater(site_force[0], 0.0)

        # Force every survivor at site 0 to rupture; other sites remain bound.
        u_off[:] = 1.0
        u_off[0, :] = 0.0
        _, site_force, _, _, failed = _shared_step(
            self.cfg, state, substrate, u_on, u_off
        )
        self.assertEqual(state.bound[0].sum(), 0)
        self.assertTrue(failed[0])
        self.assertEqual(site_force[0], 0.0)
        self.assertEqual(state.site_extension[0], 0.0)

    def test_counter_stream_is_addressable_and_repeatable(self):
        a = counter_uniforms(self.cfg, 7, 1)
        b = counter_uniforms(self.cfg, 7, 1)
        c = counter_uniforms(self.cfg, 8, 1)
        np.testing.assert_array_equal(a, b)
        self.assertFalse(np.array_equal(a, c))
        self.assertTrue(np.all((a >= 0.0) & (a < 1.0)))

    def test_dt_refinement_preserves_elastic_trajectory_direction(self):
        coarse_cfg = replace(
            self.cfg, duration=2.0, sample_interval=2.0,
            metric_sample_interval=2.0, force_ramp_time=0.5,
            total_pull_force=24.0, dt=0.05,
        )
        fine_cfg = replace(coarse_cfg, dt=0.025)
        coarse = run_elastic(coarse_cfg, spec=self.spec)
        fine = run_elastic(fine_cfg, spec=self.spec)
        dc = coarse["network"].r - coarse["network"].r0
        df = fine["network"].r - fine["network"].r0
        self.assertGreater(float(np.vdot(dc, df)), 0.0)
        relative_difference = float(np.linalg.norm(dc - df) / max(np.linalg.norm(df), 1e-12))
        self.assertLess(relative_difference, 0.05)

    def test_small_shared_cluster_ensemble_reports_real_episodes(self):
        out = shared_cluster_ensemble(self.cfg, trials=4, duration=20.0, seed0=30)
        self.assertEqual(out["trials"], 4)
        self.assertEqual(out["first_failures"].shape, (4,))
        self.assertGreaterEqual(out["representative_seed"], 30)


if __name__ == "__main__":
    unittest.main()
