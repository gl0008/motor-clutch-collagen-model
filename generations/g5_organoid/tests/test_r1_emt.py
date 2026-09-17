"""G5-R1 tests: discrete partial-EMT phenotype changes ONLY cell-cell adhesion.

Run: python -m pytest generations/g5_organoid/tests/test_r1_emt.py

Invariants: (1) emt_fraction=0 reproduces the homogeneous R0 baseline; (2) EMT does
NOT touch the clutch -- an ISOLATED single cell is identical for any emt_fraction (its
only phenotype channel, adhesion, has no neighbour to act on); (3) cohesion metrics are
sane; (4) the min / geomean adhesion pair rules both run and differ.  Personal testing,
NOT confirmed findings.
"""

from pathlib import Path
import sys
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.consistency import (  # noqa: E402
    r1_config,
    run_r0_invasion,
    run_r1_invasion,
    emt_phenotype,
    radius_of_gyration,
    detached_fraction,
    largest_connected_component_fraction,
)

SMALL = dict(n_fibers=90, n_corona_fibers=32, organoid_radius=32.0, domain_size=280.0,
             boundary_width=5.0, generation_attempts=15, dt=0.05, duration=20.0,
             sample_interval=4.0, total_pull_force=20.0)
# an ISOLATED single cell (organoid_radius < cell_spacing -> hex_centers yields 1 centre)
SINGLE = dict(n_fibers=70, n_corona_fibers=24, organoid_radius=1.0, domain_size=200.0,
              boundary_width=5.0, generation_attempts=15, dt=0.05, duration=20.0,
              sample_interval=4.0, total_pull_force=20.0)


class Phenotype(unittest.TestCase):
    def test_emt_fraction_zero_is_all_epithelial(self):
        cfg = r1_config(emt_fraction=0.0, **SMALL)
        _, centers, _, _ = __import__(
            "generations.g5_organoid.model", fromlist=["make_organoid"]).make_organoid(cfg, seed=23)
        ph = emt_phenotype(centers, cfg)
        self.assertTrue(np.allclose(ph, 1.0))

    def test_random_vs_boundary_assignment_pick_expected_count(self):
        cfg = r1_config(emt_fraction=0.25, emt_adhesion_factor=0.3, **SMALL)
        from generations.g5_organoid.model import make_organoid
        _, centers, _, _ = make_organoid(cfg, seed=23)
        for mode in ("random", "boundary"):
            ph = emt_phenotype(centers, r1_config(emt_fraction=0.25, emt_adhesion_factor=0.3,
                                                  emt_assignment=mode, **SMALL))
            self.assertEqual(int(np.sum(ph < 1.0)), int(round(0.25 * len(centers))))


class HomogeneousBaseline(unittest.TestCase):
    def test_emt_zero_reproduces_r0(self):
        """emt_fraction=0 -> all-ones adhesion -> identical to the R0 default path."""
        cfg = r1_config(emt_fraction=0.0, **SMALL)
        r1 = run_r1_invasion(cfg, seed=23)
        r0 = run_r0_invasion(cfg, seed=23)   # adhesion_scale=None -> leader_adhesion_scale (all ones)
        self.assertTrue(np.allclose(r1["centers_final"], r0["centers_final"]))
        self.assertEqual(r1["cumulative_slips"], r0["cumulative_slips"])


class ClutchIsPhenotypeIndependent(unittest.TestCase):
    def test_isolated_cell_identical_for_any_emt(self):
        """EMT changes ONLY adhesion: a lone cell (no neighbour) is invariant to EMT, so
        its clutch traction + motility are identical for emt_fraction 0 vs 1."""
        epithelial = run_r1_invasion(r1_config(emt_fraction=0.0, **SINGLE), seed=23)
        full_emt = run_r1_invasion(r1_config(emt_fraction=1.0, emt_adhesion_factor=0.25, **SINGLE), seed=23)
        self.assertEqual(epithelial["n_cells"], 1)
        self.assertTrue(np.allclose(epithelial["centers_final"], full_emt["centers_final"]))
        # clutch traction identical frame-by-frame (phenotype never entered the clutch)
        t_epi = [f["total_traction"] for f in epithelial["frames"]]
        t_emt = [f["total_traction"] for f in full_emt["frames"]]
        self.assertTrue(np.allclose(t_epi, t_emt))


class CohesionMetrics(unittest.TestCase):
    def test_metrics_are_sane_and_cohesive_at_high_adhesion(self):
        cfg = r1_config(emt_fraction=0.0, cc_adhesion=8.0, **SMALL)
        out = run_r1_invasion(cfg, seed=23)
        f = out["frames"][-1]
        self.assertGreater(f["radius_of_gyration"], 0.0)
        self.assertGreaterEqual(f["detached_fraction"], 0.0)
        self.assertLessEqual(f["detached_fraction"], 1.0)
        self.assertGreater(f["lcc_fraction"], 0.0)
        self.assertLessEqual(f["lcc_fraction"], 1.0)
        # a cohesive short run: nobody detaches, one connected cluster
        self.assertEqual(f["detached_fraction"], 0.0)
        self.assertEqual(f["lcc_fraction"], 1.0)

    def test_metric_helpers_direct(self):
        cfg = r1_config(**SMALL)
        from generations.g5_organoid.model import make_organoid
        _, centers, _, _ = make_organoid(cfg, seed=23)
        self.assertGreater(radius_of_gyration(centers), 0.0)
        self.assertEqual(detached_fraction(centers, cfg), 0.0)          # packed at rest
        self.assertEqual(largest_connected_component_fraction(centers, cfg), 1.0)


class PairRule(unittest.TestCase):
    def test_min_and_geomean_run_and_differ(self):
        """min vs geomean adhesion pair rule give different cohesion for mixed EMT."""
        common = dict(emt_fraction=0.5, emt_adhesion_factor=0.25, duration=40.0,
                      total_pull_force=25.0, **{k: v for k, v in SMALL.items()
                                                if k not in ("duration", "total_pull_force")})
        a = run_r1_invasion(r1_config(emt_pair_rule="min", **common), seed=23)
        b = run_r1_invasion(r1_config(emt_pair_rule="geomean", **common), seed=23)
        self.assertTrue(np.all(np.isfinite(a["centers_final"])))
        self.assertTrue(np.all(np.isfinite(b["centers_final"])))
        self.assertGreater(float(np.max(np.abs(a["centers_final"] - b["centers_final"]))), 1e-6)


if __name__ == "__main__":
    unittest.main(verbosity=2)
