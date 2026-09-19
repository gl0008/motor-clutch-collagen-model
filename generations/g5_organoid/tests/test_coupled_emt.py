"""Coupled EMT-leader variant tests: emt_traction_factor couples ↑traction to the ↓adhesion
partial-EMT cells (BIOLOGY_AUDIT.md B).  factor=1 must reduce exactly to R1 (adhesion-only)."""

from pathlib import Path
import sys
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.consistency import (  # noqa: E402
    r1_config, run_coupled_emt_invasion, run_r1_invasion,
    emt_phenotype, emt_cell_ids, emt_site_stall, make_random_organoid,
)

SMALL = dict(n_fibers=200, organoid_radius=28.0, domain_size=240.0, boundary_width=6.0,
             generation_attempts=20, dt=0.05, duration=100.0, sample_interval=25.0,
             cc_adhesion=6.0, emt_fraction=0.5, emt_adhesion_factor=0.15, emt_assignment="boundary")


class CoupledEMT(unittest.TestCase):
    def test_cell_ids_match_phenotype(self):
        cfg = r1_config(**SMALL)
        net, centers, _g, _r = make_random_organoid(cfg, seed=23)
        ids = set(int(x) for x in emt_cell_ids(centers, cfg))
        nz = set(int(x) for x in np.flatnonzero(emt_phenotype(centers, cfg) != 1.0))
        self.assertEqual(ids, nz)

    def test_factor_one_reduces_to_r1(self):
        """emt_traction_factor=1.0 must be byte-for-byte run_r1_invasion (decoupled)."""
        a = run_coupled_emt_invasion(r1_config(emt_traction_factor=1.0, **SMALL), seed=23)
        b = run_r1_invasion(r1_config(**SMALL), seed=23)
        self.assertTrue(np.allclose(a["centers_final"], b["centers_final"]))

    def test_site_stall_boosts_only_emt_cells(self):
        cfg = r1_config(emt_traction_factor=5.0, **SMALL)
        net, centers, _g, _r = make_random_organoid(cfg, seed=23)
        stall = emt_site_stall(centers, cfg)
        n_sec = cfg.n_contact_sectors
        emt = set(int(x) for x in emt_cell_ids(centers, cfg))
        for c in range(len(centers)):
            block = stall[c * n_sec:(c + 1) * n_sec]
            expect = cfg.motor_stall_per_site * (5.0 if c in emt else 1.0)
            self.assertTrue(np.allclose(block, expect))

    def test_coupled_keeps_force_pair_zero(self):
        out = run_coupled_emt_invasion(r1_config(emt_traction_factor=5.0, **SMALL), seed=23)
        self.assertLess(out["max_force_pair_residual"], 1e-9)
        self.assertEqual(len(out["emt_cell_ids"]), int(round(0.5 * out["n_cells"])))


if __name__ == "__main__":
    unittest.main(verbosity=2)
