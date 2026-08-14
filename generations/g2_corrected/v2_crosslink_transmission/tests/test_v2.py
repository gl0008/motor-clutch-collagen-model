from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from common.model import (  # noqa: E402
    CollagenConfig,
    Network,
    active_forces,
    connectivity_report,
    make_network_spec,
    run_fixed_pull,
)


class CorrectedV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = CollagenConfig(duration=2.0, sample_interval=0.5)
        cls.spec = make_network_spec(cls.cfg)

    def test_only_outer_boundary_beads_are_fixed(self):
        net = Network(self.spec, self.cfg)
        half = self.cfg.domain_size / 2
        distance = half - np.max(np.abs(net.r[net.fixed]), axis=1)
        self.assertTrue(np.all(distance <= self.cfg.boundary_width + 1e-9))
        interior = half - np.max(np.abs(net.r[~net.fixed]), axis=1)
        self.assertTrue(np.all(interior > self.cfg.boundary_width - 1e-9))

    def test_percolation_gate(self):
        report = connectivity_report(Network(self.spec, self.cfg))
        self.assertTrue(report["contact_fibers_connected"])
        self.assertGreaterEqual(
            report["connected_fraction"], self.cfg.required_connected_fraction
        )

    def test_gaussian_force_is_conserved(self):
        net = Network(self.spec, self.cfg)
        nodal, patches, vectors = active_forces(net, np.zeros(2), 5.0, 0.0)
        self.assertGreater(len(patches), 0)
        self.assertAlmostEqual(sum(x.weight for x in patches), 1.0, places=12)
        np.testing.assert_allclose(nodal.sum(axis=0), np.sum(vectors, axis=0), atol=1e-10)

    def test_frame_zero_contains_contacts_and_no_penetration(self):
        out = run_fixed_pull(self.cfg, spec=self.spec, include_crosslinks=True)
        self.assertGreater(len(out["contact_points"][0]), 0)
        self.assertGreaterEqual(out["min_cell_gap"].min(), -0.02)

    def test_crosslinks_move_noncontact_fibres(self):
        free = run_fixed_pull(self.cfg, spec=self.spec, include_crosslinks=False)
        linked = run_fixed_pull(self.cfg, spec=self.spec, include_crosslinks=True)
        contact = set(linked["contact_fibers"][0])
        remote_ids = [
            bead
            for fid, fibre in enumerate(linked["fibers"])
            if fid not in contact
            for bead in fibre
        ]
        free_d = np.linalg.norm(
            free["positions"][-1, remote_ids] - free["initial_positions"][remote_ids],
            axis=1,
        ).mean()
        linked_d = np.linalg.norm(
            linked["positions"][-1, remote_ids] - linked["initial_positions"][remote_ids],
            axis=1,
        ).mean()
        self.assertGreater(linked_d, free_d + 1e-10)


if __name__ == "__main__":
    unittest.main()
