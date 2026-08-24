import math
from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))

from model import (  # noqa: E402
    G4Config,
    _assemble_g4_spec,
    active_from_patches,
    bell_off_rate,
    bell_summary,
    build_g4_network,
    direct_contact_patches,
    graph_distance_from_direct,
)


class G4MechanicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = G4Config(n_fibers=32, bead_spacing=1.5, duration=0.1)
        cls.spec = _assemble_g4_spec(cls.cfg, cls.cfg.seed)

    def test_probability_sets_are_nested(self):
        low, _, _ = build_g4_network(replace(self.cfg, crosslink_probability=0.2), spec=self.spec)
        high, _, _ = build_g4_network(replace(self.cfg, crosslink_probability=0.7), spec=self.spec)
        low_ids = {(x.edge_a, x.alpha_a, x.edge_b, x.alpha_b) for x in low.crosslinks}
        high_ids = {(x.edge_a, x.alpha_a, x.edge_b, x.alpha_b) for x in high.crosslinks}
        self.assertTrue(low_ids <= high_ids)

    def test_only_direct_patch_beads_receive_active_force(self):
        network, _, _ = build_g4_network(self.cfg, spec=self.spec)
        patches = direct_contact_patches(network, np.zeros(2))
        active, _, _ = active_from_patches(network, patches, 10.0, np.zeros(2))
        forced = set(np.flatnonzero(np.linalg.norm(active, axis=1) > 0.0))
        eligible = set()
        for patch in patches:
            eligible.update(network.edges[patch.edge].tolist())
        self.assertTrue(forced)
        self.assertTrue(forced <= eligible)
        direct = {p.fiber for p in patches}
        distance = graph_distance_from_direct(network, direct)
        self.assertTrue(all(distance[fid] == 0 for fid in direct))

    def test_gaussian_distribution_conserves_scalar_force(self):
        network, _, _ = build_g4_network(self.cfg, spec=self.spec)
        patches = direct_contact_patches(network, np.zeros(2))
        _, _, vectors = active_from_patches(network, patches, 10.0, np.zeros(2))
        self.assertAlmostEqual(float(np.linalg.norm(vectors, axis=1).sum()), 10.0, places=10)

    def test_bell_landmarks_are_continuous_not_thresholds(self):
        summary = bell_summary(self.cfg)
        self.assertAlmostEqual(summary["extensions"][1], self.cfg.bell_force / self.cfg.clutch_stiffness)
        self.assertGreater(bell_off_rate(1.5001, self.cfg), bell_off_rate(1.4999, self.cfg))
        self.assertAlmostEqual(
            summary["median_lifetimes"][0], math.log(2.0) / self.cfg.clutch_off_rate0
        )

    def test_probability_zero_has_no_indirect_graph_path(self):
        network, _, _ = build_g4_network(replace(self.cfg, crosslink_probability=0.0), spec=self.spec)
        patches = direct_contact_patches(network, np.zeros(2))
        distance = graph_distance_from_direct(network, [p.fiber for p in patches])
        self.assertEqual(len(network.crosslinks), 0)
        self.assertTrue(np.all(np.isin(distance, [-1, 0])))


if __name__ == "__main__":
    unittest.main()
