from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from model import (  # noqa: E402
    PlasticityConfig,
    form_new_weak_links,
    run_v4_pair,
)
sys.path.insert(0, str(HERE.parent))
from common.model import Network, make_network_spec  # noqa: E402


class ContactPlasticityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = PlasticityConfig(
            duration=1.0,
            sample_interval=0.5,
            weak_link_capture_radius=1.5,
            weak_link_alignment_deg=90.0,
        )
        cls.spec = make_network_spec(cls.cfg)

    def test_v4_is_explicitly_gated(self):
        with self.assertRaises(RuntimeError):
            run_v4_pair(self.cfg, gates_passed=False)

    def test_new_links_join_different_fibres_and_are_stress_free(self):
        net = Network(self.spec, self.cfg)
        created = form_new_weak_links(
            net, self.cfg, time=1.0, limit=6, require_new_contact=False
        )
        self.assertGreater(len(created), 0)
        for link in created:
            fa = net.edge_fiber[link.edge_a]
            fb = net.edge_fiber[link.edge_b]
            self.assertNotEqual(fa, fb)
            pa, pb = net.crosslink_endpoints(link)
            np.testing.assert_allclose((pb - pa) - link.rest_vector, 0.0, atol=1e-12)
            self.assertEqual(link.kind, "new-weak")


if __name__ == "__main__":
    unittest.main()
