"""G5-R4 tests: FOLLOWER contact-guidance along the leader-remodelled aligned collagen track.

Run: python -m pytest generations/g5_organoid/tests/test_r4_guidance.py

Asserts the R4 invariants: guidance is ZERO on an isotropic matrix (S~0) and on epithelial
(non-guided) cells; guidance is an external CELL force only, so the force-pair stays 0 with
guidance on; ``follower_guidance=False`` reproduces the R2 path byte-for-byte; and the honest
strand metric (`connected_components`/`strand_report`) is self-consistent.  Personal testing,
NOT confirmed findings.
"""

from pathlib import Path
import sys
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.consistency import (  # noqa: E402
    r1_config, run_r2_invasion, run_r4_invasion,
    contact_guidance_forces, make_random_organoid, make_cued_organoid,
    connected_components, strand_report, leader_ids,
)

SMALL = dict(n_fibers=250, organoid_radius=28.0, domain_size=240.0, boundary_width=6.0,
             generation_attempts=20, dt=0.05, duration=100.0, sample_interval=50.0)


class TestR4Guidance(unittest.TestCase):
    def test_guidance_off_matches_r2(self):
        """follower_guidance=False -> run_r4_invasion == run_r2_invasion (same centres)."""
        over = dict(SMALL, radial_cue=True, n_leaders=3, leader_location="cue_front",
                    leader_stall_factor=6.0, budget_mode="fixed_per_leader")
        r2 = run_r2_invasion(r1_config(**over), seed=23)
        r4 = run_r4_invasion(r1_config(**dict(over, follower_guidance=False)), seed=23)
        np.testing.assert_allclose(r2["centers_final"], r4["centers_final"], atol=1e-9)

    def test_guidance_zero_when_disabled_or_zero_strength(self):
        cfg = r1_config(**dict(SMALL, follower_guidance=False, guidance_strength=20.0))
        net, centers, _gap, _rep = make_random_organoid(cfg, seed=23)
        g = contact_guidance_forces(centers, net, centers, cfg)
        self.assertEqual(float(np.abs(g).max()), 0.0)  # off -> no force
        cfg2 = r1_config(**dict(SMALL, follower_guidance=True, guidance_strength=0.0))
        g2 = contact_guidance_forces(centers, net, centers, cfg2)
        self.assertEqual(float(np.abs(g2).max()), 0.0)  # zero strength -> no force

    def test_guidance_high_threshold_gives_zero(self):
        """align_min above the achievable nematic order S<=1 -> no cell fires."""
        cfg = r1_config(**dict(SMALL, radial_cue=True, follower_guidance=True,
                               guidance_strength=20.0, guidance_align_min=1.01))
        net, centers, _gap, _rep = make_cued_organoid(cfg, seed=23)
        g = contact_guidance_forces(centers, net, centers, cfg)
        self.assertEqual(float(np.abs(g).max()), 0.0)

    def test_guidance_skips_leaders(self):
        """Leaders build the track -> they do not feel guidance (guidance_leaders=False)."""
        cfg = r1_config(**dict(SMALL, radial_cue=True, n_leaders=3, leader_location="cue_front",
                               follower_guidance=True, guidance_strength=20.0,
                               guidance_align_min=0.0))
        net, centers, _gap, _rep = make_cued_organoid(cfg, seed=23)
        g = contact_guidance_forces(centers, net, centers, cfg)
        for c in leader_ids(centers, cfg):
            self.assertEqual(float(np.linalg.norm(g[int(c)])), 0.0)

    def test_force_pair_zero_with_guidance_on(self):
        over = dict(SMALL, radial_cue=True, n_leaders=3, leader_location="cue_front",
                    leader_stall_factor=6.0, budget_mode="fixed_per_leader", cc_adhesion=1.5,
                    follower_guidance=True, guidance_strength=20.0, guidance_align_min=0.30)
        out = run_r4_invasion(r1_config(**over), seed=23)
        self.assertLess(out["max_force_pair_residual"], 1e-9)

    def test_strand_report_partitions_cells(self):
        """connected_components covers every cell exactly once; strand_report fields consistent."""
        cfg = r1_config(**dict(SMALL, radial_cue=True, n_leaders=3, leader_location="cue_front"))
        net, centers, _gap, _rep = make_cued_organoid(cfg, seed=23)
        comps = connected_components(centers, cfg)
        flat = sorted(i for comp in comps for i in comp)
        self.assertEqual(flat, list(range(len(centers))))
        rep = strand_report(centers, centers, cfg)
        self.assertEqual(rep["n_components"], len(comps))
        self.assertGreaterEqual(rep["leader_comp_size"], 1)

    def test_guidance_not_harmful_and_force_consistent(self):
        """Robust INVARIANT (not the scale-dependent benefit): with guidance on, the mechanics
        stay force-consistent and guidance does not make detachment WORSE than the control.

        The connectivity BENEFIT (leaders keep followers instead of escaping alone) is an
        EMERGENT, scale/time-dependent result -- it appears at the 2 h / radius-40 headline once
        the unguided control fragments (guided 10 followers vs unguided 0; see R4_findings.md).
        At this small/short fixture NEITHER run has fragmented yet, so the benefit is absent (and
        can be marginally negative from isotropic drift) -- asserting it here would be asserting a
        scale-fragile outcome as a universal invariant.  So we assert only what IS robust."""
        base = dict(SMALL, duration=1500.0, sample_interval=300.0, radial_cue=True,
                    n_leaders=3, leader_location="cue_front", leader_stall_factor=6.0,
                    budget_mode="fixed_per_leader", cc_adhesion=1.5)
        off = run_r4_invasion(r1_config(**dict(base, follower_guidance=False)), seed=23)
        on = run_r4_invasion(r1_config(**dict(base, follower_guidance=True,
                                              guidance_strength=20.0, guidance_align_min=0.30)), seed=23)
        self.assertLess(on["max_force_pair_residual"], 1e-9)         # mechanics stay consistent
        # guidance must not CREATE detachment the control did not have
        self.assertFalse(on["leaders_detached"] and not off["leaders_detached"])


if __name__ == "__main__":
    unittest.main()
