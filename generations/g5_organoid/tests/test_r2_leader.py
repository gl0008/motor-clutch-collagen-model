"""G5-R2 tests: static leader NUMBER + LOCATION (functional role, decoupled from EMT).

Run: python -m pytest generations/g5_organoid/tests/test_r2_leader.py

Asserts the R2 invariants: the imposed one-sided cue is length-preserving and raises the
CUE-SIDE radial order only; n_leaders=0 reproduces the R1 baseline; leaders are localized
(cue_front) vs spread (perimeter); the leader motor budget is matched_total (constant vs N_L)
or fixed_per_leader (grows); leader traction EMERGES from per-site stall (not a multiply); and
the force-pair stays 0 with leaders on.  Personal testing, NOT confirmed findings.
"""

from pathlib import Path
import sys
import math
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.model import hex_centers  # noqa: E402
from generations.g5_organoid.consistency import (  # noqa: E402
    r1_config, run_r1_invasion, run_r2_invasion,
    make_cued_organoid, leader_ids, leader_site_stall,
    _random_isotropic_spec, _cue_positions,
)
from common.model import Network, NetworkSpec  # noqa: E402

# small fast fixture (isotropic-random needs ~250 fibres for grip; force-pair holds anyway).
SMALL = dict(n_fibers=250, organoid_radius=28.0, domain_size=240.0, boundary_width=6.0,
             generation_attempts=20, dt=0.05, duration=100.0, sample_interval=50.0)


def _sector_radial_order(network, center, angle, half, band):
    """Mean radial order (2cos^2-1) of near-field segments whose midpoint azimuth is within
    ``half`` of ``angle`` and radius in the near-field band."""
    i, j = network.edges.T
    seg = network.r[j] - network.r[i]
    tan = seg / np.maximum(np.linalg.norm(seg, axis=1), 1e-12)[:, None]
    mid = 0.5 * (network.r[i] + network.r[j])
    rad = mid - center
    r = np.linalg.norm(rad, axis=1)
    er = rad / np.maximum(r, 1e-12)[:, None]
    order = 2.0 * np.square(np.sum(tan * er, axis=1)) - 1.0
    az = np.arctan2(rad[:, 1], rad[:, 0])
    dang = np.abs(np.arctan2(np.sin(az - angle), np.cos(az - angle)))
    m = (dang <= half) & (r >= band[0]) & (r < band[1])
    return float(np.mean(order[m])) if np.any(m) else float("nan")


class ImposedCue(unittest.TestCase):
    def test_cue_is_length_preserving(self):
        cfg = r1_config(radial_cue=True, cue_angle=0.0, **SMALL)
        net, centers, gap, report = make_cued_organoid(cfg, seed=23)
        self.assertGreater(report["cue"]["tract_fibers"], 0)          # the cue did something
        self.assertLess(report["cue"]["max_contour_length_change"], 1e-6)  # rigid rotation only

    def test_cue_raises_radial_order_on_cue_side_only(self):
        """Isolate the cue: on the SAME base network, the cue raises the cue-sector radial
        order while the opposite sector is ~unchanged (before vs after avoids random per-
        sector baseline variance)."""
        cfg = r1_config(radial_cue=True, cue_angle=0.0, cue_half_width=0.6, cue_band=45.0, **SMALL)
        centers = hex_centers(cfg.organoid_radius, cfg.cell_spacing)
        outer = float(np.max(np.linalg.norm(centers, axis=1))) + cfg.cell_radius
        band = (outer, outer + cfg.cue_band)
        spec = _random_isotropic_spec(cfg, centers, 23)
        net_r = Network(spec, cfg, crosslinks=[])
        cue_before = _sector_radial_order(net_r, np.zeros(2), 0.0, cfg.cue_half_width, band)
        opp_before = _sector_radial_order(net_r, np.zeros(2), math.pi, cfg.cue_half_width, band)
        pts, tract, rep = _cue_positions(spec.positions, spec.fibers, spec.fixed,
                                         centers, cfg, outer)
        self.assertGreater(len(tract), 0)
        net_c = Network(NetworkSpec(pts, spec.fibers, spec.fixed, spec.contact_fibers,
                                    spec.seed_used), cfg, crosslinks=[])
        cue_after = _sector_radial_order(net_c, np.zeros(2), 0.0, cfg.cue_half_width, band)
        opp_after = _sector_radial_order(net_c, np.zeros(2), math.pi, cfg.cue_half_width, band)
        self.assertGreater(cue_after, cue_before + 0.1)   # cue raised the cue-side order
        self.assertLess(abs(opp_after - opp_before), 0.05)  # opposite side ~unchanged


class LeaderSelection(unittest.TestCase):
    def test_zero_leaders_is_r1_baseline(self):
        cfg = r1_config(n_leaders=0, emt_fraction=0.0, **SMALL)
        a = run_r2_invasion(cfg, seed=23)
        b = run_r1_invasion(cfg, seed=23)
        self.assertTrue(np.allclose(a["centers_final"], b["centers_final"]))

    def test_cue_front_is_localized_and_differs_from_perimeter(self):
        centers = hex_centers(28.0, 18.0)
        front = leader_ids(centers, r1_config(n_leaders=3, leader_location="cue_front",
                                              cue_angle=0.0, organoid_radius=28.0))
        perim = leader_ids(centers, r1_config(n_leaders=3, leader_location="perimeter",
                                              organoid_radius=28.0))
        self.assertEqual(len(front), 3)
        # cue_front leaders sit on the +x side (outward azimuth near cue_angle=0)
        az = np.arctan2(centers[front, 1], centers[front, 0])
        self.assertLess(float(np.max(np.abs(az))), 1.2)      # clustered near angle 0
        self.assertFalse(np.array_equal(np.sort(front), np.sort(perim)))  # differs from spread


class MotorBudget(unittest.TestCase):
    def _leader_stall_sum(self, nl, mode):
        cfg = r1_config(n_leaders=nl, leader_stall_factor=3.0, budget_mode=mode,
                        organoid_radius=28.0)
        centers = hex_centers(cfg.organoid_radius, cfg.cell_spacing)
        ss = leader_site_stall(centers, cfg)
        lead = leader_ids(centers, cfg)
        n_sec = cfg.n_contact_sectors
        sites = np.concatenate([np.arange(c * n_sec, (c + 1) * n_sec) for c in lead])
        return float(ss[sites].sum())

    def test_matched_total_is_constant_across_leader_count(self):
        sums = [self._leader_stall_sum(nl, "matched_total") for nl in (1, 2, 3, 5)]
        self.assertTrue(np.allclose(sums, sums[0]))          # total front budget fixed

    def test_fixed_per_leader_grows_with_count(self):
        sums = [self._leader_stall_sum(nl, "fixed_per_leader") for nl in (1, 2, 3, 5)]
        self.assertTrue(all(sums[k + 1] > sums[k] for k in range(len(sums) - 1)))


class StallTraction(unittest.TestCase):
    def test_leader_traction_emerges_from_stall(self):
        """fixed_per_leader stall>1 raises a leader's clutch traction; factor=1 removes it."""
        base = dict(n_leaders=3, leader_location="cue_front", budget_mode="fixed_per_leader", **SMALL)
        out1 = run_r2_invasion(r1_config(leader_stall_factor=1.0, **base), seed=23)
        out3 = run_r2_invasion(r1_config(leader_stall_factor=3.0, **base), seed=23)
        lead = np.array(out1["leader_ids"])
        f1 = out1["clutch_per_cell_final"][lead].mean()
        f3 = out3["clutch_per_cell_final"][lead].mean()
        self.assertGreater(f3, f1 + 1.0)                     # higher stall -> higher traction
        self.assertLess(out3["max_force_pair_residual"], 1e-9)   # force-pair still exact
        self.assertLess(out1["max_force_pair_residual"], 1e-9)


if __name__ == "__main__":
    unittest.main(verbosity=2)
