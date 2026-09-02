"""Generation 5 tests: Stage-A scaffold balance and Stage-B radial alignment.

Run: python -m pytest generations/g5_organoid/tests/test_g5.py
"""

from pathlib import Path
import sys
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.model import (  # noqa: E402
    OrganoidConfig,
    hex_centers,
    make_organoid,
    cell_cell_forces,
    multi_cell_repulsion,
    organoid_active_forces,
    radial_alignment_profile,
    run_organoid_pull,
)

# A small, fast fixture so the suite runs in a reasonable time.  Softened collagen
# (3 MPa) is numerically stable at dt=0.05, so the pull test stays short.
SMALL = OrganoidConfig(
    n_fibers=90,
    n_corona_fibers=32,
    organoid_radius=32.0,
    domain_size=280.0,
    boundary_width=5.0,
    generation_attempts=15,
    dt=0.05,
    duration=24.0,
    sample_interval=4.0,
    total_pull_force=15.0,
)


class StageAScaffold(unittest.TestCase):
    def test_hex_packing_is_a_disk(self):
        c = hex_centers(SMALL.organoid_radius, SMALL.cell_spacing)
        self.assertGreaterEqual(len(c), 7)
        self.assertTrue(np.all(np.linalg.norm(c, axis=1) <= SMALL.organoid_radius + 1e-6))

    def test_cell_cell_balanced_at_rest(self):
        """Hex lattice at pitch == cell_spacing is ~force-balanced (Stage-A gate)."""
        c = hex_centers(SMALL.organoid_radius, SMALL.cell_spacing)
        f = cell_cell_forces(c, SMALL)
        # interior cells cancel exactly; edge cells feel the free boundary only.
        rest = np.linalg.norm(f, axis=1)
        self.assertLess(float(np.median(rest)), 1e-6)

    def test_cell_cell_is_cohesive(self):
        """A displaced cell feels a restoring (non-zero) force -> cohesion."""
        c = hex_centers(SMALL.organoid_radius, SMALL.cell_spacing)
        moved = c.copy()
        moved[1] = moved[1] + np.array([3.0, 0.0])
        f = cell_cell_forces(moved, SMALL)
        self.assertGreater(float(np.linalg.norm(f[1])), 1e-3)

    def test_network_builds_and_percolates(self):
        net, centers, gap_r, report = make_organoid(SMALL, seed=23)
        self.assertGreater(len(net.r), 500)
        self.assertEqual(len(net.fibers), SMALL.n_fibers)
        self.assertGreater(report["connected_fraction"], 0.5)
        # fibre-free gap: no bead inside the organoid disk
        self.assertGreaterEqual(float(np.min(np.linalg.norm(net.r, axis=1))), gap_r - 1.0)

    def test_repulsion_keeps_beads_outside_cells(self):
        net, centers, gap_r, _ = make_organoid(SMALL, seed=23)
        f = multi_cell_repulsion(net, centers)
        # beads sit outside the gap, so no cell should be penetrated at t=0
        self.assertEqual(float(np.max(np.abs(f))), 0.0)


class StageBAlignment(unittest.TestCase):
    def test_surface_cells_grip_fibers(self):
        net, centers, gap_r, _ = make_organoid(SMALL, seed=23)
        nodal, per_cell = organoid_active_forces(net, centers, SMALL.total_pull_force)
        gripping = sum(1 for p in per_cell if p)
        self.assertGreater(gripping, 0)
        self.assertGreater(float(np.max(np.abs(nodal))), 0.0)

    def test_radial_alignment_rises_under_pull(self):
        """Hongbo's minimal criterion: contractile pull -> collagen turns radial."""
        out = run_organoid_pull(SMALL, seed=23)
        order0 = out["frames"][0]["global_radial_order"]
        order1 = out["frames"][-1]["global_radial_order"]
        # net displacement occurred and alignment moved toward radial (+)
        self.assertGreater(out["frames"][-1]["rms_bead_disp"], 1e-3)
        self.assertGreater(order1, order0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
