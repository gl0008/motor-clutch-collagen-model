"""G5-R0 mechanical-consistency gate tests (mirror the G4D force-consistency invariants).

Run: python -m pytest generations/g5_organoid/tests/test_r0_consistency.py

These assert the R0 gate: force-pair (Newton 3rd law) residual ~0, symmetric reaction
=> no net translation, passive (no clutch) => ECM/cells inert, event-driven relocation
resets clutch state, and rigid co-moving translation is not misread as clutch loading.
All own-sim output is personal testing, NOT confirmed findings.
"""

from pathlib import Path
import sys
import types
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from generations.g5_organoid.model import (  # noqa: E402
    make_organoid,
    organoid_clutch_patches,
    cell_candidate_fibers,
    _new_clutch_state,
)
from generations.g5_organoid.consistency import (  # noqa: E402
    r0_config,
    run_r0_invasion,
    per_cell_clutch_reaction,
    _reset_site_state,
    _relative_substrate_speeds_mc,
)

# Small, fast fixture with the R0 baseline (shared clutch, G4D drag/speed) applied.
SMALL = r0_config(
    n_fibers=90,
    n_corona_fibers=32,
    organoid_radius=32.0,
    domain_size=280.0,
    boundary_width=5.0,
    generation_attempts=15,
    dt=0.05,
    duration=20.0,
    sample_interval=4.0,
    total_pull_force=15.0,
)


class ForcePair(unittest.TestCase):
    def test_force_pair_residual_is_numerically_zero(self):
        """Per-cell reaction == -(clutch traction applied to ECM), to ~machine precision."""
        out = run_r0_invasion(SMALL, seed=23)
        self.assertLess(out["max_force_pair_residual"], 1e-9)
        self.assertTrue(np.all(np.isfinite(out["centers_final"])))

    def test_symmetric_reaction_does_not_translate(self):
        """Equal-and-opposite site tractions on one cell -> zero net reaction (G4D)."""
        # two sectors of one cell (n_sec=2, M=1), inward normals +x and -x, equal force
        patches = [types.SimpleNamespace(normal_in=np.array([1.0, 0.0])),
                   types.SimpleNamespace(normal_in=np.array([-1.0, 0.0]))]
        reaction = per_cell_clutch_reaction(patches, np.array([2.5, 2.5]), n_sec=2, n_cells=1)
        self.assertLess(float(np.linalg.norm(reaction[0])), 1e-12)

    def test_symmetric_organoid_has_no_world_axis_drift(self):
        """A symmetric organoid's centroid does not systematically translate (soft check)."""
        out = run_r0_invasion(SMALL, seed=23)
        drift = float(np.linalg.norm(out["centers_final"].mean(0) - out["centers0"].mean(0)))
        self.assertLess(drift, 0.5)   # um; motion is speed-capped + symmetric -> tiny centroid drift


class PassiveControl(unittest.TestCase):
    def test_passive_no_clutch_leaves_ecm_inert(self):
        """clutch_on_rate=0 -> no traction -> the ECM does not remodel."""
        cfg = r0_config(SMALL, clutch_on_rate=0.0)
        out = run_r0_invasion(cfg, seed=23)
        moved = float(np.max(np.abs(out["final_positions"] - out["initial_positions"])))
        self.assertLess(moved, 1e-6)
        self.assertLess(out["frames"][-1]["max_cell_disp"], 1e-6)
        self.assertEqual(out["cumulative_site_failures"], 0)


class Relocation(unittest.TestCase):
    def test_reset_clears_site_history(self):
        """A relocated site carries no loading history to the new fibre."""
        state = _new_clutch_state(6, SMALL)
        state.bound[2, :] = True
        state.extension[2, :] = 5.0
        state.site_extension[2] = 5.0
        _reset_site_state(state, 2)
        self.assertFalse(bool(state.bound[2].any()))
        self.assertEqual(float(np.max(np.abs(state.extension[2]))), 0.0)
        self.assertEqual(float(state.site_extension[2]), 0.0)

    def test_every_site_failure_triggers_one_relocation(self):
        """Event-driven rule: relocations == site failures (no periodic teleporting)."""
        cfg = r0_config(SMALL, duration=60.0, bell_force=0.5,
                        clutch_off_rate0=0.05, total_pull_force=20.0)
        out = run_r0_invasion(cfg, seed=23)
        self.assertEqual(out["n_relocations"], out["cumulative_site_failures"])
        self.assertGreater(out["cumulative_site_failures"], 0)   # the failure path is exercised


class RelativeSubstrate(unittest.TestCase):
    def test_rigid_co_moving_translation_is_not_clutch_loading(self):
        """v_material == v_cell -> relative substrate speed ~ 0 (G4D cell-frame actin)."""
        net, centers, _gap, _rep = make_organoid(SMALL, seed=23)
        reach = SMALL.cell_radius + SMALL.contact_width + 2.0
        cand = cell_candidate_fibers(net, centers, reach)
        patches, _sc = organoid_clutch_patches(net, centers, SMALL, cand)
        s = next(i for i, p in enumerate(patches) if p is not None)
        v = np.array([0.013, -0.021])
        bead_vel = np.zeros_like(net.r)
        i, j = net.edges[patches[s].edge]
        bead_vel[i] = v
        bead_vel[j] = v
        site_centers = np.zeros((len(patches), 2))          # inward direction is irrelevant here
        v_cell_per_site = np.tile(v, (len(patches), 1))
        rel = _relative_substrate_speeds_mc(net, bead_vel, patches, site_centers, v_cell_per_site)
        self.assertLess(abs(float(rel[s])), 1e-12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
