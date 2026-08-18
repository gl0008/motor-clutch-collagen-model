from dataclasses import replace

import numpy as np
import pytest

from g3.config import G3Config
from g3.fixtures import build_fixture, rotated_fixture


def test_g3_defaults_are_si_and_explicitly_stable():
    cfg = G3Config()
    cfg.validate()
    assert cfg.clutch_stiffness == pytest.approx(5.0e-3)  # 5 pN/nm
    assert cfg.cell_drag == pytest.approx(0.3)            # 300 nN s/um
    assert cfg.rotational_drag == pytest.approx(cfg.cell_drag * cfg.cell_radius**2)
    assert cfg.dt / cfg.ecm_substeps < cfg.bead_drag / (2.0 * cfg.kappa_s_f)


def test_invalid_explicit_timestep_is_rejected():
    with pytest.raises(ValueError, match="stability"):
        replace(G3Config(), dt=1.0).validate()


def test_contact_stiffness_is_positive_and_included_in_stability_bound():
    with pytest.raises(ValueError, match="contact_stiffness"):
        replace(G3Config(), contact_stiffness=0.0).validate()
    with pytest.raises(ValueError, match="stability"):
        replace(G3Config(), contact_stiffness=0.1).validate()


@pytest.mark.parametrize("name,n_fibres", [
    ("empty", 0),
    ("single_fibre", 1),
    ("balanced_8", 8),
    ("isotropic_random_8", 8),
    ("aligned_8", 8),
    ("aligned_8_rotated_30", 8),
    ("asymmetric_torque", 3),
])
def test_minimal_fixture_shapes_and_far_end_anchors(name, n_fibres):
    cfg = G3Config()
    fixture = build_fixture(name, cfg, np.random.default_rng(3))
    assert fixture.network.n_fibers == n_fibres
    assert fixture.fixed_mask.sum() == 2 * n_fibres
    if n_fibres:
        assert fixture.network.n_beads == n_fibres * cfg.beads_per_fibre
        assert np.all(np.linalg.norm(fixture.initial_positions, axis=1) >= cfg.cell_radius)


def test_rotated_fixture_preserves_topology_and_distances():
    cfg = G3Config()
    fixture = build_fixture("aligned_8", cfg, np.random.default_rng(0))
    rotated = rotated_fixture(fixture, np.radians(30.0))
    assert np.array_equal(rotated.network.fiber_bonds, fixture.network.fiber_bonds)
    assert np.allclose(np.linalg.norm(rotated.initial_positions, axis=1),
                       np.linalg.norm(fixture.initial_positions, axis=1))
    assert rotated.director_angle == pytest.approx(np.radians(30.0))
