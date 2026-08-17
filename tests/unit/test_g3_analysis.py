import numpy as np
import pytest

from g3.analysis import FOI_RANDOM, displacement_statistics, fibre_orientation_index, nam_plasticity_index
from g3.config import G3Config
from g3.fixtures import build_fixture


def test_balanced_radial_fixture_has_unit_foi():
    cfg = G3Config()
    fixture = build_fixture("balanced_8", cfg, np.random.default_rng(0))
    assert fibre_orientation_index(fixture.initial_positions, fixture, np.zeros(2)) == pytest.approx(1.0)


def test_nam_kappa_is_zero_for_random_post_state_and_one_for_retained_state():
    pre = FOI_RANDOM + 0.2
    assert nam_plasticity_index(pre, FOI_RANDOM) == pytest.approx(0.0)
    assert nam_plasticity_index(pre, pre) == pytest.approx(1.0)


def test_displacement_participation_distinguishes_local_from_broad_motion():
    initial = np.zeros((10, 2))
    local = initial.copy()
    local[0, 0] = 1.0
    broad = initial.copy()
    broad[:, 0] = 1.0
    assert displacement_statistics(initial, local)["participation_ratio"] == pytest.approx(0.1)
    assert displacement_statistics(initial, broad)["participation_ratio"] == pytest.approx(1.0)
