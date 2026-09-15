from dataclasses import replace
import math

import numpy as np

import generations.g4_v4e_contact_guided_single_cell.model as model


def small_config(**changes):
    cfg = model.G4V4EConfig(
        n_fibers=60,
        geometry_attempt_factor=80,
        migration_duration=6.0,
        migration_sample_interval=1.0,
        telemetry_interval=1.0,
        event_interval=0.5,
        event_duration=2.0,
        front_maturation=1.0,
        front_loss_time=1.0,
        radial_tract_outer_surface=70.0,
        radial_cue_half_width=math.pi / 2.0,
    )
    return replace(cfg, **changes)


def test_radial_tract_preserves_bead_count_centroids_and_lengths():
    cfg = small_config()
    base = model.make_random_void_spec(cfg)
    cue, selected, report = model.make_radial_tract_spec(cfg, base)
    assert len(cue.positions) == len(base.positions)
    assert len(cue.fibers) == len(base.fibers)
    assert report["max_contour_length_change"] < 1e-10
    for fid in selected:
        ids = np.asarray(base.fibers[fid])
        assert np.allclose(base.positions[ids].mean(axis=0), cue.positions[ids].mean(axis=0))


def test_probability_is_normalized_and_front_excludes_rear():
    cfg = small_config()
    matched = model.make_matched_specs(cfg)
    network = matched["random_network"]
    candidates = model.probing_contact_patches(network, np.zeros(2), cfg.probing_reach)
    probability, components = model.protrusion_probabilities(
        network, candidates, np.zeros(2), np.asarray([1.0, 0.0]), cfg
    )
    assert np.isclose(probability.sum(), 1.0)
    assert np.all(probability >= 0.0)
    for patch, value in zip(candidates, probability):
        point = model.patch_point(network, patch)
        if point[0] < 0.0:
            assert value == 0.0
    assert components.shape == (len(candidates), 5)


def test_all_four_additive_conditions_run_with_constant_radius():
    cfg = small_config()
    matched = model.make_matched_specs(cfg)
    for mode in model.MODES:
        result = model.run_v4e_condition(cfg, mode, matched=matched)
        assert result["metrics"]["cell_radius_constant"]
        assert result["metrics"]["max_force_balance_error"] < 1e-12
        assert np.all(np.isfinite(result["positions"]))


def test_guided_condition_establishes_front_without_world_axis_parameter():
    cfg = small_config()
    matched = model.make_matched_specs(cfg)
    result = model.run_v4e_condition(
        cfg, "protrusion_guidance", matched=matched, duration=6.0
    )
    assert result["metrics"]["front_establishments"] >= 1
    direction = result["front_events"][0]["direction"]
    assert np.isclose(np.linalg.norm(direction), 1.0)


def test_zero_guidance_gains_reduce_selection_equation_to_distance_weight():
    cfg = small_config(beta_alignment=0.0, beta_memory=0.0)
    matched = model.make_matched_specs(cfg)
    network = matched["random_network"]
    candidates = model.probing_contact_patches(network, np.zeros(2), cfg.probing_reach)
    probability, components = model.protrusion_probabilities(
        network, candidates, np.zeros(2), np.zeros(2), cfg
    )
    expected = components[:, 0] / components[:, 0].sum()
    assert np.allclose(probability, expected)


def test_cue_rotation_rotates_the_constructed_tract():
    cfg0 = small_config(seed=44, radial_cue_angle=0.0)
    cfg90 = replace(cfg0, radial_cue_angle=math.pi / 2.0)
    base = model.make_random_void_spec(cfg0)
    cue0, selected0, _ = model.make_radial_tract_spec(cfg0, base)
    cue90, selected90, _ = model.make_radial_tract_spec(cfg90, base)
    assert selected0 or selected90
    if selected0:
        centroids = [cue0.positions[np.asarray(cue0.fibers[f])].mean(axis=0) for f in selected0]
        assert np.mean([c[0] for c in centroids]) > 0.0
    if selected90:
        centroids = [cue90.positions[np.asarray(cue90.fibers[f])].mean(axis=0) for f in selected90]
        assert np.mean([c[1] for c in centroids]) > 0.0
