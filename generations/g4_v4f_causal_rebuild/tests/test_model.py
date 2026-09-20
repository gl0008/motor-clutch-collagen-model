from dataclasses import replace

import numpy as np

import generations.g4_v4f_causal_rebuild.model as model


def test_stage_configuration_is_cumulative_and_physics_starts_last():
    stages = [model.config_for_stage(stage) for stage in model.MICROSTAGES]
    assert not stages[0].independent_rng
    assert stages[1].independent_rng
    assert stages[2].same_time_level_update
    assert stages[3].all_site_front_selection
    assert stages[4].matched_initial_contacts
    assert stages[5].shared_probing_reach
    for cfg in stages:
        assert not cfg.steric_reaction_to_cell
        assert not cfg.dynamic_contact_search
        assert not cfg.memoryless_controls
    full = model.config_for_stage("F_FULL")
    assert full.steric_reaction_to_cell
    assert full.dynamic_contact_search
    assert full.memoryless_controls


def test_factorial_arms_change_only_the_three_declared_switches():
    baseline = model.config_for_stage("N5_SHARED_REACH")
    ignored = {
        "steric_reaction_to_cell", "dynamic_contact_search", "memoryless_controls"
    }
    reference = {
        key: value for key, value in vars(baseline).items() if key not in ignored
    }
    for s in (False, True):
        for c in (False, True):
            for m in (False, True):
                arm = model.mechanics_arm(s, c, m)
                observed = {
                    key: value for key, value in vars(arm).items() if key not in ignored
                }
                assert observed == reference
                assert arm.steric_reaction_to_cell is s
                assert arm.dynamic_contact_search is c
                assert arm.memoryless_controls is m


def test_factorial_math_recovers_main_effect_and_interaction():
    values = {}
    for s in (0, 1):
        for c in (0, 1):
            for m in (0, 1):
                values[(s, c, m)] = np.asarray([2 * s + 3 * c + 5 * s * c + 7 * m])
    assert np.allclose(model._factorial_effect(values, 2), [7.0])
    assert np.allclose(model._factorial_interaction(values, 0, 1), [5.0])
    assert np.allclose(model._factorial_interaction(values, 0, 2), [0.0])


def test_exact_e_replay_gate_short_clock():
    report = model.validate_exact_e_replay(duration=0.1)
    assert report["passed"]
    for comparison in report["comparisons"].values():
        assert max(comparison.values()) < report["tolerance"]


def test_steric_pair_is_equal_and_opposite():
    r = np.asarray([[9.0, 0.0], [-9.5, 0.0], [0.0, 12.0]])
    values = model._steric_force_pair(r, np.zeros(2), 10.0, 0.5, 4.0)
    assert np.allclose(np.asarray(values[:2]) + np.asarray(values[2:]), 0.0)


def test_independent_counter_stream_breaks_old_seed_item_alias():
    seed = model.G4CausalConfig().counter_seed
    historical_a = model._uniform(seed + 41, 1, 0, 10, False)
    historical_b = model._uniform(seed + 42, 1, 0, 9, False)
    repaired_a = model._uniform(seed + 41, 1, 0, 10, True)
    repaired_b = model._uniform(seed + 42, 1, 0, 9, True)
    assert historical_a == historical_b
    assert repaired_a != repaired_b
