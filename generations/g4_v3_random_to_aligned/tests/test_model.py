"""Small deterministic checks for G4 v3; no long production run required."""

from dataclasses import replace

import numpy as np

from generations.g4_v3_random_to_aligned.model import (
    G4V3Config,
    _inner_boundary_beads,
    build_random_network,
    make_random_void_spec,
    run_fixed_cell_passive_check,
    run_fixed_material_traction,
    run_prescribed_contraction,
)


def small_config(**changes):
    base = G4V3Config(
        n_fibers=120,
        benchmark_ramp_time=0.10,
        benchmark_relax_time=0.20,
        benchmark_sample_interval=0.10,
        traction_duration=0.20,
        traction_sample_interval=0.10,
        force_ramp_time=0.10,
        dt=0.01,
    )
    return replace(base, **changes)


def test_random_geometry_is_not_handcrafted_radial():
    cfg = small_config()
    network, _, report = build_random_network(cfg)
    assert report["initial_global_nematic_order"] < 0.15
    assert report["contact_fragments"] >= 4
    radius = np.linalg.norm(network.r0, axis=1)
    assert np.min(radius) >= cfg.cell_radius + cfg.cell_clearance - 1e-8


def test_probability_adds_links_without_reshuffling():
    cfg = small_config()
    spec = make_random_void_spec(cfg)
    low, _, _ = build_random_network(cfg, spec=spec, crosslink_probability=0.20)
    high, _, _ = build_random_network(cfg, spec=spec, crosslink_probability=0.60)

    def identity(link):
        return (link.edge_a, link.alpha_a, link.edge_b, link.alpha_b)

    assert {identity(link) for link in low.crosslinks} <= {
        identity(link) for link in high.crosslinks
    }
    assert np.array_equal(low.r0, high.r0)


def test_prescribed_boundary_reaches_target_without_stochastic_terms():
    cfg = small_config()
    spec = make_random_void_spec(cfg)
    result = run_prescribed_contraction(
        cfg, contraction=0.01, spec=spec, crosslink_probability=1.0
    )
    inner = _inner_boundary_beads(result["network"])
    expected = 0.99 * result["network"].r0[inner]
    assert np.allclose(result["network"].r[inner], expected)
    assert np.isfinite(result["final_stretch_to_bend_ratio"])
    assert result["frames"][-1]["required_inward_boundary_force"] > 0.0


def test_fixed_passive_cell_does_not_create_alignment():
    result = run_fixed_cell_passive_check(small_config())
    assert result["initial_max_repulsion"] < 1e-8
    assert result["max_displacement"] < 1e-8
    assert np.max(np.abs(result["delta_radial_order_by_shell"])) < 1e-8


def test_hybrid_contact_weights_and_reaction_ledger():
    result = run_fixed_material_traction(small_config(), total_force=2.0)
    assert len(result["patches"]) >= 4
    assert np.isclose(sum(patch.weight for patch in result["patches"]), 1.0)
    final = result["frames"][-1]
    assert np.isclose(final["distributed_force_magnitude"], 2.0)
    assert final["action_reaction_error"] < 1e-12
