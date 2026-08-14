"""Run and save the Generation-2 baseline gates and small sensitivity checks."""

from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from common.model import (  # noqa: E402
    CollagenConfig,
    Network,
    active_forces,
    connectivity_report,
    make_network_spec,
    run_fixed_pull,
)
from v3_two_sided_migration.model import (  # noqa: E402
    MigrationConfig,
    run_speed_ensemble,
)


def _energy_error(result: dict) -> float:
    stored = (
        result["stretch_energy"]
        + result["bend_energy"]
        + result["crosslink_energy"]
    )
    balance_rate = result["active_work_rate"] - result["drag_power"]
    integrated = float(np.trapezoid(balance_rate, result["time"]))
    delta = float(stored[-1] - stored[0])
    return abs(delta - integrated) / max(abs(delta), abs(integrated), 1e-10)


def _short_variant(cfg: CollagenConfig, **changes) -> dict:
    local = replace(
        cfg,
        duration=10.0,
        sample_interval=2.0,
        **changes,
    )
    spec = make_network_spec(local)
    result = run_fixed_pull(local, spec=spec, include_crosslinks=True)
    displacement = np.linalg.norm(
        result["positions"][-1] - result["initial_positions"], axis=1
    )
    return {
        "near_displacement_um": float(result["displacement_by_shell"][-1, 0]),
        "intermediate_displacement_um": float(result["displacement_by_shell"][-1, 1]),
        "max_displacement_um": float(np.max(displacement)),
        "p99_strain": float(np.quantile(np.abs(result["bond_strain"][-1]), 0.99)),
        "connected_fraction": float(result["connectivity"]["connected_fraction"]),
    }


def main() -> dict:
    cfg = CollagenConfig()
    spec = make_network_spec(cfg)
    network = Network(spec, cfg)
    connectivity = connectivity_report(network)
    nodal, patches, vectors = active_forces(network, np.zeros(2), 5.0, 0.0)
    free = run_fixed_pull(cfg, spec=spec, include_crosslinks=False)
    linked = run_fixed_pull(cfg, spec=spec, include_crosslinks=True)

    final_strain = np.abs(linked["bond_strain"][-1])
    intermediate_ratio = float(
        linked["displacement_by_shell"][-1, 1]
        / max(free["displacement_by_shell"][-1, 1], 1e-12)
    )
    mechanics = {
        "contact_count_frame0": len(linked["contact_points"][0]),
        "contact_weight_sum": float(sum(linked["contact_weights"][0])),
        "distributed_force_error_nN": float(
            np.linalg.norm(nodal.sum(axis=0) - np.sum(vectors, axis=0))
        ),
        "min_cell_gap_um": float(np.min(linked["min_cell_gap"])),
        "p99_abs_bond_strain": float(np.quantile(final_strain, 0.99)),
        "max_abs_bond_strain": float(np.max(final_strain)),
        "crosslinked_to_free_intermediate_displacement_ratio": intermediate_ratio,
        "overdamped_energy_balance_relative_error": _energy_error(linked),
        "near_displacement_um": float(linked["displacement_by_shell"][-1, 0]),
        "intermediate_displacement_um": float(linked["displacement_by_shell"][-1, 1]),
        "far_displacement_um": float(linked["displacement_by_shell"][-1, 2]),
    }

    migration_cfg = MigrationConfig()
    speeds = run_speed_ensemble(migration_cfg, trials=20)
    mobility = {
        "cell_line": "MDA-MB-231",
        "trials": 20,
        "median_speed_um_per_min": float(np.median(speeds)),
        "min_speed_um_per_min": float(np.min(speeds)),
        "max_speed_um_per_min": float(np.max(speeds)),
        "target_um_per_min": [0.2, 0.4],
    }

    force_sensitivity = {
        str(force): _short_variant(cfg, total_pull_force=force)
        for force in (2.5, 5.0, 10.0)
    }
    compression_sensitivity = {
        str(ratio): _short_variant(cfg, compression_ratio=ratio)
        for ratio in (0.0, 0.1, 1.0)
    }
    spacing_convergence = {
        str(spacing): _short_variant(cfg, bead_spacing=spacing)
        for spacing in (0.5, 0.75, 1.0)
    }
    domain_convergence = {
        "140": _short_variant(cfg, domain_size=140.0, n_fibers=60),
        "180": _short_variant(cfg, domain_size=180.0, n_fibers=99),
    }

    gates = {
        "v2_network": bool(
            connectivity["contact_fibers_connected"]
            and connectivity["connected_fraction"] >= cfg.required_connected_fraction
            and np.all(
                cfg.domain_size / 2
                - np.max(np.abs(network.r[network.fixed]), axis=1)
                <= cfg.boundary_width + 1e-9
            )
        ),
        "v2_mechanics": bool(
            len(linked["contact_points"][0]) > 0
            and abs(sum(linked["contact_weights"][0]) - 1.0) < 1e-10
            and mechanics["distributed_force_error_nN"] < 1e-9
            and mechanics["min_cell_gap_um"] >= -0.02
            and mechanics["p99_abs_bond_strain"] < 0.05
            and mechanics["max_abs_bond_strain"] < 0.15
            and intermediate_ratio > 1.05
            and mechanics["overdamped_energy_balance_relative_error"] < 0.35
        ),
        "v3_mobility": bool(0.2 <= np.median(speeds) <= 0.4),
    }
    summary = {
        "gates": gates,
        "network": {
            **connectivity,
            "bead_count": len(network.r),
            "bond_count": len(network.edges),
            "crosslink_count": len(network.crosslinks),
            "fixed_bead_count": int(network.fixed.sum()),
            "seed_used": network.seed_used,
        },
        "mechanics": mechanics,
        "mobility": mobility,
        "sensitivity": {
            "force_nN": force_sensitivity,
            "compression_ratio": compression_sensitivity,
            "bead_spacing_um": spacing_convergence,
            "domain_size_um": domain_convergence,
        },
        "interpretation": {
            "sensitivity_runs": "10 s directional checks, not calibration fits",
            "v4_allowed": all(gates.values()),
        },
    }
    (HERE / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    main()
