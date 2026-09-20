"""Predeclared numerical and coordinate-bias checks for G4 v4F.

The validation run is intentionally separate from the demonstration builder so
an unfinished sensitivity cannot silently change a published trajectory.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
import math
from pathlib import Path

import numpy as np

from .model import G4V4FConfig, make_matched_specs, run_v4f_condition


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "validation_summary.json"
VALIDATION_JS = HERE.parent.parent / "docs" / "g4-v4f-validation.js"


def _run(cfg, *, duration, matched=None):
    matched = make_matched_specs(cfg) if matched is None else matched
    return run_v4f_condition(
        cfg, "combined", cfg.seed, matched=matched,
        duration=duration, store_positions=False,
    )


def run_validation(*, duration=600.0):
    base = G4V4FConfig()
    base_matched = make_matched_specs(base)

    # The E2 module must be dormant before a front can mature.  This verifies
    # that the corrected baseline, initial contacts and random stream are
    # genuinely shared rather than merely described as shared.
    prefront_duration = min(60.0, 0.5 * base.front_maturation)
    prefront_matched = base_matched
    prefront_e0 = run_v4f_condition(
        base, "d_control", matched=prefront_matched,
        duration=prefront_duration, store_positions=False,
    )
    prefront_e2 = run_v4f_condition(
        base, "protrusion_guidance", matched=prefront_matched,
        duration=prefront_duration, store_positions=False,
    )
    prefront_match = bool(np.allclose(prefront_e0["centers"], prefront_e2["centers"]))
    baseline_mechanics = {
        "prefront_duration": prefront_duration,
        "E0_E2_centers_identical": prefront_match,
        "E0_front_events": len(prefront_e0["front_events"]),
        "force_pair_residual": prefront_e0["metrics"]["max_force_balance_error"],
        "max_relocation_gap": prefront_e0["metrics"]["max_relocation_gap"],
        "contact_search_limit": prefront_e0["metrics"]["contact_search_limit"],
    }

    timestep = []
    for dt in (0.05, 0.025):
        result = _run(replace(base, dt=dt), duration=duration, matched=base_matched)
        timestep.append({
            "dt": dt,
            "D_cue": result["metrics"]["D_cue"],
            "persistence": result["metrics"]["persistence"],
            "front_event_order": [event["kind"] for event in result["front_events"][:4]],
        })
    timestep_direction = np.sign(timestep[0]["D_cue"]) == np.sign(timestep[1]["D_cue"])
    timestep_events = timestep[0]["front_event_order"] == timestep[1]["front_event_order"]

    radius = []
    for value in base.radius_check_values:
        result = _run(
            replace(base, cell_radius=value),
            duration=duration,
            matched=base_matched if value == base.cell_radius else None,
        )
        radius.append({
            "radius": value,
            "D_cue": result["metrics"]["D_cue"],
            "D_cue_over_R": result["metrics"]["D_cue"] / value,
            "persistence": result["metrics"]["persistence"],
        })

    guidance_sensitivity = {
        "tract_angular_sd": [],
        "front_maturation": [],
        "beta_alignment_and_memory": [],
        "contact_search_interval": [],
    }
    for degrees in (15.0, 30.0):
        radians = math.radians(degrees)
        result = _run(
            replace(base, radial_tract_angular_sd=radians),
            duration=duration,
            matched=(
                base_matched
                if abs(radians - base.radial_tract_angular_sd) < 1e-12
                else None
            ),
        )
        guidance_sensitivity["tract_angular_sd"].append({
            "value_deg": degrees,
            "D_cue": result["metrics"]["D_cue"],
            "persistence": result["metrics"]["persistence"],
        })
    for seconds in (60.0, 120.0, 180.0):
        result = _run(
            replace(base, front_maturation=seconds),
            duration=duration,
            matched=base_matched,
        )
        guidance_sensitivity["front_maturation"].append({
            "value_seconds": seconds,
            "D_cue": result["metrics"]["D_cue"],
            "persistence": result["metrics"]["persistence"],
        })
    for gain in (0.0, 1.0, 2.0, 4.0):
        result = _run(
            replace(base, beta_alignment=gain, beta_memory=gain),
            duration=duration,
            matched=base_matched,
        )
        guidance_sensitivity["beta_alignment_and_memory"].append({
            "value": gain,
            "D_cue": result["metrics"]["D_cue"],
            "persistence": result["metrics"]["persistence"],
        })
    # The search clock is a shared-baseline numerical assumption, not a
    # guidance factor.  All conditions keep the same selected value.  This
    # sensitivity checks whether the 5 s default changes the qualitative
    # combined-condition result relative to more frequent searches.
    for seconds in (0.5, 2.0, 5.0):
        result = _run(
            replace(base, contact_search_interval=seconds),
            duration=duration,
            matched=base_matched,
        )
        guidance_sensitivity["contact_search_interval"].append({
            "value_seconds": seconds,
            "D_cue": result["metrics"]["D_cue"],
            "persistence": result["metrics"]["persistence"],
            "successful_relocations": result["metrics"]["successful_relocations"],
        })

    rotations = []
    for angle in (0.0, math.pi / 2.0, math.pi):
        result = _run(
            replace(base, radial_cue_angle=angle),
            duration=duration,
            matched=base_matched if abs(angle - base.radial_cue_angle) < 1e-12 else None,
        )
        displacement = result["centers"][-1] - result["centers"][0]
        axis = np.asarray([math.cos(angle), math.sin(angle)])
        rotations.append({
            "cue_angle_deg": math.degrees(angle),
            "projected_displacement": float(displacement @ axis),
            "world_displacement": displacement.tolist(),
        })

    # A moderate field increase preserves source-fibre areal density.  This is
    # the expensive outer-boundary audit, not a fitted model variation.
    domain = []
    density = base.n_fibers / base.domain_size**2
    for size in (180.0, 240.0):
        cfg = replace(base, domain_size=size, n_fibers=int(round(density * size**2)))
        result = _run(
            cfg,
            duration=duration,
            matched=base_matched if size == base.domain_size else None,
        )
        near = float(result["delta_radial_order_by_shell"][-1][0])
        domain.append({
            "domain_size": size,
            "n_fibers": cfg.n_fibers,
            "near_delta_radial_order": near,
            "D_cue": result["metrics"]["D_cue"],
        })
    relative = abs(domain[1]["near_delta_radial_order"] - domain[0]["near_delta_radial_order"]) / max(
        abs(domain[0]["near_delta_radial_order"]),
        abs(domain[1]["near_delta_radial_order"]), 1e-12,
    )

    summary = {
        "duration_seconds": duration,
        "status": "full predeclared audit" if duration >= base.migration_duration else "short audit; full-clock gates remain exploratory",
        "baseline_mechanics": baseline_mechanics,
        "timestep": timestep,
        "radius": radius,
        "guidance_sensitivity": guidance_sensitivity,
        "cue_rotation": rotations,
        "domain": domain,
        "gates": {
            "prefront_E0_E2_identical": prefront_match,
            "baseline_has_no_front_events": bool(len(prefront_e0["front_events"]) == 0),
            "pair_force_residual_below_tolerance": bool(
                prefront_e0["metrics"]["max_force_balance_error"] < 1e-12
            ),
            "dynamic_contacts_within_reach": bool(
                prefront_e0["metrics"]["max_relocation_gap"]
                <= prefront_e0["metrics"]["contact_search_limit"] + 1e-9
            ),
            "timestep_preserves_direction": bool(timestep_direction),
            "timestep_preserves_front_event_order": bool(timestep_events),
            "cue_rotation_all_positive_projection": bool(all(row["projected_displacement"] > 0 for row in rotations)),
            "domain_near_cell_change_below_five_percent": bool(relative < 0.05),
        },
        "domain_relative_change": float(relative),
    }
    OUTPUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    VALIDATION_JS.write_text(
        "window.G4V4F_VALIDATION = "
        + json.dumps(summary, separators=(",", ":"))
        + ";\n",
        encoding="utf-8",
    )
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--duration", type=float, default=600.0)
    args = parser.parse_args()
    print(json.dumps(run_validation(duration=args.duration), indent=2))


if __name__ == "__main__":
    main()
