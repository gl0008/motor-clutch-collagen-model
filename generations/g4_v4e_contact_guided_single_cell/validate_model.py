"""Predeclared numerical and coordinate-bias checks for G4 v4E.

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

from .model import G4V4EConfig, make_matched_specs, run_v4e_condition


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "validation_summary.json"
VALIDATION_JS = HERE.parent.parent / "docs" / "g4-v4e-validation.js"


def _run(cfg, *, duration):
    matched = make_matched_specs(cfg)
    return run_v4e_condition(
        cfg, "combined", cfg.seed, matched=matched,
        duration=duration, store_positions=False,
    )


def run_validation(*, duration=600.0):
    base = G4V4EConfig()

    timestep = []
    for dt in (0.05, 0.025):
        result = _run(replace(base, dt=dt), duration=duration)
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
        result = _run(replace(base, cell_radius=value), duration=duration)
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
    }
    for degrees in (15.0, 30.0):
        result = _run(
            replace(base, radial_tract_angular_sd=math.radians(degrees)),
            duration=duration,
        )
        guidance_sensitivity["tract_angular_sd"].append({
            "value_deg": degrees,
            "D_cue": result["metrics"]["D_cue"],
            "persistence": result["metrics"]["persistence"],
        })
    for seconds in (60.0, 120.0, 180.0):
        result = _run(replace(base, front_maturation=seconds), duration=duration)
        guidance_sensitivity["front_maturation"].append({
            "value_seconds": seconds,
            "D_cue": result["metrics"]["D_cue"],
            "persistence": result["metrics"]["persistence"],
        })
    for gain in (0.0, 1.0, 2.0, 4.0):
        result = _run(
            replace(base, beta_alignment=gain, beta_memory=gain),
            duration=duration,
        )
        guidance_sensitivity["beta_alignment_and_memory"].append({
            "value": gain,
            "D_cue": result["metrics"]["D_cue"],
            "persistence": result["metrics"]["persistence"],
        })

    rotations = []
    for angle in (0.0, math.pi / 2.0, math.pi):
        result = _run(replace(base, radial_cue_angle=angle), duration=duration)
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
        result = _run(cfg, duration=duration)
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
        "timestep": timestep,
        "radius": radius,
        "guidance_sensitivity": guidance_sensitivity,
        "cue_rotation": rotations,
        "domain": domain,
        "gates": {
            "timestep_preserves_direction": bool(timestep_direction),
            "timestep_preserves_front_event_order": bool(timestep_events),
            "cue_rotation_all_positive_projection": bool(all(row["projected_displacement"] > 0 for row in rotations)),
            "domain_near_cell_change_below_five_percent": bool(relative < 0.05),
        },
        "domain_relative_change": float(relative),
    }
    OUTPUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    VALIDATION_JS.write_text(
        "window.G4V4E_VALIDATION = "
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
