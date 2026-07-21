"""Run the paired 1D experiment and build data products for the interactive lab."""

import argparse
import csv
import json
import math
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any, Dict

import numpy as np

from collagen_model.single_protrusion import (
    SingleProtrusionConfig,
    run_lifetime_ensemble,
    run_mechanism_sweep,
    run_single_protrusion_pair,
)


def _rounded(values: np.ndarray, decimals: int = 5) -> list:
    return np.round(np.asarray(values, dtype=float), decimals).tolist()


def _finite_or_none(value: float) -> Any:
    return float(value) if math.isfinite(value) else None


def _lab_condition(result: Dict[str, object]) -> Dict[str, object]:
    return {
        "tau": result["tau"],
        "time": _rounded(result["time"], 3),
        "beads": _rounded(result["bead_positions"], 4),
        "q": _rounded(np.asarray(result["sls_q"])[:, 0], 4),
        "bound": np.asarray(result["clutch_bound"], dtype=np.uint8).tolist(),
        "clutchForce": _rounded(result["clutch_force"], 4),
        "traction": _rounded(result["total_traction"], 4),
        "velocity": _rounded(result["actin_velocity"], 4),
        "loadingRate": _rounded(result["loading_rate"], 4),
        "phase": np.asarray(result["phase"], dtype=str).tolist(),
        "episodes": result["episodes"],
    }


def _write_pair_csv(path: Path, pair: Dict[str, object]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "condition",
                "time",
                "terminal_displacement",
                "sls_q",
                "bound_clutches",
                "total_traction",
                "actin_velocity",
                "loading_rate",
                "phase",
            ]
        )
        for label in ("fast", "slow"):
            result = pair[label]
            terminal_displacement = (
                np.asarray(result["bead_positions"])[:, -1]
                - pair["config"]["fibre_rest_length"]
            )
            bound_count = np.asarray(result["clutch_bound"]).sum(axis=1)
            for row in zip(
                result["time"],
                terminal_displacement,
                np.asarray(result["sls_q"])[:, 0],
                bound_count,
                result["total_traction"],
                result["actin_velocity"],
                result["loading_rate"],
                result["phase"],
            ):
                writer.writerow([label, *row])


def _write_ensemble_csv(path: Path, ensemble: Dict[str, object]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "trial",
                "fast_loading_rate",
                "slow_loading_rate",
                "fast_completed_clusters",
                "slow_completed_clusters",
            ]
        )
        fast = ensemble["fast"]
        slow = ensemble["slow"]
        for trial in range(ensemble["trials"]):
            writer.writerow(
                [
                    trial,
                    fast["trial_loading_rate"][trial],
                    slow["trial_loading_rate"][trial],
                    fast["trial_completed_clusters"][trial],
                    slow["trial_completed_clusters"][trial],
                ]
            )


def _write_episode_csv(path: Path, ensemble: Dict[str, object]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["condition", "duration", "failure_observed"])
        for label in ("fast", "slow"):
            for duration, observed in zip(
                ensemble[label]["episode_durations"],
                ensemble[label]["episode_observed"],
            ):
                writer.writerow([label, duration, int(observed)])


def _write_sweep_csv(path: Path, records: list) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "de",
            "k0_chain_over_kc",
            "kinf_over_k0",
            "episodes",
            "observed_failures",
            "km_median_lifetime",
            "median_loading_rate",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            row = dict(record)
            if not math.isfinite(row["km_median_lifetime"]):
                row["km_median_lifetime"] = ""
            writer.writerow(row)


def _summary(ensemble: Dict[str, object]) -> Dict[str, object]:
    result = {"trials": ensemble["trials"], "seed": ensemble["seed"]}
    for label in ("fast", "slow"):
        condition = ensemble[label]
        result[label] = {
            "episodes": int(len(condition["episode_durations"])),
            "observed_failures": int(np.count_nonzero(condition["episode_observed"])),
            "km_median_lifetime": _finite_or_none(condition["km_median_lifetime"]),
            "median_lifetime_reached": math.isfinite(condition["km_median_lifetime"]),
            "median_positive_loading_rate": float(condition["median_loading_rate"]),
        }
    result["acceptance"] = {
        "fast_has_lower_loading_rate": (
            result["fast"]["median_positive_loading_rate"]
            < result["slow"]["median_positive_loading_rate"]
        ),
        "fast_has_longer_cluster_lifetime": (
            result["fast"]["median_lifetime_reached"]
            and result["slow"]["median_lifetime_reached"]
            and result["fast"]["km_median_lifetime"]
            > result["slow"]["km_median_lifetime"]
        ),
    }
    return result


def _write_lab(path: Path, template_path: Path, pair: Dict[str, object], summary: Dict[str, object]) -> None:
    payload = {
        "config": pair["config"],
        "seed": pair["seed"],
        "summary": summary,
        "fast": _lab_condition(pair["fast"]),
        "slow": _lab_condition(pair["slow"]),
    }
    template = template_path.read_text(encoding="utf-8")
    marker = "__SINGLE_PROTRUSION_DATA__"
    if template.count(marker) != 1:
        raise ValueError(f"template must contain exactly one {marker} marker")
    fragment = template.replace(marker, json.dumps(payload, separators=(",", ":")))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(fragment, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/single_protrusion"))
    parser.add_argument("--lab", type=Path)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--ensemble-seed", type=int, default=1000)
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--sweep-trials", type=int, default=10)
    parser.add_argument("--sweep-only", action="store_true")
    parser.add_argument("--duration", type=float, default=30.0)
    parser.add_argument("--dt", type=float, default=1.0e-3)
    parser.add_argument("--sample-interval", type=float, default=0.02)
    args = parser.parse_args()

    config = replace(
        SingleProtrusionConfig(),
        duration=args.duration,
        dt=args.dt,
        sample_interval=args.sample_interval,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    if args.sweep_only:
        sweep = run_mechanism_sweep(
            config,
            trials=args.sweep_trials,
            seed=args.ensemble_seed + args.trials,
            workers=args.workers,
        )
        _write_sweep_csv(args.output / "mechanism_sweep.csv", sweep)
        print(json.dumps({"points": len(sweep), "trials_per_point": args.sweep_trials}, indent=2))
        return

    pair = run_single_protrusion_pair(config, args.seed)
    ensemble = run_lifetime_ensemble(
        config,
        trials=args.trials,
        seed=args.ensemble_seed,
        workers=args.workers,
    )
    sweep = run_mechanism_sweep(
        config,
        trials=args.sweep_trials,
        seed=args.ensemble_seed + args.trials,
        workers=args.workers,
    )
    summary = _summary(ensemble)

    _write_pair_csv(args.output / "paired_timeseries.csv", pair)
    _write_ensemble_csv(args.output / "ensemble_trials.csv", ensemble)
    _write_episode_csv(args.output / "cluster_episodes.csv", ensemble)
    _write_sweep_csv(args.output / "mechanism_sweep.csv", sweep)
    (args.output / "ensemble_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    (args.output / "run_metadata.json").write_text(
        json.dumps(
            {
                "config": asdict(config),
                "representative_seed": args.seed,
                "ensemble_seed": args.ensemble_seed,
                "trials": args.trials,
                "sweep_trials_per_point": args.sweep_trials,
                "sweep_seed": args.ensemble_seed + args.trials,
                "interpretation": (
                    "Single fixed-cell protrusion; cluster lifetimes are right-censored "
                    "at the simulation duration."
                ),
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    if args.lab is not None:
        template = Path(__file__).parent / "visualization" / "single-protrusion-lab.template.html"
        _write_lab(args.lab, template, pair, summary)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
