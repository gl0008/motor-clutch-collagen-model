"""Run the v4E-to-v4F causal reconstruction and save compact JSON evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .model import (
    G4CausalConfig,
    run_mechanics_factorial,
    run_microstage_ladder,
    validate_exact_e_replay,
)


HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE / "causal_results.json")
    parser.add_argument(
        "--docs-output",
        type=Path,
        default=REPO / "docs" / "g4-v4f-causal-results.json",
    )
    args = parser.parse_args()

    cfg = G4CausalConfig()
    if args.quick:
        replay_duration = 2.0
        ladder_duration = 2.0
        factorial_duration = 2.0
        seeds = (41, 42)
    else:
        replay_duration = 30.0
        ladder_duration = 300.0
        factorial_duration = 1800.0
        seeds = tuple(range(41, 61))

    payload = {
        "analysis_level": "quick smoke test" if args.quick else "predeclared inference",
        "exact_e_replay": validate_exact_e_replay(duration=replay_duration),
        "microstage_ladder": run_microstage_ladder(
            cfg, duration=ladder_duration
        ),
        "mechanics_factorial": run_mechanics_factorial(
            cfg, seeds=seeds, duration=factorial_duration
        ),
    }
    encoded = json.dumps(_plain(payload), indent=2) + "\n"
    args.output.write_text(encoded)
    args.docs_output.write_text(encoded)
    print(args.output)
    print(args.docs_output)


if __name__ == "__main__":
    main()
