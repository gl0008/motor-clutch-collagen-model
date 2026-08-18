"""Command-line entry point for reproducible G3A/G3B/G3C runs."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import yaml

from .config import DEFAULT_CONFIG_PATH, G3Config
from .simulation import run_g3, run_load_unload
from .visualization import make_stage_animation, make_summary_figure


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _git_commit():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def save_run(result, output: Path, make_gif=True, make_mp4=False):
    output.mkdir(parents=True, exist_ok=True)
    with open(output / "resolved_config.yaml", "w", encoding="utf-8") as handle:
        yaml.safe_dump({"g3": result.config.to_dict()}, handle, sort_keys=False)
    with open(output / "metrics.json", "w", encoding="utf-8") as handle:
        json.dump(_json_ready(result.summary), handle, indent=2, ensure_ascii=False)
    manifest = {
        "stage": result.stage,
        "fixture": result.fixture_name,
        "seed": result.seed,
        "git_commit": _git_commit(),
        "status": result.status,
        "interpretation": result.summary["interpretation"],
    }
    with open(output / "manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
    np.savez_compressed(output / "traces.npz", **result.traces)
    np.savez_compressed(
        output / "frames.npz",
        time=np.asarray([snapshot.time for snapshot in result.snapshots]),
        positions=np.asarray([snapshot.positions for snapshot in result.snapshots]),
        cell_center=np.asarray([snapshot.cell_center for snapshot in result.snapshots]),
        cell_angle=np.asarray([snapshot.cell_angle for snapshot in result.snapshots]),
        active_sectors=np.asarray([
            np.isin(np.arange(result.config.n_sectors), snapshot.active_sectors)
            for snapshot in result.snapshots
        ]),
    )
    make_summary_figure(result, output / f"{result.stage}_summary.png")
    if make_gif:
        make_stage_animation(result, output / f"{result.stage}_{result.fixture_name}.gif")
    if make_mp4:
        make_stage_animation(result, output / f"{result.stage}_{result.fixture_name}.mp4")


def main(argv=None):
    parser = argparse.ArgumentParser(description="G3 emergent cell-collagen guidance")
    parser.add_argument("--stage", choices=("g3a", "g3b", "g3c"), required=True)
    parser.add_argument("--fixture", default="single_fibre")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--duration", type=float, default=None, help="override simulated seconds")
    parser.add_argument("--output", required=True)
    parser.add_argument("--feedback-off", action="store_true")
    parser.add_argument("--no-gif", action="store_true")
    parser.add_argument("--mp4", action="store_true")
    parser.add_argument("--load-unload", action="store_true",
                        help="also write the G3A elastic FOI/kappa recovery diagnostic")
    args = parser.parse_args(argv)

    cfg = G3Config.from_yaml(args.config)
    result = run_g3(args.stage, args.fixture, cfg, args.seed, args.duration,
                    feedback_enabled=not args.feedback_off)
    output = Path(args.output)
    save_run(result, output, make_gif=not args.no_gif, make_mp4=args.mp4)
    if args.load_unload:
        diagnostic = run_load_unload(cfg, args.fixture, args.seed)
        with open(output / "load_unload_kappa.json", "w", encoding="utf-8") as handle:
            json.dump(_json_ready(diagnostic), handle, indent=2, ensure_ascii=False)
    print(json.dumps(_json_ready(result.summary), indent=2, ensure_ascii=False))
    return 0 if result.status == "complete" else 2


if __name__ == "__main__":
    raise SystemExit(main())
