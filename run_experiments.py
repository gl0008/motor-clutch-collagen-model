#!/usr/bin/env python3
"""Run verification and fast/slow relaxation ensembles without plotting dependencies."""

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from collagen_model.clutch import MotorClutchParameters
from collagen_model.network import make_diluted_triangular_network
from collagen_model.simulation import SimulationConfig, run_condition
from collagen_model.sls import StandardLinearSolid


def write_csv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def sls_relaxation(dt: float = 0.05, duration: float = 120.0) -> List[Tuple[float, float, float]]:
    k0, kinf, tau, extension = 18.0, 4.5, 12.0, 1.0
    bond = StandardLinearSolid(k0=k0, kinf=kinf, tau=tau)
    bond.set_extension(extension, instantaneous=True)
    rows = []
    for step in range(int(round(duration / dt)) + 1):
        time = step * dt
        analytic = (kinf + (k0 - kinf) * math.exp(-time / tau)) * extension
        rows.append((time, bond.force, analytic))
        if step * dt < duration:
            bond.set_extension(extension, dt=dt)
    return rows


def network_relaxation(
    taus: Sequence[float], dt: float = 0.05, duration: float = 120.0, strain: float = 0.05
) -> List[Tuple[str, float, float, float]]:
    """Virtual step-strain test of the whole bead network.

    The network is deformed affinely and then held fixed. The reported force is
    the mean magnitude of all axial segment forces, normalized by its value at
    the instant of deformation. This verifies that distributing SLS elements
    over a network retains the requested stress-relaxation behavior.
    """

    rows: List[Tuple[str, float, float, float]] = []
    for tau in taus:
        network = make_diluted_triangular_network(seed=77, tau=tau)
        low, _ = network.bounds
        network.positions[:, 0] = low[0] + (1.0 + strain) * (network.positions[:, 0] - low[0])
        network.instantaneous_sync()
        initial_force = float(np.mean(network.axial_force_magnitudes()))
        for step in range(int(round(duration / dt)) + 1):
            time = step * dt
            normalized = float(np.mean(network.axial_force_magnitudes()) / initial_force)
            # Because all axial elements have the same SLS parameters and the
            # affine deformation is held, this is also the analytical curve.
            analytic = network.kinf / network.k0 + (
                1.0 - network.kinf / network.k0
            ) * math.exp(-time / tau)
            rows.append(("tau={:g}".format(tau), time, normalized, analytic))
            if step * dt < duration:
                network.internal_forces(dt, zero_fixed=False)
    return rows


def polyline(points: Sequence[Tuple[float, float]], color: str, width: float = 2.0) -> str:
    coords = " ".join("{:.2f},{:.2f}".format(x, y) for x, y in points)
    return '<polyline fill="none" stroke="{}" stroke-width="{}" points="{}"/>'.format(color, width, coords)


def line_panel(
    x: np.ndarray,
    series: Sequence[Tuple[str, np.ndarray, str]],
    left: float,
    top: float,
    width: float,
    height: float,
    title: str,
    x_label: str,
    y_label: str,
) -> List[str]:
    all_y = np.concatenate([values for _, values, _ in series])
    xmin, xmax = float(np.min(x)), float(np.max(x))
    ymin, ymax = float(np.min(all_y)), float(np.max(all_y))
    if xmax == xmin:
        xmax = xmin + 1.0
    if ymax == ymin:
        ymax = ymin + 1.0
    pad = 36.0
    px0, py0 = left + pad, top + height - pad
    plot_w, plot_h = width - 2 * pad, height - 2 * pad

    def transform(xx: float, yy: float) -> Tuple[float, float]:
        return (
            px0 + plot_w * (xx - xmin) / (xmax - xmin),
            py0 - plot_h * (yy - ymin) / (ymax - ymin),
        )

    items = [
        '<rect x="{}" y="{}" width="{}" height="{}" fill="white" stroke="#c8c8c8"/>'.format(left, top, width, height),
        '<text x="{}" y="{}" font-size="15" font-weight="600">{}</text>'.format(left + 8, top + 20, title),
        '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#555"/>'.format(px0, py0, px0 + plot_w, py0),
        '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#555"/>'.format(px0, py0, px0, py0 - plot_h),
        '<text x="{}" y="{}" font-size="11" text-anchor="middle">{}</text>'.format(px0 + plot_w / 2, top + height - 5, x_label),
        '<text x="{}" y="{}" font-size="11" transform="rotate(-90 {} {})" text-anchor="middle">{}</text>'.format(left + 11, top + height / 2, left + 11, top + height / 2, y_label),
        '<text x="{}" y="{}" font-size="9">{:.3g}</text>'.format(px0 - 4, py0 + 12, xmin),
        '<text x="{}" y="{}" font-size="9" text-anchor="end">{:.3g}</text>'.format(px0 + plot_w, py0 + 12, xmax),
        '<text x="{}" y="{}" font-size="9" text-anchor="end">{:.3g}</text>'.format(px0 - 5, py0, ymin),
        '<text x="{}" y="{}" font-size="9" text-anchor="end">{:.3g}</text>'.format(px0 - 5, py0 - plot_h + 4, ymax),
    ]
    legend_x, legend_y = px0 + 8, top + 38
    for idx, (label, values, color) in enumerate(series):
        points = [transform(float(xx), float(yy)) for xx, yy in zip(x, values)]
        items.append(polyline(points, color))
        items.append('<text x="{}" y="{}" font-size="10" fill="{}">{}</text>'.format(legend_x, legend_y + 13 * idx, color, label))
    return items


def trajectory_panel(
    runs: Dict[str, List[Dict[str, np.ndarray]]],
    left: float,
    top: float,
    width: float,
    height: float,
) -> List[str]:
    colors = {"fast": "#087e8b", "slow": "#d1495b"}
    all_x = np.concatenate([run["x"] - run["x"][0] for condition in runs.values() for run in condition])
    all_y = np.concatenate([run["y"] - run["y"][0] for condition in runs.values() for run in condition])
    limit = max(0.5, float(np.max(np.abs(np.concatenate([all_x, all_y])))))
    pad = 36.0
    px0, py0 = left + pad, top + height - pad
    plot_w, plot_h = width - 2 * pad, height - 2 * pad

    def transform(xx: float, yy: float) -> Tuple[float, float]:
        return (px0 + plot_w * (xx + limit) / (2 * limit), py0 - plot_h * (yy + limit) / (2 * limit))

    items = [
        '<rect x="{}" y="{}" width="{}" height="{}" fill="white" stroke="#c8c8c8"/>'.format(left, top, width, height),
        '<text x="{}" y="{}" font-size="15" font-weight="600">Cell trajectories</text>'.format(left + 8, top + 20),
        '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#ddd"/>'.format(*transform(-limit, 0), *transform(limit, 0)),
        '<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#ddd"/>'.format(*transform(0, -limit), *transform(0, limit)),
        '<text x="{}" y="{}" font-size="11" text-anchor="middle">Δx</text>'.format(px0 + plot_w / 2, top + height - 5),
        '<text x="{}" y="{}" font-size="11" transform="rotate(-90 {} {})" text-anchor="middle">Δy</text>'.format(left + 11, top + height / 2, left + 11, top + height / 2),
    ]
    for condition, condition_runs in runs.items():
        for run in condition_runs:
            x = run["x"] - run["x"][0]
            y = run["y"] - run["y"][0]
            items.append(polyline([transform(float(xx), float(yy)) for xx, yy in zip(x, y)], colors[condition], 1.2))
    items.append('<text x="{}" y="{}" font-size="10" fill="{}">fast</text>'.format(px0 + 8, top + 38, colors["fast"]))
    items.append('<text x="{}" y="{}" font-size="10" fill="{}">slow</text>'.format(px0 + 8, top + 52, colors["slow"]))
    return items


def write_svg(path: Path, relaxation_rows, runs, summaries) -> None:
    width, height = 1080, 720
    items = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" viewBox="0 0 {} {}">'.format(width, height, width, height),
        '<rect width="100%" height="100%" fill="#f5f7f8"/>',
        '<text x="30" y="32" font-size="22" font-weight="700">SLS bead–collagen motor–clutch prototype</text>',
        '<text x="30" y="52" font-size="12" fill="#555">Verification plus exploratory fast/slow ensembles; parameters are not yet fitted.</text>',
    ]
    relaxation = np.asarray(relaxation_rows)
    items += line_panel(
        relaxation[:, 0],
        [("numerical", relaxation[:, 1], "#087e8b"), ("analytic", relaxation[:, 2], "#222222")],
        25, 75, 500, 280, "Single SLS step-strain test", "time", "force"
    )
    items += trajectory_panel(runs, 555, 75, 500, 280)
    time = summaries["fast"]["time"]
    items += line_panel(
        time,
        [("fast τ", summaries["fast"]["msd"], "#087e8b"), ("slow τ", summaries["slow"]["msd"], "#d1495b")],
        25, 385, 500, 280, "Ensemble mean squared displacement", "time", "MSD"
    )
    items += line_panel(
        time,
        [("fast τ", summaries["fast"]["bound"], "#087e8b"), ("slow τ", summaries["slow"]["bound"], "#d1495b")],
        555, 385, 500, 280, "Mean bound clutches", "time", "count"
    )
    items.append("</svg>")
    path.write_text("\n".join(items), encoding="utf-8")


def write_network_relaxation_svg(path: Path, rows) -> None:
    grouped = {}
    for label, time, numerical, analytic in rows:
        grouped.setdefault(label, {"time": [], "numerical": [], "analytic": []})
        grouped[label]["time"].append(time)
        grouped[label]["numerical"].append(numerical)
        grouped[label]["analytic"].append(analytic)
    colors = ["#087e8b", "#d1495b", "#7a5195", "#ef9b20"]
    first = next(iter(grouped.values()))
    series = []
    for color, (label, data) in zip(colors, grouped.items()):
        series.append((label + " numerical", np.asarray(data["numerical"]), color))
    # A single analytic curve per tau would sit directly under its numerical
    # curve; the CSV retains both for exact error checking.
    items = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="720" height="430" viewBox="0 0 720 430">',
        '<rect width="100%" height="100%" fill="#f5f7f8"/>',
        '<text x="25" y="30" font-size="20" font-weight="700">Network-level step-strain relaxation</text>',
        '<text x="25" y="50" font-size="12" fill="#555">All collagen segments are SLS bonds; the deformed network is held fixed.</text>',
    ]
    items += line_panel(
        np.asarray(first["time"]), series, 25, 70, 670, 330,
        "Normalized mean axial force", "time", "F(t) / F(0+)"
    )
    items.append("</svg>")
    path.write_text("\n".join(items), encoding="utf-8")


def write_network_snapshot_svg(path: Path, run: Dict[str, np.ndarray]) -> None:
    reference = run["reference_network_positions"]
    final = run["final_network_positions"]
    edges = run["edges"]
    low = np.minimum(reference.min(axis=0), final.min(axis=0))
    high = np.maximum(reference.max(axis=0), final.max(axis=0))
    span = np.maximum(high - low, 1.0e-9)
    width, height, pad = 900.0, 520.0, 45.0

    def transform(point: np.ndarray) -> Tuple[float, float]:
        return (
            pad + (width - 2 * pad) * float((point[0] - low[0]) / span[0]),
            height - pad - (height - 2 * pad) * float((point[1] - low[1]) / span[1]),
        )

    items = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="520" viewBox="0 0 900 520">',
        '<rect width="100%" height="100%" fill="#f7f8f9"/>',
        '<text x="24" y="28" font-size="20" font-weight="700">Representative coupled network realization</text>',
        '<text x="24" y="47" font-size="12" fill="#555">gray: reference network; blue: final network; red: cell trajectory</text>',
    ]
    for a, b in edges:
        x1, y1 = transform(reference[a])
        x2, y2 = transform(reference[b])
        items.append('<line x1="{:.2f}" y1="{:.2f}" x2="{:.2f}" y2="{:.2f}" stroke="#c8cdd2" stroke-width="0.7"/>'.format(x1, y1, x2, y2))
    for a, b in edges:
        x1, y1 = transform(final[a])
        x2, y2 = transform(final[b])
        items.append('<line x1="{:.2f}" y1="{:.2f}" x2="{:.2f}" y2="{:.2f}" stroke="#087e8b" stroke-width="0.8" opacity="0.8"/>'.format(x1, y1, x2, y2))
    trajectory = np.column_stack((run["x"], run["y"]))
    items.append(polyline([transform(point) for point in trajectory], "#d1495b", 2.5))
    start_x, start_y = transform(trajectory[0])
    end_x, end_y = transform(trajectory[-1])
    items.append('<circle cx="{:.2f}" cy="{:.2f}" r="4" fill="#222"/>'.format(start_x, start_y))
    items.append('<circle cx="{:.2f}" cy="{:.2f}" r="5" fill="#d1495b"/>'.format(end_x, end_y))
    items.append("</svg>")
    path.write_text("\n".join(items), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--duration", type=float, default=90.0)
    parser.add_argument("--dt", type=float, default=0.02)
    parser.add_argument("--replicates", type=int, default=6)
    parser.add_argument("--fast-tau", type=float, default=5.0)
    parser.add_argument("--slow-tau", type=float, default=50.0)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    config = SimulationConfig(duration=args.duration, dt=args.dt)
    clutch_parameters = MotorClutchParameters()
    relaxation_rows = sls_relaxation()
    write_csv(args.output / "sls_relaxation.csv", ["time", "numerical_force", "analytic_force"], relaxation_rows)
    network_relaxation_rows = network_relaxation([args.fast_tau, args.slow_tau])
    write_csv(
        args.output / "network_relaxation.csv",
        ["condition", "time", "normalized_mean_axial_force", "analytic_normalized_force"],
        network_relaxation_rows,
    )
    write_network_relaxation_svg(args.output / "network_relaxation.svg", network_relaxation_rows)

    conditions = {"fast": args.fast_tau, "slow": args.slow_tau}
    runs: Dict[str, List[Dict[str, np.ndarray]]] = {name: [] for name in conditions}
    trajectory_rows = []
    for condition, tau in conditions.items():
        for replicate in range(args.replicates):
            run = run_condition(
                tau=tau,
                seed=1000 + replicate,
                config=config,
                clutch_parameters=clutch_parameters,
            )
            runs[condition].append(run)
            for idx, time in enumerate(run["time"]):
                trajectory_rows.append([
                    condition, tau, replicate, time, run["x"][idx], run["y"][idx],
                    run["squared_displacement"][idx], run["bound_clutches"][idx],
                    run["mean_clutch_force"][idx], run["completed_lifetime_mean"][idx],
                    run["mean_network_axial_force"][idx],
                ])
    write_csv(
        args.output / "trajectories.csv",
        ["condition", "tau", "replicate", "time", "x", "y", "squared_displacement", "bound_clutches", "mean_clutch_force", "completed_lifetime_mean", "mean_network_axial_force"],
        trajectory_rows,
    )

    summaries = {}
    summary_rows = []
    for condition, condition_runs in runs.items():
        summaries[condition] = {
            "time": condition_runs[0]["time"],
            "msd": np.mean([run["squared_displacement"] for run in condition_runs], axis=0),
            "bound": np.mean([run["bound_clutches"] for run in condition_runs], axis=0),
            "force": np.mean([run["mean_clutch_force"] for run in condition_runs], axis=0),
            "lifetime": np.mean([run["completed_lifetime_mean"] for run in condition_runs], axis=0),
        }
        for idx, time in enumerate(summaries[condition]["time"]):
            summary_rows.append([
                condition, conditions[condition], time, summaries[condition]["msd"][idx],
                summaries[condition]["bound"][idx], summaries[condition]["force"][idx],
                summaries[condition]["lifetime"][idx],
            ])
    write_csv(
        args.output / "condition_summary.csv",
        ["condition", "tau", "time", "mean_squared_displacement", "mean_bound_clutches", "mean_clutch_force", "mean_completed_lifetime"],
        summary_rows,
    )
    write_svg(args.output / "summary.svg", relaxation_rows, runs, summaries)

    representative = runs["fast"][0]
    write_csv(
        args.output / "representative_nodes.csv",
        ["bead", "reference_x", "reference_y", "final_x", "final_y"],
        [
            [idx, representative["reference_network_positions"][idx, 0], representative["reference_network_positions"][idx, 1], representative["final_network_positions"][idx, 0], representative["final_network_positions"][idx, 1]]
            for idx in range(len(representative["reference_network_positions"]))
        ],
    )
    write_csv(
        args.output / "representative_edges.csv",
        ["edge", "bead_i", "bead_j"],
        [[idx, edge[0], edge[1]] for idx, edge in enumerate(representative["edges"])],
    )
    write_network_snapshot_svg(args.output / "representative_network.svg", representative)

    metadata = {
        "status": "research prototype; not calibrated or a completed paper reproduction",
        "simulation": config.__dict__,
        "motor_clutch": clutch_parameters.__dict__,
        "conditions": conditions,
        "replicates": args.replicates,
        "notes": [
            "The same random seeds/network realizations are paired across fast and slow conditions.",
            "SLS validation is expected to agree with the analytic solution to numerical precision.",
            "Migration trends are exploratory until bulk network and clutch parameters are fitted.",
        ],
    }
    (args.output / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    final_fast = float(summaries["fast"]["msd"][-1])
    final_slow = float(summaries["slow"]["msd"][-1])
    max_sls_error = max(abs(row[1] - row[2]) for row in relaxation_rows)
    max_network_error = max(abs(row[2] - row[3]) for row in network_relaxation_rows)
    print("Wrote results to {}".format(args.output.resolve()))
    print("SLS max analytic error: {:.3e}".format(max_sls_error))
    print("Network relaxation max analytic error: {:.3e}".format(max_network_error))
    print("Final ensemble MSD: fast={:.5g}, slow={:.5g}".format(final_fast, final_slow))
    print("Exploratory ratio fast/slow: {:.5g}".format(final_fast / final_slow if final_slow else float("inf")))


if __name__ == "__main__":
    main()
