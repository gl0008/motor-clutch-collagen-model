"""Reproducible G3 mechanism GIFs and compact summary figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle

from .simulation import G3RunResult


DISCLAIMER = "Personal simulation / mechanism validation; not yet a realistic 3D tumor-migration prediction."


def _bounds(result: G3RunResult):
    arrays = [snapshot.positions for snapshot in result.snapshots if snapshot.positions.size]
    centers = np.asarray([snapshot.cell_center for snapshot in result.snapshots])
    if arrays:
        points = np.vstack(arrays + [centers]) * 1.0e6
    else:
        points = centers * 1.0e6
    radius = result.config.cell_radius * 1.0e6
    if points.size == 0:
        return (-2 * radius, 2 * radius, -2 * radius, 2 * radius)
    low = points.min(axis=0)
    high = points.max(axis=0)
    span = max(float(np.max(high - low)), 4.0 * radius)
    center = 0.5 * (low + high)
    half = 0.56 * span
    return center[0] - half, center[0] + half, center[1] - half, center[1] + half


def _draw_frame(ax, result: G3RunResult, snapshot_index: int, xlim, ylim,
                show_disclaimer: bool = True):
    snapshot = result.snapshots[snapshot_index]
    ax.clear()
    positions_um = snapshot.positions * 1.0e6
    network = result.fixture.network
    for fibre in range(network.n_fibers):
        points = positions_um[network.fiber_id == fibre]
        ax.plot(points[:, 0], points[:, 1], color="#607d8b", lw=1.4, alpha=0.9)
        ax.scatter(points[:, 0], points[:, 1], s=4, color="#607d8b", alpha=0.65)

    center_um = snapshot.cell_center * 1.0e6
    radius_um = result.config.cell_radius * 1.0e6
    ax.add_patch(Circle(center_um, radius_um, facecolor="#f2b84b", edgecolor="#7a4b00",
                        lw=1.5, alpha=0.35))
    body = center_um + 0.75 * radius_um * np.array([
        np.cos(snapshot.cell_angle), np.sin(snapshot.cell_angle)])
    ax.plot([center_um[0], body[0]], [center_um[1], body[1]], color="#7a4b00", lw=2.0)

    if snapshot.bound_points.size:
        bound = snapshot.bound_points * 1.0e6
        motor = snapshot.motor_points * 1.0e6
        segments = np.stack((bound, motor), axis=1)
        ax.add_collection(LineCollection(segments, colors="#d1495b", linewidths=0.7, alpha=0.5))
        ax.scatter(bound[:, 0], bound[:, 1], s=12, color="#d1495b", marker="o",
                   label="bound material points")

    sector_angles = snapshot.cell_angle + result.protrusions.sector_angles
    sector_radius = 1.25 * radius_um
    sector_xy = center_um + sector_radius * np.column_stack((np.cos(sector_angles),
                                                              np.sin(sector_angles)))
    score = snapshot.geometry_scores
    sizes = 8.0 + 35.0 * score
    ax.scatter(sector_xy[:, 0], sector_xy[:, 1], c=score, s=sizes, cmap="viridis",
               vmin=0.0, vmax=max(1.0, float(np.max(score, initial=0.0))), alpha=0.55,
               edgecolors="none")
    active = snapshot.active_sectors
    if active.size:
        active_xy = sector_xy[active]
        ax.scatter(active_xy[:, 0], active_xy[:, 1], s=70, facecolors="none",
                   edgecolors="#c62828", linewidths=1.8, label="active protrusions")

    trajectory = np.asarray([item.cell_center for item in result.snapshots[:snapshot_index + 1]]) * 1.0e6
    if trajectory.shape[0] > 1:
        ax.plot(trajectory[:, 0], trajectory[:, 1], color="#1565c0", lw=1.5,
                label="cell trajectory")

    force = snapshot.cell_force
    magnitude = np.linalg.norm(force)
    if magnitude > 0.0:
        direction = force / magnitude
        arrow = center_um + 0.9 * radius_um * direction
        ax.annotate("", xy=arrow, xytext=center_um,
                    arrowprops={"arrowstyle": "->", "color": "#1565c0", "lw": 1.8})

    ax.set(xlim=xlim, ylim=ylim, aspect="equal", xlabel="x (µm)", ylabel="y (µm)")
    ax.set_title(
        f"{result.stage.upper()} · {result.fixture_name} · t={snapshot.time:.1f} s\n"
        f"bound={snapshot.bound_points.shape[0]}  FOI={snapshot.foi:.3f}  "
        f"torque={snapshot.cell_torque:.2e} N·m",
        fontsize=10,
    )
    ax.grid(alpha=0.15)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        unique = dict(zip(labels, handles))
        ax.legend(unique.values(), unique.keys(), loc="upper right", fontsize=7)
    if show_disclaimer:
        ax.text(0.5, -0.14, DISCLAIMER, transform=ax.transAxes, ha="center", va="top",
                fontsize=7, color="#555555")


def make_stage_animation(result: G3RunResult, path, fps: int = 10, max_frames: int = 180):
    """Write a GIF or MP4 selected by the file extension."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    count = len(result.snapshots)
    indices = np.unique(np.linspace(0, count - 1, min(count, max_frames)).round().astype(int))
    x0, x1, y0, y1 = _bounds(result)
    fig, ax = plt.subplots(figsize=(7.0, 6.3), constrained_layout=True)

    def update(frame):
        _draw_frame(ax, result, int(indices[frame]), (x0, x1), (y0, y1))

    animation = FuncAnimation(fig, update, frames=indices.size, interval=1000 / fps, repeat=True)
    if path.suffix.lower() == ".mp4":
        animation.save(path, writer=FFMpegWriter(fps=fps, bitrate=1800), dpi=120)
    else:
        animation.save(path, writer=PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    return path


def make_comparison_animation(results, path, fps: int = 10, max_frames: int = 180):
    """Write a synchronized side-by-side control animation (e.g. isotropic/aligned/rotated)."""
    results = list(results)
    if not results:
        raise ValueError("comparison requires at least one result")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n_frames = min(max(len(result.snapshots) for result in results), max_frames)
    frame_fraction = np.linspace(0.0, 1.0, n_frames)
    bounds = [_bounds(result) for result in results]
    fig, axes = plt.subplots(1, len(results), figsize=(5.4 * len(results), 5.3))
    axes = np.atleast_1d(axes)

    def update(frame):
        for ax, result, (x0, x1, y0, y1) in zip(axes, results, bounds):
            index = int(round(frame_fraction[frame] * (len(result.snapshots) - 1)))
            _draw_frame(ax, result, index, (x0, x1), (y0, y1), show_disclaimer=False)
        fig.suptitle("G3C direction controls", fontsize=11, y=0.98)
        fig.text(0.5, 0.015, DISCLAIMER, ha="center", va="bottom", fontsize=7,
                 color="#555555")
        fig.tight_layout(rect=(0.0, 0.05, 1.0, 0.94))

    animation = FuncAnimation(fig, update, frames=n_frames, interval=1000 / fps, repeat=True)
    if path.suffix.lower() == ".mp4":
        animation.save(path, writer=FFMpegWriter(fps=fps, bitrate=2400), dpi=110)
    else:
        animation.save(path, writer=PillowWriter(fps=fps), dpi=90)
    plt.close(fig)
    return path


def make_summary_figure(result: G3RunResult, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    traces = result.traces
    time = traces["time"]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8))

    axes[0].plot(time, traces["bound_count"], color="#d1495b")
    axes[0].set(xlabel="time (s)", ylabel="bound clutches", title="Attachment kinetics")

    axes[1].plot(time, traces["foi"], color="#2e7d32")
    axes[1].axhline(2.0 / np.pi, color="#777777", lw=1.0, ls="--", label="random 2/π")
    axes[1].set(xlabel="time (s)", ylabel="FOI", title="Collagen orientation")
    axes[1].legend(fontsize=7)

    x = traces["cell_x"] * 1.0e6
    y = traces["cell_y"] * 1.0e6
    axes[2].plot(x, y, color="#1565c0")
    axes[2].scatter(x[:1], y[:1], color="#2e7d32", s=25, label="start")
    axes[2].scatter(x[-1:], y[-1:], color="#c62828", s=25, label="end")
    axes[2].set(xlabel="x (µm)", ylabel="y (µm)", title="Rigid-cell trajectory", aspect="equal")
    if max(float(np.ptp(x)), float(np.ptp(y)), 0.0) < 0.01:
        axes[2].ticklabel_format(axis="both", style="sci", scilimits=(-2, 2))
    axes[2].legend(fontsize=7)

    for ax in axes:
        ax.grid(alpha=0.18)
    fig.suptitle(f"{result.stage.upper()} · {result.fixture_name}", fontsize=11, y=0.98)
    fig.text(0.5, 0.015, DISCLAIMER, ha="center", va="bottom", fontsize=7, color="#555555")
    fig.tight_layout(rect=(0.0, 0.07, 1.0, 0.90))
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path
