"""Render a G5 before/after figure: collagen reorganised into a radial pattern.

Loads a saved run (initial/final bead positions + edges + cell centres) and draws
each collagen segment coloured by its radial-alignment order (blue = tangential,
red = radial), with the organoid cell disks overlaid.

Usage:
    python generations/g5_organoid/visualize.py output/g5_fullscale_seed23.npz
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


def _order(pos, edges, center):
    seg = pos[edges[:, 1]] - pos[edges[:, 0]]
    tan = seg / np.maximum(np.linalg.norm(seg, axis=1), 1e-12)[:, None]
    mid = 0.5 * (pos[edges[:, 0]] + pos[edges[:, 1]])
    rad = mid - center
    er = rad / np.maximum(np.linalg.norm(rad, axis=1), 1e-12)[:, None]
    return 2.0 * np.square(np.sum(tan * er, axis=1)) - 1.0  # +1 radial, -1 tangential


def _panel(ax, pos, edges, centers, cell_radius, title):
    center = np.zeros(2)
    order = _order(pos, edges, center)
    segs = np.stack([pos[edges[:, 0]], pos[edges[:, 1]]], axis=1)
    lc = LineCollection(segs, cmap="coolwarm", norm=plt.Normalize(-1, 1), linewidths=0.6)
    lc.set_array(order)
    ax.add_collection(lc)
    for c in centers:
        ax.add_patch(plt.Circle(c, cell_radius, color="0.25", alpha=0.55, lw=0))
    span = np.max(np.abs(pos)) * 1.02
    ax.set_xlim(-span, span); ax.set_ylim(-span, span)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=11)
    return lc


def render(npz_path: str, out_path: str | None = None, cell_radius: float = 9.0) -> str:
    data = np.load(npz_path)
    initial, final = data["initial"], data["final"]
    edges, centers = data["edges"], data["centers"]
    out_path = out_path or str(Path(npz_path).with_suffix(".png"))

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.6))
    _panel(axes[0], initial, edges, centers, cell_radius, "t = 0  (disordered)")
    lc = _panel(axes[1], final, edges, centers, cell_radius, "after pull  (near-field radial)")
    cbar = fig.colorbar(lc, ax=axes, fraction=0.035, pad=0.02)
    cbar.set_label("radial order   (-1 tangential  →  +1 radial)")
    fig.suptitle("G5 organoid contracts → collagen turns radial (personal testing, not validated)",
                 fontsize=12)
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _crosslink_points(pos, edges, links):
    """Material-point pairs (a, b) for each crosslink at the given positions."""
    if not len(links):
        return np.empty((0, 2)), np.empty((0, 2))
    links = np.asarray(links, dtype=float)
    ea = links[:, 0].astype(int); aa = links[:, 1][:, None]
    eb = links[:, 2].astype(int); ab = links[:, 3][:, None]
    pa = (1 - aa) * pos[edges[ea, 0]] + aa * pos[edges[ea, 1]]
    pb = (1 - ab) * pos[edges[eb, 0]] + ab * pos[edges[eb, 1]]
    return pa, pb


def animate(snapshots, edges, centers, out_path, *, crosslinks=None, cell_radius=9.0,
            span=None, fps=6, title="G5 organoid contracts -> collagen turns radial",
            show_links=True):
    """Render a GIF of the collagen reorganising, coloured by radial order.

    Fibres are coloured blue (tangential) -> red (radial); crosslinks are drawn as
    small yellow ties so their role is visible; cell disks overlaid.
    """
    from matplotlib.animation import FuncAnimation, PillowWriter

    snapshots = np.asarray(snapshots)
    center = np.zeros(2)
    if span is None:
        span = float(np.max(np.abs(snapshots[0]))) * 1.02
    initial = snapshots[0]
    ghost = np.stack([initial[edges[:, 0]], initial[edges[:, 1]]], axis=1)  # t=0 outline
    fig, ax = plt.subplots(figsize=(6.6, 6.6))

    def draw(k):
        ax.clear()
        pos = snapshots[k]
        # faint t=0 network so displacement from the initial geometry is visible
        ax.add_collection(LineCollection(ghost, colors="#c9c2b4", linewidths=0.5, alpha=0.5))
        # bead-spring collagen: each segment IS a spring, coloured by radial order
        order = _order(pos, edges, center)
        segs = np.stack([pos[edges[:, 0]], pos[edges[:, 1]]], axis=1)
        lc = LineCollection(segs, cmap="coolwarm", norm=plt.Normalize(-1, 1), linewidths=0.9)
        lc.set_array(order)
        ax.add_collection(lc)
        # beads as dots (so it reads as bead-and-spring)
        ax.scatter(pos[:, 0], pos[:, 1], s=1.2, c="#40484d", alpha=0.55, linewidths=0)
        # crosslinks as gold ties between the two fibres they join
        if show_links and crosslinks is not None and len(crosslinks):
            pa, pb = _crosslink_points(pos, edges, crosslinks)
            ls = np.stack([pa, pb], axis=1)
            ax.add_collection(LineCollection(ls, colors="#d8a329", linewidths=1.3, alpha=0.95))
            mid = 0.5 * (pa + pb)
            ax.scatter(mid[:, 0], mid[:, 1], s=6, marker="D", c="#d8a329", alpha=0.9, linewidths=0)
        for c in centers:
            ax.add_patch(plt.Circle(c, cell_radius, color="0.30", alpha=0.55, lw=0))
        ax.set_xlim(-span, span); ax.set_ylim(-span, span)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("%s\nframe %d/%d   grey=t0 outline · dots=beads · gold=crosslinks · blue→red=tangential→radial"
                     % (title, k + 1, len(snapshots)), fontsize=8)

    anim = FuncAnimation(fig, draw, frames=len(snapshots), interval=1000 / fps)
    anim.save(out_path, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return out_path


def animate_invasion(bead_snaps, cell_snaps, edges, out_path, *, cell_radius=9.0,
                     span=None, fps=6, crosslinks=None,
                     title="G5 Stage D: organoid cells invade collagen"):
    """Clearer Stage-D GIF: moving cells (bright) + trails + collagen by radial order.

    The headline is CELL MOTION: disks move outward (invasion) leaving trails, while
    the collagen they drag is coloured blue (tangential) -> red (radial).  A plain
    readout reports how far the cells have invaded.
    """
    from matplotlib.animation import FuncAnimation, PillowWriter

    bead_snaps = np.asarray(bead_snaps)
    cell_snaps = np.asarray(cell_snaps)
    c0 = cell_snaps[0]
    center = np.zeros(2)
    if span is None:
        span = float(np.max(np.abs(bead_snaps[0]))) * 1.02
    fig, ax = plt.subplots(figsize=(6.6, 6.6))

    def draw(k):
        ax.clear()
        pos = bead_snaps[k]
        order = _order(pos, edges, center)
        segs = np.stack([pos[edges[:, 0]], pos[edges[:, 1]]], axis=1)
        lc = LineCollection(segs, cmap="coolwarm", norm=plt.Normalize(-1, 1),
                            linewidths=0.6, alpha=0.85)
        lc.set_array(order)
        ax.add_collection(lc)
        # cell trails (start -> current)
        cells = cell_snaps[k]
        trails = np.stack([c0, cells], axis=1)
        ax.add_collection(LineCollection(trails, colors="#111", linewidths=0.6, alpha=0.35))
        for c in cells:
            ax.add_patch(plt.Circle(c, cell_radius, color="#2c7fb8", alpha=0.85, lw=0))
        invaded = float(np.mean(np.linalg.norm(cells, axis=1) - np.linalg.norm(c0, axis=1)))
        ax.set_xlim(-span, span); ax.set_ylim(-span, span)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("%s\nframe %d/%d   cells invaded outward: %+.1f um   (blue tangential -> red radial)"
                     % (title, k + 1, len(bead_snaps), invaded), fontsize=9)

    anim = FuncAnimation(fig, draw, frames=len(bead_snaps), interval=1000 / fps)
    anim.save(out_path, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return out_path


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "output/g5_fullscale_seed23.npz"
    print("wrote", render(src))
