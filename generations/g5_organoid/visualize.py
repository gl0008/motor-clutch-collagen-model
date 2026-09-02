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


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "output/g5_fullscale_seed23.npz"
    print("wrote", render(src))
