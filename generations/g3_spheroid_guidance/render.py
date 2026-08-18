"""Render a saved G3 spheroid run as GIFs in Gloria's G2 visual grammar.

Two view modes (like the G2-v3 notebook toggle):
  * ``full``   -- the whole collagen field at true 1x scale;
  * ``follow`` -- a near-cell zoom that tracks the migrating spheroid.

Visual vocabulary matches the corrected G2 animations:
  collagen springs coloured by strain (red tension / blue compression / grey slack),
  gold diamonds for permanent crosslinks, black squares for fixed boundary beads,
  a peach spheroid body, green protrusion shafts with tip markers, magenta captured
  material points joined to the tips by red clutch spokes, a membrane heat-map for the
  emergent polarity, a faint ghost of the initial network, and a scale bar.  Geometry is
  always true 1x; nothing is display-magnified.
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle

# --- G2 palette ------------------------------------------------------------------
PAPER = "#f4f1ea"
GHOST = "#c7cdc9"
FIBER = "#5b6f79"
TENSION = "#b93a35"
COMPRESS = "#347bb0"
GOLD = "#b17d28"
INK = "#17232b"
PEACH = "#f7e3d4"
ORANGE = "#d45f32"
GREEN = "#2e7c61"
MAGENTA = "#b5379b"


def _material_points(pos, edges, links):
    a = pos[edges[links[:, 0].astype(int), 0]]
    b = pos[edges[links[:, 0].astype(int), 1]]
    pa = (1 - links[:, 1])[:, None] * a + links[:, 1][:, None] * b
    return pa


def _strain_colors(strain, scale=0.006):
    """Diverging blue-grey-red mapping, matching the G2 tension/compression cue."""
    t = np.clip(strain / scale, -1.0, 1.0)
    c = np.empty((len(strain), 4))
    base = np.array(matplotlib.colors.to_rgb(FIBER))
    red = np.array(matplotlib.colors.to_rgb(TENSION))
    blue = np.array(matplotlib.colors.to_rgb(COMPRESS))
    for k in range(3):
        c[:, k] = np.where(t >= 0, base[k] + t * (red[k] - base[k]),
                           base[k] + (-t) * (blue[k] - base[k]))
    c[:, 3] = 0.35 + 0.55 * np.abs(t)
    return c


def render(frames, meta, out_path, mode="full", fps=12, max_frames=90,
           title="G3 - spheroid collagen remodelling"):
    if len(frames) > max_frames:
        keep = np.unique(np.linspace(0, len(frames) - 1, max_frames).round().astype(int))
        frames = [frames[i] for i in keep]
    edges = np.asarray(meta["edges"])
    fixed = np.asarray(meta["fixed"], dtype=bool)
    links = np.asarray(meta["crosslinks"], dtype=float) if len(meta["crosslinks"]) else np.zeros((0, 4))
    r0 = np.asarray(meta["initial_positions"])
    R = meta["cell_radius"]
    half = meta["domain_size"] / 2.0

    seg0 = np.stack([r0[edges[:, 0]], r0[edges[:, 1]]], axis=1)
    centers = np.array([f["center"] for f in frames])
    act_vmax = max(float(np.max([f["activity"].max() for f in frames])), 1e-6)

    fig, ax = plt.subplots(figsize=(7.6, 7.6))
    fig.patch.set_facecolor(PAPER)
    ax.set_facecolor(PAPER)
    ax.set_aspect("equal")
    ax.axis("off")

    def draw(idx):
        ax.clear()
        ax.set_facecolor(PAPER)
        ax.set_aspect("equal")
        ax.axis("off")
        fr = frames[idx]
        pos = fr["positions"]
        center = fr["center"]

        # view window
        if mode == "full":
            ax.set_xlim(-half, half); ax.set_ylim(-half, half)
            win = half
        else:  # follow-cell zoom
            win = 2.6 * R
            ax.set_xlim(center[0] - win, center[0] + win)
            ax.set_ylim(center[1] - win, center[1] + win)

        # ghost of the initial network
        ax.add_collection(LineCollection(seg0, colors=GHOST, linewidths=0.5, alpha=0.45, zorder=1))

        # current springs coloured by strain
        seg = np.stack([pos[edges[:, 0]], pos[edges[:, 1]]], axis=1)
        ax.add_collection(LineCollection(seg, colors=_strain_colors(fr["strain"]),
                                         linewidths=1.0, zorder=2))

        # collagen beads drawn as circles (G2 signature: every circle is a bead)
        bm = (np.abs(pos[:, 0] - center[0]) < win) & (np.abs(pos[:, 1] - center[1]) < win)
        bead_s = 4.5 if mode != "full" else 1.8
        ax.scatter(pos[bm, 0], pos[bm, 1], s=bead_s, facecolors="#fbfdfb",
                   edgecolors=FIBER, linewidths=0.35, zorder=3)

        # fixed boundary anchors (only when they are in view)
        fb = r0[fixed]
        m = (np.abs(fb[:, 0] - center[0]) < win) & (np.abs(fb[:, 1] - center[1]) < win)
        ax.scatter(fb[m, 0], fb[m, 1], marker="s", s=9, c=INK, zorder=3)

        # permanent crosslinks as gold diamonds
        if len(links):
            mp = _material_points(pos, edges, links)
            m = (np.abs(mp[:, 0] - center[0]) < win) & (np.abs(mp[:, 1] - center[1]) < win)
            ax.scatter(mp[m, 0], mp[m, 1], marker="D", s=8, c=GOLD, alpha=0.85, zorder=4)

        # spheroid body
        ax.add_patch(Circle(center, R, facecolor=PEACH, edgecolor=ORANGE, lw=2.0, zorder=5))

        # trajectory so far
        if idx > 0:
            ax.plot(centers[:idx + 1, 0], centers[:idx + 1, 1], color=COMPRESS, lw=1.4, alpha=0.8, zorder=5)

        # protrusions, captured points, clutch spokes
        act = fr["activity"]; length = fr["length"]; engaged = fr["engaged"]
        tips = fr["tips"]; epts = fr["engaged_points"]
        angles = np.linspace(0.0, 2 * np.pi, len(act), endpoint=False)
        normals = np.column_stack([np.cos(angles), np.sin(angles)])
        bases = center[None, :] + R * normals
        for k in range(len(act)):
            # protrusion shaft: green while probing, orange once it grips a fibre
            col = ORANGE if engaged[k] else GREEN
            lw = 2.4 if engaged[k] else 1.4
            ax.plot([bases[k, 0], tips[k, 0]], [bases[k, 1], tips[k, 1]],
                    color=col, lw=lw, alpha=(1.0 if engaged[k] else 0.75),
                    solid_capstyle="round", zorder=6)
            ax.scatter(tips[k, 0], tips[k, 1], s=(20 if engaged[k] else 9),
                       c=col, alpha=(1.0 if engaged[k] else 0.75), zorder=7)
            # clutch bundle: the red spoke from the tip to the gripped collagen point
            if engaged[k] and np.isfinite(epts[k, 0]):
                ax.plot([tips[k, 0], epts[k, 0]], [tips[k, 1], epts[k, 1]],
                        color=TENSION, lw=2.0, zorder=8)
                ax.scatter(epts[k, 0], epts[k, 1], s=22, c=MAGENTA,
                           edgecolors=TENSION, linewidths=0.5, zorder=9)

        # scale bar
        bar = 20.0
        x0 = ax.get_xlim()[0] + 0.08 * (2 * win)
        y0 = ax.get_ylim()[0] + 0.09 * (2 * win)
        ax.plot([x0, x0 + bar], [y0, y0], color=INK, lw=3)
        ax.text(x0 + bar / 2, y0 + 0.015 * (2 * win), "20 um", ha="center", va="bottom",
                fontsize=8, color=INK)

        # near-shell radial order for this frame (+1 radial / 0 random / -1 tangential)
        s2 = pos[edges[:, 1]] - pos[edges[:, 0]]
        tg = s2 / np.maximum(np.linalg.norm(s2, axis=1), 1e-12)[:, None]
        midp = 0.5 * (pos[edges[:, 0]] + pos[edges[:, 1]]) - center
        rr = np.linalg.norm(midp, axis=1)
        er = midp / np.maximum(rr, 1e-12)[:, None]
        ordr = 2.0 * np.square(np.sum(tg * er, axis=1)) - 1.0
        nsh = (rr - R >= 0.0) & (rr - R < 20.0)
        sr = float(np.mean(ordr[nsh])) if np.any(nsh) else 0.0

        # header + metric strip
        t_min = fr["time"] / 60.0
        ax.set_title(f"{title}\n({'full 180 um field' if mode=='full' else 'follow-cell zoom'})",
                     fontsize=12, color=INK, loc="left")
        info = (f"t = {t_min:4.1f} min    grip = {int(engaged.sum())}    "
                f"traction = {fr['traction'].sum():4.1f} nN    radial order = {sr:+.2f}")
        ax.text(0.02, 0.025, info, transform=ax.transAxes, fontsize=9, color=INK, va="bottom",
                bbox=dict(boxstyle="round", fc="white", ec=GHOST, alpha=0.9))
        ax.text(0.985, 0.985, "Mechanism demo - not a 3D tumour-migration prediction",
                transform=ax.transAxes, fontsize=6.5, color="#6b7378", ha="right", va="top")

    anim = FuncAnimation(fig, draw, frames=len(frames), interval=1000 / fps)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(out_path), writer=PillowWriter(fps=fps))
    plt.close(fig)
    print("wrote", out_path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="directory with frames.pkl + meta.pkl")
    ap.add_argument("--out", required=True, help="output directory for the GIFs")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--stem", default="g3_spheroid")
    args = ap.parse_args()

    run = Path(args.run)
    with open(run / "frames.pkl", "rb") as fh:
        frames = pickle.load(fh)
    with open(run / "meta.pkl", "rb") as fh:
        meta = pickle.load(fh)

    out = Path(args.out)
    render(frames, meta, out / f"{args.stem}_full.gif", mode="full", fps=args.fps)
    render(frames, meta, out / f"{args.stem}_follow.gif", mode="follow", fps=args.fps)


if __name__ == "__main__":
    main()
