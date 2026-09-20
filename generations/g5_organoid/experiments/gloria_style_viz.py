"""Gloria g4-v4e-style renderer for the G5 organoid runs (visual consistency layer).

Matches the look of https://gl0008.github.io/motor-clutch-collagen-model/g4-v4e.html so our
organoid movies read like Gloria's single-cell viewer instead of the old coolwarm heat-map:

  * neutral GREY collagen fibres (true scale), GOLD crosslink dots that track the beads;
  * cells drawn as open circles (light fill), partial-EMT / leader cells in RED;
  * capped TRACTION ARROWS at each gripping site (length proportional to per-site clutch
    force, clipped to ``force_cap`` exactly like her "capped display scale");
  * RED x at sites that fully failed (all clutches ruptured) since the previous frame;
  * white background, 50 um scale bar, title + wall-clock stamp.

Plus an EVENT MICROSCOPE (``render_event_microscope``): it zooms one gripping site over a
~5-minute window and shows the site-level clutch cycle -- GREEN ring on (re)binding, the
loading arrow growing, RED x at the complete site failure (slip) -- with a force-vs-time
trace, mirroring her "five-minute pull-slip-rebind window".

Consumes the rich per-site snapshots added to ``consistency.run_r0_invasion`` (snapshots=True):
``bead_snapshots, cell_snapshots, edges, crosslink_edge/alpha, site_force_snapshots,
site_point_snapshots, site_normal_snapshots, site_failed_snapshots``.  model.py / visualize.py
are untouched.  Personal-testing visuals (CLAUDE.md 7.5); no swirling claim is implied.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection

# Gloria-style palette
FIBRE = "#7a7a7a"
XLINK = "#e8a800"
CELL_FILL = "#cfe0ee"
CELL_EDGE = "#2c6fa6"
EMT_FILL = "#f4b6b6"
EMT_EDGE = "#c0392b"
ARROW = "#12664f"      # traction arrow (dark teal)
FAIL = "#d62728"       # red x on complete site failure
BIND = "#27ae60"       # green ring on (re)bind


def _xl_points(edges, xl_edge, xl_alpha, beads):
    """Crosslink material points interpolated on the CURRENT bead positions."""
    if xl_edge is None or len(xl_edge) == 0:
        return np.empty((0, 2))
    i = edges[xl_edge, 0]
    j = edges[xl_edge, 1]
    a = xl_alpha[:, None]
    return (1.0 - a) * beads[i] + a * beads[j]


def _scale_bar(ax, span, length=50.0, label="50 um"):
    x0, y0 = -span + 0.10 * span, -span + 0.07 * span
    ax.plot([x0, x0 + length], [y0, y0], color="k", lw=3, solid_capstyle="butt", zorder=6)
    ax.text(x0 + 0.5 * length, y0 + 0.03 * span, label, ha="center", va="bottom", fontsize=8)


def run_out_from_npz(path):
    """Reconstruct the run_out-like dict the renderers need from a saved full npz.

    Lets the sim (expensive, 2 h) and the rendering (cheap, re-runnable) be decoupled so a
    render tweak or bug never forces a re-simulation.  Returns (run_out, meta) where meta
    carries red_ids / red_per_frame / title / red_label / cell_radius / summary scalars.
    """
    z = np.load(path, allow_pickle=True)
    times = z["times"]
    run_out = {
        "bead_snapshots": z["bead_snapshots"],
        "cell_snapshots": z["cell_snapshots"],
        "edges": z["edges"],
        "site_force_snapshots": z["site_force_snapshots"],
        "site_point_snapshots": z["site_point_snapshots"],
        "site_normal_snapshots": z["site_normal_snapshots"],
        "site_failed_snapshots": z["site_failed_snapshots"],
        "crosslink_edge": z["crosslink_edge"],
        "crosslink_alpha": z["crosslink_alpha"],
        "n_contact_sectors": int(z["n_contact_sectors"]),
        "frames": [{"time": float(t)} for t in times],
    }
    meta = {
        "red_ids": [int(x) for x in z["red_ids"]] if "red_ids" in z.files else [],
        "red_per_frame": ([list(r) for r in z["red_per_frame"]] if "red_per_frame" in z.files else None),
        "title": str(z["title"]) if "title" in z.files else "",
        "red_label": str(z["red_label"]) if "red_label" in z.files else "EMT/leader",
        "cell_radius": float(z["cell_radius"]) if "cell_radius" in z.files else 9.0,
        "summary": {k: (float(z[k]) if z[k].dtype.kind == "f" else int(z[k]))
                    for k in ("mean_inv", "max_inv", "detached", "rg", "lcc", "fpair", "grip_fails")
                    if k in z.files},
    }
    return run_out, meta


def _unpack(run_out):
    return dict(
        beads=np.asarray(run_out["bead_snapshots"]),
        cells=np.asarray(run_out["cell_snapshots"]),
        edges=np.asarray(run_out["edges"]),
        sf=np.asarray(run_out["site_force_snapshots"]),
        spt=np.asarray(run_out["site_point_snapshots"]),
        snrm=np.asarray(run_out["site_normal_snapshots"]),
        sfail=np.asarray(run_out["site_failed_snapshots"]),
        xl_edge=run_out.get("crosslink_edge"),
        xl_alpha=run_out.get("crosslink_alpha"),
        n_sec=int(run_out.get("n_contact_sectors", 12)),
        times=[fr["time"] for fr in run_out["frames"]],
    )


def _net_cell_traction(fk, nk, n_sec, n_cells):
    """Per-cell NET traction vector (outward = invasion direction): -sum_s force_s * inward_s.

    One resultant arrow per cell instead of one per grip site -- far less clutter, and it points
    the way the cell is actually pulling / advancing.
    """
    net = np.zeros((n_cells, 2))
    live = np.isfinite(nk[:, 0]) & (fk > 1e-9)
    for s in np.flatnonzero(live):
        net[s // n_sec] -= fk[s] * nk[s]
    return net


def _draw_main(ax, k, d, *, span, force_cap, red, cell_radius, title, red_label,
               arrow_um_per_nN, red_per_frame=None, arrows="net", show_crosslinks=False,
               show_failures=False, show_trails=True):
    """Draw one clean main-panel frame (shared by gif + png).

    ``arrows``: ``"net"`` = one resultant traction arrow per cell (default, readable),
    ``"site"`` = one per grip site (dense), ``"none"``.  ``red`` / ``red_per_frame`` mark the
    partial-EMT / leader cells (per-frame for R3's baton).  Trails show each cell's path from t0.
    """
    if red_per_frame is not None:
        red = set(int(i) for i in red_per_frame[k])
    ax.clear(); ax.set_facecolor("white")
    beads, cells, edges = d["beads"], d["cells"], d["edges"]
    F = len(beads)
    pos = beads[k]
    segs = np.stack([pos[edges[:, 0]], pos[edges[:, 1]]], axis=1)
    ax.add_collection(LineCollection(segs, colors=FIBRE, linewidths=0.45, alpha=0.55, zorder=1))
    if show_crosslinks:
        xlp = _xl_points(edges, d["xl_edge"], d["xl_alpha"], pos)
        if len(xlp):
            ax.scatter(xlp[:, 0], xlp[:, 1], s=2.0, c=XLINK, zorder=2, linewidths=0, alpha=0.5)
    c0 = cells[0]
    # faint displacement trails (t0 -> now): makes "who moved / who escaped" obvious
    if show_trails:
        trails = np.stack([c0, cells[k]], axis=1)
        ax.add_collection(LineCollection(trails, colors="#9a9a9a", linewidths=0.7, alpha=0.5, zorder=2))
    # traction arrows
    fk, pk, nk = d["sf"][k], d["spt"][k], d["snrm"][k]
    if arrows == "site":
        live = np.isfinite(pk[:, 0]) & (fk > 1e-9)
        if np.any(live):
            mag = np.minimum(fk[live], force_cap) * arrow_um_per_nN
            ax.quiver(pk[live, 0], pk[live, 1], nk[live, 0] * mag, nk[live, 1] * mag,
                      angles="xy", scale_units="xy", scale=1.0, color=ARROW, width=0.004,
                      headwidth=4, headlength=5, alpha=0.9, zorder=4)
    elif arrows == "net":
        # ONE fixed-length DIRECTIONAL arrow per cell (which way it is pulling / advancing) --
        # magnitude is deliberately not encoded, so the field stays clean and comparable.
        net = _net_cell_traction(fk, nk, d["n_sec"], len(cells[k]))
        mags = np.linalg.norm(net, axis=1)
        keep = mags > 1e-9
        if np.any(keep):
            L = 1.5 * cell_radius
            u = net[keep, 0] / mags[keep] * L
            v = net[keep, 1] / mags[keep] * L
            ax.quiver(cells[k][keep, 0], cells[k][keep, 1], u, v, angles="xy", scale_units="xy",
                      scale=1.0, color=ARROW, width=0.007, headwidth=4, headlength=5, alpha=0.85, zorder=4)
    if show_failures:
        fail = np.isfinite(pk[:, 0]) & d["sfail"][k]
        if np.any(fail):
            ax.scatter(pk[fail, 0], pk[fail, 1], s=42, c=FAIL, marker="x", linewidths=1.6, zorder=5)
    for i, c in enumerate(cells[k]):
        is_red = i in red
        ax.add_patch(plt.Circle(c, cell_radius, facecolor=(EMT_FILL if is_red else CELL_FILL),
                                edgecolor=(EMT_EDGE if is_red else CELL_EDGE), lw=1.1,
                                alpha=0.95, zorder=3))
    invaded = float(np.mean(np.linalg.norm(cells[k], axis=1) - np.linalg.norm(c0, axis=1)))
    ax.set_xlim(-span, span); ax.set_ylim(-span, span); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    _scale_bar(ax, span)
    arrow_note = {"net": "arrow = each cell's net pull", "site": "arrows = per-site traction",
                  "none": ""}[arrows]
    ax.set_title("%s\nt = %.0f min   red = %s   mean invasion %+.1f um   %s"
                 % (title, d["times"][k] / 60.0, red_label, invaded, arrow_note), fontsize=8.6)


def _auto_cap(sf):
    return float(np.nanpercentile(sf[sf > 0], 90)) if np.any(sf > 0) else 1.0


def render_gloria_style(run_out, out_gif, *, red_ids=(), red_label="EMT/leader",
                        cell_radius=9.0, span=None, fps=8, title="", force_cap=None,
                        arrow_um_per_nN=1.6, red_per_frame=None, arrows="net",
                        show_crosslinks=False, show_failures=False, show_trails=True):
    """Clean main-panel animation of a full run (snapshots=True required).

    Defaults are the readable style: light fibres, one NET traction arrow per cell, cell paths,
    no crosslink dots / no per-site failure marks.  ``red_per_frame`` moves the red circle for
    R3's baton.  Set ``arrows='site'`` / ``show_crosslinks=True`` for the dense diagnostic look.
    """
    d = _unpack(run_out)
    red = set(int(i) for i in red_ids)
    if span is None:
        span = float(np.max(np.abs(d["beads"][0]))) * 1.02
    if force_cap is None:
        force_cap = _auto_cap(d["sf"])
    fig, ax = plt.subplots(figsize=(6.8, 7.0))
    FuncAnimation(fig, lambda k: _draw_main(
        ax, k, d, span=span, force_cap=force_cap, red=red, cell_radius=cell_radius, title=title,
        red_label=red_label, arrow_um_per_nN=arrow_um_per_nN, red_per_frame=red_per_frame,
        arrows=arrows, show_crosslinks=show_crosslinks, show_failures=show_failures,
        show_trails=show_trails),
        frames=len(d["beads"]), interval=1000 / fps).save(out_gif, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return out_gif


def render_frame_png(run_out, out_png, *, red_ids=(), red_label="EMT/leader", cell_radius=9.0,
                     span=None, title="", force_cap=None, frame=0, arrow_um_per_nN=1.6,
                     red_per_frame=None, arrows="net", show_crosslinks=False,
                     show_failures=False, show_trails=True):
    """A single clean PNG (quick look / write-up hero)."""
    d = _unpack(run_out)
    red = set(int(i) for i in red_ids)
    if span is None:
        span = float(np.max(np.abs(d["beads"][0]))) * 1.02
    if force_cap is None:
        force_cap = _auto_cap(d["sf"])
    if frame < 0:
        frame = len(d["beads"]) + frame
    fig, ax = plt.subplots(figsize=(6.8, 7.0))
    _draw_main(ax, frame, d, span=span, force_cap=force_cap, red=red, cell_radius=cell_radius,
               title=title, red_label=red_label, arrow_um_per_nN=arrow_um_per_nN,
               red_per_frame=red_per_frame, arrows=arrows, show_crosslinks=show_crosslinks,
               show_failures=show_failures, show_trails=show_trails)
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_png


# ---------------------------------------------------------------------------
# Event microscope: one site's bind -> load -> slip(fail) -> rebind, ~5-min window
# ---------------------------------------------------------------------------
def pick_event_site(sf, sfail, *, half_window=10):
    """Choose the site with the clearest load -> complete-failure -> REBIND cycle.

    Score rewards a high pre-failure peak AND a visible post-failure recovery, and
    prefers events whose full +/-``half_window`` window fits inside the run (so the
    rebind is not truncated).  Returns (site, fail_frame) or (None, None).
    """
    F, S = sf.shape
    best = (-1.0, None, None)
    for s in range(S):
        for kf in np.flatnonzero(sfail[:, s]):
            pre = sf[max(0, kf - 8):kf, s]
            post = sf[kf + 1:min(F, kf + half_window + 2), s]
            if pre.size == 0 or float(pre.max()) <= 1e-9:
                continue
            if post.size == 0 or float(post.max()) <= 1e-9:
                continue                       # require a rebind after the slip
            peak = float(pre.max())
            recovery = float(post.max()) / peak
            fits = 1.0 if (kf - half_window >= 0 and kf + half_window < F) else 0.35
            score = peak * (1.0 + recovery) * fits
            if score > best[0]:
                best = (score, int(s), int(kf))
    return best[1], best[2]


def render_event_microscope(run_out, out_gif, *, cell_radius=9.0, half_window=10,
                            fps=6, title="", zoom_um=14.0, arrow_um_per_nN=1.6):
    """Zoom one gripping site over a pull-slip-rebind window (site-level clutch cycle)."""
    sf = np.asarray(run_out["site_force_snapshots"])
    spt = np.asarray(run_out["site_point_snapshots"])
    snrm = np.asarray(run_out["site_normal_snapshots"])
    sfail = np.asarray(run_out["site_failed_snapshots"])
    beads = np.asarray(run_out["bead_snapshots"])
    cells = np.asarray(run_out["cell_snapshots"])
    edges = np.asarray(run_out["edges"])
    n_sec = int(run_out.get("n_contact_sectors", 12))
    times = [fr["time"] for fr in run_out["frames"]]

    s, kf = pick_event_site(sf, sfail, half_window=half_window)
    if s is None:
        return None
    F = len(beads)
    k0 = max(0, kf - half_window)
    k1 = min(F, kf + half_window + 1)
    frames = list(range(k0, k1))
    cell = s // n_sec
    force_cap = float(np.nanpercentile(sf[sf > 0], 95)) if np.any(sf > 0) else 1.0
    ftrace = sf[k0:k1, s]

    fig, (axz, axt) = plt.subplots(1, 2, figsize=(11.0, 5.4), gridspec_kw={"width_ratios": [1.15, 1.0]})

    def draw(idx):
        k = frames[idx]
        axz.clear(); axt.clear()
        p = spt[k, s]
        if not np.isfinite(p[0]):
            # site empty this frame; center on the owning cell
            p = cells[k, cell]
        axz.set_facecolor("white")
        pos = beads[k]
        segs = np.stack([pos[edges[:, 0]], pos[edges[:, 1]]], axis=1)
        axz.add_collection(LineCollection(segs, colors=FIBRE, linewidths=1.1, alpha=0.9, zorder=1))
        # owning cell edge
        axz.add_patch(plt.Circle(cells[k, cell], cell_radius, facecolor="none",
                                 edgecolor=CELL_EDGE, lw=1.3, zorder=2))
        f = float(sf[k, s])
        bound = f > 1e-9
        prev_bound = (k > 0) and (sf[k - 1, s] > 1e-9)
        if bound:
            nrm = snrm[k, s]
            mag = min(f, force_cap) * arrow_um_per_nN
            axz.quiver(p[0], p[1], nrm[0] * mag, nrm[1] * mag, angles="xy", scale_units="xy",
                       scale=1.0, color=ARROW, width=0.012, headwidth=4, headlength=5, zorder=5)
            axz.scatter([p[0]], [p[1]], s=60, facecolors="none", edgecolors=BIND,
                        linewidths=2.2 if (bound and not prev_bound) else 1.0, zorder=6)
        if sfail[k, s]:
            axz.scatter([p[0]], [p[1]], s=160, c=FAIL, marker="x", linewidths=2.6, zorder=7)
        axz.set_xlim(p[0] - zoom_um, p[0] + zoom_um); axz.set_ylim(p[1] - zoom_um, p[1] + zoom_um)
        axz.set_aspect("equal"); axz.set_xticks([]); axz.set_yticks([])
        state = "SLIP (grip failed)" if sfail[k, s] else ("bound, loading" if bound else "unbound")
        axz.set_title("Event microscope: one grip site (cell %d)\nt = %.1f min   %s"
                      % (cell, times[k] / 60.0, state), fontsize=9)
        # force trace
        tt = np.asarray(times[k0:k1]) / 60.0
        axt.plot(tt, ftrace, "-", color=ARROW, lw=1.6)
        axt.plot(tt[idx], ftrace[idx], "o", color=ARROW, ms=7)
        for j in range(k0, k1):
            if sfail[j, s]:
                axt.axvline(times[j] / 60.0, color=FAIL, ls="--", lw=1.0, alpha=0.7)
        axt.set_xlabel("time (min)"); axt.set_ylabel("site clutch force (nN)")
        axt.set_title("pull -> slip -> rebind (red dashed = complete failure)", fontsize=9)
        axt.grid(alpha=0.25)
        fig.suptitle(title, fontsize=10)

    FuncAnimation(fig, draw, frames=len(frames), interval=1000 / fps).save(out_gif, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return out_gif
