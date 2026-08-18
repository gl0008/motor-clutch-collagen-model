"""Reproducible G3 mechanism GIFs and compact summary figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle, FancyBboxPatch

from .simulation import G3RunResult


DISCLAIMER = "Personal simulation / mechanism validation; not yet a realistic 3D tumor-migration prediction."
PAPER = "#f4f0e8"
PANEL = "#fffdf8"
INK = "#172a35"
FIBRE = "#48636f"
GHOST = "#a9b4b8"
GOLD = "#b17d28"
ORANGE = "#d45f32"
GREEN = "#2e7c61"
RED = "#b93a35"


def _add_scale_bar(ax, xlim, ylim, length_um: float = 10.0) -> None:
    """G1/G2-style physical scale bar; replaces numerical axes and grid lines."""
    x_span = xlim[1] - xlim[0]
    y_span = ylim[1] - ylim[0]
    length = min(length_um, 0.28 * x_span)
    x0 = xlim[0] + 0.045 * x_span
    y0 = ylim[0] + 0.055 * y_span
    ax.plot([x0, x0 + length], [y0, y0], color=INK, lw=1.25,
            solid_capstyle="butt", zorder=30)
    ax.text(x0, y0 + 0.018 * y_span, f"{length:g} µm", ha="left", va="bottom",
            fontsize=6.5, color=INK, zorder=30)


def _add_metric_cards(fig, panel_left: float, count: int = 4):
    """Create fixed-width metric cards matching the G1/G2 browser panels."""
    panel_width = 0.477
    inner = 0.012
    gap = 0.006
    width = (panel_width - 2.0 * inner - gap * (count - 1)) / count
    artists = []
    for slot in range(count):
        x = panel_left + inner + slot * (width + gap)
        fig.add_artist(FancyBboxPatch(
            (x, 0.035), width, 0.105, transform=fig.transFigure,
            boxstyle="round,pad=0.002,rounding_size=0.010",
            facecolor="#f0eee7", edgecolor="none", zorder=1,
        ))
        artists.append(fig.text(
            x + 0.007, 0.087, "", ha="left", va="center", fontsize=7.5,
            color=INK, fontweight="bold", family="sans-serif", zorder=2,
        ))
    return artists


def _contact_fibre_ids(result: G3RunResult) -> set[int]:
    external = result.fixture.network.external_network
    if external is not None:
        return set(map(int, external.contact_fibers))
    return set(range(min(result.config.fibre_count, result.fixture.network.n_fibers)))


def _fibre_reorientation_deg(result: G3RunResult, snapshot, fibre_ids=None) -> float:
    """Maximum true geometric segment rotation relative to frame zero."""
    bonds = result.fixture.network.fiber_bonds
    if bonds.size == 0:
        return 0.0
    edge_fibre = result.fixture.network.fiber_id[bonds[:, 0]]
    chosen = set(map(int, fibre_ids or _contact_fibre_ids(result)))
    mask = np.isin(edge_fibre, list(chosen))
    if not np.any(mask):
        return 0.0
    initial = result.initial_positions
    before = initial[bonds[mask, 1]] - initial[bonds[mask, 0]]
    after = snapshot.positions[bonds[mask, 1]] - snapshot.positions[bonds[mask, 0]]
    before /= np.maximum(np.linalg.norm(before, axis=1)[:, None], 1.0e-30)
    after /= np.maximum(np.linalg.norm(after, axis=1)[:, None], 1.0e-30)
    # Collagen is nematic: reversing a tangent is the same physical direction.
    cosine = np.clip(np.abs(np.sum(before * after, axis=1)), 0.0, 1.0)
    return float(np.degrees(np.arccos(cosine)).max(initial=0.0))


def _protrusion_switches(result: G3RunResult, snapshot_index: int) -> int:
    changes = 0
    previous = None
    for snapshot in result.snapshots[:snapshot_index + 1]:
        current = tuple(sorted(map(int, snapshot.active_sectors)))
        if previous is not None and current != previous:
            changes += 1
        previous = current
    return changes


def _mechanism_phase(snapshot) -> str:
    if snapshot.bound_points.size == 0:
        if snapshot.protrusion_lengths.max(initial=0.0) < 0.05e-6:
            return "POLARITY FORMS"
        return "PROTRUSION GROWS / SEARCHES"
    force = float(np.linalg.norm(snapshot.clutch_forces, axis=1).sum())
    return "ATTACHED" if force < 0.005e-9 else "CLUTCH LOADS COLLAGEN"


def _draw_gloria_frame(
    ax,
    result: G3RunResult,
    snapshot_index: int,
    xlim,
    ylim,
    *,
    show_displacement_vectors: bool = False,
    displacement_scale: float = 1000.0,
    full_field: bool = False,
):
    """Draw one readable G3 frame using the established G2 visual grammar."""
    snapshot = result.snapshots[snapshot_index]
    ax.clear()
    ax.set_facecolor("#fbfaf6")
    positions_um = snapshot.positions * 1.0e6
    initial_um = result.initial_positions * 1.0e6
    network = result.fixture.network
    contact_ids = _contact_fibre_ids(result)
    attached_ids = set(map(int, snapshot.bound_fibre_ids))

    for fibre in range(network.n_fibers):
        ids = network.fiber_id == fibre
        initial = initial_um[ids]
        current = positions_um[ids]
        is_contact = fibre in contact_ids
        is_attached = fibre in attached_ids
        if is_attached and result.stage != "g3s":
            colour, width, alpha, bead_size = ORANGE, 1.55, 1.0, 5.0
        elif is_contact and result.stage != "g3s":
            colour, width, alpha, bead_size = "#287f7b", 1.05, 0.92, 3.0
        else:
            colour = FIBRE
            width = 0.38 if full_field else 0.45
            alpha = 0.42 if full_field else 0.24
            bead_size = 0.65 if full_field else 0.8
        if show_displacement_vectors or is_attached:
            ax.plot(initial[:, 0], initial[:, 1], color=GHOST, lw=0.55,
                    alpha=0.52, ls="--", zorder=1)
        ax.plot(current[:, 0], current[:, 1], color=colour, lw=width,
                alpha=alpha, solid_capstyle="round", zorder=2 if not is_attached else 5)
        ax.scatter(current[:, 0], current[:, 1], s=bead_size, color=colour,
                   alpha=alpha, linewidths=0, zorder=3 if is_contact else 2)

    # G2-style material-point crosslink connectors and diamonds. Links touching
    # an attached fibre are emphasized to expose the force-transmission path.
    if network.xl_edge_a is not None and network.xl_edge_a.size:
        ea = network.fiber_bonds[network.xl_edge_a]
        eb = network.fiber_bonds[network.xl_edge_b]
        aa = network.xl_alpha_a[:, None]
        ab = network.xl_alpha_b[:, None]
        pa = (1.0 - aa) * positions_um[ea[:, 0]] + aa * positions_um[ea[:, 1]]
        pb = (1.0 - ab) * positions_um[eb[:, 0]] + ab * positions_um[eb[:, 1]]
        edge_fibre_a = network.fiber_id[ea[:, 0]]
        edge_fibre_b = network.fiber_id[eb[:, 0]]
        transmitting = np.isin(edge_fibre_a, list(attached_ids)) | np.isin(
            edge_fibre_b, list(attached_ids)
        )
        # For a many-patch spheroid, highlighting every link which touches a
        # loaded fibre paints an artificial orange/gold halo.  The surface
        # patch markers and numeric engagement count already convey adhesion;
        # keep the matrix itself visually readable.
        if result.stage == "g3s":
            transmitting[:] = False
        segments = np.stack((pa, pb), axis=1)
        ax.add_collection(LineCollection(
            segments[~transmitting], colors=GOLD, linewidths=0.35,
            alpha=0.26 if full_field else 0.18, zorder=3,
        ))
        if np.any(transmitting):
            ax.add_collection(LineCollection(
                segments[transmitting], colors=GOLD, linewidths=1.0,
                alpha=0.95, zorder=6,
            ))
        midpoint = 0.5 * (pa + pb)
        ax.scatter(midpoint[:, 0], midpoint[:, 1], marker="D",
                   s=np.where(transmitting, 13.0, 3.2), color=GOLD,
                   alpha=np.where(transmitting, 0.95, 0.34), linewidths=0, zorder=6)

    fixed = result.fixture.fixed_mask
    if np.any(fixed) and full_field:
        ax.scatter(positions_um[fixed, 0], positions_um[fixed, 1], marker="s", s=6,
                   color=INK, alpha=0.9, linewidths=0, zorder=7)

    center_um = snapshot.cell_center * 1.0e6
    radius_um = result.config.cell_radius * 1.0e6
    ax.add_patch(Circle(center_um, radius_um, facecolor="#f8eee8", edgecolor=ORANGE,
                        lw=1.8, zorder=10))
    ax.scatter(center_um[0], center_um[1], s=5, color="#287f7b", zorder=11)

    # Intracellular activity stays visible but subordinate to physical
    # protrusions; it is no longer a ring of large competing symbols.
    sector_angles = snapshot.cell_angle + result.protrusions.sector_angles
    if result.stage != "g3s":
        activity_xy = center_um + 1.055 * radius_um * np.column_stack(
            (np.cos(sector_angles), np.sin(sector_angles))
        )
        activity = snapshot.polarity_activity
        activity_size = 2.0 + 18.0 * activity / max(float(activity.max(initial=0.0)), 1.0e-12)
        ax.scatter(activity_xy[:, 0], activity_xy[:, 1], s=activity_size,
                   color="#d78773", alpha=0.32, linewidths=0, zorder=11)

    lengths = snapshot.protrusion_lengths * 1.0e6
    lab_angles = snapshot.cell_angle + result.protrusions.sector_angles
    bases = center_um + radius_um * np.column_stack((np.cos(lab_angles), np.sin(lab_angles)))
    tips = snapshot.protrusion_tips * 1.0e6
    attached_sectors = set(map(int, snapshot.bound_sector_ids))
    active_set = set(map(int, snapshot.active_sectors))
    if result.stage == "g3s":
        # A spheroid is represented by surface adhesion patches, not a ring of
        # cartoon protrusions. A patch turns orange only after physical capture.
        for sector in snapshot.active_sectors:
            sector = int(sector)
            attached = sector in attached_sectors
            ax.scatter(bases[sector, 0], bases[sector, 1], s=28 if attached else 13,
                       facecolor=ORANGE if attached else "#7aa99a", edgecolor="#fffefa",
                       linewidth=0.7, zorder=14)
    # A sector turnover is a physical transition, not teleportation: keep the
    # old shaft visible while it retracts and let the new one grow from zero.
    for sector in ([] if result.stage == "g3s" else np.flatnonzero(lengths > 0.02)):
        sector = int(sector)
        if sector in active_set:
            continue
        ax.plot([bases[sector, 0], tips[sector, 0]],
                [bases[sector, 1], tips[sector, 1]], color="#8fc4b0",
                lw=1.2, alpha=0.38, solid_capstyle="round", zorder=12)
        ax.scatter(tips[sector, 0], tips[sector, 1], s=9,
                   color="#8fc4b0", alpha=0.38, linewidths=0, zorder=12)
    for sector in ([] if result.stage == "g3s" else snapshot.active_sectors):
        sector = int(sector)
        if lengths[sector] <= 0.01:
            continue
        ax.plot([bases[sector, 0], tips[sector, 0]],
                [bases[sector, 1], tips[sector, 1]], color=GREEN,
                lw=3.0, solid_capstyle="round", zorder=13)
        if sector in attached_sectors:
            ax.scatter(tips[sector, 0], tips[sector, 1], s=45,
                       facecolor="#fffefa", edgecolor=RED, linewidth=1.7, zorder=15)
        else:
            ax.scatter(tips[sector, 0], tips[sector, 1], s=25,
                       color=GREEN, edgecolor="#fffefa", linewidth=0.6, zorder=14)

    if snapshot.bound_points.size and result.stage != "g3s":
        bound = snapshot.bound_points * 1.0e6
        motor = snapshot.motor_points * 1.0e6
        ax.add_collection(LineCollection(
            np.stack((bound, motor), axis=1), colors=RED, linewidths=1.0,
            alpha=0.62, zorder=14,
        ))
        ax.scatter(bound[:, 0], bound[:, 1], s=18, facecolor="#fffefa",
                   edgecolor=RED, linewidth=1.2, zorder=15)
        force_nN = snapshot.clutch_forces * 1.0e9
        nonzero = np.linalg.norm(force_nN, axis=1) > 0.0
        if np.any(nonzero):
            ax.quiver(
                bound[nonzero, 0], bound[nonzero, 1],
                1.8 * force_nN[nonzero, 0], 1.8 * force_nN[nonzero, 1],
                angles="xy", scale_units="xy", scale=1.0, color=ORANGE,
                width=0.004, headwidth=4.0, headlength=5.0, alpha=0.88, zorder=16,
            )

    if show_displacement_vectors:
        if result.stage == "g3s":
            # A fixed, sparse set of the genuinely most-displaced material
            # points communicates deformation without a distracting arrow field.
            magnitude = np.linalg.norm(positions_um - initial_um, axis=1)
            radius = np.linalg.norm(initial_um - snapshot.cell_center[None, :] * 1e6, axis=1)
            magnitude[radius < result.config.cell_radius * 1e6 + 1.0] = -np.inf
            selected = np.argsort(magnitude)[-20:]
        else:
            bead_mask = np.isin(network.fiber_id, list(contact_ids | attached_ids))
            selected = np.flatnonzero(bead_mask)[::2]
        delta_um = positions_um[selected] - initial_um[selected]
        visible = np.linalg.norm(delta_um, axis=1) > 1.0e-7
        selected = selected[visible]
        delta_um = delta_um[visible]
        if selected.size:
            starts = initial_um[selected]
            ends = starts + displacement_scale * delta_um
            ax.add_collection(LineCollection(
                np.stack((starts, ends), axis=1), colors=GREEN,
                linewidths=0.75, alpha=0.72, zorder=8,
            ))

    trajectory = np.asarray(
        [item.cell_center for item in result.snapshots[:snapshot_index + 1]]
    ) * 1.0e6
    if trajectory.shape[0] > 1 and np.ptp(trajectory, axis=0).max(initial=0.0) > 0.0:
        ax.plot(trajectory[:, 0], trajectory[:, 1], color="#347bb0", lw=1.8, zorder=12)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    ax.grid(False)
    _add_scale_bar(ax, xlim, ylim)
    for spine in ax.spines.values():
        spine.set_color("#d8d2c5")
        spine.set_linewidth(0.8)


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
    initial_um = result.initial_positions * 1.0e6
    ax.set_facecolor(PANEL)
    for fibre in range(network.n_fibers):
        initial = initial_um[network.fiber_id == fibre]
        points = positions_um[network.fiber_id == fibre]
        ax.plot(initial[:, 0], initial[:, 1], color="#c8d0d3", lw=0.55, alpha=0.34)
        ax.plot(points[:, 0], points[:, 1], color="#647983", lw=0.82, alpha=0.88)
        ax.scatter(points[:, 0], points[:, 1], s=1.8, color="#647983", alpha=0.72)

    # Exact material-point crosslink locations for the G2-derived network.
    if network.xl_edge_a is not None and network.xl_edge_a.size:
        ea = network.fiber_bonds[network.xl_edge_a]
        eb = network.fiber_bonds[network.xl_edge_b]
        aa = network.xl_alpha_a[:, None]
        ab = network.xl_alpha_b[:, None]
        pa = (1.0 - aa) * positions_um[ea[:, 0]] + aa * positions_um[ea[:, 1]]
        pb = (1.0 - ab) * positions_um[eb[:, 0]] + ab * positions_um[eb[:, 1]]
        midpoint = 0.5 * (pa + pb)
        ax.scatter(midpoint[:, 0], midpoint[:, 1], marker="D", s=5.0,
                   color="#b89435", alpha=0.62, linewidths=0)

    fixed = result.fixture.fixed_mask
    if np.any(fixed):
        ax.scatter(positions_um[fixed, 0], positions_um[fixed, 1], marker="s", s=8,
                   color="#17242b", alpha=0.9, zorder=5)

    center_um = snapshot.cell_center * 1.0e6
    radius_um = result.config.cell_radius * 1.0e6
    ax.add_patch(Circle(center_um, radius_um, facecolor="#f6ded0", edgecolor="#c87453",
                        lw=1.7, alpha=0.96, zorder=6))
    body = center_um + 0.75 * radius_um * np.array([
        np.cos(snapshot.cell_angle), np.sin(snapshot.cell_angle)])
    ax.plot([center_um[0], body[0]], [center_um[1], body[1]], color="#7a4b00", lw=2.0)

    if snapshot.bound_points.size:
        bound = snapshot.bound_points * 1.0e6
        motor = snapshot.motor_points * 1.0e6
        segments = np.stack((bound, motor), axis=1)
        ax.add_collection(LineCollection(segments, colors="#d1495b", linewidths=0.8, alpha=0.56,
                                         zorder=9))
        ax.scatter(bound[:, 0], bound[:, 1], s=11, color="#d1495b", marker="o",
                   label="bound material points")

        # Display-scaled fibre force vectors, matching G2's explicit visual convention.
        magnitudes = np.linalg.norm(snapshot.clutch_forces, axis=1)
        nonzero = magnitudes > 0.0
        if np.any(nonzero):
            force_nN = snapshot.clutch_forces[nonzero] * 1.0e9
            starts = bound[nonzero]
            ends = starts + 1.8 * force_nN
            force_segments = np.stack((starts, ends), axis=1)
            ax.add_collection(LineCollection(force_segments, colors="#e67e22", linewidths=0.8,
                                             alpha=0.58, zorder=10))

    sector_angles = snapshot.cell_angle + result.protrusions.sector_angles
    sector_radius = 1.08 * radius_um
    sector_xy = center_um + sector_radius * np.column_stack((np.cos(sector_angles),
                                                              np.sin(sector_angles)))
    activity = snapshot.polarity_activity
    sizes = 8.0 + 75.0 * activity / max(float(activity.max(initial=0.0)), 1.0e-12)
    ax.scatter(sector_xy[:, 0], sector_xy[:, 1], c=activity, s=sizes, cmap="magma",
               vmin=0.0, vmax=max(float(activity.max(initial=0.0)), 1.0e-12), alpha=0.82,
               edgecolors="none", zorder=8)

    # Explicit protrusion shafts and tips. Inactive/retracting probes stay faint;
    # active probes are green, and contacted tips are linked by red clutch spokes.
    lengths = snapshot.protrusion_lengths * 1.0e6
    lab_angles = snapshot.cell_angle + result.protrusions.sector_angles
    bases = center_um + radius_um * np.column_stack((np.cos(lab_angles), np.sin(lab_angles)))
    tips = snapshot.protrusion_tips * 1.0e6
    visible = lengths > 0.02
    active_mask = np.zeros(lengths.size, dtype=bool)
    active_mask[snapshot.active_sectors] = True
    for sector in np.flatnonzero(visible):
        active_colour = "#2a9d78" if active_mask[sector] else "#9ac9b8"
        active_alpha = 0.98 if active_mask[sector] else 0.35
        ax.plot([bases[sector, 0], tips[sector, 0]], [bases[sector, 1], tips[sector, 1]],
                color=active_colour, lw=2.2 if active_mask[sector] else 1.0,
                alpha=active_alpha, zorder=8)
        ax.scatter(tips[sector, 0], tips[sector, 1], s=18 if active_mask[sector] else 8,
                   color=active_colour, alpha=active_alpha, zorder=9)
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
        fontsize=10, color=INK, fontweight="bold",
    )
    ax.grid(alpha=0.12, color="#7f8c8d")
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        unique = dict(zip(labels, handles))
        ax.legend(unique.values(), unique.keys(), loc="upper right", fontsize=7)
    if show_disclaimer:
        ax.text(0.5, -0.14, DISCLAIMER, transform=ax.transAxes, ha="center", va="top",
                fontsize=7, color="#555555")


def make_stage_animation(result: G3RunResult, path, fps: int = 10, max_frames: int = 180):
    """Write an attachment-to-network-response animation in G1/G2 card grammar."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    count = len(result.snapshots)
    indices = np.unique(np.linspace(0, count - 1, min(count, max_frames)).round().astype(int))
    # Slow the early causal sequence so growth and first capture are readable
    # instead of flashing past in a few tenths of a second.
    first_attached = next(
        (i for i, snapshot in enumerate(result.snapshots) if snapshot.bound_points.size),
        min(4, indices.size - 1),
    )
    early_cut = min(first_attached + 2, indices.size)
    indices = np.concatenate((np.repeat(indices[:early_cut], 3), indices[early_cut:]))
    # Match Gloria's G1/G2 panel grammar exactly: two white condition cards,
    # serif titles, border-only scenes, no numerical axes/grid, four metrics per card.
    fig = plt.figure(figsize=(13.1, 7.45))
    fig.patch.set_facecolor(PAPER)
    panel_lefts = (0.012, 0.508)
    for left in panel_lefts:
        fig.add_artist(FancyBboxPatch(
            (left, 0.012), 0.477, 0.976, transform=fig.transFigure,
            boxstyle="round,pad=0.008,rounding_size=0.018",
            facecolor="#fffefa", edgecolor="#d8d2c5", linewidth=0.9, zorder=-1,
        ))
    axes = np.asarray([
        fig.add_axes((0.027, 0.175, 0.447, 0.735)),
        fig.add_axes((0.523, 0.175, 0.447, 0.735)),
    ])
    center0 = result.snapshots[0].cell_center * 1.0e6
    radius_um = result.config.cell_radius * 1.0e6
    near_half = 1.85 * radius_um
    response_half = 3.5 * radius_um
    spheroid = result.stage == "g3s"
    question = (
        "Question: can a coarse surface shell load an irregular collagen mesh?"
        if spheroid else
        "Question: can one growing tip attach, then transmit a declared 5 nN load?"
    )
    left_title = "Surface-shell growth and attachment" if spheroid else "Protrusion growth and attachment"
    right_title = "G2-scale shell-driven fibre response" if spheroid else "G2-scale force-driven fibre response"
    # At G2's 5 nN protocol the physical response is tens of nanometres over
    # this short causal movie.  Show its direction legibly without pretending
    # the fibre geometry itself has moved 25 times farther; every metric card
    # remains in true units.
    display_scale = 25.0 if result.config.traction_target > 0.0 else 1000.0
    fig.text(0.5, 0.993, question, ha="center", va="top", fontsize=8.2,
             color="#5f6c74")
    fig.text(0.027, 0.958, left_title, ha="left", va="top",
             fontsize=18, fontweight="bold", family="serif", color=INK)
    fig.text(0.523, 0.958, right_title, ha="left", va="top",
             fontsize=18, fontweight="bold", family="serif", color=INK)
    fig.text(0.967, 0.927, f"displacement vectors ×{display_scale:g}" +
             " (visual aid; cards are true scale)",
             ha="right", va="top",
             fontsize=7.2, color="#5f6c74")
    cards = [_add_metric_cards(fig, left) for left in panel_lefts]

    def update(frame):
        index = int(indices[frame])
        snapshot = result.snapshots[index]
        center = snapshot.cell_center * 1.0e6
        _draw_gloria_frame(
            axes[0], result, index,
            (center[0] - near_half, center[0] + near_half),
            (center[1] - near_half, center[1] + near_half),
        )
        _draw_gloria_frame(
            axes[1], result, index,
            (center0[0] - response_half, center0[0] + response_half),
            (center0[1] - response_half, center0[1] + response_half),
            show_displacement_vectors=True, displacement_scale=display_scale,
        )
        force_nN = float(np.linalg.norm(snapshot.clutch_forces, axis=1).sum() * 1.0e9)
        attached = set(map(int, snapshot.bound_fibre_ids))
        rotation_mdeg = 1000.0 * _fibre_reorientation_deg(
            result, snapshot, attached if attached else _contact_fibre_ids(result)
        )
        max_displacement_um = float(np.linalg.norm(
            snapshot.positions - result.initial_positions, axis=1
        ).max(initial=0.0) * 1.0e6)
        left_values = [
            f"{_mechanism_phase(snapshot)}\nmechanism state",
            f"{snapshot.time:4.1f} s\ntime",
            f"{snapshot.protrusion_lengths.max(initial=0.0) * 1e6:5.2f} µm\n{'surface extension' if spheroid else 'protrusion length'}",
            f"{snapshot.bound_points.shape[0]}\nbound clutches",
        ]
        right_values = [
            f"{force_nN:6.3f} nN\nclutch load",
            f"{rotation_mdeg:7.3f} mdeg\ntrue max Δ angle",
            f"{max_displacement_um:.5f} µm\nmax bead Δr",
            f"{len(attached)}\nloaded fibres",
        ]
        for artist, value in zip(cards[0], left_values):
            artist.set_text(value)
        for artist, value in zip(cards[1], right_values):
            artist.set_text(value)

    animation = FuncAnimation(fig, update, frames=indices.size, interval=1000 / fps, repeat=True)
    if path.suffix.lower() == ".mp4":
        animation.save(path, writer=FFMpegWriter(fps=fps, bitrate=1800), dpi=120)
    else:
        animation.save(path, writer=PillowWriter(fps=fps), dpi=100)
    plt.close(fig)
    return path


def _spheroid_shell_metrics(result: G3RunResult, snapshot):
    """True displacement plus local radial collagen order in two annuli."""
    bonds = result.fixture.network.fiber_bonds
    a, b = snapshot.positions[bonds[:, 0]], snapshot.positions[bonds[:, 1]]
    segment = b - a
    length = np.linalg.norm(segment, axis=1)
    midpoint = 0.5 * (a + b)
    radius = np.linalg.norm(midpoint - snapshot.cell_center, axis=1)
    radial = (midpoint - snapshot.cell_center) / np.maximum(radius[:, None], 1e-30)
    tangent = segment / np.maximum(length[:, None], 1e-30)
    order = 2.0 * np.square(np.sum(tangent * radial, axis=1)) - 1.0
    shell_radius = result.config.cell_radius
    near = (radius >= shell_radius) & (radius < shell_radius + 15e-6)
    outer = (radius >= shell_radius + 15e-6) & (radius < shell_radius + 40e-6)
    value = lambda mask: float(np.mean(order[mask])) if np.any(mask) else float("nan")
    displacement = float(np.linalg.norm(snapshot.positions - result.initial_positions, axis=1).max(initial=0.0) * 1e6)
    return displacement, value(near), value(outer)


def make_spheroid_remodelling_animation(result: G3RunResult, path, fps: int = 8, max_frames: int = 24):
    """One full-field, one-question collective collagen-remodelling movie."""
    if result.stage != "g3s":
        raise ValueError("spheroid remodelling animation requires a g3s result")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    count = len(result.snapshots)
    indices = np.unique(np.linspace(0, count - 1, min(count, max_frames)).round().astype(int))
    indices = np.concatenate((np.repeat(indices[:min(3, len(indices))], 3), indices[3:]))
    fig = plt.figure(figsize=(12.8, 7.3), facecolor=PAPER)
    ax = fig.add_axes((0.035, 0.10, 0.67, 0.82))
    side = fig.add_axes((0.735, 0.10, 0.235, 0.82))
    fig.text(0.035, 0.975, "Collective spheroid traction remodels a random collagen field", ha="left", va="top",
             fontsize=17, fontweight="bold", family="serif", color=INK)
    fig.text(0.035, 0.942, "Question: does distributed surface traction create measurable collagen displacement and radial alignment?", ha="left", va="top", fontsize=8.5, color="#5f6c74")

    def update(frame):
        index = int(indices[frame])
        snapshot = result.snapshots[index]
        half = result.config.scaled_domain_size * 0.5e6
        _draw_gloria_frame(ax, result, index, (-half, half), (-half, half),
                            show_displacement_vectors=True, displacement_scale=5.0, full_field=True)
        side.clear(); side.set_axis_off(); side.set_facecolor("#fffefa")
        force = float(np.linalg.norm(snapshot.clutch_forces, axis=1).sum() * 1e9)
        engaged = len(set(map(int, snapshot.bound_sector_ids)))
        displacement, near_order, outer_order = _spheroid_shell_metrics(result, snapshot)
        rows = [
            (f"t = {snapshot.time:.1f} s", "time"),
            (f"{engaged} / {result.config.spheroid_surface_patches}", "engaged surface patches"),
            (f"{force:.1f} nN", "applied collective traction"),
            (f"{displacement:.3f} µm", "true maximum bead displacement"),
            (f"{near_order:+.3f}", "near-surface radial order"),
            (f"{outer_order:+.3f}", "outer-shell radial order"),
        ]
        side.text(0.0, 0.97, "Protocol", va="top", fontsize=13, fontweight="bold", family="serif", color=INK)
        side.text(0.0, 0.89, f"{result.config.spheroid_surface_patches} surface patches\n{result.config.cell_radius*1e6:.0f} µm radius spheroid\n{result.config.traction_target*1e9:.0f} nN target (sensitivity)", va="top", fontsize=8.5, color="#43535b")
        y = 0.69
        for value, label in rows:
            side.text(0.0, y, value, va="top", fontsize=12, fontweight="bold", color=INK)
            side.text(0.0, y - 0.038, label, va="top", fontsize=7.5, color="#5f6c74")
            y -= 0.105
        side.text(0.0, 0.02, "Dashed = initial mesh\nGreen vectors = 5× display aid\nAll values above are true scale.", va="bottom", fontsize=7.4, color="#5f6c74")

    animation = FuncAnimation(fig, update, frames=len(indices), interval=1000 / fps, repeat=True)
    animation.save(path, writer=PillowWriter(fps=fps), dpi=95)
    plt.close(fig)
    return path


def make_comparison_animation(
    results,
    path,
    fps: int = 10,
    max_frames: int = 180,
    panel_labels=None,
    title: str = "G3-R synchronized control · same scale and seed",
):
    """Write a synchronized Gloria-style G3B or G3C two-card comparison."""
    results = list(results)
    if not results:
        raise ValueError("comparison requires at least one result")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n_frames = min(max(len(result.snapshots) for result in results), max_frames)
    frame_fraction = np.linspace(0.0, 1.0, n_frames)
    is_motion = any(result.stage == "g3c" for result in results)
    if is_motion:
        common_bounds = _bounds(results[0])
        bounds = [common_bounds] * len(results)
    else:
        radius = results[0].config.cell_radius * 1.0e6
        half = 2.75 * radius
        bounds = []
        for result in results:
            centers = np.asarray([snapshot.cell_center for snapshot in result.snapshots]) * 1.0e6
            center = centers[0]
            bounds.append((center[0] - half, center[0] + half,
                           center[1] - half, center[1] + half))
    if len(results) != 2:
        raise ValueError("Gloria-format G3 comparisons require exactly two panels")
    fig = plt.figure(figsize=(13.1, 7.45))
    fig.patch.set_facecolor(PAPER)
    panel_lefts = (0.012, 0.508)
    for left in panel_lefts:
        fig.add_artist(FancyBboxPatch(
            (left, 0.012), 0.477, 0.976,
            transform=fig.transFigure, boxstyle="round,pad=0.008,rounding_size=0.018",
            facecolor="#fffefa", edgecolor="#d8d2c5", linewidth=0.9, zorder=-1,
        ))
    axes = np.asarray([
        fig.add_axes((0.027, 0.175, 0.447, 0.735)),
        fig.add_axes((0.523, 0.175, 0.447, 0.735)),
    ])
    fig.text(0.5, 0.993, title, ha="center", va="top", fontsize=8.2,
             color="#5f6c74")
    panel_labels = list(panel_labels) if panel_labels is not None else [None] * len(results)
    display_scale = 25.0 if any(result.config.traction_target > 0.0 for result in results) else 1000.0
    fig.text(0.967, 0.927, f"displacement vectors ×{display_scale:g}" +
             " (visual aid; cards are true scale)",
             ha="right", va="top", fontsize=7.2, color="#5f6c74")
    for panel, label in enumerate(panel_labels):
        if label:
            fig.text(panel_lefts[panel] + 0.015, 0.958, label, ha="left", va="top",
                     fontsize=18, fontweight="bold", family="serif", color=INK)
    cards = [_add_metric_cards(fig, left) for left in panel_lefts]

    def update(frame):
        for panel, (ax, result, bounds_one) in enumerate(zip(axes, results, bounds)):
            x0, x1, y0, y1 = bounds_one
            index = int(round(frame_fraction[frame] * (len(result.snapshots) - 1)))
            snapshot = result.snapshots[index]
            _draw_gloria_frame(
                ax, result, index, (x0, x1), (y0, y1), full_field=is_motion,
                show_displacement_vectors=True, displacement_scale=display_scale,
            )
            force_nN = float(np.linalg.norm(snapshot.clutch_forces, axis=1).sum() * 1.0e9)
            if is_motion:
                displacement = np.linalg.norm(
                    snapshot.cell_center - result.snapshots[0].cell_center
                ) * 1.0e6
                rotation = np.degrees(snapshot.cell_angle - result.snapshots[0].cell_angle)
                values = [
                    f"{displacement:.4f} µm\ncell Δr",
                    f"{rotation:.4f}°\ncell rotation",
                    f"{snapshot.bound_points.shape[0]}\nbound",
                    f"{force_nN:.3f} nN\nclutch load",
                ]
            else:
                active_degrees = np.degrees(
                    snapshot.cell_angle + result.protrusions.sector_angles[snapshot.active_sectors]
                )
                angle_text = ", ".join(f"{angle % 360:.0f}°" for angle in active_degrees)
                values = [
                    f"{angle_text}\nactive θ",
                    f"{_protrusion_switches(result, index)}\nsector switches",
                    f"{snapshot.bound_points.shape[0]}\nbound",
                    f"{force_nN:.3f} nN\nclutch load",
                ]
            for artist, value in zip(cards[panel], values):
                artist.set_text(value)

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
