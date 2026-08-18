"""Re-render saved G3 G2-scale trajectories without changing their physics."""

from __future__ import annotations

from pathlib import Path
import shutil

from g3.run import load_saved_run
from g3.visualization import make_comparison_animation, make_stage_animation


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "g3_g2scale"
ASSETS = ROOT / "docs" / "assets" / "g3_g2scale"


def main():
    runs = {
        name: load_saved_run(RESULTS / name)
        for name in (
            "g3a_5nN", "g3b_feedback_off_5nN", "g3b_feedback_on_5nN",
            "g3c_fixed_5nN", "g3c_released_5nN", "g3s_surface_shell_20nN",
        )
    }
    ASSETS.mkdir(parents=True, exist_ok=True)
    renders = {
        "g3a_g2scale_loading.gif": lambda path: make_stage_animation(runs["g3a_5nN"], path, fps=8, max_frames=18),
        "g3b_g2scale_feedback.gif": lambda path: make_comparison_animation(
            [runs["g3b_feedback_off_5nN"], runs["g3b_feedback_on_5nN"]], path, fps=8, max_frames=18,
            panel_labels=["Intrinsic polarity only", "+ adhesion feedback"],
            title="G3B diagnostic · no OFF/ON separation in this short matched run",
        ),
        "g3c_g2scale_released.gif": lambda path: make_comparison_animation(
            [runs["g3c_fixed_5nN"], runs["g3c_released_5nN"]], path, fps=8, max_frames=18,
            panel_labels=["Rigid cell (fixed)", "Released rigid cell"],
            title="G3C question · does releasing the cell change reaction-driven motion?",
        ),
        "g3s_surface_shell_20nN.gif": lambda path: make_stage_animation(
            runs["g3s_surface_shell_20nN"], path, fps=8, max_frames=18),
    }
    for filename, render in renders.items():
        path = RESULTS / filename
        render(path)
        shutil.copy2(path, ASSETS / filename)


if __name__ == "__main__":
    main()
