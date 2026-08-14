"""Precompute the corrected V2 no-crosslink / crosslinked comparison."""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from common.model import CollagenConfig, run_fixed_pull_pair  # noqa: E402
from common.serialization import write_data_js  # noqa: E402


if __name__ == "__main__":
    result = run_fixed_pull_pair(CollagenConfig())
    write_data_js(
        HERE / "demo" / "data.js",
        result,
        quantize_keys=("without_crosslinks", "with_crosslinks"),
    )
    free = result["without_crosslinks"]
    linked = result["with_crosslinks"]
    print(
        "Generation 2 / V2:",
        linked["derived"]["bead_count"],
        "beads,",
        linked["derived"]["crosslink_count"],
        "permanent crosslinks",
    )
    print("near/intermediate/far displacement without links:", free["displacement_by_shell"][-1])
    print("near/intermediate/far displacement with links:   ", linked["displacement_by_shell"][-1])
