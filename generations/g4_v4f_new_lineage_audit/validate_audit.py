"""Validate the audit record without rerunning v4E or v4F simulations."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUDIT = json.loads((HERE / "lineage_audit.json").read_text())

EXPECTED_CHANGE_IDS = {
    "FNEW-REQ-001",
    "FNEW-REQ-002",
    "FNEW-REQ-003",
    "FNEW-HARM-001",
    "FNEW-NUM-001",
    "FNEW-NUM-002",
    "FNEW-INIT-001",
    "FNEW-NUM-003",
    "FNEW-PRES-001",
}


def manifest(path: Path) -> dict:
    text = path.read_text().strip()
    return json.loads(text.split("=", 1)[1].rstrip(";"))


def main() -> None:
    assert AUDIT["audit_type"] == "documentation_only"
    assert AUDIT["new_physics"] is False
    assert AUDIT["simulation_rerun"] is False
    assert AUDIT["comparison_verdict"]["parent_child_one_factor"] is False

    changes = AUDIT["changes"]
    ids = [item["id"] for item in changes]
    assert len(ids) == len(set(ids)), "duplicate audit change ID"
    assert set(ids) == EXPECTED_CHANGE_IDS, "missing or unexpected change entry"
    assert all(item["causal_boundary"] for item in changes)

    parent = AUDIT["lineage"]["parent_implementation"]
    child = AUDIT["lineage"]["child_implementation"]
    merge_base = subprocess.check_output(
        ["git", "merge-base", parent, child], cwd=ROOT, text=True
    ).strip()
    assert merge_base == parent, "recorded parent is not the child merge base"

    e = manifest(ROOT / "docs" / "g4-v4e-manifest.js")
    f = manifest(ROOT / "docs" / "g4-v4f-manifest.js")
    assert e["representativeSeed"] == 45
    assert f["representativeSeed"] == 51

    shared = set(e["config"]) & set(f["config"])
    differences = {
        key: (e["config"][key], f["config"][key])
        for key in shared
        if e["config"][key] != f["config"][key]
    }
    assert differences == {}, f"unexpected inherited config differences: {differences}"
    assert set(f["config"]) - set(e["config"]) == {"contactSearchInterval"}

    source = (ROOT / "generations" / "g4_v4f_corrected_ofat_guidance" / "model.py").read_text()
    required_markers = (
        "def _steric_force_pair",
        "def _dynamic_relocate_site",
        "max_gap = contact_width + probe_reach",
        "patches = _even_initial_sites",
        "Every sufficiently long-lived contact may become the front",
        "Explicit Euler must advance both members of the steric pair",
    )
    for marker in required_markers:
        assert marker in source, f"audit marker missing from v4F source: {marker}"

    print("G4 v4F_new lineage audit: PASS")
    print(f"  {len(changes)} changes classified")
    print("  ancestry, representative seeds and inherited manifest parameters verified")


if __name__ == "__main__":
    main()
