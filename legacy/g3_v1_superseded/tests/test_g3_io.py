from dataclasses import replace

import json
import numpy as np

from g3.config import G3Config
from g3.run import save_run
from g3.simulation import run_g3


def test_save_run_writes_replayable_manifest_metrics_and_arrays(tmp_path):
    cfg = replace(G3Config(), n_clutches=10, n_motors=10, bind_rate=10.0,
                  metrics_interval=0.01, frame_interval=0.01)
    result = run_g3("g3a", "single_fibre", cfg, seed=5, duration=0.02)
    save_run(result, tmp_path, make_gif=False)
    for name in ("resolved_config.yaml", "metrics.json", "manifest.json",
                 "traces.npz", "frames.npz", "g3a_summary.png"):
        assert (tmp_path / name).exists()
    metrics = json.loads((tmp_path / "metrics.json").read_text(encoding="utf-8"))
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    traces = np.load(tmp_path / "traces.npz")
    assert metrics["seed"] == 5
    assert manifest["stage"] == "g3a"
    assert traces["time"][-1] == 0.02
