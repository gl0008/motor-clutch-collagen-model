"""Compact JavaScript export for precomputed browser animations."""

from __future__ import annotations

import base64
import json
from pathlib import Path
from typing import Any

import numpy as np


def _plain(value: Any):
    if isinstance(value, np.ndarray):
        if np.issubdtype(value.dtype, np.floating):
            return np.round(value, 5).tolist()
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


def quantize_positions(result: dict, quantum: float = 0.0005) -> dict:
    """Replace frame positions by base64 int16 displacements from frame zero."""

    result = dict(result)
    positions = np.asarray(result.pop("positions"), dtype=float)
    initial = np.asarray(result["initial_positions"], dtype=float)
    delta = np.rint((positions - initial[None, :, :]) / quantum)
    maximum = float(np.max(np.abs(delta))) if delta.size else 0.0
    if maximum > np.iinfo(np.int16).max:
        raise ValueError(
            f"position quantisation overflow ({maximum}); increase quantum"
        )
    packed = np.asarray(delta, dtype="<i2").tobytes(order="C")
    result["position_encoding"] = {
        "dtype": "int16-le",
        "quantum_um": quantum,
        "shape": list(delta.shape),
        "base64": base64.b64encode(packed).decode("ascii"),
    }
    return result


def write_data_js(path: Path, data: dict, *, quantize_keys=()) -> None:
    payload = dict(data)
    for key in quantize_keys:
        payload[key] = quantize_positions(payload[key])
    text = "window.MODEL_DATA=" + json.dumps(
        _plain(payload), separators=(",", ":"), allow_nan=False
    ) + ";\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
