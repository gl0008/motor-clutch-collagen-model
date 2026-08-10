"""Generate V2 simulation data consumed by demo/index.html."""
import json
from pathlib import Path
import numpy as np
from model import Config, run


def serial(value):
    if isinstance(value, np.ndarray):
        return np.round(value,5).tolist() if np.issubdtype(value.dtype,np.floating) else value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: serial(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [serial(v) for v in value]
    return value


if __name__ == "__main__":
    here=Path(__file__).resolve().parent
    result=run(Config())
    (here/"demo"/"data.js").write_text("window.MODEL_DATA="+json.dumps(serial(result),separators=(",",":"))+";\n")
    print(f"V2: {len(result['fibers'])} fibers, {len(result['crosslinks'])} permanent crosslinks")
    print(f"alignment {result['alignment'][0]:.3f} -> {result['alignment'][-1]:.3f}")
