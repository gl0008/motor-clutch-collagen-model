"""Generate synchronized fixed-versus-moving V3 data."""
import json
from pathlib import Path
import numpy as np
from model import Config, run_fixed_moving_pair

def serial(v):
    if isinstance(v,np.ndarray):return np.round(v,5).tolist() if np.issubdtype(v.dtype,np.floating) else v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {k:serial(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [serial(x) for x in v]
    return v

if __name__=="__main__":
    here=Path(__file__).resolve().parent;data=run_fixed_moving_pair(Config())
    (here/"demo"/"data.js").write_text("window.MODEL_DATA="+json.dumps(serial(data),separators=(",",":"))+";\n")
    dx=data["moving"]["cell_center"][-1,0]-data["moving"]["cell_center"][0,0]
    print(f"V3 final moving-cell dx = {dx:.5f} um; fixed-cell dx = 0")
