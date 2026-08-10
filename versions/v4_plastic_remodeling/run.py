"""Generate V4 permanent-versus-plastic load/unload comparison."""
import json
from pathlib import Path
import numpy as np
from model import Config,run_remodel_pair
def serial(v):
    if isinstance(v,np.ndarray):return np.round(v,5).tolist() if np.issubdtype(v.dtype,np.floating) else v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {k:serial(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [serial(x) for x in v]
    return v
if __name__=="__main__":
    here=Path(__file__).resolve().parent;data=run_remodel_pair(Config())
    (here/"demo"/"data.js").write_text("window.MODEL_DATA="+json.dumps(serial(data),separators=(",",":"))+";\n")
    p=data["plastic"];print(f"V4 broken={p['broken_total'][-1]}, reformed={p['reformed_total'][-1]}, residual dS={p['residual_alignment_change']:.5f}")
