import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from model import Config,gaussian_forces,make_network
class TestV1(unittest.TestCase):
    def test_weights_and_fixed_ends(self):
        cfg=Config(duration=.04);net=make_network(cfg);anchors=net.r[net.fixed].copy()
        active,w=gaussian_forces(net,cfg,1.0)
        self.assertAlmostEqual(float(w.sum()),1.0)
        net.advance(active,cfg.dt)
        self.assertTrue(np.allclose(net.r[net.fixed],anchors))
    def test_no_sls_or_crosslinks(self):
        net=make_network(Config())
        self.assertFalse(hasattr(net,"crosslinks"));self.assertFalse(hasattr(net,"maxwell_force"))
if __name__=="__main__":unittest.main()
