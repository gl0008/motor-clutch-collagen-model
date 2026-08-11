import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from model import Config,run_remodel_pair
class TestV4(unittest.TestCase):
    def test_control_never_breaks(self):
        cfg=Config(duration=.32,load_duration=.16,crosslink_force_threshold=0,crosslink_force_scale=.001,crosslink_off_rate0=20,crosslink_on_rate=.01)
        pair=run_remodel_pair(cfg)
        self.assertEqual(int(pair["elastic"]["broken_total"][-1]),0)
        self.assertGreater(int(pair["plastic"]["broken_total"][-1]),0)
    def test_same_initial_state(self):
        cfg=Config(duration=.08,load_duration=.04)
        pair=run_remodel_pair(cfg)
        self.assertTrue(np.allclose(pair["elastic"]["positions"][0],pair["plastic"]["positions"][0]))
        self.assertAlmostEqual(pair["elastic"]["alignment"][0],pair["plastic"]["alignment"][0])
if __name__=="__main__":unittest.main()

