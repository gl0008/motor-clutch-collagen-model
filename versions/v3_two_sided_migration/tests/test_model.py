import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from model import Config,counter_uniform,run_fixed_moving_pair

class TestV3(unittest.TestCase):
    def test_counter_stream(self):
        self.assertEqual(counter_uniform(7,10,1,3,0),counter_uniform(7,10,1,3,0))
        self.assertNotEqual(counter_uniform(7,10,1,3,0),counter_uniform(7,10,1,3,1))

    def test_fixed_and_moving_conditions(self):
        cfg=Config(duration=.3)
        pair=run_fixed_moving_pair(cfg)
        self.assertTrue(np.all(pair["fixed"]["cell_center"]==0))
        self.assertGreater(np.max(np.abs(pair["fixed"]["reaction"])),0)
        # Same stream means the initial zero-force stochastic histories agree.
        self.assertTrue(np.array_equal(pair["fixed"]["bound"][0],pair["moving"]["bound"][0]))

if __name__=="__main__":unittest.main()

