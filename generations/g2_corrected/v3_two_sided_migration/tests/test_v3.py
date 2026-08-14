from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from model import (  # noqa: E402
    MigrationConfig,
    _random_stream,
    run_migration,
    run_speed_ensemble,
)
sys.path.insert(0, str(HERE.parent))
from common.model import make_network_spec  # noqa: E402


class MigrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = MigrationConfig(duration=12.0, sample_interval=2.0)
        cls.spec = make_network_spec(cls.cfg)
        cls.stream = _random_stream(cls.cfg, 1234)

    def test_fixed_cell_remains_exactly_fixed(self):
        result = run_migration(
            self.cfg, spec=self.spec, moving=False, random_stream=self.stream
        )
        np.testing.assert_array_equal(result["cell_center"], 0.0)

    def test_common_stream_is_reproducible(self):
        a = _random_stream(self.cfg, 99)
        b = _random_stream(self.cfg, 99)
        np.testing.assert_array_equal(a[0], b[0])
        np.testing.assert_array_equal(a[1], b[1])

    def test_released_v3_cell_stays_on_protrusion_axis(self):
        result = run_migration(
            self.cfg, spec=self.spec, moving=True, random_stream=self.stream
        )
        np.testing.assert_array_equal(result["cell_center"][:, 1], 0.0)

    def test_mobility_calibration_has_expected_scale(self):
        short = MigrationConfig(duration=240.0, dt=0.1, sample_interval=6.0)
        speeds = run_speed_ensemble(short, trials=8)
        self.assertGreater(np.median(speeds), 0.05)
        self.assertLess(np.median(speeds), 0.8)


if __name__ == "__main__":
    unittest.main()
