import sys, unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from model import Config, active_forces, make_network, run


class TestV2(unittest.TestCase):
    def test_scale_and_contact_hybrid(self):
        cfg=Config(duration=.12)
        net=make_network(cfg)
        self.assertAlmostEqual(cfg.cell_radius,9.0)
        force,patches=active_forces(net,np.zeros(2),3.0)
        self.assertGreaterEqual(len(patches),3)
        self.assertAlmostEqual(sum(p.weight for p in patches),1.0,places=12)
        contacted=set(p.fiber for p in patches)
        direct=set(net.edge_fiber[np.flatnonzero(np.linalg.norm(force,axis=1)>0)])
        self.assertTrue(direct.issubset(contacted))
        self.assertTrue(all(0 <= p.surface_distance <= cfg.contact_width for p in patches))

    def test_crosslinks_are_permanent_and_anchors_fixed(self):
        cfg=Config(duration=.12)
        net=make_network(cfg); n=len(net.crosslinks); anchors=net.r[net.fixed].copy()
        result=run(cfg)
        self.assertGreater(n,0)
        self.assertEqual(len(result["crosslinks"]),n)
        self.assertTrue(np.allclose(result["positions"][-1][result["fixed"]],anchors))

    def test_no_active_force_no_motion(self):
        cfg=Config(total_pull_force=0,duration=.06)
        result=run(cfg)
        self.assertLess(np.max(np.abs(result["positions"][-1]-result["positions"][0])),1e-10)

    def test_crosslink_transmits_force_to_other_fiber(self):
        net=make_network(Config(duration=.06))
        link=net.crosslinks[0]
        a0,a1=net.edges[link.edge_a];b0,b1=net.edges[link.edge_b]
        net.r[a0] += np.array([.1,0])
        _,_,fx,_=net.elastic_forces(True)
        self.assertGreater(np.linalg.norm(fx[b0])+np.linalg.norm(fx[b1]),0)


if __name__=="__main__": unittest.main()
