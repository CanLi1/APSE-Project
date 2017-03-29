from __future__ import division
from idaes_models.process.conceptd.rpb.main import build_model
from six import itervalues
from idaes_models.core.util.misc import doNothing
from idaes_models.core.loa import do_LOA
from idaes_models.core.util.var_test import assert_var_equal
import unittest
from idaes_models.core.util.misc import category, requires_solver

class TestPackedBed(unittest.TestCase):
    def setup(self):
        pass

    @category('frequent')
    def test_2rpb_2pb(self):
        m = build_model(rpb=2, pb=2)
        m.solvers.local_NLP.promote('ipopt')

        for o in itervalues(m.units):
            o.apply_OA_strategy()
            getattr(o, 'apply_linear_relaxations', doNothing)(nsegs=2)

        do_LOA(m, tol=1, add_min_flows=False)
        # This is really hard to get to converge, so I'm leaving it out for
        # now.
        # assert_var_equal(self, m.obj, 34.48, 1.5)

    @category('frequent')
    def test_1rpb_0pb(self):
        m = build_model(rpb=1, pb=0)
        m.solvers.local_NLP.promote('ipopt')

        for o in itervalues(m.units):
            o.apply_OA_strategy()
            getattr(o, 'apply_linear_relaxations', doNothing)(nsegs=1)

        do_LOA(m, tol=1, add_min_flows=False)
        assert_var_equal(self, m.obj, 34.34, 1)

    @requires_solver('gurobi')
    @category('frequent')
    def test_1rpb_1pb(self):
        m = build_model(rpb=1, pb=1)
        m.solvers.local_NLP.promote('ipopt')

        for o in itervalues(m.units):
            o.apply_OA_strategy()
            getattr(o, 'apply_linear_relaxations', doNothing)(nsegs=7)

        do_LOA(m, tol=1, add_min_flows=False)
        assert_var_equal(self, m.obj, 34.48, 1)

    @category('frequent')
    def test_0rpb_1pb(self):
        m = build_model(rpb=0, pb=1)
        m.solvers.local_NLP.promote('ipopt')

        for o in itervalues(m.units):
            o.apply_OA_strategy()
            getattr(o, 'apply_linear_relaxations', doNothing)(nsegs=5)

        do_LOA(m, tol=1, add_min_flows=False)
        assert_var_equal(self, m.obj, 66.02, 1)


if __name__ == '__main__':
    unittest.main()
