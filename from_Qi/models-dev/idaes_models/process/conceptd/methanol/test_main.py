from __future__ import division
from idaes_models.process.conceptd.methanol.main import do_LOA, build_model, manual_init_LOA
from six import itervalues
from idaes_models.core.util.misc import doNothing
from idaes_models.core.util.var_test import assert_var_equal
import unittest
from idaes_models.core.util.misc import category

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

class TestMethanol(unittest.TestCase):
    @category('frequent')
    def test_LOA(self):
        m = build_model()
        for o in itervalues(m.units):
            # print(o.name, o.type)
            getattr(o, 'apply_linear_relaxations', doNothing)()
            # do not use OA on interconnection nodes
            o.apply_OA_strategy(oa_ports=False)

        # m.solvers.local_NLP.options["outlev"] = 3  # use for troubleshooting CONOPT
        m.solvers.local_NLP.promote('ipopt')
        m.solvers.local_NLP.ipopt.tol = 1E-4
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")
        do_LOA(m)
        # LOA depends too heavily on choice of initialization. Testing this
        # with a different case.
        # assert_var_equal(self, m.local_LB, -2333, 1)
        # assert_var_equal(self, m.global_UB, -1727, 1)

    @category('frequent')
    def test_fixed_init_LOA(self):
        m = build_model()
        for o in itervalues(m.units):
            # print(o.name, o.type)
            getattr(o, 'apply_linear_relaxations', doNothing)()
            # do not use OA on interconnection nodes
            o.apply_OA_strategy(oa_ports=False)

        # m.solvers.local_NLP.options["outlev"] = 3  # use for troubleshooting CONOPT
        m.solvers.local_NLP.promote('ipopt')
        m.solvers.local_NLP.ipopt.tol = 1E-4
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")
        do_LOA(m, custom_init=manual_init_LOA)
        assert_var_equal(self, m.local_LB, -2744, 5)
        assert_var_equal(self, m.global_UB, -1575, 1)

    # def test_GLOA(self):
    #     pass
    #     # don't run GLOA for now until I figure out how to selectively run tests


if __name__ == '__main__':
    unittest.main()
