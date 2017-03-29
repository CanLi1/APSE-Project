from __future__ import division
from idaes_models.process.conceptd.water.main import do_LOA, build_model, do_GLOA
from six import itervalues
from idaes_models.core.util.misc import doNothing
from idaes_models.core.util.var_test import assert_var_equal
import unittest
from idaes_models.core.util.misc import category


class TestWater(unittest.TestCase):
    @category('frequent')
    def test_default_LOA(self):
        m = build_model(trial=None)
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")
        for o in itervalues(m.units):
            getattr(o, 'apply_linear_relaxations', doNothing)()
            getattr(o, 'apply_OA_strategy', doNothing)(oa_ports=True)

        # m.solvers.local_NLP.options["outlev"] = 3  # use for troubleshooting CONOPT
        m.solvers.local_NLP.promote('conopt')
        m.solvers.local_NLP.ipopt.tol = 1E-4
        do_LOA(m, add_min_flows=False)

        assert_var_equal(self, m.local_LB, 301255, 100)
        assert_var_equal(self, m.global_UB, 349525, 100)

    @category('frequent')
    def test_trial01_LOA(self):
        m = build_model(trial='01')
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")
        for o in itervalues(m.units):
            getattr(o, 'apply_linear_relaxations', doNothing)()
            getattr(o, 'apply_OA_strategy', doNothing)(oa_ports=True)

        m.solvers.local_NLP.promote('conopt')
        do_LOA(m)

    @category('skip-broken')
    def test_GLOA_one_iter(self):
        m = build_model(trial=None)
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")

        for o in itervalues(m.units):
            o.apply_linear_relaxations()
            o.apply_OA_strategy()

        m.solvers.local_NLP.promote('conopt')
        do_LOA(m, tol=100)
        do_GLOA(m, tol=100, iterlim=1, contract_bounds=True, do_self_proj=True)

    @category('nightly')
    def test_GLOA_bounding(self):
        m = build_model(trial=None)
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")

        for o in itervalues(m.units):
            o.apply_linear_relaxations()
            o.apply_OA_strategy()

        m.solvers.local_NLP.promote('conopt')
        do_LOA(m, tol=100)
        do_GLOA(m, tol=100, iterlim=32, contract_bounds=True, do_self_proj=False)
        assert_var_equal(self, m.global_UB, 348337, 150)
        self.assertTrue(m.global_LB + 100 >= m.global_UB,
                        msg="Lower bound did not converge")

    @category('skip-broken')
    def test_GLOA_full(self):
        m = build_model(trial=None)
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")

        for o in itervalues(m.units):
            o.apply_linear_relaxations()
            o.apply_OA_strategy()

        m.solvers.local_NLP.promote('conopt')
        do_LOA(m, tol=100)
        do_GLOA(m, tol=100, iterlim=32, contract_bounds=True, do_self_proj=True)


if __name__ == '__main__':
    unittest.main()
