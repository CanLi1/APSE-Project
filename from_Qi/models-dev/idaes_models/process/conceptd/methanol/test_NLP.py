from __future__ import division
from idaes_models.process.conceptd.methanol.main import build_model, solve_init_NLP
from idaes_models.core.util.var_test import assert_var_equal
import unittest
from idaes_models.core.util.misc import category

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

class TestMethanolNLP(unittest.TestCase):
    @category('frequent')
    def test_NLP_1(self):
        m = build_model()

        m.solvers.local_NLP.promote('ipopt')
        # m.solvers.local_NLP.outlev = 3  # use for troubleshooting CONOPT
        # m.solvers.local_NLP.max_iter = 3000
        m.solvers.local_NLP.tol = 1E-4
        # m.solvers.local_NLP.ipopt.linear_solver = 'ma57'
        solve_init_NLP(m, init_index=1, tee=False)

        assert_var_equal(self, m.global_UB, 859, 1)

    @category('frequent')
    def test_NLP_2(self):
        m = build_model()

        m.solvers.local_NLP.promote('ipopt')
        # m.solvers.local_NLP.outlev = 3  # use for troubleshooting CONOPT
        # m.solvers.local_NLP.max_iter = 5000
        m.solvers.local_NLP.tol = 1E-4
        solve_init_NLP(m, init_index=-1, tee=False)

        assert_var_equal(self, m.global_UB, -1575, 1)


if __name__ == '__main__':
    unittest.main()
