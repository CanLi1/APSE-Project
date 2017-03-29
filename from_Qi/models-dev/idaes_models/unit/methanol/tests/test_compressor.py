from __future__ import division
from idaes_models.unit.methanol.compressor import Compressor
from idaes_models.process.conceptd.methanol.flowsheet import MethanolModel
from six import itervalues
from idaes_models.core.util.misc import doNothing
from idaes_models.core.util.var_test import assert_var_equal
from pyomo.opt import SolverFactory
import unittest
from idaes_models.core.util.misc import category

class TestCompressor(unittest.TestCase):
    @category('frequent')
    def test_case_1(self):
        m = MethanolModel()
        m.comp_set = ['H2', 'CO', 'MeOH', 'CH4']
        m.max_flow = 20
        m.solvers.local_NLP.promote('ipopt')

        m.add_unit(
            Compressor(name='ms_feed_compr1', parent=m, compr_cost=25, bounds={'flow': (0, 5)},
                       block_bounds={'work': (1.7723, 6.6116),
                                     'P_ratio': (1.4165, 1.4341),
                                     'Pin': (None, 10),
                                     'Pout': (None, 115),
                                     'Tin': (None, 300),
                                     'Tout': (None, 500)}))
        m.build_units()
        m.units.ms_feed_compr1.equip_exists.fix(1)
        m._activate_standard_objective()
        for o in itervalues(m.units):
            getattr(o, 'apply_NLP', doNothing)()

        compressor = m.units.ms_feed_compr1

        compressor.flow['H2'].fix(2.8781)
        compressor.flow['CO'].fix(1.1992)
        compressor.flow['MeOH'].fix(0)
        compressor.flow['CH4'].fix(0.7195)
        compressor.Tin.fix(300)
        compressor.Pin.fix(10)
        compressor.Pout.fix(47.694)
        compressor.equip_exists.fix(1)

        m.dual.activate()
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")
        m.solve(using='local_NLP', tee=True, skip_trivial_constraints=True)

        tol = 1E-4

        assert_var_equal(self, compressor.Tout, 430.22, 1E-3)
        assert_var_equal(self, compressor.work, 2.5985, tol)
        assert_var_equal(self, compressor.P_ratio, 1.4341, tol)

    @category('frequent')
    def test_case_2(self):
        m = MethanolModel()
        m.comp_set = ['H2', 'CO', 'MeOH', 'CH4']
        m.max_flow = 20
        m.solvers.local_NLP.promote('ipopt')

        m.add_unit(
            Compressor(name='ms_feed_compr1', parent=m,
                       bounds={'flow': (0, 5), 'Pout': (25, None)},
                       block_bounds={'work': (3.3556, 3.7407),
                                     'P_ratio': (1.7400, 1.7400),
                                     'Pin': (None, 10),
                                     'Pout': (None, 115),
                                     'Tin': (None, 300),
                                     'Tout': (None, 600)
                                     }))
        m.build_units()
        m.units.ms_feed_compr1.equip_exists.fix(1)
        m._activate_standard_objective()
        for o in itervalues(m.units):
            getattr(o, 'apply_NLP', doNothing)()

        compressor = m.units.ms_feed_compr1

        compressor.flow['H2'].fix(2.3687)
        compressor.flow['CO'].fix(1.0932)
        compressor.flow['MeOH'].fix(0)
        compressor.flow['CH4'].fix(0.1822)
        compressor.Tin.fix(300)
        compressor.Pin.fix(10)
        compressor.Pout.fix(110.249)
        compressor.equip_exists.fix(1)

        m.dual.activate()
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")
        m.solve(using='local_NLP', tee=True, skip_trivial_constraints=True)

        tol = 1E-4

        assert_var_equal(self, compressor.Tout, 522.00, 1E-3)
        assert_var_equal(self, compressor.work, 3.3654, tol)
        assert_var_equal(self, compressor.P_ratio, 1.7400, tol)
