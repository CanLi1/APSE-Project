from __future__ import division
from idaes_models.process.conceptd.methanol.flowsheet import MethanolModel
from idaes_models.unit.methanol.flash import Flash
from six import itervalues
from idaes_models.core.util.misc import doNothing
from idaes_models.core.util.var_test import assert_var_equal
import unittest
from pyomo.core.base.connector import ConnectorExpander
from idaes_models.core.util.misc import category, requires_solver

class TestFlash(unittest.TestCase):
    @requires_solver('baron')
    @category('frequent')
    def test_case_1(self):
        m = MethanolModel()
        m.comp_set = ['H2', 'CO', 'MeOH', 'CH4']
        m.max_flow = 20
        m.min_T = 300
        m.max_T = 900
        m.max_P = 150
        antoine = {
            'A': {'H2': 13.6333, 'CO': 14.3686, 'MeOH': 18.5875, 'CH4': 15.2243},
            'B': {'H2': 164.9, 'CO': 530.22, 'MeOH': 3626.55, 'CH4': 897.84},
            'C': {'H2': 3.19, 'CO': -13.15, 'MeOH': -34.29, 'CH4': -7.16}}

        m.add_unit(
            Flash(name='flash', parent=m, antoine=antoine, key_comp='H2',
                  bounds={'vap_recover[MeOH]': (None, 0.9),
                          'flow_liq[MeOH]': (0.09, 0.9)}))

        m.build_units()

        m._activate_standard_objective()
        for o in itervalues(m.units):
            getattr(o, 'apply_NLP', doNothing)()
        m.dual.deactivate()

        flash = m.units.flash
        flash.flow_in['H2'].fix(5.6888)
        flash.flow_in['CO'].fix(0.9082)
        flash.flow_in['MeOH'].fix(2.4398)
        flash.flow_in['CH4'].fix(5.5818)
        flash.T.fix(326.39)
        flash.P.fix(53.885)

        xfrm = ConnectorExpander()
        xfrm.apply(instance=m)

        # m.solvers.local_NLP.promote('ipopt')
        # m.solvers.local_NLP.options["outlev"] = 3  # use for troubleshooting CONOPT
        m.solve(using='global_NLP')

        tol = 1E-4

        print("Nonlinear flash tests, case 1")
        assert_var_equal(self, flash.flow_vap['H2'], 5.6553, tol)
        assert_var_equal(self, flash.flow_vap['H2'], 5.6553, tol)
        assert_var_equal(self, flash.flow_vap['CO'], 0.9009, tol)
        assert_var_equal(self, flash.flow_vap['MeOH'], 1.5398, 1E-3)
        assert_var_equal(self, flash.flow_vap['CH4'], 5.5226, tol)
        assert_var_equal(self, flash.flow_liq['H2'], 0.0335, tol)
        assert_var_equal(self, flash.flow_liq['CO'], 0.0073, tol)
        assert_var_equal(self, flash.flow_liq['MeOH'], 0.9000, 1E-3)
        assert_var_equal(self, flash.flow_liq['CH4'], 0.0591, tol)
        assert_var_equal(self, flash.total_flow_liq, 1.0000, 1E-3)

    @requires_solver('baron')
    @category('frequent')
    def test_case_2(self):
        m = MethanolModel()
        m.comp_set = ['H2', 'CO', 'MeOH', 'CH4']
        m.max_flow = 20
        m.min_T = 300
        m.max_T = 900
        m.max_P = 150
        antoine = {
            'A': {'H2': 13.6333, 'CO': 14.3686, 'MeOH': 18.5875, 'CH4': 15.2243},
            'B': {'H2': 164.9, 'CO': 530.22, 'MeOH': 3626.55, 'CH4': 897.84},
            'C': {'H2': 3.19, 'CO': -13.15, 'MeOH': -34.29, 'CH4': -7.16}}

        m.add_unit(
            Flash(name='flash', parent=m, antoine=antoine, key_comp='H2',
                  bounds={'vap_recover[MeOH]': (None, 0.9),
                          'flow_liq[MeOH]': (0.09, 0.9)}))

        m.build_units()

        m._activate_standard_objective()
        for o in itervalues(m.units):
            getattr(o, 'apply_NLP', doNothing)()
        m.dual.deactivate()

        flash = m.units.flash
        flash.flow_in['H2'].fix(4.13842899891)
        flash.flow_in['CO'].fix(0.817359559226)
        flash.flow_in['MeOH'].fix(2.96407656406)
        flash.flow_in['CH4'].fix(2.30922044107)
        flash.T.fix(369.053318103)
        flash.P.fix(78.9202114431)

        xfrm = ConnectorExpander()
        xfrm.apply(instance=m)

        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")
        m.solve(using='global_NLP')

        tol = 1E-4

        print("Nonlinear flash tests, case 2")
        assert_var_equal(self, flash.flow_vap['H2'], 4.08904071212, tol)
        assert_var_equal(self, flash.flow_vap['CO'], 0.805462144237, tol)
        assert_var_equal(self, flash.flow_vap['MeOH'], 2.06407656512, 1E-3)
        assert_var_equal(self, flash.flow_vap['CH4'], 2.27050613187, tol)
        assert_var_equal(self, flash.flow_liq['H2'], 0.049388286785, tol)
        assert_var_equal(self, flash.flow_liq['CO'], 0.0118974149893, tol)
        assert_var_equal(self, flash.flow_liq['MeOH'], 0.899999998939, 1E-3)
        assert_var_equal(self, flash.flow_liq['CH4'], 0.0387143091988, tol)
        assert_var_equal(self, flash.total_flow_liq, 1.0000, 1E-3)
