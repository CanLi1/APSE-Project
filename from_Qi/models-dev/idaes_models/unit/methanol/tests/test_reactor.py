from __future__ import division
from idaes_models.process.conceptd.methanol.flowsheet import MethanolModel
from idaes_models.unit.methanol.reactor import Reactor
from six import itervalues
from idaes_models.core.util.misc import doNothing
from idaes_models.core.util.var_test import assert_var_equal
import unittest
from idaes_models.core.util.misc import category

class TestReactor(unittest.TestCase):
    @category('frequent')
    def test_case_1(self):
        m = MethanolModel()
        m.comp_set = ['H2', 'CO', 'MeOH', 'CH4']
        m.max_flow = 20
        rxn_stoich = {'H2': -1.0, 'CO': -0.5, 'MeOH': 0.5, 'CH4': 0.0}

        m.add_unit(
            Reactor(name='expensive_rxn', parent=m, stoich=rxn_stoich, key_comp='H2', v_param=0.1,
                    vol_cost=10, rxn_cost=250, bounds={
                        'conv': (0, 0.415), 'eq_conv': (0, 0.415), 'extent': (0, 5),
                        'flow_in': (0, 10), 'Pin': (0, 150)},
                    block_bounds={'flow_out[CO]': (0.1, None)}
                    ))

        m.build_units()

        m.units.expensive_rxn.equip_exists.fix(1)
        m._activate_standard_objective()
        for o in itervalues(m.units):
            getattr(o, 'apply_NLP', doNothing)()

        erxn = m.units.expensive_rxn

        erxn.flow_in['H2'].fix(7.8571)
        erxn.flow_in['CO'].fix(1.9924)
        erxn.flow_in['MeOH'].fix(1.3557)
        erxn.flow_in['CH4'].fix(5.5818)
        # model.reactor.Tin.fix(455.39)  # needed to round this original number for feasibility
        erxn.Tin.fix(455.393)
        erxn.Pin.fix(150)

        m.dual.activate()
        # m.solvers.local_NLP.promote('ipopt')
        # m.solvers.local_NLP.outlev = 3  # use for troubleshooting
        m.solvers.local_NLP.halt_on_ampl_error = "yes"
        m.solve(using='local_NLP', tee=True, skip_trivial_constraints=True)

        tol = 1E-4

        assert_var_equal(self, erxn.Pout, 135, tol),
        assert_var_equal(self, erxn.Tout, 523, 1E-3),
        assert_var_equal(self, erxn.flow_out['H2'], 5.6888, tol),
        assert_var_equal(self, erxn.flow_out['CO'], 0.9082, tol),
        assert_var_equal(self, erxn.flow_out['MeOH'], 2.4398, tol),
        assert_var_equal(self, erxn.flow_out['CH4'], 5.5818, tol),
        assert_var_equal(self, erxn.total_flow_in, 16.7869, tol),
        assert_var_equal(self, erxn.total_flow_out, 14.6187, tol),
        assert_var_equal(self, erxn.eq_conv, 0.4135, tol),
        assert_var_equal(self, erxn.conv, 0.2760, tol),
        assert_var_equal(self, erxn.extent, 2.1683, tol)
