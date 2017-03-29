from __future__ import division
from idaes_models.unit.clc.steam_rxn_cocur import SteamRxnCocur
from idaes_models.unit.clc.steam_rxn_counter import SteamRxnCounter
from idaes_models.process.conceptd.clc.flowsheet import CLCFlowsheet
import unittest
from idaes_models.core.util.misc import category
from pyomo.opt import TerminationCondition
# from idaes_models.core.util.var_test import assert_var_equal

class TestSteamRxnCocur(unittest.TestCase):
    """Tests for steam reactors
    """
    @category('frequent')
    def test_cocurrent(self):
        """Tests the co-current steam reactor
        """
        m = CLCFlowsheet()
        m.comp_set = ['C', 'CH4', 'CO', 'H2', 'CO2', 'H2O', 'Fe2O3',
                      'Fe3O4', 'FeO', 'Fe', 'O2', 'N2']
        m.max_flow = 300
        m.add_unit(SteamRxnCocur(name="steam_reactor", parent=m))

        # TODO finish test. Fix relevant co-current stream reactor variables
        # here in order to specify inputs

        results = m.units.steam_reactor.solve(tee=True, options={'outlev': 5, 'tol': 1e-8})
        for s in m.units.steam_reactor.k:
            getattr(m.units.steam_reactor, 'y_' + s).pprint()
        m.units.steam_reactor.T.pprint()
        m.units.steam_reactor.F.pprint()
        m.units.steam_reactor.Hc.pprint()
        m.units.steam_reactor.H.pprint()
        m.units.steam_reactor.Q.pprint()
        self.assertIs(results.solver.termination_condition, TerminationCondition.optimal)

        # TODO add assert statements. Add assert statements here in order to
        # check that the model variables are at the expected values. Example
        # below:

        # assert_var_equal(self, m.units.air_reactor.F['SolidOut'], 122, 1)

    @category('frequent')
    def test_countercurrent(self):
        """Tests the counter-current steam reactor
        """
        m = CLCFlowsheet()
        m.comp_set = ['C', 'CH4', 'CO', 'H2', 'CO2', 'H2O', 'Fe2O3',
                      'Fe3O4', 'FeO', 'Fe', 'O2', 'N2']
        m.max_flow = 300
        m.add_unit(SteamRxnCounter(name="steam_reactor", parent=m))

        # TODO finish test. Fix relevant co-current stream reactor variables
        # here in order to specify inputs

        results = m.units.steam_reactor.solve(tee=True, options={'outlev': 5, 'tol': 1e-8})
        for s in m.units.steam_reactor.k:
            getattr(m.units.steam_reactor, 'y_' + s).pprint()
        m.units.steam_reactor.T.pprint()
        m.units.steam_reactor.F.pprint()
        m.units.steam_reactor.Hc.pprint()
        m.units.steam_reactor.H.pprint()
        m.units.steam_reactor.Q.pprint()
        self.assertIs(results.solver.termination_condition, TerminationCondition.optimal)

        # TODO add assert statements. Add assert statements here in order to
        # check that the model variables are at the expected values. Example
        # below:

        # assert_var_equal(self, m.units.air_reactor.F['SolidOut'], 122, 1)


if __name__ == "__main__":
    unittest.main()
