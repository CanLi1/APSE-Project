from __future__ import division
from idaes_models.process.conceptd.water.main import build_model
from six import itervalues
from idaes_models.core.util.misc import doNothing
from idaes_models.core.util.var_test import assert_var_equal
import unittest
from idaes_models.core.util.misc import category, requires_solver
from idaes_models.core.plugins.propagate_fixed import propagate_var_fix, reset_propagated_var_fix
from pyomo.environ import value
from pyomo.opt import TerminationCondition


class TestWater(unittest.TestCase):
    @requires_solver('gurobi')
    @category('frequent')
    def test_LOA_NLP(self):
        m = build_model(trial=None)
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")
        for o in itervalues(m.units):
            getattr(o, 'apply_linear_relaxations', doNothing)()
            getattr(o, 'apply_OA_strategy', doNothing)(oa_ports=True)

        # m.solvers.local_NLP.options["outlev"] = 3  # use for troubleshooting CONOPT
        m.solvers.local_NLP.promote('ipopt')
        # m.solvers.local_NLP.outlev = 3  # use for troubleshooting CONOPT
        # m.solvers.local_NLP.halt_on_ampl_error = 'yes'
        m.solvers.local_NLP.ipopt.tol = 1E-4
        m.solvers.local_NLP.max_iter = 1000
        # do_LOA(m)
        m.solve_local_MIP()
        # m.units.in1.outlet.display_flows()
        # m.units.in2.outlet.display_flows()
        # m.units.in3.outlet.display_flows()
        # m.units.in4.outlet.display_flows()
        # m.units.in5.outlet.display_flows()
        # m.units.tru2.outlet.out_streams.display()
        # m.units.tru2.outlet.lin_cuts.mass_balance.display()
        m._activate_standard_objective()
        for o in itervalues(m.units):
            o.apply_NLP()
        m.dual.activate()
        propagate_var_fix(m, tmp=True)
        for o in itervalues(m.units):
            o.introspect_flows()
        for o in itervalues(m.units):
            o.deactivate_trivial_constraints()
        for o in itervalues(m.units):
            o.set_min_flows()
        # m.units.tru2.outlet.display()
        # print(m.units.tru2.outlet.equip.flow_comps_unit.active)
        # m.units.in1.outlet.equip.flow_comps_out.pprint()
        # print(m.units.in1.outlet.equip.flow_comps_out['tru1', '1'].active)
        # m.units.in1.outlet.fc_tru4.display()
        # m.units.tru4.inlet.fc_in1.display()
        # m.units.tru4.inlet.fc_unit.display()
        # m.units.tru4.flow_in.display()
        # getattr(m.units.tru4.links, 'var_link_with_inlet.expanded').pprint()
        results = m.solve(using='local_NLP', tee=False, skip_trivial_constraints=True)
        self.assertIs(results.solver.termination_condition, TerminationCondition.optimal)
        m._global_UB = min([m.global_UB, value(m.obj.expr)])
        # print(results)
        # print(m._tmp_propagate_fixed)
        reset_propagated_var_fix(m)
        for o in itervalues(m.units):
            o.reset_introspect_fixed()
            o.reset_trivial_constraints()
            o.reset_min_flows()
        m.add_oa_cut()
        m.dual.deactivate()
        # m.add_integer_cut(tmp=True)
        m._activate_standard_objective()
        m._activate_penalized_oa_objective()
        for o in itervalues(m.units):
            o.apply_MIP()
        m.solve_local_MIP()
        # m.units.in1.outlet.display_flows()
        # m.units.in2.outlet.display_flows()
        # m.units.in3.outlet.display_flows()
        # m.units.in4.outlet.frac.display()
        # m.units.in4.outlet.lin_cuts.mc_exit_flows.display()
        # m.units.in5.outlet.display_flows()

        # assert_var_equal(self, m.local_LB, 301255, 100)
        # assert_var_equal(self, m.global_UB, 349525, 100)


if __name__ == '__main__':
    unittest.main()
