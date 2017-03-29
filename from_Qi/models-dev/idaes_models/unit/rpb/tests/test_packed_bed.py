from __future__ import division
from idaes_models.process.conceptd.rpb.flowsheet import PackedBedModel
from idaes_models.unit.rpb.packed_bed import RotatingPackedBed, PackedBed
from six import itervalues
from idaes_models.core.util.misc import doNothing
from idaes_models.core.util.var_test import assert_var_equal
import unittest

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

class TestPackedBed(unittest.TestCase):
    def test_case_rpb(self):
        m = PackedBedModel()
        m.comp_set = ['A', 'water']
        m.max_flow = 10
        m.solvers.local_NLP.promote('ipopt')
        m.solvers.local_NLP.outlev = 3  # use for troubleshooting
        m.solvers.local_NLP.halt_on_ampl_error = "yes"
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")
        rpb = m.add_unit(
            RotatingPackedBed(name='rpb', parent=m, bounds={'frac_vap_in[A]': (1E-5, 1E-2), 'frac_vap_out[A]': (1E-5, 1E-2)}))
        m.build_units()
        m._activate_standard_objective()
        for o in itervalues(m.units):
            getattr(o, 'apply_NLP', doNothing)()

        rpb.Rout.fix(0.4542)
        rpb.Rin.fix(0.1514)
        rpb.total_flow_vap_in.fix(2.7709)
        rpb.frac_vap_in['A'].fix(0.01)
        rpb.h.fix(0.3431)
        rpb.speed.fix(13.3333)
        rpb.frac_vap_out.display()

        m.solve(using='local_NLP', tee=True, skip_trivial_constraints=True)

        rpb.display_variables()

        assert_var_equal(self, rpb.u, 0.01709, 5E-4)
        assert_var_equal(self, rpb.V, 0.198, 1E-3)
        assert_var_equal(self, rpb.total_flow_liq_in, 0.014, 5E-4)
        assert_var_equal(self, rpb.t, 0.19675, 1E-3)
        assert_var_equal(self, rpb.kL, 0.157, 1E-3)
        assert_var_equal(self, rpb.exp_base, 4.4797, 1E-3)
        assert_var_equal(self, rpb.frac_vap_out['A'], 0.000113, 1E-5)
        assert_var_equal(self, rpb.rotate_cost, 4.429, 1E-3)
        assert_var_equal(self, rpb.equip_cost, 33.33, 1E-2)

    def test_case_pb(self):
        m = PackedBedModel()
        m.comp_set = ['A', 'water']
        m.max_flow = 10
        m.solvers.local_NLP.promote('ipopt')
        m.solvers.local_NLP.outlev = 3  # use for troubleshootingCONOPT
        m.solvers.local_NLP.halt_on_ampl_error = "yes"
        # if not m.check_solvers_available():
        #     self.skipTest("Not all solvers available for superstructure synthesis.")
        pb = m.add_unit(
            PackedBed(name='pb', parent=m, bounds={'frac_vap_in[A]': (1E-5, 1E-2), 'frac_vap_out[A]': (1E-5, 1E-2)}))
        m.build_units()
        m._activate_standard_objective()
        for o in itervalues(m.units):
            getattr(o, 'apply_NLP', doNothing)()

        pb.R.fix(0.1)
        pb.total_flow_vap_in.fix(0.013)
        pb.frac_vap_in['A'].fix(0.008)
        pb.h.fix(1.515)
        pb.frac_vap_out.display()

        m.solve(using='local_NLP', tee=True, skip_trivial_constraints=True)

        pb.display_variables()

        assert_var_equal(self, pb.V, 0.048, 1E-3)
        assert_var_equal(self, pb.total_flow_liq_in, 1.3E-4, 1E-5)
        assert_var_equal(self, pb.exp_base, 4.68629, 1E-3)
        assert_var_equal(self, pb.frac_vap_out['A'], 7.368E-5, 5E-6)
        assert_var_equal(self, pb.equip_cost, 0.31488, 1E-2)


if __name__ == '__main__':
    unittest.main()
