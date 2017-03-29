from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import Var, Constraint, NonNegativeReals, Param, Block
from idaes_models.core.util.var import none_if_empty as ne
from idaes_models.core.ports import InPort, OutPort


@ProcBlock("Valve")
class _Valve(UnitModel):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('type', [])
        kwargs['type'].append('Valve')
        self._gamma = kwargs.pop('gamma')
        self.linearize = kwargs.pop('linearize', 'OA')
        super(_Valve, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        super(_Valve, self)._build_unit_sets()

    def _build_unit_params(self):
        super(_Valve, self)._build_unit_params()
        b = self
        b.gamma = Param(initialize=self._gamma)

    def _build_unit_vars(self):
        super(_Valve, self)._build_unit_vars()
        b = self

        b.Tin = Var(domain=NonNegativeReals, initialize=500, bounds=(300, 900))
        b.Tout = Var(domain=NonNegativeReals, initialize=500, bounds=(300, 900))
        b.Pin = Var(domain=NonNegativeReals, initialize=15, bounds=(1, 150))
        b.Pout = Var(domain=NonNegativeReals, initialize=9.5, bounds=(1, 150))
        b.flow = Var(self.comps, domain=NonNegativeReals, initialize=1.0, bounds=(0, 20))
        b.P_ratio = Var(domain=NonNegativeReals)
        b.P_rat = Var(domain=NonNegativeReals)

        # Calculated bounds
        b.P_ratio.setlb((b.Pout.lb / b.Pin.ub) ** b.gamma)
        b.P_ratio.setub(min(b.Pout.ub / b.Pin.lb, 1) ** b.gamma)

        # Calculated initialization
        b.P_ratio.value = (b.Pout.value / b.Pin.value) ** b.gamma

    def _build_unit_constraints(self):
        super(_Valve, self)._build_unit_constraints()
        b = self

        equip = b.equip

        def exit_P(equip, *vs):
            b = equip.parent_block()
            return b.Pout[ne(vs)] <= b.Pin[ne(vs)]
        equip.exit_P = Constraint(rule=exit_P)

        # -------------------------
        # - Nonlinear constraints -
        # -------------------------

        def pressure_ratio(equip, *vs):
            b = equip.parent_block()
            return b.Pin[ne(vs)] ** b.gamma * b.P_ratio[ne(vs)] == b.Pout[ne(vs)] ** b.gamma
        equip.pressure_ratio = Constraint(rule=pressure_ratio)

        def exit_T(equip, *vs):
            b = equip.parent_block()
            return b.Tout[ne(vs)] == b.Tin[ne(vs)] * b.P_ratio[ne(vs)]
        equip.exit_T = Constraint(rule=exit_T)

    def _build_unit_ports(self):
        super(_Valve, self)._build_unit_ports()
        b = self
        self.add_unit(InPort(name='inlet', parent=self, fc=b.flow, T=b.Tin, P=b.Pin))
        self.add_unit(OutPort(name='outlet', parent=self, fc=b.flow, T=b.Tout, P=b.Pout))

    def get_flow_vars(self):
        yield self.flow

    def apply_linear_relaxations(self, nsegs=1):
        b = self
        lin_cuts = b.lin_cuts
        from idaes_models.core.util.mccormick import add_mccormick_relaxation
        from idaes_models.core.util.concave import add_concave_linear_underest
        super(_Valve, self).apply_linear_relaxations(nsegs=nsegs)

        lin_cuts.mc_Prat = Block()
        b.P_rat.setlb(b.Pout.lb / b.Pin.ub)
        b.P_rat.setub(min(b.Pout.ub / b.Pin.lb, 1))
        add_mccormick_relaxation(lin_cuts.mc_Prat, b.Pout, b.Pin, b.P_rat, nsegs, None, 1.0)

        def concave_Pratio(P_rat):
            return P_rat ** b.gamma
        lin_cuts.env_Pratio = Block()
        add_concave_linear_underest(lin_cuts.env_Pratio, 'concave_Pratio', 1, b.P_rat, b.P_ratio, concave_Pratio)
        # also need overestimators
        f_lb = b.P_rat.lb ** b.gamma
        f_ub = b.P_rat.ub ** b.gamma
        df_lb = b.gamma * b.P_rat.lb ** (b.gamma - 1)
        df_ub = b.gamma * b.P_rat.ub ** (b.gamma - 1)
        lin_cuts.env_Pratio.lb_overest = Constraint(
            expr=b.P_ratio <= df_lb * (b.P_rat - b.P_rat.lb) + f_lb)
        lin_cuts.env_Pratio.ub_overest = Constraint(
            expr=b.P_ratio <= df_ub * (b.P_rat - b.P_rat.ub) + f_ub)

        lin_cuts.mc_exit_T = Block()
        add_mccormick_relaxation(lin_cuts.mc_exit_T, b.Tout, b.Tin, b.P_ratio, nsegs, None, 1.0)

    def get_vars_to_bound(self):
        yield self.block.Pout
        yield self.block.Pin
        yield self.block.P_rat
        yield self.block.P_ratio
