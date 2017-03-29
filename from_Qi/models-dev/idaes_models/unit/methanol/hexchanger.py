from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import Var, Param, Constraint, NonNegativeReals, Binary, Block, Reals
from idaes_models.core.util.var import none_if_empty as ne
from idaes_models.core.ports import InPort, OutPort

class BaseHeatExchanger(UnitModel):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('type', [])
        kwargs['type'].append('HEX')
        self._Cp = kwargs.pop('Cp', 35.0)
        self._cost_Q = kwargs.pop('cost_Q', 0)
        self._P_init = kwargs.pop('P_init', 1.0)
        super(BaseHeatExchanger, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        super(BaseHeatExchanger, self)._build_unit_sets()

    def _build_unit_params(self):
        super(BaseHeatExchanger, self)._build_unit_params()
        b = self
        # subclasses should declare Q_sign (-1: cooler, 1: heater)
        b.cost_Q = Param(initialize=self._cost_Q)

    def _build_unit_vars(self):
        super(BaseHeatExchanger, self)._build_unit_vars()
        b = self
        b.Tin = Var(domain=NonNegativeReals, initialize=305, bounds=(300, 900))
        b.Tout = Var(domain=NonNegativeReals, initialize=301, bounds=(300, 900))
        b.P = Var(domain=NonNegativeReals, initialize=self._P_init, bounds=(0, 150))
        b.flow = Var(self.comps, domain=NonNegativeReals, initialize=1, bounds=(0, 20))
        b.total_flow = Var(domain=NonNegativeReals, initialize=1, bounds=(0, 20))
        b.equip_cost = Var(within=NonNegativeReals, initialize=0, bounds=(0, 20000))
        b.equip_exists = Var(domain=Binary, initialize=1)
        b.Q = Var(within=NonNegativeReals, initialize=1, bounds=(0, 50))
        b.Cp = Param(initialize=self._Cp)
        # dT = Tout - Tin
        b.dT = Var(domain=Reals, initialize=0, bounds=(-1800, 1800))
        # FdT = total_flow * dT
        b.FdT = Var(domain=Reals, initialize=0, bounds=(-36000, 36000))

    def _build_unit_constraints(self):
        super(BaseHeatExchanger, self)._build_unit_constraints()
        b = self

        equip = b.equip

        def total_flow_calc(equip, *vs):
            b = equip.parent_block()
            return b.total_flow[ne(vs)] == sum(b.flow[c, vs] for c in self.comps)
        equip.total_flow_calc = Constraint(rule=total_flow_calc)

        def cost_calc(equip, *vs):
            b = equip.parent_block()
            return b.equip_cost[ne(vs)] >= b.Q[ne(vs)] * b.cost_Q - (b.Q[ne(vs)].ub * b.cost_Q * (1 - b.equip_exists))
        equip.cost_calc = Constraint(rule=cost_calc)

        def flow_ub(equip, c, *vs):
            b = equip.parent_block()
            return b.flow[c, vs] <= b.flow[c, vs].ub * b.equip_exists
        equip.flow_ub = Constraint(self.comps, rule=flow_ub)

        def cost_ub(equip, *vs):
            b = equip.parent_block()
            return b.equip_cost[ne(vs)] <= b.equip_cost[ne(vs)].ub * b.equip_exists
        equip.cost_ub = Constraint(rule=cost_ub)

        def Q_ub(equip, *vs):
            b = equip.parent_block()
            return b.Q[ne(vs)] <= b.Q[ne(vs)].ub * b.equip_exists
        equip.no_Q = Constraint(rule=Q_ub)

        def dT_defn(equip, *vs):
            b = equip.parent_block()
            return b.dT[ne(vs)] == b.Tout[ne(vs)] - b.Tin[ne(vs)]
        equip.dT_defn = Constraint(rule=dT_defn)

        def duty(equip, *vs):
            b = equip.parent_block()
            return b.Q_sign * b.Q[ne(vs)] == b.FdT[ne(vs)] * b.Cp * (3600 * 8500 * 1.0E-12)
        equip.duty = Constraint(rule=duty)

        def FdT(nl_con, *vs):
            b = equip.parent_block()
            return b.FdT[ne(vs)] == b.total_flow[ne(vs)] * b.dT[ne(vs)]
        equip.FdT = Constraint(rule=FdT)

    def _build_unit_ports(self):
        super(BaseHeatExchanger, self)._build_unit_ports()
        b = self

        self.add_unit(InPort(name='inlet', parent=self, fc=b.flow, T=b.Tin, P=b.P))
        self.add_unit(OutPort(name='outlet', parent=self, fc=b.flow, T=b.Tout, P=b.P))

    def get_flow_vars(self):
        yield self.flow

    def apply_linear_relaxations(self, nsegs=1):
        b = self
        lin_cuts = b.lin_cuts
        from idaes_models.core.util.mccormick import add_mccormick_relaxation
        from idaes_models.core.util.var import max_ub, min_lb
        super(BaseHeatExchanger, self).apply_linear_relaxations(nsegs=nsegs)

        lin_cuts.mc_FdT = Block()
        b.FdT.setub(max_ub(b.dT, b.total_flow))
        b.FdT.setlb(min_lb(b.dT, b.total_flow))
        add_mccormick_relaxation(lin_cuts.mc_FdT, b.FdT, b.total_flow, b.dT, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

    def get_vars_to_bound(self):
        yield self.total_flow, self.block_bounds
        yield self.Tin, self.block_bounds
        yield self.Tout, self.block_bounds
        yield self.dT, self.block_bounds


@ProcBlock("Cooler")
class _Cooler(BaseHeatExchanger):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('type', [])
        kwargs['type'].append('Cooler')
        kwargs.setdefault('cost_Q', 700)  # Cooling cost $700 per 1E9 kJ
        super(_Cooler, self).__init__(*args, **kwargs)

    def _build_unit_params(self):
        super(_Cooler, self)._build_unit_params()
        b = self
        b.Q_sign = Param(initialize=-1)

    def _build_unit_constraints(self):
        super(_Cooler, self)._build_unit_constraints()
        b = self
        b.equip.cooler_defn = Constraint(expr=b.Tin >= b.Tout)


@ProcBlock("Heater")
class _Heater(BaseHeatExchanger):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('type', [])
        kwargs['type'].append('Heater')
        kwargs.setdefault('cost_Q', 8000)  # Heating cost $8000 per 1E9 kJ
        super(_Heater, self).__init__(*args, **kwargs)

    def _build_unit_params(self):
        super(_Heater, self)._build_unit_params()
        b = self
        b.Q_sign = Param(initialize=1)

    def _build_unit_constraints(self):
        super(_Heater, self)._build_unit_constraints()
        b = self
        b.equip.heater_defn = Constraint(expr=b.Tin <= b.Tout)
