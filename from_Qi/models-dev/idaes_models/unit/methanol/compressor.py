from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import Var, Constraint, Block, Param, NonNegativeReals, Binary
from idaes_models.core.util.var import none_if_empty as ne
from idaes_models.core.ports import InPort, OutPort

@ProcBlock("Compressor")
class _Compressor(UnitModel):
    """ Simple compressor
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('type', [])
        kwargs['type'].append('Compressor')
        self._alpha = kwargs.pop('alpha', 0.720)
        self._gamma = kwargs.pop('gamma', 0.23077)
        self._compr_eff = kwargs.pop('compr_eff', 0.750)
        self._cost_elec = kwargs.pop('cost_elec', 0.255)
        self._compr_cost = kwargs.pop('compr_cost', 0)
        super(_Compressor, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        super(_Compressor, self)._build_unit_sets()

    def _build_unit_params(self):
        super(_Compressor, self)._build_unit_params()
        b = self
        # Compressor factor
        alpha = b.alpha = Param(initialize=self._alpha)
        # Cp/Cv
        gamma = b.gamma = Param(initialize=self._gamma)
        # compressor efficiency
        compr_eff = b.compr_eff = Param(initialize=self._compr_eff)
        # work coefficient
        b.work_coeff = Param(initialize=alpha / (1000 * compr_eff * gamma))
        # electricity cost in $1000s per kW.yr
        b.cost_elec = Param(initialize=self._cost_elec)
        b.compr_cost = Param(initialize=self._compr_cost)

    def _build_unit_vars(self):
        super(_Compressor, self)._build_unit_vars()
        b = self
        b.flow = Var(self.comps, domain=NonNegativeReals, initialize=1.0, bounds=(0, 20))
        b.Pin = Var(domain=NonNegativeReals, initialize=1.0, bounds=(1, 150))
        b.Pout = Var(domain=NonNegativeReals, initialize=1.5, bounds=(1, 150))
        b.Tin = Var(domain=NonNegativeReals, initialize=300, bounds=(300, 900))
        b.Tout = Var(domain=NonNegativeReals, initialize=300, bounds=(300, 900))
        b.equip_cost = Var(domain=NonNegativeReals, initialize=0, bounds=(0, 20000))
        # work in kW
        b.work = Var(domain=NonNegativeReals, bounds=(0, 50))
        b.equip_exists = Var(domain=Binary, initialize=1)
        b.total_flow = Var(domain=NonNegativeReals, bounds=(0, 20))
        # Pressure ratio between outlet and inlet, to the gamma power
        b.P_ratio = Var(domain=NonNegativeReals, bounds=(1, 1.74))
        # b.P_ratio_m1 = Var(domain=NonNegativeReals, bounds=(0, 0.74), initialize=0.5)
        # Pressure ratio between outlet and inlet, used only for linear terms
        b.P_rat = Var(domain=NonNegativeReals, bounds=(1, 11.03))
        b.FTin = Var(domain=NonNegativeReals, initialize=0)
        b.PratFTin = Var(domain=NonNegativeReals, initialize=0)

        # Calculated initialization
        b.P_ratio.value = (b.Pout.value / b.Pin.value) ** b.gamma
        b.total_flow.value = sum(b.flow[c].value for c in self.comps)
        b.work.value = b.work_coeff * (b.P_ratio.value - 1.0) * b.Tin.value * b.total_flow.value

    def _build_unit_constraints(self):
        super(_Compressor, self)._build_unit_constraints()
        b = self
        equip = b.equip

        def total_flow_calc(equip, *vs):
            b = equip.parent_block()
            return b.total_flow[ne(vs)] == sum(b.flow[c, vs] for c in self.comps)
        equip.total_flow_calc = Constraint(rule=total_flow_calc)

        def cost_calc(equip, *vs):
            b = equip.parent_block()
            return b.equip_cost[ne(vs)] >= b.work[ne(vs)] * 0.175 * 1000 - (b.work[ne(vs)].ub * 0.175 * 1000 * (1 - b.equip_exists)) + b.work[ne(vs)] * b.cost_elec - (b.work[ne(vs)].ub * b.cost_elec * (1 - b.equip_exists)) + b.compr_cost * b.equip_exists
        equip.cost_calc = Constraint(rule=cost_calc)

        def flow_ub(equip, c, *vs):
            b = equip.parent_block()
            return b.flow[c, vs] <= b.flow[c, vs].ub * b.equip_exists
        equip.flow_ub = Constraint(self.comps, rule=flow_ub)

        def cost_ub(equip, *vs):
            b = equip.parent_block()
            return b.equip_cost[ne(vs)] <= b.equip_cost[ne(vs)].ub * b.equip_exists
        equip.cost_ub = Constraint(rule=cost_ub)

        def work_ub(equip, *vs):
            b = equip.parent_block()
            return b.work[ne(vs)] <= b.work[ne(vs)].ub * b.equip_exists
        equip.work_ub = Constraint(rule=work_ub)

        def pressure_ratio(equip, *vs):
            b = equip.parent_block()
            return b.P_ratio[ne(vs)] == (b.Pout[ne(vs)] / b.Pin[ne(vs)]) ** b.gamma
        equip.pressure_ratio = Constraint(rule=pressure_ratio)

        def work_calc(equip, *vs):
            b = equip.parent_block()
            return b.work[ne(vs)] ==\
                b.work_coeff * (b.P_ratio[ne(vs)] - 1.0) * b.Tin[ne(vs)] * b.total_flow[ne(vs)]
            """ can also be written as:
            work = coeff * P_ratio * Tin * flow - coeff * Tin * flow
            """
        equip.work_calc = Constraint(rule=work_calc)

        def exit_temperature(equip, *vs):
            b = equip.parent_block()
            return b.Tout[ne(vs)] == b.P_ratio[ne(vs)] * b.Tin[ne(vs)]
        equip.exit_temperature = Constraint(rule=exit_temperature)

    def _build_unit_ports(self):
        super(_Compressor, self)._build_unit_ports()
        b = self

        self.add_unit(InPort(name='inlet', parent=self, fc=b.flow, T=b.Tin, P=b.Pin))
        self.add_unit(OutPort(name='outlet', parent=self, fc=b.flow, T=b.Tout, P=b.Pout))

    def apply_linear_relaxations(self, nsegs=1):
        b = self
        lin_cuts = b.lin_cuts
        from idaes_models.core.util.mccormick import add_mccormick_relaxation
        from idaes_models.core.util.var import min_lb, max_ub
        from idaes_models.core.util.concave import add_concave_linear_underest
        super(_Compressor, self).apply_linear_relaxations(nsegs=nsegs)

        def concave_Pratio(P_rat):
            return P_rat ** b.gamma
        lin_cuts.env_Pratio = Block()
        add_concave_linear_underest(lin_cuts.env_Pratio, 'concave_Pratio', 1, b.P_rat, b.P_ratio, concave_Pratio, exists=b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_P_rat = Block()
        add_mccormick_relaxation(lin_cuts.mc_P_rat, b.Pout, b.Pin, b.P_rat, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_TinFlow = Block()
        b.FTin.setub(max_ub(b.Tin, b.total_flow))
        b.FTin.setlb(min_lb(b.Tin, b.total_flow))
        add_mccormick_relaxation(lin_cuts.mc_TinFlow, b.FTin, b.total_flow, b.Tin, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_PratFTin = Block()
        b.PratFTin.setub(max_ub(b.FTin, b.P_ratio))
        b.PratFTin.setlb(min_lb(b.FTin, b.P_ratio))
        add_mccormick_relaxation(lin_cuts.mc_PratFTin, b.PratFTin, b.P_ratio, b.FTin, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.linear_work = Constraint(expr=b.work == b.work_coeff * b.PratFTin - b.work_coeff * b.FTin)

        lin_cuts.mc_exit_T = Block()
        add_mccormick_relaxation(lin_cuts.mc_exit_T, b.Tout, b.P_ratio, b.Tin, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

    def get_vars_to_bound(self):
        yield self.Pin
        yield self.Pout
        yield self.P_rat
        yield self.P_ratio
        yield self.total_flow
        yield self.Tin
        yield self.FTin

    def get_flow_vars(self):
        yield self.flow
