from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import Param, Var, Constraint, NonNegativeReals, Binary, Block, NonPositiveReals
from idaes_models.core.util.nonlinear_fcns import exp
from idaes_models.core.util.var import none_if_empty as ne
from idaes_models.core.ports import InPort, OutPort
__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


@ProcBlock("Reactor")
class _Reactor(UnitModel):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('type', [])
        kwargs['type'].append('Reactor')
        self.key_comp = kwargs.pop('key_comp')
        self._stoich = kwargs.pop('stoich')
        self._v_param = kwargs.pop('v_param')
        self._vol_cost = kwargs.pop('vol_cost')
        self._rxn_cost = kwargs.pop('rxn_cost')
        self._Cp = kwargs.pop('Cp', 35.0)
        self._heat_rxn = kwargs.pop('heat_rxn', -15.0)
        super(_Reactor, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        super(_Reactor, self)._build_unit_sets()

    def _build_unit_params(self):
        super(_Reactor, self)._build_unit_params()
        b = self
        b.stoich = Param(self.comps, initialize=self._stoich)
        b.Cp = Param(initialize=self._Cp)
        b.heat_rxn = Param(initialize=self._heat_rxn)
        b.v_param = Param(initialize=self._v_param)
        b.vol_cost = Param(initialize=self._vol_cost)
        b.rxn_fixed_cost = Param(initialize=self._rxn_cost)
        b.activ_E = Param(initialize=-1800)

    def _build_unit_vars(self):
        super(_Reactor, self)._build_unit_vars()
        b = self
        b.Tin = Var(domain=NonNegativeReals, initialize=500, bounds=(423, 873))
        b.Tout = Var(domain=NonNegativeReals, initialize=535, bounds=(523, 873))
        b.Pin = Var(domain=NonNegativeReals, initialize=25, bounds=(25, 150))
        b.Pout = Var(domain=NonNegativeReals, initialize=25, bounds=(25, 150))
        b.flow_in = Var(self.comps, domain=NonNegativeReals, initialize=1.0, bounds=(0, 20))
        b.flow_out = Var(self.comps, domain=NonNegativeReals, initialize=1.0, bounds=(0, 20))
        b.equip_cost = Var(domain=NonNegativeReals, initialize=0, bounds=(0, 15000))
        b.equip_exists = Var(domain=Binary, initialize=1)
        b.V = Var(domain=NonNegativeReals, initialize=100, bounds=(0, 100))
        b.total_flow_in = Var(domain=NonNegativeReals, bounds=(0, 20))
        b.total_flow_out = Var(domain=NonNegativeReals, bounds=(0, 20))
        # equilibrium conversion
        b.eq_conv = Var(domain=NonNegativeReals, bounds=(0, 0.42))
        b.conv = Var(domain=NonNegativeReals, bounds=(0, 0.42))
        b.extent = Var(domain=NonNegativeReals, bounds=(0, 5))
        b.FTin = Var(domain=NonNegativeReals)
        b.FTout = Var(domain=NonNegativeReals)
        b.InvTout = Var(domain=NonPositiveReals)
        b.Psq = Var(domain=NonNegativeReals)
        b.RxnExp = Var(domain=NonNegativeReals)
        b.InhibFactor = Var(domain=NonNegativeReals)

        self.init_calculated_values()

    def _build_unit_constraints(self):
        super(_Reactor, self)._build_unit_constraints()
        b = self

        equip = b.equip
        equip.equip = Block()

        # if block is active, fix volume to 100 m^3
        def react_vol(equip, *vs):
            b = equip.parent_block()
            return b.V[ne(vs)] >= 100 * b.equip_exists
        equip.react_vol = Constraint(rule=react_vol)

        def mass_balance(equip, c, *vs):
            b = equip.parent_block()
            return b.flow_out[c, vs] == b.flow_in[c, vs] + b.stoich[c] * b.extent[ne(vs)]
        equip.mass_balance = Constraint(self.comps, rule=mass_balance)

        def flow_out_lb(equip, c, *vs):
            b = equip.parent_block()
            return 0.01 * b.equip_exists <= b.flow_out[c, vs]
        equip.flow_out_lb = Constraint(self.comps, rule=flow_out_lb)

        def total_flow_in_calc(equip, *vs):
            b = equip.parent_block()
            return b.total_flow_in[ne(vs)] == sum(b.flow_in[c, vs] for c in self.comps)
        equip.total_flow_in_calc = Constraint(rule=total_flow_in_calc)

        def total_flow_out_calc(equip, *vs):
            b = equip.parent_block()
            return b.total_flow_out[ne(vs)] == sum(b.flow_out[c, vs] for c in self.comps)
        equip.total_flow_out_calc = Constraint(rule=total_flow_out_calc)

        def pressure_drop(equip, *vs):
            b = equip.parent_block()
            return b.Pout[ne(vs)] == 0.9 * b.Pin[ne(vs)]
        equip.pressure_drop = Constraint(rule=pressure_drop)

        def cost_calc(equip, *vs):
            b = equip.parent_block()
            return b.equip_cost[ne(vs)] >= b.vol_cost * b.V[ne(vs)] - (b.vol_cost * b.V[ne(vs)].ub * (1 - b.equip_exists)) + b.rxn_fixed_cost * b.equip_exists
        equip.cost_calc = Constraint(rule=cost_calc)

        def flow_in_ub(equip, c, *vs):
            b = equip.parent_block()
            return b.flow_in[c, vs] <= b.flow_in[c, vs].ub * b.equip_exists
        equip.flow_in_ub = Constraint(self.comps, rule=flow_in_ub)

        def flow_out_ub(equip, c, *vs):
            b = equip.parent_block()
            return b.flow_out[c, vs] <= b.flow_out[c, vs].ub * b.equip_exists
        equip.flow_out_ub = Constraint(self.comps, rule=flow_out_ub)

        def cost_ub(equip, *vs):
            b = equip.parent_block()
            return b.equip_cost[ne(vs)] <= b.equip_cost[ne(vs)].ub * b.equip_exists
        equip.cost_ub = Constraint(rule=cost_ub)

        def rxn_conv_relation(equip, *vs):
            # Conversion <= equilibrium conversion
            b = equip.parent_block()
            return b.conv[ne(vs)] <= b.eq_conv[ne(vs)]
        equip.rxn_conv_relation = Constraint(rule=rxn_conv_relation)

        self.__add_nonlinear_constraints()

    def _build_unit_ports(self):
        super(_Reactor, self)._build_unit_ports()
        b = self
        self.add_unit(InPort(name='inlet', parent=self, fc=b.flow_in, T=b.Tin, P=b.Pin))
        self.add_unit(OutPort(name='outlet', parent=self, fc=b.flow_out, T=b.Tout, P=b.Pout))

    def init_calculated_values(self):
        b = self
        # Intialize calculated model variable values
        b.total_flow_in.value = sum(b.flow_in[c].value for c in self.comps)
        b.eq_conv.value = 0.415 * (1 - (26.25 * exp(-1800 / b.Tout.value)) / (b.Pin.value / 10) ** 2.0)
        b.conv.value = b.eq_conv.value * (1 - exp(-b.v_param * b.V.value)) * (b.flow_in['H2'].value + b.flow_in['CO'].value + b.flow_in['MeOH'].value) / b.total_flow_in.value
        b.extent.value = b.conv.value * b.flow_in[self.key_comp].value
        for c in self.comps:
            b.flow_out[c].value = b.flow_in[c].value + b.stoich[c] * b.extent.value
        b.total_flow_out.value = sum(b.flow_out[c].value for c in self.comps)

    def __add_nonlinear_constraints(self):
        b = self
        equip = b.equip

        def rxn_extent(equip, *vs):
            b = equip.parent_block()
            return b.extent[ne(vs)] == b.conv[ne(vs)] * b.flow_in[self.key_comp, vs]
        equip.rxn_extent = Constraint(rule=rxn_extent)

        def heat_balance(equip, *vs):
            b = equip.parent_block()
            return (b.total_flow_in[ne(vs)] * b.Tin[ne(vs)] - b.total_flow_out[ne(vs)] * b.Tout[ne(vs)]) * b.Cp == b.heat_rxn * b.extent[ne(vs)]
        equip.heat_balance = Constraint(rule=heat_balance)

        def equil_conversion(equip, *vs):
            b = equip.parent_block()
            return b.eq_conv[ne(vs)] == 0.415 * (1 - (26.25 * exp(b.activ_E / b.Tout[ne(vs)])) / (b.Pin[ne(vs)] / 10.0) ** 2.0)
        equip.equil_conversion = Constraint(rule=equil_conversion)

        def conversion(equip, *vs):
            b = equip.parent_block()
            return b.conv[ne(vs)] == b.eq_conv[ne(vs)] * (1 - exp(-b.v_param * b.V[ne(vs)])) * (b.flow_in['H2', vs] + b.flow_in['CO', vs] + b.flow_in['MeOH', vs]) / b.total_flow_in[ne(vs)]
        equip.conversion = Constraint(rule=conversion)

    def get_flow_vars(self):
        yield self.flow_in
        yield self.flow_out

    def apply_linear_relaxations(self, nsegs=1):
        b = self
        lin_cuts = b.lin_cuts
        from idaes_models.core.util.mccormick import add_mccormick_relaxation
        from idaes_models.core.util.var import min_lb, max_ub, ub, lb
        super(_Reactor, self).apply_linear_relaxations(nsegs=nsegs)

        lin_cuts.mc_rxn_extent = Block()
        # extent = flow_in_of_key * conversion
        add_mccormick_relaxation(lin_cuts.mc_rxn_extent, b.extent, b.flow_in[self.key_comp], b.conv, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_FTin = Block()
        # FTin = total_flow_in * Tin
        b.FTin.setub(max_ub(b.total_flow_in, b.Tin))
        b.FTin.setlb(min_lb(b.total_flow_in, b.Tin))
        add_mccormick_relaxation(lin_cuts.mc_FTin, b.FTin, b.total_flow_in, b.Tin, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_FTout = Block()
        # FTout = total_flow_out * Tout
        b.FTout.setub(max_ub(b.total_flow_out, b.Tout))
        b.FTout.setlb(min_lb(b.total_flow_out, b.Tout))
        add_mccormick_relaxation(lin_cuts.mc_FTout, b.FTout, b.total_flow_out, b.Tout, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.linear_heat_balance = Constraint(expr=(b.FTin - b.FTout) * b.Cp == b.heat_rxn * b.extent)

        lin_cuts.mc_InvTout = Block()
        # InvTout = -1800 / Tout
        b.InvTout.setlb(min(b.activ_E / lb(b.Tout), b.activ_E / ub(b.Tout)))
        b.InvTout.setub(max(b.activ_E / lb(b.Tout), b.activ_E / ub(b.Tout)))
        add_mccormick_relaxation(lin_cuts.mc_InvTout, b.activ_E, b.Tout, b.InvTout, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        # exponential function is actually convex, need to linearize
        # TODO add support for equip_exists
        # TODO add support for block bounds
        lin_cuts.env_RxnExp = Block()
        # RxnExp = 26.25 * exp(InvTout)
        f_lb = 26.25 * exp(b.InvTout.lb)
        f_ub = 26.25 * exp(b.InvTout.ub)
        df_lb = f_lb
        df_ub = f_ub
        b.RxnExp.setub(f_ub)
        b.RxnExp.setlb(f_lb)
        lin_cuts.env_RxnExp.lin_overest = Constraint(
            expr=b.RxnExp <= (f_ub - f_lb) / (b.InvTout.ub - b.InvTout.lb) *
            (b.InvTout - b.InvTout.lb) + f_lb)
        lin_cuts.env_RxnExp.lb_underest = Constraint(
            expr=b.RxnExp >= df_lb * (b.InvTout - b.InvTout.lb) + f_lb)
        lin_cuts.env_RxnExp.ub_underest = Constraint(
            expr=b.RxnExp >= df_ub * (b.InvTout - b.InvTout.ub) + f_ub)

        # P squared is convex, need to linearize
        # TODO add support for equip_exists
        # TODO add support for block bounds
        lin_cuts.env_Psq = Block()
        # Psq = (P / 10) ^ 2
        b.Psq.setub((b.Pin.ub / 10) ** 2)
        b.Psq.setlb((b.Pin.lb / 10) ** 2)
        f_lb = (b.Pin.lb / 10) ** 2
        f_ub = (b.Pin.ub / 10) ** 2
        df_lb = (2 / 10) * (b.Pin.lb / 10)
        df_ub = (2 / 10) * (b.Pin.ub / 10)
        lin_cuts.env_Psq.lin_overest = Constraint(
            expr=b.Psq <= (f_ub - f_lb) / (b.Pin.ub - b.Pin.lb) * (b.Pin - b.Pin.lb) + f_lb)
        lin_cuts.env_Psq.lb_underest = Constraint(
            expr=b.Psq >= df_lb * (b.Pin - b.Pin.lb) + f_lb)
        lin_cuts.env_Psq.ub_underest = Constraint(
            expr=b.Psq >= df_ub * (b.Pin - b.Pin.ub) + f_ub)

        lin_cuts.mc_InhibFactor = Block()
        # InhibFactor = RxnExp / Psq
        # set bounds on InhibFactor
        b.InhibFactor.setub(max(lb(b.RxnExp) / lb(b.Psq), lb(b.RxnExp) / ub(b.Psq), ub(b.RxnExp) / lb(b.Psq), ub(b.RxnExp) / ub(b.Psq)))
        b.InhibFactor.setlb(min(lb(b.RxnExp) / lb(b.Psq), lb(b.RxnExp) / ub(b.Psq), ub(b.RxnExp) / lb(b.Psq), ub(b.RxnExp) / ub(b.Psq)))
        add_mccormick_relaxation(lin_cuts.mc_InhibFactor, b.RxnExp, b.InhibFactor, b.Psq, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

    def get_vars_to_bound(self):
        yield self.block.flow_in, self.block_bounds
        yield self.block.conv, self.block_bounds
        yield self.block.total_flow_in, self.block_bounds
        yield self.block.Tin, self.block_bounds
        yield self.block.total_flow_out, self.block_bounds
        yield self.block.Tout, self.block_bounds
        yield self.block.InvTout, self.block_bounds
        yield self.block.Pin, self.block_bounds
        yield self.block.InhibFactor, self.block_bounds
