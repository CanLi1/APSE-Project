from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import Var, Constraint, Param, Set, NonNegativeReals, Block, Reals
from idaes_models.core.util.nonlinear_fcns import log, exp
from idaes_models.core.util.var import none_if_empty as ne
from idaes_models.core.ports import InPort, OutPort


@ProcBlock("Flash")
class _Flash(UnitModel):
    """ Flash unit
    """
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('type', [])
        kwargs['type'].append('flash')
        self._antoine = kwargs.pop('antoine')
        self.key_comp = kwargs.pop('key_comp')
        super(_Flash, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        super(_Flash, self)._build_unit_sets()
        b = self
        b.antoine_coefficients = Set(initialize=['A', 'B', 'C'])

    def _build_unit_params(self):
        super(_Flash, self)._build_unit_params()
        b = self

        def ant_init(b, a, c): return self._antoine[a][c]
        b.antoine = Param(b.antoine_coefficients, self.comps, initialize=ant_init)

    def _build_unit_vars(self):
        super(_Flash, self)._build_unit_vars()
        b = self
        b.T = Var(domain=NonNegativeReals, initialize=400, bounds=(300, 500))
        b.P = Var(domain=NonNegativeReals, initialize=25, bounds=(2.5, 150))
        b.flow_in = Var(self.comps, domain=NonNegativeReals, initialize=4, bounds=(0, 20))
        b.flow_vap = Var(self.comps, domain=NonNegativeReals, initialize=3.5, bounds=(0, 20))
        b.flow_liq = Var(self.comps, domain=NonNegativeReals, initialize=0.2, bounds=(0, 1))
        b.total_flow_liq = Var(domain=NonNegativeReals, bounds=(0.1, 1.0))
        b.Pvap = Var(self.comps, domain=NonNegativeReals)
        b.vap_recover = Var(self.comps, domain=NonNegativeReals, initialize=0.5, bounds=(0.01, 0.9999))
        b.lnPvap = Var(self.comps, domain=Reals)
        b.Tinv = Var(self.comps, domain=Reals)
        b.nonkey_comps = self.comps - self.key_comp
        b.recov_keyPvap = Var(self.comps - self.key_comp, domain=NonNegativeReals)
        b.keyRecov_Pvap = Var(self.comps - self.key_comp, domain=NonNegativeReals)
        b.recovs_keyPvap = Var(self.comps - self.key_comp, domain=NonNegativeReals)
        b.recovs_Pvap = Var(self.comps - self.key_comp, domain=NonNegativeReals)
        b.PF = Var(domain=NonNegativeReals)
        b.PvapF = Var(self.comps, domain=NonNegativeReals)

        # Initialize variables that require calculations
        for c in self.comps:
            b.Pvap[c].setlb((10.0 / 7500.6168) * exp(b.antoine['A', c] - b.antoine['B', c] / (b.T.lb - b.antoine['C', c])))
            b.Pvap[c].value = (10.0 / 7500.6168) * exp(b.antoine['A', c] - b.antoine['B', c] / (b.T.value - b.antoine['C', c]))
            b.Pvap[c].setub((10.0 / 7500.6168) * exp(b.antoine['A', c] - b.antoine['B', c] / (b.T.ub - b.antoine['C', c])))
        b.vap_recover['H2'].value = 0.995
        for c in self.comps:
            b.vap_recover[c].value = (b.Pvap[c].value * b.vap_recover['H2'].value / b.Pvap['H2'].value) / (1.0 + b.vap_recover['H2'].value * (b.Pvap[c].value - b.Pvap['H2'].value) / b.Pvap['H2'].value)
            b.flow_vap[c].value = b.vap_recover[c].value * b.flow_in[c].value
            b.flow_liq[c].value = b.flow_in[c].value - b.flow_vap[c].value
        b.total_flow_liq.value = sum(b.flow_liq[c].value for c in self.comps)

    def _build_unit_constraints(self):
        super(_Flash, self)._build_unit_constraints()
        b = self

        equip = b.equip

        def mass_balance(equip, c, *vs):
            b = equip.parent_block()
            return b.flow_in[c, vs] == b.flow_vap[c, vs] + b.flow_liq[c, vs]
        equip.mass_balance = Constraint(self.comps, rule=mass_balance)

        def total_liq_flow_calc(equip, *vs):
            b = equip.parent_block()
            return b.total_flow_liq[ne(vs)] == sum(b.flow_liq[c, vs] for c in self.comps)
        equip.total_liq_flow_calc = Constraint(rule=total_liq_flow_calc)

        # -------------------------
        # - Nonlinear constraints -
        # -------------------------

        def Pvap_calc(equip, c, *vs):
            b = equip.parent_block()
            return log(b.Pvap[c, vs] / 10 * 7500.6168) == b.antoine['A', c] - b.antoine['B', c] / (b.T[ne(vs)] - b.antoine['C', c])
        equip.Pvap_calc = Constraint(self.comps, rule=Pvap_calc)

        def recovery_calc(equip, c, *vs):
            b = equip.parent_block()
            return b.vap_recover[self.key_comp, vs] * (b.vap_recover[c, vs] * b.Pvap[self.key_comp, vs] + (1 - b.vap_recover[c, vs]) * b.Pvap[c, vs]) == b.Pvap[self.key_comp, vs] * b.vap_recover[c, vs]
        equip.recovery_calc = Constraint(b.nonkey_comps, rule=recovery_calc)

        def component_recovery(equip, c, *vs):
            b = equip.parent_block()
            return b.flow_vap[c, vs] == b.vap_recover[c, vs] * b.flow_in[c, vs]
        equip.component_recovery = Constraint(self.comps, rule=component_recovery)

        def flash_pressure(equip, *vs):
            b = equip.parent_block()
            return b.P[ne(vs)] * b.total_flow_liq[ne(vs)] == sum(b.Pvap[c, vs] * b.flow_liq[c, vs] for c in self.comps)
        equip.flash_pressure = Constraint(rule=flash_pressure)

    def _build_unit_ports(self):
        super(_Flash, self)._build_unit_ports()
        b = self
        self.add_unit(InPort(name='inlet', parent=self, fc=b.flow_in, T=b.T, P=b.P))
        self.add_unit(OutPort(name='vap_out', parent=self, fc=b.flow_vap, T=b.T, P=b.P))
        self.add_unit(OutPort(name='liq_out', parent=self, fc=b.flow_liq, T=b.T, P=b.P))

    def get_flow_vars(self):
        yield self.flow_in
        yield self.flow_vap
        yield self.flow_liq

    def apply_linear_relaxations(self, nsegs=1):
        b = self
        lin_cuts = b.lin_cuts
        from idaes_models.core.util.concave import add_concave_linear_underest
        from idaes_models.core.util.mccormick import add_mccormick_relaxation
        from idaes_models.core.util.var import max_ub, min_lb
        super(_Flash, self).apply_linear_relaxations(nsegs=nsegs)

        def concave_lnPvap(c, Pvap):
            return log(Pvap / 10 * 7500.6168)
        lin_cuts.env_lnPvap = Block()
        add_concave_linear_underest(lin_cuts.env_lnPvap, 'concave_lnPvap', 1, b.Pvap, b.lnPvap, concave_lnPvap, self.comps)
        # also need overestimators
        f_lb = {}
        f_ub = {}
        df_lb = {}
        df_ub = {}
        for c in self.comps:
            f_lb[c] = log(b.Pvap[c].lb / 10 * 7500.6168)
            f_ub[c] = log(b.Pvap[c].ub / 10 * 7500.6168)
            df_lb[c] = 1.0 / b.Pvap[c].lb
            df_ub[c] = 1.0 / b.Pvap[c].ub

        def lnPvap_lb_overest(lin_cuts, c):
            return b.lnPvap[c] <= df_lb[c] * (b.Pvap[c] - b.Pvap[c].lb) + f_lb[c]
        lin_cuts.env_lnPvap.lnPvap_lb_overest = Constraint(self.comps, rule=lnPvap_lb_overest)

        def lnPvap_ub_overest(lin_cuts, c):
            return b.lnPvap[c] <= df_ub[c] * (b.Pvap[c] - b.Pvap[c].ub) + f_ub[c]
        lin_cuts.env_lnPvap.lnPvap_ub_overest = Constraint(self.comps, rule=lnPvap_ub_overest)

        lin_cuts.mc_Tinv = Block()
        for c in self.comps:
            b.Tinv[c].setub(max(b.antoine['B', c] / (b.T.ub - b.antoine['C', c]), b.antoine['B', c] / (b.T.lb - b.antoine['C', c])))
            b.Tinv[c].setlb(min(b.antoine['B', c] / (b.T.ub - b.antoine['C', c]), b.antoine['B', c] / (b.T.lb - b.antoine['C', c])))
        for c in self.comps:
            add_mccormick_relaxation(lin_cuts.mc_Tinv, b.antoine['B', c], b.T - b.antoine['C', c], b.Tinv[c], nsegs, c, 1.0)

        def linear_Pvap(lin_cuts, c):
            return b.lnPvap[c] == b.antoine['A', c] - b.Tinv[c]
        lin_cuts.linear_Pvap = Constraint(self.comps, rule=linear_Pvap)

        # recovery_calc
        for c in b.nonkey_comps:
            b.recov_keyPvap[c].setub(max_ub(b.vap_recover[c], b.Pvap[self.key_comp]))
            b.recov_keyPvap[c].setlb(min_lb(b.vap_recover[c], b.Pvap[self.key_comp]))
            b.keyRecov_Pvap[c].setub(max_ub(b.vap_recover[self.key_comp], b.Pvap[c]))
            b.keyRecov_Pvap[c].setlb(min_lb(b.vap_recover[self.key_comp], b.Pvap[c]))
            b.recovs_keyPvap[c].setub(max_ub(b.vap_recover[self.key_comp], b.recov_keyPvap[c]))
            b.recovs_keyPvap[c].setlb(min_lb(b.vap_recover[self.key_comp], b.recov_keyPvap[c]))
            b.recovs_Pvap[c].setub(max_ub(b.vap_recover[c], b.keyRecov_Pvap[c]))
            b.recovs_Pvap[c].setlb(min_lb(b.vap_recover[c], b.keyRecov_Pvap[c]))

        lin_cuts.mc_recov_keyPvap = Block()
        for c in b.nonkey_comps:
            add_mccormick_relaxation(lin_cuts.mc_recov_keyPvap, b.recov_keyPvap[c], b.vap_recover[c], b.Pvap[self.key_comp], nsegs, c, 1.0)

        lin_cuts.mc_keyRecov_Pvap = Block()
        for c in b.nonkey_comps:
            add_mccormick_relaxation(lin_cuts.mc_keyRecov_Pvap, b.keyRecov_Pvap[c], b.vap_recover[self.key_comp], b.Pvap[c], nsegs, c, 1.0)

        lin_cuts.mc_recovs_keyPvap = Block()
        for c in b.nonkey_comps:
            add_mccormick_relaxation(lin_cuts.mc_recovs_keyPvap, b.recovs_keyPvap[c], b.vap_recover[self.key_comp], b.recov_keyPvap[c], nsegs, c, 1.0)

        lin_cuts.mc_recovs_Pvap = Block()
        for c in b.nonkey_comps:
            add_mccormick_relaxation(lin_cuts.mc_recovs_Pvap, b.recovs_Pvap[c], b.vap_recover[c], b.keyRecov_Pvap[c], nsegs, c, 1.0)

        def linear_recovery_calc(lin_cuts, c):
            return b.recovs_keyPvap[c] + b.keyRecov_Pvap[c] - b.recovs_Pvap[c] == b.recov_keyPvap[c]
        lin_cuts.linear_recovery_calc = Constraint(b.nonkey_comps, rule=linear_recovery_calc)

        lin_cuts.mc_component_recovery = Block()
        for c in self.comps:
            add_mccormick_relaxation(lin_cuts.mc_component_recovery, b.flow_vap[c], b.flow_in[c], b.vap_recover[c], nsegs, c, 1.0)

        # flash pressure
        lin_cuts.mc_PF = Block()
        b.PF.setub(max_ub(b.P, b.total_flow_liq))
        b.PF.setlb(min_lb(b.P, b.total_flow_liq))
        add_mccormick_relaxation(lin_cuts.mc_PF, b.PF, b.total_flow_liq, b.P, nsegs, None, 1.0)

        lin_cuts.mc_PvapF = Block()
        for c in self.comps:
            b.PvapF[c].setub(max_ub(b.Pvap[c], b.flow_liq[c]))
            b.PvapF[c].setlb(min_lb(b.Pvap[c], b.flow_liq[c]))
        for c in self.comps:
            add_mccormick_relaxation(lin_cuts.mc_PvapF, b.PvapF[c], b.flow_liq[c], b.Pvap[c], nsegs, c, 1.0)

        lin_cuts.linear_flash_pressure = Constraint(expr=b.PF == sum(b.PvapF[c] for c in self.comps))

    def get_vars_to_bound(self):
        yield self.block.flow_in
        yield self.block.flow_liq
        yield self.block.Pvap
        yield self.block.T
        yield self.block.vap_recover
        yield self.block.recov_keyPvap
        yield self.block.keyRecov_Pvap
