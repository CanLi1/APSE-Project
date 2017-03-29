from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import Var, Constraint, Block, NonNegativeReals, Param
from idaes_models.core.util.var import none_if_empty as ne
from idaes_models.core.ports import InPort

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


@ProcBlock("Sink")
class _Sink(UnitModel):
    """ Simple sink node
    """

    def __init__(self, *args, **kwargs):
        """ Initialize sink """
        kwargs.setdefault('type', [])
        kwargs['type'].append('Sink')
        self._price = kwargs.pop('price', None)
        self.fix_T = kwargs.pop('T', None)
        self.fix_P = kwargs.pop('P', None)
        self.min_frac = kwargs.pop('min_frac', {})
        super(_Sink, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        super(_Sink, self)._build_unit_sets()

    def _build_unit_params(self):
        super(_Sink, self)._build_unit_params()

    def _build_unit_vars(self):
        super(_Sink, self)._build_unit_vars()
        b = self

        # Temperature
        T = b.T = Var(domain=NonNegativeReals, initialize=300, bounds=(300, 900))
        if self.fix_T is not None:
            T.fix(self.fix_T)
        # Pressure
        P = b.P = Var(domain=NonNegativeReals, initialize=1, bounds=(0, 150))
        if self.fix_P is not None:
            P.fix(self.fix_P)
        # Flow
        b.flow = Var(self.comps, domain=NonNegativeReals, initialize=1, bounds=(0, 20))
        b.total_flow_in = Var(domain=NonNegativeReals, initialize=1, bounds=(0, 20))
        b.price = Param(initialize=self._price)
        b.equip_profit = Var(domain=NonNegativeReals)
        if self.min_frac:
            frac = b.frac =\
                Var(self.comps, domain=NonNegativeReals, bounds=(0, 1))
            for c in self.comps:
                if c in self.min_frac:
                    frac[c].setlb(self.min_frac[c])

    def _build_unit_constraints(self):
        super(_Sink, self)._build_unit_constraints()
        b = self

        equip = b.equip

        def total_flow_in(equip, *vs):
            b = equip.parent_block()
            return b.total_flow_in[ne(vs)] == sum(b.flow[c, vs] for c in self.comps)
        equip.total_flow_in = Constraint(rule=total_flow_in)

        if self.min_frac:
            def frac_flow(equip, c, *vs):
                b = equip.parent_block()
                return b.flow[c, vs] == b.total_flow_in[ne(vs)] * b.frac[c, vs]
            equip.frac_flow = Constraint(self.comps, rule=frac_flow)

        def profit_calc(equip, *vs):
            b = equip.parent_block()
            return b.equip_profit[ne(vs)] == b.total_flow_in[ne(vs)] * b.price
        equip.profit_calc = Constraint(rule=profit_calc)

    def _build_unit_ports(self):
        super(_Sink, self)._build_unit_ports()
        b = self
        self.add_unit(InPort(name='inlet', parent=self, fc=b.flow, T=b.T, P=b.P))

    def get_flow_vars(self):
        yield self.flow

    def generate_cut_gen_problem(self):
        if self.min_frac:
            return super(_Sink, self).generate_cut_gen_problem()

    def apply_linear_relaxations(self, nsegs=1):
        super(_Sink, self).apply_linear_relaxations(nsegs=nsegs)
        if self.min_frac:
            b = self
            lin_cuts = b.lin_cuts
            from idaes_models.core.util.mccormick import add_mccormick_relaxation
            from idaes_models.core.util.var import min_lb, max_ub

            lin_cuts.mc_frac = Block()
            for c in self.comps:
                b.flow[c].setub(max_ub(b.total_flow_in, b.frac[c]))
                b.flow[c].setlb(min_lb(b.total_flow_in, b.frac[c]))
            for c in self.comps:
                add_mccormick_relaxation(lin_cuts.mc_frac, b.flow[c], b.total_flow_in, b.frac[c], nsegs, c, 1.0)

    def get_vars_to_bound(self):
        yield self.total_flow_in
        if self.min_frac:
            yield self.frac
