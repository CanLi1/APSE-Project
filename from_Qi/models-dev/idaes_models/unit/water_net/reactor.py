from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import Var, Block, Constraint, NonNegativeReals, Param, Binary, RangeSet, ConcreteModel, Objective, minimize, value
from idaes_models.core.util.concave import add_concave_linear_underest
from idaes_models.core.util.var import SliceVar
from idaes_models.core.util.var import none_if_empty as ne
from idaes_models.core.util.cut_gen\
    import count_vars, get_sum_sq_diff, clone_block
from idaes_models.core.util.oa import add_oa_constraints
from idaes_models.core.ports import InPort, OutPort

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


@ProcBlock("Reactor")
class _Reactor(UnitModel):
    """ Treatment unit """

    def __init__(self, *args, **kwargs):
        self._min_flow = kwargs.pop('min_flow', 0)
        self._theta = kwargs.pop('theta', 0)
        self._beta = kwargs.pop('beta', 0)
        self._gamma = kwargs.pop('gamma', 0)
        self._perf = kwargs.pop('perf', None)
        super(_Reactor, self).__init__(*args, **kwargs)
        if self._perf is None:
            self._perf = {c: 0 for c in self.comps}

    def _build_unit_params(self):
        super(_Reactor, self)._build_unit_params()
        b = self

        b.performance = Param(self.comps, initialize=self._perf)
        b.min_flow = Param(initialize=self._min_flow)
        b.theta = Param(initialize=self._theta)
        b.beta = Param(initialize=self._beta)
        b.gamma = Param(initialize=self._gamma)

    def _build_unit_vars(self):
        super(_Reactor, self)._build_unit_vars()
        b = self

        b.flow_in = Var(
            self.comps, domain=NonNegativeReals, bounds=(0, self._max_flow), initialize=1)
        b.flow_out = Var(
            self.comps, domain=NonNegativeReals, bounds=(0, self._max_flow), initialize=1)
        b.equip_exists = Var(domain=Binary, initialize=1)

        b.lCTU = Var(domain=NonNegativeReals, initialize=0)
        b.nCTU = Var(domain=NonNegativeReals, initialize=0)
        b.CTU = b.lCTU + b.nCTU

    def _build_unit_constraints(self):
        super(_Reactor, self)._build_unit_constraints()
        b = self

        equip = b.equip

        def mass_balance(equip, c, *vs):
            b = equip.parent_block()
            return b.flow_in[c, vs] * (1 - b.performance[c]) ==\
                b.flow_out[c, vs]
        equip.mass_balance = Constraint(self.comps, rule=mass_balance)

        if b.min_flow > 0:
            # If there is a minimum flow, apply it
            def min_flow(equip, *vs):
                b = equip.parent_block()
                return b.flow_in['water', vs] >= b.min_flow * b.equip_exists
            equip.min_flow = Constraint(rule=min_flow)

        def lin_unit_cost(equip, *vs):
            b = equip.parent_block()
            return b.lCTU[ne(vs)] >= b.beta * b.flow_in['water', vs] +\
                b.gamma * b.equip_exists
        equip.lin_unit_cost = Constraint(rule=lin_unit_cost)

        def nl_unit_cost(equip, *vs):
            b = equip.parent_block()
            return b.nCTU[ne(vs)] >= b.theta * b.flow_in['water', vs] ** 0.7
        equip.nonlin_unit_cost = Constraint(rule=nl_unit_cost)

        def max_flow(equip, c, *vs):
            b = equip.parent_block()
            return b.flow_in[c, vs] <= b.flow_in[c, vs].ub * b.equip_exists
        equip.max_flow = Constraint(self.comps, rule=max_flow)

    def _build_unit_ports(self):
        super(_Reactor, self)._build_unit_ports()
        b = self
        self.add_unit(InPort(name='inlet', parent=self, fc=b.flow_in))
        self.add_unit(OutPort(name='outlet', parent=self, fc=b.flow_out))

    def apply_linear_relaxations(self, nsegs=1, **kwargs):
        b = self
        lin_cuts = b.lin_cuts
        super(_Reactor, self).apply_linear_relaxations(nsegs=nsegs, **kwargs)

        def concave_unit_cost(flow_water):
            return b.theta * flow_water ** 0.7
        lin_cuts.env_unit_cost = Block()
        add_concave_linear_underest(lin_cuts.env_unit_cost, 'unit_cost', 1, SliceVar(b.flow_in, {1: 'water'}), b.nCTU, concave_unit_cost, exists=b.equip_exists, block_bounds=self.block_bounds)

    def reconstruct_envelopes(self):
        super(_Reactor, self).reconstruct_envelopes()
        # TODO below is a hack and should probably be handled with a block bound
        self.equip.max_flow.reconstruct()

    def add_oa_cut(self, iter_num):
        if abs(self.equip_exists.value) <= 1E-6:
            return  # skip if unit was not selected in this iteration
        add_oa_constraints(
            iter_num, self.oa_cuts, self.jacs,
            self.equip_exists)

    def get_flow_vars(self):
        yield self.flow_in
        yield self.flow_out

    def get_vars_to_bound(self):
        for c in self.comps - 'water':
            yield self.flow_in[c]
        yield self.flow_in['water'], self.block_bounds
        # yield self.flow_in['water']

    def generate_cut_gen_problem(self):
        if abs(self.equip_exists.value) <= 1E-3:
            # Do not solve cut generation problem for inactive units
            return None

        self.apply_NLP()
        num_vars = count_vars(self)

        b = ConcreteModel(name=self.local_name)
        var_set = b.var_set = RangeSet(num_vars)
        b.lbda = Var(b.var_set, domain=NonNegativeReals, bounds=(0, 1), initialize=0)
        b.lbda_sum = Constraint(expr=sum(b.lbda[vs] for vs in b.var_set) == 1)

        clone_block(self, b, var_set, b.lbda)

        sum_sq_diff = get_sum_sq_diff(self, b)
        b.obj = Objective(expr=sum_sq_diff, sense=minimize)
        return b

    def apply_self_proj_cut(self, cg_prob):
        if value(cg_prob.obj.expr) >= 1E-3:
            # generate cut
            print(self.unit_name, value(cg_prob.obj.expr))
            # self.cg_prob = cg_prob
            print('Adding self-projection cut for ' + self.local_name)
            self.lin_cuts.self_proj_cuts.add(
                2 * sum((cg_prob._flow_in_clone[c].value - self.flow_in[c].value) * (self.flow_in[c] - cg_prob._flow_in_clone[c].value) for c in self.comps) +
                2 * sum((cg_prob._flow_out_clone[c].value - self.flow_out[c].value) * (self.flow_out[c] - cg_prob._flow_out_clone[c].value) for c in self.comps) +
                2 * ((cg_prob._lCTU_clone.value - self.lCTU.value) * (self.lCTU - cg_prob._lCTU_clone.value)) +
                2 * ((cg_prob._nCTU_clone.value - self.nCTU.value) * (self.nCTU - cg_prob._nCTU_clone.value)) >= 0
            )
        else:
            print('Self-projection cut for ' +
                  self.local_name + ' redundant.')
