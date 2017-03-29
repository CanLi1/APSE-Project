from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import PositiveReals, NonNegativeReals, Param, Var, Binary, Constraint
from idaes_models.core.ports import OutPort


@ProcBlock("Feed")
class _Feed(UnitModel):
    def __init__(self, *args, **kwargs):
        """ Initialize feed """
        kwargs.setdefault('type', [])
        kwargs['type'].append('Feed')
        self._price = kwargs.pop('price', None)
        self.fix_T = kwargs.pop('T', None)
        self.fix_P = kwargs.pop('P', None)
        self.fix_frac = kwargs.pop('frac', {})
        super(_Feed, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        super(_Feed, self)._build_unit_sets()

    def _build_unit_params(self):
        super(_Feed, self)._build_unit_params()
        # Feed cost
        b = self
        b.price = Param(initialize=self._price)

    def _build_unit_vars(self):
        super(_Feed, self)._build_unit_vars()
        b = self
        # Temperature
        T = b.T = Var(domain=PositiveReals, initialize=300, bounds=(300, 900))
        if self.fix_T is not None:
            T.fix(self.fix_T)
        # Pressure
        P = b.P = Var(domain=PositiveReals, initialize=1, bounds=(0, 150))
        if self.fix_P is not None:
            P.fix(self.fix_P)
        # Flow
        b.flow = Var(self.comps, domain=NonNegativeReals,
                     initialize=1, bounds=(0, 5))
        # Cost
        b.equip_cost = Var(domain=NonNegativeReals,
                           initialize=0, bounds=(0, 6000))
        # Composition
        b.frac = Var(self.comps, domain=NonNegativeReals, bounds=(0, 1))
        for c in self.comps:
            if c in self.fix_frac:
                b.frac[c].fix(self.fix_frac[c])
            else:
                b.frac[c].value = 1.0 / len(self.comps)
        # Total flow
        b.total_flow_out = Var(domain=NonNegativeReals,
                               initialize=1, bounds=(0, 20))
        # Feed existence
        b.equip_exists = Var(domain=Binary, initialize=1)

    def _build_unit_constraints(self):
        super(_Feed, self)._build_unit_constraints()
        b = self

        equip = b.equip

        def total_flow_out_defn(equip):
            return b.total_flow_out == sum(b.flow[c] for c in self.comps)
        equip.total_flow_out_defn = Constraint(rule=total_flow_out_defn)

        def total_flow_out_lb(equip):
            return 0.5 * b.equip_exists <= b.total_flow_out
        equip.total_flow_out_lb = Constraint(rule=total_flow_out_lb)

        def frac_flow(equip, c):
            return b.flow[c] == b.total_flow_out * b.frac[c]
        equip.frac_flow = Constraint(self.comps, rule=frac_flow)

        def cost_calc(equip):
            return b.equip_cost == b.price * b.total_flow_out
        equip.cost_calc = Constraint(rule=cost_calc)

        def flow_ub(equip, c):
            return b.flow[c] <= b.flow[c].ub * b.equip_exists
        equip.flow_ub = Constraint(self.comps, rule=flow_ub)

        def cost_ub(equip):
            return b.equip_cost <= b.equip_cost.ub * b.equip_exists
        equip.cost_ub = Constraint(rule=cost_ub)

    def _build_unit_ports(self):
        super(_Feed, self)._build_unit_ports()
        b = self
        self.add_unit(OutPort(name='outlet', parent=self,
                              fc=b.flow, T=b.T, P=b.P))

    def get_flow_vars(self):
        yield self.flow

    def generate_cut_gen_problem(self):
        pass  # empty placeholder
