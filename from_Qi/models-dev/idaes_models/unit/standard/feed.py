from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import NonNegativeReals, Param, Var, Binary, Constraint
from idaes_models.core.ports import OutPort


@ProcBlock("Feed")
class _Feed(UnitModel):
    """Indicates a feed stream to the process.

    TOOD: this unit still has some relics of its promotion from a specific unit
    model to core that need to be resolved.

    For example:
    - cost_ub

    TODO: There should be some thought about how Feed and Sink will behave when
    nesting models withing each other.

    This unit model also does not support self projection cuts

    Attributes:
        fix_flow (dict): Dictionary enabling flow rates of specified components to be held constant
        fix_frac (dict): Dictionary enabling mole fractions of specified components to be held constant
        fix_total_flow (float): Allows the total flow rate to be held constant
        price (float): Assigns a cost to the total molar flowrate
    """
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('type', [])
        kwargs['type'].append('Feed')
        self._price = kwargs.pop('price', 0)
        self.fix_frac = kwargs.pop('frac', {})
        self.fix_total_flow = kwargs.pop('total_flow', None)
        self.fix_flow = kwargs.pop('flow', {})
        super(_Feed, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        super(_Feed, self)._build_unit_sets()

    def _build_unit_params(self):
        super(_Feed, self)._build_unit_params()
        b = self
        # Feed cost
        b.price = Param(initialize=self._price)

    def _build_unit_vars(self):
        super(_Feed, self)._build_unit_vars()
        b = self
        # Flow
        b.flow = Var(self.comps, domain=NonNegativeReals, initialize=1, bounds=(0, self._max_flow))
        for c in self.fix_flow:
            b.flow[c].fix(self.fix_flow[c])
        # Cost
        b.equip_cost = Var(domain=NonNegativeReals, initialize=0, bounds=(0, 6000))
        # Composition
        b.frac = Var(self.comps, domain=NonNegativeReals, bounds=(0, 1))
        for c in self.comps:
            if c in self.fix_frac:
                b.frac[c].fix(self.fix_frac[c])
            else:
                b.frac[c].value = 1.0 / len(self.comps)
        # Total flow
        b.total_flow_out = Var(domain=NonNegativeReals, initialize=1, bounds=(0, self._max_flow))
        if all(c in self.fix_flow for c in self.comps):
            # fix_flow contains all components
            self.fix_total_flow = sum(self.fix_flow[c] for c in self.comps)
        if self.fix_total_flow is not None:
            b.total_flow_out.fix(self.fix_total_flow)
        # Feed existence
        b.equip_exists = Var(domain=Binary, initialize=1)

    def _build_unit_constraints(self):
        super(_Feed, self)._build_unit_constraints()
        b = self

        equip = b.equip

        def total_flow_out_defn(equip):
            return b.total_flow_out == sum(b.flow[c] for c in self.comps)
        equip.total_flow_out_defn = Constraint(rule=total_flow_out_defn)

        if self.fix_total_flow is None and \
                len(self.fix_frac) < len(self.comps) - 1:
            raise NotImplementedError('Nonlinear feed not supported')

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
        self.add_unit(OutPort(name='outlet', parent=self, fc=b.flow))

    def __call__(self, b, indx=0):
        super(_Feed, self).__call__(b, indx)

    def get_flow_vars(self):
        yield self.flow

    def generate_cut_gen_problem(self):
        pass  # empty placeholder
