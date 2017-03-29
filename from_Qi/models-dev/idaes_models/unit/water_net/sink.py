from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import Var, Constraint, NonNegativeReals
from idaes_models.core.ports import InPort

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


@ProcBlock("Sink")
class _Sink(UnitModel):
    """ Simple sink node
    """

    def __init__(self, *args, **kwargs):
        """ Initialize sink """
        self.max_cont = kwargs.pop('max_cont', None)
        super(_Sink, self).__init__(*args, **kwargs)
        if self.max_cont is None:
            self.max_cont = {c: 0 for c in self.comps if c != 'water'}

    def _build_unit_sets(self):
        super(_Sink, self)._build_unit_sets()

    def _build_unit_params(self):
        super(_Sink, self)._build_unit_params()

    def _build_unit_vars(self):
        super(_Sink, self)._build_unit_vars()
        b = self

        b.flow = Var(self.comps, domain=NonNegativeReals, initialize=1, bounds=(0, self._max_flow))
        b.total_flow_out = sum(b.flow[c] for c in self.comps)

    def _build_unit_constraints(self):
        super(_Sink, self)._build_unit_constraints()
        b = self

        def max_contaminant(equip, c, *vs):
            b = equip.parent_block()
            return b.flow[c, vs] <= self.max_cont[c] if c != 'water' else Constraint.Skip
        b.equip.max_contaminant = Constraint(self.comps, rule=max_contaminant)

    def _build_unit_ports(self):
        super(_Sink, self)._build_unit_ports()
        b = self
        self.add_unit(InPort(name='inlet', parent=self, fc=b.flow))

    def get_flow_vars(self):
        yield self.flow

class ConcSink(UnitModel):
    """ Simple sink node
    """

    def __init__(self, comp, *args, **kwargs):
        """ Initialize sink """
        super(ConcSink, self).__init__(*args, **kwargs)
        self.comps = comp

    def __call__(self, b, indx=0):
        super(ConcSink, self).__call__(b, indx)
        b.conc = Var(self.comps - 'water', domain=NonNegativeReals, initialize=1, bounds=(0, 1))
        b.total_flow = Var(domain=NonNegativeReals)
        b.total_flow.fix(60)
        b.flow = Var(self.comps - 'water', domain=NonNegativeReals, bounds=(0, self._max_flow))

        # --------------
        # - Connectors -
        # --------------

        from connectors import conc_connector
        b.in_port = conc_connector(b.total_flow, b.conc)

        max_cont = {1: 3, 2: 3, 3: 3, 4: 3}

        def max_contaminant(_, c):
            return b.flow[c] <= max_cont[int(c)] if c != 'water' else Constraint.Skip
        b.max_contaminant = Constraint(self.comps, rule=max_contaminant)

        equip = b.equip

        def flow_calc(equip, c):
            return b.flow[c] == b.total_flow * b.conc[c]
        equip.flow_calc = Constraint(self.comps - 'water', rule=flow_calc)
