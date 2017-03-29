from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import NonNegativeReals, Var
from idaes_models.core.ports import InPort

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


@ProcBlock("Sink")
class _Sink(UnitModel):
    """ Simple sink node
    """

    def __init__(self, *args, **kwargs):
        """ Initialize sink """
        super(_Sink, self).__init__(*args, **kwargs)

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

    def _build_unit_ports(self):
        super(_Sink, self)._build_unit_ports()
        b = self
        self.add_unit(InPort(name='inlet', parent=self, fc=b.flow))

    def __call__(self, b, indx=0):
        super(_Sink, self).__call__(b, indx)

    def get_flow_vars(self):
        yield self.flow
