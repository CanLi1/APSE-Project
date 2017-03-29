from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import Var, NonNegativeReals
from six import iteritems
from idaes_models.core.ports import OutPort

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


@ProcBlock("Feed")
class _Feed(UnitModel):
    """ Simple feed node
    """

    def __init__(self, *args, **kwargs):
        """ Initialize feed """
        self.fix_flow = kwargs.pop('flow', None)
        self.fix_exists = kwargs.pop('exists', None)
        super(_Feed, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        super(_Feed, self)._build_unit_sets()

    def _build_unit_params(self):
        super(_Feed, self)._build_unit_params()

    def _build_unit_vars(self):
        super(_Feed, self)._build_unit_vars()
        b = self

        b.flow = Var(self.comps, domain=NonNegativeReals, initialize=1, bounds=(0, self.max_flow))
        if self.fix_flow is not None:
            for k, v in iteritems(self.fix_flow):
                b.flow[k].fix(v)
        b.total_flow_out = sum(b.flow[c] for c in self.comps)

    def _build_unit_constraints(self):
        super(_Feed, self)._build_unit_constraints()

    def _build_unit_ports(self):
        super(_Feed, self)._build_unit_ports()
        # for port in itervalues(self.unit_ports):
        #     port.construct_unit()
        # b.ports.links.add(expr=b.out_port == self.unit_ports.outlet.get_unit_link())
        b = self

        self.add_unit(OutPort(name='outlet', parent=self, fc=b.flow))

    def get_flow_vars(self):
        yield self.flow
