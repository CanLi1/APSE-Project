from __future__ import division
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from idaes_models.core.ports import InPort, OutPort
from pyomo.environ import Var, NonNegativeReals

__author__ = "Qi Chen <andrew.cmu.edu>"

@ProcBlock("SingleChoiceNode")
class _SingleChoiceNode(UnitModel):
    """Superstructure node that functions as a SingleChoiceMixer and/or SingleChoiceSplitter
    """
    def __init__(self, *args, **kwargs):
        self.single_choice_inlet = kwargs.pop('single_choice_inlet', True)
        self.single_choice_outlet = kwargs.pop('single_choice_outlet', True)
        self._flow_scheme = 'fc'
        super(_SingleChoiceNode, self).__init__(*args, **kwargs)

    def _build_unit_vars(self):
        super(_SingleChoiceNode, self)._build_unit_vars()
        self.P = Var(domain=NonNegativeReals, initialize=1, bounds=(0, self.max_P))
        self.T = Var(domain=NonNegativeReals, initialize=300, bounds=(self.min_T, self.max_T))
        self.fc = Var(self.comps, domain=NonNegativeReals, initialize=1, bounds=(0, self.max_flow))
        self.F = Var(domain=NonNegativeReals, initialize=1.0, bounds=(0, self.max_flow))
        self.y = Var(self.comps, domain=NonNegativeReals, initialize=0.5, bounds=(0, 1))
        self.vf = Var(domain=NonNegativeReals, initialize=0.5, bounds=(0, 1))

    def _build_unit_ports(self):
        super(_SingleChoiceNode, self)._build_unit_ports()
        if self._flow_scheme == 'fc':
            # individual component flows
            self.add_unit(InPort(name='inlet', parent=self, single_choice=self.single_choice_inlet, P=self.P, T=self.T, fc=self.fc, vf=self.vf))
            self.add_unit(OutPort(name='outlet', parent=self, single_choice=self.single_choice_outlet, P=self.P, T=self.T, fc=self.fc, vf=self.vf))
        elif self._flow_scheme == 'Fy':
            # total flows and compositions
            self.add_unit(InPort(name='inlet', parent=self, single_choice=self.single_choice_inlet, P=self.P, T=self.T, F=self.F, y=self.y, vf=self.vf))
            self.add_unit(OutPort(name='outlet', parent=self, single_choice=self.single_choice_outlet, P=self.P, T=self.T, F=self.F, y=self.y, vf=self.vf))
        else:
            raise ValueError('Unknown flow scheme: ' + self._flow_scheme)

    def fix_flows(self):
        pass

    def unfix_flows(self):
        pass

    def get_flow_vars(self):
        return iter([])


@ProcBlock("SingleChoiceMixer")
class _SingleChoiceMixer(_SingleChoiceNode):
    """Superstructure node where only a single inlet will be active
    """
    def __init__(self, *args, **kwargs):
        kwargs['single_choice_inlet'] = True
        kwargs['single_choice_outlet'] = False
        super(_SingleChoiceMixer, self).__init__(*args, **kwargs)

@ProcBlock("SingleChoiceSplitter")
class _SingleChoiceSplitter(_SingleChoiceNode):
    """Superstructure node where only a single outlet will be active
    """
    def __init__(self, *args, **kwargs):
        kwargs['single_choice_inlet'] = False
        kwargs['single_choice_outlet'] = True
        super(_SingleChoiceSplitter, self).__init__(*args, **kwargs)
