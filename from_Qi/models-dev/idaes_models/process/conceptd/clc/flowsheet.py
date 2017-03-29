from __future__ import division
from idaes_models.core.flowsheet_model import FlowsheetModel
from idaes_models.core import ProcBlock
from pyomo.environ import Objective

@ProcBlock("CLCFlowsheet")
class _CLCFlowsheet(FlowsheetModel):

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('name', 'clc_flowsheet')
        super(_CLCFlowsheet, self).__init__(*args, **kwargs)

    def _activate_standard_objective(self):
        self.del_component('obj')
        self.obj = Objective(expr=-1 * self.units.fuel_reactor.gas_out.F_unit)
        for obj in self.component_objects(
                ctype=Objective, active=True):
            if obj is not self.obj:
                obj.deactivate()

    def _activate_penalized_oa_objective(self):
        self.del_component('oa_obj')
        self.oa_obj = Objective(expr=-1 * self.units.fuel_reactor.gas_out.F_unit + sum(var for var in self._get_slack_vars()))
        for obj in self.component_objects(
                ctype=Objective, active=True):
            if obj is not self.oa_obj:
                obj.deactivate()
