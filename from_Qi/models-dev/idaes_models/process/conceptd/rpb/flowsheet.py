"""
Flow-based model for Packed Bed process
"""
from __future__ import division
from idaes_models.core.flowsheet_model import FlowsheetModel
from idaes_models.core import ProcBlock
import pyomo.environ as pe
from six import itervalues
from idaes_models.unit.standard.feed import Feed
from idaes_models.unit.rpb.packed_bed import PackedBed
from idaes_models.unit.rpb.packed_bed import RotatingPackedBed
from idaes_models.unit.rpb.sink import Sink
from six.moves import range
from itertools import product
from idaes_models.core.plugins.propagate_fixed import propagate_var_fix

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


@ProcBlock("PackedBedModel")
class _PackedBedModel(FlowsheetModel):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('name', 'methanol_flowsheet')
        super(_PackedBedModel, self).__init__(*args, **kwargs)

    def _activate_standard_objective(self):
        if not hasattr(self, 'obj'):
            equip_costs = sum(o.equip_cost for o in itervalues(self.units) if hasattr(o, 'equip_cost'))
            self.obj = pe.Objective(expr=equip_costs, sense=pe.minimize)
        else:
            self.obj.activate()
        for obj in self.component_objects(
                ctype=pe.Objective, active=True):
            if obj is not self.obj:
                obj.deactivate()

    def _activate_penalized_oa_objective(self):
        equip_costs = sum(o.equip_cost for o in itervalues(self.units) if hasattr(o, 'equip_cost'))
        penalty_cost = sum(var for var in self._get_slack_vars())
        try:
            self.del_component('oa_obj')
        except:
            pass
        self.oa_obj = pe.Objective(expr=equip_costs + 1000 * penalty_cost, sense=pe.minimize)
        for obj in self.component_objects(
                ctype=pe.Objective, active=True):
            if obj is not self.oa_obj:
                obj.deactivate()

def build_model(rpb=2, pb=2, recycle=True, crossover=False):
    m = PackedBedModel()
    m.comp_set = ['A', 'water']
    m.max_flow = 10

    vap_in = m.add_unit(
        Feed(name='vapor_in', parent=m, flow={'A': 0.0278, 'water': 2.78 - 0.0278}))
    liq_in = m.add_unit(
        Feed(name='liquid_in', parent=m, frac={'A': 0, 'water': 1}, bounds={'flow': (0, 0.0278)}))

    # outlet sinks
    m.add_unit(Sink(name='vapor_out', parent=m, bounds={'flow[A]': (None, 0.000278)}))
    m.add_unit(Sink(name='liquid_out', parent=m))

    for i in range(pb):
        n = i + 1
        pb = m.add_unit(PackedBed(name='packed_bed' + str(n), parent=m))
        m.connect(m.units.vapor_in, pb, to_port='vap_in')
        m.connect(m.units.liquid_in, pb, to_port='liq_in')
        if recycle:
            m.connect(pb, pb, from_port='vap_out', to_port='vap_in')
            m.connect(pb, pb, from_port='liq_out', to_port='liq_in')
        m.connect(pb, m.units.vapor_out, from_port='vap_out')
        m.connect(pb, m.units.liquid_out, from_port='liq_out')
    for j in range(rpb):
        n = j + 1
        rpb = m.add_unit(RotatingPackedBed(name='rotating_packed_bed' + str(n), parent=m, block_bounds={'total_flow_liq_in': (1E-5, None)}))
        m.connect(m.units.vapor_in, rpb, to_port='vap_in')
        m.connect(m.units.liquid_in, rpb, to_port='liq_in')
        if recycle:
            m.connect(rpb, rpb, from_port='vap_out', to_port='vap_in')
            m.connect(rpb, rpb, from_port='liq_out', to_port='liq_in')
        m.connect(rpb, m.units.vapor_out, from_port='vap_out')
        m.connect(rpb, m.units.liquid_out, from_port='liq_out')

    if crossover:
        for i, j in product(range(pb + rpb), repeat=2):
            if i == j:
                continue
            if i < pb:
                from_unit = 'packed_bed' + str(i + 1)
            else:
                from_unit = 'rotating_packed_bed' + str(i + 1 - pb)
            if j < pb:
                to_unit = 'packed_bed' + str(j + 1)
            else:
                to_unit = 'rotating_packed_bed' + str(j + 1 - pb)
            m.connect(m.units[from_unit], m.units[to_unit], from_port='vap_out', to_port='vap_in')
            m.connect(m.units[from_unit], m.units[to_unit], from_port='liq_out', to_port='liq_in')

    m.register_logic_constraint(all_of=(vap_in, liq_in))

    m.build_units()
    m.build_logic()
    m.build_links()
    m.expand_connectors()
    propagate_var_fix(m)

    return m
