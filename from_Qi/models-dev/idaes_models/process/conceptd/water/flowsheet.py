"""
Flow-based model for Ruiz water network synthesis model
"""
from __future__ import division
from idaes_models.core.flowsheet_model import FlowsheetModel
from idaes_models.core import ProcBlock
from pyomo.environ import Objective, minimize
import os
import pandas as pd
from idaes_models.unit.water_net.feed import Feed
from idaes_models.unit.water_net.sink import Sink
from idaes_models.unit.water_net.reactor import Reactor
from six import itervalues, iteritems
from idaes_models.core.plugins.propagate_fixed import propagate_var_fix

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


@ProcBlock("WaterModel")
class _WaterModel(FlowsheetModel):
    """Represents the superstructure synthesis problem for water treatment.
    Holds the model objects for the problem.

    Attributes:
        oa_obj (TYPE): Description
        obj (TYPE): Description
    """
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('name', 'water_treatment_flowsheet')
        super(_WaterModel, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        super(_WaterModel, self).__call__(*args, **kwargs)

    def _activate_standard_objective(self):
        if not hasattr(self, 'obj'):
            total_cost = sum(o.CTU for o in itervalues(self.units) if isinstance(o, Reactor))
            self.obj = Objective(expr=total_cost, sense=minimize)
        else:
            self.obj.activate()
        for obj in self.component_objects(ctype=Objective, active=True):
            if obj.local_name != 'obj':
                obj.deactivate()

    def _activate_penalized_oa_objective(self):
        total_cost = sum(o.CTU for o in itervalues(self.units) if isinstance(o, Reactor))
        penalty_cost = sum(var for var in self._get_slack_vars())
        try:
            self.del_component('oa_obj')
        except:
            pass
        self.oa_obj = Objective(expr=total_cost + 1000 * penalty_cost)
        for obj in itervalues(
                self.component_map(ctype=Objective, active=True)):
            if obj.local_name != 'oa_obj':
                obj.deactivate()


def build_model(trial=None):
    """Builds the water treatment flowsheet model

    Args:
        trial (str, optional): trial number of case to be tested. For example,
        '01' or '02'. If none provided, default data case will be used.

    Returns:
        WaterModel: water treatment network model
    """
    # Load in input data from CSV files
    cur_path = os.path.dirname(os.path.realpath(__file__)) + '/'
    if trial is None:
        datDec = ''  # data decorator
        datPath = 'dat'  # data path
    else:
        datDec = trial
        datPath = 'dat/Trial' + datDec
    # inlet flows
    dat_inlet_flows = pd.read_csv(
        '{}{}/f{}.csv'.format(cur_path, datPath, datDec),
        header=0, index_col=0)
    in_streams = list(dat_inlet_flows.index)
    comps = list(dat_inlet_flows.columns.values)
    # unit performance data
    dat_tru_perf = pd.read_csv(
        '{}{}/cp{}.csv'.format(cur_path, datPath, datDec),
        header=0, index_col=0)
    trus = list(dat_tru_perf.index)
    # make sure that data for components is consistent
    assert comps == list(dat_tru_perf.columns.values), \
        "Components do not match up between inlet flow data file" \
        " and treatment unit performance file"
    # load treatment unit data
    dat_trus = pd.read_csv(
        '{}{}/trus{}.csv'.format(cur_path, datPath, datDec),
        header=0, index_col=0)
    assert trus == list(dat_trus.index), \
        "Treatment units do not match up between performance file and info file"
    # load max contaminant outflow data
    dat_max_cont = pd.read_csv(
        '{}{}/mc{}.csv'.format(cur_path, datPath, datDec),
        header=0, index_col=0)
    # treatment unit minimum flows
    min_unit_flows = {k: float(v) for k, v in iteritems(dict(dat_trus['min_flow']))}
    # treatment unit cost parameters
    theta = {k: float(v) for k, v in iteritems(dict(dat_trus['theta']))}
    beta = {k: float(v) for k, v in iteritems(dict(dat_trus['beta']))}
    gamma = {k: float(v) for k, v in iteritems(dict(dat_trus['gamma']))}

    m = WaterModel()
    m.comp_set = comps
    m.max_flow = 300

    # Define inlets and outlets, with corresponding splitters and mixer
    for i in in_streams:
        m.add_unit(Feed(name='in' + str(i), parent=m, flow=dat_inlet_flows.ix[i].to_dict()))
    m.add_unit(Sink(name='out', parent=m, max_cont=dat_max_cont.iloc[0].to_dict()))

    # Define treatment units, with corresponding mixers and splitters
    for i in trus:
        m.add_unit(Reactor(name='tru' + str(i), parent=m, comps=m.comps, perf=dat_tru_perf.ix[i].to_dict(), min_flow=min_unit_flows[i], theta=theta[i], beta=beta[i], gamma=gamma[i]))

    # Define connections between splitters and mixers
    for i in in_streams:
        for j in trus:
            m.connect(m.units['in' + str(i)], m.units['tru' + str(j)])
        m.connect(m.units['in' + str(i)], m.units.out)
    for i in trus:
        for j in trus:
            m.connect(m.units['tru' + str(i)], m.units['tru' + str(j)])
        m.connect(m.units['tru' + str(i)], m.units.out)

    m.build_units()
    m.build_links()
    m.expand_connectors()
    propagate_var_fix(m)

    return m
