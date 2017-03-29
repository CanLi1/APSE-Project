"""
Flow-based model for Methanol process
"""
from __future__ import division
from idaes_models.core.flowsheet_model import FlowsheetModel
from idaes_models.core import ProcBlock
import pyomo.environ as pe
from six import itervalues
from bunch import Bunch
from idaes_models.unit.methanol.feed import Feed
from idaes_models.unit.methanol.sink import Sink
from idaes_models.unit.methanol.compressor import Compressor
from idaes_models.unit.methanol.hexchanger import Cooler, Heater
from idaes_models.unit.methanol.reactor import Reactor
from idaes_models.unit.methanol.valve import Valve
from idaes_models.unit.methanol.flash import Flash
from idaes_models.unit.standard.single_choice import SingleChoiceMixer, SingleChoiceNode, SingleChoiceSplitter
from idaes_models.core.plugins.propagate_fixed import propagate_var_fix

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


@ProcBlock("MethanolModel")
class _MethanolModel(FlowsheetModel):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('name', 'methanol_flowsheet')
        super(_MethanolModel, self).__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        super(_MethanolModel, self).__call__(*args, **kwargs)

    def _activate_standard_objective(self):
        if not hasattr(self, 'obj'):
            equip_costs = sum(o.equip_cost for o in itervalues(self.units) if hasattr(o, 'equip_cost'))
            equip_profit = sum(o.equip_profit for o in itervalues(self.units) if hasattr(o, 'equip_profit'))
            self.obj = pe.Objective(expr=equip_costs - equip_profit, sense=pe.minimize)
        else:
            self.obj.activate()
        for obj in itervalues(
                self.component_map(ctype=pe.Objective, active=True)):
            if obj.local_name != 'obj':
                obj.deactivate()

    def _activate_penalized_oa_objective(self):
        equip_costs = sum(o.equip_cost for o in itervalues(self.units) if hasattr(o, 'equip_cost'))
        equip_profit = sum(o.equip_profit for o in itervalues(self.units) if hasattr(o, 'equip_profit'))
        penalty_cost = sum(var for var in self._get_slack_vars())
        try:
            self.del_component('oa_obj')
        except:
            pass
        self.oa_obj = pe.Objective(expr=equip_costs - equip_profit + 1000 * penalty_cost, sense=pe.minimize)
        for obj in self.component_objects(
                ctype=pe.Objective, active=True):
            if obj.local_name != 'oa_obj':
                obj.deactivate()

    def display_costs(self):
        for o in itervalues(self.units):
            cost = getattr(o, 'equip_cost', None)
            exists = getattr(o, 'equip_exists', Bunch(value=1))
            if cost is not None:
                if abs(exists.value - 1.0) <= 1E-6:
                    print(o.name, cost.value)

    def display_slacks(self):
        def compare(x, y):
            if x.value > y.value:
                return 1
            elif x.value < y.value:
                return -1
            else:
                return 0
        print('Active slack penalties')
        for v in sorted(filter(lambda x: x.value > 0, self._get_slack_vars()), cmp=compare):
            print(v.name, v.value)


def build_model():
    """Constructs and populates the methanol model

    Returns:
        FlowsheetModel: methanol flowsheet model
    """
    m = MethanolModel()
    m.comp_set = ['H2', 'CO', 'MeOH', 'CH4']
    m.max_flow = 20
    m.min_T = 300
    m.max_T = 900
    m.max_P = 150
    antoine = {
        'A': {'H2': 13.6333, 'CO': 14.3686, 'MeOH': 18.5875, 'CH4': 15.2243},
        'B': {'H2': 164.9, 'CO': 530.22, 'MeOH': 3626.55, 'CH4': 897.84},
        'C': {'H2': 3.19, 'CO': -13.15, 'MeOH': -34.29, 'CH4': -7.16}}
    rxn_stoich = {'H2': -1.0, 'CO': -0.5, 'MeOH': 0.5, 'CH4': 0.0}

    # Inlets
    cfeed = m.add_unit(
        Feed(name='cheap_feed', parent=m, price=795.6, T=300, P=10,
             frac={'H2': .6, 'CO': .25, 'MeOH': 0, 'CH4': 0.15}))
    efeed = m.add_unit(
        Feed(name='expensive_feed', parent=m, price=1009.8, T=300, P=10,
             frac={'H2': .65, 'CO': .3, 'MeOH': 0, 'CH4': 0.05}))
    # Outlets
    prod = m.add_unit(
        Sink(name='product', parent=m, price=7650, T=400, min_frac={'MeOH': 0.9}))
    purge = m.add_unit(Sink(name='purge', parent=m, price=642.6, T=400))
    # Add feed mixer
    feed_mix = m.add_unit(SingleChoiceNode(name='feed_to_compress', parent=m))
    m.connect(cfeed, feed_mix)
    m.connect(efeed, feed_mix)
    # Feed compressors
    ss_feed_compr = m.add_unit(
        Compressor(name='ss_feed_compr', parent=m,
                   bounds={'flow': (0, 5), 'Pout': (25, None)},
                   block_bounds={'work': (3.3556, 3.7407),
                                 'P_ratio': (1.7400, 1.7400),
                                 'Pin': (None, 10),
                                 'Pout': (None, 115),
                                 'Tin': (None, 300),
                                 'Tout': (None, 600)
                                 }))
    ms_feed_compr1 = m.add_unit(
        Compressor(name='ms_feed_compr1', parent=m, compr_cost=25, bounds={'flow': (0, 5)},
                   block_bounds={'work': (1.7723, 6.6116),
                                 'P_ratio': (1.4165, 1.4341),
                                 'Pin': (None, 10),
                                 'Pout': (None, 115),
                                 'Tin': (None, 300),
                                 'Tout': (None, 500)}))
    ms_feed_cooler = m.add_unit(
        Cooler(name='ms_feed_cooler', parent=m, bounds={'Q': (0, 10), 'flow': (0, 5), 'P': (0, 115),
               'Tout': (300, 500)}))
    ms_feed_compr2 = m.add_unit(
        Compressor(name='ms_feed_compr2', parent=m, compr_cost=25, bounds={'flow': (0, 5)},
                   block_bounds={'work': (1.7723, 6.6116),
                                 'P_ratio': (1.2873, 1.3027),
                                 'Pin': (None, 115),
                                 'Pout': (None, 150),
                                 'Tin': (None, 500),
                                 'Tout': (None, 600)}))
    m.connect(feed_mix, ss_feed_compr)
    m.connect(feed_mix, ms_feed_compr1)
    m.connect(ms_feed_compr1, ms_feed_cooler)
    m.connect(ms_feed_cooler, ms_feed_compr2)
    # Feed compr mixer
    feed_compr_mix = m.add_unit(SingleChoiceMixer(name='feed_compr_mixer', parent=m))
    m.connect(ss_feed_compr, feed_compr_mix)
    m.connect(ms_feed_compr2, feed_compr_mix)
    # Recycle mixer
    # preheat_split = m.add_unit(
    #     'preheat_split',
    #     SingleChoiceSplitter(comps=m.comps,
    #                          bounds={'T_feed_compr_mixer': (350, None),
    #                          'P': (25, None)}))
    # Reaction preheater
    rxn_preheat = m.add_unit(
        Heater(name='rxn_preheat', parent=m, bounds={'Tout': (0, 500), 'Q': (0, 10)}))
    m.connect(feed_compr_mix, rxn_preheat)
    # Preheat to
    # Reaction splitter
    rxn_split = m.add_unit(SingleChoiceSplitter(name='rxn_split', parent=m))
    m.connect(rxn_preheat, rxn_split)
    # Reactors
    crxn = m.add_unit(
        Reactor(name='cheap_rxn', parent=m, stoich=rxn_stoich, key_comp='H2', v_param=0.05,
                vol_cost=5, rxn_cost=100, bounds={
                    'conv': (0, 0.415), 'eq_conv': (0, 0.415), 'extent': (0, 5),
                    'flow_in': (0, 10), 'Pin': (0, 150)},
                block_bounds={'flow_out[CO]': (0.1, None)}
                ))
    erxn = m.add_unit(
        Reactor(name='expensive_rxn', parent=m, stoich=rxn_stoich, key_comp='H2', v_param=0.1,
                vol_cost=10, rxn_cost=250, bounds={
                    'conv': (0, 0.415), 'eq_conv': (0, 0.415), 'extent': (0, 5),
                    'flow_in': (0, 10), 'Pin': (0, 150)},
                block_bounds={'flow_out[CO]': (0.1, None)}
                ))
    m.connect(rxn_split, crxn)
    m.connect(rxn_split, erxn)
    # Reactor mix
    rxn_mix = m.add_unit(SingleChoiceMixer(name='rxn_mix', parent=m))
    m.connect(crxn, rxn_mix)
    m.connect(erxn, rxn_mix)
    # Valve
    valve = m.add_unit(
        Valve(name='valve', parent=m, gamma=0.23077, bounds={'Tin': (523, None)}))
    m.connect(rxn_mix, valve)
    # Flash cooler
    flash_cool = m.add_unit(
        Cooler(name='flash_cooler', parent=m, bounds={'Q': (0, 1.5185)},
               block_bounds={'Q': (None, 1.5)}))  # this was 1.47
    m.connect(valve, flash_cool)
    # Flash
    flash = m.add_unit(
        Flash(name='flash', parent=m, antoine=antoine, key_comp='H2',
              bounds={'vap_recover[MeOH]': (None, 0.9),
                      'flow_liq[MeOH]': (0.09, 0.9)}))
    m.connect(flash_cool, flash)
    # Product Heater
    prod_heat = m.add_unit(Heater(name='prod_heater', parent=m, bounds={'Q': (0, 10)}))
    m.connect(flash, prod_heat, from_port='liq_out')
    m.connect(prod_heat, prod)
    # Purge heater
    purge_heat = m.add_unit(Heater(name='purge_heater', parent=m, bounds={'Q': (0, 10)}))
    m.connect(flash, purge_heat, from_port='vap_out')
    m.connect(purge_heat, purge)
    # Recycle compressor split
    recy_compr_split = m.add_unit(SingleChoiceSplitter(name='recy_compr_split', parent=m))
    m.connect(flash, recy_compr_split, from_port='vap_out')
    # Recycle compressors
    ss_recy_compr = m.add_unit(
        Compressor(name='ss_recy_compr', parent=m, bounds={'flow': (0, 5)},
                   block_bounds={'work': (1.0434, 4.3386),
                                 'P_ratio': (1.0779, 50),
                                 'Pin': (None, 135),
                                 'Pout': (None, 150),
                                 'Tin': (None, 400),
                                 'Tout': (None, 500)}))
    ms_recy_compr1 = m.add_unit(
        Compressor(name='ms_recy_compr1', parent=m, compr_cost=25,
                   bounds={'flow': (0, 5), 'Pout': (0, 135)},
                   block_bounds={'work': (0.9467, 2.1692),
                                 'P_ratio': (1.0429, 25),
                                 'Pin': (None, 135),
                                 'Pout': (None, 150),
                                 'Tin': (None, 400),
                                 'Tout': (None, 450)}))
    ms_recy_cooler = m.add_unit(
        Cooler(name='ms_recy_cooler', parent=m, bounds={'Q': (0, 10), 'flow': (0, 5), 'P': (0, 130)}))
    ms_recy_compr2 = m.add_unit(
        Compressor(name='ms_recy_compr2', parent=m, compr_cost=25,
                   bounds={'flow': (0, 5)},
                   block_bounds={'work': (0.9467, 2.1692),
                                 'P_ratio': (1.0411, 25),
                                 'Pin': (None, 135),
                                 'Pout': (None, 150),
                                 'Tin': (None, 450),
                                 'Tout': (None, 500)}))
    m.connect(recy_compr_split, ss_recy_compr)
    m.connect(recy_compr_split, ms_recy_compr1)
    m.connect(ms_recy_compr1, ms_recy_cooler)
    m.connect(ms_recy_cooler, ms_recy_compr2)
    # Recycle compressor mixer
    recy_compr_mix = m.add_unit(SingleChoiceMixer(name='recy_compr_mix', parent=m))
    m.connect(ss_recy_compr, recy_compr_mix)
    m.connect(ms_recy_compr2, recy_compr_mix)
    m.connect(recy_compr_mix, rxn_preheat)

    # cheap or expensive feed
    m.register_logic_constraint(one_of=(cfeed, efeed))
    # single or multi-stage feed compression
    m.register_logic_constraint(one_of=(ss_feed_compr, ms_feed_compr1))
    m.register_logic_constraint(one_of=(ss_feed_compr, ms_feed_cooler))
    m.register_logic_constraint(one_of=(ss_feed_compr, ms_feed_compr2))
    m.register_logic_constraint(one_of=(crxn, erxn))  # cheap or expensive reactor
    # single or multi-stage recycle compression
    m.register_logic_constraint(one_of=(ss_recy_compr, ms_recy_compr1))
    m.register_logic_constraint(one_of=(ss_recy_compr, ms_recy_cooler))
    m.register_logic_constraint(one_of=(ss_recy_compr, ms_recy_compr2))

    m.build_units()
    m.build_logic()
    m.build_links()
    m.expand_connectors()
    propagate_var_fix(m)

    return m
