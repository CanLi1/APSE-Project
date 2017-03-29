"""
Global logic-based outer approximation
"""
from __future__ import division
from idaes_models.process.conceptd.methanol.flowsheet import build_model
# import pandas as pd
from six import itervalues, iteritems
from pprint import pprint
from idaes_models.core.loa import do_LOA, do_GLOA

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def main():
    """Main function to execute run"""
    m = build_model()

    m.solvers.global_MINLP.MaxTime = 43200  # 12 hrs
    # Uncomment this line to solve the full space problem in BARON
    # m.solve_global_MINLP()

    for o in itervalues(m.units):
        o.apply_linear_relaxations()
        o.apply_OA_strategy(oa_ports=False)

    # m.solvers.local_NLP.options["outlev"] = 3  # use for troubleshooting CONOPT
    m.solvers.local_NLP.promote('ipopt')
    m.solvers.local_NLP.ipopt.tol = 1E-4
    # do_LOA(m, tol=10, custom_init=manual_init_LOA)
    do_LOA(m, tol=10)

    # Run rigorous GLOA algorithm
    m.solvers.global_NLP.MaxTime = 1200
    m.solvers.self_proj_cut_gen.MaxTime = 30
    do_GLOA(m, tol=10, iterlim=17)

    pprint(dict(m.solver_time))
    print("Total CPU time", m.total_solver_time)

    return m

def manual_init_LOA(m, init_index=1):
    """Manual set covering initialization for LOA"""

    if not solve_init_NLP(m, init_index):
        raise RuntimeError('Unable to solve initialization local NLP 1')
    m.add_oa_cut()
    m.add_integer_cut(tmp=True)

    if not solve_init_NLP(m, -1 * init_index):
        raise RuntimeError('Unable to solve initialization local NLP 2')
    m.add_oa_cut()
    m.add_integer_cut(tmp=True)

def solve_init_NLP(m, init_index=1, tee=False):
    keyEquip = set(['cheap_feed', 'expensive_feed', 'cheap_rxn', 'expensive_rxn', 'ms_feed_compr1', 'ms_feed_compr2', 'ms_feed_cooler', 'ss_feed_compr', 'ms_recy_compr1', 'ms_recy_compr2', 'ms_recy_cooler', 'ss_recy_compr'])
    print('Begin LOA initialization')
    nlpInit = {
        1: set(['cheap_feed', 'expensive_rxn', 'ms_feed_compr1', 'ms_feed_compr2', 'ms_feed_cooler', 'ms_recy_compr1', 'ms_recy_compr2', 'ms_recy_cooler']),
        2: set(['cheap_feed', 'cheap_rxn', 'ss_feed_compr', 'ss_recy_compr']),
        3: set(['cheap_feed', 'cheap_rxn', 'ss_feed_compr', 'ms_recy_compr1', 'ms_recy_compr2', 'ms_recy_cooler']),
        4: set(['cheap_feed', 'cheap_rxn', 'ms_feed_compr1', 'ms_feed_compr2', 'ms_feed_cooler', 'ss_recy_compr']),
        5: set(['cheap_feed', 'cheap_rxn', 'ms_feed_compr1', 'ms_feed_compr2', 'ms_feed_cooler', 'ms_recy_compr1', 'ms_recy_compr2', 'ms_recy_cooler']),
        6: set(['cheap_feed', 'expensive_rxn', 'ss_feed_compr', 'ss_recy_compr']),
        7: set(['cheap_feed', 'expensive_rxn', 'ss_feed_compr', 'ms_recy_compr1', 'ms_recy_compr2', 'ms_recy_cooler']),
        8: set(['cheap_feed', 'expensive_rxn', 'ms_feed_compr1', 'ms_feed_compr2', 'ms_feed_cooler', 'ss_recy_compr'])}
    if init_index > 0:
        for name, o in iteritems(m.units):
            if name not in (keyEquip - nlpInit[init_index]):
                exists = getattr(o, 'equip_exists', None)
                if exists is not None:
                    exists.value = 1
            else:
                exists = getattr(o, 'equip_exists', None)
                if exists is not None:
                    exists.value = 0
    else:
        init_index = -1 * init_index
        for name, o in iteritems(m.units):
            if name not in keyEquip or name in (keyEquip - nlpInit[init_index]):
                exists = getattr(o, 'equip_exists', None)
                if exists is not None:
                    exists.value = 1
            else:
                exists = getattr(o, 'equip_exists', None)
                if exists is not None:
                    exists.value = 0
    return m.solve_local_NLP(tee=tee)


if __name__ == '__main__':
    m = main()
