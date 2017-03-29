"""
GLOA for packed bed vs rotating packed bed system

Original PB/RPB model developed by Qian, Zhi <zhiqian@andrew.cmu.edu>.
Manuscript in preparation.
"""
from __future__ import division
from six import itervalues
from idaes_models.process.conceptd.rpb.flowsheet import build_model

from idaes_models.core.util.misc import doNothing
from idaes_models.core.loa import do_LOA, do_GLOA

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def main():
    m = build_model(rpb=1, pb=1)

    for o in itervalues(m.units):
        o.apply_OA_strategy()
        getattr(o, 'apply_linear_relaxations', doNothing)(nsegs=3)

    # m.solvers.global_MINLP.MaxTime = 30
    # m.solve_global_MINLP()
    m.solvers.local_NLP.promote('ipopt')
    do_LOA(m, tol=1)
    m.solvers.global_NLP.MaxTime = 30
    m.solvers.self_proj_cut_gen.MaxTime = 30
    do_GLOA(m, tol=1, iterlim=4, contract_bounds=True)
    return m


if __name__ == '__main__':
    m = main()
