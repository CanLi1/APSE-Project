"""
Global logic-based outer approximation
"""
from __future__ import division
from pprint import pprint
from idaes_models.process.conceptd.water.flowsheet import build_model
from six import itervalues
from idaes_models.core.loa import do_LOA, do_GLOA

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def main():
    """Builds a flowsheet model and solves it

    Returns:
        WaterModel: water treatment network model
    """
    m = build_model()

    for o in itervalues(m.units):
        o.apply_linear_relaxations()
        o.apply_OA_strategy(oa_ports=True)

    m.solvers.global_MINLP.MaxTime = 14400
    m.solvers.global_MINLP.EpsA = 100
    # m.solve_global_MINLP()

    # Do non-rigorous LOA in order to initialize GLOA
    m.solvers.local_NLP.promote('conopt')
    # m.solvers.local_NLP.outlev = 3  # use for troubleshooting CONOPT
    # m.solvers.local_NLP.halt_on_ampl_error = 'yes'
    # m.solvers.local_NLP.ipopt.tol = 1E-4
    do_LOA(m, tol=100)

    # Run rigorous GLOA algorithm
    do_GLOA(m, tol=100, iterlim=32, contract_bounds=True, do_self_proj=False)

    pprint(dict(m.solver_time))
    print("Total CPU time", m.total_solver_time)

    return m


if __name__ == '__main__':
    m = main()
