"""
Very simple demonstartion of using a Pyomo ExternalFunction.
"""
from __future__ import division
from __future__ import print_function
__author__ = "John Eslick"
__version__ = "1.0.0"

from pyomo.environ import *
import sys
model = ConcreteModel()
cset = list(range(100000))

model.T = Var(cset)
model.P = Var(initialize=101325)
model.V = Var(cset)
model.V_func = ExternalFunction(library="/home/jovyan/models/src/prop_lib/phys_prop_C_auto.so", function="V_vap_py_first_five_args")

model.y = Var(cset)
"""
    T = ra[at[0]]
    P = ra[at[1]]
    y_CO2 = y0 = ra[at[2]]
    y_H2O = y1 = ra[at[3]]
    y_N2 = y2 = ra[at[4]]
    y_MEA = y3 = ra[at[5]]
    y_O2 = y4 = ra[at[6]]
"""
#model.cp_ideal_func = ExternalFunction(library="/home/jovyan/models/models-fork/models/src/prop_lib/phys_prop_C_auto.so", function="cp_ideal")

for i in cset:
    model.T[i].fix(300 + i/1000.0)

def V_rule(model, i):
    return model.V[i]==model.V_func(model.T[i], model.P)

model.eq = Constraint(cset, rule=V_rule)
model.o = Objective(expr=model.V[0])
solver = SolverFactory('ipopt')
results = solver.solve(model, tee=True)
print(results)