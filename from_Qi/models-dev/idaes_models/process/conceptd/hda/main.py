from __future__ import division
from pyomo.environ import *
from pyomo.opt import SolverFactory

nlpsolver = SolverFactory('conopt')
mipsolver = SolverFactory('cplex')

m = HDAModel()
