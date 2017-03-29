"""
Model constructor and holder for the Hydrodealkylation of toluene process.
This will draw inspiration from John Eslick's process_base and flowsheet_model classes as well as my own (Qi's) model_flow classes in the methanol and water modules. In the future, this model will be used as a basis for developing a unified process flowsheet base object.
"""
from __future__ import division
from pyomo.environ import *

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

class HDAModel():
	def __init__(self, *args, **kwargs):

		raise NotImplementedError('This needs to be updated. It is woefully out of date')

		m.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
		m.alpha = Param(initialize=0.3665)  # compressor coefficient
		m.compeff = Param(initialize=0.750)  # compressor efficiency
		m.gamma = Param(initialize=1.300)  # ratio of cp to cv
		m.abseff = Param(initialize=0.333)  # absorber efficiency
		m.disteff = Param(initialize=0.5)  # column tray efficiency
		m.flow_ub = Param(initialize=50)  # flow upper bound
		m.P_ub = Param(initialize=4.0)  # pressure upper bound
		m.T_ub = Param(initialize=7.0)  # temperature upper bound

		m.cost_elec = Param(initialize=0.340)  # electricity cost ($340. per kw-yr)
		m.cost_Qc = Param(initialize=0.700)  # cooling cost ($700 per 1E9 kj)
		m.cost_Qh = Param(initialize=8.000)  # heating cost ($8000 per 1E9 kj)
		m.cost_fuel = Param(initialize=4.0)  # fuel cost furnace ($4 per 1e6 btu)

		# chemical components: hydrogen, methane, benzene, toluene, diphenyl
		m.comp_set = Set(initialize=['H2', 'CH4', 'BEN', 'TOL', 'DIP'])
