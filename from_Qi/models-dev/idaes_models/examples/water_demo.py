"""
Small demonstration example for water treatment network problem. This is a flat
model useful primarily for training purposes, getting acquainted with Pyomo.
"""
from __future__ import division
from pyomo.environ import ConcreteModel, Set, Constraint, NonNegativeReals, Var, Binary, minimize, Objective, Block
from pyomo.opt import SolverFactory
from idaes_models.core.util.mccormick import add_mccormick_relaxation

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"

def main():
	m = ConcreteModel('water treatment')
	comp = m.Comps = Set(initialize=['A', 'B', 'C'])
	tru = m.treatment_units = Set(initialize=['unit1', 'unit2'])
	outlet = m.outlet = Set(initialize=['out'])
	inlet = m.inlet = Set(initialize=['in'])
	in_flow = {'A': 10, 'B': 20, 'C': 100}
	m.in_split_frac = Var(tru | outlet, domain=NonNegativeReals, bounds=(0, 1))
	m.flow = Var(inlet | tru, outlet | tru, comp, domain=NonNegativeReals, bounds=(0, 500))
	m.flow_into_tru = Var(tru, comp, domain=NonNegativeReals, bounds=(0, 500))
	m.flow_from_tru = Var(tru, comp, domain=NonNegativeReals, bounds=(0, 500))
	m.tru_split_frac = Var(tru, tru | outlet, domain=NonNegativeReals, bounds=(0, 1))
	m.flow_out = Var(comp, domain=NonNegativeReals, bounds=(0, 500))
	m.unit_exists = Var(tru, domain=Binary)

	def inlet_mass_bal(m, sink, comp):
		return in_flow[comp] * m.in_split_frac[sink] == m.flow['in', sink, comp]
	m.inlet_mass_bal = Constraint(tru | outlet, comp, rule=inlet_mass_bal)
	# Sum of split fractions are equal to 1
	m.inlet_frac_sum = Constraint(expr=sum(m.in_split_frac[sink] for sink in tru | outlet) == 1)

	def tru_mix_bal(m, t, comp):
		return sum(m.flow[source, t, comp] for source in inlet | tru) == m.flow_into_tru[t, comp]
	m.tru_mix_bal = Constraint(tru, comp, rule=tru_mix_bal)

	perf = {
		'unit1': {'A': 0.2, 'B': 0.35, 'C': 0},
		'unit2': {'A': 0.3, 'B': 0.2, 'C': 0}}
	def tru_performance(m, tru, comp):
		return m.flow_from_tru[tru, comp] == m.flow_into_tru[tru, comp] * (1 - perf[tru][comp])
	m.tru_performance = Constraint(tru, comp, rule=tru_performance)

	def tru_split_bal(m, tru, sink, comp):
		return m.flow_from_tru[tru, comp] * m.tru_split_frac[tru, sink] == m.flow[tru, sink, comp]
	m.tru_split_bal = Constraint(tru, tru | outlet, comp, rule=tru_split_bal)

	def tru_split_sum(m, t):
		return sum(m.tru_split_frac[t, sink] for sink in tru | outlet) == 1
	m.tru_split_sum = Constraint(tru, rule=tru_split_sum)

	def outlet_mix_bal(m, comp):
		return m.flow_out[comp] == sum(m.flow[source, 'out', comp] for source in inlet | tru)
	m.outlet_mix_bal = Constraint(comp, rule=outlet_mix_bal)

	outlet_limit = {'A': 3, 'B': 10}
	def outlet_mix_limit(m, comp):
		return m.flow_out[comp] <= outlet_limit[comp] if comp in outlet_limit else Constraint.NoConstraint
	m.outlet_mix_limit = Constraint(comp, rule=outlet_mix_limit)

	def unit_exist_flow_limit(m, tru):
		return m.flow_into_tru[tru, 'C'] <= m.flow_into_tru[tru, 'C'].ub * m.unit_exists[tru]
	m.unit_exist_flow_limit = Constraint(tru, rule=unit_exist_flow_limit)

	m.objective = Objective(expr=sum(m.flow_into_tru[t, 'C'] for t in tru), sense=minimize)

	globalsolver = SolverFactory('baron')
	globalsolver.options['CplexLibName'] = "/opt/ibm/ILOG/CPLEX_Studio1263/cplex/bin/x86-64_linux/libcplex1263.so"
	globalsolver.options['MaxTime'] = 180
	globalsolver.options['allowipopt'] = 0
	globalsolver.options['EpsA'] = 0.01
	globalsolver.options['EpsR'] = 0.0001
	results = globalsolver.solve(m, tee=True)

	print(results)
	m.display()

	# Block for holding linearized equations
	m.lin = Block()
	nsegs = 1  # number of McCormick segments
	for (source, sink, c) in tru * (tru | outlet) * comp:
		add_mccormick_relaxation(m.lin, m.flow[source, sink, c], m.flow_from_tru[source, c], m.tru_split_frac[source, sink], nsegs, (source, sink, c), 1.0)

	m.tru_split_bal.deactivate()

	mipsolver = SolverFactory('gurobi')
	results = mipsolver.solve(m, tee=True)

	print(results)
	# m.display()


if __name__ == '__main__':
	main()
