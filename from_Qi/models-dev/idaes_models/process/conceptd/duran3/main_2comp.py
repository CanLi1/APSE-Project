from __future__ import division
from pyomo.environ import *
from pyomo.environ import exp as pyo_exp
from casadi import exp as cas_exp
from casadi import SX, SXFunction
from pyomo.opt import SolverFactory
from math import copysign
from collections import Iterable, OrderedDict

"""
Created on Fri 2016-03-11 21:48

@author: qtothec

Re-implementation of Duran example 3 superstructure synthesis (the "eight process problem") in Pyomo with outer approximation cuts generated automatically by automatic differentiation. This variation is meant to be a test bed for handling of indexed variables.
In this test, we assume that all flows are composed of a mixture of component A and component B.

Ref:
	SELECT OPTIMAL PROCESS FROM WITHIN GIVEN SUPERSTRUCTURE.
	MARCO DURAN , PH.D. THESIS (EX3) , 1984.
	CARNEGIE-MELLON UNIVERSITY , PITTSBURGH , PA.

(original problem, my implementation may vary)
	Problem type:	convex MINLP
			size:	 8  binary variables
					26  continuous variables
					32  constraints

Pictoral representation can be found on
page 969 of Turkay & Grossmann, 1996.
http://dx.doi.org/10.1016/0098-1354(95)00219-7

Used result of set covering algorithm from Turkay paper cited above for initial NLP selection
"""

"""
Solver Prep
"""
nlpsolver = SolverFactory('minos')
# nlpsolver.options["halt_on_ampl_error"] = "yes"
# nlpsolver.options["max_iter"] = 5000

mipsolver = SolverFactory('cplex')

baron = SolverFactory('baron')

"""
	Declare Pyomo model
"""
m = ConcreteModel(name='DuranEx3')

"""
	Set declarations
"""
I = m.I = Set(initialize=xrange(2, 25+1))		# process streams
J = m.J = Set(initialize=xrange(1, 8+1))		# process units
C = m.C = Set(initialize=['A', 'B'])			# chemical components
PI = m.PI = Set(initialize=xrange(1, 4+1))		# set of purely integer constraints
DS = m.DS = Set(initialize=xrange(1, 4+1))		# set of design specifications
"""		1: Unit 8
		2: Unit 8
		3: Unit 4
		4: Unit 4
"""
MB = m.MB = Set(initialize=xrange(1, 7+1))		# set of mass balances
"""		1: 4-6-7
		2: 3-5-8
		3: 4-5
		4: 1-2
		5: 1-2-3
		6: 6-7-4
		7: 6-7
"""

"""
	Parameter and initial point declarations
"""
# FIXED COST INVESTMENT COEFF FOR PROCESS UNITS
# Format: process #: cost
fixed_cost = {1: 5, 2: 8, 3: 6, 4: 10, 5: 6, 6: 7, 7: 4, 8: 5}
CF = m.CF = Param(J, initialize=fixed_cost)

# VARIABLE COST COEFF FOR PROCESS UNITS - STREAMS
# Format: stream #: cost
variable_cost = {3: -10, 5: -15, 9: -40, 19: 25, 21: 35, 25: -35, 17: 80, 14: 15, 10: 15, 2: 1, 4: 1, 18: -65, 20: -60, 22: -80}
CV = m.CV = Param(I, initialize=variable_cost, default=0)

# initial point information for equipment selection (for each NLP subproblem)
initY = {
	'sub1': {1: 1, 2: 0, 3: 1, 4: 1, 5: 0, 6: 0, 7: 1, 8: 1},
	'sub2': {1: 0, 2: 1, 3: 1, 4: 1, 5: 0, 6: 1, 7: 0, 8: 1},
	'sub3': {1: 1, 2: 0, 3: 1, 4: 0, 5: 1, 6: 0, 7: 0, 8: 1}
	}
# initial point information for stream flows
initX = {2: 2, 3: 1.5, 6: 0.75, 7: 0.5, 8: 0.5, 9: 0.75, 11: 1.5, 12: 1.34, 13: 2, 14: 2.5, 17: 2, 18: 0.75, 19: 2, 20: 1.5, 23: 1.7, 24: 1.5, 25: 0.5}
initXc = {k: v / len(C) for k, v in initX.iteritems()}

"""
	Variable declarations
"""
# BINARY VARIABLE DENOTING EXISTENCE-NONEXISTENCE
Y = m.Y = Var(J, domain=Binary, initialize=initY['sub1'])
# FLOWRATES OF PROCESS STREAMS
X = m.X = Var(I, domain=NonNegativeReals, initialize=initX)
Xc = m.Xc = Var(I, C, domain=NonNegativeReals, initialize=initXc)
# OBJECTIVE FUNCTION CONSTANT TERM
CONSTANT = m.constant = Param(initialize=122.0)

"""
	Set up data structures and functions for automatic differentiation
"""
nl = m.nonlinear_eqns = Block()
jacs = {}  # Dictionary data structure to hold jacobian information for nonlinear constraints

# Helper functions to generate equality constraints
def q_exp(x):
	# Apply a different exponential function depending on whether a Pyomo or CasADi object is being passed
	if isinstance(x, SX):
		return cas_exp(x)
	else:
		return pyo_exp(x)

def create_nonlinear_eq_constr(name, m, pyo_vars, eq_expr, *sets):
	def eq_rule(b, *sets):
		return eq_expr(*sets, **pyo_vars) == 0
	new_constr = Constraint(*sets, rule=eq_rule)
	setattr(m.nonlinear_eqns, name, new_constr)
	cas_vars = {}
	flat_vars = {}
	for var_name, pyo_var in pyo_vars.iteritems():
		if pyo_var.is_indexed():
			# if the pyomo variable is indexed, generate a dictionary of CasADi variables
			cas_vars[var_name] = {indx: SX.sym(var_name + str(indx)) for indx in pyo_var._index}
			flat_vars.update({var_name + str(indx): (cas_vars[var_name][indx], pyo_var[indx]) for indx in pyo_var._index})
		else:
			cas_vars[var_name] = SX.sym(var_name)
			flat_vars[var_name] = (cas_vars[var_name], pyo_var)
	cas_var_name_list, var_list = zip(*flat_vars.iteritems())
	cas_var_list, pyo_var_list = zip(*var_list)
	if len(sets) == 0:
		f = eq_expr(**cas_vars)
		opt = {'input_scheme': cas_var_name_list, 'output_scheme': ['f']}
		eqn = SXFunction(name, cas_var_list, [f], opt)
		jacs[new_constr] = OrderedDict([(k, {'jac': eqn.jacobian(k, 'f'), 'var': flat_vars[k][1]}) for k in cas_var_name_list])
	else:
		indices = reduce(lambda x, y: x * y, sets)  # Cartesian product of sets
		for i in indices:
			# for each index, generate function. If statement below checks to see if i is a tuple or just a single index and acts appropriately.
			f = eq_expr(*i, **cas_vars) if isinstance(i, Iterable) else eq_expr(i, **cas_vars)
			opt = {'input_scheme': cas_var_name_list, 'output_scheme': ['f']}
			eqn = SXFunction(name, cas_var_list, [f], opt)
			jacs[new_constr[i]] = OrderedDict([(k, {'jac': eqn.jacobian(k, 'f'), 'var': flat_vars[k][1]}) for k in cas_var_name_list])

def get_jacobian(eqn, var_name):
	return float(jacs[eqn][var_name]['jac'](map(lambda x: x['var'].value, jacs[eqn].values()))[0])

"""
	Constraint definitions
"""

# INPUT-OUTPUT RELATIONS FOR process units 1 through 8
# commented lines indicate what the traditional Pyomo code would be

# nl.inout1 = Constraint(expr=exp(X[3]) - 1 == X[2])
def rule_inout1(c, Xc):
	return q_exp(Xc[3, c]) - 1 - Xc[2, c]
create_nonlinear_eq_constr('inout1', m, {'Xc': Xc}, rule_inout1, C)
# def rule_inout1(i, X):
# 	return q_exp(X[3]) - 1 - X[2]
# create_nonlinear_eq_constr('inout1', m, {'X': X}, rule_inout1, I)
# def rule_inout1(x3, x2):
# 	return q_exp(x3) - 1 - x2
# create_nonlinear_eq_constr('inout1', m, {'x2': X[2], 'x3': X[3]}, rule_inout1)
# for c in C:
# 	create_nonlinear_eq_constr('inout1' + c, m, {'x2': Xc[2, c], 'x3': Xc[3, c]}, rule_inout1)

# nl.inout2 = Constraint(expr=exp(X[5] / 1.2) - 1 == X[4])
def rule_inout2(x5, x4):
	return q_exp(x5 / 1.2) - 1 - x4
create_nonlinear_eq_constr('inout2', m, {'x5': X[5], 'x4': X[4]}, rule_inout2)

def rule_inout3(b, c): return 1.5 * Xc[9, c] + Xc[10, c] == Xc[8, c]
m.inout3 = Constraint(C, rule=rule_inout3)

def rule_inout4(b, c): return 1.25 * (Xc[12, c] + Xc[14, c]) == Xc[13, c]
m.inout4 = Constraint(C, rule=rule_inout4)

def rule_inout5(b, c): return Xc[15, c] == 2 * Xc[16, c]
m.inout5 = Constraint(C, rule=rule_inout5)

# nl.inout6 = Constraint(expr=exp(X[20] / 1.5) - 1 == X[19])
def rule_inout6(x20, x19):
	return q_exp(x20 / 1.5) - 1 - x19
create_nonlinear_eq_constr('inout6', m, {'x20': X[20], 'x19': X[19]}, rule_inout6)

# nl.inout7 = Constraint(expr=exp(X[22]) - 1 == X[21])
def rule_inout7(x22, x21):
	return q_exp(x22) - 1 - x21
create_nonlinear_eq_constr('inout7', m, {'x22': X[22], 'x21': X[21]}, rule_inout7)

# nl.inout8 = Constraint(expr=exp(X[18]) - 1 == X[10] + X[17])
def rule_inout8(x18, x10, x17):
	return q_exp(x18) - 1 - x10 - x17
create_nonlinear_eq_constr('inout8', m, {'x18': X[18], 'x10': X[10], 'x17': X[17]}, rule_inout8)

# Mass balance equations
def rule_massbal1(b, c): return Xc[13, c] == Xc[19, c] + Xc[21, c]
m.massbal1 = Constraint(C, rule=rule_massbal1)
def rule_massbal2(b, c): return Xc[17, c] == Xc[9, c] + Xc[16, c] + Xc[25, c]
m.massbal2 = Constraint(C, rule=rule_massbal2)
def rule_massbal3(b, c): return Xc[11, c] == Xc[12, c] + Xc[15, c]
m.massbal3 = Constraint(C, rule=rule_massbal3)
def rule_massbal4(b, c): return Xc[3, c] + Xc[5, c] == Xc[6, c] + Xc[11, c]
m.massbal4 = Constraint(C, rule=rule_massbal4)
def rule_massbal5(b, c): return Xc[6, c] == Xc[7, c] + Xc[8, c]
m.massbal5 = Constraint(C, rule=rule_massbal5)
def rule_massbal6(b, c): return Xc[23, c] == Xc[20, c] + Xc[22, c]
m.massbal6 = Constraint(C, rule=rule_massbal6)
def rule_massbal7(b, c): return Xc[23, c] == Xc[14, c] + Xc[24, c]
m.massbal7 = Constraint(C, rule=rule_massbal7)
def rule_compbal(b, i): return X[i] == sum(Xc[i, c] for c in C)
m.compbal = Constraint(I, rule=rule_compbal)

# process specifications
def rule_specs1(b, c): return Xc[10, c] <= 0.8 * Xc[17, c]
m.specs1 = Constraint(C, rule=rule_specs1)
def rule_specs2(b, c): return Xc[10, c] >= 0.4 * Xc[17, c]
m.specs2 = Constraint(C, rule=rule_specs2)
def rule_specs3(b, c): return Xc[12, c] <= 5 * Xc[14, c]
m.specs3 = Constraint(C, rule=rule_specs3)
def rule_specs4(b, c): return Xc[12, c] >= 2 * Xc[14, c]
m.specs4 = Constraint(C, rule=rule_specs4)

# Logical constraints (big-M) for each process
# These allow for iff unit j exists
m.logical1 = Constraint(expr=X[2] <= 10 * Y[1])
m.logical2 = Constraint(expr=X[4] <= 10 * Y[2])
m.logical3 = Constraint(expr=X[9] <= 10 * Y[3])
m.logical4 = Constraint(expr=X[12] + X[14] <= 10 * Y[4])
m.logical5 = Constraint(expr=X[15] <= 10 * Y[5])
m.logical6 = Constraint(expr=X[19] <= 10 * Y[6])
m.logical7 = Constraint(expr=X[21] <= 10 * Y[7])
m.logical8 = Constraint(expr=X[10] + X[17] <= 10 * Y[8])

# pure integer constraints
m.pureint1 = Constraint(expr=Y[1] + Y[2] == 1)
m.pureint2 = Constraint(expr=Y[4] + Y[5] <= 1)
m.pureint3 = Constraint(expr=Y[6] + Y[7] - Y[4] == 0)
m.pureint4 = Constraint(expr=Y[3] - Y[8] <= 0)

"""
	Profit (objective) function definition
"""
m.profit = Objective(expr=(-1) * sum(Y[j] * CF[j] for j in J) - sum(X[i] * CV[i] for i in I) - CONSTANT, sense=maximize)

"""
	Bound definitions
"""
# x (flow) upper bounds
x_ubs = {3: 2, 5: 2, 9: 2, 10: 1, 14: 1, 17: 2, 19: 2, 21: 2, 25: 3}
for i, x_ub in x_ubs.iteritems():
	X[i].setub(x_ub)
	for c in C:
		Xc[i, c].setub(x_ub)

"""
	Other problem setup
"""
# fix integer variables
Y.fix()

# Set up duals information reporting
duals = m.dual = Suffix(direction=Suffix.IMPORT_EXPORT)

"""
	Set up results reporting
"""
lower_bound = 0
upper_bound = float('inf')
oa_iter = 1
oa_master_iter = 0
iter_results = []

def record_result(prob_type):
	new_result = {
		'Mstr': oa_master_iter if prob_type == 'master' else '',
		'Sub': oa_iter if prob_type == 'sub' else '',
		'obj': round(value(m.profit.expr)),
		'UB': round(upper_bound) if prob_type == 'master' else '',
		'LB': round(lower_bound) if prob_type == 'sub' else '',
		'status': str(results.solver.termination_condition)
	}
	iter_results.append(new_result)
# Ordered list of field names for the CSV output headers
fieldnames = ['Mstr', 'Sub', 'obj', 'LB', 'UB', 'status']

"""
	First NLP subproblem
"""
# Solve first NLP subproblem
results = nlpsolver.solve(m)

# record results
lower_bound = value(m.profit.expr)
record_result('sub')
print('NLP Subproblem ' + str(oa_iter) + ': ' + str(round(value(m.profit.expr), 1)) + ' L: ' + str(round(lower_bound, 1)) + ' U: ' + str(round(upper_bound, 1)))

"""
	Set up OA data structures and cut generation function
"""
lin = m.linear_cuts = Block()  # model data structure to hold generated linear constraints and related additional parameters and variables
lin.deactivate()  # keep this deactivated for now, since we are still solving NLPs
oa_cuts = lin.oa_cuts = ConstraintList()  # list of outer approximation cuts
int_cuts = lin.integer_cuts = ConstraintList()  # list of integer cuts

def add_oa_cuts():
	# For each entry in the dictionary of nonlinear functions that we have calculated jacobians for, add an outer approximation cut.
	for constr, jac_entry in jacs.iteritems():
		oa_cuts.add(expr=copysign(1, duals[constr]) * sum(get_jacobian(constr, var_name) * (var_entry['var'] - var_entry['var'].value) for var_name, var_entry in jac_entry.iteritems()) <= 0)

# adds OA and integer cuts for the first subproblem
add_oa_cuts()
int_cuts.add((1, sum(1 - Y[j] for j in J if round(Y[j].value, 1) == 1) + sum(Y[j] for j in J if round(Y[j].value, 1) == 0), None))

"""
	Second NLP subproblem
"""
# Fix binaries for second NLP subproblem
for j in J:
	Y[j].fix(initY['sub2'][j])

# Solve second NLP subproblem
results = nlpsolver.solve(m)
oa_iter += 1

# record results and add OA and integer cuts
lower_bound = max([lower_bound, value(m.profit.expr)])
record_result('sub')
add_oa_cuts()
int_cuts.add((1, sum(1 - Y[j] for j in J if round(Y[j].value, 1) == 1) + sum(Y[j] for j in J if round(Y[j].value, 1) == 0), None))
print('NLP Subproblem ' + str(oa_iter) + ': ' + str(round(value(m.profit.expr), 1)) + ' L: ' + str(round(lower_bound, 1)) + ' U: ' + str(round(upper_bound, 1)))

"""
	Third NLP subproblem
"""
# Fix binaries for third NLP subproblem
for j in J:
	Y[j].fix(initY['sub3'][j])

# Solve third NLP subproblem
results = nlpsolver.solve(m)
oa_iter += 1

# record results and add OA and integer cuts
lower_bound = max([lower_bound, value(m.profit.expr)])
record_result('sub')
add_oa_cuts()
int_cuts.add((1, sum(1 - Y[j] for j in J if round(Y[j].value, 1) == 1) + sum(Y[j] for j in J if round(Y[j].value, 1) == 0), None))
print('NLP Subproblem ' + str(oa_iter) + ': ' + str(round(value(m.profit.expr), 1)) + ' L: ' + str(round(lower_bound, 1)) + ' U: ' + str(round(upper_bound, 1)))

"""
	Begin loop to iterate between master problems and subproblems
"""
while lower_bound < upper_bound and oa_iter < 10:
	# MIP master problem
	Y.unfix()  # unfix binaries
	nl.deactivate()  # deactivate nonlinear equations
	lin.activate()  # activate linear cuts
	results = mipsolver.solve(m, warmstart=True)  # solve MIP
	oa_master_iter += 1  # increment master iteration counter
	upper_bound = min([upper_bound, value(m.profit.expr)])  # update UB
	record_result('master')  # record iteration result
	print('MIP Master ' + str(oa_master_iter) + ': ' + str(round(value(m.profit.expr), 1)) + ' L: ' + str(round(lower_bound, 1)) + ' U: ' + str(round(upper_bound, 1)))  # print results to screen

	# Terminate loop if bounds have closed
	if lower_bound >= upper_bound:
		break

	# NLP subproblem
	Y.fix()  # fix binary variables
	nl.activate()  # activate nonlinear equations
	lin.deactivate()  # deactivate linear equations
	results = nlpsolver.solve(m)  # solve NLP
	oa_iter += 1  # increment iteration counter
	lower_bound = max([lower_bound, value(m.profit.expr)])  # update LB
	record_result('sub')  # record iteration result
	add_oa_cuts()  # adds OA cut
	int_cuts.add((1, sum(1 - Y[j] for j in J if round(Y[j].value, 1) == 1) + sum(Y[j] for j in J if round(Y[j].value, 1) == 0), None))  # add integer cut
	print('NLP Subproblem ' + str(oa_iter) + ': ' + str(round(value(m.profit.expr), 1)) + ' L: ' + str(round(lower_bound, 1)) + ' U: ' + str(round(upper_bound, 1)))  # print results to screen

"""
	Write results to a CSV file
"""
import csv

with open('iter_results.csv', 'wb') as output_file:
	dict_writer = csv.DictWriter(output_file, fieldnames=fieldnames)
	dict_writer.writeheader()
	dict_writer.writerows(iter_results)
