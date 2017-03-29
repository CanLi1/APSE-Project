from __future__ import division
from __future__ import print_function
from pyomo.environ import *
from pyomo.opt import SolverFactory

import csv

def resi(a,b,c,z):
    return z**3 + a*z**2 + b*z + c

data = []
with open('testing1.txt', 'rb') as data_file:
    rdr = csv.reader(data_file, delimiter=' ')
    for row in rdr:
        row = map(float, row)
        data.append(row)

model = ConcreteModel()
opt = SolverFactory('ipopt')

a = model.a = Var(domain=Reals, initialize=-0.996637378369)
b = model.b = Var(domain=Reals, initialize=1)
c = model.c = Var(domain=Reals, initialize=1)
zl = model.zl = Var(domain=Reals, initialize=1)
zv = model.zv = Var(domain=Reals, initialize=1)

model.f_zl = ExternalFunction(
    library="./phys_prop.so", function="ceos_z_liq")

model.f_zv = ExternalFunction(
    library="./phys_prop.so", function="ceos_z_vap")

model.eq_zl = Constraint(expr=zl == model.f_zl(a,b,c))
model.eq_zv = Constraint(expr=zv == model.f_zv(a,b,c))

for i,r in enumerate(data):
    a.fix(data[i][2])
    b.fix(data[i][3])
    c.fix(data[i][4])

    res = opt.solve(model)
    print("z_liq = {:+0.5e}, z_vap = {:+0.5e}, res_liq = {:+0.5e}, res_vap = {:+0.5e}".format(
        zl.value - data[i][5], 
        zv.value - data[i][6],
        resi(a.value, b.value, c.value, zl.value),
        resi(a.value, b.value, c.value, zv.value)
        ))
        

