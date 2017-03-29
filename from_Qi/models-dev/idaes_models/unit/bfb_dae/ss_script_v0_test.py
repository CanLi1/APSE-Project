from ref_ss_2 import bfb_ss
from pyomo.opt import SolverFactory
from pyomo.core import Var
"""
written by David M Thierry
based on ss_script_v0.py
test bfb_ss write matlab file
"""
__author__ = 'David M Thierry'
solver = SolverFactory('ipopt')
solver.options['linear_solver'] = 'ma57'
solver.options['max_iter'] = 100
solver.options['print_level'] = 5


nfe_x = 5
kord_x = 3

mod = bfb_ss(nfe_x, kord_x, 0, 0)

instance = mod.create_instance('dat_ss_1.dat')
results = solver.solve(instance, tee=True, logfile='loghere.log', keepfiles=True)
instance.solutions.load_from(results)


file_dat = open('dmatlab.m', 'w')
file_dat.write('% matlab file written \n')
print('writing matlab file..')
spl = ['h', 'c', 'n']
dummy1 = 0
for v in instance.component_objects(Var, active=True):
    vobj = getattr(instance, str(v))
    file_dat.write('\n\n')
    # if dummy1 > 0:

    name = str(v)

    for i in vobj._index:
        if i == None:
            k = 0
        else:
            k = len(i)
    if k == 0:
        continue
    if k == 1:
        print(k, str(v))
    elif k == 2:
        file_dat.write(name + ' = [' + '\n')
        for ix in range(1, nfe_x + 1):
            for jx in range(1, kord_x + 1):
                idx = (ix, jx)
                valx = vobj[idx].value
                stidx = ''
                for elm in idx:
                    stidx += str(elm) + ' '
                # stidx = str(ix) + ' ' + str(jx) + '\t'
                file_dat.write('\t' + str(valx) + '\n')
                dummy1 += 1
        file_dat.write('];\n')
    elif k == 3:
        for c in spl:
            file_dat.write(name + '_' + c + ' = [' + '\n')
            for it in range(1, nfe_x + 1):
                for jt in range(1, kord_x + 1):
                    idx = (it, jt, c)
                    valx = vobj[idx].value
                    stidx = ''
                    for elm in idx:
                        stidx += str(elm) + ' '
                    # stidx = str(ix) + ' ' + str(jx) + '\t'
                    file_dat.write('\t' + str(valx) + '\n')
                    dummy1 += 1
            file_dat.write('];\n')
import time
file_dat.write('% written ' + str(dummy1) + ' vars')
file_dat.write("% Current time " + time.strftime("%X") + '\n\n\n')

file_dat.write('lx' + ' = [' + '\n')
for ix in range(1, nfe_x + 1):
    for jx in range(1, kord_x + 1):
        idx = (ix, jx)
        valx = instance.l[idx]
        stidx = ''
        # stidx = str(ix) + ' ' + str(jx) + '\t'
        file_dat.write('\t' + str(valx) + '\n')
file_dat.write('];\n')

file_dat.write('plot(lx, Gb)\n')

file_dat.close()