from __future__ import division
from idaes_models.unit.clc.air_rxn_counter import AirRxnCounter
from idaes_models.unit.clc.steam_rxn_counter import SteamRxnCounter
from idaes_models.unit.clc.fuel_rxn_counter import FuelRxnCounter
from idaes_models.process.conceptd.clc.flowsheet import CLCFlowsheet
from idaes_models.core.plugins.propagate_fixed import propagate_var_fix
from six import itervalues
from idaes_models.core.loa import do_LOA
from math import log
from pyomo.environ import value
from pyomo.environ import Constraint

def main():
    m = CLCFlowsheet()
    m.max_flow = 3000
    m.comp_set = ['C', 'CH4', 'CO', 'H2', 'CO2', 'H2O', 'Fe2O3',
                  'Fe3O4', 'FeO', 'Fe', 'O2', 'N2']
    air = m.add_unit(AirRxnCounter(name="air_reactor", parent=m))
    steam = m.add_unit(SteamRxnCounter(name="steam_reactor", parent=m))
    fuel = m.add_unit(FuelRxnCounter(name="fuel_reactor", parent=m))

    #m.connect(air, fuel, from_port='solid_out', to_port='solid_in')
    #m.connect(fuel, steam, from_port='solid_out', to_port='solid_in')
    #m.connect(steam, air, from_port='solid_out', to_port='solid_in')

    m.register_logic_constraint(all_of=(air,fuel))

    m.build_units()
    m.build_logic()
    m.build_links()
    m.expand_connectors()
    propagate_var_fix(m)

    #getattr(fuel.links, 'var_link_with_solid_out.expanded').deactivate()
    #getattr(steam.links, 'var_link_with_solid_out.expanded').deactivate()
    #getattr(air.links, 'var_link_with_solid_out.expanded').deactivate()

    for o in itervalues(m.units):
        o.apply_linear_relaxations()
        o.apply_OA_strategy(oa_ports=False)

    m.solvers.local_NLP.tol = 1E-7
    m.solvers.local_NLP.outlev = 5
    m.solvers.local_NLP.bound_push = 1e-8
    m.solvers.local_NLP.mu_init = 1e-8
    m.solvers.local_NLP.max_iter = 15000
    m.solvers.local_NLP.promote('ipopt')
    # m.solvers.local_NLP.halt_on_ampl_error = 'yes'

    #m._activate_standard_objective()
    for o in itervalues(m.units):
        o.apply_NLP()
    m.dual.activate()

    # Air reactor fixes
    m.units.air_reactor.T.setlb(1100)
    m.units.air_reactor.T.setub(1100)
    m.units.air_reactor.Pt.setlb(100000)
    m.units.air_reactor.Pt.setub(100000)
    '''m.units.air_reactor.F['GasIn'].setlb(1167.085 * 1000 / 3600 )
    m.units.air_reactor.F['GasIn'].setub(1167.085 * 1000 / 3600 )
    m.units.air_reactor.F['GasIn'].fix(1167.085 * 1000 / 3600 )'''
    
    m.units.air_reactor.y_GasIn['CH4'].fix(0)
    m.units.air_reactor.y_GasIn['CO2'].fix(0)
    m.units.air_reactor.y_GasIn['CO'].fix(0)
    m.units.air_reactor.y_GasIn['H2'].fix(1e-5)
    m.units.air_reactor.y_GasIn['H2O'].fix(0.010-1e-5)
    m.units.air_reactor.y_GasIn['O2'].fix(0.209)
    m.units.air_reactor.y_GasIn['N2'].fix(0.781)

    # m.units.air_reactor.F['SolidIn'].setlb(361.3)
    # m.units.air_reactor.F['SolidIn'].setub(361.3)
    # m.units.air_reactor.F['SolidIn'].fix(361.3)
    m.units.air_reactor.y_SolidIn['Fe2O3'].setlb(0.004705231)
    m.units.air_reactor.y_SolidIn['Fe2O3'].setub(0.004705231)
    m.units.air_reactor.y_SolidIn['Fe3O4'].setlb(0.189593136)
    m.units.air_reactor.y_SolidIn['Fe3O4'].setub(0.189593136)
    m.units.air_reactor.y_SolidIn['FeO'].setlb(0.77193468)
    m.units.air_reactor.y_SolidIn['FeO'].setub(0.77193468)
    m.units.air_reactor.y_SolidIn['Fe'].setlb(0.033490174)
    m.units.air_reactor.y_SolidIn['Fe'].setub(0.033490174)
    m.units.air_reactor.y_SolidIn['C'].setlb(1 - 0.004705231 - 0.189593136 - 0.77193468 - 0.033490174)
    m.units.air_reactor.y_SolidIn['C'].setub(1 - 0.004705231 - 0.189593136 - 0.77193468 - 0.033490174)

    # Fuel reactor fixes
    '''m.units.fuel_reactor.T.setlb(1100)
    m.units.fuel_reactor.T.setub(1100)
    m.units.fuel_reactor.Pt.setlb(100000)
    m.units.fuel_reactor.Pt.setub(100000)
    m.units.fuel_reactor.F['GasIn'].setlb(188.732 * 1000 / 3600)
    m.units.fuel_reactor.F['GasIn'].setub(188.732 * 1000 / 3600)
    m.units.fuel_reactor.F['GasIn'].fix(188.732 * 1000 / 3600)
    
    m.units.fuel_reactor.y_GasIn['CH4'].fix(0.973 - 3e-5)
    m.units.fuel_reactor.y_GasIn['CO2'].fix(0.013)
    m.units.fuel_reactor.y_GasIn['CO'].fix(1e-5)
    m.units.fuel_reactor.y_GasIn['H2'].fix(1e-5)
    m.units.fuel_reactor.y_GasIn['H2O'].fix(1e-5)
    m.units.fuel_reactor.y_GasIn['O2'].fix(0)
    m.units.fuel_reactor.y_GasIn['N2'].fix(0.014)'''

    '''m.units.fuel_reactor.y_GasOut['O2'].setlb(0)
    m.units.fuel_reactor.y_GasOut['O2'].setub(0)'''
    '''m.units.fuel_reactor.y_GasOut['O2'].fix(0)'''

    # Steam reactor fixes
    '''m.units.steam_reactor.T.setlb(1100)
    m.units.steam_reactor.T.setub(1100)
    m.units.steam_reactor.Pt.setlb(100000)
    m.units.steam_reactor.Pt.setub(100000)
    m.units.steam_reactor.F['GasIn'].setlb(736.584 * 1000 / 3600)
    m.units.steam_reactor.F['GasIn'].setub(736.584 * 1000 / 3600)
    m.units.steam_reactor.F['GasIn'].fix(736.584 * 1000 / 3600)

    m.units.steam_reactor.y_GasIn['CH4'].fix(0)
    m.units.steam_reactor.y_GasIn['CO2'].fix(0)
    m.units.steam_reactor.y_GasIn['CO'].fix(0)
    m.units.steam_reactor.y_GasIn['H2'].fix(1e-5)
    m.units.steam_reactor.y_GasIn['H2O'].fix(1-1e-5)
    m.units.steam_reactor.y_GasIn['O2'].fix(0)
    m.units.steam_reactor.y_GasIn['N2'].fix(0)
    
    m.units.steam_reactor.y_GasOut['O2'].setlb(0)
    m.units.steam_reactor.y_GasOut['O2'].setub(0)'''


    # Constraint that fixes the total Fe in circulation
    m.units.air_reactor.SolidsCirc1 = Constraint(expr = \
            m.units.air_reactor.F['SolidIn'] *\
            (m.units.air_reactor.y_SolidIn['Fe'] +\
            2*m.units.air_reactor.y_SolidIn['Fe2O3'] +\
            3*m.units.air_reactor.y_SolidIn['Fe3O4'] +\
            m.units.air_reactor.y_SolidIn['FeO']) >= 495)

    m.units.air_reactor.SolidsCirc2 = Constraint(expr = \
            m.units.air_reactor.F['SolidIn'] *\
            (m.units.air_reactor.y_SolidIn['Fe'] +\
            2*m.units.air_reactor.y_SolidIn['Fe2O3'] +\
            3*m.units.air_reactor.y_SolidIn['Fe3O4'] +\
            m.units.air_reactor.y_SolidIn['FeO']) <= 505)

    '''m.units.fuel_reactor.SolidsCirc = Constraint(expr = \
            m.units.fuel_reactor.F['SolidIn'] *\
            (m.units.fuel_reactor.y_SolidIn['Fe'] +\
            2*m.units.fuel_reactor.y_SolidIn['Fe2O3'] +\
            3*m.units.fuel_reactor.y_SolidIn['Fe3O4'] +\
            m.units.fuel_reactor.y_SolidIn['FeO']) == 500)'''

    '''m.units.steam_reactor.SolidsCirc = Constraint(expr = \
            m.units.steam_reactor.F['SolidIn'] *\
            (m.units.steam_reactor.y_SolidIn['Fe'] +\
            2*m.units.steam_reactor.y_SolidIn['Fe2O3'] +\
            3*m.units.steam_reactor.y_SolidIn['Fe3O4'] +\
            m.units.steam_reactor.y_SolidIn['FeO']) == 500)'''

    air.do_init()
    #fuel.do_init()
    #steam.do_init()

    propagate_var_fix(m, tmp=True)
    for o in itervalues(m.units):
        o.introspect_flows()
    for o in itervalues(m.units):
        o.deactivate_trivial_constraints()

    '''m.solve(
            using='local_NLP', tee=True, skip_trivial_constraints=True)

    #m.units.fuel_reactor.y_GasOut['O2'].setlb(0)
    #m.units.fuel_reactor.y_GasOut['O2'].setub(0)
    #m.units.fuel_reactor.y_GasOut['O2'].fix(0)

    # No O2 reaction in fuel reactor
    m.units.fuel_reactor.No_O2_Rxn = Constraint(expr = 0 == \
    m.units.fuel_reactor.y_GasIn['O2']*m.units.fuel_reactor.F['GasIn'] -\
    m.units.fuel_reactor.y_GasOut['O2']*m.units.fuel_reactor.F['GasOut'])'''
    
    m.solve(
        using='local_NLP', tee=True, skip_trivial_constraints=True)                                                   


    print""
    print "Molar Component Flow Rates"
    print""
    
    for k in sorted(m.units.air_reactor.k):
        for j in m.units.air_reactor.comps:
            if value(getattr(m.units.air_reactor, 'y_' + k)[j]):
                print 'Air',k,j,value(m.units.air_reactor.F[k]) * value(getattr(m.units.air_reactor, 'y_' + k)[j]) 
    print""
    
    '''for k in sorted(m.units.fuel_reactor.k):
        for j in m.units.fuel_reactor.comps:
            if value(getattr(m.units.fuel_reactor, 'y_' + k)[j]):
                print 'Fuel',k,j,value(m.units.fuel_reactor.F[k]) * value(getattr(m.units.fuel_reactor, 'y_' + k)[j]) 
    print""
    
    for k in sorted(m.units.steam_reactor.k):
        for j in m.units.steam_reactor.comps:
            if value(getattr(m.units.steam_reactor, 'y_' + k)[j]):
                print 'Steam',k,j,value(m.units.steam_reactor.F[k]) * value(getattr(m.units.steam_reactor, 'y_' + k)[j]) 
    print""'''

    m.units.air_reactor.T.pprint()
    #m.units.fuel_reactor.T.pprint()
    #m.units.steam_reactor.T.pprint()

    m.units.air_reactor.Q.pprint()
    #m.units.fuel_reactor.Q.pprint()
    #m.units.steam_reactor.Q.pprint()

    m.units.air_reactor.F.pprint()
    #m.units.fuel_reactor.F.pprint()
    #m.units.steam_reactor.F.pprint()

    m.units.air_reactor.Pt.pprint()
    #m.units.fuel_reactor.Pt.pprint()
    #m.units.steam_reactor.Pt.pprint()

    for j in m.units.fuel_reactor.comps:
        if j in m.units.fuel_reactor.solid:
            print(j, m.units.fuel_reactor.y_SolidOut[j].value)

    '''for j in m.units.steam_reactor.comps:
        if j in m.units.steam_reactor.solid:
            print(j, m.units.steam_reactor.y_SolidOut[j].value)    '''

    return m


if __name__ == '__main__':
    m = main()
