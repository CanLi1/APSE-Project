"""
A simple test bed for basic unit models


            
"""
from __future__ import division
from __future__ import print_function

__author__ = "Andrew Lee"
__version__ = "1.0.0"

from pyomo.environ import *

import idaes_models.core.flowsheet_model as fs
import BFB as BFB
#import BFB_210217 as BFB

class Flowsheet(fs.FlowsheetModel):
    def __init__(self, *args, **kwargs):
        """
        Create a flowsheet model.
        """
        fs.FlowsheetModel.__init__(self, *args, **kwargs)
        
        # Create set of discretisation points
        nfe = 45
        fe_a = 0.51
        fe_b = 0.22
        fe_set = {0:0.0}
        for i in range(1,nfe+1):
            if i < nfe*fe_a:
                fe_set[float(i)] = i*fe_b/(nfe*fe_a)
            else:
                fe_set[float(i)] = fe_b + (i-nfe*fe_a)*(1-fe_b)/(nfe*(1-fe_a))

        # Create unit models
        self.add_unit(o=BFB.BFB(name = 'BFB',
                                parent = self,
                                fe_set = fe_set,
                                s_inlet = "Bottom",
                                s_outlet = "Overflow",
                                gas_prop_lib = 'gas_prop_test',
                                sol_prop_lib = 'solid_prop_test',
                                hx_prop_lib = 'hx_prop_test',
                                vb_method = 'Davidson'))

def setInputs(fs_obj):
    # Change Parameters
    fs_obj.BFB.Kd = 1
    
    #Fix some variables
    fs_obj.BFB.Dt.fix(15)
    fs_obj.BFB.Lb.fix(1.88)
    fs_obj.BFB.nor.fix(2500)
    fs_obj.BFB.dx.fix(0.1)
    fs_obj.BFB.lhx.fix(0.15)

    fs_obj.BFB.Gas_In_F.fix(4337.1)   # mol/s
    fs_obj.BFB.Gas_In_P.fix(122200)   # Pa
    fs_obj.BFB.Gas_In_T.fix(338.89)        # K
    fs_obj.BFB.Gas_In_y['CO2'].fix(0.1285)
    fs_obj.BFB.Gas_In_y['H2O'].fix(0.0785)
    fs_obj.BFB.Gas_In_y['N2'].fix(0.793)

    fs_obj.BFB.Solid_In_F.fix(370.37) # kg/s
    fs_obj.BFB.Solid_In_T.fix(355.56)      # K
    fs_obj.BFB.Solid_In_x['Car'].fix(1.27865)
    fs_obj.BFB.Solid_In_x['H2O'].fix(0.213835)
    
    fs_obj.BFB.HX_In_F.fix(27800)    # mol/s
    fs_obj.BFB.HX_In_T.fix(305.37)         # s
    fs_obj.BFB.HX_In_P.fix(112000)     # Pa
    fs_obj.BFB.HX_In_y['H2O'].fix(1)

def print_summary(fs_obj):  
    
    fs_obj.BFB.Gas_Out_F.display()
    fs_obj.BFB.Gas_Out_T.display()
    fs_obj.BFB.Gas_Out_P.display()
    fs_obj.BFB.Gas_Out_y.display()
    fs_obj.BFB.Solid_Out_F.display()
    fs_obj.BFB.Solid_Out_T.display()
    fs_obj.BFB.Solid_Out_x.display()
    fs_obj.BFB.HX_Out_T.display()
    fs_obj.BFB.HX_Out_P.display()
    
    removal = {}
    mbal_tol = {}
    mbal_tol['Sorb'] = 1 - fs_obj.BFB.Solid_Out_F.value \
                            / fs_obj.BFB.Solid_In_F.value
    for j in fs_obj.BFB.GasList:
        if j == 'CO2':
            removal[j] = 1 - value(fs_obj.BFB.Gas_Out_F) \
                                * value(fs_obj.BFB.Gas_Out_y['CO2']) \
                            / (value(fs_obj.BFB.Gas_In_F) \
                                * value(fs_obj.BFB.Gas_In_y['CO2']))
            mbal_tol[j] = 1 - (value(fs_obj.BFB.Gas_Out_F) \
                                * value(fs_obj.BFB.Gas_Out_y['CO2']) \
                            + value(fs_obj.BFB.Solid_Out_F.value) \
                                * value(fs_obj.BFB.Solid_Out_x['Car'])) \
                            / (value(fs_obj.BFB.Gas_In_F) \
                                    * value(fs_obj.BFB.Gas_In_y['CO2']) \
                                + value(fs_obj.BFB.Solid_In_F) \
                                    * value(fs_obj.BFB.Solid_In_x['Car']))
        elif j == 'H2O':
            removal[j] = 1 - value(fs_obj.BFB.Gas_Out_F) \
                                * value(fs_obj.BFB.Gas_Out_y['H2O']) \
                            / (value(fs_obj.BFB.Gas_In_F) \
                                * value(fs_obj.BFB.Gas_In_y['H2O']))
            mbal_tol[j] = 1 - (fs_obj.BFB.Gas_Out_F.value \
                                * fs_obj.BFB.Gas_Out_y['H2O'].value \
                            + fs_obj.BFB.Solid_Out_F.value \
                                * fs_obj.BFB.Solid_Out_x['H2O'].value) \
                            / (fs_obj.BFB.Gas_In_F.value \
                                    * fs_obj.BFB.Gas_In_y['H2O'].value \
                                + fs_obj.BFB.Solid_In_F.value \
                                    * fs_obj.BFB.Solid_In_x['H2O'].value)
        else:
            mbal_tol[j] = (fs_obj.BFB.Gas_In_F.value \
                                * fs_obj.BFB.Gas_In_y['N2'].value \
                            - fs_obj.BFB.Gas_Out_F.value \
                                * fs_obj.BFB.Gas_Out_y['N2'].value) \
                            / (fs_obj.BFB.Gas_In_F.value \
                                * fs_obj.BFB.Gas_In_y['N2'].value)
    ebal_gas = value(fs_obj.BFB.Gas_In_F) \
                    * value(fs_obj.BFB.gas_prop_b[0].h_vap) \
                - value(fs_obj.BFB.Gb[fs_obj.BFB.nfe]) \
                    * value(fs_obj.BFB.gas_prop_b[fs_obj.BFB.nfe].h_vap) \
                - value(fs_obj.BFB.Ge[fs_obj.BFB.nfe]) \
                    * value(fs_obj.BFB.gas_prop_e[fs_obj.BFB.nfe].h_vap)
    # Modify next equation for boundary conditions
    # Overflow conditions
    ebal_sol = value(fs_obj.BFB.Solid_In_F)*value(fs_obj.BFB.sol_prop_f.h_sol) \
                - value(fs_obj.BFB.Solid_Out_F) \
                    * value(fs_obj.BFB.sol_prop_e[fs_obj.BFB.feom[fs_obj.BFB.nfe]].h_sol)
    # Underflow conditions
    '''ebal_sol = value(fs_obj.BFB.Solid_In_F)*value(fs_obj.BFB.sol_prop_f.h_sol) \
                - value(fs_obj.BFB.Solid_Out_F) \
                    * value(fs_obj.BFB.sol_prop_e[1].h_sol)'''
    ebal_tol = ebal_gas + ebal_sol + value(fs_obj.BFB.Q)

    print('Removal:',removal)
    print('Mass Balance Tolerance:',mbal_tol)
    print('Energy Balance Tolerance:','%.5f'%(ebal_tol),
                          "(Gas:",'%.1f'%(ebal_gas),
                          "Solid:",'%.1f'%(ebal_sol),
                          "HX:",'%.1f'%(value(fs_obj.BFB.Q)),")")
    
    #fs_obj.pprint()

def main():
    """
    Make the flowsheet object, fix some variables, and solve the problem
    """
    # Creare flowsheet
    fs_obj = Flowsheet()
    
    # Fix variables
    setInputs(fs_obj)
    
    # Expand connectors
    #fs_obj.expand_connectors()
    
    # Initialize model
    fs_obj.BFB._initialize(outlvl=1,
                           optarg={"tol"            : 1e-8,
                                   "max_cpu_time"   : 300,
                                   "print_level"    : 5})#,
                                   #"halt_on_ampl_error": 'yes'})
    '''fs_obj.load_json(fname="file_name.json")'''

    # Solve flowhseet model
    sopts = {"tol"            : 1e-8,
             "max_cpu_time"   : 300,
             "print_level"    : 5}

    fs_obj.solve(tee=False,options=sopts)

    # Print results to:
        #  1) check right number of variables and equations
        #  2) make sure the problem solved successfully.
    fs_obj.display_results()
    print_summary(fs_obj)
    
    '''for i in fs_obj.BFB.l:
        if i != 0:
            print(i,", r_c:",value(fs_obj.BFB.sol_prop_c[i].r[2]),", r_e:",value(fs_obj.BFB.sol_prop_e[i].r[2]))
    for i in fs_obj.BFB.l:
        if i != 0:
            print(i,", rs_c:",value(fs_obj.BFB.rsc['Car',i]),", rs_e:",value(fs_obj.BFB.rse['Car',i]))'''

    '''fs_obj.save_json(fname="file_name.json")'''
    
if __name__ == "__main__":
    main()

'''Gas_Out_F : Gas Flowrate at Gas Outlet [mol/s]
    Size=1, Index=None
    Key  : Lower : Value         : Upper : Fixed : Stale : Domain
    None :  None : 4092.63889789 :  None : False : False :  Reals
Gas_Out_T : Gas Temeprature at Gas Outlet [K]
    Size=1, Index=None
    Key  : Lower : Value         : Upper : Fixed : Stale : Domain
    None :  None : 346.915000928 :  None : False : False :  Reals
Gas_Out_P : Pressure at Gas Outlet [Pa]
    Size=1, Index=None
    Key  : Lower : Value         : Upper : Fixed : Stale : Domain
    None :  None : 116814.445594 :  None : False : False :  Reals
Gas_Out_y : Gas Mole Fractions at Gas Outlet [mol/mol]
    Size=3, Index=BFB.Gas_Out_y_index
    Key : Lower : Value           : Upper : Fixed : Stale : Domain
    CO2 :  None :  0.092209342401 :  None : False : False :  Reals
    H2O :  None : 0.0674232599843 :  None : False : False :  Reals
     N2 :  None :  0.840367397615 :  None : False : False :  Reals
Solid_Out_F : Solid Flowrate at Solid Outlet [kg/s]
    Size=1, Index=None
    Key  : Lower : Value  : Upper : Fixed : Stale : Domain
    None :  None : 370.37 :  None : False : False :  Reals
Solid_Out_T : Solid Temperature at Solid Outlet [K]
    Size=1, Index=None
    Key  : Lower : Value         : Upper : Fixed : Stale : Domain
    None :  None : 346.176907353 :  None : False : False :  Reals
Solid_Out_x : Solid Active Mass at Solid Outlet [kg/kg]
    Size=2, Index=BFB.Solid_Out_x_index
    Key : Lower : Value          : Upper : Fixed : Stale : Domain
    Car :  None :  1.76448256889 :  None : False : False :  Reals
    H2O :  None : 0.388048066839 :  None : False : False :  Reals
HX_Out_T : HX Fluid Temperature at Outlet [K]
    Size=1, Index=None
    Key  : Lower : Value         : Upper : Fixed : Stale : Domain
    None :  None : 318.515142335 :  None : False : False :  Reals
HX_Out_P : HX Fluid Pressure at Outlet [Pa]
    Size=1, Index=None
    Key  : Lower : Value         : Upper : Fixed : Stale : Domain
    None :  None : 130183.482239 :  None : False : False :  Reals
Removal: {'H2O': 0.18951667802714833, 'CO2': 0.3228641788031533}
Mass Balance Tolerance: {'H2O': -2.220446049250313e-16, 'CO2': 0.0, 'N2': -2.6444024471141236e-16, 'Sorb': 0.0}
Energy Balance Tolerance: -0.00000 (Gas: -688609.6 Solid: 28193144.7 HX: -27504535.1 )'''
