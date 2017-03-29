"""
A flowsheet for BFB
"""
from __future__ import division
from __future__ import print_function

__author__ = "Andrew Lee"  
__version__ = "1.0.1"
# stochastics added by David Woodruff, Feb 2017

from pyomo.environ import *

import idaes_models.core.flowsheet_model as fs
import BFB as BFB

# Flags to control what we are debugging.
# The original (pre-stochastic) behavior is matched when one uses True True
# As of Feb 28 stochastic is not converging for the EF for target = False; i.e. False False fails

usetarget = False # set to False to use cost as the objective function.
deterministic =True # to False to use stochastic
print ("usetarget="+str(usetarget), \
       "deterministic="+str(deterministic))

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
                fe_set[i] = i*fe_b/(nfe*fe_a)
            else:
                fe_set[i] = fe_b + (i-nfe*fe_a)*(1-fe_b)/(nfe*(1-fe_a))

        # Create unit models
        self.add_unit(o=BFB.BFB(name = 'BFB',
                                parent = self,
                                fe_set = fe_set,
                                s_inlet = "Bottom",
                                s_outlet = "Overflow",
                                gas_prop_lib = 'gas_prop_test',
                                sol_prop_lib = 'solid_prop_test',
                                hx_prop_lib = 'hx_prop_test'))

#==========
def create_objective_fct(fs_obj):
    # Create the parameters and functions for an objective function.
    # Input: fs_obj: the flowsheet
    #        usetarge: if true, use a target, otherwise use cost

    if usetarget:
        # switched to param to make it easier to do quick tests using default
        fs_obj.Lb = Param(default = 4.0, mutable=True)
        fs_obj.Dt = Param(default = 10.0, mutable=True)
        fs_obj.dx = Param(default = 0.2, mutable=True)
        fs_obj.lhx = Param(default = 0.25, mutable=True)

        def _e(fs_obj):
            return (fs_obj.BFB.Lb - fs_obj.Lb)**2 \
                     + (fs_obj.BFB.Dt - fs_obj.Dt)**2 \
                     + (fs_obj.BFB.dx - fs_obj.dx)**2 \
                     + (fs_obj.BFB.lhx - fs_obj.lhx)**2
        fs_obj.obj = Objective(rule = _e)
        if not deterministic:
            # sillyness just for debugging
            fs_obj.FirstStageCost = Expression(rule = _e)
            fs_obj.SecondStageCost = Expression(rule = _e)
            

    else:
        # For stochastic, we are going to use a cost based objective.
        # Do so in a way that supports stochastic programming.
        # But this also works for a deterministic program.
        # NOTE: some stochastic programming programmers prefer the
        #       the stage costs in a list, particularly for multi-stage.

        # first stage cost Params
        fs_obj.SurfaceCostRate = Param(default = 1.0)
        fs_obj.PerTubeCost = Param(default = 1.0)
        # second stage cost Params
        fs_obj.SolidRateCost = Param(default = 1.0)  # factor in the NPV coeff
        fs_obj.CoolantRateCost = Param(default = 0.1) # factor in the NPV coeff

        def _e1(fs_obj):
            # approximate cost
            return fs_obj.SurfaceCostRate * 3.14 * fs_obj.BFB.Dt * fs_obj.BFB.Lb \
                + fs_obj.PerTubeCost * 3.14 * fs_obj.BFB.Dt / fs_obj.BFB.lhx
        fs_obj.FirstStageCost = Expression(rule = _e1)

        def _e2(fs_obj):
            # approximate cost
            return fs_obj.SolidRateCost * fs_obj.BFB.Solid_In_F \
                + fs_obj.SolidRateCost * fs_obj.BFB.HX_In_F 
        fs_obj.SecondStageCost = Expression(rule = _e2)

        fs_obj.CostObjective = Objective(expr = \
                               fs_obj.FirstStageCost + fs_obj.SecondStageCost)

    
#==========
   
def create_bfb_instance():
    """ 
    Create the flowsheet object for a given scenario_name.
    In general, the flowsheets should be implmented to have one of these.
    """
    # Create flowsheet
    fs_obj = Flowsheet()

    # Fix variables (on the BFB block; some will be unfixed below)
    setInputs(fs_obj)

    # Expand connectors
    #fs_obj.expand_connectors()

    # Feb 28, 2017 hack by DLW to make BFB happy (and if you undo it after init using numeric=False, the nl writer won't be happy)
    # The next line effectively "undoes" the last few lines of setInputs
    do_Suffix_driven_fixing(fs_obj, numeric = True)  

    # Initialize model
    fs_obj.BFB._initialize(outlvl=1,
                           optarg={"tol"            : 1e-8,
                                   "max_cpu_time"   : 300,
                                   "print_level"    : 5})#,
                                   #"halt_on_ampl_error": 'yes'})
    '''fs_obj.load_json(fname="file_name.json")'''

    # Solve flowhseet model
    print('Begin Optimisation')

    # Create an Objective function
    create_objective_fct(fs_obj)
    
    # Add bounds to state variables
    for i in fs_obj.BFB.l:
        fs_obj.BFB.vg[i].setlb(10*fs_obj.BFB.sol_prop_e[1].v_mf)

        fs_obj.BFB.P[i].setub(150000)
        fs_obj.BFB.P[i].setlb(100000)
        
        fs_obj.BFB.Tgb[i].setub(380)
        fs_obj.BFB.Tgb[i].setlb(300)
        fs_obj.BFB.Tgc[i].setub(380)
        fs_obj.BFB.Tgc[i].setlb(300)
        fs_obj.BFB.Tge[i].setub(380)
        fs_obj.BFB.Tge[i].setlb(300)
        fs_obj.BFB.Tsc[i].setub(380)
        fs_obj.BFB.Tsc[i].setlb(300)
        fs_obj.BFB.Tse[i].setub(380)
        fs_obj.BFB.Tse[i].setlb(300)

        for j in fs_obj.BFB.GasList:
            if j == 'N2':
                fs_obj.BFB.yb[j,i].setub(1.0)
                fs_obj.BFB.yb[j,i].setlb(0.0)
                fs_obj.BFB.yc[j,i].setub(1.0)
                fs_obj.BFB.yc[j,i].setlb(0.0)
                fs_obj.BFB.ye[j,i].setub(1.0)
                fs_obj.BFB.ye[j,i].setlb(0.0)
            else:
                fs_obj.BFB.yb[j,i].setub(0.2)
                fs_obj.BFB.yb[j,i].setlb(0.0)
                fs_obj.BFB.yc[j,i].setub(0.2)
                fs_obj.BFB.yc[j,i].setlb(0.0)
                fs_obj.BFB.ye[j,i].setub(0.2)
                fs_obj.BFB.ye[j,i].setlb(0.0)

        for j in fs_obj.BFB.SolidList:
            fs_obj.BFB.xc[j,i].setub(2.0)
            fs_obj.BFB.xc[j,i].setlb(0.0)
            fs_obj.BFB.xe[j,i].setub(2.0)
            fs_obj.BFB.xe[j,i].setlb(0.0)

    # Adding bounds to manipulated variables
    fs_obj.BFB.Lb.setub(5.0)
    fs_obj.BFB.Lb.setlb(0.5)
    
    fs_obj.BFB.Dt.setub(15.0)
    fs_obj.BFB.Dt.setlb(1.0)

    fs_obj.BFB.dx.setub(0.2)
    fs_obj.BFB.dx.setlb(0.02)
    
    fs_obj.BFB.lhx.setub(1.0)
    fs_obj.BFB.lhx.setlb(0.1)

    fs_obj.BFB.Gas_In_F.setub(6500.0)
    fs_obj.BFB.Gas_In_F.setlb(2500.0)

    fs_obj.BFB.Gas_In_T.setub(340.0)
    fs_obj.BFB.Gas_In_T.setlb(300.0)

    for j in fs_obj.BFB.GasList:
        if j == 'N2':
            fs_obj.BFB.Gas_In_y[j].setub(1.0)
            fs_obj.BFB.Gas_In_y[j].setlb(0.0)
        elif j == 'CO2':
            fs_obj.BFB.Gas_In_y[j].setub(0.2)
            fs_obj.BFB.Gas_In_y[j].setlb(0.08)
        else:
            fs_obj.BFB.Gas_In_y[j].setub(0.15)
            fs_obj.BFB.Gas_In_y[j].setlb(0.03)

    # Unfix variables
    fs_obj.BFB.Lb.unfix()
    fs_obj.BFB.Dt.unfix()
    fs_obj.BFB.dx.unfix()
    fs_obj.BFB.lhx.unfix()

    return fs_obj

#===========================
# begin stoch
import os
import json
def pysp_instance_creation_callback(scenario_name, node_names):
    """
    This callback is used by PySP to instantiate scenario instances.
    """
    print ("Debug: In callback for scenario=",str(scenario_name))
    fs_obj = create_bfb_instance()
    dir = "stoch_dir" ## st.PySP_dir()
    fname = os.path.join(dir, scenario_name + ".json")
    print ("Debug: in scenario callback; fname=",fname)
    with open(fname, "r") as scenfile:
        ScenIn = json.load(scenfile)

    # DLW: as of Feb 2017 the following only works for singleton Params.
    #    When you fix this, I hope you delete both of these comment lines!
    for i in ScenIn:
        if not hasattr(fs_obj.BFB, i):
            raise RuntimeError("Found "+i+" in "+fname+" but not the model.")
        p = getattr(fs_obj.BFB, i)
        p.store_values(ScenIn[i])

    return fs_obj
# end stoch
#========================================================

#=========
def setup_Var_suffixes(fs_obj):
    # WORK IN PROGRESS
    # NOTE as of Feb 28, 2017 this gets called but then "undone" and the literals are used to fix Vars, really.
    """ Put Param suffixes on some Vars so the values will appear
    in Pyomo expressions as Params that can easily be changed, e.g.,
    for stochastic programming, where the code does not want to
    be hard-wired for fixing Vars, but rather changes any Param
    it is fed (see pysp_instance_creation_callback)
    """
    
    # DLW: As of Feb 2017 this is the "by-hand" version. Cryptic note:
    #      The number of lines could be cut almost in half by using gettr
    #      and function with a string argument to set up the suffix.

    fs_obj.fixed_values = Suffix(direction=Suffix.LOCAL)

    # Declare Params for the values… and associate with Suffixes
    fs_obj.BFB.Gas_In_F__fix = Param(initialize=4337.1, mutable=True)
    fs_obj.fixed_values[fs_obj.BFB.Gas_In_F] = fs_obj.BFB.Gas_In_F__fix
    fs_obj.BFB.Gas_In_P__fix = Param(initialize=122200.0, mutable=True)
    fs_obj.fixed_values[fs_obj.BFB.Gas_In_P] = fs_obj.BFB.Gas_In_P__fix
    fs_obj.BFB.Gas_In_T__fix = Param(initialize=338.89, mutable=True)
    fs_obj.fixed_values[fs_obj.BFB.Gas_In_T] = fs_obj.BFB.Gas_In_T__fix
    # NOTE to scenario data creators: the Params for y are singletons.
    fs_obj.BFB.Gas_In_y_CO2__fix = Param(initialize=0.1285, mutable=True)
    fs_obj.fixed_values[fs_obj.BFB.Gas_In_y['CO2']] = \
                        fs_obj.BFB.Gas_In_y_CO2__fix
    fs_obj.BFB.Gas_In_y_H2O__fix = Param(initialize=0.0785, mutable=True)
    fs_obj.fixed_values[fs_obj.BFB.Gas_In_y['H2O']] = \
                        fs_obj.BFB.Gas_In_y_H2O__fix
    fs_obj.BFB.Gas_In_y_N2__fix = Param(initialize=0.793, mutable=True)
    fs_obj.fixed_values[fs_obj.BFB.Gas_In_y['N2']] = \
                        fs_obj.BFB.Gas_In_y_N2__fix

    
def do_Suffix_driven_fixing(fs_obj, numeric = False):
    # For Vars that have a fixed_value suffix: fix at the suffix.
    # numeric = True is a hack to make other code happy by using the value.
    import six
    for comp, parm in six.iteritems(fs_obj.fixed_values):
        if not numeric:
            comp.fix(parm)
        else:
            comp.fix(value(parm))

#==========
        
def setInputs(fs_obj):
    # Change Parameters
    fs_obj.BFB.Kd = 1
    
    #Fix some variables
    fs_obj.BFB.Dt.fix(15)
    fs_obj.BFB.Lb.fix(1.88)
    fs_obj.BFB.nor.fix(2500)
    fs_obj.BFB.dx.fix(0.1)
    fs_obj.BFB.lhx.fix(0.15)

    """ A suffix is used for these to help stochastics:
    fs_obj.BFB.Gas_In_F.fix(4337.1)   # mol/s
    fs_obj.BFB.Gas_In_P.fix(122200)   # Pa
    fs_obj.BFB.Gas_In_T.fix(338.89)        # K
    fs_obj.BFB.Gas_In_y['CO2'].fix(0.1285)
    fs_obj.BFB.Gas_In_y['H2O'].fix(0.0785)
    fs_obj.BFB.Gas_In_y['N2'].fix(0.793)
    """

    fs_obj.BFB.Solid_In_F.fix(370.37) # kg/s
    fs_obj.BFB.Solid_In_T.fix(385.56)      # K
    fs_obj.BFB.Solid_In_x['Car'].fix(1.27865)
    fs_obj.BFB.Solid_In_x['H2O'].fix(0.213835)
    
    fs_obj.BFB.HX_In_F.fix(27800)    # mol/s
    fs_obj.BFB.HX_In_T.fix(305.37)         # s
    fs_obj.BFB.HX_In_P.fix(112000)     # Pa
    fs_obj.BFB.HX_In_y['H2O'].fix(1)

    # One can imagine setInputs as just the Param Fix and these 2 calls.
    setup_Var_suffixes(fs_obj)
    do_Suffix_driven_fixing(fs_obj)  

def print_summary(fs_obj):  
    
    fs_obj.BFB.Gas_Out_F.display()
    fs_obj.BFB.Gas_Out_T.display()
    fs_obj.BFB.Gas_Out_P.display()
    fs_obj.BFB.Gas_Out_y.display()
    fs_obj.BFB.Solid_Out_F.display()
    fs_obj.BFB.Solid_Out_T.display()
    fs_obj.BFB.Solid_Out_x.display()
    
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
                    * value(fs_obj.BFB.sol_prop_e[fs_obj.BFB.nfe].h_sol)
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

    #fs_obj.BFB.pprint()
    fs_obj.BFB.Lb.display()
    fs_obj.BFB.Dt.display()
    fs_obj.BFB.dx.display()
    fs_obj.BFB.lhx.display()

def main():

    sopts = {"tol"            : 1e-8,
             "max_cpu_time"   : 300,
             "print_level"    : 5}

    if deterministic:
        # Make the flowsheet object, fix some variables, and solve the problem.

        # the next line is the "original"
        fs_obj = create_bfb_instance()

        # the next line is for debugging scenario data
        # fs_obj =  pysp_instance_creation_callback("Scenario_2", None)

        if usetarget:
            Lb = [4.0,2.8]
            Dt = [10.0,6.0]
            dx = [0.2,0.05]
            lhx = [0.25,0.12]
            runcount = len(Lb)
        else:
            runcount = 1

        for i in range(runcount):

            if usetarget:
                fs_obj.Lb = Lb[i]
                fs_obj.Dt = Dt[i]
                fs_obj.dx = dx[i]
                fs_obj.lhx = lhx[i]

            fs_obj.solve(tee=False, options=sopts)

            # Print results to:
                #  1) check right number of variables and equations
                #  2) make sure the problem solved successfully.
            fs_obj.display_results()
            #print_summary(fs_obj)
            print("Run "+str(i+1)+":",value(fs_obj.BFB.Lb),value(fs_obj.BFB.Dt),value(fs_obj.BFB.dx),value(fs_obj.BFB.lhx))
            '''fs_obj.save_json(fname="file_name.json")'''

    else:  # must be stochastic
        # one target is enough, so use defaults if using target
        import sys
        import pyomo.pysp.daps.basicclasses as bc
        import pyomo.pysp.daps.stoch_solver as st

        # stoch_dir needs to have the scenario tree node files with daps names
        # and needs the scenario tree template file 'bfb_stoch_tree_input.json'
        tree_model = bc.Tree_2Stage_json_dir('stoch_dir',
                     'bfb_stoch_tree_input.json')
        
        fsfilename = __file__
        solver = st.StochSolver(fsfilename, tree_model)
        ef_instance = solver.solve_ef('ipopt', sopts)
        print ("We have a stochastic solution... Now do something with it.")
        # ef_instance.pprint() # not a great thing to do with it, but something.
        
if __name__ == "__main__":
    main()
