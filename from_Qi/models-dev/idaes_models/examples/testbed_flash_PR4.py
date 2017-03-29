"""
A simple test bed for basic unit models


            
"""
from __future__ import division
from __future__ import print_function

__author__ = "Anthony Burgard"
__version__ = "1.0.0"

from pyomo.environ import *

import idaes_models.core.unit_model as unit
import idaes_models.core.flowsheet_model as fs
from idaes_models.unit.standard.flash_PR4 import flash
#from properties import *

class Flowsheet(fs.FlowsheetModel):
    """
    Create the flowsheet class.  Contains all the unit models and 
    connections between them.  Also contains the main Pyomo model and
    solver.
    """
    def __init__(self, comp, name="Flowsheet",solver="ipopt"):
        """
        Create a flowsheet model.
        """
        fs.FlowsheetModel.__init__(self, name, solver)
        self.comp = comp
        model = self.model
        #External root finding functions
        plib = '../properties/physical/cubic_eos/phys_prop.so'
        self.model.ceos_z_liq = ExternalFunction(
            library=plib, function="ceos_z_liq")
        self.model.ceos_z_vap = ExternalFunction(
            library=plib, function="ceos_z_vap")
        self.add_unit("Fl1",flash(comp,Nout=2,main_model=model))
        
        
        '''# Create property calls
        self._import_prop_calls(model)'''
        # Create unit models
        self.add_unit("Fl1",flash(comp,Nout=2,main_model=model))
                     
    '''def _import_prop_calls(self, model, p_lib="phys_prop.so"):
        """
        Load external functions for property calls. These will be used
        by unit and other sub models.
        """
        model.p_h_vap = ExternalFunction(library=p_lib, function="h_vap")'''

def setInputs(fs_obj):

    fs_obj.model.Fl1.In_F.fix(1)
    fs_obj.model.Fl1.In_T.fix(95)
    fs_obj.model.Fl1.In_P.fix(7*101325)
    fs_obj.model.Fl1.In_y['N2'].fix(0.60)
    fs_obj.model.Fl1.In_y['O2'].fix(0.35)
    fs_obj.model.Fl1.In_y['Ar'].fix(0.05)
    
    for j in fs_obj.comp:
        fs_obj.model.Fl1.Out_y[1,j] = 1/3
        fs_obj.model.Fl1.Out_y[2,j] = 1/3
        
    fs_obj.model.Fl1.T.fix(104.34) 
    fs_obj.model.Fl1.P.fix(101325*7)
    
 
    
def print_summary(fs_obj):
    """
    Print some key results from the model.
    """
    print("---")
    print("Flash F_in = {}".format(fs_obj.model.Fl1.In_F.value))
    print("Flash T_in = {}".format(fs_obj.model.Fl1.In_T.value))
    print("Flash P_in = {}".format(fs_obj.model.Fl1.In_P.value))  
    for j in fs_obj.comp:
        print("Flash y[{}]_out = {}".format(
          j,fs_obj.model.Fl1.In_y[j].value))            
    print("---")

    # Outlet streams for flash
    Nout = [1,2]
   
    print("Flash Q = {}".format(fs_obj.model.Fl1.Q.value))
    for i in Nout:
        print("Flash F[{}]_Out = {}".format(i,fs_obj.model.Fl1.Out_F[i].value))
        print("Flash T[{}]_Out = {}".format(i,fs_obj.model.Fl1.Out_T[i].value))
        print("Flash P[{}]_Out = {}".format(i,fs_obj.model.Fl1.Out_P[i].value))  
        for j in fs_obj.comp:
            print("Flash y[{},{}]_out = {}".format(
              i,j,fs_obj.model.Fl1.Out_y[i,j].value))            
        print("---")
   
def main():
    """
    Make the flowsheet object, fix some variables, and solve the problem
    """
    # Component List
    comp = ['N2','O2','Ar']
    fs_obj = Flowsheet(comp)

    #Fix some variables
    setInputs(fs_obj)

    #Expand connectors (consider comining with sovle in flowsheet?)
    fs_obj.expand_connectors()

    #Solve flowhseet model
    opt = SolverFactory('ipopt') 
    results = opt.solve(fs_obj.model,tee=True,keepfiles=False,options={\
            'outlev':5,'bound_push':1e-8,'mu_init':1e-8}) 
    
    #fs_obj.solve(tee=True,keepfiles=False)
    # skip_trivial_constraints=True

    # Print results to:
    #  1) check right number of variables and equations
    #  2) make sure the problem solved successfully.

    print_summary(fs_obj)
 
    fs_obj.model.display()
    
    #Try out load and save from json
    #fs_obj.save_json("test.json")
    #fs_obj.model.Heater.Q.value = 2.0
    #fs_obj.load_json("test.json")
    #print_summary(fs_obj)
    
    
if __name__ == "__main__":
    main()
