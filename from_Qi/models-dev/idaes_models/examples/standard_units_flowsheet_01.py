"""
A simple test bed for basic unit models


            
"""
from __future__ import division
from __future__ import print_function

__author__ = "Andrew Lee"
__version__ = "1.0.0"

from pyomo.environ import *

import idaes_models.core.unit_model as unit
import idaes_models.core.flowsheet_model as fs
from idaes_models.unit.standard.heater import Heater
from idaes_models.unit.standard.mixer import Mixer
import os

class Flowsheet(fs.FlowsheetModel):
    """
    Create the flowsheet class.  Contains all the unit models and 
    connections between them.  Also contains the main Pyomo model and
    solver.
    """
    def __init__(self, comp, name="Flowsheet", solver="ipopt"):
        """
        Create a flowsheet model.
        """
        p_lib = "../properties/physical/MEA_Simple/phys_prop.so"
        fs.FlowsheetModel.__init__(
            self, 
            name, 
            solver,
            prop_funcs={
                "h_vap":("h_vap",p_lib)})
        self.comp = comp
        
        model = self.model
        # Create unit models
        self.add_unit("Heater", 
                      Heater(comp, prop_dict=
                        {"h_vap":model.phys_prop.h_vap}))
        self.add_unit("Mixer",  
                      Mixer(comp, NI=2, prop_dict=
                        {"h_vap":model.phys_prop.h_vap}))
        # Create streams
        self.connect("S3", model.Mixer.Outlet, model.Heater.Inlet)
                     
def setInputs(fs_obj):
    fs_obj.model.Mixer.In_F[1].fix(100)
    fs_obj.model.Mixer.In_T[1].fix(298.15)
    fs_obj.model.Mixer.In_P[1].fix(101325)
    fs_obj.model.Mixer.In_y[1, 'CO2'].fix(0)
    fs_obj.model.Mixer.In_y[1, 'H2O'].fix(0)
    fs_obj.model.Mixer.In_y[1, 'N2'].fix(1)
    fs_obj.model.Mixer.In_y[1, 'MEA'].fix(0)
    fs_obj.model.Mixer.In_y[1, 'O2'].fix(0)
    
    fs_obj.model.Mixer.In_F[2].fix(10)
    fs_obj.model.Mixer.In_T[2].fix(330)
    fs_obj.model.Mixer.In_P[2].fix(2e5)
    fs_obj.model.Mixer.In_y[2, 'CO2'].fix(1)
    fs_obj.model.Mixer.In_y[2, 'H2O'].fix(0)
    fs_obj.model.Mixer.In_y[2, 'N2'].fix(0)
    fs_obj.model.Mixer.In_y[2, 'MEA'].fix(0)
    fs_obj.model.Mixer.In_y[2, 'O2'].fix(0)
    
    fs_obj.model.Heater.deltaT.fix(10)
    
def print_summary(fs_obj):
    """
    Print some key results from the model.
    """
    print("Results:")
    print("---")
    print("---")
    print("Mixer F_Out = {}".format(fs_obj.model.Mixer.Out_F.value))
    print("Mixer T_Out = {}".format(fs_obj.model.Mixer.Out_T.value))
    print("Mixer P_Out = {}".format(fs_obj.model.Mixer.Out_P.value))
    for i in fs_obj.comp:
        print("Mixer y[{}]_out = {}".format(
              i, fs_obj.model.Mixer.Out_y[i].value))
    print("---")
    print("Heater F_Out = {}".format(fs_obj.model.Heater.Out_F.value))
    print("Heater T_Out = {}".format(fs_obj.model.Heater.Out_T.value))
    print("Heater P_Out = {}".format(fs_obj.model.Heater.Out_P.value))
    for i in fs_obj.comp:
        print("Heater y[{}]_out = {}".format(
              i, fs_obj.model.Heater.In_y[i].value))
    print("---")
    print("Heater deltaT = {}".format(fs_obj.model.Heater.deltaT.value))
    print("Heater Q = {}".format(fs_obj.model.Heater.Q.value))
    
def main():
    """
    Make the flowsheet object, fix some variables, and solve the problem
    """
    # Component List
    comp = ['CO2', 'H2O', 'N2', 'MEA', 'O2']
    fs_obj = Flowsheet(comp)
    #Fix some variables
    setInputs(fs_obj)
    #Expand connectors (consider comining with sovle in flowsheet?)
    fs_obj.expand_connectors()
    #Solve flowhseet model
    fs_obj.solve()
    
    # Print results to:
    #  1) check right number of variables and equations
    #  2) make sure the problem solved successfully.
    fs_obj.display_results()
    print_summary(fs_obj)
    
    #Try out load and save from json
    #fs_obj.save_json("test.json")
    #fs_obj.model.Heater.Q.value = 2.0
    #fs_obj.load_json("test.json")
    #print_summary(fs_obj)
    
    
if __name__ == "__main__":
    main()
