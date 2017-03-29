"""
A simple test bed for basic unit models


            
"""
from __future__ import division
from __future__ import print_function

__author__ = "Andrew Lee"
__version__ = "1.0.0"

from pyomo.environ import *

import idaes_models.core.flowsheet_model as fs
from idaes_models.properties.physical.PR_ASU import PropPack as prop

class Flowsheet(fs.FlowsheetModel):
    """
    Create the flowsheet class.  Contains all the unit models and 
    connections between them.  Also contains the main Pyomo model and
    solver.
    """
    def __init__(self, *args, **kwargs):
        """
        Create a flowsheet model.
        """
        fs.FlowsheetModel.__init__(self, *args, **kwargs)

        # Create unit models
        self.add_unit(o=prop(name="P1",parent=self,phases='VL',mix_props=1))

def setInputs(fs_obj):
    fs_obj.P1.T.fix(104)
    fs_obj.P1.P.fix(709275)
    fs_obj.P1.z['N2'].fix(0.60)
    fs_obj.P1.z['O2'].fix(0.35)
    fs_obj.P1.z['Ar'].fix(0.05)
    
def print_summary(fs_obj):
    """
    Print some key results from the model.
    """
    print("Results:")
    print("T = {}".format(fs_obj.P1.T.value))
    print("P = {}".format(fs_obj.P1.P.value))
    print()
    print("L = {}".format(fs_obj.P1.L.value))
    print("V = {}".format(fs_obj.P1.V.value))
    print("Z_liq_e = {}".format(fs_obj.P1.Z_liq_e.value))
    print("Z_vap_e = {}".format(fs_obj.P1.Z_vap_e.value))
    print("beta = {}".format(fs_obj.P1.beta.value))
    print("Pe = {}".format(fs_obj.P1.Pe.value))
    for i in fs_obj.P1.comp:
        print("P1 phi_liq[{}] = {}".format(
              i, value(fs_obj.P1.phi_liq[i])))
    for i in fs_obj.P1.comp:
        print("P1 phi_vap[{}] = {}".format(
              i, value(fs_obj.P1.phi_vap[i])))
    
    print()
    for i in fs_obj.P1.comp:
        print("P1 x[{}] = {}".format(
              i, value(fs_obj.P1.x[i])))
    for i in fs_obj.P1.comp:
        print("P1 y[{}] = {}".format(
              i, value(fs_obj.P1.y[i])))

    print()
    print("Pl = {}".format(fs_obj.P1.Pl.value))
    print("Pv = {}".format(fs_obj.P1.Pv.value))
    print("V_liq = {}".format(fs_obj.P1.V_liq.value))
    print("V_vap = {}".format(fs_obj.P1.V_vap.value))
    print("h_liq = {}".format(fs_obj.P1.h_liq.value))
    print("h_vap = {}".format(fs_obj.P1.h_vap.value))
    print("s_liq = {}".format(fs_obj.P1.s_liq.value))
    print("s_vap = {}".format(fs_obj.P1.s_vap.value))
    print("g_liq = {}".format(fs_obj.P1.g_liq.value))
    print("g_vap = {}".format(fs_obj.P1.g_vap.value))
    print()
    print("h_mix = {}".format(fs_obj.P1.h_mix.value))
    print()
    print("mw_liq = {}".format(fs_obj.P1.mw_liq.value))
    print("mw_vap = {}".format(fs_obj.P1.mw_vap.value))
    print()
    print("mw_mix = {}".format(fs_obj.P1.mw_mix.value))
    
def main():
    """
    Make the flowsheet object, fix some variables, and solve the problem
    """
    # Creare flowsheet
    fs_obj = Flowsheet()
    
    # Fix variables
    setInputs(fs_obj)
    
    # Expand connectors
    fs_obj.expand_connectors()
    
    # Initialize properties
    fs_obj.P1._initialize(T=value(fs_obj.P1.T),
                                     P=value(fs_obj.P1.P),
                                     z={'Ar':value(fs_obj.P1.z['Ar']),
                                        'N2':value(fs_obj.P1.z['N2']),
                                        'O2':value(fs_obj.P1.z['O2'])},
                                    outlvl=1)

    # Solve flowhseet model
    fs_obj.solve(tee=True,options={'tol':1e-10})
    
    # Print results to:
    #  1) check right number of variables and equations
    #  2) make sure the problem solved successfully.
    fs_obj.display_results()
    print_summary(fs_obj)

if __name__ == "__main__":
    main()
