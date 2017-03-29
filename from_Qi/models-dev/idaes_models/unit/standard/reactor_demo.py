"""
A simple test bed for basic unit models


            
"""
from __future__ import division
from __future__ import print_function

__author__ = "Andrew Lee"
__version__ = "1.0.0"

from pyomo.environ import *

import idaes_models.core.flowsheet_model as fs
from idaes_models.unit.standard.rstoic import UnitModel as rstoic
from idaes_models.unit.standard.rcstr import UnitModel as rcstr
from idaes_models.unit.standard.requil import UnitModel as requil
from idaes_models.unit.standard.rpfr import UnitModel as rpfr
from idaes_models.unit.standard.rgibbs import UnitModel as rgibbs

class Flowsheet(fs.FlowsheetModel):
    """
    Create the flowsheet class.  Contains all the unit models and 
    connections between them.  Also contains the main Pyomo model and
    solver.
    """
    def __init__(self, name="Flowsheet", solver="ipopt"):
        """
        Create a flowsheet model.
        """
        fs.FlowsheetModel.__init__(self, name, solver)
        model = self.model

        # Create unit models
        self.add_unit("R1", rstoic(main_model=model))
        self.add_unit("R2", rcstr(main_model=model))
        self.add_unit("R3", requil(main_model=model))
        self.add_unit("R4", rpfr(main_model=model))
        self.add_unit("R5", rgibbs(main_model=model))

def setInputs(fs_obj):
    # Stoichiometric reactor
    fs_obj.model.R1.In_F.fix(100)
    fs_obj.model.R1.In_T.fix(298.15)
    fs_obj.model.R1.In_P.fix(101325)
    fs_obj.model.R1.In_y['a'].fix(0.5)
    fs_obj.model.R1.In_y['b'].fix(0.5)
    fs_obj.model.R1.In_y['c'].fix(0)
    fs_obj.model.R1.In_y['d'].fix(0)
    fs_obj.model.R1.In_y['e'].fix(0)
    fs_obj.model.R1.In_y['f'].fix(0)

    fs_obj.model.R1.Out_P.fix(101325)
    
    fs_obj.model.R1.x_rxn[1].fix(20)
    fs_obj.model.R1.x_rxn[2].fix(5)
    fs_obj.model.R1.x_rxn[3].fix(5)

    fs_obj.model.R1.Q.fix(-1500000)
    
    # CSTR
    fs_obj.model.R2.V.fix(100)

    fs_obj.model.R2.In_F.fix(100)
    fs_obj.model.R2.In_T.fix(298.15)
    fs_obj.model.R2.In_P.fix(101325)
    fs_obj.model.R2.In_y['a'].fix(0.5)
    fs_obj.model.R2.In_y['b'].fix(0.5)
    fs_obj.model.R2.In_y['c'].fix(0)
    fs_obj.model.R2.In_y['d'].fix(0)
    fs_obj.model.R2.In_y['e'].fix(0)
    fs_obj.model.R2.In_y['f'].fix(0)

    fs_obj.model.R2.Out_P.fix(101325)

    fs_obj.model.R2.Q.fix(-1500000)

    # K Equilibrium Reactor
    fs_obj.model.R3.In_F.fix(100)
    fs_obj.model.R3.In_T.fix(298.15)
    fs_obj.model.R3.In_P.fix(101325)
    fs_obj.model.R3.In_y['a'].fix(0.5)
    fs_obj.model.R3.In_y['b'].fix(0.5)
    fs_obj.model.R3.In_y['c'].fix(0)
    fs_obj.model.R3.In_y['d'].fix(0)
    fs_obj.model.R3.In_y['e'].fix(0)
    fs_obj.model.R3.In_y['f'].fix(0)

    fs_obj.model.R3.Out_P.fix(101325)

    fs_obj.model.R3.Q.fix(-1500000)

    # PFR
    fs_obj.model.R4.Lr.fix(2)
    fs_obj.model.R4.A.fix(0.5)    

    fs_obj.model.R4.In_F.fix(100)
    fs_obj.model.R4.In_T.fix(298.15)
    fs_obj.model.R4.In_P.fix(101325)
    fs_obj.model.R4.In_y['a'].fix(0.5)
    fs_obj.model.R4.In_y['b'].fix(0.5)
    fs_obj.model.R4.In_y['c'].fix(0)
    fs_obj.model.R4.In_y['d'].fix(0)
    fs_obj.model.R4.In_y['e'].fix(0)
    fs_obj.model.R4.In_y['f'].fix(0)
    
    # Gibbs reactor
    fs_obj.model.R5.In_F.fix(100)
    fs_obj.model.R5.In_T.fix(298.15)
    fs_obj.model.R5.In_P.fix(101325)
    fs_obj.model.R5.In_y['H2'].fix(0.1)
    fs_obj.model.R5.In_y['N2'].fix(0.474)
    fs_obj.model.R5.In_y['O2'].fix(0.126)
    fs_obj.model.R5.In_y['CO2'].fix(0)
    fs_obj.model.R5.In_y['CH4'].fix(0.3)
    fs_obj.model.R5.In_y['CO'].fix(0)
    fs_obj.model.R5.In_y['H2O'].fix(0)
    fs_obj.model.R5.In_y['NH3'].fix(0)

    fs_obj.model.R5.Out_P.fix(101325)

    fs_obj.model.R5.Q.fix(0)
    
def print_summary(fs_obj):
    """
    Print some key results from the model.
    """
    print("Results:")
    print("---Stoichiometric Reactor---")
    print("R1 F_Out = {}".format(fs_obj.model.R1.Out_F.value))
    print("R1 T_Out = {}".format(fs_obj.model.R1.Out_T.value))
    print("R1 P_Out = {}".format(fs_obj.model.R1.Out_P.value))
    for i in fs_obj.model.R1.comp:
        print("R1 y[{}]_out = {}".format(
              i, fs_obj.model.R1.Out_y[i].value))
    print("---CSTR---")
    print("R2 F_Out = {}".format(fs_obj.model.R2.Out_F.value))
    print("R2 T_Out = {}".format(fs_obj.model.R2.Out_T.value))
    print("R2 P_Out = {}".format(fs_obj.model.R2.Out_P.value))
    for i in fs_obj.model.R2.comp:
        print("R2 y[{}]_out = {}".format(
              i, fs_obj.model.R2.Out_y[i].value))
    print("---Equilibrium Reactor---")
    print("R3 F_Out = {}".format(fs_obj.model.R3.Out_F.value))
    print("R3 T_Out = {}".format(fs_obj.model.R3.Out_T.value))
    print("R3 P_Out = {}".format(fs_obj.model.R3.Out_P.value))
    for i in fs_obj.model.R3.comp:
        print("R3 y[{}]_out = {}".format(
              i, fs_obj.model.R3.Out_y[i].value))
    print("---PFR---")
    print("R4 F_Out = {}".format(fs_obj.model.R4.Out_F.value))
    print("R4 T_Out = {}".format(fs_obj.model.R4.Out_T.value))
    print("R4 P_Out = {}".format(fs_obj.model.R4.Out_P.value))
    for i in fs_obj.model.R4.comp:
        print("R4 y[{}]_out = {}".format(
              i, fs_obj.model.R4.Out_y[i].value))
    print("---Gibbs reactor---")
    print("R5 F_Out = {}".format(fs_obj.model.R5.Out_F.value))
    print("R5 T_Out = {}".format(fs_obj.model.R5.Out_T.value))
    print("R5 P_Out = {}".format(fs_obj.model.R5.Out_P.value))
    for i in fs_obj.model.R5.comp:
        print("R5 y[{}]_out = {}".format(
              i, fs_obj.model.R5.Out_y[i].value))

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
    
    rstoic._initialize(fs_obj.model.R1)
    rcstr._initialize(fs_obj.model.R2,T=400,y={'a':0.3,
                                                'b':0.1,
                                                'c':0.2,
                                                'd':0.2,
                                                'e':0.0,
                                                'f':0.2})
    requil._initialize(fs_obj.model.R3,T=400,y={'a':0.2,
                                                'b':0.1,
                                                'c':0.2,
                                                'd':0.2,
                                                'e':0.0,
                                                'f':0.2})
    rpfr._initialize(fs_obj.model.R4,T=400,y={'a':0.2,
                                                'b':0.1,
                                                'c':0.2,
                                                'd':0.2,
                                                'e':0.0,
                                                'f':0.2})
    rgibbs._initialize(fs_obj.model.R5,T=300,y={'H2':0.1,
                                                'N2':0.474,
                                                'O2':0.126,
                                                'CO2':0.0,
                                                'CH4':0.3,
                                                'CO':0.0,
                                                'H2O':0.0,
                                                'NH3':0.0})

    # Solve flowhseet model
    fs_obj.solve(tee=True)
    
    # Print results to:
    #  1) check right number of variables and equations
    #  2) make sure the problem solved successfully.
    fs_obj.display_results()
    print_summary(fs_obj)

if __name__ == "__main__":
    main()
