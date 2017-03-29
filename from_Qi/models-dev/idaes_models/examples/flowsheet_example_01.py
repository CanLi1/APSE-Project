"""
This is a simple example of using the standard Flowsheet class. Below is a diagram of the model::

                                    A->B       Separate A from B
        Feed --->---[Mixer]--->---[Reactor]---->----[Separator]----> Product
         (A)           |                                  |         (Mostly B)
                       |                                  |
                       |        Recycle (Mostly A)        |
                       +-----------<---------<------------+

As you can see from the diagram, the :class:`Flowsheet` object connects
together three unit models: :class:`Mixer`, :class:`Reactor`,
and :class:`Separator`.
"""
from __future__ import division
from __future__ import print_function

__author__ = "John Eslick"
__version__ = "1.0.0"

import sys
#from pyomo.environ import *
from pyomo.environ import Var, Param, Constraint, Set, NonNegativeReals
import idaes_models.core.unit_model as unit
import idaes_models.core.flowsheet_model as fs

# To keep this simple it is all in one file.  First do classes for the
# mixer, seperator, and reactor unit models.  The do the flowsheet model
# class.  Then in the main function create a flowsheet object and solve.

class Mixer(unit.UnitModel):
    """
    Mixer unit class
    """
    def __init__(self, comp):
        """
        Create a mixer object with the component set comp.
        """
        unit.UnitModel.__init__(self) #call base class constructor
        self.comp = comp
    
    def __call__(self, b, index=0):
        """
        Lets object be used as a Pyomo rule for creating a model block

        Args:
            b: the pyomo block to set up.
            index: the index if creating an indexed block
        """
        
        self.pyomo_block = b
        
        # First make model variables
        b.f_in1 = Var(
            self.comp, domain=NonNegativeReals, initialize=1.0,
            doc="Inlet 1 component flow rates")
        b.f_in2 = Var(
            self.comp, domain=NonNegativeReals, initialize=1.0,
            doc="Inlet 2 component flow rates")
        b.f_out = Var(
            self.comp, domain=NonNegativeReals, initialize=1.0,
            doc="Outlet component flow rates")
        # The temperature and pressure variables may seem a little weird
        # because I am not really using them, but they are part of the
        # standard stream connectors, so need to include them           
        b.T = Var(domain=NonNegativeReals, initialize=300, 
                  doc="Inlet 1 and outlet temperature")
        b.T2 = Var(domain=NonNegativeReals, initialize=300,
                   doc="Inlet 2 temperture")
        b.P = Var(domain=NonNegativeReals, initialize=101325,
                  doc="Inlet 1 and outlet pressure")
        b.P2 = Var(domain=NonNegativeReals, initialize=101325,
                   doc="Inlet 2 pressure")
        
        #Some mass balance equations
        b.mass_balance = Constraint(self.comp, noruleinit=True)
        for i in self.comp:
            b.mass_balance.add(i, expr=b.f_out[i] == 
                               b.f_in1[i] + b.f_in2[i])
        
        # Connectors for 2 inlet streams and 1 outlet
        b.inlet1 = self.one_phase_fluid_flow_connector(b.T, b.P, b.f_in1)
        b.inlet2 = self.one_phase_fluid_flow_connector(b.T2, b.P2, b.f_in2)
        b.outlet = self.one_phase_fluid_flow_connector(b.T, b.P, b.f_out)
        return b
        
class Separator(unit.UnitModel):
    """
    Separator unit class
    """
    def __init__(self, comp):
        """
        Create a mixer object with the component set comp.
        """
        unit.UnitModel.__init__(self) #call base class constructor
        self.comp = comp
    
    def __call__(self, b, index=0):
        """
        Lets object be used as a Pyomo rule for creating a model block
        """
        self.pyomo_block = b
        # Variables T and P are not really used, just for standard ports
        T = b.T = Var(domain=NonNegativeReals, initialize=300,
                      doc="Temperature")
        P = b.P = Var(domain=NonNegativeReals, initialize=101325,
                      doc="Pressure")
        f_in = b.f_in = Var(self.comp, domain=NonNegativeReals, 
                            initialize=1.0, doc="Inlet component flows")
        frac1 = b.frac1 = Var(self.comp, domain=NonNegativeReals, 
                              initialize=1.0, bounds=(0,1),
                              doc="Mole fraction of each component in "\
                                  "outlet 1")
        f_out1 = b.f_out1 = Var(self.comp, domain=NonNegativeReals, 
                                initialize=1.0, doc="Outlet 1 flows")
        f_out2 = b.f_out2 = Var(self.comp, domain=NonNegativeReals, 
                                initialize=1.0, doc="Outlet 2 flows")
        #Constraints
        b.mass_bal = Constraint(self.comp, noruleinit=True)
        b.split = Constraint(self.comp, noruleinit=True)
        for i in self.comp:
            b.split.add(i, expr=f_out1[i] == f_in[i]*frac1[i])
            b.mass_bal.add(i, expr=f_out2[i] == f_in[i] - f_out1[i])
        #Connectors
        b.inlet = self.one_phase_fluid_flow_connector(T, P, f_in)
        b.outlet1 = self.one_phase_fluid_flow_connector(T, P, f_out1)
        b.outlet2 = self.one_phase_fluid_flow_connector(T, P, f_out2)
        return b
        
        
class Reactor(unit.UnitModel):
    """
    Reactor unit class
    """
    def __init__(self, comp,  stoich, key_comp):
        """
        Create a mixer object with the component set comp.
        """
        unit.UnitModel.__init__(self) #call base class constructor
        self.comp = comp
        self.key_comp = key_comp
        self.stoich = stoich
        
    def __call__(self, b, index=0):
        """
        Lets object be used as a Pyomo rule for creating a model block
        """
        self.pyomo_block = b
        #Variables
        b.stoich = Param(self.comp, initialize=self.stoich)
        b.T = Var(domain=NonNegativeReals, initialize=300, 
                  doc="Temperature")
        b.P = Var(domain=NonNegativeReals, initialize=101325,
                  doc="Pressure")
        b.f_in = Var(self.comp, domain=NonNegativeReals, initialize=1.0,
                     doc="Inlet component molar flow rates")
        b.f_out = Var(self.comp, domain=NonNegativeReals, initialize=1.0,
                      doc="Outlet molar flow rates")
        b.conv = Var(domain=NonNegativeReals, bounds=(0,1),
                     doc="Conversion")
        #Constraints
        key = self.key_comp
        b.mass_bal = Constraint(self.comp, noruleinit=True)
        for i in self.comp:
            b.mass_bal.add(i, expr=b.f_out[i] == b.f_in[i] - 
                           b.stoich[i]*b.f_in[key]*b.conv/b.stoich[key])
        #Connectors
        b.inlet = self.one_phase_fluid_flow_connector(b.T, b.P, b.f_in)
        b.outlet = self.one_phase_fluid_flow_connector(b.T, b.P, b.f_out)             
        return b

class Flowsheet(fs.FlowsheetModel):
    def __init__(self, comp, stoich, key_comp, name="F", solver="ipopt"):
        """
        Create a flowsheet model.
        """
        fs.FlowsheetModel.__init__(self, name, solver)
        # Create unit models
        self.add_unit("mix", Mixer(comp))
        self.add_unit("sep", Separator(comp))
        self.add_unit("react", Reactor(comp, stoich, key_comp))
        # Connect unit models
        self.connect("s_mr", self.model.mix.outlet, self.model.react.inlet)
        self.connect("s_rs", self.model.react.outlet, self.model.sep.inlet)
        self.connect("s_sm", self.model.sep.outlet2, self.model.mix.inlet2)
        # Not using T and P here just need them for the standard 
        # connectors so jut fix the mixer T and P, through the 
        # connectors all the other T and P will also be set.
        self.model.mix.T.fix(300)
        self.model.mix.P.fix(101325)

def main():
    """
    Make the flowsheet object, fix some variables, and solve the problem
    """
    comp = Set(initialize=['A', 'B'])
    fs_obj = Flowsheet(comp, {'A':-1,'B':1}, 'A')
    #Fix some variables
    fs_obj.model.mix.f_in1['A'].fix(100)
    fs_obj.model.mix.f_in1['B'].fix(0)
    fs_obj.model.react.conv.fix(0.8)
    fs_obj.model.sep.frac1['A'].fix(0.1)
    fs_obj.model.sep.frac1['B'].fix(0.9)
    
    fs_obj.expand_connectors()
    fs_obj.solve()
    # Print results to 1) check right number of variables and equations
    # 2) make sure the problem solved successfully.
    fs_obj.display_results()
    
    print("Results:")
    print("Product A: {} mol/s".format(fs_obj.model.sep.f_out1['A'].value))
    print("Product B: {} mol/s".format(fs_obj.model.sep.f_out1['B'].value))
    
    
if __name__ == "__main__":
    main()
