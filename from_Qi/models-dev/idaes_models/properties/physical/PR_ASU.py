# -*- coding: utf-8 -*-
"""
Package containing necessary parameters for using the Peng-Robinson EoS for
mixtures of Ar, N2 and O2.

Package inherits from Cubic_EoS to solve the phase equilibrium and
thermodynamic properties of the mixture.

@author: alee
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Some more inforation about this module
__author__ = "Andrew Lee"
__version__ = "0.1"

# Import Pyomo
from pyomo.environ import *

# Import IDAES cores
import idaes_models.properties.physical.Cubic_EoS as base_prop

class PropPack(base_prop.PropPack):
    """
    This package contains the necessary parameters and methods for using the
    IDAES general cubic EoS model for systems containing Ar, N2 and O2.
    """
    # List of all chemical species included in package
    comp = ['N2','O2','Ar']

    def _thermo_params(self):
        # Create package parameters
        """ This section of for the users entered parameters relating to
            their system and the species of interest.
            Users must provide values for all the necessary parameters for
            their EoS.
            
            Users need to be sure they provide the correct parameters for the
            desired EoS, most notably the correct set of kappa parameters.
        """
        blk = self

        # Set Equation of State to use
        blk.EoS_ID = 'PR'

        # List of all chemical elements that constitute the chemical species 
        blk.elem = ['N', 'O', 'Ar']

        # Elemental composition of all species
        blk.elem_comp = {'N2'  : {'N':2, 'O':0, 'Ar':0},
                         'O2'  : {'N':0, 'O':2, 'Ar':0},
                         'Ar'  : {'N':0, 'O':0, 'Ar':1},
        }
        # Critical properties - temperature in K, pressure in Pa
        blk.Tc = {'N2' : 126.20,
                  'O2' : 154.58,
                  'Ar' : 150.86
        }
        blk.Pc = {'N2' : 3394387.5,
                  'O2' : 5045985.0,
                  'Ar' : 4873732.5
        }
        # Pitzer accentricity factor (from Prop. Gases & Liquids)
        blk.omega = {'N2' : 0.040,
                     'O2' : 0.021,
                     'Ar' :-0.004
        }
        # Peng-Robinson binary interaction parameters
        blk.kappa = {'N2' : {'N2': 0.0000, 'O2':-0.0119, 'Ar':-0.0026},
                     'O2' : {'N2':-0.0119, 'O2': 0.0000, 'Ar': 0.0104},
                     'Ar' : {'N2':-0.0026, 'O2': 0.0104, 'Ar': 0.0000},
        }
        # Antoine coefficients for ideal vapour (units: bar, K)
        blk.Antoine = {'Ar' : {1 : 3.29555, 2 : 215.240, 3 : -22.233},
                       'N2' : {1 : 3.73620, 2 : 264.651, 3 :  -6.788},
                       'O2' : {1 : 3.95230, 2 : 340.024, 3 :  -4.144}
        }
    def _ideal_pc_props(self):
        # Ideal gas properties
        """ In this section, the user needs to define the methods used to
            calcuate the pure component, ideal gas properties of their species.
        """
        blk = self

        # Ideal gas heat capacity parameters (from Prop. Gases & Liquids)
        blk.cp_ig_param = {
            'N2' : {1:31.128960, 2:-1.356e-2, 3:2.678e-5, 4:-1.167e-8, 5:0},
            'O2' : {1:28.087192, 2:-3.678e-6, 3:1.745e-5, 4:-1.064e-8, 5:0},
            'Ar' : {1:20.790296, 2:-3.209e-5, 3:5.163e-8, 4: 0,        5:0}
        }
        def rule_cp_ig_pc(b, j):
            return b.cp_ig_pc[j] == (b.cp_ig_param[j][5])*b.T**4 \
                    + (b.cp_ig_param[j][4])*b.T**3 \
                    + (b.cp_ig_param[j][3])*b.T**2 \
                    + (b.cp_ig_param[j][2])*b.T \
                    + b.cp_ig_param[j][1]
        blk.eq_cp_ig_pc = Constraint(blk.comp, rule=rule_cp_ig_pc)

        def rule_h_ig_pc(b, j):
            return b.h_ig_pc[j] == (b.cp_ig_param[j][5]/5)*(b.T**5-b.Tref**5) \
                    + (b.cp_ig_param[j][4]/4)*(b.T**4-b.Tref**4) \
                    + (b.cp_ig_param[j][3]/3)*(b.T**3-b.Tref**3) \
                    + (b.cp_ig_param[j][2]/2)*(b.T**2-b.Tref**2) \
                    + b.cp_ig_param[j][1]*(b.T-b.Tref)
        blk.eq_h_ig_pc = Constraint(blk.comp, rule=rule_h_ig_pc)
        
        def rule_s_ig_pc(b, j):
            return b.s_ig_pc[j] == (b.cp_ig_param[j][4]/3)*(b.T**3-b.Tref**3) \
                    + (b.cp_ig_param[j][3]/2)*(b.T**2-b.Tref**2) \
                    + b.cp_ig_param[j][2]*(b.T-b.Tref) \
                    + b.cp_ig_param[j][1]*log(b.T/b.Tref)
        blk.eq_s_ig_pc = Constraint(blk.comp, rule=rule_s_ig_pc)

    def _trans_params(self):
        # Parameters for transport property models
        """ In this section, the user needs to provide the necessary parameters
            for the transport property models.
        """
        blk = self

        # Component molecular weights (mol/kg)
        blk.mw_pc = {'Ar' : 0.04401,
                     'N2' : 0.028014,
                     'O2' : 0.03199806
        }
    def __init__(self, *args, **kwargs):
        """
        Calls the general cubic EoS package to build the property models.
        
        Args
        ====
        plib: an external library containing external function calls
                    for returning the roots of a cubic EoS.
        """
        # Get build arguments
        self.plib = kwargs.pop("plib", 'phys_prop.so')
        self.flag_phases = kwargs.pop("phases", 'VL')
        self.flag_mprop = kwargs.pop("mix_props", 1)

        # Call general cubic EoS constructor
        base_prop.PropPack.__init__(self, *args, **kwargs)
