"""
IDAES air reactor model

Rxn 4-11 regressed from Fig 2 of Huang, et al, J Sust Bioenergy Syst, 2013,3,33-39
Rxn 7a-8a regressed from Fig 3 of Xin, et al, J Nat Gas Chem, 2015, 14(4), 248-253
Rxn ox1-3 regressed from Fig 5 of Fan, et al, AiChE J, 2015, 61, 1, 2-22
Rxn Fe1 and Fe2 are regressed from Fig 4 and 6 of Xin, et al.

Rxn 4  : CO + 3 Fe2O3 <-> CO2 + 2 Fe3O4
Rxn 5  : CO + Fe2O3 <-> CO2 + 2 FeO
Rxn 6  : H2 + 3 Fe2O3 <-> H2O + 2 Fe3O4
Rxn 7  : H2 + Fe2O3 <-> H2O + 2 FeO
Rxn 8  : CH4 + 3 Fe2O3 <-> 2 H2 + CO + 2 Fe3O4
Rxn 9  : CH4 + 4 Fe2O3 <-> 2 H2O + CO2 + 8 FeO
Rxn 10 : C + 3 Fe2O3 <-> CO + 2 Fe3O4
Rxn 11 : C + 2 Fe2O3 <-> CO2 + 4 FeO
Rxn 7a : CH4 + Fe3O4 <-> 3 Fe + CO2 + 2 H2O
Rxn 8a : CH4 + 4 FeO <-> 4 Fe + CO2 + 2 H2O
Rxn ox1: 2 Fe + O2 <-> 2 FeO
Rxn ox2: 6 FeO + O2 <-> 2 Fe3O4
Rxn ox3: 4 Fe3O4 + O2 <-> 6 Fe2O3
Rxn Fe1: Fe3O4 + Fe <-> FeO
Rxn Fe2: Fe2O3 + FeO <-> Fe3O4
"""

from __future__ import division
from __future__ import print_function

from idaes_models.core.unit_model import UnitModel
from idaes_models.core.ports import InPort, OutPort
from idaes_models.core import ProcBlock
from idaes_models.process.conceptd.clc.flowsheet import CLCFlowsheet

from pyomo.environ import *

__author__ = "Tony Burgard, "\
             "Qi Chen <qichen@andrew.cmu.edu>"
__version__ = "1.0.0"


@ProcBlock("AirRxnCocur")
class _AirRxnCocur(UnitModel):
    """
    Equilibrium Air Reactor Model

    """

    def __init__(self, *args, **kwargs):
        """
        Create a reactor object with the component set comp and two outlets.
        Args
        ====
        comp: a Pyomo Set containing a list of components
        Nout: the number of outlets
        """
        # call base class constructor
        super(_AirRxnCocur, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        """Make reactor sets
        """
        b = self
        # Reaction
        b.i = Set(initialize=['4', '5', '6', '7', '8', '9', '10', '11', '7a',
                              '8a', 'ox1', 'ox2', 'ox3', 'Fe1', 'Fe2'])
        # Streams
        b.k = Set(initialize=['GasIn', 'SolidIn', 'GasOut', 'SolidOut'])
        # Elements
        b.e = Set(initialize=['C', 'O', 'H', 'Fe', 'N'])
        # Compounds to ignore for equilibrium calculations
        b.ignore = Set(initialize=['CH4', 'CO', 'CO2'])
        # Gas components
        b.gas = Set(initialize=['CH4', 'CO', 'H2', 'CO2', 'H2O', 'O2', 'N2'])
        # Solid components
        b.solid = Set(initialize=['C', 'Fe2O3', 'Fe3O4', 'FeO', 'Fe'])
        # Gas streams
        b.gas_streams = Set(initialize=['GasIn', 'GasOut'])

    def _build_unit_params(self):
        """
        Make reactor variables
        """
        b = self
        # Equilibrium constants
        b.EqConst = {
            ('4', 'A'): 7.515672,
            ('4', 'B'): 4830.116457,
            ('5', 'A'): 3.295034,
            ('5', 'B'): 818.8999531,
            ('6', 'A'): 11.14578,
            ('6', 'B'): 870.2915707,
            ('7', 'A'): 7.651818,
            ('7', 'B'): -3891.91264,
            ('8', 'A'): 41.68227,
            ('8', 'B'): -26442.10669,
            ('9', 'A'): 55.84036,
            ('9', 'B'): -37059.91603,
            ('10', 'A'): 28.49731,
            ('10', 'B'): -15572.83047,
            ('11', 'A'): 28.11654,
            ('11', 'B'): -19250.97118,
            ('7a', 'A'): 36.52848,
            ('7a', 'B'): -36362.43546,
            ('8a', 'A'): 30.85680,
            ('8a', 'B'): -30715.77616,
            ('ox1', 'A'): -24.0225,
            ('ox1', 'B'): 70418.51336,
            ('ox2', 'A'): -36.0141,
            ('ox2', 'B'): 59803.18789,
            ('ox3', 'A'): -15.9728,
            ('ox3', 'B'): 63988.6020,
            ('Fe1', 'A'): 8.272653,
            ('Fe1', 'B'): -4573.15,
            ('Fe2', 'A'): -0.63963,
            ('Fe2', 'B'): 2697.507
        }
        # Number of elements in each component
        b.ElemComp = {
            ('CH4', 'C'): 1,
            ('CO2', 'C'): 1,
            ('CO', 'C'): 1,
            ('C', 'C'): 1,
            ('CH4', 'H'): 4,
            ('H2', 'H'): 2,
            ('H2O', 'H'): 2,
            ('CO2', 'O'): 2,
            ('CO', 'O'): 1,
            ('H2O', 'O'): 1,
            ('Fe2O3', 'O'): 3,
            ('FeO', 'O'): 1,
            ('Fe3O4', 'O'): 4,
            ('O2', 'O'): 2,
            ('Fe2O3', 'Fe'): 2,
            ('FeO', 'Fe'): 1,
            ('Fe3O4', 'Fe'): 3,
            ('Fe', 'Fe'): 1,
            ('N2', 'N'): 1
        }
        # Standard heat of formation in J/mol (from NIST)
        b.DfHo = {
            ('Fe2O3'): -825500,
            ('Fe3O4'): -1120890,
            ('FeO'): -272040,
            ('Fe'): 0,
            ('C'): 0,
            ('CH4'): -74870,
            ('CO'): -110530,
            ('CO2'): -393510,
            ('H2'): 0,
            ('H2O'): -241826,
            ('N2'): 0,
            ('O2'): 0
        }
        # Shomate equation constants (from NIST)
        b.Sh = {
            ('Fe2O3', 'A'): 110.9362000,
            ('Fe2O3', 'B'): 32.0471400,
            ('Fe2O3', 'C'): -9.1923330,
            ('Fe2O3', 'D'): 0.9015060,
            ('Fe2O3', 'E'): 5.4336770,
            ('Fe2O3', 'F'): -843.1471000,
            ('Fe2O3', 'G'): 228.3548000,
            ('Fe2O3', 'H'): -825.5032000,
            ('Fe3O4', 'A'): 200.8320000,
            ('Fe3O4', 'B'): 0.0000002,
            ('Fe3O4', 'C'): -0.0000001,
            ('Fe3O4', 'D'): 0.0000000,
            ('Fe3O4', 'E'): 0.0000000,
            ('Fe3O4', 'F'): -1174.1350000,
            ('Fe3O4', 'G'): 388.0790000,
            ('Fe3O4', 'H'): -1120.8940000,
            ('FeO', 'A'): 45.7512000,
            ('FeO', 'B'): 18.7855300,
            ('FeO', 'C'): -5.9522010,
            ('FeO', 'D'): 0.8527790,
            ('FeO', 'E'): -0.0812650,
            ('FeO', 'F'): -286.7429000,
            ('FeO', 'G'): 110.3120000,
            ('FeO', 'H'): -272.0441000,
            ('Fe', 'A'): 23.9744900,
            ('Fe', 'B'): 8.3677500,
            ('Fe', 'C'): 0.0002770,
            ('Fe', 'D'): -0.0000860,
            ('Fe', 'E'): -0.0000050,
            ('Fe', 'F'): 0.2680270,
            ('Fe', 'G'): 62.0633600,
            ('Fe', 'H'): 7.7880150,
            ('C', 'A'): 21.1751000,
            ('C', 'B'): -0.8124280,
            ('C', 'C'): 0.4485370,
            ('C', 'D'): -0.0432560,
            ('C', 'E'): -0.0131030,
            ('C', 'F'): 710.3470000,
            ('C', 'G'): 183.8734000,
            ('C', 'H'): 716.6690000,
            ('CH4', 'A'): -0.7030290,
            ('CH4', 'B'): 108.4773000,
            ('CH4', 'C'): -42.5215700,
            ('CH4', 'D'): 5.8627880,
            ('CH4', 'E'): 0.6785650,
            ('CH4', 'F'): -76.8437600,
            ('CH4', 'G'): 158.7163000,
            ('CH4', 'H'): -74.8731000,
            ('CO', 'A'): 25.5675900,
            ('CO', 'B'): 6.0961300,
            ('CO', 'C'): 4.0546560,
            ('CO', 'D'): -2.6713010,
            ('CO', 'E'): 0.1310210,
            ('CO', 'F'): -118.0089000,
            ('CO', 'G'): 227.3665000,
            ('CO', 'H'): -110.5271000,
            ('CO2', 'A'): 24.9973500,
            ('CO2', 'B'): 55.1869600,
            ('CO2', 'C'): -33.6913700,
            ('CO2', 'D'): 7.9483870,
            ('CO2', 'E'): -0.1366380,
            ('CO2', 'F'): -403.6075000,
            ('CO2', 'G'): 228.2431000,
            ('CO2', 'H'): -393.5224000,
            ('H2', 'A'): 33.0661780,
            ('H2', 'B'): -11.3634170,
            ('H2', 'C'): 11.4328160,
            ('H2', 'D'): -2.7728740,
            ('H2', 'E'): -0.1585580,
            ('H2', 'F'): -9.9807970,
            ('H2', 'G'): 172.7079740,
            ('H2', 'H'): 0.0000000,
            ('H2O', 'A'): 30.0920000,
            ('H2O', 'B'): 6.8325140,
            ('H2O', 'C'): 6.7934350,
            ('H2O', 'D'): -2.5344800,
            ('H2O', 'E'): 0.0821390,
            ('H2O', 'F'): -250.8810000,
            ('H2O', 'G'): 223.3967000,
            ('H2O', 'H'): -241.8264000,
            ('N2', 'A'): 19.5058300,
            ('N2', 'B'): 19.8870500,
            ('N2', 'C'): -8.5985350,
            ('N2', 'D'): 1.3697840,
            ('N2', 'E'): 0.5276010,
            ('N2', 'F'): -4.9352020,
            ('N2', 'G'): 212.3900000,
            ('N2', 'H'): 0.0000000,
            ('O2', 'A'): 30.0323500,
            ('O2', 'B'): 8.7729720,
            ('O2', 'C'): -3.9881330,
            ('O2', 'D'): 0.7883130,
            ('O2', 'E'): -0.7415990,
            ('O2', 'F'): -11.3246800,
            ('O2', 'G'): 236.1663000,
            ('O2', 'H'): 0.0000000
        }

    def _build_unit_vars(self):
        b = self
        b.T = Var(doc='Air reactor temperature (K)', domain=NonNegativeReals)
        for s in b.k:  # streams
            setattr(b, 'y_' + s, Var(b.comps, doc='Component mole fractions',
                                     within=NonNegativeReals, bounds=(0, 1)))
        b.logYeq = Var(b.comps,
                       doc='Log of either partial pressure (gas) or mole\
                           fraction (solid) used for equilibrium conditions',
                       bounds=(-20, None))
        b.logYeqSolidout = Var(b.comps,
                               doc='Yeq at solid outlet for countercurrent\
                                    reactors', bounds=(-20, None))
        b.F = Var(b.k, doc='Air reactor molar flow rates in mol/s',
                  domain=NonNegativeReals)
        b.p = Var(b.k, b.comps, doc='Partial pressures',
                  domain=NonNegativeReals)
        b.Pt = Var(doc='Air reactor total pressure', domain=NonNegativeReals)
        b.Hc = Var(
            b.comps, doc='Standard molar enthalpy of component in J/mol', bounds=(-1E12, 1E12))
        b.H = Var(b.k, doc='Enthalpy of stream in J/s', bounds=(-1E12, 1E12))
        b.Q = Var(doc='Air reactor heat duty in J/s')
        b.equip_exists = Var(domain=Binary, initialize=1)

    def _build_unit_constraints(self):

        b = self
        equip = b.equip

        equip.EQ_rxn4 = Constraint(
            expr=b.EqConst['4', 'A'] + b.EqConst['4', 'B'] / b.T ==
            b.logYeq['CO2'] + 2 * b.logYeq['Fe3O4'] - b.logYeq['CO'] - 3 * b.logYeq['Fe2O3'])

        equip.EQ_rxn5 = Constraint(
            expr=b.EqConst['5', 'A'] + b.EqConst['5', 'B'] / b.T ==
            b.logYeq['CO2'] + 2 * b.logYeq['FeO'] - b.logYeq['CO'] - b.logYeq['Fe2O3'])

        equip.EQ_rxn6 = Constraint(
            expr=b.EqConst['6', 'A'] + b.EqConst['6', 'B'] / b.T ==
            b.logYeq['H2O'] + 2 * b.logYeq['Fe3O4'] - b.logYeq['H2'] - 3 * b.logYeq['Fe2O3'])

        equip.EQ_rxn7 = Constraint(
            expr=b.EqConst['7', 'A'] + b.EqConst['7', 'B'] / b.T ==
            b.logYeq['H2O'] + 2 * b.logYeq['FeO'] - b.logYeq['H2'] - b.logYeq['Fe2O3'])

        equip.EQ_rxn8 = Constraint(
            expr=b.EqConst['8', 'A'] + b.EqConst['8', 'B'] / b.T ==
            2 * b.logYeq['H2'] + b.logYeq['CO'] + 2 * b.logYeq['Fe3O4'] - b.logYeq['CH4'] - 3 * b.logYeq['Fe2O3'])

        equip.EQ_rxn9 = Constraint(
            expr=b.EqConst['9', 'A'] + b.EqConst['9', 'B'] / b.T ==
            2 * b.logYeq['H2O'] + b.logYeq['CO2'] + 8 * b.logYeq['FeO'] - b.logYeq['CH4'] - 4 * b.logYeq['Fe2O3'])

        equip.EQ_rxn10 = Constraint(
            expr=b.EqConst['10', 'A'] + b.EqConst['10', 'B'] / b.T ==
            b.logYeq['CO'] + 2 * b.logYeq['Fe3O4'] - b.logYeq['C'] - 3 * b.logYeq['Fe2O3'])

        equip.EQ_rxn11 = Constraint(
            expr=b.EqConst['11', 'A'] + b.EqConst['11', 'B'] / b.T ==
            b.logYeq['CO2'] + 4 * b.logYeq['FeO'] - b.logYeq['C'] - 2 * b.logYeq['Fe2O3'])

        equip.EQ_rxn7a = Constraint(
            expr=b.EqConst['7a', 'A'] + b.EqConst['7a', 'B'] / b.T ==
            3 * b.logYeq['Fe'] + b.logYeq['CO2'] + 2 * b.logYeq['H2O'] - b.logYeq['Fe3O4'] - b.logYeq['CH4'])

        equip.EQ_rxn8a = Constraint(
            expr=b.EqConst['8a', 'A'] + b.EqConst['8a', 'B'] / b.T ==
            4 * b.logYeq['Fe'] + b.logYeq['CO2'] + 2 * b.logYeq['H2O'] - 4 * b.logYeq['FeO'] - b.logYeq['CH4'])

        equip.EQ_rxn10_counter = Constraint(
            expr=b.EqConst['10', 'A'] + b.EqConst['10', 'B'] / b.T ==
            b.logYeq['CO'] + 2 * b.logYeq['Fe3O4'] - b.logYeqSolidout['C'] - 3 * b.logYeq['Fe2O3'])

        equip.EQ_fe1 = Constraint(
            expr=b.EqConst['Fe1', 'A'] + b.EqConst['Fe1', 'B'] / b.T ==
            4 * b.logYeqSolidout['FeO'] - b.logYeqSolidout['Fe3O4'] - b.logYeqSolidout['Fe'])

        equip.EQ_fe2 = Constraint(
            expr=b.EqConst['Fe2', 'A'] + b.EqConst['Fe2', 'B'] / b.T ==
            b.logYeqSolidout['Fe3O4'] - b.logYeqSolidout['Fe2O3'] - b.logYeqSolidout['FeO'])

        equip.EQ_ox1 = Constraint(
            expr=b.EqConst['ox1', 'A'] + b.EqConst['ox1', 'B'] / b.T ==
            2 * b.logYeq['Fe3O4'] - 6 * b.logYeq['FeO'] - b.logYeq['O2'])

        equip.EQ_ox2 = Constraint(
            expr=b.EqConst['ox2', 'A'] + b.EqConst['ox2', 'B'] / b.T ==
            6 * b.logYeq['Fe2O3'] - 4 * b.logYeq['Fe3O4'] - b.logYeq['O2'])

        equip.EQ_ox3 = Constraint(
            expr=b.EqConst['ox3', 'A'] + b.EqConst['ox3', 'B'] / b.T ==
            2 * b.logYeq['FeO'] - 2 * b.logYeq['Fe'] - b.logYeq['O2'])

        @equip.Constraint(self.gas - self.ignore)
        def EqYGasIn(equip, j):
            return b.logYeq[j] + log(b.Pt) == log(b.p['GasIn', j])

        @equip.Constraint(self.gas - self.ignore)
        def EqYGasOut(equip, j):
            return b.logYeq[j] + log(b.Pt) == log(b.p['GasOut', j])

        @equip.Constraint(self.solid - self.ignore)
        def EqYSolidIn(equip, j):
            return b.logYeq[j] + log(b.F['SolidIn']) == log(b.y_SolidIn[j] * b.F['SolidIn'])

        @equip.Constraint(self.solid - self.ignore)
        def EqYSolidOut1(equip, j):
            return b.logYeq[j] + log(b.F['SolidOut']) == log(b.y_SolidOut[j] * b.F['SolidOut'])

        @equip.Constraint(self.solid - self.ignore)
        def EqYSolidOut2(equip, j):
            return b.logYeqSolidout[j] + log(b.F['SolidOut']) == log(b.y_SolidOut[j] * b.F['SolidOut'])

        equip.CH4conversion = Constraint(
            expr=b.F['GasOut'] * b.y_GasOut['CH4'] == 0.05 * b.F['GasIn'] * b.y_GasIn['CH4'])
        equip.COtoH2ratio = Constraint(
            expr=b.y_GasOut['CO'] == b.y_GasOut['H2'] * 1)
        equip.Fullyreduced = Constraint(expr=b.y_SolidOut['Fe3O4'] == 0)
        equip.No_C_Rxn = Constraint(
            expr=b.y_SolidIn['C'] * b.F['SolidIn'] == b.y_SolidOut['C'] * b.F['SolidOut'])

        @equip.Constraint(b.gas_streams, b.gas)
        def EqPartialpressures(equip, k, j):
            return getattr(b, 'y_' + k)[j] * b.Pt == b.p[k, j]

        equip.EQ_rxn4.deactivate()
        equip.EQ_rxn5.deactivate()
        # equip.EQ_rxn6.deactivate()
        equip.EQ_rxn7.deactivate()
        equip.EQ_rxn8.deactivate()
        equip.EQ_rxn9.deactivate()
        equip.EQ_rxn10_counter.deactivate()
        equip.EQ_rxn10.deactivate()
        equip.EQ_rxn11.deactivate()
        equip.EQ_rxn7a.deactivate()
        equip.EQ_rxn8a.deactivate()
        equip.EQ_fe1.deactivate()
        equip.EQ_fe2.deactivate()

        equip.EqYGasIn.deactivate()
        equip.EqYSolidIn.deactivate()
        equip.EqYSolidOut2.deactivate()

        equip.CH4conversion.deactivate()
        equip.COtoH2ratio.deactivate()
        equip.Fullyreduced.deactivate()

        self._build_mass_balance()
        self._build_energy_balance()

    def _build_mass_balance(self):
        """
        Make total and component mass balance equations.
        """
        b = self
        equip = b.equip

        equip.SumFrac_GasOut = Constraint(
            expr=1 == sum(b.y_GasOut[j] for j in b.comps))
        equip.SumFrac_SolidOut = Constraint(
            expr=1 == sum(b.y_SolidOut[j] for j in b.comps))

        def ElemBal_rule(equip, e):
            relevant_comps = [j for j in b.comps if (j, e) in b.ElemComp]
            return b.F['SolidIn'] * sum(b.y_SolidIn[j] * b.ElemComp[j, e] for j in relevant_comps) +\
                b.F['GasIn'] * sum(b.y_GasIn[j] * b.ElemComp[j, e] for j in relevant_comps) ==\
                b.F['SolidOut'] * sum(b.y_SolidOut[j] * b.ElemComp[j, e] for j in relevant_comps) +\
                b.F['GasOut'] * sum(b.y_GasOut[j] *
                                    b.ElemComp[j, e] for j in relevant_comps)
        equip.EqElemBal = Constraint(self.e, rule=ElemBal_rule)

        # enforce sum(y) = 1 at gas inlet
        # @equip.Constraint()
        # def y_sum(equip):
        #     return sum(b.y_GasIn[j] for j in b.comps) == 1

    def _build_energy_balance(self):
        """
        Make energy balance equations.
        """
        b = self
        equip = b.equip

        def EQ_Hc_rule(equip, j):
            return b.Hc[j] == b.DfHo[j] + \
                (b.Sh[j, 'A'] * (b.T / 1000) +
                 b.Sh[j, 'B'] * (b.T / 1000)**2 / 2 +
                 b.Sh[j, 'C'] * (b.T / 1000)**3 / 3 +
                 b.Sh[j, 'D'] * (b.T / 1000)**4 / 4 -
                 b.Sh[j, 'E'] / (b.T / 1000) +
                 b.Sh[j, 'F'] - b.Sh[j, 'H']) * 1000
        equip.EQ_Hc = Constraint(self.comps, rule=EQ_Hc_rule)

        def EQ_H_rule(equip, k):
            return b.H[k] == \
                sum(b.F[k] * getattr(b, 'y_' + k)[j] * b.Hc[j]
                    for j in b.comps)
        equip.EQ_H = Constraint(self.k, rule=EQ_H_rule)

        equip.EQ_ebal = Constraint(expr=b.H['GasIn'] + b.H['SolidIn'] + b.Q ==
                                   b.H['GasOut'] + b.H['SolidOut'])

    def do_init(self):
        '''
        Initialize variable values
        '''
        b = self

        # Initialize properties
        b.T = 1100
        b.Pt = 100000

        for k in b.k:
            b.F[k] = 0
            for j in b.comps:
                getattr(b, 'y_' + k)[j] = 0

        for j in b.comps:
            b.logYeq[j] = 0

        b.F['GasIn'] = 1167.085 * 1000 / 3600
        b.F['SolidIn'] = 736.16 * 1000 / 3600

        # Example inlet conditions for air reactor
        b.y_GasIn['CH4'] = 0
        b.y_GasIn['CO2'] = 0
        b.y_GasIn['CO'] = 0
        b.y_GasIn['H2'] = 1e-5
        b.y_GasIn['H2O'] = 0.010 - 1e-5
        b.y_GasIn['O2'] = 0.209
        b.y_GasIn['N2'] = 0.781

        b.y_SolidIn['Fe2O3'] = 1e-5
        b.y_SolidIn['Fe3O4'] = 0.103
        b.y_SolidIn['FeO'] = 0.897 - 3e-5
        b.y_SolidIn['Fe'] = 1e-5
        b.y_SolidIn['C'] = 1e-5

        for j in b.comps:
            if j in b.gas:
                b.y_SolidIn[j].fix(0)
                b.y_SolidOut[j].fix(0)
            if j in b.solid:
                b.y_GasIn[j].fix(0)
                b.y_GasOut[j].fix(0)
        # b.y_GasOut['O2'].fix(0)

        for j in b.comps:
            b.y_SolidOut[j] = value(b.y_SolidIn[j])
            b.y_GasOut[j] = value(b.y_GasIn[j])

        b.F['SolidOut'] = value(b.F['SolidIn'])
        b.F['GasOut'] = value(b.F['GasIn'])

        for k in b.k:
            for j in b.comps:
                b.p[k, j] = value(getattr(b, 'y_' + k)[j]) * value(b.Pt)

        for j in b.comps:
            if j not in b.ignore:
                b.logYeq[j] = log(
                    (value(b.p['GasIn', j]) + value(b.p['SolidIn', j])))

        '''for k in b.k:
        for j in b.comps:
            if j in b.ignore:
                if j in b.gas:
                    b.p[k,j].fix(0)'''

        for j in b.comps:
            b.Hc[j] = b.DfHo[j] + \
                (b.Sh[j, 'A'] * (value(b.T) / 1000) +
                 b.Sh[j, 'B'] * (value(b.T) / 1000)**2 / 2 +
                 b.Sh[j, 'C'] * (value(b.T) / 1000)**3 / 3 +
                 b.Sh[j, 'D'] * (value(b.T) / 1000)**4 / 4 -
                 b.Sh[j, 'E'] / (value(b.T) / 1000) +
                 b.Sh[j, 'F'] - b.Sh[j, 'H']) * 1000

        for k in b.k:
            b.H[k] = sum(value(b.F[k]) * value(getattr(b, 'y_' + k)
                                               [j]) * value(b.Hc[j]) for j in b.comps)

        b.Q = -value(b.H['GasIn']) - value(b.H['SolidIn']) +\
            value(b.H['GasOut']) + value(b.H['SolidOut'])

    def _build_unit_ports(self):
        """
        Make inlet and outlet stream connectors.
        """
        b = self
        self.add_unit(InPort(name='gas_in', parent=self,
                             F=b.F['GasIn'], y=b.y_GasIn))
        self.add_unit(InPort(name='solid_in', parent=self,
                             F=b.F['SolidIn'], y=b.y_SolidIn))
        self.add_unit(OutPort(name='gas_out', parent=self,
                              F=b.F['GasOut'], y=b.y_GasOut))
        self.add_unit(OutPort(name='solid_out', parent=self,
                              F=b.F['SolidOut'], y=b.y_SolidOut))

    def get_flow_vars(self):
        b = self
        yield b.F


if __name__ == '__main__':
    m = CLCFlowsheet()
    m.max_flow = 20
    m.comp_set = ['C', 'CH4', 'CO', 'H2', 'CO2', 'H2O', 'Fe2O3',
                  'Fe3O4', 'FeO', 'Fe', 'O2', 'N2']
    m.add_unit(AirRxnCocur(name="air_reactor", parent=m))

    m.units.air_reactor.solve(tee=True, options={
                              'outlev': 5, 'tol': 1e-8, 'bound_push': 1e-8, 'mu_init': 1e-8})
    for s in m.units.air_reactor.k:
        getattr(m.units.air_reactor, 'y_' + s).pprint()
    m.units.air_reactor.T.pprint()
    m.units.air_reactor.F.pprint()
    m.units.air_reactor.Hc.pprint()
    m.units.air_reactor.H.pprint()
    m.units.air_reactor.Q.pprint()
#    sr.display()
