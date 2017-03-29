from __future__ import division
import math
from idaes_models.core.unit_model import UnitModel
from idaes_models.core import ProcBlock
from pyomo.environ import Var, NonNegativeReals, Block, Constraint, Binary, Param, NonPositiveReals
from idaes_models.core.ports import InPort, OutPort
from idaes_models.core.util.nonlinear_fcns import exp, log

__author__ = "Qi Chen <qichen@andrew.cmu.edu>"


@ProcBlock("PackedBed")
class _PackedBed(UnitModel):
    def __init__(self, *args, **kwargs):
        super(_PackedBed, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        super(_PackedBed, self)._build_unit_sets()

    def _build_unit_params(self):
        super(_PackedBed, self)._build_unit_params()
        b = self
        # Flow cost factor: multiplied with the total liquid flow into unit
        b.flow_factor = Param(initialize=360 * 1.6)
        # Size cost factor: multipled with volume
        b.size_factor = Param(initialize=5)
        # Mass transfer coefficient
        b.kL = Param(initialize=0.0032)
        # Some mass transfer term
        b.a = Param(initialize=400)
        # liquid film mean lifetime (seconds)
        b.t = Param(initialize=2)

    def _build_unit_vars(self):
        super(_PackedBed, self)._build_unit_vars()
        b = self
        b.flow_vap_in = Var(self.comps, domain=NonNegativeReals, initialize=0.1, bounds=(0, self._max_flow))
        b.flow_vap_out = Var(self.comps, domain=NonNegativeReals, initialize=0.1, bounds=(0, self._max_flow))
        b.flow_liq_in = Var(self.comps, domain=NonNegativeReals, initialize=0.1, bounds=(0, self._max_flow))
        b.flow_liq_out = Var(self.comps, domain=NonNegativeReals, initialize=0.1, bounds=(0, self._max_flow))
        b.total_flow_vap_in = Var(domain=NonNegativeReals, initialize=0.1, bounds=(0, self._max_flow))
        b.total_flow_vap_out = Var(domain=NonNegativeReals, initialize=0.1, bounds=(0, self._max_flow))
        b.total_flow_liq_in = Var(domain=NonNegativeReals, initialize=0.1, bounds=(0, self._max_flow))
        b.total_flow_liq_out = Var(domain=NonNegativeReals, initialize=0.1, bounds=(0, self._max_flow))
        b.frac_vap_in = Var(self.comps, domain=NonNegativeReals, initialize=0.1, bounds=(1E-6, 1))
        b.frac_vap_out = Var(self.comps, domain=NonNegativeReals, initialize=0.1, bounds=(1E-6, 1))
        b.equip_exists = Var(domain=Binary, initialize=1)
        b.equip_cost = Var(domain=NonNegativeReals, initialize=1, bounds=(0, 1000))
        # volume of packed bed
        b.V = Var(domain=NonNegativeReals, initialize=1, bounds=(0.01, 40))
        # neight of packed bed
        b.h = Var(domain=NonNegativeReals, initialize=1, bounds=(0.1, 20))
        # radius of packed bed
        b.R = Var(domain=NonNegativeReals, initialize=0.557, bounds=(0.1, 1))
        # ratio of radial thickness to axial thickness of packing
        b.RH = Var(domain=NonNegativeReals, initialize=1, bounds=(1 / 30, 1 / 15))
        b.exp_base = Var(domain=NonNegativeReals, initialize=1, bounds=(0, 23))
        # linear approx of R^2
        b.R_sq = Var(domain=NonNegativeReals, initialize=1, bounds=(0, 1))
        # linear approx of height * R^2
        b.HR_sq = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of total_flow_vap_in * exp_base
        b.FvapinExpbase = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of exp(exp_base) + 1
        b.exp_val = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of frac_vap_out[A] * (1 - frac_vap_in[A])
        b.frac_vap_inout = Var(domain=NonNegativeReals, initialize=1)

    def _build_unit_constraints(self):
        super(_PackedBed, self)._build_unit_constraints()
        b = self
        equip = b.equip

        def total_flow_vap_in_calc(equip):
            b = equip.parent_block()
            return b.total_flow_vap_in == \
                sum(b.flow_vap_in[c] for c in self.comps)
        equip.total_flow_vap_in_calc = Constraint(rule=total_flow_vap_in_calc)

        def total_flow_liq_in_calc(equip):
            b = equip.parent_block()
            return b.total_flow_liq_in == \
                sum(b.flow_liq_in[c] for c in self.comps)
        equip.total_flow_liq_in_calc = Constraint(rule=total_flow_liq_in_calc)

        def total_flow_vap_out_calc(equip):
            b = equip.parent_block()
            return b.total_flow_vap_out == \
                sum(b.flow_vap_out[c] for c in self.comps)
        equip.total_flow_vap_out_calc = Constraint(rule=total_flow_vap_out_calc)

        def total_flow_liq_out_calc(equip):
            b = equip.parent_block()
            return b.total_flow_liq_out == \
                sum(b.flow_liq_out[c] for c in self.comps)
        equip.total_flow_liq_out_calc = Constraint(rule=total_flow_liq_out_calc)

        def cost(equip):
            b = equip.parent_block()
            return b.equip_cost == b.flow_factor * b.total_flow_liq_in + b.size_factor * b.V
        equip.cost = Constraint(rule=cost)

        def flow_in_out_vap(equip):
            b = equip.parent_block()
            return b.total_flow_vap_in == b.total_flow_vap_out
        equip.flow_in_out_vap = Constraint(rule=flow_in_out_vap)

        def flow_in_out_liq(equip):
            b = equip.parent_block()
            return b.total_flow_liq_in == b.total_flow_liq_out
        equip.flow_in_out_liq = Constraint(rule=flow_in_out_liq)

        def pb_liq_flow(equip):
            b = equip.parent_block()
            return b.total_flow_liq_in * 100 == b.total_flow_vap_in
        equip.pb_liq_flow = Constraint(rule=pb_liq_flow)

        def frac_vap_in_calc(equip, c):
            b = equip.parent_block()
            return b.frac_vap_in[c] * b.total_flow_vap_in == b.flow_vap_in[c]
        equip.frac_vap_in_calc = Constraint(self.comps, rule=frac_vap_in_calc)

        def frac_vap_out_calc(equip, c):
            b = equip.parent_block()
            return b.frac_vap_out[c] * b.total_flow_vap_out == b.flow_vap_out[c]
        equip.frac_vap_out_calc = Constraint(self.comps, rule=frac_vap_out_calc)

        def exp_base_calc(equip):
            b = equip.parent_block()
            return b.total_flow_vap_in * b.exp_base == b.kL * b.a * b.h * math.pi * b.R ** 2
        equip.exp_base_calc = Constraint(rule=exp_base_calc)

        def mass_transfer(equip):
            b = equip.parent_block()
            return 10 * b.frac_vap_in['A'] == 10 * b.frac_vap_out['A'] * (exp(b.exp_base) + 1) * (1 - b.frac_vap_in['A'])
        equip.mass_transfer = Constraint(rule=mass_transfer)

        def volume_calc(equip):
            b = equip.parent_block()
            return b.V == math.pi * b.R ** 2 * b.h
        equip.volume_calc = Constraint(rule=volume_calc)

        def RH_calc(equip):
            b = equip.parent_block()
            return b.RH * b.h == b.R
        equip.RH_calc = Constraint(rule=RH_calc)

    def _build_unit_ports(self):
        super(_PackedBed, self)._build_unit_ports()
        b = self
        self.add_unit(InPort(name='vap_in', parent=self, fc=b.flow_vap_in))
        self.add_unit(OutPort(name='vap_out', parent=self, fc=b.flow_vap_out))
        self.add_unit(InPort(name='liq_in', parent=self, fc=b.flow_liq_in))
        self.add_unit(OutPort(name='liq_out', parent=self, fc=b.flow_liq_out))

    def get_flow_vars(self):
        yield self.flow_vap_in
        yield self.flow_vap_out
        yield self.flow_liq_in
        yield self.flow_liq_out

    def apply_linear_relaxations(self, nsegs=1, **kwargs):
        b = self
        lin_cuts = b.lin_cuts
        self.lin_cuts_nsegs = nsegs
        super(_PackedBed, self).apply_linear_relaxations(nsegs=nsegs, **kwargs)

        from idaes_models.core.util.mccormick import add_mccormick_relaxation
        from idaes_models.core.util.var import tighten_mc_var
        from functools import partial
        tighten_mc = partial(tighten_mc_var, block_bounds=self.block_bounds)
        from idaes_models.core.util.convex import add_convex_relaxation

        lin_cuts.mc_FvapinX = Block()
        # setup_mccormick_cuts(lin_cuts.mc_FvapinX, 'Fvapin', 1, self.comps)
        for c in self.comps:
            tighten_mc(b.flow_vap_in[c], b.total_flow_vap_in, b.frac_vap_in[c])
            add_mccormick_relaxation(lin_cuts.mc_FvapinX, b.flow_vap_in[c], b.total_flow_vap_in, b.frac_vap_in[c], nsegs, c, b.equip_exists, block_bounds=self.block_bounds)
            # add_mccormick_cut(lin_cuts.mc_FvapinX, 'Fvapin', c, b.flow_vap_in[c], b.total_flow_vap_in, b.frac_vap_in[c], b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_FvapoutX = Block()
        # setup_mccormick_cuts(lin_cuts.mc_FvapoutX, 'Fvapout', 1, self.comps)
        for c in self.comps:
            tighten_mc(b.flow_vap_out[c], b.total_flow_vap_out, b.frac_vap_out[c])
            add_mccormick_relaxation(lin_cuts.mc_FvapoutX, b.flow_vap_out[c], b.total_flow_vap_out, b.frac_vap_out[c], nsegs, c, b.equip_exists, block_bounds=self.block_bounds)
            # add_mccormick_cut(lin_cuts.mc_FvapoutX, 'Fvapout', c, b.flow_vap_out[c], b.total_flow_vap_out, b.frac_vap_out[c], b.equip_exists, block_bounds=self.block_bounds)

        def R_sq(R):
            return R ** 2

        def dR_sq(R):
            return 2 * R
        lin_cuts.env_R_sq = Block()
        # add_convex_relaxation(lin_cuts.env_R_sq, b.R_sq, b.R, R_sq, dR_sq, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_convex_relaxation(lin_cuts.env_R_sq, b.R_sq, b.R, R_sq, dR_sq, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')

        lin_cuts.mc_HR_sq = Block()
        tighten_mc(b.HR_sq, b.h, b.R_sq)
        # setup_mccormick_cuts(lin_cuts.mc_HR_sq, 'HR_sq', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_HR_sq, 'HR_sq', None, b.HR_sq, b.h, b.R_sq, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_HR_sq, b.HR_sq, b.h, b.R_sq, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_FvapinExpbase = Block()
        tighten_mc(b.FvapinExpbase, b.total_flow_vap_in, b.exp_base)
        # setup_mccormick_cuts(lin_cuts.mc_FvapinExpbase, 'FvapinExpbase', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_FvapinExpbase, 'FvapinExpbase', None, b.FvapinExpbase, b.total_flow_vap_in, b.exp_base, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_FvapinExpbase, b.FvapinExpbase, b.total_flow_vap_in, b.exp_base, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.linear_exp_base_calc = Constraint(expr=b.FvapinExpbase == b.kL * b.a * math.pi * b.HR_sq)

        def exp_val(x):
            return exp(x) + 1

        def dexp_val(x):
            return exp(x)
        lin_cuts.env_exp_val = Block()
        # add_convex_relaxation(lin_cuts.env_exp_val, b.exp_val, b.exp_base, exp_val, dexp_val, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_convex_relaxation(lin_cuts.env_exp_val, b.exp_val, b.exp_base, exp_val, dexp_val, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')

        lin_cuts.mc_frac_vap_inout = Block()
        tighten_mc(b.frac_vap_inout, b.frac_vap_out['A'], 1 - b.frac_vap_in['A'])
        # setup_mccormick_cuts(lin_cuts.mc_frac_vap_inout, 'frac_vap_inout', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_frac_vap_inout, 'frac_vap_inout', None, b.frac_vap_inout, b.frac_vap_out['A'], 1 - b.frac_vap_in['A'], b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_frac_vap_inout, b.frac_vap_inout, b.frac_vap_out['A'], 1 - b.frac_vap_in['A'], nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_mass_transfer = Block()
        tighten_mc(b.frac_vap_in['A'], b.frac_vap_inout, b.exp_val)
        # setup_mccormick_cuts(lin_cuts.mc_mass_transfer, 'mass_transfer', nsegs=1)
        # # frac_vap_inout * exp_val = frac_vap_in[A]
        # add_mccormick_cut(lin_cuts.mc_mass_transfer, 'mass_transfer', None, b.frac_vap_in['A'], b.frac_vap_inout, b.exp_val, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_mass_transfer, b.frac_vap_in['A'], b.frac_vap_inout, b.exp_val, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.linear_volume_calc = Constraint(expr=b.V == math.pi * b.HR_sq)

        lin_cuts.mc_RH = Block()
        tighten_mc(b.R, b.RH, b.h)
        # setup_mccormick_cuts(lin_cuts.mc_RH, 'RH_calc', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_RH, 'RH_calc', None, b.R, b.RH, b.h, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_RH, b.R, b.RH, b.h, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

    def get_vars_to_bound(self):
        yield self.R, self.block_bounds
        yield self.h
        yield self.R_sq
        yield self.total_flow_vap_in, self.block_bounds
        yield self.exp_base
        yield self.frac_vap_out, self.block_bounds
        yield self.frac_vap_in, self.block_bounds
        yield self.frac_vap_inout
        yield self.RH

    def generate_cut_gen_problem(self):
        raise NotImplementedError()

    def apply_self_proj_cut(self):
        raise NotImplementedError()


@ProcBlock("RotatingPackedBed")
class _RotatingPackedBed(UnitModel):
    def __init__(self, *args, **kwargs):
        super(_RotatingPackedBed, self).__init__(*args, **kwargs)

    def _build_unit_sets(self):
        super(_RotatingPackedBed, self)._build_unit_sets()

    def _build_unit_params(self):
        super(_RotatingPackedBed, self)._build_unit_params()
        b = self
        b.rotate_factor = Param(initialize=2.2)
        # Flow cost factor: multiplied with the total liquid flow into unit
        b.flow_factor = Param(initialize=360)
        # Size cost factor: multipled with volume
        b.size_factor = Param(initialize=121)
        # Some mass transfer term
        b.a = Param(initialize=400)
        # density of 1
        b.rho = Param(initialize=1)
        b.k_ov = Param(initialize=1)
        b.D = Param(initialize=0.003364)
        b.Ns = Param(initialize=90)

    def _build_unit_vars(self):
        super(_RotatingPackedBed, self)._build_unit_vars()
        b = self
        b.flow_vap_in = Var(self.comps, domain=NonNegativeReals, initialize=1, bounds=(0, self._max_flow))
        b.flow_vap_out = Var(self.comps, domain=NonNegativeReals, initialize=1, bounds=(0, self._max_flow))
        b.flow_liq_in = Var(self.comps, domain=NonNegativeReals, initialize=1, bounds=(0, self._max_flow))
        b.flow_liq_out = Var(self.comps, domain=NonNegativeReals, initialize=1, bounds=(0, self._max_flow))
        b.total_flow_vap_in = Var(domain=NonNegativeReals, initialize=1, bounds=(0, self._max_flow))
        b.total_flow_vap_out = Var(domain=NonNegativeReals, initialize=1, bounds=(0, self._max_flow))
        b.total_flow_liq_in = Var(domain=NonNegativeReals, initialize=1, bounds=(0, self._max_flow))
        b.total_flow_liq_out = Var(domain=NonNegativeReals, initialize=1, bounds=(0, self._max_flow))
        b.frac_vap_in = Var(self.comps, domain=NonNegativeReals, initialize=0.01, bounds=(1E-6, 1 - 1E-6))
        b.frac_vap_out = Var(self.comps, domain=NonNegativeReals, initialize=0.001, bounds=(1E-6, 1 - 1E-6))
        b.equip_exists = Var(domain=Binary, initialize=1)
        b.equip_cost = Var(domain=NonNegativeReals, initialize=1, bounds=(0, 1000))
        b.rotate_cost = Var(domain=NonNegativeReals, initialize=1, bounds=(0, 1000))
        # volume of rotating packed bed
        b.V = Var(domain=NonNegativeReals, initialize=1, bounds=(0.005, 1))
        # neight of rotating packed bed
        b.h = Var(domain=NonNegativeReals, initialize=0.266, bounds=(0.01, 0.4))
        # radius of rotating packed bed (inner and outer radius)
        b.Rin = Var(domain=NonNegativeReals, initialize=0.147, bounds=(0.05, 0.26))
        b.Rout = Var(domain=NonNegativeReals, initialize=0.441, bounds=(0.2, 2))
        # R = Rout - Rin
        b.R = Var(domain=NonNegativeReals, initialize=1, bounds=(0.01, 2))
        b.speed = Var(domain=NonNegativeReals, initialize=10, bounds=(10, 40))
        b.erf = Var(domain=NonNegativeReals, initialize=1, bounds=(0, 1))
        # b.erft = Var(domain=NonNegativeReals, initialize=1, bounds=(0, 1))
        b.t = Var(domain=NonNegativeReals, initialize=0.104, bounds=(0.001, 1.44))
        b.kL = Var(domain=NonNegativeReals, initialize=0.095, bounds=(0, 10))
        # velocity of liquid stream entering RPB
        b.u = Var(domain=NonNegativeReals, initialize=0.031, bounds=(0.0001, 1))
        # ratio of radial thickness to axial thickness of packing
        b.RH = Var(domain=NonNegativeReals, initialize=20, bounds=(0.88, 1.2))
        b.exp_base = Var(domain=NonNegativeReals, initialize=1, bounds=(0, 100))
        # omega (speed) squared
        b.wSq = Var(domain=NonNegativeReals, initialize=1)
        # w^2 * h
        b.wSqh = Var(domain=NonNegativeReals, initialize=1)
        # Rout^4
        b.RoutFourth = Var(domain=NonNegativeReals, initialize=1)
        # Rin^4
        b.RinFourth = Var(domain=NonNegativeReals, initialize=1)
        # Rout^4 - Rin^4
        b.RFourthDiff = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of omega (speed) * R
        b.wR = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of 0.02107 * (w*R)^0.5488
        b.wRexp = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of total_flow_liq_in ** 0.2279
        b.FLiqInExp = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of t^0.5
        b.t1 = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of t^1.5
        b.t2 = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of t^2.5
        b.t3 = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of t^3.5
        b.t4 = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of t * erf
        b.terf = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of exp(-k_ov * t)
        b.expt = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of t^0.5 * exp(-k_ov * t)
        b.t1expt = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of t * kL
        b.tkL = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of Rout^2
        b.RoutSq = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of Rin^2
        b.RinSq = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of Rout^2 - Rin^2
        b.RSqDiff = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of h * (Rout^2 - Rin^2)
        b.hRSqDiff = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of total_flow_vap_in * exp_base
        b.FvapinExpbase = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of kL * h * (Rout^2 - Rin^2)
        b.kLhRSqDiff = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of (1 - frac_vap_in[A]) * frac_vap_out[A]
        b.frac_vap_inout = Var(domain=NonNegativeReals, initialize=1)
        # linear approx of frac_vap_in[A] - frac_vap_inout
        b.frac_diff = Var(domain=NonNegativeReals, initialize=0.001, bounds=(1E-9, 1))
        # linear approx of ln(frac_vap_in[A] - frac_vap_inout)
        b.lnFracDiff = Var(domain=NonPositiveReals, initialize=-1)
        # linear approx of ln(frac_vap_inout)
        b.lnFrac = Var(domain=NonPositiveReals, initialize=-1)

    def _build_unit_constraints(self):
        super(_RotatingPackedBed, self)._build_unit_constraints()
        b = self
        equip = b.equip
        nl_con = equip.nl_con = Block()

        def total_flow_vap_in_calc(equip):
            b = equip.parent_block()
            return b.total_flow_vap_in == \
                sum(b.flow_vap_in[c] for c in self.comps)
        equip.total_flow_vap_in_calc = Constraint(rule=total_flow_vap_in_calc)

        def total_flow_liq_in_calc(equip):
            b = equip.parent_block()
            return b.total_flow_liq_in == \
                sum(b.flow_liq_in[c] for c in self.comps)
        equip.total_flow_liq_in_calc = Constraint(rule=total_flow_liq_in_calc)

        def total_flow_vap_out_calc(equip):
            b = equip.parent_block()
            return b.total_flow_vap_out == \
                sum(b.flow_vap_out[c] for c in self.comps)
        equip.total_flow_vap_out_calc = Constraint(rule=total_flow_vap_out_calc)

        def total_flow_liq_out_calc(equip):
            b = equip.parent_block()
            return b.total_flow_liq_out == \
                sum(b.flow_liq_out[c] for c in self.comps)
        equip.total_flow_liq_out_calc = Constraint(rule=total_flow_liq_out_calc)

        def cost(equip):
            b = equip.parent_block()
            return b.equip_cost == b.rotate_cost + b.flow_factor * b.total_flow_liq_in + b.size_factor * b.V
        equip.cost = Constraint(rule=cost)

        def R_defn(equip):
            b = equip.parent_block()
            return b.R == b.Rout - b.Rin
        equip.R_defn = Constraint(rule=R_defn)

        def R_out_min(equip):
            b = equip.parent_block()
            return b.Rout >= 2 * b.Rin
        equip.R_out_min = Constraint(rule=R_out_min)

        def R_out_max(equip):
            b = equip.parent_block()
            return b.Rout <= 3 * b.Rin
        equip.R_out_max = Constraint(rule=R_out_max)

        def flow_in_out_vap(equip):
            b = equip.parent_block()
            # return b.flow_vap_in['water'] == b.flow_vap_out['water']
            return b.total_flow_vap_in == b.total_flow_vap_out
        equip.flow_in_out_vap = Constraint(rule=flow_in_out_vap)

        def flow_in_out_liq(equip):
            b = equip.parent_block()
            return b.total_flow_liq_in == b.total_flow_liq_out
        equip.flow_in_out_liq = Constraint(rule=flow_in_out_liq)

        def rpb_liq_flow(equip):
            b = equip.parent_block()
            return b.total_flow_liq_in * 200 == b.total_flow_vap_in
        equip.rpb_liq_flow = Constraint(rule=rpb_liq_flow)

        def rotate_cost(nl_con):
            b = nl_con.parent_block().parent_block()
            return b.rotate_cost == b.rotate_factor * (b.Rout ** 4 - b.Rin ** 4) / 4 * b.rho * math.pi * b.h * b.speed ** 2
        nl_con.rotate_cost = Constraint(rule=rotate_cost)

        def frac_vap_in_calc(nl_con, c):
            b = nl_con.parent_block().parent_block()
            return b.frac_vap_in[c] * b.total_flow_vap_in == b.flow_vap_in[c]
        nl_con.frac_vap_in_calc = Constraint(self.comps, rule=frac_vap_in_calc)

        def frac_vap_out_calc(nl_con, c):
            b = nl_con.parent_block().parent_block()
            return b.frac_vap_out[c] * b.total_flow_vap_out == b.flow_vap_out[c]
        nl_con.frac_vap_out_calc = Constraint(self.comps, rule=frac_vap_out_calc)

        def u_defn(nl_con):
            b = nl_con.parent_block().parent_block()
            return b.u == 0.02107 * b.total_flow_liq_in ** 0.2279 * (b.speed * b.R) ** 0.5488
        nl_con.u_defn = Constraint(rule=u_defn)

        def t_defn(nl_con):
            b = nl_con.parent_block().parent_block()
            return b.t * b.u * b.Ns == b.R
        nl_con.t_defn = Constraint(rule=t_defn)

        def erf_calc(nl_con):
            """Error function calculation

            We use the Taylor series expansion for the error function. This
            requires the value of sqrt(t) to be between 0 and ~1.2 for
            reasonable accuracy.
            """
            b = nl_con.parent_block().parent_block()
            return b.erf == 2 / (math.pi) ** 0.5 * (b.t ** 0.5 - b.t ** 1.5 / 3 + b.t ** 2.5 / 10 - b.t ** 3.5 / 42)
        nl_con.erf_calc = Constraint(rule=erf_calc)

        def kL_calc(nl_con):
            b = nl_con.parent_block().parent_block()
            return b.t * b.kL == (b.k_ov * b.D) ** 0.5 * (b.t * b.erf + (b.t / (math.pi * b.k_ov)) ** 0.5 * exp(-b.k_ov * b.t) + 1 / (2 * b.k_ov) * b.erf)
        nl_con.kL_calc = Constraint(rule=kL_calc)

        def exp_base_calc(nl_con):
            b = nl_con.parent_block().parent_block()
            return b.total_flow_vap_in * b.exp_base == b.kL * b.a * b.h * math.pi * (b.Rout ** 2 - b.Rin ** 2)
        nl_con.exp_base_calc = Constraint(rule=exp_base_calc)

        def mass_transfer(nl_con):
            b = nl_con.parent_block().parent_block()
            # return b.frac_vap_in['A'] == b.frac_vap_out['A'] * (exp(b.exp_base) + 1) * (1 - b.frac_vap_in['A'])
            return log((b.frac_vap_in['A'] - (1 - b.frac_vap_in['A']) * b.frac_vap_out['A']) ** 2) / 2 - log(1 - b.frac_vap_in['A']) - log(b.frac_vap_out['A']) == b.exp_base
        nl_con.mass_transfer = Constraint(rule=mass_transfer)

        def volume_calc(nl_con):
            b = nl_con.parent_block().parent_block()
            return b.V == math.pi * (b.Rout ** 2 - b.Rin ** 2) * b.h
        nl_con.volume_calc = Constraint(rule=volume_calc)

        def RH_calc(nl_con):
            b = nl_con.parent_block().parent_block()
            return b.RH * b.h == b.R
        nl_con.RH_calc = Constraint(rule=RH_calc)

    def _build_unit_ports(self):
        super(_RotatingPackedBed, self)._build_unit_ports()
        b = self
        self.add_unit(InPort(name='vap_in', parent=self, fc=b.flow_vap_in))
        self.add_unit(OutPort(name='vap_out', parent=self, fc=b.flow_vap_out))
        self.add_unit(InPort(name='liq_in', parent=self, fc=b.flow_liq_in))
        self.add_unit(OutPort(name='liq_out', parent=self, fc=b.flow_liq_out))

    def get_flow_vars(self):
        yield self.flow_vap_in
        yield self.flow_vap_out
        yield self.flow_liq_in
        yield self.flow_liq_out

    def apply_linear_relaxations(self, nsegs=1):
        b = self
        lin_cuts = b.lin_cuts
        self.lin_cuts_nsegs = nsegs
        super(_RotatingPackedBed, self).apply_linear_relaxations(nsegs=nsegs)

        from idaes_models.core.util.mccormick import add_mccormick_relaxation
        # from idaes_models.core.util.var import max_ub, min_lb
        from idaes_models.core.util.var import \
            tighten_block_bound as tbb, lb, ub, \
            tighten_var_bound as tighten_vb, tighten_mc_var, is_fixed_by_bounds
        from functools import partial
        tighten_bb = partial(tbb, block_bounds=self.block_bounds)
        tighten_mc = partial(tighten_mc_var, block_bounds=self.block_bounds)
        lbb = partial(lb, block_bounds=self.block_bounds)
        ubb = partial(ub, block_bounds=self.block_bounds)
        from idaes_models.core.util.concave import add_concave_relaxation
        from idaes_models.core.util.convex import add_convex_relaxation

        lin_cuts.mc_FvapinX = Block()
        # setup_mccormick_cuts(lin_cuts.mc_FvapinX, 'Fvapin', 1, self.comps)
        for c in self.comps:
            tighten_mc(b.flow_vap_in[c], b.total_flow_vap_in, b.frac_vap_in[c])
            # add_mccormick_cut(lin_cuts.mc_FvapinX, 'Fvapin', c, b.flow_vap_in[c], b.total_flow_vap_in, b.frac_vap_in[c], b.equip_exists, block_bounds=self.block_bounds)
            add_mccormick_relaxation(lin_cuts.mc_FvapinX, b.flow_vap_in[c], b.total_flow_vap_in, b.frac_vap_in[c], nsegs, c, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_FvapoutX = Block()
        # setup_mccormick_cuts(lin_cuts.mc_FvapoutX, 'Fvapout', 1, self.comps)
        for c in self.comps:
            tighten_mc(b.flow_vap_out[c], b.total_flow_vap_out, b.frac_vap_out[c])
            # add_mccormick_cut(lin_cuts.mc_FvapoutX, 'Fvapout', c, b.flow_vap_out[c], b.total_flow_vap_out, b.frac_vap_out[c], b.equip_exists, block_bounds=self.block_bounds)
            add_mccormick_relaxation(lin_cuts.mc_FvapoutX, b.flow_vap_out[c], b.total_flow_vap_out, b.frac_vap_out[c], nsegs, c, b.equip_exists, block_bounds=self.block_bounds)

        def wSq(speed):
            return speed ** 2

        def dwSq(speed):
            return 2 * speed
        lin_cuts.env_wSq = Block()
        # add_convex_relaxation(lin_cuts.env_wSq, b.wSq, b.speed, wSq, dwSq, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_convex_relaxation(lin_cuts.env_wSq, b.wSq, b.speed, wSq, dwSq, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')

        # wSq * h
        lin_cuts.mc_wSqh = Block()
        # b.wSqh.setlb(min_lb(b.wSq, b.h))
        # b.wSqh.setub(max_ub(b.wSq, b.h))
        tighten_mc(b.wSqh, b.wSq, b.h)
        # setup_mccormick_cuts(lin_cuts.mc_wSqh, 'wSqh', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_wSqh, 'wSqh', None, b.wSqh, b.wSq, b.h, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_wSqh, b.wSqh, b.wSq, b.h, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        # Rout^4
        def RoutFourth(Rout):
            return Rout ** 4

        def dRoutFourth(Rout):
            return 4 * Rout ** 3
        lin_cuts.env_RoutFourth = Block()
        # add_convex_relaxation(lin_cuts.env_RoutFourth, b.RoutFourth, b.Rout, RoutFourth, dRoutFourth, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_convex_relaxation(lin_cuts.env_RoutFourth, b.RoutFourth, b.Rout, RoutFourth, dRoutFourth, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')

        # Rin^4
        def RinFourth(Rin):
            return Rin ** 4

        def dRinFourth(Rin):
            return 4 * Rin ** 3
        lin_cuts.env_RinFourth = Block()
        # add_convex_relaxation(lin_cuts.env_RinFourth, b.RinFourth, b.Rin, RinFourth, dRinFourth, b.equip_exists, self.block_bounds)
        add_convex_relaxation(lin_cuts.env_RinFourth, b.RinFourth, b.Rin, RinFourth, dRinFourth, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')

        # Rout^4 - Rin^4
        lin_cuts.linear_RFourthDiff = Constraint(expr=b.RFourthDiff == b.RoutFourth - b.RinFourth)
        # b.RFourthDiff.setlb(0)
        # b.RFourthDiff.setub(b.RoutFourth.ub - b.RinFourth.lb)
        tighten_vb(b.RFourthDiff, (0, ub(b.RoutFourth) - lb(b.RinFourth)))
        tighten_bb(b.RFourthDiff, (0, ubb(b.RoutFourth) - lbb(b.RinFourth)))

        # linear_rotate_cost
        lin_cuts.mc_rotate_cost = Block()
        tighten_mc(b.rotate_cost, b.rotate_factor / 4 * b.RFourthDiff,
                   b.rho * math.pi * b.wSqh)
        # setup_mccormick_cuts(lin_cuts.mc_rotate_cost, 'rotate_cost', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_rotate_cost, 'rotate_cost', None, b.rotate_cost, b.rotate_factor / 4 * b.RFourthDiff, b.rho * math.pi * b.wSqh, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_rotate_cost, b.rotate_cost, b.rotate_factor / 4 * b.RFourthDiff, b.rho * math.pi * b.wSqh, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_wR = Block()
        # b.wR.setlb(min_lb(b.speed, b.R))
        # b.wR.setub(max_ub(b.speed, b.R))
        tighten_mc(b.wR, b.speed, b.R)
        # setup_mccormick_cuts(lin_cuts.mc_wR, 'wR', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_wR, 'wR', None, b.wR, b.speed, b.R, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_wR, b.wR, b.speed, b.R, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        def wRexp(wR):
            return 0.02107 * wR ** 0.5488

        def dwRexp(wR):
            return 0.02107 * 0.5488 * wR ** (0.5488 - 1)
        lin_cuts.env_wRexp = Block()
        # add_concave_relaxation(lin_cuts.env_wRexp, b.wRexp, b.wR, wRexp, dwRexp, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_concave_relaxation(lin_cuts.env_wRexp, b.wRexp, b.wR, wRexp, dwRexp, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        # # b.wRexp.setlb(0.02107 * b.wR.lb ** 0.5488)
        # # b.wRexp.setub(0.02107 * b.wR.ub ** 0.5488)
        # tighten_vb(b.wRexp,
        #            (0.02107 * lb(b.wR) ** 0.5488,
        #             0.02107 * ub(b.wR) ** 0.5488))
        # tighten_bb(b.wRexp,
        #            (0.02107 * lbb(b.wR) ** 0.5488,
        #             0.02107 * ubb(b.wR) ** 0.5488))
        # add_concave_linear_underest(lin_cuts.env_wRexp, 'wRexp', 1, b.wR, b.wRexp, wRexp, exists=b.equip_exists, block_bounds=self.block_bounds)
        # f_lb = 0.02107 * b.wR.lb ** 0.5488
        # f_ub = 0.02107 * b.wR.ub ** 0.5488
        # df_lb = 0.02107 * 0.5488 * b.wR.lb ** (0.5488 - 1)
        # df_ub = 0.02107 * 0.5488 * b.wR.ub ** (0.5488 - 1)
        # lin_cuts.env_wRexp.overest_lb = Constraint(expr=b.wRexp <= df_lb * (b.wR - b.wR.lb) + f_lb)
        # lin_cuts.env_wRexp.underest_lb = Constraint(expr=b.wRexp <= df_ub * (b.wR - b.wR.ub) + f_ub)

        def FLiqInExp(total_flow_liq_in):
            return total_flow_liq_in ** 0.2279

        def dFLiqInExp(total_flow_liq_in):
            return 0.2279 * total_flow_liq_in ** (0.2279 - 1)
        lin_cuts.env_FLiqInExp = Block()
        # b.FLiqInExp.setlb(b.total_flow_liq_in.lb ** 0.2279)  # disabled for now due to division by zero
        # tighten_vb(b.FLiqInExp,
        #            (None,
        #             ub(b.total_flow_liq_in) ** 0.2279))
        # tighten_bb(b.FLiqInExp,
        #            (None,
        #             ubb(b.total_flow_liq_in) ** 0.2279))
        # add_concave_linear_underest(lin_cuts.env_FLiqInExp, 'FLiqInExp', 1, b.total_flow_liq_in, b.FLiqInExp, FLiqInExp, exists=b.equip_exists, block_bounds=self.block_bounds)
        # # f_lb = b.total_flow_liq_in.lb ** 0.2279
        # f_ub = b.total_flow_liq_in.ub ** 0.2279
        # # df_lb = 0.2279 * b.total_flow_liq_in.lb ** (0.2279 - 1)
        # df_ub = 0.2279 * b.total_flow_liq_in.ub ** (0.2279 - 1)
        # # lin_cuts.env_FLiqInExp.overest_lb = Constraint(expr=b.FLiqInExp <= df_lb * (b.total_flow_liq_in - b.total_flow_liq_in.lb) + f_lb)
        # lin_cuts.env_FLiqInExp.overest_ub = Constraint(expr=b.FLiqInExp <= df_ub * (b.total_flow_liq_in - b.total_flow_liq_in.ub) + f_ub)
        add_concave_relaxation(lin_cuts.env_FLiqInExp, b.FLiqInExp, b.total_flow_liq_in, FLiqInExp, dFLiqInExp, nsegs, None, b.equip_exists, block_bounds=self.block_bounds, bound_contract='monotonic_increase')

        lin_cuts.mc_u_defn = Block()
        tighten_mc(b.u, b.FLiqInExp, b.wRexp)
        # setup_mccormick_cuts(lin_cuts.mc_u_defn, 'u_defn', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_u_defn, 'u_defn', None, b.u, b.FLiqInExp, b.wRexp, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_u_defn, b.u, b.FLiqInExp, b.wRexp, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_t_defn = Block()
        tighten_mc(b.R, b.Ns * b.t, b.u)
        # setup_mccormick_cuts(lin_cuts.mc_t_defn, 't_defn', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_t_defn, 't_defn', None, b.R, b.Ns * b.t, b.u, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_t_defn, b.R, b.Ns * b.t, b.u, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)
        if is_fixed_by_bounds(b.t):
            b.t.fix()
            pass

        def t1(t):
            return t ** 0.5

        def dt1(t):
            return 0.5 * t ** -0.5
        lin_cuts.env_t1 = Block()
        # add_concave_relaxation(lin_cuts.env_t1, b.t1, b.t, t1, dt1, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_concave_relaxation(lin_cuts.env_t1, b.t1, b.t, t1, dt1, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        # f_lb = b.t.lb ** 0.5
        # f_ub = b.t.ub ** 0.5
        # df_lb = 0.5 * b.t.lb ** -0.5
        # df_ub = 0.5 * b.t.ub ** -0.5
        # tighten_vb(b.t1, (lb(b.t) ** 0.5, ub(b.t) ** 0.5))
        # # b.t1.setlb(f_lb)
        # # b.t1.setub(f_ub)
        # add_concave_linear_underest(lin_cuts.env_t1, 't1', 1, b.t, b.t1, t1, exists=b.equip_exists, block_bounds=self.block_bounds)
        # if not is_fixed_by_bounds(b.t):
        #     lin_cuts.env_t1.overest_lb = Constraint(expr=b.t1 <= df_lb * (b.t - b.t.lb) + f_lb)
        #     lin_cuts.env_t1.overest_ub = Constraint(expr=b.t1 <= df_ub * (b.t - b.t.ub) + f_ub)
        #     pass
        # else:
        #     b.t1.fix(f_lb)
        #     pass

        def t2(t):
            return t ** 1.5

        def dt2(t):
            return 1.5 * t ** 0.5
        lin_cuts.env_t2 = Block()
        # we know that this is convex because t > 0
        # add_convex_relaxation(lin_cuts.env_t2, b.t2, b.t, t2, dt2, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_convex_relaxation(lin_cuts.env_t2, b.t2, b.t, t2, dt2, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')

        def t3(t):
            return t ** 2.5

        def dt3(t):
            return 2.5 * t ** 1.5
        lin_cuts.env_t3 = Block()
        # we know that this is convex because t > 0
        # add_convex_relaxation(lin_cuts.env_t3, b.t3, b.t, t3, dt3, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_convex_relaxation(lin_cuts.env_t3, b.t3, b.t, t3, dt3, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')

        def t4(t):
            return t ** 3.5

        def dt4(t):
            return 3.5 * t ** 2.5
        lin_cuts.env_t4 = Block()
        # we know that this is convex because t > 0
        # add_convex_relaxation(lin_cuts.env_t4, b.t4, b.t, t4, dt4, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_convex_relaxation(lin_cuts.env_t4, b.t4, b.t, t4, dt4, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')

        lin_cuts.linear_erf = Constraint(expr=b.erf == 2 / math.pi ** 0.5 * (b.t1 - b.t2 / 3 + b.t3 / 10 - b.t4 / 42))

        lin_cuts.mc_terf = Block()
        # b.terf.setlb(min_lb(b.t, b.erf))
        # b.terf.setub(max_ub(b.t, b.erf))
        tighten_mc(b.terf, b.t, b.erf)
        # setup_mccormick_cuts(lin_cuts.mc_terf, 'terf', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_terf, 'terf', None, b.terf, b.t, b.erf, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_terf, b.terf, b.t, b.erf, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        def expt(t):
            return exp(-b.k_ov * t)

        def dexpt(t):
            return -b.k_ov * exp(-b.k_ov * t)
        lin_cuts.env_expt = Block()
        # add_convex_relaxation(lin_cuts.env_expt, b.expt, b.t, expt, dexpt, b.equip_exists, self.block_bounds, bound_contract='monotonic_decrease')
        add_convex_relaxation(lin_cuts.env_expt, b.expt, b.t, expt, dexpt, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_decrease')

        lin_cuts.mc_t1expt = Block()
        # b.t1expt.setlb(min_lb(b.t1, b.expt))
        # b.t1expt.setub(max_ub(b.t1, b.expt))
        tighten_mc(b.t1expt, b.t1, b.expt)
        # setup_mccormick_cuts(lin_cuts.mc_t1expt, 't1expt', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_t1expt, 't1expt', None, b.t1expt, b.t1, b.expt, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_t1expt, b.t1expt, b.t1, b.expt, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_tkL = Block()
        # b.tkL.setlb(min_lb(b.t, b.kL))
        # b.tkL.setub(max_ub(b.t, b.kL))
        tighten_mc(b.tkL, b.t, b.kL)
        # setup_mccormick_cuts(lin_cuts.mc_tkL, 'tkL', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_tkL, 'tkL', None, b.tkL, b.t, b.kL, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_tkL, b.tkL, b.t, b.kL, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.linear_kL_defn = Constraint(expr=b.tkL == (b.k_ov * b.D) ** 0.5 * (b.terf + 1 / (math.pi * b.k_ov) ** 0.5 * b.t1expt + 1 / (2 * b.k_ov) * b.erf))

        def RoutSq(Rout):
            return Rout ** 2

        def dRoutSq(Rout):
            return 2 * Rout
        lin_cuts.env_RoutSq = Block()
        # add_convex_relaxation(lin_cuts.env_RoutSq, b.RoutSq, b.Rout, RoutSq, dRoutSq, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_convex_relaxation(lin_cuts.env_RoutSq, b.RoutSq, b.Rout, RoutSq, dRoutSq, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        # f_lb = b.Rout.lb ** 2
        # f_ub = b.Rout.ub ** 2
        # df_lb = 2 * b.Rout.lb
        # df_ub = 2 * b.Rout.ub
        # tighten_vb(b.RoutSq, (lb(b.Rout) ** 2, ub(b.Rout) ** 2))
        # # b.RoutSq.setlb(f_lb)
        # # b.RoutSq.setub(f_ub)
        # lin_cuts.env_RoutSq.underest_lb = Constraint(expr=b.RoutSq >= df_lb * (b.Rout - b.Rout.lb) + f_lb)
        # lin_cuts.env_RoutSq.underest_ub = Constraint(expr=b.RoutSq >= df_ub * (b.Rout - b.Rout.ub) + f_ub)
        # lin_cuts.env_RoutSq.overest = Constraint(expr=b.RoutSq <= (f_ub - f_lb) / (b.Rout.ub - b.Rout.lb) * (b.Rout - b.Rout.lb) + f_lb)

        def RinSq(Rin):
            return Rin ** 2

        def dRinSq(Rin):
            return 2 * Rin
        lin_cuts.env_RinSq = Block()
        # add_convex_relaxation(lin_cuts.env_RinSq, b.RinSq, b.Rin, RinSq, dRinSq, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_convex_relaxation(lin_cuts.env_RinSq, b.RinSq, b.Rin, RinSq, dRinSq, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        # f_lb = b.Rin.lb ** 2
        # f_ub = b.Rin.ub ** 2
        # df_lb = 2 * b.Rin.lb
        # df_ub = 2 * b.Rin.ub
        # tighten_vb(b.RinSq, (lb(b.Rin) ** 2, ub(b.Rin) ** 2))
        # # b.RinSq.setlb(f_lb)
        # # b.RinSq.setub(f_ub)
        # lin_cuts.env_RinSq.underest_lb = Constraint(expr=b.RinSq >= df_lb * (b.Rin - b.Rin.lb) + f_lb)
        # lin_cuts.env_RinSq.underest_ub = Constraint(expr=b.RinSq >= df_ub * (b.Rin - b.Rin.ub) + f_ub)
        # lin_cuts.env_RinSq.overest = Constraint(expr=b.RinSq <= (b.RinSq.ub - b.RinSq.lb) / (b.Rin.ub - b.Rin.lb) * (b.Rin - b.Rin.lb) + f_lb)

        lin_cuts.linear_RSqDiff = Constraint(expr=b.RSqDiff == b.RoutSq - b.RinSq)
        # b.RSqDiff.setlb(0)  # Rout >= Rin
        # b.RSqDiff.setub(b.RoutSq.ub - b.RinSq.lb)
        tighten_vb(b.RSqDiff, (0, ub(b.RoutSq) - lb(b.RinSq)))
        tighten_bb(b.RSqDiff, (lbb(b.RoutSq) - ubb(b.RoutSq), ubb(b.RoutSq) - lbb(b.RinSq)))

        lin_cuts.mc_FvapinExpbase = Block()
        tighten_mc(b.FvapinExpbase, b.total_flow_vap_in, b.exp_base)
        # b.FvapinExpbase.setlb(min_lb(b.total_flow_vap_in, b.exp_base))
        # b.FvapinExpbase.setub(max_ub(b.total_flow_vap_in, b.exp_base))
        # setup_mccormick_cuts(lin_cuts.mc_FvapinExpbase, 'FvapinExpbase', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_FvapinExpbase, 'FvapinExpbase', None, b.FvapinExpbase, b.total_flow_vap_in, b.exp_base, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_FvapinExpbase, b.FvapinExpbase, b.total_flow_vap_in, b.exp_base, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_hRSqDiff = Block()
        # b.hRSqDiff.setlb(min_lb(b.h, b.RSqDiff))
        # b.hRSqDiff.setub(max_ub(b.h, b.RSqDiff))
        tighten_mc(b.hRSqDiff, b.h, b.RSqDiff)
        # setup_mccormick_cuts(lin_cuts.mc_hRSqDiff, 'hRSqDiff', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_hRSqDiff, 'hRSqDiff', None, b.hRSqDiff, b.h, b.RSqDiff, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_hRSqDiff, b.hRSqDiff, b.h, b.RSqDiff, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.mc_kLhRSqDiff = Block()
        # b.kLhRSqDiff.setlb(min_lb(b.kL, b.hRSqDiff))
        # b.kLhRSqDiff.setub(max_ub(b.kL, b.hRSqDiff))
        tighten_mc(b.kLhRSqDiff, b.kL, b.hRSqDiff)
        # setup_mccormick_cuts(lin_cuts.mc_kLhRSqDiff, 'kLhRSqDiff', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_kLhRSqDiff, 'kLhRSqDiff', None, b.kLhRSqDiff, b.kL, b.hRSqDiff, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_kLhRSqDiff, b.kLhRSqDiff, b.kL, b.hRSqDiff, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.linear_exp_base_calc = Constraint(expr=b.FvapinExpbase == b.kLhRSqDiff * b.a * math.pi)

        # mass transfer equation
        lin_cuts.mc_frac_vap_inout = Block()
        # b.frac_vap_inout.setlb(min_lb(b.frac_vap_out['A'], 1 - b.frac_vap_in['A']))
        # b.frac_vap_inout.setub(max_ub(b.frac_vap_out['A'], 1 - b.frac_vap_in['A']))
        tighten_mc(b.frac_vap_inout, b.frac_vap_out['A'], 1 - b.frac_vap_in['A'])
        # setup_mccormick_cuts(lin_cuts.mc_frac_vap_inout, 'frac_vap_inout', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_frac_vap_inout, 'frac_vap_inout', None, b.frac_vap_inout, b.frac_vap_out['A'], 1 - b.frac_vap_in['A'], b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_frac_vap_inout, b.frac_vap_inout, b.frac_vap_out['A'], 1 - b.frac_vap_in['A'], nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

        lin_cuts.linear_frac_diff = Constraint(expr=b.frac_diff == b.frac_vap_in['A'] - b.frac_vap_inout)
        tighten_vb(b.frac_diff, (None, ub(b.frac_vap_in['A']) - lb(b.frac_vap_inout)))
        tighten_bb(b.frac_diff, (None, ubb(b.frac_vap_in['A']) - lbb(b.frac_vap_inout)))
        # b.frac_diff.setub(b.frac_vap_in['A'].ub - b.frac_vap_inout.lb)

        def lnFracDiff(frac_diff):
            return log(frac_diff)

        def dlnFracDiff(frac_diff):
            return 1 / frac_diff
        lin_cuts.env_lnFracDiff = Block()
        # add_concave_relaxation(lin_cuts.env_lnFracDiff, b.lnFracDiff, b.frac_diff, lnFracDiff, dlnFracDiff, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_concave_relaxation(lin_cuts.env_lnFracDiff, b.lnFracDiff, b.frac_diff, lnFracDiff, dlnFracDiff, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        # f_lb = log(b.frac_diff.lb)
        # f_ub = log(b.frac_diff.ub)
        # df_lb = 1 / b.frac_diff.lb
        # df_ub = 1 / b.frac_diff.ub
        # tighten_vb(b.lnFracDiff, (log(lb(b.frac_diff)), log(ub(b.frac_diff))))
        # # b.lnFracDiff.setlb(f_lb)
        # # b.lnFracDiff.setub(f_ub)
        # add_concave_linear_underest(lin_cuts.env_lnFracDiff, 'lnFracDiff', 1, b.frac_diff, b.lnFracDiff, lnFracDiff, exists=b.equip_exists, block_bounds=self.block_bounds)
        # lin_cuts.env_lnFracDiff.overest_lb = Constraint(expr=b.lnFracDiff <= df_lb * (b.frac_diff - b.frac_diff.lb) + f_lb)
        # lin_cuts.env_lnFracDiff.overest_ub = Constraint(expr=b.lnFracDiff <= df_ub * (b.frac_diff - b.frac_diff.ub) + f_ub)

        def lnFrac(frac_vap_inout):
            return log(frac_vap_inout)

        def dlnFrac(frac_vap_inout):
            return 1 / frac_vap_inout
        lin_cuts.env_lnFrac = Block()
        # add_concave_relaxation(lin_cuts.env_lnFrac, b.lnFrac, b.frac_vap_inout, lnFrac, dlnFrac, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        add_concave_relaxation(lin_cuts.env_lnFrac, b.lnFrac, b.frac_vap_inout, lnFrac, dlnFrac, nsegs, None, b.equip_exists, self.block_bounds, bound_contract='monotonic_increase')
        # f_lb = log(b.frac_vap_inout.lb)
        # f_ub = log(b.frac_vap_inout.ub)
        # df_lb = 1 / b.frac_vap_inout.lb
        # df_ub = 1 / b.frac_vap_inout.ub
        # tighten_vb(b.lnFrac, (log(lb(b.frac_vap_inout)), log(ub(b.frac_vap_inout))))
        # # b.lnFrac.setlb(f_lb)
        # # b.lnFrac.setub(f_ub)
        # add_concave_linear_underest(lin_cuts.env_lnFrac, 'lnFrac', 1, b.frac_vap_inout, b.lnFrac, lnFrac, exists=b.equip_exists, block_bounds=self.block_bounds)
        # lin_cuts.env_lnFrac.overest_lb = Constraint(expr=b.lnFrac <= df_lb * (b.frac_vap_inout - b.frac_vap_inout.lb) + f_lb)
        # lin_cuts.env_lnFrac.overest_ub = Constraint(expr=b.lnFrac <= df_ub * (b.frac_vap_inout - b.frac_vap_inout.ub) + f_ub)

        lin_cuts.linear_mass_transfer = Constraint(expr=b.lnFracDiff - b.lnFrac == b.exp_base)

        lin_cuts.linear_volume_calc = Constraint(expr=b.V == math.pi * b.hRSqDiff)

        lin_cuts.mc_RH = Block()
        tighten_mc(b.R, b.RH, b.h)
        # setup_mccormick_cuts(lin_cuts.mc_RH, 'RH_calc', nsegs=1)
        # add_mccormick_cut(lin_cuts.mc_RH, 'RH_calc', None, b.R, b.RH, b.h, b.equip_exists, block_bounds=self.block_bounds)
        add_mccormick_relaxation(lin_cuts.mc_RH, b.R, b.RH, b.h, nsegs, None, b.equip_exists, block_bounds=self.block_bounds)

    def get_vars_to_bound(self):
        yield self.total_flow_vap_in, self.block_bounds
        yield self.frac_vap_in, self.block_bounds
        yield self.total_flow_vap_out, self.block_bounds
        yield self.frac_vap_out, self.block_bounds
        yield self.speed
        yield self.wSq
        yield self.h
        yield self.Rout
        yield self.Rin
        yield self.RFourthDiff
        yield self.wSqh
        yield self.R
        yield self.wR
        yield self.total_flow_liq_in, self.block_bounds
        yield self.FLiqInExp
        yield self.wRexp
        yield self.t
        yield self.u
        yield self.t1
        yield self.expt
        yield self.kL
        yield self.exp_base
        yield self.RSqDiff
        yield self.hRSqDiff
        yield self.frac_diff
        yield self.frac_vap_inout
        yield self.RH

    def generate_cut_gen_problem(self):
        raise NotImplementedError()

    def apply_self_proj_cut(self):
        raise NotImplementedError()
