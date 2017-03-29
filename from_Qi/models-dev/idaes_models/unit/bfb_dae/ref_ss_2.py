# !/usr/env python
from __future__ import print_function
from __future__ import division
from pyomo.environ import *
from cpoinsc import collptsgen
import math
"""bfb_ss - bvp function
Orthogonal collocation (OC)
Model generation function.
written by David M Thierry, July 2016

based on a heavily modified version of FD by Mingzhao Yu; OC implemented
based on ref_ss_1.py
intent: steady state simulation
written for release
note: you may require numpy to run this.
use: given initial guess (100 points).
"""

__author__ = 'David M Thierry'


def bfb_ss(nfe_x, kord_x, _alp_gauB_x, _bet_gauB_x):
    """returns model for a given nfe_x finite elements, and kord_x number of c-points
    if Radau is selected order of approx = 2kord_x - 1
    ir Legendre collocation is selected order = 2kord_x
    recommended: use Legendre collocation
    Recommended use MA57 as linear solver
    Do not use asymmetric discretization with i=1 and j=0 calculated"""

    # length of the bed
    _L = 5.
    # time horizon
    _t = 1.
    # alpha and beta of gauss orthogonal polynomials
    # _alp_gauB_x = 0
    # _bet_gauB_x = 0

    _zi0 = -3.
    _maxh = 0.5
    # tolerance
    _epsi = 1e-04
    _taunc = 0.5
    _hs = 0.

    # differential variables
    # cbin
    # cein
    # ebin
    # ecwin
    # eein
    # hxh
    # P
    # Phx
    # ccwin


    # Lagrange polynomials
    def lgr_x(j, tau, kord, alp, bet):
        tauk = collptsgen(kord, alp, bet)
        tauk.reverse()
        tauk.append(0.)
        tauk.reverse()
        # tauk = [
        #     0.,
        #     0.155051,
        #     0.644949,
        #     1.000000]
        out = 1
        for k in range(0, kord + 1):
            if j != k:
                out *= (tau - tauk[k]) / (tauk[j] - tauk[k])
        return out

    def lgry(j, tau, kord, alp, bet):
        tauk = collptsgen(kord, alp, bet)
        tauk.reverse()
        tauk.append(0.)
        tauk.reverse()
        out = 1
        # for legendre [0, K-1]
        if j == 0:
            return 0
        else:
            for k in range(1, kord + 1):
                if j != k:
                    out *= (tau - tauk[k]) / (tauk[j] - tauk[k])
            return out

    def lgrdot_x(j, tau, kord, alp, bet):
        tauk = collptsgen(kord, alp, bet)
        tauk.reverse()
        tauk.append(0.)
        tauk.reverse()
        out1 = 1
        for k in range(0, kord + 1):
            if k != j:
                out1 *= 1 / (tauk[j] - tauk[k])
        out2 = 1
        out3 = 0
        for m in range(0, kord + 1):
            if m != j:
                out2 = 1  # initialize multiplication
                for n in range(0, kord + 1):
                    if n != m and n != j:
                        out2 *= tau - tauk[n]
                        # elif n == j:
                        # print ("we've got a problem here")
                out3 += out2
        out = out3 * out1

        return out

    def lgrydot(j, tau, kord, alp, bet):
        tauk = collptsgen(kord, alp, bet)
        tauk.reverse()
        tauk.append(0.)
        tauk.reverse()
        out1 = 1
        for k in range(1, kord + 1):
            if k != j:
                out1 *= 1 / (tauk[j] - tauk[k])
        out2 = 1
        out3 = 0
        for m in range(1, kord + 1):
            if m != j:
                out2 = 1  # initialize multiplication
                for n in range(1, kord + 1):
                    if n != m and n != j:
                        out2 *= tau - tauk[n]
                        # elif n == j:
                        # print ("we've got a problem here")
                out3 += out2
        out = out3 * out1

        return out

    # mod
    mod = AbstractModel()

    mod.nfe_x = Param(initialize=nfe_x)
    mod.ncp_x = Param(initialize=kord_x)


    tau_x = collptsgen(kord_x, _alp_gauB_x, _bet_gauB_x)

    # start at zero
    tau_i_x = {0: 0.}
    # create a list
    for ii in range(1, kord_x + 1):
        tau_i_x[ii] = tau_x[ii - 1]

    print('taui', tau_i_x)

    # ======= SETS fe cp ======= #
    # For finite element = 1 .. nfe_x

    fe_x_list = [ii for ii in range(1, nfe_x + 1)]
    cp_x_list = [ii for ii in range(0, kord_x + 1)]

    mod.fe_x = Set(initialize=fe_x_list)
    mod.fe_x_i = Set(initialize=fe_x_list)
    mod.cp_x_i = Set(initialize=cp_x_list)

    # collocation points

    mod.cp_x = Set(initialize=cp_x_list)

    # components
    mod.sp = Set(initialize=['c', 'h', 'n'])

    # create collocation param
    mod.taucp_x = Param(mod.cp_x, initialize=tau_i_x)


    # ldot def
    def __ldoti_x(m, j, k):
        return lgrdot_x(j, m.taucp_x[k], kord_x, _alp_gauB_x, _bet_gauB_x)

    def __ldotyi(m, j, k):
        if j > 0:
            return lgrydot(j, m.taucp_x[k], kord_x, _alp_gauB_x, _bet_gauB_x)
        else:
            return 0.0

    def __lj1_x(m, j):
        return lgr_x(j, 1, kord_x, _alp_gauB_x, _bet_gauB_x)

    def __ljy1(m, j):
        if j > 0:
            return lgry(j, 1, kord_x, _alp_gauB_x, _bet_gauB_x)
        else:
            return 0.0

    mod.ldot_x = Param(mod.cp_x, mod.cp_x, initialize=__ldoti_x)
    mod.lydot = Param(mod.cp_x, mod.cp_x, initialize=__ldotyi)
    mod.l1_x = Param(mod.cp_x, initialize=__lj1_x)
    mod.l1y = Param(mod.cp_x, initialize=__ljy1)
    #

    mod.A1 = Param(initialize=56200.0)
    mod.A2 = Param(initialize=2.62)
    mod.A3 = Param(initialize=98.9)
    mod.ah = Param(initialize=0.8)
    mod.Ao = Param(initialize=4e-4)
    mod.ap = Param(initialize=90.49773755656109)
    mod.Ax = Param(initialize=63.32160110335595)
    mod.cpg_mol = Param(initialize=30.0912)

    _cpgcsb = {'c': 38.4916, 'h': 33.6967}
    _cpgcst = {'c': 38.4913, 'h': 33.6968}
    _cpgcsc = {'c': 39.4617, 'h': 29.1616}
    _cpgcgc = {'c': 39.2805, 'h': 33.7988, 'n': 29.1578}
    _cpgcge = {'c': 39.3473, 'h': 33.8074, 'n': 29.1592}
    _cpgcse = {'c': 39.3486, 'h': 33.8076}

    mod.cpgcsb = Param(mod.sp, initialize=_cpgcsb)
    mod.cpgcst = Param(mod.sp, initialize=_cpgcst)
    mod.cpgcgc = Param(mod.sp, initialize=_cpgcgc)
    mod.cpgcge = Param(mod.sp, initialize=_cpgcge)
    mod.cpgcsc = Param(mod.sp, initialize=_cpgcsc)
    mod.cpgcse = Param(mod.sp, initialize=_cpgcse)
    mod.cps = Param(initialize=1.13)
    mod.Cr = Param(initialize=1.0)

    mod.dH1 = Param(initialize=-52100)
    mod.dH2 = Param(initialize=-36300)
    mod.dH3 = Param(initialize=-64700)
    mod.dp = Param(initialize=1.5e-4)
    mod.dPhx = Param(initialize=0.01)
    mod.dS1 = Param(initialize=-78.5)
    mod.dS2 = Param(initialize=-88.1)
    mod.dS3 = Param(initialize=-175)
    mod.Dt = Param(initialize=9.0)
    mod.Dte = Param(initialize=2.897869210295575)
    mod.dx = Param(initialize=0.02)
    mod.E1 = Param(initialize=28200)
    mod.E2 = Param(initialize=58200)
    mod.E3 = Param(initialize=57700)
    mod.emf = Param(initialize=0.5)
    mod.fw = Param(initialize=0.2)
    mod.hw = Param(initialize=1.5)
    mod.K_d = Param(initialize=100)
    mod.kg = Param(initialize=0.0280596)
    mod.kp = Param(initialize=1.36)
    mod.Lb = Param(initialize=5.0)
    mod.m1 = Param(initialize=1.17)
    mod.mug = Param(initialize=1.92403e-5)
    mod.nv = Param(initialize=2350.0)
    mod.Nx = Param(initialize=941.0835981537471)
    mod.phis = Param(initialize=1.0)
    mod.Pr = Param(initialize=0.7385161078322956)
    mod.rhohx = Param(initialize=985.9393497324021)
    mod.rhos = Param(initialize=442.0)

    mod.pi = Param(initialize=3.14)
    mod.R = Param(initialize=8.314472)
    mod.gc = Param(initialize=9.81)

    # heat exchanger input condition
    mod.HXIn_P = Param(initialize=1.12)
    mod.HXIn_T = Param(initialize=33)
    mod.HXIn_F = Param(initialize=60000)

    # Gas input condition
    mod.GasIn_T = Param(initialize=40)
    _GasIn_z = {'c': 0.13, 'h': 0.06, 'n': 0.81}
    mod.GasIn_z = Param(mod.sp, initialize=_GasIn_z)
    mod.flue_gas_P = Param(initialize=1.8)

    # Solid input condition
    _nin = {'c': 0.01, 'h':0.7, 'n':0.7}
    mod.nin = Param(mod.sp, initialize=_nin)
    # mol/kg
    mod.SolidIn_T = Param(initialize=50)
    mod.sorbent_P = Param(initialize=1.5)

    # atmosphere pressure
    mod.Out2_P = Param(initialize=1.0)

    # input gas valve
    mod.CV_1 = Param(initialize=5.696665718420114)
    mod.per_opening1 = Param(initialize=85.0)

    # output gas valve
    mod.CV_2 = Param(initialize=12.97878700936543)
    mod.per_opening2 = Param(initialize=50.0)

    # input solid valve
    mod.CV_3 = Param(initialize=17483.58063173724)
    mod.per_opening3 = Param(initialize=50)
    mod.eavg = Param(initialize=0.591951)

    # output solid valve
    mod.CV_4 = Param(initialize=11187.66019532553)
    mod.per_opening4 = Param(initialize=50)

    # ######################################################################################################################
    # By abstract model init.

    mod.HXIn_h_ix = Param()
    mod.GasIn_P_ix = Param()
    mod.GasIn_F_ix = Param()
    mod.GasOut_P_ix = Param()
    mod.GasOut_F_ix = Param()
    mod.GasOut_T_ix = Param()
    mod.GasOut_z_ix = Param(mod.sp)
    mod.SolidIn_Fm_ix = Param()
    mod.SolidOut_Fm_ix = Param()
    mod.SolidIn_P_ix = Param()
    mod.SolidOut_P_ix = Param()
    mod.SolidOut_T_ix = Param()
    mod.SorbOut_F_ix = Param()
    mod.rhog_in_ix = Param()
    mod.rhog_out_ix = Param()
    mod.DownOut_P_ix = Param()
    mod.h_downcomer_ix = Param()

    mod.cb_ix = Param(mod.fe_x, mod.sp)
    mod.cbin_ix = Param(mod.fe_x, mod.sp)
    mod.cc_ix = Param(mod.fe_x, mod.sp)
    mod.ccwin_ix = Param(mod.fe_x, mod.sp)
    mod.ce_ix = Param(mod.fe_x, mod.sp)
    mod.cein_ix = Param(mod.fe_x, mod.sp)
    mod.D_ix = Param(mod.fe_x, mod.sp)
    mod.dl_ix = Param()
    mod.g1_ix = Param()
    mod.Kbc_ix = Param(mod.fe_x, mod.sp)
    mod.Kce_ix = Param(mod.fe_x, mod.sp)
    mod.Kgbulk_ix = Param(mod.fe_x, mod.sp)
    mod.Ksbulk_ix = Param(mod.fe_x, mod.sp)
    mod.nc_ix = Param(mod.fe_x, mod.sp)
    mod.ne_ix = Param(mod.fe_x, mod.sp)
    mod.rgc_ix = Param(mod.fe_x, mod.sp)
    mod.rge_ix = Param(mod.fe_x, mod.sp)
    mod.rsc_ix = Param(mod.fe_x, mod.sp)
    mod.rse_ix = Param(mod.fe_x, mod.sp)
    mod.Sit_ix = Param()
    mod.Sot_ix = Param()
    mod.yb_ix = Param(mod.fe_x, mod.sp)
    mod.yc_ix = Param(mod.fe_x, mod.sp)
    mod.ye_ix = Param(mod.fe_x, mod.sp)

    mod.hsinb_ix = Param()
    mod.hsint_ix = Param()
    mod.vmf_ix = Param()
    mod.db0_ix = Param()
    #
    # /////////////\\\\\\\\\\\\\\\\ #
    lst2 = [ii for ii in range(1, 102)]
    # Initial guess
    mod.fe2 = Set(initialize=lst2)
    mod.ar_ix = Param(mod.fe2)
    mod.cbt_ix = Param(mod.fe2)
    mod.cct_ix = Param(mod.fe2)
    mod.cet_ix = Param(mod.fe2)
    mod.db_ix = Param(mod.fe2)
    mod.dbe_ix = Param(mod.fe2)
    mod.dbm_ix = Param(mod.fe2)
    mod.dbu_ix = Param(mod.fe2)
    mod.delta_ix = Param(mod.fe2)
    mod.dthx_ix = Param(mod.fe2)
    mod.e_ix = Param(mod.fe2)
    mod.ebin_ix = Param(mod.fe2)
    mod.ecwin_ix = Param(mod.fe2)
    mod.ed_ix = Param(mod.fe2)
    mod.eein_ix = Param(mod.fe2)
    mod.fb_ix = Param(mod.fe2)
    mod.fc_ix = Param(mod.fe2)
    mod.fcw_ix = Param(mod.fe2)
    mod.fn_ix = Param(mod.fe2)
    mod.g2_ix = Param(mod.fe2)
    mod.g3_ix = Param(mod.fe2)
    mod.gb_ix = Param(mod.fe2)
    mod.hbc_ix = Param(mod.fe2)
    mod.hce_ix = Param(mod.fe2)
    mod.hd_ix = Param(mod.fe2)
    mod.hgbulk_ix = Param(mod.fe2)
    mod.hl_ix = Param(mod.fe2)
    mod.hp_ix = Param(mod.fe2)
    mod.hsbulk_ix = Param(mod.fe2)
    mod.hsc_ix = Param(mod.fe2)
    mod.hse_ix = Param(mod.fe2)
    mod.ht_ix = Param(mod.fe2)
    mod.hxh_ix = Param(mod.fe2)
    mod.jc_ix = Param(mod.fe2)
    mod.je_ix = Param(mod.fe2)
    mod.k1c_ix = Param(mod.fe2)
    mod.k1e_ix = Param(mod.fe2)
    mod.k2c_ix = Param(mod.fe2)
    mod.k2e_ix = Param(mod.fe2)
    mod.k3c_ix = Param(mod.fe2)
    mod.k3e_ix = Param(mod.fe2)
    mod.kcebs_ix = Param(mod.fe2)
    mod.ke1c_ix = Param(mod.fe2)
    mod.ke1e_ix = Param(mod.fe2)
    mod.ke2c_ix = Param(mod.fe2)
    mod.ke2e_ix = Param(mod.fe2)
    mod.ke3c_ix = Param(mod.fe2)
    mod.ke3e_ix = Param(mod.fe2)
    mod.kpa_ix = Param(mod.fe2)
    mod.nup_ix = Param(mod.fe2)
    mod.p_ix = Param(mod.fe2)
    mod.phx_ix = Param(mod.fe2)
    mod.r1c_ix = Param(mod.fe2)
    mod.r1e_ix = Param(mod.fe2)
    mod.r2c_ix = Param(mod.fe2)
    mod.r2e_ix = Param(mod.fe2)
    mod.r3c_ix = Param(mod.fe2)
    mod.r3e_ix = Param(mod.fe2)
    mod.red_ix = Param(mod.fe2)
    mod.rhog_ix = Param(mod.fe2)
    mod.tau_ix = Param(mod.fe2)
    mod.tgb_ix = Param(mod.fe2)
    mod.tgc_ix = Param(mod.fe2)
    mod.tge_ix = Param(mod.fe2)
    mod.thx_ix = Param(mod.fe2)
    mod.tsc_ix = Param(mod.fe2)
    mod.tse_ix = Param(mod.fe2)
    mod.ttube_ix = Param(mod.fe2)
    mod.vb_ix = Param(mod.fe2)
    mod.vbr_ix = Param(mod.fe2)
    mod.ve_ix = Param(mod.fe2)
    mod.vg_ix = Param(mod.fe2)

    # /////////////\\\\\\\\\\\\\\\\ #
    mod.cb_ix_c = Param(mod.fe2)
    mod.cb_ix_h = Param(mod.fe2)
    mod.cb_ix_n = Param(mod.fe2)
    mod.cbin_ix_c = Param(mod.fe2)
    mod.cbin_ix_h = Param(mod.fe2)
    mod.cbin_ix_n = Param(mod.fe2)
    mod.cc_ix_c = Param(mod.fe2)
    mod.cc_ix_h = Param(mod.fe2)
    mod.cc_ix_n = Param(mod.fe2)
    mod.ccwin_ix_c = Param(mod.fe2)
    mod.ccwin_ix_h = Param(mod.fe2)
    mod.ccwin_ix_n = Param(mod.fe2)
    mod.ce_ix_c = Param(mod.fe2)
    mod.ce_ix_h = Param(mod.fe2)
    mod.ce_ix_n = Param(mod.fe2)
    mod.cein_ix_c = Param(mod.fe2)
    mod.cein_ix_h = Param(mod.fe2)
    mod.cein_ix_n = Param(mod.fe2)
    mod.d_ix_c = Param(mod.fe2)
    mod.d_ix_h = Param(mod.fe2)
    mod.d_ix_n = Param(mod.fe2)
    mod.kbc_ix_c = Param(mod.fe2)
    mod.kbc_ix_h = Param(mod.fe2)
    mod.kbc_ix_n = Param(mod.fe2)
    mod.kce_ix_c = Param(mod.fe2)
    mod.kce_ix_h = Param(mod.fe2)
    mod.kce_ix_n = Param(mod.fe2)
    mod.kgbulk_ix_c = Param(mod.fe2)
    mod.kgbulk_ix_h = Param(mod.fe2)
    mod.kgbulk_ix_n = Param(mod.fe2)
    mod.ksbulk_ix_c = Param(mod.fe2)
    mod.ksbulk_ix_h = Param(mod.fe2)
    mod.ksbulk_ix_n = Param(mod.fe2)
    mod.nc_ix_c = Param(mod.fe2)
    mod.nc_ix_h = Param(mod.fe2)
    mod.nc_ix_n = Param(mod.fe2)
    mod.ne_ix_c = Param(mod.fe2)
    mod.ne_ix_h = Param(mod.fe2)
    mod.ne_ix_n = Param(mod.fe2)
    mod.rsc_ix_c = Param(mod.fe2)
    mod.rsc_ix_h = Param(mod.fe2)
    mod.rsc_ix_n = Param(mod.fe2)
    mod.rse_ix_c = Param(mod.fe2)
    mod.rse_ix_h = Param(mod.fe2)
    mod.rse_ix_n = Param(mod.fe2)
    mod.rgc_ix_c = Param(mod.fe2)
    mod.rgc_ix_h = Param(mod.fe2)
    mod.rgc_ix_n = Param(mod.fe2)
    mod.rge_ix_c = Param(mod.fe2)
    mod.rge_ix_h = Param(mod.fe2)
    mod.rge_ix_n = Param(mod.fe2)
    mod.yb_ix_c = Param(mod.fe2)
    mod.yb_ix_h = Param(mod.fe2)
    mod.yb_ix_n = Param(mod.fe2)
    mod.yc_ix_c = Param(mod.fe2)
    mod.yc_ix_h = Param(mod.fe2)
    mod.yc_ix_n = Param(mod.fe2)
    mod.ye_ix_c = Param(mod.fe2)
    mod.ye_ix_h = Param(mod.fe2)
    mod.ye_ix_n = Param(mod.fe2)
    mod.llast = Param(initialize=5.)
    #
    # variable finite element lenght
    # mod.hi = Var(within=NonNegativeReals0, _L), initialize=(_L / 100))
    mod.lenleft = Param()

    def __ir_hi(m, i):
        reg = m.lenleft / m.nfe_x
        return reg

    mod.hi_x = Param(mod.fe_x, initialize=__ir_hi)

    # mod.lenleft.display()


    def __l_irule(m, j, k):
        # type: (model, fe_x, cp) -> current_x
        h0 = sum(m.hi_x[i] for i in range(1, j))
        # return float(m.hi[j] * taui[k] + m.hi[j] * (float(j) - 1.))
        return float(m.hi_x[j] * tau_i_x[k] + h0)

    # print taui
    # length of the bed
    mod.l = Param(mod.fe_x, mod.cp_x, initialize=__l_irule)


    def __ini_cp(i, y, k, taucp):
        dy = y[i + 1] - y[i]
        if i == 1 and k == 1:
            yx = y[i]
            # yx = dy * taucp[k] + y[i]
        else:
            yx = dy * taucp[k] + y[i]
        return yx

    def __ini_cp_dv(i, y, k, taucp):
        dy = y[i + 1] - y[i]
        yx = dy * taucp[k] + y[i]
        return yx


    def __fdp_x(x, y, l):
        # simple int function
        n = len(y)
        dx = l / (n)
        i1 = math.floor(x / dx + 1)
        i1 = int(i1)
        i2 = i1 + 1
        # print(n, dx, x, i1, i2)
        if i2 > n:
            return y[n]
        else:
            dy = y[i2] - y[i1]
            out = (dy) * x + y[i1]
            return out

    # /////////////\\\\\\\\\\\\\\\\ #


    def __ir_ar(m, ix, k):
        h = __l_irule(m, ix, k)
        # print(h)
        return __fdp_x(h, m.ar_ix, m.lenleft)
        # return __ini_cp(ix, m.ar_ix, k, taui)

    def __ir_cbt(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.cbt_ix, m.lenleft)
        # return __ini_cp(i, m.cbt_ix, k, taui)

    def __ir_cct(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.cct_ix, m.lenleft)
        # return __ini_cp(i, m.cct_ix, k, taui)

    def __ir_cet(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.cet_ix, m.lenleft)
        # return __ini_cp(i, m.cet_ix, k, taui)

    def __ir_db(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.db_ix, m.lenleft)
        # return __ini_cp(i, m.db_ix, k, taui)

    def __ir_dbe(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.dbe_ix, m.lenleft)
        # return __ini_cp(i, m.dbe_ix, k, taui)

    def __ir_dbm(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.dbm_ix, m.lenleft)
        # return __ini_cp(i, m.dbm_ix, k, taui)

    def __ir_dbu(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.dbu_ix, m.lenleft)
        # return __ini_cp(i, m.dbu_ix, k, taui)

    def __ir_delta(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.delta_ix, m.lenleft)
        # return __ini_cp(i, m.delta_ix, k, taui)

    def __ir_dthx(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.dthx_ix, m.lenleft)
        # return __ini_cp(i, m.dthx_ix, k, taui)

    def __ir_e(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.e_ix, m.lenleft)
        # return __ini_cp(i, m.e_ix, k, taui)

    def __ir_ebin(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ebin_ix, m.lenleft)
        # return __ini_cp_dv(i, m.ebin_ix, k, taui)

    def __ir_ecwin(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ecwin_ix, m.lenleft)
        # return __ini_cp_dv(i, m.ecwin_ix, k, taui)

    def __ir_ed(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ed_ix, m.lenleft)
        # return __ini_cp(i, m.ed_ix, k, taui)

    def __ir_eein(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.eein_ix, m.lenleft)
        # return __ini_cp_dv(i, m.eein_ix, k, taui)

    def __ir_fb(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.fb_ix, m.lenleft)
        # return __ini_cp(i, m.fb_ix, k, taui)

    def __ir_fc(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.fc_ix, m.lenleft)
        # return __ini_cp(i, m.fc_ix, k, taui)

    def __ir_fcw(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.fcw_ix, m.lenleft)
        # return __ini_cp(i, m.fcw_ix, k, taui)

    def __ir_fn(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.fn_ix, m.lenleft)
        # return __ini_cp(i, m.fn_ix, k, taui)

    def __ir_g2(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.g2_ix, m.lenleft)
        # return __ini_cp(i, m.g2_ix, k, taui)

    def __ir_g3(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.g3_ix, m.lenleft)
        # return __ini_cp(i, m.g3_ix, k, taui)

    def __ir_gb(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.gb_ix, m.lenleft)
        # return __ini_cp(i, m.gb_ix, k, taui)

    def __ir_hbc(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.hbc_ix, m.lenleft)
        # return __ini_cp(i, m.hbc_ix, k, taui)

    def __ir_hce(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.hce_ix, m.lenleft)
        # return __ini_cp(i, m.hce_ix, k, taui)

    def __ir_hd(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.hd_ix, m.lenleft)
        # return __ini_cp(i, m.hd_ix, k, taui)

    def __ir_hgbulk(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.hgbulk_ix, m.lenleft)
        # return __ini_cp(i, m.hgbulk_ix, k, taui)

    def __ir_hl(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.hl_ix, m.lenleft)
        # return __ini_cp(i, m.hl_ix, k, taui)

    def __ir_hp(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.hp_ix, m.lenleft)
        # return __ini_cp(i, m.hp_ix, k, taui)

    def __ir_hsbulk(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.hsbulk_ix, m.lenleft)
        # return __ini_cp(i, m.hsbulk_ix, k, taui)

    def __ir_hsc(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.hsc_ix, m.lenleft)
        # return __ini_cp(i, m.hsc_ix, k, taui)

    def __ir_hse(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.hse_ix, m.lenleft)
        # return __ini_cp(i, m.hse_ix, k, taui)

    def __ir_ht(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ht_ix, m.lenleft)
        # return __ini_cp(i, m.ht_ix, k, taui)

    def __ir_hxh(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.hxh_ix, m.lenleft)
        # return __ini_cp_dv(i, m.hxh_ix, k, taui)

    def __ir_jc(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.jc_ix, m.lenleft)
        # return __ini_cp_dv(i, m.jc_ix, k, taui)

    def __ir_je(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.je_ix, m.lenleft)
        # return __ini_cp_dv(i, m.je_ix, k, taui)

    def __ir_k1c(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.k1c_ix, m.lenleft)
        # return __ini_cp(i, m.k1c_ix, k, taui)

    def __ir_k1e(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.k1e_ix, m.lenleft)
        # return __ini_cp(i, m.k1e_ix, k, taui)

    def __ir_k2c(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.k2c_ix, m.lenleft)
        # return __ini_cp(i, m.k2c_ix, k, taui)

    def __ir_k2e(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.k2e_ix, m.lenleft)
        # return __ini_cp(i, m.k2e_ix, k, taui)

    def __ir_k3c(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.k3c_ix, m.lenleft)
        # return __ini_cp(i, m.k3c_ix, k, taui)

    def __ir_k3e(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.k3e_ix, m.lenleft)
        # return __ini_cp(i, m.k3e_ix, k, taui)

    def __ir_kcebs(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.kcebs_ix, m.lenleft)
        # return __ini_cp(i, m.kcebs_ix, k, taui)

    def __ir_ke1c(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ke1c_ix, m.lenleft)
        # return __ini_cp(i, m.ke1c_ix, k, taui)

    def __ir_ke1e(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ke1e_ix, m.lenleft)
        # return __ini_cp(i, m.ke1e_ix, k, taui)

    def __ir_ke2c(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ke2c_ix, m.lenleft)
        # return __ini_cp(i, m.ke2c_ix, k, taui)

    def __ir_ke2e(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ke2e_ix, m.lenleft)
        # return __ini_cp(i, m.ke2e_ix, k, taui)

    def __ir_ke3e(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ke3c_ix, m.lenleft)
        # return __ini_cp(i, m.ke3c_ix, k, taui)

    def __ir_ke3c(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ke3e_ix, m.lenleft)
        # return __ini_cp(i, m.ke3e_ix, k, taui)

    def __ir_kpa(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.kpa_ix, m.lenleft)
        # return __ini_cp(i, m.kpa_ix, k, taui)

    def __ir_nup(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.nup_ix, m.lenleft)
        # return __ini_cp(i, m.nup_ix, k, taui)

    def __ir_p(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.p_ix, m.lenleft)
        # return __ini_cp_dv(i, m.p_ix, k, taui)

    def __ir_phx(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.phx_ix, m.lenleft)
        # return __ini_cp_dv(i, m.phx_ix, k, taui)

    def __ir_r1c(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.r1c_ix, m.lenleft)
        # return __ini_cp(i, m.r1c_ix, k, taui)

    def __ir_r1e(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.r1e_ix, m.lenleft)
        # return __ini_cp(i, m.r1e_ix, k, taui)

    def __ir_r2c(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.r2c_ix, m.lenleft)
        # return __ini_cp(i, m.r2c_ix, k, taui)

    def __ir_r2e(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.r2e_ix, m.lenleft)
        # return __ini_cp(i, m.r2e_ix, k, taui)

    def __ir_r3c(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.r3c_ix, m.lenleft)
        # return __ini_cp(i, m.r3c_ix, k, taui)

    def __ir_r3e(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.r3e_ix, m.lenleft)
        # return __ini_cp(i, m.r3e_ix, k, taui)

    def __ir_red(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.red_ix, m.lenleft)
        # return __ini_cp(i, m.red_ix, k, taui)

    def __ir_rhog(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.rhog_ix, m.lenleft)
        # return __ini_cp(i, m.rhog_ix, k, taui)

    def __ir_tau(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.tau_ix, m.lenleft)
        # return __ini_cp(i, m.tau_ix, k, taui)

    def __ir_tgb(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.tgb_ix, m.lenleft)
        # return __ini_cp(i, m.tgb_ix, k, taui)

    def __ir_tgc(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.tgc_ix, m.lenleft)
        # return __ini_cp(i, m.tgc_ix, k, taui)

    def __ir_tge(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.tge_ix, m.lenleft)
        # return __ini_cp(i, m.tge_ix, k, taui)

    def __ir_thx(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.thx_ix, m.lenleft)
        # return __ini_cp(i, m.thx_ix, k, taui)

    def __ir_tsc(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.tsc_ix, m.lenleft)
        # return __ini_cp(i, m.tsc_ix, k, taui)

    def __ir_tse(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.tse_ix, m.lenleft)
        # return __ini_cp(i, m.tse_ix, k, taui)

    def __ir_ttube(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ttube_ix, m.lenleft)
        # return __ini_cp(i, m.ttube_ix, k, taui)

    def __ir_vb(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.vb_ix, m.lenleft)
        # return __ini_cp(i, m.vb_ix, k, taui)

    def __ir_vbr(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.vbr_ix, m.lenleft)
        # return __ini_cp(i, m.vbr_ix, k, taui)

    def __ir_ve(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.ve_ix, m.lenleft)
        # return __ini_cp(i, m.ve_ix, k, taui)

    def __ir_vg(m, ix, k):
        h = __l_irule(m, ix, k)
        return __fdp_x(h, m.vg_ix, m.lenleft)
        # return __ini_cp(i, m.vg_ix, k, taui)

    # /////////////\\\\\\\\\\\\\\\\ #
    # double index variables


    def __ir_cb(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.cb_ix_c, m.lenleft)
            # return __ini_cp(i, m.cb_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.cb_ix_h, m.lenleft)
            # return __ini_cp(i, m.cb_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.cb_ix_n, m.lenleft)
            # return __ini_cp(i, m.cb_ix_n, k, taui)

    def __ir_cbin(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.cbin_ix_c, m.lenleft)
            # return __ini_cp_dv(ix, m.cbin_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.cbin_ix_h, m.lenleft)
            # return __ini_cp_dv(ix, m.cbin_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.cbin_ix_n, m.lenleft)
            # return __ini_cp_dv(ix, m.cbin_ix_n, k, taui)

    def __ir_cc(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.cc_ix_c, m.lenleft)
            # return __ini_cp(ix, m.cc_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.cc_ix_h, m.lenleft)
            # return __ini_cp(ix, m.cc_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.cc_ix_n, m.lenleft)
            # return __ini_cp(ix, m.cc_ix_n, k, taui)


    def __ir_ccwin(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.ccwin_ix_c, m.lenleft)
            # return __ini_cp_dv(ix, m.ccwin_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.ccwin_ix_h, m.lenleft)
            # return __ini_cp_dv(ix, m.ccwin_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.ccwin_ix_n, m.lenleft)
            # return __ini_cp_dv(ix, m.ccwin_ix_n, k, taui)

    # def __if_ccwin_h(m, ix, k):
    # def __if_ccwin_n(m, ix, k):

    def __ir_ce(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.ce_ix_c, m.lenleft)
            # return __ini_cp(ix, m.ce_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.ce_ix_h, m.lenleft)
            # return __ini_cp(ix, m.ce_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.ce_ix_n, m.lenleft)
            # return __ini_cp(ix, m.ce_ix_n, k, taui)

    # def __if_ce_h(m, ix, k):
    # def __if_ce_n(m, ix, k):


    def __ir_cein(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.cein_ix_c, m.lenleft)
            # return __ini_cp_dv(ix, m.cein_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.cein_ix_h, m.lenleft)
            # return __ini_cp_dv(ix, m.cein_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.cein_ix_n, m.lenleft)
            # return __ini_cp_dv(ix, m.cein_ix_n, k, taui)

    def __ir_d(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.d_ix_c, m.lenleft)
            # return __ini_cp(ix, m.d_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.d_ix_h, m.lenleft)
            # return __ini_cp(ix, m.d_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.d_ix_n, m.lenleft)
            # return __ini_cp(ix, m.d_ix_n, k, taui)

    def __ir_kbc(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.kbc_ix_c, m.lenleft)
            # return __ini_cp(ix, m.kbc_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.kbc_ix_h, m.lenleft)
            # return __ini_cp(ix, m.kbc_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.kbc_ix_n, m.lenleft)
            # return __ini_cp(ix, m.kbc_ix_n, k, taui)

    def __ir_kce(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.kce_ix_c, m.lenleft)
            # return __ini_cp(ix, m.kce_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.kce_ix_h, m.lenleft)
            # return __ini_cp(ix, m.kce_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.kce_ix_n, m.lenleft)
            # return __ini_cp(ix, m.kce_ix_n, k, taui)

    def __ir_kgbulk(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.kgbulk_ix_c, m.lenleft)
            # return __ini_cp(ix, m.kgbulk_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.kgbulk_ix_h, m.lenleft)
            # return __ini_cp(ix, m.kgbulk_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.kgbulk_ix_n, m.lenleft)
            # return __ini_cp(ix, m.kgbulk_ix_n, k, taui)

    def __ir_ksbulk(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.ksbulk_ix_c, m.lenleft)
            # return __ini_cp(ix, m.ksbulk_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.ksbulk_ix_h, m.lenleft)
            # return __ini_cp(ix, m.ksbulk_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.ksbulk_ix_n, m.lenleft)
            # return __ini_cp(ix, m.ksbulk_ix_n, k, taui)

    def __ir_nc_c(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.nc_ix_c, m.lenleft)
            # return __ini_cp(ix, m.nc_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.nc_ix_h, m.lenleft)
            # return __ini_cp(ix, m.nc_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.nc_ix_n, m.lenleft)
            # return __ini_cp(ix, m.nc_ix_n, k, taui)

    def __ir_ne_c(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.ne_ix_c, m.lenleft)
            # return __ini_cp(ix, m.ne_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.ne_ix_h, m.lenleft)
            # return __ini_cp(ix, m.ne_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.ne_ix_n, m.lenleft)
            # return __ini_cp(ix, m.ne_ix_n, k, taui)

    def __ir_rsc_c(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.rsc_ix_c, m.lenleft)
            # return __ini_cp(ix, m.rsc_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.rsc_ix_h, m.lenleft)
            # return __ini_cp(ix, m.rsc_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.rsc_ix_n, m.lenleft)
            # return __ini_cp(ix, m.rsc_ix_n, k, taui)

    def __ir_rse_c(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.rse_ix_c, m.lenleft)
            # return __ini_cp(ix, m.rse_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.rse_ix_h, m.lenleft)
            # return __ini_cp(ix, m.rse_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.rse_ix_n, m.lenleft)
            # return __ini_cp(ix, m.rse_ix_n, k, taui)

    def __ir_rgc_c(m, ix, k, j):
        h = __l_irule(m, ix, k)
        if j == 'c':
            return __fdp_x(h, m.rgc_ix_c, m.lenleft)
            # return __ini_cp(ix, m.rgc_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.rgc_ix_h, m.lenleft)
            # return __ini_cp(ix, m.rgc_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.rgc_ix_n, m.lenleft)
            # return __ini_cp(ix, m.rgc_ix_n, k, taui)

    def __ir_rge_c(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.rge_ix_c, m.lenleft)
            # return __ini_cp(ix, m.rge_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.rge_ix_h, m.lenleft)
            # return __ini_cp(ix, m.rge_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.rge_ix_n, m.lenleft)
            # return __ini_cp(ix, m.rge_ix_n, k, taui)

    def __ir_yb_c(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.yb_ix_c, m.lenleft)
            # return __ini_cp(ix, m.yb_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.yb_ix_h, m.lenleft)
            # return __ini_cp(ix, m.yb_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.yb_ix_n, m.lenleft)
            # return __ini_cp(ix, m.yb_ix_n, k, taui)

    def __ir_yc_c(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.yc_ix_c, m.lenleft)
            # return __ini_cp(ix, m.yc_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.yb_ix_h, m.lenleft)
            # return __ini_cp(ix, m.yb_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.yb_ix_n, m.lenleft)
            # return __ini_cp(ix, m.yb_ix_n, k, taui)

    def __ir_ye_c(m, ix, k, j):
        h = __l_irule(m, ix, k)

        if j == 'c':
            return __fdp_x(h, m.ye_ix_c, m.lenleft)
            # return __ini_cp(ix, m.ye_ix_c, k, taui)
        elif j == 'h':
            return __fdp_x(h, m.ye_ix_h, m.lenleft)
            # return __ini_cp(ix, m.ye_ix_h, k, taui)
        elif j == 'n':
            return __fdp_x(h, m.ye_ix_n, m.lenleft)
            # return __ini_cp(ix, m.ye_ix_n, k, taui)

    def gasout_zi_rule(m, i):
        return m.GasOut_z_ix[i]

    mod.HXIn_h = Var(within=Reals, initialize=mod.HXIn_h_ix)
    mod.GasIn_P = Var(within=NonNegativeReals, initialize=mod.GasIn_P_ix)
    mod.GasIn_F = Var(within=NonNegativeReals, initialize=mod.GasIn_F_ix)
    mod.GasOut_P = Var(within=NonNegativeReals, initialize=mod.GasOut_P_ix)
    mod.GasOut_F = Var(within=NonNegativeReals, initialize=mod.GasOut_F_ix)
    mod.GasOut_T = Var(within=NonNegativeReals, initialize=mod.GasOut_T_ix)
    mod.GasOut_z = Var(mod.sp, within=NonNegativeReals, initialize=gasout_zi_rule)
    mod.SolidIn_Fm = Var(within=NonNegativeReals, initialize=mod.SolidIn_Fm_ix)
    # mod.SolidOut_Fm = Var(within=NonNegativeReals0, 10000000), initialize=mod.SolidOut_Fm_ix)
    mod.SolidIn_P = Var(within=NonNegativeReals, initialize=mod.SolidIn_P_ix)
    mod.SolidOut_P = Var(within=NonNegativeReals, initialize=mod.SolidOut_P_ix)
    # mod.SolidOut_T = Var(within=NonNegativeReals, initialize=mod.SolidOut_T_ix)
    # mod.SorbOut_F = Var(within=NonNegativeReals, initialize=mod.SorbOut_F_ix)
    mod.rhog_in = Var(within=NonNegativeReals, initialize=mod.rhog_in_ix)
    mod.rhog_out = Var(within=NonNegativeReals, initialize=mod.rhog_out_ix)
    mod.DownOut_P = Var(within=NonNegativeReals, initialize=mod.DownOut_P_ix)
    mod.h_downcomer = Var(within=NonNegativeReals, initialize=mod.h_downcomer_ix)
    # mod.hsinb = Var(within=Reals-100, initialize=mod.hsinb_ix)
    mod.hsint = Var(within=Reals, initialize=mod.hsint_ix)
    mod.vmf = Var(within=NonNegativeReals, initialize=mod.vmf_ix)
    mod.db0 = Var(within=NonNegativeReals, initialize=mod.db0_ix)
    mod.Sit = Var(within=NonNegativeReals, initialize=mod.Sit_ix)
    mod.Sot = Var(within=NonNegativeReals, initialize=mod.Sot_ix)
    mod.g1 = Var(within=NonNegativeReals, initialize=mod.g1_ix)
    mod.Ar = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ar)
    mod.cbt = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_cbt)
    mod.cct = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_cct)
    mod.cet = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_cet)
    mod.cb = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_cb)
    mod.cbin = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_cbin)
    mod.cc = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_cc)
    mod.ccwin = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_ccwin)
    mod.ce = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_ce)
    mod.cein = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_cein)
    mod.D = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_d)
    mod.db = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_db)
    mod.dbe = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_dbe)
    mod.dbm = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_dbm)
    mod.dbu = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_dbu)
    mod.delta = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_delta)
    mod.dThx = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_dthx)
    mod.e = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_e)
    mod.ebin = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ebin)
    mod.ecwin = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_ecwin)
    mod.ed = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ed)
    mod.eein = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_eein)
    mod.fb = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_fb)
    mod.fc = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_fc)
    mod.fcw = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_fcw)
    mod.fn = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_fn)
    mod.g2 = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_g2)
    mod.g3 = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_g3)
    mod.Gb = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_gb)
    mod.Hbc = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_hbc)
    mod.Hce = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_hce)
    mod.hd = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_hd)
    mod.Hgbulk = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_hgbulk)
    mod.hl = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_hl)
    mod.hp = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_hp)
    mod.Hsbulk = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_hsbulk)
    mod.hsc = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_hsc)
    mod.hse = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_hse)
    mod.ht = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ht)
    mod.hxh = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_hxh)
    mod.Jc = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_jc)
    mod.Je = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_je)
    mod.k1c = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_k1c)
    mod.k1e = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_k1e)
    mod.k2c = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_k2c)
    mod.k2e = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_k2e)
    mod.k3c = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_k3c)
    mod.k3e = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_k3e)
    mod.Kbc = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_kbc)
    mod.Kce = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_kce)
    mod.Kcebs = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_kcebs)
    mod.Ke1c = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ke1c)
    mod.Ke1e = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ke1e)
    mod.Ke2c = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ke2c)
    mod.Ke2e = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ke2e)
    mod.Ke3c = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ke3c)
    mod.Ke3e = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ke3e)
    mod.Kgbulk = Var(mod.fe_x, mod.cp_x, mod.sp, within=Reals, initialize=__ir_kgbulk)
    mod.kpa = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_kpa)
    mod.Ksbulk = Var(mod.fe_x, mod.cp_x, mod.sp, within=Reals, initialize=__ir_ksbulk)
    mod.nc = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_nc_c)
    mod.ne = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_ne_c)
    mod.Nup = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_nup)
    mod.P = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_p)
    mod.Phx = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_phx)
    mod.r1c = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_r1c)
    mod.r1e = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_r1e)
    mod.r2c = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_r2c)
    mod.r2e = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_r2e)
    mod.r3c = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_r3c)
    mod.r3e = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=__ir_r3e)
    mod.Red = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_red)
    mod.rgc = Var(mod.fe_x, mod.cp_x, mod.sp, within=Reals, initialize=__ir_rgc_c)
    mod.rge = Var(mod.fe_x, mod.cp_x, mod.sp, within=Reals, initialize=__ir_rge_c)
    mod.rhog = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_rhog)
    mod.rsc = Var(mod.fe_x, mod.cp_x, mod.sp, within=Reals, initialize=__ir_rsc_c)
    mod.rse = Var(mod.fe_x, mod.cp_x, mod.sp, within=Reals, initialize=__ir_rse_c)
    mod.tau = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_tau)
    mod.Tgb = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_tgb)
    mod.Tgc = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_tgc)
    mod.Tge = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_tge)
    mod.Thx = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_thx)
    mod.Tsc = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_tsc)
    mod.Tse = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_tse)
    mod.Ttube = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ttube)
    mod.vb = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_vb)
    mod.vbr = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_vbr)
    mod.ve = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_ve)
    mod.vg = Var(mod.fe_x, mod.cp_x, within=NonNegativeReals, initialize=__ir_vg)
    mod.yb = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_yb_c)
    mod.yc = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_yc_c)
    mod.ye = Var(mod.fe_x, mod.cp_x, mod.sp, within=NonNegativeReals, initialize=__ir_ye_c)
    mod.z = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=0)

    # ////dv\\\\ #
    mod.dcbin_dx = Var(mod.fe_x, mod.cp_x, mod.sp, within=Reals, initialize=0)
    mod.dcein_dx = Var(mod.fe_x, mod.cp_x, mod.sp, initialize=0)
    mod.debin_dx = Var(mod.fe_x, mod.cp_x, initialize=0)
    mod.decwin_dx = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=0)
    mod.deein_dx = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=0)
    mod.dhxh_dx = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=0)
    mod.dP_dx = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=0)
    mod.dz_dx = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=0)
    mod.dPhx_dx = Var(mod.fe_x, mod.cp_x, within=Reals, initialize=0)
    mod.dccwin_dx = Var(mod.fe_x, mod.cp_x, mod.sp, initialize=0)

    def __icbl(m, j):
        if j == 'c':
            return m.cbin_ix_c[100]
        elif j == 'h':
            return m.cbin_ix_h[100]
        elif j == 'n':
            return m.cbin_ix_n[100]

    def __iceinl(m, j):
        if j == 'c':
            return m.cein_ix_c[100]
        elif j == 'h':
            return m.cein_ix_h[100]
        elif j == 'n':
            return m.cein_ix_n[100]

    def __ebinl(m):
        return m.ebin_ix[100]

    def __ecwinl(m):
        return m.ecwin_ix[100]

    def __eeinl(m):
        return m.eein_ix[100]

    def __hxhl(m):
        return m.hxh_ix[100]

    def __pl(m):
        return m.p_ix[100]

    def __phxl(m):
        return m.phx_ix[100]

    def __ccwinl(m, j):
        if j == 'c':
            return m.ccwin_ix_c[100]
        elif j == 'h':
            return m.ccwin_ix_h[100]
        elif j == 'n':
            return m.ccwin_ix_n[100]

    def __ir_hse_l(m):
        return m.hse_ix[100]

    def __ir_ne_c_l(m, j):
        if j == 'c':
            return m.ne_ix_c[100]
        elif j == 'h':
            return m.ne_ix_h[100]
        elif j == 'n':
            return m.ne_ix_n[100]

    def __gbl(m):
        return m.gb_ix[100]

    def __tgbl(m):
        return m.tgb_ix[100]

    def __ybl(m, j):
        if j == 'c':
            return m.yb_ix_c[100]
        elif j == 'h':
            return m.yb_ix_h[100]
        elif j == 'n':
            return m.yb_ix_n[100]

    # last cp last fe
    mod.cbin_l = Var(mod.sp, within=NonNegativeReals, initialize=__icbl)
    mod.cein_l = Var(mod.sp, within=NonNegativeReals, initialize=__iceinl)
    mod.ebin_l = Var(within=NonNegativeReals, initialize=__ebinl)
    mod.ecwin_l = Var(within=Reals, initialize=__ecwinl)
    mod.eein_l = Var(within=Reals, initialize=__eeinl)
    mod.hxh_l = Var(within=Reals, initialize=__hxhl)
    mod.P_l = Var(within=NonNegativeReals, initialize=__pl)
    mod.Phx_l = Var(within=NonNegativeReals, initialize=__phxl)
    mod.ccwin_l = Var(mod.sp, within=NonNegativeReals, initialize=__ccwinl)
    mod.hse_l = Var(within=Reals, initialize=__ir_hse_l)
    mod.ne_l = Var(mod.sp, within=NonNegativeReals, initialize=__ir_ne_c_l)
    mod.Gb_l = Var(within=NonNegativeReals, initialize=__gbl)
    mod.Tgb_l = Var(within=NonNegativeReals, initialize=__tgbl)
    mod.yb_l = Var(mod.sp, within=NonNegativeReals, initialize=__ybl)
    mod.z_l = Var(within=Reals, initialize=0)
    #

    def a1_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.vg[i, j] * m.Ax * m.cbt[i, j] * 3600 == m.Gb[i, j]
        else:
            return Constraint.Skip

    mod.a1 = Constraint(mod.fe_x, mod.cp_x, rule=a1_rule)

    def a3_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.ebin[i, j] == (m.Gb[i, j] / 3600) * m.cpg_mol * m.Tgb[i, j]
        elif j == 0 and i == 1:
            return m.ebin[i, j] == (m.Gb[i, j] / 3600) * m.cpg_mol * m.Tgb[i, j]
        else:
            return Constraint.Skip

    mod.a3 = Constraint(mod.fe_x, mod.cp_x, rule=a3_rule)

    def a4_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.ecwin[i, j] == m.Jc[i, j] * m.hsc[i, j]
        else:

            return Constraint.Skip

    mod.a4 = Constraint(mod.fe_x, mod.cp_x, rule=a4_rule)

    def a5_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.eein[i, j] == m.Je[i, j] * m.hse[i, j]
        else:
            return Constraint.Skip

    mod.a5 = Constraint(mod.fe_x, mod.cp_x, rule=a5_rule)

    def a7_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return m.cbin[i, j, k] == m.yb[i, j, k] * m.Gb[i, j] / 3600
        elif j == 0 and i == 1:
            return m.cbin[i, j, k] == m.yb[i, j, k] * m.Gb[i, j] / 3600
        else:
            return Constraint.Skip

    mod.a7 = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=a7_rule)

    def a8_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return m.ccwin[i, j, k] == m.Jc[i, j] * m.nc[i, j, k]
        else:
            return Constraint.Skip

    mod.a8 = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=a8_rule)

    def a9_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return m.cein[i, j, k] == m.Je[i, j] * m.ne[i, j, k]
        else:
            return Constraint.Skip

    mod.a9 = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=a9_rule)

    def a11_rule_2(m, i, j):
        if 0 < j <= kord_x:
            return m.z[i, j] == m.Je[i, j] - m.Jc[i, j]
        else:
            return Constraint.Skip

    mod.a11_2 = Constraint(mod.fe_x, mod.cp_x, rule=a11_rule_2)

    # dveqn P
    def a12_rule(m, i, j):
        if 0 < j <= kord_x:
            return (m.dP_dx[i, j]) * 100000 == -m.hi_x[i] * (1 - m.e[i, j]) * m.rhos * m.gc
        else:
            return Constraint.Skip

    mod.a12 = Constraint(mod.fe_x, mod.cp_x, rule=a12_rule)

    # ae for Gb
    def a13_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Gb[i, j] / 3600 == m.vb[i, j] * m.Ax * m.delta[i, j] * m.cbt[i, j]
        else:
            return Constraint.Skip

    mod.a13 = Constraint(mod.fe_x, mod.cp_x, rule=a13_rule)

    # ae for JC
    def a14_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Jc[i, j] == m.fw * m.delta[i, j] * m.rhos * (1 - m.ed[i, j]) * m.vb[i, j]
        else:
            return Constraint.Skip

    mod.a14 = Constraint(mod.fe_x, mod.cp_x, rule=a14_rule)

    # bubble, cw, emulsion gas molefrac, a15, a16, a17
    def a15_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return m.cb[i, j, k] == m.yb[i, j, k] * m.cbt[i, j]
        else:
            return Constraint.Skip

    mod.a15 = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=a15_rule)

    def a16_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return m.cc[i, j, k] == m.yc[i, j, k] * m.cct[i, j]
        else:
            return Constraint.Skip

    mod.a16 = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=a16_rule)

    def a17_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return m.ce[i, j, k] == m.ye[i, j, k] * m.cet[i, j]
        else:
            return Constraint.Skip

    mod.a17 = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=a17_rule)

    # total concentration emulsion, cw, bubble a18, a19, a20
    def a18_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.cet[i, j] == sum(m.ce[i, j, k] for k in m.sp)
        else:
            return Constraint.Skip

    mod.a18 = Constraint(mod.fe_x, mod.cp_x, rule=a18_rule)

    def a19_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.cct[i, j] == sum(m.cc[i, j, k] for k in m.sp)
        else:
            return Constraint.Skip

    mod.a19 = Constraint(mod.fe_x, mod.cp_x, rule=a19_rule)

    def a20_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.cbt[i, j] == sum(m.cb[i, j, k] for k in m.sp)
        else:
            return Constraint.Skip

    mod.a20 = Constraint(mod.fe_x, mod.cp_x, rule=a20_rule)

    # IG equation
    def a21_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.cbt[i, j] == m.P[i, j] * 100 / (8.314 * (m.Tgb[i, j] + 273.16))
        else:
            return Constraint.Skip

    mod.a21 = Constraint(mod.fe_x, mod.cp_x, rule=a21_rule)

    # Diffusivity a22, a23, a24
    def a22_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.D[i, j, 'c'] == (0.1593 - 0.1282 * (m.P[i, j] - 1.4) + 0.001 * (m.Tge[i, j] - 60) + 0.0964 * (
                (m.P[i, j] - 1.4) ** 2) - 0.0006921 * ((m.P[i, j] - 1.4) * (m.Tge[i, j] - 60)) -
                                      3.3532e-06 * (m.Tge[i, j] - 60) ** 2) * m.ye[i, j, 'h'] / (
                                         m.ye[i, j, 'h'] + m.ye[i, j, 'n']) + \
                                     (0.1495 - 0.1204 * (m.P[i, j] - 1.4) + 0.0008896 * (m.Tge[i, j] - 60) + 0.0906 * (
                                         (m.P[i, j] - 1.4) ** 2) -
                                      0.0005857 * (m.P[i, j] - 1.4) * (m.Tge[i, j] - 60) -
                                      3.559e-06 * (m.Tge[i, j] - 60) ** 2) * m.ye[i, j, 'n'] / (
                                         m.ye[i, j, 'h'] + m.ye[i, j, 'n'])
        else:
            return Constraint.Skip

    mod.a22 = Constraint(mod.fe_x, mod.cp_x, rule=a22_rule)

    def a23_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.D[i, j, 'h'] == (0.1593 - 0.1282 * (m.P[i, j] - 1.4) + 0.001 * (m.Tge[i, j] - 60) +
                                      0.0964 * ((m.P[i, j] - 1.4) ** 2) - 0.0006921 * (
                                          (m.P[i, j] - 1.4) * (m.Tge[i, j] - 60)) -
                                      3.3532e-06 * (m.Tge[i, j] - 60) ** 2) * m.ye[i, j, 'c'] / (
                                         m.ye[i, j, 'c'] + m.ye[i, j, 'n']) + \
                                     (0.2165 - 0.1743 * (m.P[i, j] - 1.4) + 0.001377 * (m.Tge[i, j] - 60) + 0.13109 * (
                                         (m.P[i, j] - 1.4) ** 2) -
                                      0.0009115 * (m.P[i, j] - 1.4) * (m.Tge[i, j] - 60) -
                                      4.8394e-06 * (m.Tge[i, j] - 60) ** 2) * m.ye[i, j, 'n'] / (
                                         m.ye[i, j, 'c'] + m.ye[i, j, 'n'])
        else:
            return Constraint.Skip

    mod.a23 = Constraint(mod.fe_x, mod.cp_x, rule=a23_rule)

    def a24_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.D[i, j, 'n'] == (0.1495 - 0.1204 * (m.P[i, j] - 1.4) + 0.0008896 * (m.Tge[i, j] - 60) + 0.0906 * (
                (m.P[i, j] - 1.4) ** 2) -
                                      0.0005857 * (m.P[i, j] - 1.4) * (m.Tge[i, j] - 60) -
                                      3.559e-06 * (m.Tge[i, j] - 60) ** 2) * m.ye[i, j, 'c'] / (
                                         m.ye[i, j, 'h'] + m.ye[i, j, 'c']) + \
                                     (0.2165 - 0.1743 * (m.P[i, j] - 1.4) + 0.001377 * (m.Tge[i, j] - 60) + 0.13109 * (
                                         (m.P[i, j] - 1.4) ** 2) -
                                      0.0009115 * (m.P[i, j] - 1.4) * (m.Tge[i, j] - 60) -
                                      4.8394e-06 * (m.Tge[i, j] - 60) ** 2) * m.ye[i, j, 'h'] / (
                                         m.ye[i, j, 'h'] + m.ye[i, j, 'c'])
        else:
            return Constraint.Skip

    mod.a24 = Constraint(mod.fe_x, mod.cp_x, rule=a24_rule)

    # density
    def a25_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rhog[i, j] == m.P[i, j] * 100 * (
                m.ye[i, j, 'c'] * 44.01 + m.ye[i, j, 'n'] * 28.01 + m.ye[i, j, 'h'] * 18.02) \
                                   / (8.314 * (m.Tge[i, j] + 273.16))
        else:
            return Constraint.Skip

    mod.a25 = Constraint(mod.fe_x, mod.cp_x, rule=a25_rule)

    def a26_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Ar[i, j] == (m.dp ** 3) * m.rhog[i, j] * (m.rhos - m.rhog[i, j]) * m.gc / (m.mug ** 2)
        else:
            return Constraint.Skip

    mod.a26 = Constraint(mod.fe_x, mod.cp_x, rule=a26_rule)

    def a27_rule(m, i, j):
        if 0 < j <= kord_x:
            return (1 - m.e[i, j]) == (1 - m.ed[i, j]) * (1 - m.delta[i, j])
        else:
            return Constraint.Skip

    mod.a27 = Constraint(mod.fe_x, mod.cp_x, rule=a27_rule)

    def a28_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.vbr[i, j] == 0.711 * sqrt(m.gc * m.db[i, j])
        else:
            return Constraint.Skip

    mod.a28 = Constraint(mod.fe_x, mod.cp_x, rule=a28_rule)

    # point [1,0] undefined! skip to [1,1]
    def a29_rule(m):
        return m.db0 == 1.38 * (m.gc ** (-0.2)) * ((m.vg[1, 1] - m.ve[1, 1]) * m.Ao) ** 0.4

    mod.a29 = Constraint(rule=a29_rule)

    def a30_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.dbe[i, j] == (m.Dt / 4) * (-m.g1 + m.g3[i, j]) ** 2
        else:
            return Constraint.Skip

    mod.a30 = Constraint(mod.fe_x, mod.cp_x, rule=a30_rule)

    def a31_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.dbm[i, j] == 2.59 * (m.gc ** (-0.2)) * ((m.vg[i, j] - m.ve[i, j]) * m.Ax) ** 0.4
        else:
            return Constraint.Skip

    mod.a31 = Constraint(mod.fe_x, mod.cp_x, rule=a31_rule)

    def a32_rule(m):
        return m.g1 == 2.56E-2 * sqrt(m.Dt / m.gc) / m.vmf

    mod.a32 = Constraint(rule=a32_rule)

    def a33_rule(m, i, j):
        if 0 < j <= kord_x:
            return 4 * m.g2[i, j] == m.Dt * (m.g1 + m.g3[i, j]) ** 2
        else:
            return Constraint.Skip

    mod.a33 = Constraint(mod.fe_x, mod.cp_x, rule=a33_rule)

    def a34_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.g3[i, j] == sqrt(m.g1 ** 2 + 4 * m.dbm[i, j] / m.Dt)
        else:
            return Constraint.Skip

    mod.a34 = Constraint(mod.fe_x, mod.cp_x, rule=a34_rule)

    # x included?
    def a35_rule(m, i, j):
        if 0 < j <= kord_x:
            return (((sqrt(m.dbu[i, j]) - sqrt(m.dbe[i, j])) / (sqrt(m.db0) - sqrt(m.dbe[i, j]))) ** (
                1 - m.g1 / m.g3[i, j])) * \
                   (((sqrt(m.dbu[i, j]) - sqrt(m.g2[i, j])) / (sqrt(m.db0) - sqrt(m.g2[i, j]))) ** (
                       1 + m.g1 / m.g3[i, j])) == \
                   exp(-0.3 * (m.l[i, j]) / m.Dt)
        else:
            return Constraint.Skip

    mod.a35 = Constraint(mod.fe_x, mod.cp_x, rule=a35_rule)

    def a36_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.fc[i, j] == 3. * (m.vmf / m.emf) / (m.vbr[i, j] - (m.vmf / m.emf))
        else:
            return Constraint.Skip

    mod.a36 = Constraint(mod.fe_x, mod.cp_x, rule=a36_rule)

    def a37_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.fcw[i, j] == m.fc[i, j] + m.fw
        else:
            return Constraint.Skip

    mod.a37 = Constraint(mod.fe_x, mod.cp_x, rule=a37_rule)

    def a38_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return m.Kbc[i, j, k] == 1.32 * 4.5 * (m.vmf / m.db[i, j]) + 5.85 * (
                ((m.D[i, j, k] * 1E-4) ** 0.5) * (m.gc ** 0.25) / (m.db[i, j] ** (5 / 4)))
        else:
            return Constraint.Skip

    mod.a38 = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=a38_rule)

    def a39_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return m.Kce[i, j, k] == 6.77 * sqrt(m.ed[i, j] * (m.D[i, j, k] * 1E-4) * m.vbr[i, j] / (m.db[i, j] ** 3))
        else:
            return Constraint.Skip

    mod.a39 = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=a39_rule)

    def a40_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Kcebs[i, j] == 3 * (1 - m.ed[i, j]) / ((1 - m.delta[i, j]) * m.ed[i, j]) * (
                m.ve[i, j] / m.db[i, j])
        else:
            return Constraint.Skip

    mod.a40 = Constraint(mod.fe_x, mod.cp_x, rule=a40_rule)

    def a41_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Hbc[i, j] == 1.32 * 4.5 * m.vmf * m.cbt[i, j] * m.cpg_mol / m.db[i, j] + \
                                  5.85 * sqrt((m.kg / 1000) * m.cbt[i, j] * m.cpg_mol) * (m.gc ** 0.25) / (
                                      m.db[i, j] ** (5 / 4))
        else:
            return Constraint.Skip

    mod.a41 = Constraint(mod.fe_x, mod.cp_x, rule=a41_rule)

    def a42_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Hce[i, j] == 6.78 * sqrt(
                m.ed[i, j] * m.vb[i, j] * (m.kg / 1000) * m.cct[i, j] * m.cpg_mol / (m.db[i, j] ** 3))
        else:
            return Constraint.Skip

    mod.a42 = Constraint(mod.fe_x, mod.cp_x, rule=a42_rule)

    def a43_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Nup[i, j] == 1000 * m.hp[i, j] * m.dp / m.kg
        else:
            return Constraint.Skip

    mod.a43 = Constraint(mod.fe_x, mod.cp_x, rule=a43_rule)

    def a44_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Red[i, j] == m.ve[i, j] * m.dp * m.rhog[i, j] / m.mug
        else:
            return Constraint.Skip

    mod.a44 = Constraint(mod.fe_x, mod.cp_x, rule=a44_rule)

    def a45_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Nup[i, j] == 0.03 * (m.Red[i, j] ** 1.3)
        else:
            return Constraint.Skip

    mod.a45 = Constraint(mod.fe_x, mod.cp_x, rule=a45_rule)

    def a46_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.kpa[i, j] == (3.58 - 2.5 * m.ed[i, j]) * m.kg * ((m.kp / m.kg) ** (0.46 - 0.46 * m.ed[i, j]))
        else:
            return Constraint.Skip

    mod.a46 = Constraint(mod.fe_x, mod.cp_x, rule=a46_rule)

    def a47_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.fn[i, j] == m.vg[i, j] / m.vmf
        else:
            return Constraint.Skip

    mod.a47 = Constraint(mod.fe_x, mod.cp_x, rule=a47_rule)

    def a48_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.tau[i, j] == 0.44 * ((m.dp * m.gc / ((m.vmf ** 2) * ((m.fn[i, j] - m.ah) ** 2))) ** 0.14) * (
                (m.dp / m.dx) ** 0.225)
        else:
            return Constraint.Skip

    mod.a48 = Constraint(mod.fe_x, mod.cp_x, rule=a48_rule)

    def a49_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.fb[i, j] == 0.33 * (((m.vmf ** 2) * ((m.fn[i, j] - m.ah) ** 2) / (m.dp * m.gc)) ** 0.14)
        else:
            return Constraint.Skip

    mod.a49 = Constraint(mod.fe_x, mod.cp_x, rule=a49_rule)

    def a50_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.hd[i, j] == 2 * sqrt(
                (m.kpa[i, j] / 1000) * m.rhos * m.cps * (1 - m.ed[i, j]) / (m.pi * m.tau[i, j]))
        else:
            return Constraint.Skip

    mod.a50 = Constraint(mod.fe_x, mod.cp_x, rule=a50_rule)

    def a51_rule(m, i, j):
        if 0 < j <= kord_x:
            return 1000 * m.hl[i, j] * m.dp / m.kg == 0.009 * (m.Ar[i, j] ** 0.5) * (m.Pr ** 0.33)
        else:
            return Constraint.Skip

    mod.a51 = Constraint(mod.fe_x, mod.cp_x, rule=a51_rule)

    def a52_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.ht[i, j] == m.fb[i, j] * m.hd[i, j] + (1 - m.fb[i, j]) * m.hl[i, j]
        else:
            return Constraint.Skip

    mod.a52 = Constraint(mod.fe_x, mod.cp_x, rule=a52_rule)

    # pde Pressure in hex
    def a53_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.dPhx_dx[i, j] == m.hi_x[i] * m.dPhx + m.hi_x[i] * m.rhohx * 1E-5
        else:
            return Constraint.Skip

    mod.a53 = Constraint(mod.fe_x, mod.cp_x, rule=a53_rule)

    def a54_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.dThx[i, j] == m.Ttube[i, j] - m.Tse[i, j]
        else:
            return Constraint.Skip

    mod.a54 = Constraint(mod.fe_x, mod.cp_x, rule=a54_rule)

    def a55_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.ht[i, j] * m.dThx[i, j] * m.Cr == m.hw * (m.Thx[i, j] - m.Ttube[i, j])
        else:
            return Constraint.Skip

    mod.a55 = Constraint(mod.fe_x, mod.cp_x, rule=a55_rule)

    def a56_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Thx[i, j] == 33.2104 + 14170.15 * (m.hxh[i, j] + 0.285)
        else:
            return Constraint.Skip

    mod.a56 = Constraint(mod.fe_x, mod.cp_x, rule=a56_rule)

    def a57_rule(m):
        return 10 * 1.75 / (m.phis * m.emf ** 3) * (m.dp * m.vmf * m.rhog[1, 1] / m.mug) ** 2 + \
               10 * 150 * (1 - m.emf) / ((m.phis ** 2) * (m.emf ** 3)) * (m.dp * m.vmf * m.rhog[1, 1] / m.mug) == \
               10 * m.dp ** 3 * m.rhog[1, 1] * (m.rhos - m.rhog[1, 1]) * m.gc / m.mug ** 2

    mod.a57 = Constraint(rule=a57_rule)

    def a58_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.k1c[i, j] == m.A1 * (m.Tsc[i, j] + 273.15) * exp(-m.E1 / (m.R * (m.Tsc[i, j] + 273.15)))
        else:
            return Constraint.Skip

    mod.a58 = Constraint(mod.fe_x, mod.cp_x, rule=a58_rule)

    def a59_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.k2c[i, j] == m.A2 * (m.Tsc[i, j] + 273.15) * exp(-m.E2 / (m.R * (m.Tsc[i, j] + 273.15)))
        else:
            return Constraint.Skip

    mod.a59 = Constraint(mod.fe_x, mod.cp_x, rule=a59_rule)

    def a60_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.k3c[i, j] == m.A3 * (m.Tsc[i, j] + 273.15) * exp(-m.E3 / (m.R * (m.Tsc[i, j] + 273.15)))
        else:
            return Constraint.Skip

    mod.a60 = Constraint(mod.fe_x, mod.cp_x, rule=a60_rule)

    def a61_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.k1e[i, j] == m.A1 * (m.Tse[i, j] + 273.15) * exp(-m.E1 / (m.R * (m.Tse[i, j] + 273.15)))
        else:
            return Constraint.Skip

    mod.a61 = Constraint(mod.fe_x, mod.cp_x, rule=a61_rule)

    def a62_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.k2e[i, j] == m.A2 * (m.Tse[i, j] + 273.15) * exp(-m.E2 / (m.R * (m.Tse[i, j] + 273.15)))
        else:
            return Constraint.Skip

    mod.a62 = Constraint(mod.fe_x, mod.cp_x, rule=a62_rule)

    def a63_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.k3e[i, j] == m.A3 * (m.Tse[i, j] + 273.15) * exp(-m.E3 / (m.R * (m.Tse[i, j] + 273.15)))
        else:
            return Constraint.Skip

    mod.a63 = Constraint(mod.fe_x, mod.cp_x, rule=a63_rule)

    def a64_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Ke1c[i, j] * m.P[i, j] * 1E5 == exp(-m.dH1 / (m.R * (m.Tsc[i, j] + 273.15)) + m.dS1 / m.R)
        else:
            return Constraint.Skip

    mod.a64 = Constraint(mod.fe_x, mod.cp_x, rule=a64_rule)

    def a65_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Ke2c[i, j] * m.P[i, j] * 1E5 == exp(-m.dH2 / (m.R * (m.Tsc[i, j] + 273.15)) + m.dS2 / m.R)
        else:
            return Constraint.Skip

    mod.a65 = Constraint(mod.fe_x, mod.cp_x, rule=a65_rule)

    def a66_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Ke3c[i, j] * m.P[i, j] * 1E5 == exp(-m.dH3 / (m.R * (m.Tsc[i, j] + 273.15)) + m.dS3 / m.R)
        else:
            return Constraint.Skip

    mod.a66 = Constraint(mod.fe_x, mod.cp_x, rule=a66_rule)

    def a67_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Ke1e[i, j] * m.P[i, j] * 1E5 == exp(-m.dH1 / (m.R * (m.Tse[i, j] + 273.15)) + m.dS1 / m.R)
        else:
            return Constraint.Skip

    mod.a67 = Constraint(mod.fe_x, mod.cp_x, rule=a67_rule)

    def a68_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Ke2e[i, j] * m.P[i, j] * 1E5 == exp(-m.dH2 / (m.R * (m.Tse[i, j] + 273.15)) + m.dS2 / m.R)
        else:
            return Constraint.Skip

    mod.a68 = Constraint(mod.fe_x, mod.cp_x, rule=a68_rule)

    def a69_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Ke3e[i, j] * m.P[i, j] * 1E5 == exp(-m.dH3 / (m.R * (m.Tse[i, j] + 273.15)) + m.dS3 / m.R)
        else:
            return Constraint.Skip

    mod.a69 = Constraint(mod.fe_x, mod.cp_x, rule=a69_rule)

    def a70_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.r1c[i, j] == m.k1c[i, j] * (
                (m.P[i, j] * m.yc[i, j, 'h'] * 1E5) - (m.nc[i, j, 'h'] * m.rhos / m.Ke1c[i, j]))
        else:
            return Constraint.Skip

    mod.a70 = Constraint(mod.fe_x, mod.cp_x, rule=a70_rule)

    def a71_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.r2c[i, j] == m.k2c[i, j] * ((1 - 2 * (m.nc[i, j, 'n'] * m.rhos / m.nv) -
                                                  (m.nc[i, j, 'c'] * m.rhos / m.nv)) * m.nc[i, j, 'h'] * m.rhos * m.P[
                                                     i, j] *
                                                 m.yc[i, j, 'c'] * 1E5 -
                                                 (
                                                     ((m.nc[i, j, 'n'] * m.rhos / m.nv) + (
                                                         m.nc[i, j, 'c'] * m.rhos / m.nv)) *
                                                     m.nc[
                                                         i, j, 'c'] * m.rhos /
                                                     m.Ke2c[i, j]))
        else:
            return Constraint.Skip

    mod.a71 = Constraint(mod.fe_x, mod.cp_x, rule=a71_rule)

    def a72_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.r3c[i, j] == m.k3c[i, j] * (((1 - 2 * (m.nc[i, j, 'n'] * m.rhos / m.nv) -
                                                   (m.nc[i, j, 'c'] * m.rhos / m.nv)) ** 2) * (
                                                     (m.P[i, j] * m.yc[i, j, 'c'] * 1E5) ** m.m1) -
                                                 ((m.nc[i, j, 'n'] * m.rhos / m.nv) * (
                                                     (m.nc[i, j, 'n'] * m.rhos / m.nv) + (
                                                         m.nc[i, j, 'c'] * m.rhos / m.nv)) /
                                                  m.Ke3c[i, j]))
        else:
            return Constraint.Skip

    mod.a72 = Constraint(mod.fe_x, mod.cp_x, rule=a72_rule)

    def a73_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.r1e[i, j] == m.k1e[i, j] * (
                (m.P[i, j] * m.ye[i, j, 'h'] * 1E5) - (m.ne[i, j, 'h'] * m.rhos / m.Ke1e[i, j]))
        else:
            return Constraint.Skip

    mod.a73 = Constraint(mod.fe_x, mod.cp_x, rule=a73_rule)

    def a74_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.r2e[i, j] == m.k2e[i, j] * ((1. - 2. * (m.ne[i, j, 'n'] * m.rhos / m.nv) -
                                                  (m.ne[i, j, 'c'] * m.rhos / m.nv)) * m.ne[i, j, 'h'] * m.rhos * (
                                                     m.P[i, j] * m.ye[i, j, 'c'] * 1E5) -
                                                 (((m.ne[i, j, 'n'] * m.rhos / m.nv) +
                                                   (m.ne[i, j, 'c'] * m.rhos / m.nv)) * m.ne[i, j, 'c'] * m.rhos /
                                                  m.Ke2e[
                                                      i, j])
                                                 )
        else:
            return Constraint.Skip

    mod.a74 = Constraint(mod.fe_x, mod.cp_x, rule=a74_rule)

    def a75_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.r3e[i, j] == \
                   m.k3e[i, j] * (
                       ((1. - 2. * (m.ne[i, j, 'n'] * m.rhos / m.nv) -
                         (m.ne[i, j, 'c'] * m.rhos / m.nv)) ** 2) * ((m.P[i, j] * m.ye[i, j, 'c'] * 1E5) ** m.m1) -
                       ((m.ne[i, j, 'n'] * m.rhos / m.nv) * (
                           (m.ne[i, j, 'n'] * m.rhos / m.nv) + (m.ne[i, j, 'c'] * m.rhos / m.nv)) / m.Ke3e[i, j]))
        else:
            return Constraint.Skip

    mod.a75 = Constraint(mod.fe_x, mod.cp_x, rule=a75_rule)

    def a76_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rgc[i, j, 'c'] == (m.nv * m.r3c[i, j] + m.r2c[i, j]) / 1000.
        else:
            return Constraint.Skip

    mod.a76 = Constraint(mod.fe_x, mod.cp_x, rule=a76_rule)

    def a77_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rge[i, j, 'c'] == (m.nv * m.r3e[i, j] + m.r2e[i, j]) / 1000.
        else:
            return Constraint.Skip

    mod.a77 = Constraint(mod.fe_x, mod.cp_x, rule=a77_rule)

    def a78_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rsc[i, j, 'c'] == m.r2c[i, j]
        else:
            return Constraint.Skip

    mod.a78 = Constraint(mod.fe_x, mod.cp_x, rule=a78_rule)

    def a79_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rse[i, j, 'c'] == m.r2e[i, j]
        else:
            return Constraint.Skip

    mod.a79 = Constraint(mod.fe_x, mod.cp_x, rule=a79_rule)

    def a80_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rgc[i, j, 'h'] == m.r1c[i, j] / 1000
        else:
            return Constraint.Skip

    mod.a80 = Constraint(mod.fe_x, mod.cp_x, rule=a80_rule)

    def a81_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rge[i, j, 'h'] == m.r1e[i, j] / 1000
        else:
            return Constraint.Skip

    mod.a81 = Constraint(mod.fe_x, mod.cp_x, rule=a81_rule)

    def a82_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rsc[i, j, 'h'] == m.r1c[i, j] - m.r2c[i, j]
        else:
            return Constraint.Skip

    mod.a82 = Constraint(mod.fe_x, mod.cp_x, rule=a82_rule)

    def a83_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rse[i, j, 'h'] == m.r1e[i, j] - m.r2e[i, j]
        else:
            return Constraint.Skip

    mod.a83 = Constraint(mod.fe_x, mod.cp_x, rule=a83_rule)

    def a84_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rgc[i, j, 'n'] == 0
        else:
            return Constraint.Skip

    mod.a84 = Constraint(mod.fe_x, mod.cp_x, rule=a84_rule)

    def a85_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rge[i, j, 'n'] == 0
        else:
            return Constraint.Skip

    mod.a85 = Constraint(mod.fe_x, mod.cp_x, rule=a85_rule)

    def a86_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rsc[i, j, 'n'] == m.nv * m.r3c[i, j]
        else:
            return Constraint.Skip

    mod.a86 = Constraint(mod.fe_x, mod.cp_x, rule=a86_rule)

    def a87_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.rse[i, j, 'n'] == m.nv * m.r3e[i, j]
        else:
            return Constraint.Skip

    mod.a87 = Constraint(mod.fe_x, mod.cp_x, rule=a87_rule)

    def a88_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.hsc[i, j] == ((m.nc[i, j, 'h'] + m.nc[i, j, 'c']) * (m.cpgcsc['h'] * m.Tsc[i, j] + m.dH1) +
                                   m.nc[i, j, 'c'] * (m.cpgcsc['c'] * m.Tsc[i, j] + m.dH2) +
                                   m.nc[i, j, 'n'] * (m.cpgcsc['c'] * m.Tsc[i, j] + m.dH3)) * 1E-3 + m.cps * m.Tsc[i, j]
        else:
            return Constraint.Skip

    mod.a88 = Constraint(mod.fe_x, mod.cp_x, rule=a88_rule)

    def a89_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.hse[i, j] == ((m.ne[i, j, 'h'] + m.ne[i, j, 'c']) * (m.cpgcse['h'] * m.Tse[i, j] + m.dH1) +
                                   m.ne[i, j, 'c'] * (m.cpgcse['c'] * m.Tse[i, j] + m.dH2) +
                                   m.ne[i, j, 'n'] * (m.cpgcse['c'] * m.Tse[i, j] + m.dH3)) * 1E-3 + m.cps * m.Tse[i, j]
        else:
            return Constraint.Skip

    mod.a89 = Constraint(mod.fe_x, mod.cp_x, rule=a89_rule)

    # equation A.1 Gas phase component balance
    def _de_ngb_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return 0 == -m.dcbin_dx[i, j, k] + m.hi_x[i] * (
                -m.Ax * m.delta[i, j] * m.Kbc[i, j, k] * (m.cb[i, j, k] - m.cc[i, j, k])
            ) + m.Kgbulk[i, j, k]
        else:
            return Constraint.Skip

    mod.de_ngb = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=_de_ngb_rule)

    # equation A.2 Gas phase energy balance
    def _de_hgb_rule(m, i, j):
        if 0 < j <= kord_x:
            return 0 == -m.debin_dx[i, j] + \
                        m.hi_x[i] * (-m.Ax * m.delta[i, j] * m.Hbc[i, j] * (m.Tgb[i, j] - m.Tgc[i, j])) + m.Hgbulk[i, j]
        else:
            return Constraint.Skip

    mod.de_hgb = Constraint(mod.fe_x, mod.cp_x, rule=_de_hgb_rule)

    # equation A.3 Gas phase component balance
    def _de_ngc_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return 0 == m.Kbc[i, j, k] * (m.cb[i, j, k] - m.cc[i, j, k]) - m.Kce[i, j, k] * (
                m.cc[i, j, k] - m.ce[i, j, k]) \
                        - m.fcw[i, j] * (1. - m.ed[i, j]) * m.rgc[i, j, k]
        else:
            return Constraint.Skip

    mod.de_ngc = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=_de_ngc_rule)

    # equation A.4 Gas phase energy balance
    def _de_hgc_rule(m, i, j):
        if 0 < j <= kord_x:
            return 0 == m.Hbc[i, j] * (m.Tgb[i, j] - m.Tgc[i, j]) - m.Hce[i, j] * (m.Tgc[i, j] - m.Tge[i, j]) - \
                        m.fcw[i, j] * (1 - m.ed[i, j]) * m.rhos * m.ap * m.hp[i, j] * (m.Tgc[i, j] - m.Tsc[i, j]) - \
                        m.fcw[i, j] * (1 - m.ed[i, j]) * sum(m.rgc[i, j, k] * m.cpgcgc[k] for k in m.sp) * m.Tgc[i, j]
        else:
            return Constraint.Skip

    mod.de_hgc = Constraint(mod.fe_x, mod.cp_x, rule=_de_hgc_rule)

    # equation A.5 Solid phase adsorbed species balance
    def _de_nsc_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return 0 == -m.dccwin_dx[i, j, k] * m.Ax - m.Ksbulk[i, j, k] + \
                        m.hi_x[i] * (-
                                     m.Ax * m.delta[i, j] * m.rhos * m.Kcebs[i, j] * (m.nc[i, j, k] - m.ne[i, j, k]) +
                                     m.Ax * m.fcw[i, j] * m.delta[i, j] * (1 - m.ed[i, j]) * m.rsc[i, j, k])
        else:
            return Constraint.Skip

    mod.de_nsc = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=_de_nsc_rule)

    # equation A.6 Solid phase energy balance
    def _de_hsc_rule(m, i, j):
        if 0 < j <= kord_x:
            return 0 == -m.decwin_dx[i, j] * m.Ax - m.Hsbulk[i, j] + \
                        m.hi_x[i] * (
                            - m.Ax * m.delta[i, j] * m.rhos * m.Kcebs[i, j] * (m.hsc[i, j] - m.hse[i, j]) +
                            m.Ax * m.fcw[i, j] * m.delta[i, j] * (1 - m.ed[i, j]) * sum(
                                (m.rgc[i, j, k] * m.cpgcgc[k]) for k in m.sp) * (m.Tgc[i, j]) +
                            m.Ax * m.fcw[i, j] * m.delta[i, j] * (1 - m.ed[i, j]) * m.rhos * m.ap * m.hp[i, j] * (
                                m.Tgc[i, j] - m.Tsc[i, j]))
        else:
            return Constraint.Skip

    mod.de_hsc = Constraint(mod.fe_x, mod.cp_x, rule=_de_hsc_rule)

    # equation A.7 Gas phase component balance
    def _de_nge_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return 0 == m.hi_x[i] * m.Ax * m.delta[i, j] * m.Kce[i, j, k] * (m.cc[i, j, k] - m.ce[i, j, k]) - \
                        m.hi_x[i] * m.Ax * (1. - m.fcw[i, j] * m.delta[i, j] - m.delta[i, j]) * (1. - m.ed[i, j]) * \
                        m.rge[
                            i, j, k] - \
                        m.Kgbulk[i, j, k]
        else:
            return Constraint.Skip

    mod.de_nge = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=_de_nge_rule)

    # equation A.8 Gas phase energy balance
    def _de_hge_rule(m, i, j):
        if 0 < j <= kord_x:
            return 0 == m.hi_x[i] * m.Ax * m.delta[i, j] * m.Hce[i, j] * (m.Tgc[i, j] - m.Tge[i, j]) - \
                        m.hi_x[i] * m.Ax * (1 - m.fcw[i, j] * m.delta[i, j] - m.delta[i, j]) * (
                            1. - m.ed[i, j]) * m.rhos * m.ap * m.hp[i, j] * (m.Tge[i, j] - m.Tse[i, j]) - \
                        m.Hgbulk[i, j] - \
                        m.hi_x[i] * m.Ax * (1. - m.fcw[i, j] * m.delta[i, j] - m.delta[i, j]) * (1. - m.ed[i, j]) * \
                        sum(m.rge[i, j, k] * m.cpgcge[k] for k in m.sp) * m.Tge[i, j]
        else:
            return Constraint.Skip

    mod.de_hge = Constraint(mod.fe_x, mod.cp_x, rule=_de_hge_rule)

    # equation A.9 Solid phase adsorbed species balance
    def _de_nse_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return 0 == m.dcein_dx[i, j, k] * m.Ax + m.Ksbulk[i, j, k] + \
                        m.hi_x[i] * (
                            m.Ax * m.delta[i, j] * m.rhos * m.Kcebs[i, j] * (m.nc[i, j, k] - m.ne[i, j, k]) +
                            m.Ax * (1 - m.fcw[i, j] * m.delta[i, j] -
                                    m.delta[i, j]) * (1 - m.ed[i, j]) * m.rse[i, j, k])
        else:
            return Constraint.Skip

    mod.de_nse = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=_de_nse_rule)

    # equation A.10 Solid phase energy balance
    def _de_hse_rule(m, i, j):
        if 0 < j <= kord_x:
            return 0 == m.deein_dx[i, j] * m.Ax + m.Hsbulk[i, j] + \
                        m.hi_x[i] * (
                            m.Ax * m.delta[i, j] * m.rhos * m.Kcebs[i, j] * (m.hsc[i, j] - m.hse[i, j]) +
                            m.Ax * (1 - m.fcw[i, j] * m.delta[i, j] - m.delta[i, j]) * (1 - m.ed[i, j]) *
                            sum((m.rge[i, j, k] * m.cpgcge[k]) for k in m.sp) * m.Tge[i, j] +
                            m.Ax * (
                                1. - m.fcw[i, j] * m.delta[i, j] - m.delta[i, j]
                            ) * (1. - m.ed[i, j]) * m.rhos * m.ap * m.hp[i, j] * (m.Tge[i, j] - m.Tse[i, j])
                            + m.pi * m.dx * m.ht[i, j] * m.dThx[i, j] * m.Nx * m.Cr)
        else:
            return Constraint.Skip

    mod.de_hse = Constraint(mod.fe_x, mod.cp_x, rule=_de_hse_rule)

    # Solid flux balance
    def _de_ws_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.dz_dx[i, j] == 0.0
        else:
            return Constraint.Skip

    mod.de_ws = Constraint(mod.fe_x, mod.cp_x, rule=_de_ws_rule)

    def i1_rule(m, i, j, k):
        if 0 < j <= kord_x:
            return m.Kgbulk[i, j, k] == m.K_d * (m.cet[i, j] - m.cbt[i, j]) * m.yb[i, j, k]
        else:
            return Constraint.Skip

    mod.i1 = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=i1_rule)

    def i2_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.Hgbulk[i, j] == m.K_d * (m.cet[i, j] - m.cbt[i, j]) * m.cpg_mol * m.Tgb[i, j]
        else:
            return Constraint.Skip

    mod.i2 = Constraint(mod.fe_x, mod.cp_x, rule=i2_rule)

    # order - 1 approximation of derivative of Jc
    def i3_rule(m, i, k, c):
        if 0 < k <= kord_x:
            return m.Ksbulk[i, k, c] == \
                   -m.Ax * sum(m.lydot[j, k] * m.Jc[i, j] for j in m.cp_x if 0 < j <= kord_x) * m.ne[i, k, c]
        else:
            return Constraint.Skip

    mod.i3 = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=i3_rule)

    # order - 1 approximation of derivative of Jc
    def i4_rule(m, i, k):
        if 0 < k <= kord_x:
            return m.Hsbulk[i, k] == \
                   -m.Ax * sum(m.lydot[j, k] * m.Jc[i, j] for j in m.cp_x if 0 < j <= kord_x) * m.hse[i, k]
        else:
            return Constraint.Skip

    mod.i4 = Constraint(mod.fe_x, mod.cp_x, rule=i4_rule)

    def i5_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.db[i, j] == m.dbu[i, j]
        else:
            return Constraint.Skip

    mod.i5 = Constraint(mod.fe_x, mod.cp_x, rule=i5_rule)

    def i6_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.vb[i, j] == \
                   1.55 * ((m.vg[i, j] - m.vmf) + 14.1 * (m.db[i, j] + 0.005)) * (m.Dte ** 0.32) + m.vbr[i, j]
        else:
            return Constraint.Skip

    mod.i6 = Constraint(mod.fe_x, mod.cp_x, rule=i6_rule)

    def i7_rule(m, i, j):
        if 0 < j <= kord_x:
            return (1 - m.emf) * (
                (m.dp ** 0.1) * (m.gc ** 0.118) * 2.05 * ((m.l[i, j]) ** 0.043)) == \
                   2.54 * (m.mug ** 0.066) * (1. - m.ed[i, j])
        else:
            return Constraint.Skip

    mod.i7 = Constraint(mod.fe_x, mod.cp_x, rule=i7_rule)

    def i8_rule(m, i, j):
        if 0 < j <= kord_x:
            return m.ve[i, j] * ((m.dp ** 0.568) * (m.gc ** 0.663) * (0.08518 * (m.rhos - m.rhog[i, j]) + 19.09) *
                                 ((m.l[i, j]) ** 0.244)) == \
                   m.vmf * 188. * 1.02 * (m.mug ** 0.371)
        else:
            return Constraint.Skip

    mod.i8 = Constraint(mod.fe_x, mod.cp_x, rule=i8_rule)

    def e1_rule(m):
        return m.HXIn_h == -0.2831 - 2.9863e-6 * (m.HXIn_P - 1.3) + 7.3855e-05 * (m.HXIn_T - 60)

    mod.e1 = Constraint(rule=e1_rule)

    # Heat-Exchanger fluid energy balance
    # pde
    def e3_rule(m, i, j):
        if 0 < j <= kord_x:
            return 0 == (m.HXIn_F / 3600) * m.dhxh_dx[i, j] - \
                        m.hi_x[i] * 1E-6 * m.pi * m.dx * m.ht[i, j] * m.dThx[i, j] * m.Nx * m.Cr
        else:
            return Constraint.Skip

    mod.e3 = Constraint(mod.fe_x, mod.cp_x, rule=e3_rule)

    def e5_rule(m):
        return m.hsint == ((m.nin['h'] + m.nin['c']) * (m.cpgcst['h'] * m.SolidIn_T + m.dH1) +
                           m.nin['c'] * (m.cpgcst['c'] * m.SolidIn_T + m.dH2) +
                           m.nin['n'] * (m.cpgcst['c'] * m.SolidIn_T + m.dH3)) * 1E-3 + m.cps * m.SolidIn_T

    mod.e5 = Constraint(rule=e5_rule)

    # given by valves/specified
    # Gas inlet P
    def e7_rule(m):
        return m.GasIn_P == m.P[1, 0] + 0.034

    #
    mod.e7 = Constraint(rule=e7_rule)

    # Gas inlet F
    def e8_rule(m):
        return m.Gb[1, 0] == m.GasIn_F

    mod.e8 = Constraint(rule=e8_rule)

    # T gas b inlet
    def e9_rule(m):
        return m.Tgb[1, 0] == m.GasIn_T

    mod.e9 = Constraint(rule=e9_rule)

    # Gas In mfrac
    def e10_rule(m, k):
        return m.yb[1, 0, k] == m.GasIn_z[k]

    mod.e10 = Constraint(mod.sp, rule=e10_rule)

    # Gas out P
    def x_3_rule(m):
        return m.GasOut_P == m.P_l

    mod.x_3 = Constraint(rule=x_3_rule)

    # Gas out F
    def e12_rule(m):
        return m.Gb_l == m.GasOut_F

    mod.e12 = Constraint(rule=e12_rule)

    # Gas out T
    def e13_rule(m):
        return m.GasOut_T == m.Tgb_l

    mod.e13 = Constraint(rule=e13_rule)

    # Gas out molef
    def e14_rule(m, j):
        return m.GasOut_z[j] == m.yb_l[j]

    mod.e14 = Constraint(mod.sp, rule=e14_rule)

    def e15_rule(m):
        return m.Sit == m.SolidIn_Fm / 3600

    mod.e15 = Constraint(rule=e15_rule)

    def e16_rule(m):
        return m.SolidIn_P == m.GasOut_P

    mod.e16 = Constraint(rule=e16_rule)

    def v1_rule(m):
        return (m.GasIn_F / 3600) == \
               (m.CV_1 * (m.per_opening1 / 100) * ((m.flue_gas_P - m.GasIn_P) / m.rhog_in) ** 0.5)

    mod.v1 = Constraint(rule=v1_rule)

    def v2_rule(m):
        return m.rhog_in == \
               m.GasIn_P * 100 * (m.GasIn_z['c'] * 44.01 + m.GasIn_z['n'] * 28.01 + m.GasIn_z['h'] * 18.02) / (
                   8.314 * (m.GasIn_T + 273.16))

    mod.v2 = Constraint(rule=v2_rule)

    def v4_rule(m):
        return m.GasOut_F / 3600 == m.CV_2 * (m.per_opening2 / 100) * ((m.GasOut_P - m.Out2_P) / m.rhog_out) ** 0.5

    mod.v4 = Constraint(rule=v4_rule)

    def v5_rule(m):
        return m.rhog_out == m.GasOut_P * 100 * (
            m.GasOut_z['c'] * 44.01 + m.GasOut_z['n'] * 28.01 + m.GasOut_z['h'] * 18.02) / \
                             (8.314 * (m.GasOut_T + 273.16))

    mod.v5 = Constraint(rule=v5_rule)

    def v3_rule(m):
        return (m.SolidIn_Fm / 3600) == m.CV_3 * (m.per_opening3 / 100) * (
                                                                              (m.sorbent_P - m.SolidIn_P) / (
                                                                                  2. * m.rhos)) ** 0.5

    mod.v3 = Constraint(rule=v3_rule)

    def __dvar_x_cbin_(m, i, k, c):
        if 0 < k <= kord_x:
            return m.dcbin_dx[i, k, c] == sum(m.ldot_x[j, k] * m.cbin[i, j, c] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __dvar_x_cein_(m, i, k, c):
        if 0 < k <= kord_x:
            return m.dcein_dx[i, k, c] == sum(m.ldot_x[j, k] * m.cein[i, j, c] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __dvar_x_ebin_(m, i, k):
        if 0 < k <= kord_x:
            return m.debin_dx[i, k] == sum(m.ldot_x[j, k] * m.ebin[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __dvar_x_ecwin_(m, i, k):
        if 0 < k <= kord_x:
            return m.decwin_dx[i, k] == sum(m.ldot_x[j, k] * m.ecwin[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __dvar_x_eein_(m, i, k):
        if 0 < k <= kord_x:
            return m.deein_dx[i, k] == sum(m.ldot_x[j, k] * m.eein[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __dvar_x_hxh_(m, i, k):  #
        if 0 < k <= kord_x:
            return m.dhxh_dx[i, k] == sum(m.ldot_x[j, k] * m.hxh[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __dvar_x_p_(m, i, k):
        if 0 < k <= kord_x:
            return m.dP_dx[i, k] == sum(m.ldot_x[j, k] * m.P[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __dvar_x_phx_(m, i, k):  #
        if 0 < k <= kord_x:
            return m.dPhx_dx[i, k] == sum(m.ldot_x[j, k] * m.Phx[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __dvar_x_ccwin_(m, i, k, c):
        if 0 < k <= kord_x:
            return m.dccwin_dx[i, k, c] == sum(m.ldot_x[j, k] * m.ccwin[i, j, c] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __dvar_z_(m, i, k):
        if 0 < k <= kord_x:
            return m.dz_dx[i, k] == sum(m.ldot_x[j, k] * m.z[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    mod.dvar_c_z = Constraint(mod.fe_x, mod.cp_x, rule=__dvar_z_)
    mod.dvar_x_cbin = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=__dvar_x_cbin_)
    mod.dvar_x_cein = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=__dvar_x_cein_)
    mod.dvar_x_ebin = Constraint(mod.fe_x, mod.cp_x, rule=__dvar_x_ebin_)
    mod.dvar_x_ecwin = Constraint(mod.fe_x, mod.cp_x, rule=__dvar_x_ecwin_)
    mod.dvar_x_eein = Constraint(mod.fe_x, mod.cp_x, rule=__dvar_x_eein_)
    mod.dvar_x_hxh = Constraint(mod.fe_x, mod.cp_x, rule=__dvar_x_hxh_)
    mod.dvar_x_p = Constraint(mod.fe_x, mod.cp_x, rule=__dvar_x_p_)
    mod.dvar_x_phx = Constraint(mod.fe_x, mod.cp_x, rule=__dvar_x_phx_)
    mod.dvar_x_ccwin = Constraint(mod.fe_x, mod.cp_x, mod.sp, rule=__dvar_x_ccwin_)

    # continuity of fe eqns
    def __cp_x_cbin(m, i, c):
        if i < nfe_x:
            return m.cbin[i + 1, 0, c] == sum(m.l1_x[j] * m.cbin[i, j, c] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __cp_x_cein(m, i, c):
        if i < nfe_x:
            return m.cein[i + 1, 0, c] == sum(m.l1_x[j] * m.cein[i, j, c] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __cp_x_ebin(m, i):
        if i < nfe_x:
            return m.ebin[i + 1, 0] == sum(m.l1_x[j] * m.ebin[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __cp_x_ecwin(m, i):
        if i < nfe_x:
            return m.ecwin[i + 1, 0] == sum(m.l1_x[j] * m.ecwin[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __cp_x_eein(m, i):
        if i < nfe_x:
            return m.eein[i + 1, 0] == sum(m.l1_x[j] * m.eein[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __cp_x_hxh(m, i):
        if i < nfe_x:
            return m.hxh[i + 1, 0] == sum(m.l1_x[j] * m.hxh[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __cp_x_p(m, i):
        if i < nfe_x:
            return m.P[i + 1, 0] == sum(m.l1_x[j] * m.P[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __cp_x_phx(m, i):
        if i < nfe_x:
            return m.Phx[i + 1, 0] == sum(m.l1_x[j] * m.Phx[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    def __cp_x_ccwin(m, i, c):
        if i < nfe_x:
            return m.ccwin[i + 1, 0, c] == sum(m.l1_x[j] * m.ccwin[i, j, c] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    mod.cp1_c = Constraint(mod.fe_x, mod.sp, rule=__cp_x_cbin)
    mod.cp2_c = Constraint(mod.fe_x, mod.sp, rule=__cp_x_cein)
    mod.cp3_c = Constraint(mod.fe_x, rule=__cp_x_ebin)
    mod.cp4_c = Constraint(mod.fe_x, rule=__cp_x_ecwin)
    mod.cp5_c = Constraint(mod.fe_x, rule=__cp_x_eein)
    mod.cp6_c = Constraint(mod.fe_x, rule=__cp_x_hxh)
    mod.cp8_c = Constraint(mod.fe_x, rule=__cp_x_p)
    mod.cp9_c = Constraint(mod.fe_x, rule=__cp_x_phx)
    mod.cp10_c = Constraint(mod.fe_x, mod.sp, rule=__cp_x_ccwin)

    def __cp_z_(m, i):
        if i < nfe_x:
            return m.z[i + 1, 0] == sum(m.l1_x[j] * m.z[i, j] for j in m.cp_x if j <= kord_x)
        else:
            return Constraint.Skip

    mod.cpz_c = Constraint(mod.fe_x, rule=__cp_z_)

    # dv last fe last cp
    def __zl_cbin_(m, c):
        return m.cbin_l[c] == sum(m.l1_x[j] * m.cbin[nfe_x, j, c] for j in m.cp_x if j <= kord_x)

    def __zl_cein_(m, c):
        return m.cein_l[c] == sum(m.l1_x[j] * m.cein[nfe_x, j, c] for j in m.cp_x if j <= kord_x)

    def __zl_ebin_(m):
        return m.ebin_l == sum(m.l1_x[j] * m.ebin[nfe_x, j] for j in m.cp_x if j <= kord_x)

    def __zl_ecwin_(m):
        return m.ecwin_l == sum(m.l1_x[j] * m.ecwin[nfe_x, j] for j in m.cp_x if j <= kord_x)

    def __zl_eein_(m):
        return m.eein_l == sum(m.l1_x[j] * m.eein[nfe_x, j] for j in m.cp_x if j <= kord_x)

    def __zl_hxh_(m):
        return m.hxh_l == sum(m.l1_x[j] * m.hxh[nfe_x, j] for j in m.cp_x if j <= kord_x)

    def __zl_p_(m):
        return m.P_l == sum(m.l1_x[j] * m.P[nfe_x, j] for j in m.cp_x if j <= kord_x)

    def __zl_phx_(m):
        return m.Phx_l == sum(m.l1_x[j] * m.Phx[nfe_x, j] for j in m.cp_x if j <= kord_x)

    def __zl_ccwin_(m, c):
        return m.ccwin_l[c] == sum(m.l1_x[j] * m.ccwin[nfe_x, j, c] for j in m.cp_x if j <= kord_x)

    mod.zl_cbin = Constraint(mod.sp, rule=__zl_cbin_)
    mod.zl_cein = Constraint(mod.sp, rule=__zl_cein_)
    mod.zl_ebin = Constraint(rule=__zl_ebin_)
    mod.zl_ecwin = Constraint(rule=__zl_ecwin_)
    mod.zl_eein = Constraint(rule=__zl_eein_)
    mod.zl_hxh = Constraint(rule=__zl_hxh_)  #
    mod.zl_p = Constraint(rule=__zl_p_)
    mod.zl_phx = Constraint(rule=__zl_phx_)  #
    mod.zl_ccwin = Constraint(mod.sp, rule=__zl_ccwin_)

    def __zl_z_(m):
        return m.z_l == sum(m.l1_x[j] * m.z[nfe_x, j] for j in m.cp_x if j <= kord_x)

    mod.zl_z = Constraint(rule=__zl_z_)  #

    # ae last fe last cp
    def __yl_hse(m):
        return m.hse_l == sum(m.l1y[j] * m.hse[nfe_x, j] for j in m.cp_x if 0 < j <= kord_x)

    def __yl_ne(m, c):
        return m.ne_l[c] == sum(m.l1y[j] * m.ne[nfe_x, j, c] for j in m.cp_x if 0 < j <= kord_x)

    def __yl_gb(m):
        return m.Gb_l == sum(m.cbin_l[c] for c in m.sp) * 3600

    def __yl_tgb(m):
        return m.ebin_l == (m.Gb_l / 3600) * m.cpg_mol * m.Tgb_l

    def __yl_yb(m, c):
        return m.cbin_l[c] == m.yb_l[c] * m.Gb_l / 3600

    mod.yl_hse = Constraint(rule=__yl_hse)
    mod.yl_ne = Constraint(mod.sp, rule=__yl_ne)
    mod.yl_gb = Constraint(rule=__yl_gb)
    mod.yl_tgb = Constraint(rule=__yl_tgb)
    mod.yl_yb = Constraint(mod.sp, rule=__yl_yb)

    # some bcs
    def _bc_ccs_b(m, c):
        return m.ccwin[1, 0, c] == m.cein[1, 0, c]

    mod.bc_ccs_b = Constraint(mod.sp, rule=_bc_ccs_b)

    def _bc_ecs_b(m):
        return m.ecwin[1, 0] == m.eein[1, 0]

    mod.bc_ecs_b = Constraint(rule=_bc_ecs_b)

    def _bc_z_b(m):
        return m.z[1, 0] == 0.0

    mod.bc_z_b = Constraint(rule=_bc_z_b)

    def _sot_o_(m):
        return m.Sit == m.Sot + m.z_l * m.Ax

    mod.sot_o_ = Constraint(rule=_sot_o_)

    def _bc_ccs_t(m, j):
        return m.ccwin_l[j] * m.Ax + m.Sit * m.nin[j] == \
               m.cein_l[j] * m.Ax + m.Sot * m.ne_l[j]

    mod.bc_ccs_t = Constraint(mod.sp, rule=_bc_ccs_t)

    def _bc_ecs_t(m):
        return m.ecwin_l * m.Ax + m.Sit * m.hsint == \
               m.eein_l * m.Ax + m.Sot * m.hse_l

    mod.bc_ecs_t = Constraint(rule=_bc_ecs_t)

    # bc_hxh
    def _bc_hxin_(m):
        return 0 == m.HXIn_h - m.hxh_l

    mod.bc_hxin_ = Constraint(rule=_bc_hxin_)

    # bc_Phx
    def _bc_Phx_(m):
        return m.Phx_l == m.HXIn_P

    mod.bc_Phx_ = Constraint(rule=_bc_Phx_)

    mod.objFunc = Objective(expr=1, sense=minimize)

    return mod
