"""
ORDINARY DIFFERENTIAL EQUATIONS ODE
===========================================================================
Module containing reactor system ODEs for solving.
"""

from Conversion import c2k
import Constant_Lib as cc
import numpy as np
import Property_Lib as pl
import VLE
from scipy import integrate as integ
from simple_pid import PID
from pprint import pprint

class ODE():
    def __init__(self, scen):
        """
            Initializes instance for ordinary differential equation integration routine.

            Arguments:


            Outputs:

        """

        self.scenario = scen
        self.pid_jacket = None

        tmax = self.scenario.rxn_time * 60 * 60
        self.N = int(tmax)
        self.t = np.linspace(0, tmax, self.N)

        ic = VLE.Initial_Conditions(self.scenario)
        H2O = pl.Water(self.scenario.T0, ic.P, ic)
        H2O2 = pl.Hydrogen_Peroxide(self.scenario.T0, ic.P, ic)
        O2 = pl.Oxygen(self.scenario.T0, ic.P, ic)
        cp = pl.Common_Properties(self.scenario.T0, ic, H2O, H2O2, O2)

        self.Y0 = [self.scenario.T0, self.scenario.T0, H2O.n, H2O2.n, O2.n]

        self.data = [np.zeros((self.N, len(self.Y0))), np.zeros((self.N, len(vars(H2O)))),
                         np.zeros((self.N, len(vars(H2O2)))), np.zeros((self.N, len(vars(O2)))),
                         np.zeros((self.N, len(vars(cp))))]

        def vdir(obj):
            return obj.__dict__.keys()

        property_master = [H2O, H2O2, O2, cp]

        temp_data = []

        for i in property_master:
            attributes = vdir(i)
            value = [getattr(i, j) for j in attributes]
            temp_data.append(value)

        for i in range(len(self.data) - 1):
            self.data[i+1][0][:] = temp_data[:][i]

    def initialize(self):
        self.solver = integ.ode(self.BPR_PID)
        self.solver.set_integrator(self.scenario.integrator)
        self.solver.set_initial_value(self.Y0, self.t[0])

        self.pid_jacket = PID(Kp = 0.016, Ki = 0, Kd = 0, setpoint = self.scenario.rxn_temp,
                         output_limits = (-self.scenario.max_rate, self.scenario.max_rate),
                              auto_mode = True, sample_time = 0.01)

    def integrate(self):
        k = 1

        while self.solver.successful() and self.solver.t <= self.t[-1]:

            if self.t[k] / 3600 >= 9:
                self.pid_jacket.setpoint = self.scenario.T0

                if self.data[5][k - 1][] <= cc.Patm:
                    break

            rate[k] = pid_jacket(sol[k - 1, 0])

            solver.set_f_params(Urx, Dv, rate[k], Cv_BPR_max,
                                list(data[k - 1, i] for i in [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]))
            solver.integrate(t[k])
            sol[k] = solver.y
            data[k] = checkoutput(solver.y, list(data[k - 1, i] for i in [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]))
            dsol[k] = BPR_PID(t[k], sol[k, :], Urx, Dv, rate[k], Cv_BPR_max,
                              list(data[k, i] for i in [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]))

            P_critical[k] = P_crit(data[k, 8], sol[k, 0], data[k, 3], data[k, 4], data[k, 5])

            if P_discharge <= P_critical[k]:
                crit = 1
            else:
                crit = 0

            if data[k, 8] > P_BPR:
                n_vent[k] = ventflow(sol[k, 0], data[k, 8], P_BPR, data[k, 3], data[k, 4], data[k, 5], Cv_BPR_max,
                                     D_BPR, 0)
            else:
                n_vent[k] = 0

            if data[k, 8] >= P_RD:
                break

            k += 1

    def BPR_PID(self, t, Y):

        T, Tj, nH2O, nH2O2, nO2 = Y

        mR = nH2O *cc.MH2O + nH2O2 *cc.MH2O2 + nO2 *cc.MO2 # Total reaction mass (g)

        uc = VLE.Update_Conditions(self.scenario, T, nH2O, nH2O2, nO2, data)

        H2O = pl.Water(T, uc.P, uc)
        H2O2 = pl.Hydrogen_Peroxide(T, uc.P, uc)
        O2 = pl.Oxygen(T, uc.P, uc)
        cp = pl.Common_Properties(T, uc, H2O, H2O2, O2)
        kin = pl.Kinetics(T, self.scenario.kf)

        x = cp.nG *(H2O.y *cc.MH2O + H2O2.y *cc.MH2O2 + O2.y *cc.MO2 ) /(cp.nG *(H2O.y *cc.MH2O + H2O2.y *cc.MH2O2 + O2.y *cc.MO2) + cp.nL *
                    (H2O.x *cc.MH2O + H2O2.x *cc.MH2O2)) # Mass fraction vapour in vessel

        CpL = H2O.x *H2O.cpl + H2O2.x *H2O2.cpl # Average liquid constant pressure heat capacity (J/(g*K))

        CpG = H2O.y *H2O.cpg + H2O2.y *H2O2.cpg + O2.y *O2.cpg # Average vapour constant pressure heat capacity (J/(g*K))

        Cp = x* CpG + (1 - x) * CpL  # Average heat capacity in vessel (J/(g*K))

        pG = (1 / (cc.R * c2k(T) * 1000)) * (H2O.P * cc.MH2O + H2O2.P * cc.MH2O2 + O2.P * cc.MO2)  # Average vapour density (kg/L)

        pL = H2O.density * H2O.x + H2O2.density * H2O2.x  # Average liquid density (kg/L)

        vL = 1 / pL  # Average specific gravity (L/kg)

        vfg = (1 / pG) - (1 / pL)  # Change in specific volume upon vaporization (L/kg)

        if cc.Patm <= P_crit(P, T, yH2O, yH2O2, yO2):
            flag = 1
        else:
            flag = 0

        if cp.P > self.scenario.P_BPR:
            n_vent = ventflow(T, P, P_BPR, yH2O, yH2O2, yO2, Cv_BPR_max, D_BPR, 0)
        else:
            n_vent = 0

        # Differential equations for change in molar amount of components (mol/s)
        dnH2O_dt = kin.rate * nH2O2 - n_vent * H2O.y
        dnH2O2_dt = -kin.rate * nH2O2 - n_vent * H2O2.y
        dnO2_dt = kin.rate * nH2O2 / 2 - n_vent * O2.y

        dTj_dt = self.scenario.max_rate / 60  # Jacket Temperature (C/s)

        Qr = dnH2O2_dt * cc.MH2O2 * cc.dH_rxn  # Reaction heat flux (J/s)

        QHEx = self.scenario.Ux * Awet(VL, Dv) * (T - Tj)  # External heat exchanger heat flux (J/s)

        Q_vap = ((vL / vfg) + 1) * cp.dhvap * n_vent * (H2O.y + H2O2.y) * cc.MH2O  # Latent heat of evaporation (J/s)

        Xp = ((cp.dhvap - cp.P * vfg) / (vfg * Cp * 1000)) * (x * cp.dvGdt / 1000 + (1 - x) * cp.dvLdt / 1000)

        dT_dt = (Qr - QHEx - Q_vap) / (mR * Cp * (1 - Xp))

        return [dT_dt, dTj_dt, dnH2O_dt, dnH2O2_dt, dnO2_dt]