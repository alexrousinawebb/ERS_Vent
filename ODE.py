"""
ORDINARY DIFFERENTIAL EQUATIONS ODE
===========================================================================
Module containing reactor system ODEs for solving.
"""

from Conversion import c2k, A_wet
import Constant_Lib as cc
import numpy as np
import Property_Lib as pl
import VLE
from scipy import integrate as integ
from simple_pid import PID
import ERS
import matplotlib.pyplot as plt

class ODE(ERS.ERS):
    def __init__(self, scen):
        """
            Initializes instance for ordinary differential equation integration routine.

            Arguments:


            Outputs:

        """

        self.scenario = scen
        self.pid_jacket = None
        self.ramp_rate = None
        self.n_vent = 0
        self.critical_flow = False

        tmax = self.scenario.rxn_time * 60 * 60
        self.N = int(tmax)
        self.t = np.linspace(0, tmax, self.N)
        self.xe = 1

        ic = VLE.Initial_Conditions(self.scenario)
        self.H2O = pl.Water(self.scenario.T0, ic.P, ic)
        self.H2O2 = pl.Hydrogen_Peroxide(self.scenario.T0, ic.P, ic)
        self.O2 = pl.Oxygen(self.scenario.T0, ic.P, ic)
        self.cp = pl.Common_Properties(self.scenario.T0, ic, self.H2O, self.H2O2, self.O2)
        self.Y0 = [self.scenario.T0, self.scenario.T0, self.H2O.n, self.H2O2.n, self.O2.n]

        self.data = [np.zeros((self.N, len(self.Y0))), np.zeros((self.N, len(vars(self.H2O)))),
                         np.zeros((self.N, len(vars(self.H2O2)))), np.zeros((self.N, len(vars(self.O2)))),
                         np.zeros((self.N, len(vars(self.cp))))]

        self.store_data(0)

        self.data[0][0][:] = self.Y0

    def initialize(self):
        self.solver = integ.ode(self.rxn_vent_ode)
        self.solver.set_integrator(self.scenario.integrator)
        self.solver.set_initial_value(self.Y0, self.t[0])
        self.pid_config()

    def integrate(self):
        k = 1

        while self.solver.successful() and self.solver.t <= self.t[-1]:

            if self.t[k] / 3600 >= 9:
                self.pid_jacket.setpoint = self.scenario.T0

                if self.cp.P <= cc.Patm:
                    break

            self.ramp_rate = self.pid_jacket(self.data[0][k - 1][0])
            self.solver.set_f_params(k)

            self.solver.integrate(self.t[k])
            self.data[0][k][:] = self.solver.y
            self.store_data(k)

            if self.cp.P >= self.scenario.P_RD:
                break

            if self.scenario.plot_rt is True and k % 60 == 0:
                self.plot_realtime(k)

            k += 1

        plt.show()

    def rxn_vent_ode(self, t, Y, k):

        T, Tj, nH2O, nH2O2, nO2 = Y

        mR = nH2O *cc.MH2O + nH2O2 *cc.MH2O2 + nO2 *cc.MO2 # Total reaction mass (g)

        uc = VLE.Update_Conditions(self.scenario, T, nH2O, nH2O2, nO2, self.data, k)

        self.H2O = pl.Water(T, uc.P, uc)
        self.H2O2 = pl.Hydrogen_Peroxide(T, uc.P, uc)
        self.O2 = pl.Oxygen(T, uc.P, uc)
        self.cp = pl.Common_Properties(T, uc, self.H2O, self.H2O2, self.O2)
        kin = pl.Kinetics(T, self.scenario.kf)

        x = self.cp.nG * (self.H2O.y * cc.MH2O + self.H2O2.y * cc.MH2O2 + self.O2.y * cc.MO2) / \
            (self.cp.nG * (self.H2O.y * cc.MH2O + self.H2O2.y * cc.MH2O2 + self.O2.y * cc.MO2) +
             self.cp.nL * (self.H2O.x * cc.MH2O + self.H2O2.x * cc.MH2O2))  # Mass fraction vapour in vessel

        CpL = self.H2O.x * self.H2O.cpl + self.H2O2.x * self.H2O2.cpl  # Average liquid constant pressure heat capacity (J/(g*K))

        CpG = self.H2O.y * self.H2O.cpg + self.H2O2.y * self.H2O2.cpg + self.O2.y * self.O2.cpg  # Average vapour constant pressure heat capacity (J/(g*K))

        Cp = x * CpG + (1 - x) * CpL  # Average heat capacity in vessel (J/(g*K))

        pG = (1 / (cc.R * c2k(T) * 1000)) * (
                    self.H2O.P * cc.MH2O + self.H2O2.P * cc.MH2O2 + self.O2.P * cc.MO2)  # Average vapour density (kg/L)

        pL = self.H2O.density * self.H2O.x + self.H2O2.density * self.H2O2.x  # Average liquid density (kg/L)

        vL = 1 / pL  # Average specific gravity (L/kg)

        vfg = (1 / pG) - (1 / pL)  # Change in specific volume upon vaporization (L/kg)

        self.vent_calc(k)

        # Differential equations for change in molar amount of components (mol/s)
        dnH2O_dt = kin.rate * nH2O2 - (self.H2O.y * self.xe + self.H2O.x * (1 - self.xe)) * self.n_vent
        dnH2O2_dt = -kin.rate * nH2O2 - (self.H2O2.y * self.xe + self.H2O2.x * (1 - self.xe)) * self.n_vent
        dnO2_dt = kin.rate * nH2O2 / 2 - self.O2.y * self.xe * self.n_vent

        dTj_dt = self.ramp_rate / 60  # Jacket Temperature (C/s)

        Qr = dnH2O2_dt * cc.MH2O2 * cc.dH_rxn  # Reaction heat flux (J/s)

        QHEx = self.scenario.Ux * A_wet(self.cp.VL, self.scenario.D) * (T - Tj)  # External heat exchanger heat flux (J/s)

        Q_vap = ((vL / vfg) + 1) * self.cp.dhvap * self.n_vent * (self.H2O.y + self.H2O2.y) * cc.MH2O  # Latent heat of evaporation (J/s)

        Xp = ((self.cp.dhvap - self.cp.P * vfg) / (vfg * Cp * 1000)) * (x * self.cp.dvGdt / 1000 + (1 - x) *
                                                                        self.cp.dvLdt / 1000)

        dT_dt = (Qr - QHEx - Q_vap) / (mR * Cp * (1 - Xp))

        return [dT_dt, dTj_dt, dnH2O_dt, dnH2O2_dt, dnO2_dt]

    def vdir(self, obj):
        return obj.__dict__.keys()

    def store_data(self, index):
        property_master = [self.H2O, self.H2O2, self.O2, self.cp]

        temp_data = []

        for i in property_master:
            attributes = self.vdir(i)
            value = [getattr(i, j) for j in attributes]
            temp_data.append(value)

        for i in range(len(self.data) - 1):
            self.data[i + 1][index][:] = temp_data[:][i]

    def plot_realtime(self, k):
        plt.figure(1)
        plt.scatter(self.t[k], self.data[0][k][0], color='r')
        plt.scatter(self.t[k], self.data[0][k][1], color='b')

        plt.figure(2)
        plt.scatter(self.t[k], self.data[1][k][9], color='r')
        plt.scatter(self.t[k], self.data[2][k][9], color='b')
        plt.scatter(self.t[k], self.data[1][k][10], color='g')
        plt.scatter(self.t[k], self.data[2][k][10], color='y')
        plt.scatter(self.t[k], self.data[3][k][6], color='k')

        plt.figure(3)
        plt.scatter(self.t[k], self.data[1][k][13], color='r')
        plt.scatter(self.t[k], self.data[2][k][13], color='b')
        # plt.scatter(self.t[k], self.data[3][k][])
        plt.scatter(self.t[k], self.data[4][k][4], color='g')

        plt.pause(0.05)

    def pid_config(self):
        self.pid_jacket = PID(Kp=self.scenario.Kp, Ki=self.scenario.Ki, Kd=self.scenario.Kd,
                              setpoint=self.scenario.rxn_temp,
                              output_limits=(-self.scenario.max_rate, self.scenario.max_rate),
                              auto_mode=True, sample_time=0.01)

    def vent_calc(self, k):
        if self.cp.P > self.scenario.P_BPR:
            self.n_vent = self.ventflow(self.H2O, self.H2O2, self.O2, self.cp,
                                        self.data[0][k - 1][0], self.scenario.P_BPR)
        else:
            self.n_vent = 0
