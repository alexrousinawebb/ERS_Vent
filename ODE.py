"""
ORDINARY DIFFERENTIAL EQUATIONS ODE
===========================================================================
Module containing reactor system ODEs for solving.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate as integ
from simple_pid import PID
from pprint import pprint as pprint

import Constant_Lib as cc
import ERS
import Property_Lib as pl
import VLE
from Conversion import c2k, A_wet


class ODE(ERS.ERS):
    def __init__(self, scen):
        """
            Initializes instance for ordinary differential equation integration routine.

            Arguments:
            scen:   Object containing all relevant parameters from the user generated reaction scenario
        """

        #  Initialize common parameters
        self.scenario = scen
        self.pid_jacket = None
        self.ramp_rate = None
        self.n_vent_vap = None
        self.n_vent_TF = None
        self.venting = None
        self.n_vent = 0
        self.critical_flow = False
        self.xe = 1
        self.CpG = None
        self.CpL = None
        self.Cp = None
        self.pL = None
        self.pG = None
        self.vG = None
        self.vL = None
        self.jgx = None
        self.Ui = None

        #  Generate timespan mesh for integration
        tmax = (self.scenario.rxn_time + self.scenario.cool_time) * 60 * 60
        self.N = int(tmax)
        self.t = np.linspace(0, tmax, self.N)

        #  Initialize starting reaction conditions and compound objects
        ic = VLE.Initial_Conditions(self.scenario)
        self.H2O = pl.Water(self.scenario.T0, ic.P, ic)
        self.H2O2 = pl.Hydrogen_Peroxide(self.scenario.T0, ic.P, ic)
        self.O2 = pl.Oxygen(self.scenario.T0, ic.P, ic)
        self.cp = pl.Common_Properties(self.scenario.T0, ic, self.H2O, self.H2O2, self.O2)
        self.Y0 = [self.scenario.T0, self.scenario.T0, self.H2O.n, self.H2O2.n, self.O2.n]

        # Generate data storage list of arrays
        self.data = [np.zeros((self.N, len(self.Y0))), np.zeros((self.N, len(vars(self.H2O)))),
                         np.zeros((self.N, len(vars(self.H2O2)))), np.zeros((self.N, len(vars(self.O2)))),
                         np.zeros((self.N, len(vars(self.cp))))]

        #  Write initial conditions to data storage array
        self.store_data(0)
        self.data[0][0][:] = self.Y0

    def initialize(self):
        """
            Initializes ode solver for integration.
        """
        self.solver = integ.ode(self.rxn_vent_ode)
        self.solver.set_integrator(self.scenario.integrator)
        self.solver.set_initial_value(self.Y0, self.t[0])
        self.pid_config()

    def integrate(self):
        k = 1

        while self.solver.successful() and self.solver.t < self.t[-1] and self.cp.P <= self.scenario.MAWP:

            self.venting = False

            if self.t[k] / 3600 >= self.scenario.rxn_time:
                self.pid_jacket.setpoint = self.scenario.T0

            if self.cp.P >= self.scenario.P_RD:
                break

            self.ramp_rate = self.pid_jacket(self.data[0][k - 1][0])
            self.solver.set_f_params(k)

            self.solver.integrate(self.t[k])
            self.data[0][k][:] = self.solver.y
            self.store_data(k)

            if self.scenario.plot_rt is True and k % 60 == 0:
                self.plot_realtime(k)

            k += 1

        self.t = self.t[:k]
        self.data = [i[:k, :] for i in self.data]

    def rxn_vent_ode(self, t, Y, k):
        """
            Main ODE funtion for integration of both venting, nonventing, and BPR enabled scenarios.

            Arguments:
                k: Integration iteration number. Used for indexing current data set.
        """

        T, Tj, nH2O, nH2O2, nO2 = Y

        # Total reaction mass (g)
        mR = nH2O *cc.MH2O + nH2O2 *cc.MH2O2 + nO2 *cc.MO2

        #  Update VLE conditions in reactor system and generate compound objects/kinetic data
        uc = VLE.Update_Conditions(self.scenario, T, nH2O, nH2O2, nO2, self.data, k)
        self.H2O = pl.Water(T, uc.P, uc)
        self.H2O2 = pl.Hydrogen_Peroxide(T, uc.P, uc)
        self.O2 = pl.Oxygen(T, uc.P, uc)
        self.cp = pl.Common_Properties(T, uc, self.H2O, self.H2O2, self.O2)
        kin = pl.Kinetics(T, self.scenario.kf)

        #  Mass fraction vapour in vessel
        x = self.cp.nG * (self.H2O.y * cc.MH2O + self.H2O2.y * cc.MH2O2 + self.O2.y * cc.MO2) / \
            (self.cp.nG * (self.H2O.y * cc.MH2O + self.H2O2.y * cc.MH2O2 + self.O2.y * cc.MO2) +
             self.cp.nL * (self.H2O.x * cc.MH2O + self.H2O2.x * cc.MH2O2))

        #  Average phase and vessel properties
        self.CpL = self.H2O.x * self.H2O.cpl + self.H2O2.x * self.H2O2.cpl  # Average liquid constant pressure heat capacity (J/(g*K))

        self.CpG = self.H2O.y * self.H2O.cpg + self.H2O2.y * self.H2O2.cpg + self.O2.y * self.O2.cpg  # Average vapour constant pressure heat capacity (J/(g*K))

        self.Cp = x * self.CpG + (1 - x) * self.CpL  # Average heat capacity in vessel (J/(g*K))

        self.pG = (1 / (cc.R * c2k(T) * 1000)) * (
                    self.H2O.P * cc.MH2O + self.H2O2.P * cc.MH2O2 + self.O2.P * cc.MO2)  # Average vapour density (kg/L)

        self.pL = self.H2O.density * self.H2O.x + self.H2O2.density * self.H2O2.x  # Average liquid density (kg/L)

        self.vL = 1 / self.pL  # Average specific gravity (L/kg)

        self.vG = 1 / self.pG

        vfg = self.vG - self.vL  # Change in specific volume upon vaporization (L/kg)

        #  Calculate vent flow rate
        self.vent_calc(k)

        #  Differential equations for change in molar amount of components (mol/s)
        dnH2O_dt = kin.rate * nH2O2 - (self.H2O.y * self.xe + self.H2O.x * (1 - self.xe)) * self.n_vent
        dnH2O2_dt = -kin.rate * nH2O2 - (self.H2O2.y * self.xe + self.H2O2.x * (1 - self.xe)) * self.n_vent
        dnO2_dt = kin.rate * nH2O2 / 2 - self.O2.y * self.xe * self.n_vent

        #  Vessel thermodynamics
        dTj_dt = self.ramp_rate / 60

        Q_rxn = dnH2O2_dt * cc.MH2O2 * cc.dH_rxn

        Q_HEx = self.scenario.Ux * A_wet(self.cp.VL, self.scenario.D) * (T - Tj)

        Q_vap = ((self.vL / vfg) + 1) * self.cp.dhvap * self.n_vent * (self.H2O.y + self.H2O2.y) * cc.MH2O

        Xp = ((self.cp.dhvap - self.cp.P * vfg) / (vfg * self.Cp * 1000)) * (x * self.cp.dvGdt / 1000 + (1 - x) *
                                                                        self.cp.dvLdt / 1000)

        dT_dt = (Q_rxn - Q_HEx - Q_vap) / (mR * self.Cp * (1 - Xp))

        return [dT_dt, dTj_dt, dnH2O_dt, dnH2O2_dt, dnO2_dt]

    def vdir(self, obj):
        """
            Parse non-callable and non-dunderfunction object attributes
        """
        return obj.__dict__.keys()

    def store_data(self, index):
        """
            Store integration data into previously initialized data array.
        """
        property_master = [self.H2O, self.H2O2, self.O2, self.cp]

        temp_data = []

        for i in property_master:
            attributes = self.vdir(i)
            value = [getattr(i, j) for j in attributes]
            temp_data.append(value)

        for i in range(len(self.data) - 1):
            self.data[i + 1][index][:] = temp_data[:][i]

    def plot_realtime(self, k):
        """
            Plot integration parameters in real time during integration routine.
        """
        plt.figure(1, figsize=(15,10))

        plt.subplot(2, 2, 1)
        plt.scatter(self.t[k], self.data[0][k][0], color='r', s=1)
        plt.scatter(self.t[k], self.data[0][k][1], color='b', s=1)
        plt.xlabel('Time (s)')
        plt.ylabel('Temperature (C)')
        plt.title("Temperature Profile")
        plt.legend(['Reactor', 'Jacket'])

        plt.subplot(2, 2, 2)
        plt.scatter(self.t[k], self.data[1][k][9]*100, color='r', s=1)
        plt.scatter(self.t[k], self.data[2][k][9]*100, color='b', s=1)
        plt.scatter(self.t[k], self.data[1][k][10]*100, color='g', s=1)
        plt.scatter(self.t[k], self.data[2][k][10]*100, color='y', s=1)
        plt.scatter(self.t[k], self.data[3][k][6]*100, color='k', s=1)
        plt.xlabel('Time (s)')
        plt.ylabel('Mole Fraction (%)')
        plt.title("Reactor Composition")
        plt.legend(['Water (l)', 'Hydrogen Peroxide (l)', 'Water (v)', 'Hydrogen Peroxide (v)', 'Oxygen (v)'])

        plt.subplot(2, 2, 3)
        plt.scatter(self.t[k], self.data[1][k][13], color='r', s=1)
        plt.scatter(self.t[k], self.data[2][k][13], color='b', s=1)
        plt.scatter(self.t[k], self.data[3][k][9], color='k', s=1)
        plt.scatter(self.t[k], self.data[4][k][4], color='g', s=1)
        plt.xlabel('Time (s)')
        plt.ylabel('Pressure (kPa)')
        plt.title("Reactor Pressure")
        plt.legend(['Water', 'Hydrogen Peroxide', 'Oxygen', 'Total'])

        plt.pause(0.05)

    def plot_vals(self):
        """
            Plot integration parameters in real time during integration routine.
        """
        plt.figure(2, figsize=(15,10))

        plt.subplot(2, 2, 1)
        plt.plot(self.t[:], [i[0] for i in self.data[0]], color='r')
        plt.plot(self.t[:], [i[1] for i in self.data[0]], color='b')
        plt.xlabel('Time (s)')
        plt.ylabel('Temperature (C)')
        plt.title("Temperature Profile")
        plt.legend(['Reactor', 'Jacket'])

        plt.subplot(2, 2, 2)
        plt.plot(self.t[:], [i[9] for i in self.data[1]], color='r')
        plt.plot(self.t[:], [i[9] for i in self.data[2]], color='b')
        plt.plot(self.t[:], [i[10] for i in self.data[1]], color='g')
        plt.plot(self.t[:], [i[10] for i in self.data[2]], color='y')
        plt.plot(self.t[:], [i[6] for i in self.data[3]], color='k')
        plt.xlabel('Time (s)')
        plt.ylabel('Mole Fraction (%)')
        plt.title("Reactor Composition")
        plt.legend(['Water (l)', 'Hydrogen Peroxide (l)', 'Water (v)', 'Hydrogen Peroxide (v)', 'Oxygen (v)'])

        plt.subplot(2, 2, 3)
        plt.plot(self.t[:], [i[13] for i in self.data[1]], color='r')
        plt.plot(self.t[:], [i[13] for i in self.data[2]], color='b')
        plt.plot(self.t[:], [i[9] for i in self.data[3]], color='k')
        plt.plot(self.t[:], [i[4] for i in self.data[4]], color='g')
        plt.xlabel('Time (s)')
        plt.ylabel('Pressure (kPa)')
        plt.title("Reactor Pressure")
        plt.legend(['Water', 'Hydrogen Peroxide', 'Oxygen', 'Total'])

        plt.show()

    def pid_config(self):
        """
            Configure jacket proportional-integral-derivative controller
        """
        self.pid_jacket = PID(Kp=self.scenario.Kp, Ki=self.scenario.Ki, Kd=self.scenario.Kd,
                              setpoint=self.scenario.rxn_temp,
                              output_limits=(-self.scenario.max_rate, self.scenario.max_rate),
                              auto_mode=True, sample_time=0.01)

    def vent_calc(self, k):
        """
            Calculate reactor vent (BPR/RD/PRV) flow for various scenarios.
        """
        if self.scenario.P_RD > self.cp.P > self.scenario.P_BPR and self.scenario.BPR == True:
            self.ventflow(self.data[0][k - 1][0], self.scenario.P_BPR)
            self.n_vent = self.n_vent_vap
            self.xe = 1

        elif self.cp.P >= self.scenario.P_RD and self.scenario.PRV == True and self.venting == True:
            self.ventflow(self.data[0][k - 1][0], cc.Patm)

            if self.scenario.TF_vent is True:
                self.two_phase()

                if self.TF == False:
                    self.n_vent = self.n_vent_vap
                    self.xe = 1
                elif self.TF == True:
                    self.flow_twophase(self.data[0][k - 1][0],cc.Patm)

                if self.xe > 1:
                    self.xe = 1
                elif self.xe < 0:
                    self.xe = 0
            elif self.scenario.TF_vent is False:
                self.n_vent = self.n_vent_vap
                self.xe = 1

        else:
            self.n_vent = 0
            self.xe = 1