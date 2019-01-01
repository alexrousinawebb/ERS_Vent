"""
EMERGENCY RELIEF SYSTEM ERS
===========================================================================
Module containing DIERS technology for single phase and two-phase flow.
"""

import Constant_Lib as cc
import numpy as np
from Conversion import c2k, A_relief
from scipy import optimize as opt

class ERS():
    def critical_pressure(self, k, P):
        """
            Calculate the critical pressure for choked flow.
        """
        return (2 / (k + 1)) ** (k / (k - 1)) * P

    def ventflow(self, T, P_discharge):
        """
            Calculates the vent flow rate for all vapour venting scenarios.
        """

        Mw = self.H2O.y * cc.MH2O + self.H2O2.y * cc.MH2O2 + self.O2.y * cc.MO2

        Z = self.H2O.Z * self.H2O.y + self.H2O2.Z * self.H2O2.y + self.O2.Z * self.O2.y

        Kd = self.scenario.BPR_max_Cv / (27.66 * A_relief(self.scenario.D_BPR) * 1550)

        if self.critical_flow is True:

            C = 520 * np.sqrt(self.cp.k * ((2 / (self.cp.k + 1)) ** ((self.cp.k + 1) / (self.cp.k - 1))))

            self.n_vent_vap = ((A_relief(self.scenario.D_BPR) * 1000000 * C * Kd * self.cp.P / 13160) * np.sqrt(
                Mw / (c2k(T) * Z)) * 1000 / 3600) / Mw

        elif self.critical_flow is False:

            r = self.cp.P / P_discharge

            F2 = np.sqrt((self.cp.k / (self.cp.k - 1)) * r ** (2 / self.cp.k) *
                         ((1 - r ** ((self.cp.k - 1) / self.cp.k)) / (1 - r)))

            self.n_vent_vap = (((A_relief(self.scenario.D_BPR) * 1000000 * F2 * Kd / 17.9) * np.sqrt(
                Mw * self.cp.P * (self.cp.P - P_discharge) / (Z * c2k(T)))) * 1000 / 3600) / Mw

    def voidfrac(self, z):

        C0 = 1.5

        F = np.empty(1)

        if self.scenario.flow_regime is 'churn-turbulent':
            F[0] = ((z[0] * ((1 - z[0]) ** 2)) / ((1 - z[0] ** 3) * (1 - C0 * z[0]))) - (self.jgx / self.Ui)
        elif self.scenario.flow_regime is 'bubbly':
            F[0] = ((self.jgx / self.Ui) / (2 + C0 * (self.jgx / self.Ui))) - z[0]

        return F

    def two_phase(self):

        self.jgx = (A_relief(self.scenario.D_RD) * self.n_vent_vap) / (self.pG * A_relief(self.scenario.D * 39.3701) * 1000)

        if self.scenario.flow_regime is 'churn-turbulent':
            Ux_factor = 1.53
        elif self.scenario.flow_regime is 'bubbly':
            Ux_factor = 1.18

        self.Ui = Ux_factor * (self.cp.st * cc.g * 1000 * (self.pL - self.pG)) ** (1 / 4) / np.sqrt(1000 * self.pL)

        alpha = opt.fsolve(self.voidfrac, 0.8)

        alphaves = (self.scenario.VR - self.cp.VL) / self.scenario.VR

        if alpha <= alphaves:
            self.TF = 0
            self.Xm = np.empty(1)
            self.jgi = np.empty(1)
        else:
            self.TF = 1

            C0 = 1.5

            if self.scenario.flow_regime is 'churn-turbulent':
                self.jgi = 2 * alphaves * Ux / (1 - C0 * alphaves)
                a_m = 2 * alphaves / (1 + C0 * alphaves)
            elif self.scenario.flow_regime is 'bubbly':
                self.jgi = alphaves * (1 - alphaves) ** 2 * Ux / ((1 - alphaves ** 3) * (1 - C0 * alphaves))
                a_m = alphaves

            self.Xm = a_m * self.pG / (a_m * self.pG + (1 - a_m) * self.pL)