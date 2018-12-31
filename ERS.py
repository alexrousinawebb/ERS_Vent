"""
EMERGENCY RELIEF SYSTEM ERS
===========================================================================
Module containing DIERS technology for single phase and two-phase flow.
"""

import Constant_Lib as cc
import numpy as np
from Conversion import c2k, A_relief

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

            W = ((A_relief(self.scenario.D_BPR) * 1000000 * C * Kd * self.cp.P / 13160) * np.sqrt(
                Mw / (c2k(T) * Z)) * 1000 / 3600) / Mw

        elif self.critical_flow is False:

            r = self.cp.P / P_discharge

            F2 = np.sqrt((self.cp.k / (self.cp.k - 1)) * r ** (2 / self.cp.k) *
                         ((1 - r ** ((self.cp.k - 1) / self.cp.k)) / (1 - r)))

            W = (((A_relief(self.scenario.D_BPR) * 1000000 * F2 * Kd / 17.9) * np.sqrt(
                Mw * self.cp.P * (self.cp.P - P_discharge) / (Z * c2k(T)))) * 1000 / 3600) / Mw

        return W

    def voidfrac(self, z):

        C0 = 1.5

        F = np.empty(1)

        if self.scenario.flow_regime is 'churn-turbulent':
            F[0] = ((z[0] * ((1 - z[0]) ** 2)) / ((1 - z[0] ** 3) * (1 - C0 * z[0]))) - (self.jgx / self.Ui)
        elif self.scenario.flow_regime is 'bubbly':
            F[0] = ((self.jgx / self.Ui) / (2 + C0 * (self.jgx / self.Ui))) - z[0]

        return F

    def two_phase(self, T, P_discharge, D_RD, Dv, Kd_RD, data, flag):

        Mw = self.H2O.y * cc.MH2O + self.H2O2.y * cc.MH2O2 + self.O2.y * cc.MO2  # (g/mol)

        Z = self.H2O.Z * self.H2O.y + self.H2O2.Z * self.H2O2.y + self.O2.Z * self.O2.y

        if P_discharge > P_crit(P_discharge, T, yH2O, yH2O2, yO2):
            r = P_discharge / P

            F2 = sqrt((k / (k - 1)) * (r ** (2 / k)) * ((1 - r ** ((k - 1) / k)) / (1 - r)))

            W_vap = (F2 * Kd_RD / 17.9) * sqrt((Mw * P * (P - P_discharge)) / (Z * C2K(T))) * (
                        1000000 / 3600)  # (kg/(m**2*s))
        else:
            C = 520 * sqrt(k * ((2 / (k + 1)) ** ((k + 1) / (k - 1))))

            W_vap = ((C * P * Kd_RD / 13160) * sqrt(Mw / (C2K(T) * Z))) * (1000000 / 3600)

        jgx = (A_relief(D_RD) * W_vap) / (pG * A_relief(Dv * 39.3701) * 1000)

        if flag == 1:  # churn turbulent
            Ux_factor = 1.53
        else:  # bubbly
            Ux_factor = 1.18

        Ux = Ux_factor * (sigma * g * 1000 * (pL - pG)) ** (1 / 4) / sqrt(1000 * pL)  # Churn-turbulent regime

        alpha = fsolve(voidfrac, 0.8, args=(jgx, Ux, flag))

        alphaves = (VR - VL) / VR

        if alpha <= alphaves:
            TF = 0
            Xm = empty(1)
            jgi = empty(1)
        else:
            TF = 1

            C0 = 1.5

            if flag == 1:  # churn turbulent
                jgi = 2 * alphaves * Ux / (1 - C0 * alphaves)
                a_m = 2 * alphaves / (1 + C0 * alphaves)
            else:  # bubbly
                jgi = alphaves * (1 - alphaves) ** 2 * Ux / ((1 - alphaves ** 3) * (1 - C0 * alphaves))
                a_m = alphaves

            Xm = a_m * pG / (a_m * pG + (1 - a_m) * pL)

        n_vent = W_vap * A_relief(D_RD) * 1000 / Mw

        return [TF, n_vent, jgi, Xm]