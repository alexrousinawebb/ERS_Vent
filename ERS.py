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
        return (2 / (k + 1)) ** (k / (k - 1)) * P

    def ventflow(self, H2O, H2O2, O2, cp, T, P_discharge):

        Mw = H2O.y * cc.MH2O + H2O2.y * cc.MH2O2 + O2.y * cc.MO2

        Z = H2O.Z * H2O.y + H2O2.Z * H2O2.y + O2.Z * O2.y

        Kd = self.scenario.BPR_max_Cv / (27.66 * A_relief(self.scenario.D_BPR) * 1550)

        if self.critical_flow is True:

            C = 520 * np.sqrt(cp.k * ((2 / (cp.k + 1)) ** ((cp.k + 1) / (cp.k - 1))))

            W = ((A_relief(self.scenario.D_BPR) * 1000000 * C * Kd * cp.P / 13160) * np.sqrt(
                Mw / (c2k(T) * Z)) * 1000 / 3600) / Mw

        elif self.critical_flow is False:

            r = cp.P / P_discharge

            F2 = np.sqrt((cp.k / (cp.k - 1)) * r ** (2 / cp.k) * ((1 - r ** ((cp.k - 1) / cp.k)) / (1 - r)))

            W = (((A_relief(self.scenario.D_BPR) * 1000000 * F2 * Kd / 17.9) * np.sqrt(
                Mw * cp.P * (cp.P - P_discharge) / (Z * c2k(T)))) * 1000 / 3600) / Mw

        return W