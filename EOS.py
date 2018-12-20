"""
EOS Equations of State
===========================================================================
Module for calculating EOS properties for non-ideal gasses.

Implemented:
    RK-EOS - Redlich-Kwong
    SRK-EOS - Soave-Redlich-Kwong

"""

from Conversion import c2k
import numpy as np

class RK_EOS():
    def reduced_tempertaure(self):
        """
            Calculate the reduced temperature.
            For use in estimating compressibility factors of non-ideal vapours using RK-EOS.
        """

        return c2k(self.T) / self.Tc

    def reduced_pressure(self):
        """
            Calculate the reduced pressure.
            For use in estimating compressibility factors of non-ideal vapours using RK-EOS.
        """

        return self.P / self.Pc

    def compress_solver(self, z):
        """
            Solver for calculating compressibility factors using Redlich-Kwong equations of state.
        """

        A = 0.42748 * self.Pr / self.Tr ** 2.5
        B = 0.08664 * self.Pr / self.Tr

        F = np.empty(1)
        F[0] = z[0] ** 3 - z[0] ** 2 + (A - B - B ** 2) * z[0] - A * B

        return F