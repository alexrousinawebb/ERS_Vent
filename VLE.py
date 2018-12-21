"""
VAPOUR LIQUID EQUILIBRIUM VLE
===========================================================================
Module for performing VLE calculations.
"""

import numpy as np
from scipy import optimize as opt
from Conversion import c2k
import Constant_Lib as constant
import Property_Lib as pl

class Solvers():
    def VLE(self, z):
        """
            Solver for calculating thermodynamically stable equilibrium conditions given overall composition,
            temperature, and quantity of material.
        """

        xH2O = z[0]
        xH2O2 = z[1]
        yH2O = z[2]
        yH2O2 = z[3]
        yO2 = z[4]
        nL = z[5]
        nG = z[6]
        P = z[7]
        PH2O = z[8]
        PH2O2 = z[9]
        VL = z[10]
        ZO2 = z[11]

        H2O = pl.Water(self.T, P)
        H2O2 = pl.Hydrogen_Peroxide(self.T, P)
        O2 = pl.Oxygen(self.T, P)

        A = 0.42748 * O2.Pr / O2.Tr ** 2.5
        B = 0.08664 * O2.Pr / O2.Tr

        F = np.empty(12)
        F[0] = nL * xH2O + nG * yH2O - self.ntotal * self.zH2O
        F[1] = nL * xH2O2 + nG * yH2O2 - self.ntotal * self.zH2O2
        F[2] = nL + nG - self.ntotal
        F[3] = PH2O + PH2O2 + (ZO2 * self.nO2 * constant.R * c2k(self.T)) / (self.VR - VL) - P
        F[4] = PH2O - xH2O * H2O.Psat * H2O.gamma
        F[5] = PH2O2 - xH2O2 * H2O2.Psat * H2O2.gamma
        F[6] = yH2O - PH2O / P
        F[7] = yH2O2 - PH2O2 / P
        F[8] = yO2 - (ZO2 * self.nO2 * constant.R * c2k(self.T)) / ((self.VR - VL) * P)
        F[9] = xH2O + xH2O2 - yH2O - yH2O2 - yO2
        F[10] = VL - nL * ((xH2O * constant.MH2O / (H2O.density * 1000)) + (xH2O2 * constant.MH2O2 / (H2O2.density * 1000)))
        F[11] = ZO2 ** 3 - ZO2 ** 2 + (A - B - B ** 2) * ZO2 - A * B

        return F

class Equilibrate(Solvers):
    def equilibrate(self, data):
        """
            Works with VLE to calculate thermodynamically stable vessel conditions given data.

            Arguments:
            data: list containing best initial guess for vessel thermodynamic conditions
                contains:
                [
                xH2O (mol fraction water in liquid phase),
                xH2O2 (mol fraction hydrogen peroxide in liquid phase),
                yH2O (mol fraction water in vapour phase),
                yH2O2 (mol fraction hydrogen peroxide in vapour phase),
                yO2 (mol fraction oxygen in vapour phase),
                nL (total moles of liquid in vessel),
                nG (total moles of vapour in vessel),
                P (total vessel pressure),
                PH2O (water partial pressure),
                PH2O2 (hydrogen peroxide partial pressure),
                VL (liquid volume),
                ZO2 (oxygen compressibility factor)
                ]
        """

        initvals = np.asarray(data)

        (pl.Water.x, pl.Hydrogen_Peroxide.x, pl.Water.y, pl.Hydrogen_Peroxide.y, pl.Oxygen.y, self.nL,
         self.nG, self.P, pl.Water.P, pl.Hydrogen_Peroxide.P, self.VL, ZO2) = opt.fsolve(self.VLE, initvals)

        self.VG = self.VR - self.VL

        pl.Oxygen.x = 0

        pl.Oxygen.P = self.nO2 * constant.R * c2k(self.T) / self.VG

class Initial_Conditions(Equilibrate):
    def __init__(self, temperature, H2O2_massfraction, reactor_charge, scenario, pressure=101):
        """
            Initializes instance for starting conditions of the system.

            Arguments:
                temperature: Temperature of liquid in degrees Celsius
                pressure: Starting headspace pressure in kilopascals
                H2O2_massfraction: Starting H2O2 mass fraction in liquid phase
                reactor_charge: Total mass of starting reactants loaded to reactor (H2O + H2O2) in kg

            Outputs:
                zH2O: total mole fraction of water
                zH2O2: total mole fraction of hydrogen peroxide
                zO2: total mole fraction of oxygen
                mH2O: total mass water in system in kg
                mH2O2: total mass hydrogen peroxide in system in kg
                nH2O: total amount water in system in mol
                nH2O2: total amount hydrogen peroxide in system in mol
                nO2: total amount oxygen in system in mol
                ntotal: total amount in system in mol
                VG: headspace volume in L
        """

        self.T = temperature
        self.XH2O2 = H2O2_massfraction
        self.mR = reactor_charge
        self.VR = scenario.VR
        self.P0 = pressure

        self.ntotal = None
        self.nL = None
        self.nG = None

        self.P = None

        self.VL = None
        self.VG = None

        self.initial_conditions()

    def initial_conditions(self):
        """
            Calculates thermodynamically stable starting conditions for the reactor.
        """

        H2O = pl.Water(self.T)
        H2O2 = pl.Hydrogen_Peroxide(self.T)

        mH2O2 = self.mR * self.XH2O2
        mH2O = self.mR * (1 - self.XH2O2)

        nH2O2 = mH2O2 * 1000 / constant.MH2O2
        nH2O = mH2O * 1000 / constant.MH2O

        VG = self.VR - mH2O / H2O.density - mH2O2 / H2O2.density
        VL = mH2O / H2O.density - mH2O2 / H2O2.density

        self.nO2 = self.P0 * VG / (constant.R * c2k(self.T))
        mO2 = self.nO2*constant.MO2

        self.ntotal = nH2O + nH2O2 + self.nO2

        self.zH2O = nH2O / self.ntotal
        self.zH2O2 = nH2O2 / self.ntotal
        zO2 = self.nO2 / self.ntotal

        data = [self.zH2O, self.zH2O2, self.zH2O, self.zH2O2, zO2, self.ntotal, self.nO2, self.P0, H2O.Psat, H2O2.Psat, VL, 0.95]

        self.equilibrate(data)

        #  Set water attributes
        pl.Water.m = mH2O
        pl.Water.z = self.zH2O
        pl.Water.n = nH2O

        #  Set hydrogen peroxide attributes
        pl.Hydrogen_Peroxide.m = mH2O2
        pl.Hydrogen_Peroxide.z = self.zH2O2
        pl.Hydrogen_Peroxide.n = nH2O2

        #  Set oxygen attributes
        pl.Oxygen.m = mO2
        pl.Oxygen.z = zO2
        pl.Oxygen.n = self.nO2