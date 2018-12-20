import numpy as np
from scipy import optimize as opt
from threading import Thread

class Initial_Conditions():
    def __init__(self, temperature, pressure, H2O2_massfraction, reactor_charge):
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
        self.P = pressure
        self.XH2O2 = H2O2_massfraction
        self.mR = reactor_charge

        self.zH2O = None
        self.zH2O2 = None
        self.zO2 = None
        self.mH2O = None
        self.mH2O2 = None
        self.mO2 = None
        self.nH2O = None
        self.nH2O2 = None
        self.nO2 = None
        self.ntotal = None

    def __call__(self, *args, **kwargs):


    def initial_conditions(self):
        """
            Calculates thermodynamically stable starting conditions for the reactor.
        """

        self.mH2O2 = self.mR * self.XH2O2
        self.mH2O = self.mR * (1 - self.XH2O2)

        self.nH2O2 = self.mH2O2 * 1000 / self.MH2O2
        self.nH2O = self.mH2O * 1000 / self.MH2O

        VG = self.VR - self.mH2O / self.density - self.mH2O2 / self.density

        self.nO2 = self.P * VG / (self.R * self.c2k(self.T))
        self.mO2 = self.mO2*self.MO2

        self.ntotal = self.nH2O + self.nH2O2 + self.nO2

        self.zH2O = self.nH2O / self.ntotal
        self.zH2O2 = self.nH2O2 / self.ntotal
        self.zO2 = self.nO2 / self.ntotal

        self.xH2O = (self.mH2O / self.MH2O) / ((self.mH2O / self.MH2O) + (self.mH2O2 / self.MH2O2))

        H2O = Water(self.T, self.P, self.xH2O)
        H2O2 = Hydrogen_Peroxide(self.T, self.P, xH2O)
        O2 = Oxygen(self.T, self.P)

        H2O.runmain()
        H2O2.runmain()
        O2.runmain()

class Vessel_Parameters(Conversion):
    def __init__(self, reactor_volume, heat_transfer_coefficient=450, aspect_ratio=1.5):
        """
            Provides access to common reactor vessel properties.

            Attributes:
                D: Reactor diameter in meters
                h: Reactor height in meters
                AR: Reactor aspect ratio (h/D)
                Ux: Reactor heat transfer coefficient in (J/(s*m**2*K))
                VR: Reactor volume in gallons
        """
        self.VR = Conversion.g2l(reactor_volume)
        self.Ux = heat_transfer_coefficient
        self.AR = aspect_ratio
        self.D = 2 * ((self.VR * 0.001) / (2 * np.pi * self.AR)) ** (1 / 3)
        self.h = self.D * self.AR

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

    A = 0.42748 * O2.Pr / O2.Tr ** 2.5
    B = 0.08664 * O2.Pr / O2.Tr

    F = np.empty(12)
    F[0] = nL * xH2O + nG * yH2O - self.ntotal * self.zH2O
    F[1] = nL * xH2O2 + nG * yH2O2 - self.ntotal * self.zH2O2
    F[2] = nL + nG - self.ntotal
    F[3] = PH2O + PH2O2 + (ZO2 * self.nO2 * self.R * Conversion.c2k(self.T)) / (self.VR - VL) - P
    F[4] = PH2O - xH2O * H2O.Psat * H2O.gamma
    F[5] = PH2O2 - xH2O2 * H2O2.Psat * H2O2.gamma
    F[6] = yH2O - PH2O / P
    F[7] = yH2O2 - PH2O2 / P
    F[8] = yO2 - (ZO2 * self.nO2 * self.R * Conversion.c2k(self.T)) / ((self.VR - VL) * P)
    F[9] = xH2O + xH2O2 - yH2O - yH2O2 - yO2
    F[10] = VL - nL * ((xH2O * self.MH2O / (H2O.density * 1000)) + (xH2O2 * self.MH2O2 / (H2O2.density * 1000)))
    F[11] = ZO2 ** 3 - ZO2 ** 2 + (A - B - B ** 2) * ZO2 - A * B

    return F

def equilibrate(self, data):
        """
            Works with VLE to calculate thermodynamically stable vessel conditions given data.

            Arguments:
            data: list containing best initial guess for vessel thermodynamic conditions
                contains:
                (
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
                )
        """

        initvals = np.asarray(data)

        outputA = opt.fsolve(self.VLE, initvals,
                         args=(self.ntotal, self.zH2O, self.zH2O2, self.T, self.nO2))

        (xH2O, xH2O2, yH2O, yH2O2, yO2, nL, nG, P, PH2O, PH2O2, VL, ZO2) = outputA

        VG = self.VR - VL

        xO2 = 0

        PO2 = self.nO2 * self.R * Conversion.c2k(self.T) / VG

        return [xH2O, xH2O2, xO2, yH2O, yH2O2, yO2, nL, nG, P, PH2O, PH2O2, PO2, VL, VG, ZO2]

T = 25
P = 150
xH2O = 0.5

H2O = Water(T, P, xH2O)
H2O.runmain()

H2O2 = Hydrogen_Peroxide(T, P, xH2O)
H2O2.runmain()

O2 = Oxygen(T, P)
O2.runmain()

pp = Common_Properties(T)
pp.runmain()

print('surface tension = ' + str(pp.st))
print('saturation pressure = ' + str(H2O2.Psat))
print('enthalpy of vaporization = ' + str(pp.dhvap))
print('vapour heat capacity = ' + str(H2O2.cpg))
print('liquid heat capacity = ' + str(H2O2.cpl))
print('liquid density = ' + str(H2O2.density))
print('reduced temperature = ' + str(H2O2.Tr))
print('activity = ' + str(H2O2.gamma))
print('compressibility = ' + str(H2O2.Z))