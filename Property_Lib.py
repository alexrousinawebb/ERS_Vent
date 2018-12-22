"""
PROPERTY LIBRARY
===========================================================================
Module for creating instances of physical properties for water, hydrogen peroxide, oxygen, and common system physical properties.

Water/Hydrogen Peroxide:
    Density
    Saturation Pressure
    Liquid Heat Capacity
    Vapour Heat Capacity
    Activity Coefficient
    Compressibility Factor from RK-EOS
    Reduced Temperature
    Reduced Pressure
    
Oxygen:
    Compressibility Factor from RK-EOS
    Reduced Temperature
    Reduced Pressure
    Vapour Heat Capacity
    
Common:
    Surface Tension
    Enthalpy of Vaporization    
"""

import numpy as np
from scipy import optimize as opt
from threading import Thread as th
from Conversion import c2k
from EOS import RK_EOS
import Constant_Lib as constant
import VLE

class Water(RK_EOS):
    def __init__(self, temperature=25, pressure=101):
        """
            Initializes instance of water (H2O).

            Arguments:
            T: Temperature of liquid in degrees Celsius
            P: Pressure in headspace of reactor in kilopascals
            xH2O: Mole fraction of water in liquid phase
        """

        #  Solve for
        self.density = None
        self.cpl = None
        self.cpg = None
        self.gamma = None
        self.Z = None
        self.Psat = None
        self.Pr = None
        self.Tr = None
        self.m = None
        self.x = 1
        self.y = None
        self.z = None
        self.n = None
        self.P = None


        self.p_L(temperature)
        self.antoine(temperature)
        self.heat_capacity_L(temperature)
        self.heat_capacity_G(temperature)
        self.activity(temperature, self.x)
        self.compress(temperature, pressure)

    def p_L(self, T):
        """
            Calculate density of liquid water in kg/L.
        """

        A = 999.83952
        B = 16.945176
        C = -7.987040e-3
        D = -46.170461e-6
        E = 105.56302e-9
        F = -280.54253e-12
        G = 16.897850e-3

        self.density = ((A + B * T + C * T ** 2 + D * T ** 3 + E * T ** 4 + F * T ** 5) / (1 + G * T)) / 1000

    def antoine(self, T):
        """
            Calculate the saturation pressure of water in kPa.
        """

        if T > 99:
            A = 8.14019
            B = 1810.94
            C = 244.485
        else:
            A = 8.07131
            B = 1730.63
            C = 233.426

        self.Psat = (10 ** (A - (B / (C + T)))) * (101.325 / 760)

    def heat_capacity_L(self, T):
        """
            Calculate the constant pressure heat capacity of liquid water in J/(g*K).
        """

        A = -203.606
        B = 1523.290
        C = -3196.413
        D = 2474.455
        E = 3.855326

        Tref = c2k(T) / 1000

        self.cpl = (A + B * Tref + C * Tref ** 2 + D * Tref ** 3 + E / Tref ** 2) / constant.MH2O

    def heat_capacity_G(self, T):
        """
            Calculate the constant pressure heat capacity of water vapour in J/(g*K).
        """
        A = 30.09200
        B = 6.832514
        C = 6.793435
        D = -2.534480
        E = 0.082139

        Tref = c2k(T) / 1000

        self.cpg = (A + B * Tref + C * Tref ** 2 + D * Tref ** 3 + E / Tref ** 2) / constant.MH2O

    def activity(self, T, xH2O):
        """
            Calculate the activity coefficient for water (H2O).
        """

        Ca0 = -999.883
        Ca1 = -2499.584
        Ca2 = 8.261924
        Ca3 = 327.4487
        P10 = 17418.34
        P11 = -109.9125
        P12 = 0.1663847
        P20 = -6110.401
        P21 = 28.08669
        P22 = -0.03587408
        Ca01 = 126.7385
        Ca11 = -2558.776
        Ca21 = 12.33364
        Ca31 = 343.105
        Ca02 = 63.18354
        Ca12 = -149.9278
        Ca22 = 0.4745954
        Ca32 = 348.1642
        Ca03 = 59.42228
        Ca13 = -199.2644
        Ca23 = 0.8321514
        Ca33 = 346.2121

        T_K = c2k(T)

        if T_K > 0 and T_K <= 317.636:
            Ba = Ca0 + ((Ca1 * Ca2) / (np.pi * (Ca2 ** 2 + (T_K - Ca3) ** 2)))

        elif T_K > 317.636 and T_K <= 348.222:
            Ba = ((Ca0 + ((Ca1 * Ca2) / (np.pi * (Ca2 ** 2 + ((T_K - Ca3) ** 2))))) + (
                        P12 * T_K ** 2 + P11 * T_K + P10)) / 2

        elif T_K > 348.222 and T_K <= 391.463:
            Ba = P22 * T_K ** 2 + P21 * T_K + P20

        else:
            Ba = -612.9613

        Bb = Ca01 + ((Ca11 * Ca21) / (np.pi * (Ca21 ** 2 + ((T_K - Ca31) ** 2))))

        Bc = Ca02 + (Ca12 / (1 + np.exp(Ca22 * (T_K - Ca32))))

        Bd = Ca03 + (Ca13 / (1 + np.exp(Ca23 * (T_K - Ca33))))

        self.gamma = np.exp(((1 - xH2O ** 2) / (constant.R * T_K)) * (
                    Ba + Bb * (1 - 4 * xH2O) + Bc * (1 - 2 * xH2O) * (1 - 6 * xH2O) + Bd * ((1 - 2 * xH2O) ** 2) * (
                        1 - 8 * xH2O)))

    def compress(self, T, P):

        initval = np.asarray(0.99)

        self.Tr = self.reduced_tempertaure(T, constant.TcH2O)
        self.Pr = self.reduced_pressure(P, constant.PcH2O)

        self.Z = float(opt.fsolve(self.compress_solver, initval))

    def inherit_properties(self, equilibrium_conditions):
        self.x = equilibrium_conditions.xH2O
        self.y = equilibrium_conditions.yH2O
        self.P = equilibrium_conditions.PH2O
        self.z = equilibrium_conditions.zH2O
        self.m = equilibrium_conditions.mH2O
        self.n = equilibrium_conditions.nH2O

class Hydrogen_Peroxide(RK_EOS):
    def __init__(self, temperature=25, pressure=101):
        """
            Initializes instance of hydrogen peroxide (H2O2).

            Arguments:
            T: Temperature of liquid in degrees Celsius
            P: Pressure in headspace in kilopascals
            xH2O: Mole fraction water in liqid phase
        """

        #  Solve For
        self.density = None
        self.cpl = None
        self.cpg = None
        self.gamma = None
        self.Z = None
        self.Psat = None
        self.Pr = None
        self.Tr = None
        self.m = None
        self.x = 0
        self.y = None
        self.z = None
        self.n = None
        self.P = None

        #  Calculate physical properties based on input
        self.p_L(temperature)
        self.antoine(temperature)
        self.heat_capacity_G(temperature)
        self.heat_capacity_L(temperature)
        self.activity(temperature, self.x)
        self.compress(temperature, pressure)

    def p_L(self, T):
        """
            Calculate density of liquid hydrogen peroxide in kg/L.
        """

        Jb = 0.39763
        Jc = 0.02206
        Jd = 0.05187
        Kb = -2.8732E-3
        Kc = 3.5357E-3
        Kd = -1.9414E-3
        Lb = 3.2488E-5
        Lc = -6.0947E-5
        Ld = 3.9061E-5
        Mb = -1.6363E-7
        Mc = 3.6165E-7
        Md = -2.5500E-7

        N = Jb + Kb * T + Lb * (T ** 2) + Mb * (T ** 3)
        O = Jc + Kc * T + Lc * (T ** 2) + Mc * (T ** 3)
        P = Jd + Kd * T + Ld * (T ** 2) + Md * (T ** 3)

        if T >= 100:
            self.density = 1.2456174226244978
        else:
            H2O = Water(T)

            self.density = H2O.density + N + O ** 2 + P ** 3

    def antoine(self, T):
        """
            Calculate the saturation pressure of hydrogen peroxide in kPa.
        """

        D = 7.96917
        E = 1886.76
        F = 220.6

        self.Psat = (10 ** (D - (E / (F + T)))) * (101.325 / 760)

    def heat_capacity_L(self, T):
        """
            Calculate the constant pressure heat capacity of liquid hydrogen peroxide in J/(g*K).
        """

        A = 0.657
        B = 2.11e-4

        self.cpl = (A + B * T) * 4.184

    def heat_capacity_G(self, T):
        """
            Calculate the constant pressure heat capacity of hydrogen peroxide vapour in J/(g*K).
        """

        F = 34.25667
        G = 55.18445
        H = -35.15443
        I = 9.087440
        J = -0.422157

        Tref = c2k(T) / 1000

        self.cpg = (F + G * Tref + H * Tref ** 2 + I * Tref ** 3 + J / Tref ** 2) / constant.MH2O2

    def activity(self, T, xH2O):
        """
            Calculate the activity coefficient for hydrogen peroxide (H2O).
        """

        Ca0 = -999.883
        Ca1 = -2499.584
        Ca2 = 8.261924
        Ca3 = 327.4487
        P10 = 17418.34
        P11 = -109.9125
        P12 = 0.1663847
        P20 = -6110.401
        P21 = 28.08669
        P22 = -0.03587408
        Ca01 = 126.7385
        Ca11 = -2558.776
        Ca21 = 12.33364
        Ca31 = 343.105
        Ca02 = 63.18354
        Ca12 = -149.9278
        Ca22 = 0.4745954
        Ca32 = 348.1642
        Ca03 = 59.42228
        Ca13 = -199.2644
        Ca23 = 0.8321514
        Ca33 = 346.2121

        T_K = c2k(T)

        if T_K > 0 and T_K <= 317.636:
            Ba = Ca0 + ((Ca1 * Ca2) / (np.pi * (Ca2 ** 2 + (T_K - Ca3) ** 2)))

        elif T_K > 317.636 and T_K <= 348.222:
            Ba = ((Ca0 + ((Ca1 * Ca2) / (np.pi * (Ca2 ** 2 + ((T_K - Ca3) ** 2))))) + (
                        P12 * T_K ** 2 + P11 * T_K + P10)) / 2

        elif T_K > 348.222 and T_K <= 391.463:
            Ba = P22 * T_K ** 2 + P21 * T_K + P20

        else:
            Ba = -612.9613

        Bb = Ca01 + ((Ca11 * Ca21) / (np.pi * (Ca21 ** 2 + ((T_K - Ca31) ** 2))))

        Bc = Ca02 + (Ca12 / (1 + np.exp(Ca22 * (T_K - Ca32))))

        Bd = Ca03 + (Ca13 / (1 + np.exp(Ca23 * (T_K - Ca33))))

        self.gamma = np.exp((xH2O ** 2 / (constant.R * T_K)) * (
                    Ba + Bb * (3 - 4 * xH2O) + Bc * (1 - 2 * xH2O) * (5 - 6 * xH2O) + Bd * ((1 - 2 * xH2O) ** 2) * (
                        7 - 8 * xH2O)))

    def compress(self, T, P):

        initval = np.asarray(0.99)

        self.Tr = self.reduced_tempertaure(T, constant.TcH2O2)
        self.Pr = self.reduced_pressure(P, constant.PcH2O2)

        self.Z = float(opt.fsolve(self.compress_solver, initval))

    def inherit_properties(self, equilibrium_conditions):
        self.x = equilibrium_conditions.xH2O2
        self.y = equilibrium_conditions.yH2O2
        self.P = equilibrium_conditions.PH2O2
        self.z = equilibrium_conditions.zH2O2
        self.m = equilibrium_conditions.mH2O2
        self.n = equilibrium_conditions.nH2O2

class Oxygen(RK_EOS):
    def __init__(self, temperature=25, pressure=101):
        """
            Initializes instance of oxygen (O2).

            Arguments:
            T: Temperature of medium in degrees Celsius
            P: Pressure in headspace in kilopascals
        """

        #  Solve for
        self.cpg = None
        self.Z = None
        self.Tr = None
        self.Pr = None
        self.m = None
        self.x = None
        self.y = None
        self.z = None
        self.n = None
        self.P = None

        #  Calculate physical properties
        self.heat_capacity_G(temperature)
        self.compress(temperature, pressure)

    def heat_capacity_G(self, T):
        """
            Calculate the constant pressure heat capacity of oxygen gas in J/(g*K).
        """

        K = 31.32234
        L = -20.23531
        M = 57.86644
        N = -36.50624
        O = -0.007374

        Tref = c2k(T) / 1000

        self.cpg = (K + L * Tref + M * Tref ** 2 + N * Tref ** 3 + O / Tref ** 2) / constant.MO2

    def compress(self, T, P):

        initval = np.asarray(0.99)

        self.Tr = self.reduced_tempertaure(T, constant.TcO2)
        self.Pr = self.reduced_pressure(P, constant.PcO2)

        self.Z = float(opt.fsolve(self.compress_solver, initval))

    def inherit_properties(self, equilibrium_conditions):
        self.x = equilibrium_conditions.xO2
        self.y = equilibrium_conditions.yO2
        self.P = equilibrium_conditions.PO2
        self.z = equilibrium_conditions.zO2
        self.m = equilibrium_conditions.mO2
        self.n = equilibrium_conditions.nO2

class Common_Properties():
    def __init__(self, temperature=25):
        """
            Initializes instance for physical properties common to the system.

            Arguments:
            T: Temperature of liquid in degrees Celsius
        """

        self.st = None
        self.dhvap = None
        self.dvLdt = None
        self.dvGdt = None
        self.T = temperature
        self.TcH2O = constant.TcH2O

    def surface_tension(self):
        """
            Calculate surface tension of liquid water in N/m.
        """

        B = 235.8E-3
        b = -0.625
        u = 1.256

        self.st = B * (((self.TcH2O - c2k(self.T)) / self.TcH2O) ** u) * (1 + b * (self.TcH2O - c2k(self.T)) / self.TcH2O)

    def enthvap(self):
        """
            Calculate the enthalpy of vaporization of water in J/g.
        """

        A = -3e-5
        B = 0.0051
        C = -2.75588
        D = 2500.2

        self.dhvap = A * self.T ** 3 + B * self.T ** 2 + C * self.T + D

    def runmain(self):
        th(target = self.surface_tension).start()
        th(target=self.enthvap).start()

    # def delta_sv_G(self):
    #     self.dvGdt = constant.R * ((1 / (compress(PH2O, T, 'H2O') * constant.MH2O * PH2O)) + (
    #                 1 / (compress(PH2O2, T, 'H2O2') * constant.MH2O2 * PH2O2)) + (1 / (
    #                 compress(PO2, T, 'O2') * constant.MO2 * PO2)))  # This is not exactly correct Z and P also f(P)
    # def delta_sv_L(self):
    #     A = 999.83952
    #     B = 16.945176
    #     C = -7.987040e-3
    #     D = -46.170461e-6
    #     E = 105.56302e-9
    #     F = -280.54253e-12
    #     G = 16.897850e-3
    #
    #     self.dvLdt = -G * (A + B * self.T + C * self.T ** 2 + D * self.T ** 3 + E * self.T ** 4 + F * self.T ** 5) / (
    #             1000 * (G * self.T + 1) ** 2) + (B + 2 * C * self.T + 3 * D * self.T ** 2 + 4 * E * self.T ** 3 +
    #                                              5 * F * self.T ** 4) / (1000 * (G * self.T + 1))

class Kinetics():
    def __init__(self, temperature, kf=1):
        """
            Initializes instance of rate kinetics.

            Arguments:
            T: Temperature of medium in degrees Celsius
            kf: Contamination factor for acceleration of hydrogen peroxide decomposition.
                kf=1 is no contamination.
        """

        #  Solve for
        self.rate = None

        #  Inputs
        self.kf = kf
        self.T = temperature

    def kinetics(self):
        self.rate = constant.A_ar * self.kf * np.exp(-constant.Ea / c2k(self.T))