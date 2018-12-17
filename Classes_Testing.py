import numpy as np
from scipy import optimize as opt
from threading import Thread

def c2k(temperature):
    """
        Convert temperature from degrees Celsius to kelvin.
    """
    return temperature + 273.15

def g2l(volume_gal):
    """
        Convert volume from gallons (gal) to litres (L).
    """
    return volume_gal*3.78541

class Common_Constants:
    def __init__(self):
        # Molecular weights (g/mol)
        self.MH2O = 18.01528
        self.MO2 = 31.998
        self.MH2O2 = 34.0147

        # Universal constants
        self.R = 8.3145
        self.g = 9.81

        # System properties
        self.Patm = 101.325

        # Critical temperatures (K)
        self.TcO2 = 154.6
        self.TcH2O = 647.096
        self.TcH2O2 = 728

        # Critical pressures (kPa)
        self.PcO2 = 5050
        self.PcH2O = 22060
        self.PcH2O2 = 22000

class Vessel_Parameters:
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
        self.VR = g2l(reactor_volume)
        self.Ux = heat_transfer_coefficient
        self.AR = aspect_ratio
        self.D = 2 * ((self.VR * 0.001) / (2 * np.pi * self.AR)) ** (1 / 3)
        self.h = self.D * self.AR

class Common_Functions(Common_Constants):
    def compress_solver(self, z):
        """
            Solver for calculating compressibility factors using Redlich-Kwong equations of state.
        """

        A = 0.42748 * self.Pr / self.Tr ** 2.5
        B = 0.08664 * self.Pr / self.Tr

        F = np.empty(1)
        F[0] = z[0] ** 3 - z[0] ** 2 + (A - B - B ** 2) * z[0] - A * B

        return F

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
        F[3] = PH2O + PH2O2 + (ZO2 * self.nO2 * self.R * c2k(self.T)) / (self.VR - VL) - P
        F[4] = PH2O - xH2O * H2O.Psat * H2O.gamma
        F[5] = PH2O2 - xH2O2 * H2O2.Psat * H2O2.gamma
        F[6] = yH2O - PH2O / P
        F[7] = yH2O2 - PH2O2 / P
        F[8] = yO2 - (ZO2 * self.nO2 * self.R * c2k(self.T)) / ((self.VR - VL) * P)
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

        PO2 = self.nO2 * self.R * c2k(self.T) / VG

        return [xH2O, xH2O2, xO2, yH2O, yH2O2, yO2, nL, nG, P, PH2O, PH2O2, PO2, VL, VG, ZO2]

class Water(Common_Functions):
    def __init__(self, temperature, pressure, xH2O):
        """
            Initializes instance of water (H2O).

            Arguments:
            T: Temperature of liquid in degrees Celsius
            P: Pressure in headspace of reactor in kilopascals
            xH2O: Mole fraction of water in liquid phase
        """
        Common_Functions.__init__(self)
        self.density = None
        self.Psat = None
        self.Pr = None
        self.cpl = None
        self.cpg = None
        self.Tr = None
        self.gamma = None
        self.Z = None
        self.T = temperature
        self.P = pressure
        self.xH2O = xH2O
        self.Tref = c2k(self.T) / 1000

    def p_L(self):
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

        self.density = ((A + B * self.T + C * self.T ** 2 + D * self.T ** 3 + E * self.T ** 4 + F * self.T ** 5) / (1 + G * self.T)) / 1000

    def antoine(self):
        """
            Calculate the saturation pressure of water in kPa.
        """

        if self.T > 99:
            A = 8.14019
            B = 1810.94
            C = 244.485
        else:
            A = 8.07131
            B = 1730.63
            C = 233.426

        self.Psat = (10 ** (A - (B / (C + self.T)))) * (101.325 / 760)

    def heat_capacity_L(self):
        """
            Calculate the constant pressure heat capacity of liquid water in J/(g*K).
        """

        A = -203.606
        B = 1523.290
        C = -3196.413
        D = 2474.455
        E = 3.855326

        self.cpl = (A + B * self.Tref + C * self.Tref ** 2 + D * self.Tref ** 3 + E / self.Tref ** 2) / self.MH2O

    def heat_capacity_G(self):
        """
            Calculate the constant pressure heat capacity of water vapour in J/(g*K).
        """
        A = 30.09200
        B = 6.832514
        C = 6.793435
        D = -2.534480
        E = 0.082139

        self.cpg = (A + B * self.Tref + C * self.Tref ** 2 + D * self.Tref ** 3 + E / self.Tref ** 2) / self.MH2O

    def reduced_tempertaure(self):
        """
            Calculate the reduced temperature for water.
            For use in estimating compressibility factors of non-ideal vapours using RK-EOS.
        """

        self.Tr = c2k(self.T) / self.TcH2O

    def reduced_pressure(self):
        """
            Calculate the reduced pressure for water.
            For use in estimating compressibility factors of non-ideal vapours using RK-EOS.
        """

        self.Pr = self.P / self.PcH2O

    def activity(self):
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

        if c2k(self.T) > 0 and c2k(self.T) <= 317.636:
            Ba = Ca0 + ((Ca1 * Ca2) / (np.pi * (Ca2 ** 2 + (c2k(self.T) - Ca3) ** 2)))

        if c2k(self.T) > 317.636 and c2k(self.T) <= 348.222:
            Ba = ((Ca0 + ((Ca1 * Ca2) / (np.pi * (Ca2 ** 2 + ((c2k(self.T) - Ca3) ** 2))))) + (
                        P12 * c2k(self.T) ** 2 + P11 * c2k(self.T) + P10)) / 2

        elif c2k(self.T) > 348.222 and c2k(self.T) <= 391.463:
            Ba = P22 * c2k(self.T) ** 2 + P21 * c2k(self.T) + P20

        else:
            Ba = -612.9613

        Bb = Ca01 + ((Ca11 * Ca21) / (np.pi * (Ca21 ** 2 + ((c2k(self.T) - Ca31) ** 2))))

        Bc = Ca02 + (Ca12 / (1 + np.exp(Ca22 * (c2k(self.T) - Ca32))))

        Bd = Ca03 + (Ca13 / (1 + np.exp(Ca23 * (c2k(self.T) - Ca33))))

        self.gamma = np.exp(((1 - self.xH2O ** 2) / (self.R * c2k(self.T))) * (
                    Ba + Bb * (1 - 4 * self.xH2O) + Bc * (1 - 2 * self.xH2O) * (1 - 6 * self.xH2O) + Bd * ((1 - 2 * self.xH2O) ** 2) * (
                        1 - 8 * self.xH2O)))

    def compress(self):

        initval = np.asarray(0.99)

        self.Z = opt.fsolve(self.compress_solver, initval)

    def runmain(self):
        if __name__ == '__main__':
            Thread(target = self.antoine).start()
            Thread(target=self.heat_capacity_L).start()
            Thread(target=self.heat_capacity_G).start()
            Thread(target=self.p_L).start()
            Thread(target=self.reduced_tempertaure).start()
            Thread(target=self.reduced_pressure).start()
            Thread(target=self.activity).start()
            Thread(target=self.compress).start()

class Hydrogen_Peroxide(Common_Functions):
    def __init__(self, temperature, pressure, xH2O):
        """
            Initializes instance of hydrogen peroxide (H2O2).

            Arguments:
            T: Temperature of liquid in degrees Celsius
            P: Pressure in headspace in kilopascals
            xH2O: Mole fraction water in liqid phase
        """

        Common_Functions.__init__(self)
        self.density = None
        self.Psat = None
        self.Pr = None
        self.cpl = None
        self.cpg = None
        self.Tr = None
        self.gamma = None
        self.Z = None
        self.T = temperature
        self.P = pressure
        self.xH2O = xH2O
        self.Tref = c2k(self.T) / 1000

    def p_L(self):
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

        N = Jb + Kb * self.T + Lb * (self.T ** 2) + Mb * (self.T ** 3)
        O = Jc + Kc * self.T + Lc * (self.T ** 2) + Mc * (self.T ** 3)
        P = Jd + Kd * self.T + Ld * (self.T ** 2) + Md * (self.T ** 3)

        if self.T >= 100:
            self.density = 1.2456174226244978
        else:
            self.density = H2O.density + N + O ** 2 + P ** 3

    def antoine(self):
        """
            Calculate the saturation pressure of hydrogen peroxide in kPa.
        """

        D = 7.96917
        E = 1886.76
        F = 220.6

        self.Psat = (10 ** (D - (E / (F + self.T)))) * (101.325 / 760)

    def heat_capacity_L(self):
        """
            Calculate the constant pressure heat capacity of liquid hydrogen peroxide in J/(g*K).
        """

        A = 0.657
        B = 2.11e-4

        self.cpl = (A + B * self.T) * 4.184

    def heat_capacity_G(self):
        """
            Calculate the constant pressure heat capacity of hydrogen peroxide vapour in J/(g*K).
        """

        F = 34.25667
        G = 55.18445
        H = -35.15443
        I = 9.087440
        J = -0.422157

        self.cpg = (F + G * self.Tref + H * self.Tref ** 2 + I * self.Tref ** 3 + J / self.Tref ** 2) / self.MH2O2

    def reduced_tempertaure(self):
        """
            Calculate the reduced temperature for hydrogen peroxide.
            For use in estimating compressibility factors of non-ideal vapours using RK-EOS.
        """

        self.Tr = c2k(self.T) / self.TcH2O2

    def reduced_pressure(self):
        """
            Calculate the reduced pressure for hydrogen peroxide.
            For use in estimating compressibility factors of non-ideal vapours using RK-EOS.
        """

        self.Pr = self.P / self.PcH2O2

    def activity(self):
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

        if c2k(self.T) > 0 and c2k(self.T) <= 317.636:
            Ba = Ca0 + ((Ca1 * Ca2) / (np.pi * (Ca2 ** 2 + (c2k(self.T) - Ca3) ** 2)))

        if c2k(self.T) > 317.636 and c2k(self.T) <= 348.222:
            Ba = ((Ca0 + ((Ca1 * Ca2) / (np.pi * (Ca2 ** 2 + ((c2k(self.T) - Ca3) ** 2))))) + (
                        P12 * c2k(self.T) ** 2 + P11 * c2k(self.T) + P10)) / 2

        elif c2k(self.T) > 348.222 and c2k(self.T) <= 391.463:
            Ba = P22 * c2k(self.T) ** 2 + P21 * c2k(self.T) + P20

        else:
            Ba = -612.9613

        Bb = Ca01 + ((Ca11 * Ca21) / (np.pi * (Ca21 ** 2 + ((c2k(self.T) - Ca31) ** 2))))

        Bc = Ca02 + (Ca12 / (1 + np.exp(Ca22 * (c2k(self.T) - Ca32))))

        Bd = Ca03 + (Ca13 / (1 + np.exp(Ca23 * (c2k(self.T) - Ca33))))

        self.gamma = np.exp((self.xH2O ** 2 / (self.R * c2k(self.T))) * (
                    Ba + Bb * (3 - 4 * self.xH2O) + Bc * (1 - 2 * self.xH2O) * (5 - 6 * self.xH2O) + Bd * ((1 - 2 * self.xH2O) ** 2) * (
                        7 - 8 * self.xH2O)))

    def compress(self):

        initval = np.asarray(0.99)

        self.Z = opt.fsolve(self.compress_solver, initval)

    def runmain(self):
        if __name__ == '__main__':
            Thread(target = self.antoine).start()
            Thread(target=self.heat_capacity_L).start()
            Thread(target=self.heat_capacity_G).start()
            Thread(target=self.p_L).start()
            Thread(target=self.reduced_tempertaure).start()
            Thread(target=self.reduced_pressure).start()
            Thread(target=self.activity).start()
            Thread(target=self.compress).start()

class Oxygen(Common_Functions):
    def __init__(self, temperature, pressure):
        """
            Initializes instance of oxygen (O2).

            Arguments:
            T: Temperature of medium in degrees Celsius
            P: Pressure in headspace in kilopascals
        """

        Common_Functions.__init__(self)
        self.cpg = None
        self.Tr = None
        self.Pr = None
        self.Z = None
        self.T = temperature
        self.P = pressure
        self.Tref = c2k(self.T) / 1000

    def heat_capacity_G(self):
        """
            Calculate the constant pressure heat capacity of oxygen gas in J/(g*K).
        """

        K = 31.32234
        L = -20.23531
        M = 57.86644
        N = -36.50624
        O = -0.007374

        self.cpg = (K + L * self.Tref + M * self.Tref ** 2 + N * self.Tref ** 3 + O / self.Tref ** 2) / self.MO2

    def reduced_tempertaure(self):
        """
            Calculate the reduced temperature for oxygen.
            For use in estimating compressibility factors of non-ideal vapours using RK-EOS.
        """

        self.Tr = c2k(self.T) / self.TcO2

    def reduced_pressure(self):
        """
            Calculate the reduced pressure for oxygen.
            For use in estimating compressibility factors of non-ideal vapours using RK-EOS.
        """

        self.Pr = self.P / self.PcO2

    def compress(self):

        initval = np.asarray(0.99)

        self.Z = opt.fsolve(self.compress_solver, initval)

    def runmain(self):
        if __name__ == '__main__':
            Thread(target = self.heat_capacity_G).start()
            Thread(target=self.reduced_tempertaure).start()
            Thread(target=self.reduced_pressure).start()
            Thread(target=self.compress).start()

class Common_Properties(Common_Constants):
    def __init__(self, temperature):
        """
            Initializes instance for physical properties common to the system.

            Arguments:
            T: Temperature of liquid in degrees Celsius
        """

        Common_Constants.__init__(self)
        self.st = None
        self.dhvap = None
        self.T = temperature

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
        if __name__ == '__main__':
            Thread(target = self.surface_tension).start()
            Thread(target=self.enthvap).start()

# class Initial_Conditions(Common_Functions):
#     def __init__(self, temperature, pressure, H2O2_massfraction, reactor_charge):
#         """
#             Initializes instance for starting conditions of the system.
#
#             Arguments:
#                 temperature: Temperature of liquid in degrees Celsius
#                 pressure: Starting headspace pressure in kilopascals
#                 H2O2_massfraction: Starting H2O2 mass fraction in liquid phase
#                 reactor_charge: Total mass of starting reactants loaded to reactor (H2O + H2O2) in kg
#
#             Outputs:
#                 zH2O: total mole fraction of water
#                 zH2O2: total mole fraction of hydrogen peroxide
#                 zO2: total mole fraction of oxygen
#                 mH2O: total mass water in system in kg
#                 mH2O2: total mass hydrogen peroxide in system in kg
#                 nH2O: total amount water in system in mol
#                 nH2O2: total amount hydrogen peroxide in system in mol
#                 nO2: total amount oxygen in system in mol
#                 ntotal: total amount in system in mol
#                 VG: headspace volume in L
#         """
#
#         self.T = temperature
#         self.P = pressure
#         self.XH2O2 = H2O2_massfraction
#         self.mR = reactor_charge
#
#         self.zH2O = None
#         self.zH2O2 = None
#         self.zO2 = None
#         self.mH2O = None
#         self.mH2O2 = None
#         self.nH2O = None
#         self.nH2O2 = None
#         self.nO2 = None
#         self.ntotal = None
#         self.VG = None
#
#     def initial_conditions(self):
#         """
#             Calculates thermodynamically stable starting conditions for the reactor.
#         """
#
#         self.mH2O2 = self.mR * self.XH2O2
#         self.mH2O = self.mR * (1 - self.XH2O2)
#
#         self.nH2O2 = self.mH2O2 * 1000 / cc.MH2O2
#         self.nH2O = self.mH2O * 1000 / cc.MH2O
#
#         self.VG = vp.VR - self.mH2O / H2O.density - self.mH2O2 / H2O2.density
#
#         self.nO2 = self.P * self.VG / (cc.R * c2k(self.T))
#
#         self.ntotal = self.nH2O + self.nH2O2 + self.nO2
#
#         self.zH2O = self.nH2O / self.ntotal
#         self.zH2O2 = self.nH2O2 / self.ntotal
#         self.zO2 = self.nO2 / self.ntotal
#
#         self.xH2O = (self.mH2O / cc.MH2O) / ((self.mH2O / cc.MH2O) + (self.mH2O2 / cc.MH2O2))
#
#         H2O = Water(self.T, self.P, self.xH2O)
#         H2O2 = Hydrogen_Peroxide(self.T, self.P, xH2O)
#         O2 = Oxygen(self.T, self.P)
#
#         H2O.runmain()
#         H2O2.runmain()
#         O2.runmain()



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