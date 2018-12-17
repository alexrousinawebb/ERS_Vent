import numpy as np

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
        self. MH2O2 = 34.0147

        # Universal constants
        self.R = 8.3145
        self.g = 9.81

        # System properties
        self.Patm = 101.325

        # Critical temperatures (K)
        self.TcH2O = 647.096

class Water:
    def __init__(self, temperature):
        """
            Initializes instance of water (H2O).

            Arguments:
            T: Temperature of liquid in degrees Celsius
        """
        self.st = None
        self.density = None
        self.Psat = None
        self.dhvap = None
        self.cpl = None
        self.cpg = None
        self.Tr = None
        self.gamma = None
        self.T = temperature

    def surface_tension(self):
        """
            Calculate surface tension of liquid water in N/m.
        """

        B = 235.8E-3
        b = -0.625
        u = 1.256
        Tc = 647.15

        self.st = B * (((Tc - c2k(self.T)) / Tc) ** u) * (1 + b * (Tc - c2k(self.T)) / Tc)

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

    def enthvap(self):
        """
            Calculate the enthalpy of vaporization of water in J/g.
        """

        A = -3e-5
        B = 0.0051
        C = -2.75588
        D = 2500.2

        self.dhvap = A * self.T ** 3 + B * self.T ** 2 + C * self.T + D

    def heat_capacity_L(self):
        """
            Calculate the constant pressure heat capacity of liquid water in J/(g*K).
        """

        A = -203.606
        B = 1523.290
        C = -3196.413
        D = 2474.455
        E = 3.855326
        Tref = c2k(self.T) / 1000

        self.cpl = (A + B * Tref + C * Tref ** 2 + D * Tref ** 3 + E / Tref ** 2) / cc.MH2O

    def heat_capacity_G(self):
        """
            Calculate the constant pressure heat capacity of water vapour in J/(g*K).
        """
        A = 30.09200
        B = 6.832514
        C = 6.793435
        D = -2.534480
        E = 0.082139
        Tref = c2k(self.T) / 1000

        self.cpg = (A + B * Tref + C * Tref ** 2 + D * Tref ** 3 + E / Tref ** 2) / cc.MH2O

    def reduced_tempertaure(self):
        """
            Calculate the reduced temperature for water.
            For use in estimating compressibility factors of non-ideal vapours using RK-EOS.
        """

        self.Tr = c2k(self.T) / cc.TcH2O

    def activity(self, xH2O):
        """
            Calculate the activity coefficient for water (H2O).

            Arguments:
            xH2O: Mole fraction water in liquid phase
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

        self.gamma = np.exp(((1 - xH2O ** 2) / (cc.R * c2k(self.T))) * (
                    Ba + Bb * (1 - 4 * xH2O) + Bc * (1 - 2 * xH2O) * (1 - 6 * xH2O) + Bd * ((1 - 2 * xH2O) ** 2) * (
                        1 - 8 * xH2O)))




T = 0

cc = Common_Constants()

H2O = Water(T)
H2O.surface_tension()
H2O.antoine()
H2O.enthvap()
H2O.heat_capacity_G()
H2O.heat_capacity_L()
H2O.p_L()
H2O.reduced_tempertaure()

print('surface tension = ' + str(H2O.st))
print('saturation pressure = ' + str(H2O.Psat))
print('enthalpy of vaporization = ' + str(H2O.dhvap))
print('vapour heat capacity = ' + str(H2O.cpg))
print('liquid heat capacity = ' + str(H2O.cpl))
print('liquid density = ' + str(H2O.density))
print('reduced temperature = ' + str(H2O.Tr))