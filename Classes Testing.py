class PhysProp:
    def __init__(self):
        self.surface_tension = None
        self.pH2O = None
        self.pH2O2 = None

    def st(self, T):
        """
            Calculate surface tension of liquid water in N/m.

            Arguments:
            T: Temperature of liquid in degrees Celcius
        """

        B = 235.8E-3
        b = -0.625
        u = 1.256
        Tc = 647.15

        self.surface_tension = B * (((Tc - T) / Tc) ** u) * (1 + b * (Tc - T) / Tc)  # returned as (N/m)

    def p_L(self, T):
        """
            Calculate density of liquid water and hydrogen peroxide in kg/L.

            Arguments:
            T: Temperature of liquid in degrees Celcius
        """

        A = 999.83952
        B = 16.945176
        C = -7.987040e-3
        D = -46.170461e-6
        E = 105.56302e-9
        F = -280.54253e-12
        G = 16.897850e-3

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

        self.pH2O = ((A + B * T + C * T ** 2 + D * T ** 3 + E * T ** 4 + F * T ** 5) / (1 + G * T)) / 1000

        N = Jb + Kb * T + Lb * (T ** 2) + Mb * (T ** 3)
        O = Jc + Kc * T + Lc * (T ** 2) + Mc * (T ** 3)
        P = Jd + Kd * T + Ld * (T ** 2) + Md * (T ** 3)

        if T >= 100:
            self.pH2O2 = 1.2456174226244978
        else:
            self.pH2O2 = self.pH2O + N + O ** 2 + P ** 3

    def antoine(self, ):

T = 398.15

pp = PhysProp()
pp.st(T)
pp.p_L(T)

sigma = pp.surface_tension
pH2O = pp.pH2O
pH2O2 = pp.pH2O2
