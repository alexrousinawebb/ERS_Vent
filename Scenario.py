"""
SCENARIO
===========================================================================
Module initiates and stores a valid ERS scenario.
"""

import numpy as np
from Conversion import g2l

class Scenario():
    def __init__(self, reactor_volume, set_temp, reaction_time, H2O2_massfraction, charge_mass, D_RD=None,
                 P_RD=None, P_BPR=None, D_BPR=0.5, BPR_max_Cv=5.5, start_pressure=101.325, kf=1, start_temp=25,
                 heat_transfer_coefficient=450, aspect_ratio=1.5, MAWP=100000, max_rate_jacket=2, integrator='lsoda'):
        """
            Initializes scenario instance.

            Arguments (Required):
            reactor_volume:             Total reactor volume in gallons
            set_temp:                   Reaction temperature in degrees celsius
            reaction_time:              Total reaction time in hours
            H2O2_massfraction:          Mass fraction of hydrogen peroxide in starting mixture
            charge_mass:                Total reaction mass charged to reactor in kilograms
            D_RD:                       Rupture disk diameter in inches
            P_RD:                       Rupture disk/PRV activation pressure in kilopascals
            P_BPR:                      Backpressure regulator set point pressure in kilopascals

            Arguments (Optional):
            D_BPR:                      Backpressure regulator orifice diameter in inches
            BPR_max_Cv:                 Max discharge coefficient for backpressure regulator
            start_pressure:             Initial pressure in reactor headspace at start of reaction period in kilopascals
            kf:                         Contamination factor for hydrogen peroxide (1 - 10,000)
            start_temp:                 Starting temperature, default at ambient conditions
            heat_transfer_coefficeint:  Reactor heat transfer coefficient in (W/m**2 K)
            aspect_ratio:               Reactor aspect ratio (height/width ratio)
            max_rate_jacket:            Maximum rate of cooling (temperature change) of jacket fluid in degrees celsius
                                        per second
            integrator:                 Choice of integration method. Options include - vode, zvode, lsoda, dopri5,
                                                                                        dop853
            MAWP:                       Maximum allowable working pressure for reaction vessel in kilopascals
                                        Defines where to stop integration in the event of an overpressure scenario
        """

        #  Reactor Parameters
        self.VR = g2l(reactor_volume)
        self.Ux = heat_transfer_coefficient
        self.AR = aspect_ratio
        self.D = 2 * ((self.VR * 0.001) / (2 * np.pi * self.AR)) ** (1 / 3)
        self.h = self.D * self.AR
        self.MAWP = MAWP

        #  ERS Design Parameters
        self.D_RD = D_RD
        self.D_BPR = D_BPR
        self.BPR_max_Cv = BPR_max_Cv
        self.P_RD = P_RD
        self.P_BPR = P_BPR

        #  Chemical Properties
        self.XH2O2 = H2O2_massfraction
        self.mR = charge_mass
        self.kf = kf

        #  Reaction Conditions
        self.rxn_temp = set_temp
        self.rxn_time = reaction_time
        self.T0 = start_temp
        self.P0 = start_pressure
        self.max_rate = max_rate_jacket

        #  Miscellaneous
        self.integrator = integrator

