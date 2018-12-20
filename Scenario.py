"""
SCENARIO
===========================================================================
Module initiates and stores a valid ERS scenario.
"""

import numpy as np
from Conversion import g2l

class Scenario():
    def __init__(self, reactor_volume, heat_transfer_coefficient=450, aspect_ratio=1.5):
        """
            Initializes scenario instance.

            Arguments:
            T: Temperature of liquid in degrees Celsius
            P: Pressure in headspace of reactor in kilopascals
            xH2O: Mole fraction of water in liquid phase
        """
        self.VR = g2l(reactor_volume)
        self.Ux = heat_transfer_coefficient
        self.AR = aspect_ratio
        self.D = 2 * ((self.VR * 0.001) / (2 * np.pi * self.AR)) ** (1 / 3)
        self.h = self.D * self.AR