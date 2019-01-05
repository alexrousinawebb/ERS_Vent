"""
UNIT CONVERSION
===========================================================================
Basic module for converting between standard metric and imperial unit sets
"""

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

def A_relief(D):
    """
        Convert circular diameter in inches to area in square meters.
    """
    return np.pi * ((D * 0.0254) ** 2) / 4

def A_wet(VL, Dv):
    """
        Calculate wetted reactor area using surface area of cylinder approximation in square meters.
    """
    hwet = (VL * 0.001) / (np.pi * (Dv / 2) ** 2)

    return np.pi * (Dv / 2) ** 2 + 2 * np.pi * (Dv / 2) * hwet

def cv2kd(Cv, D):
    return Cv / (27.66 * A_relief(D) * 1550)