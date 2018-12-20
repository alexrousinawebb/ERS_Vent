"""
UNIT CONVERSION
===========================================================================
Extremely basic module for converting between standard metric and imperial unit sets
"""

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