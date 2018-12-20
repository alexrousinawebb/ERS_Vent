import numpy as np
import Property_Lib as pp
from scipy import optimize as opt

T = 25
P = 150
xH2O = 0.5

H2O = pp.Water(T, P, xH2O)
H2O.runmain()

H2O2 = pp.Hydrogen_Peroxide(T, P, xH2O)
H2O2.runmain()

O2 = pp.Oxygen(T, P)
O2.runmain()

pprop = pp.Common_Properties(T)
pprop.runmain()

H2O2 = pp.Hydrogen_Peroxide(T, P, xH2O)
H2O2.runmain()

# print('surface tension = ' + str(pprop.st))
# print('saturation pressure = ' + str(H2O2.Psat))
# print('enthalpy of vaporization = ' + str(pprop.dhvap))
# print('vapour heat capacity = ' + str(H2O2.cpg))
# print('liquid heat capacity = ' + str(H2O2.cpl))
# print('liquid density = ' + str(H2O2.density))
# print('reduced temperature = ' + str(H2O2.Tr))
# print('activity = ' + str(H2O2.gamma))
# print('compressibility = ' + str(H2O2.Z))