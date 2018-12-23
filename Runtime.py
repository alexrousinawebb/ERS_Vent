import numpy as np
import Property_Lib as pl
from scipy import optimize as opt
from Scenario import Scenario
import VLE
from pprint import pprint

T = 90
VR = 100
XH2O2 = 0.3
mR = 304

scen = Scenario(VR)

ic = VLE.Initial_Conditions(T, XH2O2, mR, scen)
ic.initial_conditions()

H2O = pl.Water(T, ic.P, ic)
H2O2 = pl.Hydrogen_Peroxide(T, ic.P, ic)
O2 = pl.Oxygen(T, ic.P, ic)
cp = pl.Common_Properties(T, ic, H2O, H2O2, O2)

print('H2O:')
pprint(vars(H2O))
print('H2O2:')
pprint(vars(H2O2))
print('O2:')
pprint(vars(O2))
print('Common Properties:')
pprint(vars(cp))