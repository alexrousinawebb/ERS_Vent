import numpy as np
import Property_Lib as pl
from scipy import optimize as opt
from Scenario import Scenario
import ODE
from pprint import pprint

T_rxn = 110
VR = 100
XH2O2 = 0.3
mR = 304
t_rxn = 6
D_RD = 2
P_RD = 1000
P_BPR = 200

scen = Scenario(VR, T_rxn, t_rxn, XH2O2, mR, D_RD, P_RD, P_BPR)

ode1 = ODE.ODE(scen)

ic = VLE.Initial_Conditions(scen)
ic.initial_conditions()

H2O = pl.Water(scen.T0, ic.P, ic)
H2O2 = pl.Hydrogen_Peroxide(scen.T0, ic.P, ic)
O2 = pl.Oxygen(scen.T0, ic.P, ic)
cp = pl.Common_Properties(scen.T0, ic, H2O, H2O2, O2)

print('H2O:')
print(np.asarray(vars(H2O)))
print('H2O2:')
pprint(vars(H2O2))
print('O2:')
pprint(vars(O2))
print('Common Properties:')
pprint(vars(cp))