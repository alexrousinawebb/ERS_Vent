import numpy as np
import Property_Lib as pl
from scipy import optimize as opt
from Scenario import Scenario
import ODE
import ERS
from pprint import pprint

T_rxn = 110
VR = 100
XH2O2 = 0.3
mR = 304
t_rxn = 12
D_RD = 2
P_RD = 1000
P_BPR = 200

scen = Scenario(VR, T_rxn, t_rxn, XH2O2, mR, D_RD, P_RD, P_BPR)

ode1 = ODE.ODE(scen)
ode1.initialize()
ode1.integrate()

