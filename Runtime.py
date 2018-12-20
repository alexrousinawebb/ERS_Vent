import numpy as np
import Property_Lib as pl
from scipy import optimize as opt
from Scenario import Scenario
import VLE

T = 25
VR = 100
XH2O2 = 0.3
mR =304

scen = Scenario(VR)

ic = VLE.Initial_Conditions(T, XH2O2, mR, scen)
ic.initial_conditions()