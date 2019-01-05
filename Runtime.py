from Scenario import Scenario
import ODE
from pprint import pprint as pprint

T_rxn = 110
VR = 100
XH2O2 = 0.3
mR = 304
t_rxn = 24
D_RD = 12
P_RD = 1000
P_BPR = 200

scen = Scenario(VR, T_rxn, t_rxn, XH2O2, mR, D_RD, P_RD, P_BPR, cooldown_time=2, kf=3000, RD=True, BPR=True,
                TF_vent=True)

ode1 = ODE.ODE(scen)
ode1.initialize_heatup()
ode1.integrate(plot_rt=True)
ode1.initialize_vent(integrator='vode')
ode1.integrate(plot_rt=True)
ode1.plot_vals()

