from Scenario import Scenario
import ODE

T_rxn = 110
VR = 100
XH2O2 = 0.3
mR = 304
t_rxn = 12
D_RD = 2
P_RD = 1000
P_BPR = 200

scen = Scenario(VR, T_rxn, t_rxn, XH2O2, mR, D_RD, P_RD, P_BPR, plot_rt=True)

ode1 = ODE.ODE(scen)
ode1.initialize()
ode1.integrate()

