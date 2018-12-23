"""
ORDINARY DIFFERENTIAL EQUATIONS ODE
===========================================================================
Module containing reactor system ODEs for solving.
"""

def BPR_PID(t, Y, Urx, Dv, rate, Cv_BPR_max, data):

    T, Tj, nH2O, nH2O2, nO2 = Y

    R = 8.3145 # Ideal gas constant (L*kPa/mol*K)

    MH2O = 18.01528 # Molecular weight water (g/mol)
    MH2O2 = 34.0147 # Molecular weight hydrogen peroxide (g/mol)
    MO2 = 31.998 # Molecular weight oxygen (g/mol)

    [pH2O, pH2O2] = densityL(T) # Liquid density (kg/L)

    [PH2Osat, PH2O2sat] = antoine(T) # Saturation pressures for volatile components (kPa)

    [CpH2OL, CpH2O2L] = CpLiq(T) # Liquid constant pressure heat capacities (J/(g*K))

    [CpH2OG, CpH2O2G, CpO2G] = CpGas(T) # Vapour constant pressure heat capacities (J/(g*K))

    dHvap = enthvap(T) # Enthalpy of vaporization (J/g)

    zH2O = nH2 O /(nH2O + nH2O2 + nO2) # Mole fraction water in vessel

    zH2O2 = nH2O 2 /(nH2O + nH2O2 + nO2) # Mole fraction hydrogen peroxide in vessel

    ntotal = nH2O + nH2O2 + nO2 # Total molar amount in vessel (mol)

    mR = nH2 O *MH2O + nH2O 2 *MH2O2 + nO 2 *MO2 # Total reaction mass (g)

    (xH2O, xH2O2, xO2, yH2O, yH2O2, yO2, nL, nG, P, PH2O, PH2O2, PO2, VL, VG, ZO2) = Equilibrate(ntotal,
                                                                                                 zH2O, zH2O2, T, PH2Osat, PH2O2sat, pH2O, pH2O2, MH2O, MH2O2, nO2, 1, data)

    x = n G *(yH2 O *MH2O + yH2O 2 *MH2O2 + yO 2 *MO2 ) /(n G *(yH2 O *MH2O + yH2O 2 *MH2O2 + yO 2 *MO2) + n L *
                (xH2 O *MH2O + xH2O 2 *MH2O2)) # Mass fraction vapour in vessel

    CpL = xH2 O *CpH2OL + xH2O 2 *CpH2O2L # Average liquid constant pressure heat capacity (J/(g*K))

    CpG = yH2 O *CpH2OG + yH2O 2 *CpH2O2G + yO 2 *CpO2G # Average vapour constant pressure heat capacity (J/(g*K))

    Cp = x* CpG + (1 - x) * CpL  # Average heat capacity in vessel (J/(g*K))

    pG = (1 / (R * C2K(T) * 1000)) * (PH2O * MH2O + PH2O2 * MH2O2 + PO2 * MO2)  # Average vapour density (kg/L)

    pL = pH2O * xH2O + pH2O2 * xH2O2  # Average liquid density (kg/L)

    vL = 1 / pL  # Average specific gravity (L/kg)

    vfg = (1 / pG) - (1 / pL)  # Change in specific volume upon vaporization (L/kg)

    if P_discharge <= P_crit(P, T, yH2O, yH2O2, yO2):
        flag = 1
    else:
        flag = 0

    if P > P_BPR:
        n_vent = ventflow(T, P, P_BPR, yH2O, yH2O2, yO2, Cv_BPR_max, D_BPR, 0)
    else:
        n_vent = 0

    # Differential equations for change in molar amount of components (mol/s)
    dnH2O_dt = krxn(T, kf) * nH2O2 - n_vent * yH2O
    dnH2O2_dt = -krxn(T, kf) * nH2O2 - n_vent * yH2O2
    dnO2_dt = krxn(T, kf) * nH2O2 / 2 - n_vent * yO2

    dTj_dt = rate / 60  # Jacket Temperature (C/s)

    Qr = dnH2O2_dt * MH2O2 * dHrxnL  # Reaction heat flux (J/s)

    QHEx = Urx * Awet(VL, Dv) * (T - Tj)  # External heat exchanger heat flux (J/s)

    Q_vap = ((vL / vfg) + 1) * dHvap * n_vent * (yH2O + yH2O2) * MH2O  # Latent heat of evaporation (J/s)

    Xp = ((dHvap - P * vfg) / (vfg * Cp * 1000)) * (x * dvGdT(PH2O, PH2O2, PO2, T) / 1000 + (1 - x) * dvLdT(T) / 1000)

    dT_dt = (Qr - QHEx - Q_vap) / (mR * Cp * (1 - Xp))

    return [dT_dt, dTj_dt, dnH2O_dt, dnH2O2_dt, dnO2_dt]