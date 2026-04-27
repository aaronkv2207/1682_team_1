from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from ambiance import Atmosphere
from conceptual_design import ureg

# TODO: Check new plane inputs (for each plane) in aero_main and update drag build up.
# I started roughly changing some things

# parameters
b = 17.75  # [m]
c = 2.54  # [m]
S = b * c
AR = b**2 / S
# x_to = 45  # [m]

m = 7500  # [kg]
W = m * 9.81  # [N]

# fuselage (roskam pt 2)
D_f = 1.6  # [m]
l_f = 15  # [m] f
lambda_f = l_f / D_f


# nacelles (roskam pt 2)
l_n = c / 2  # [m] length of cowl (est)
D_n = 1.1 * 1.116  # [m] diameter of cowl (est)
l_1 = l_n / 2
D_hl = 1.116  # [m] diameter of fan (from prop team)
D_ef = D_hl


# tail
S_h = 13.33  # [m^2]
c_h = 2.11  # [m]

S_v = 5.99  # [m^2]
c_v = 2.23  # [m]

l_v = 8.0  # [m]
vt_ar = 1.2
vt_c = np.sqrt(S_v / vt_ar)
l_h = l_v + 0.25 * vt_c  # [m]

V_h = S_h * l_h / (S * c)
V_v = S_v * l_v / (S * b)


# landing gear (still need to add!!)
tire_width = 0.284  # [m] = 11.2 [in]
tire_radius = 0.818 / 2  # [m] = 32.2 [in] diameter

# trouser_length = ...


# # TODO: strut
# c_s = ... # chord of strut
# l_s = ... # length of strut
# S_s = c_s * l_s


# wetted area
S_wet_fuse = np.pi * D_f * l_f * (1 - 2 / lambda_f) ** (2 / 3) * (1 + 1 / lambda_f**2)
S_wet_fancowl = (
    l_n
    * D_n
    * (
        2
        + 0.35 * l_1 * l_n
        + 0.8 * l_1 * D_hl / l_n * D_n
        + 1.15 * (1 - l_1 / l_n) * D_ef / D_n
    )
)
S_wet_wing = 2 * S
S_wet_htail = 2 * S_h
S_wet_vtail = 2 * S_v

S_wet_gear = tire_width * (np.pi * tire_radius * 2)
# S_wet_strut = 2*S_s #TODO: Implement
# s_wet_trousers = ...

# S_wet_total = S_wet_fuse + 8*S_wet_fancowl + S_wet_wing + S_wet_htail + S_wet_vtail + S_wet_gear


# cruise
h_cruise = 5486.4  # [m] = 18000 ft
rho_cruise = Atmosphere(h=h_cruise).density[0]  # [kg/m^3]
# print('air density cruise', rho_cruise)
mu_cruise = Atmosphere(h=h_cruise).dynamic_viscosity[0]  # [Pa * s]
print("dynamic viscosity cruise", mu_cruise)
V_cruise = 125  # [m/s]
C_L_cruise = W / (1 / 2 * rho_cruise * V_cruise**2 * S)
print("C_L cruise", C_L_cruise)
e_cruise = 0.95  # NOTE: should be using jvl e calculated in aero_main


# takeoff
h_takeoff = 0  # [m]
rho_takeoff = Atmosphere(h=h_takeoff).density[0]
# print('air density takeoff', rho_takeoff)
mu_takeoff = Atmosphere(h=h_takeoff).dynamic_viscosity[
    0
]  # [Ns/m^2] dynamic viscosity at sea level
# print('dynamic viscosity takeoff', mu_takeoff)
# # print('takeoff velocity', V_takeoff)
# C_L_takeoff = cl_max  # from wt data???
# e_takeoff = 0.9

# print(C_L_takeoff)


# landing
h_landing = h_takeoff
rho_landing = rho_takeoff
mu_landing = mu_takeoff
# V_landing = 20 # NOTE: should be using jvl velocity in aero_main
# C_L_landing = 5 # NOTE: should be using jvl CL in aero_main
# e_landing = 0.9 # NOTE: should be using jvl e in aero_main


def calc_Re_l(rho, V, mu, l):  #
    Re_l = rho * V * l / mu
    return Re_l


def calc_C_f(Re_l):  # skin friction coefficient
    C_f = 0.455 / (
        np.log10(Re_l) ** 2.58
    )  # assuming fully turbulent flow for a conservative and realistic estimate!
    return C_f


def calc_CDA(S_wet, C_f, K_f):  # drag area
    CDA = S_wet * C_f * K_f
    return CDA


# calculate form factors
def calc_K_f_airfoil(tc):
    K_f = 1 + 2.0 * tc + 60 * tc**4
    return K_f


def calc_K_f_axi(dl):
    K_f = 1 + 1.5 * dl ** (3 / 2) + 7 * dl**3
    return K_f


# calculate profile drag (dissipation summation buildup)
def calc_C_Dp(rho, V, mu):
    V_i = V_inf = V

    # S_wets = [S_wet_fuse, S_wet_fancowl, S_wet_wing, S_wet_htail, S_wet_vtail, S_wet_gear, S_wet_strut] # TODO: implement S_wet_strut
    # ls = [l_f, l_n, c, c_h, c_v, tire_radius*2, D_s] # TODO: implement S_wet_strut
    S_wets = [
        S_wet_fuse,
        S_wet_fancowl,
        S_wet_wing,
        S_wet_htail,
        S_wet_vtail,
        S_wet_gear,
    ]
    ls = [l_f, l_n, c, c_h, c_v, tire_radius * 2]
    fineness_ratios = [
        D_f / l_f,
        D_n / (2 * l_n),
        0.12,
        0.10,
        0.10,
        0.10,
        0.25,
    ]  # need to check fineness for landing gear

    Re_ls = []
    C_fs = []
    CDAs = []
    K_fs = []

    for i, s_val in enumerate(S_wets):
        Re_l = calc_Re_l(rho, V, mu, ls[i])
        Re_ls.append(Re_l)

        C_f = calc_C_f(Re_l)
        C_fs.append(C_f)

        if i <= 1:
            K_f = calc_K_f_axi(fineness_ratios[i])
            K_fs.append(K_f)

        else:
            K_f = calc_K_f_airfoil(fineness_ratios[i])
            K_fs.append(K_f)

        # K_f = calc_K_f_airfoil(fineness_ratios[i])
        # K_fs.append(K_f)

        CDA = calc_CDA(s_val, C_fs[i], K_fs[i])

        if i == 1:  # 8 total fans
            CDA *= 8

        CDAs.append(CDA)

    Dps = []
    C_Dps = []
    for CDA in CDAs:
        Dp = CDA * (1 / 2 * rho * V_i**3 / V_inf)
        Dps.append(Dp)
        C_Dp = Dp / (1 / 2 * rho * V**2 * S)
        C_Dps.append(C_Dp)

    C_Dps = np.array(C_Dps)
    print(np.round(C_Dps, 3))

    # Dp = sum(CDAs) * (1/2 * rho * V_i**3 / V_inf) #NOTE: Would need to change this logic if you V_i != V_inf
    # C_Dp = Dp / (1/2 * rho * V**2 * S)

    # # printing all the parts for debugging...
    # Re_ls = np.array(Re_ls)
    # print(np.round(Re_ls,2))

    # C_fs = np.array(C_fs)
    # print(np.round(C_fs,4))

    # normalized_S = []

    # for s in S_wets:
    #     normalized_S.append(s/S)

    # normalized_S = np.array(normalized_S)
    # print(np.round(normalized_S,2))

    return sum(Dps), sum(C_Dps)


# calculate induced drag
def calc_C_Di(C_L, rho, V, e):
    C_Di = C_L**2 / (np.pi * AR * e)
    Di = C_Di * 1 / 2 * rho * V**2 * S

    return Di, C_Di

