#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:34:04 2026

@author: brendachow
"""

import numpy as np
from ambiance import Atmosphere
from conceptual_design import V_STALL, cl_max, ureg

# parameters
AR = 8 
b = 19.92  # [m]
c = 2.49  # [m]
S = b**2 / AR
x_to = 45  # [m]

m = 7500  # [kg]
W = m * 9.81  # [N]

# fuselage (roskam pt 2)
D_f = 1.6  # [m] 
l_f = 15  # [m] 
lambda_f = l_f / D_f


# nacelles (roskam pt 2)
l_n = c/2  # [m] length of cowl (est)
D_n = 1.1*1.116  # [m] diameter of cowl (est)
l_1 = l_n / 2
D_hl = 1.116  # [m] diameter of fan (from prop team)
D_ef = D_hl


# tail
S_h = 11.98  # [m^2]
c_h = 2.45 # [m]

S_v = 7.61  # [m^2]
c_v = 3.04 # [m]

l_h = 9.06  # [m] estimated from twin otter
l_v = 9.06  # [m] estimated from twin otter

V_h = S_h * l_h / (S * c)
V_v = S_v * l_v / (S * b)


# landing gear (still need to add!!)
tire_width = ...
tire_radius = ...
trouser_length = ...


# wetted area
S_wet_fuse = np.pi * D_f * l_f * (1 - 2 / lambda_f) ** (2 / 3) * (1 + 1 / lambda_f**2)
S_wet_fancowl = (l_n * D_n * (2 + 0.35 * l_1 * l_n + 0.8 * l_1 * D_hl / l_n * D_n + 1.15 * (1 - l_1 / l_n) * D_ef / D_n))
S_wet_wing = S * 2
S_wet_htail = S_h * 2
S_wet_vtail = S_v * 2
S_wet_gear = ...

S_wet_total = S_wet_fuse + 8*S_wet_fancowl + S_wet_wing + S_wet_htail + S_wet_vtail


# cruise
h_cruise = 5486.4 # [m]
rho_cruise = Atmosphere(h=h_cruise).density[0]  # [kg/m^2]
mu_cruise = Atmosphere(h=h_cruise).dynamic_viscosity[0] # [Pa * s]
V_cruise = 125  # [m/s]
C_L_cruise = W / (1 / 2 * rho_cruise * V_cruise**2 * S)
e_cruise = 0.95


# takeoff
h_takeoff = 0 # [m]
rho_takeoff = Atmosphere(h=h_takeoff).density[0]
mu_takeoff = Atmosphere(h=h_takeoff).dynamic_viscosity[0]  # [Ns/m^2] dynamic viscosity at sea level
V_takeoff = 1.1 * V_STALL.magnitude  # [m/s]
C_L_takeoff = cl_max  # from wt data???
e_takeoff = 0.9


print("AR =", AR)
print("wing area =", round(S, 2), "m^2")
print("total wetted area =", round(S_wet_total, 2), "m^2")
print("takeoff distance =", x_to, "m")
print("MTOW =", round(W, 2), "N")

print("h tail coefficient =", round(V_h, 3))
print("v tail coefficient =", round(V_v, 3))


def calc_Re_l(rho, V, mu, l):
    Re_l = rho*V*l/mu
    return Re_l

def calc_C_f(Re_l):    
    C_f = 0.455 / (np.log10(Re_l ** 2.58)) # assuming fully turbulent flow for a conservative and realistic estimate!
    return C_f

def calc_CDA(S_wet, C_f, K_f):
    CDA = S_wet * C_f * K_f
    return CDA

def calc_K_f_airfoil(tc):
    K_f = 1 + 2.0 * tc + 60 * tc**4
    return K_f

def calc_K_f_axi(dl):
    K_f = 1 + 1.5 * dl**(3/2) + 7 * dl**3
    return K_f

# calculate profile drag (dissipation summation buildup)
def calc_C_Dp(rho, V, mu):
    V_i = V_inf = V


    S_wets = [S_wet_fuse, S_wet_fancowl, S_wet_wing, S_wet_htail, S_wet_vtail]
    ls = [l_f, l_n, c, c_h, c_v]
    fineness_ratios = [D_f/l_f, 0.12, 0.12, 0.10, 0.10]

    Re_ls = []
    C_fs = []
    CDAs = []
    K_fs = []

    for i, s_val in enumerate(S_wets):
        Re_l = calc_Re_l(rho, V, mu, ls[i])
        Re_ls.append(Re_l)

        C_f = calc_C_f(Re_l)
        C_fs.append(C_f)

        if i == 0:
            K_f = calc_K_f_axi(fineness_ratios[i])
            K_fs.append(K_f)

        else: 
            K_f = calc_K_f_airfoil(fineness_ratios[i])
            K_fs.append(K_f)

        CDA = calc_CDA(S_wets[i], C_fs[i], K_fs[i])

        if i==1: # 8 total fans
            CDA*=8

        CDAs.append(CDA)
    
    # print(CDAs)

    Dp = sum(CDAs) * (1/2 * rho * V_i**3 / V_inf) #NOTE: Would need to change this logic if you V_i != V_inf
    C_Dp = Dp / (1/2 * rho * V**2 * S)

    return Dp, C_Dp


# calculate induced drag
def calc_C_Di(C_L, rho, V, e):

    C_Di = C_L**2 / (np.pi * AR * e)
    Di = C_Di * 1/2 * rho * V**2 * S

    return Di, C_Di


Dp_cruise, C_Dp_cruise = calc_C_Dp(rho_cruise, V_cruise, mu_cruise)
Dp_takeoff, C_Dp_takeoff = calc_C_Dp(rho_takeoff, V_takeoff, mu_takeoff)

Di_cruise, C_Di_cruise = calc_C_Di(C_L_cruise, rho_cruise, V_cruise, e_cruise)
Di_takeoff, C_Di_takeoff = calc_C_Di(C_L_takeoff, rho_takeoff, V_takeoff, e_takeoff)

print('-------')
print('cruising :)')
print('profile drag coefficient =', round(C_Dp_cruise,3))
print('induced drag coefficient =', round(C_Di_cruise,3))
print('total drag coefficient =', round(C_Dp_cruise+C_Di_cruise,3))

print('-------')
print('takeoff!')
print('profile drag coefficient =', round(C_Dp_takeoff,3))
print('induced drag coefficient =', round(C_Di_takeoff,3))
print('total drag coefficient =', round(C_Dp_takeoff+C_Di_takeoff,3))
