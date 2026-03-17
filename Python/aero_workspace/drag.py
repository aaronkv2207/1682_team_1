#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:34:04 2026

@author: brendachow
"""

import numpy as np
from ambiance import Atmosphere
from conceptual_design import V_STALL, cl_max, ureg
import matplotlib.pyplot as plt

# parameters
b = 19.92  # [m]
c = 2.49  # [m]
S = b*c
AR = b**2/S
x_to = 45  # [m]

m = 7500  # [kg]
W = m * 9.81  # [N]

# fuselage (roskam pt 2)
D_f = 1.6  # [m] 
l_f = 15  # [m] f
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
tire_width = 0.284 # [m] = 11.2 [in] 
tire_radius = 0.818/2 # [m] = 32.2 [in] diameter

# trouser_length = ...


# wetted area
S_wet_fuse = np.pi * D_f * l_f * (1 - 2 / lambda_f) ** (2 / 3) * (1 + 1 / lambda_f**2)
S_wet_fancowl = (l_n * D_n * (2 + 0.35 * l_1 * l_n + 0.8 * l_1 * D_hl / l_n * D_n + 1.15 * (1 - l_1 / l_n) * D_ef / D_n))
S_wet_wing = 2*S 
S_wet_htail = 2*S_h 
S_wet_vtail = 2*S_v 

S_wet_gear = tire_width * (np.pi*tire_radius*2)
# s_wet_trousers = ...

# S_wet_total = S_wet_fuse + 8*S_wet_fancowl + S_wet_wing + S_wet_htail + S_wet_vtail + S_wet_gear


# cruise
h_cruise = 5486.4 # [m] = 18000 ft
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


# landing
h_landing = h_takeoff
rho_landing = rho_takeoff
mu_landing = mu_takeoff
V_landing = 20
C_L_landing = 1.665
e_landing = 0.9


print("AR =", AR)
print("wing area =", round(S, 2), "m^2")
# print("total wetted area =", round(S_wet_total, 2), "m^2")
print("takeoff distance =", x_to, "m")
print("MTOW =", round(W, 2), "N")

print("h tail coefficient =", round(V_h, 3))
print("v tail coefficient =", round(V_v, 3))


def calc_Re_l(rho, V, mu, l): # 
    Re_l = rho*V*l/mu
    return Re_l

def calc_C_f(Re_l): # skin friction coefficient
    C_f = 0.455 / (np.log10(Re_l) ** 2.58) # assuming fully turbulent flow for a conservative and realistic estimate!
    return C_f

def calc_CDA(S_wet, C_f, K_f): # drag area
    CDA = S_wet * C_f * K_f
    return CDA

# calculate form factors
def calc_K_f_airfoil(tc):
    K_f = 1 + 2.0 * tc + 60 * tc**4
    return K_f

def calc_K_f_axi(dl):
    K_f = 1 + 1.5 * dl**(3/2) + 7 * dl**3
    return K_f

# calculate profile drag (dissipation summation buildup)
def calc_C_Dp(rho, V, mu):
    V_i = V_inf = V

    S_wets = [S_wet_fuse, S_wet_fancowl, S_wet_wing, S_wet_htail, S_wet_vtail, S_wet_gear]
    ls = [l_f, l_n, c, c_h, c_v, tire_radius*2]
    fineness_ratios = [D_f/l_f, 0.12, 0.12, 0.10, 0.10, 0.10] # need to check fineness for landing gear

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

        CDA = calc_CDA(s_val, C_fs[i], K_fs[i])

        if i==1: # 8 total fans
            CDA*=8

        CDAs.append(CDA)

    Dp = sum(CDAs) * (1/2 * rho * V_i**3 / V_inf) #NOTE: Would need to change this logic if you V_i != V_inf
    C_Dp = Dp / (1/2 * rho * V**2 * S)

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

    return Dp, C_Dp


# calculate induced drag
def calc_C_Di(C_L, rho, V, e):

    C_Di = C_L**2 / (np.pi * AR * e)
    Di = C_Di * 1/2 * rho * V**2 * S

    return Di, C_Di

print(C_L_cruise)

# calcualte drag during cruise!
Dp_cruise, C_Dp_cruise = calc_C_Dp(rho_cruise, V_cruise, mu_cruise)
Di_cruise, C_Di_cruise = calc_C_Di(C_L_cruise, rho_cruise, V_cruise, e_cruise)

print('-------')
print('cruising :)')
print('profile drag coefficient =', round(C_Dp_cruise,3))
print('induced drag coefficient =', round(C_Di_cruise,3))
print('total drag coefficient =', round(C_Dp_cruise+C_Di_cruise,3))
         

# calculate drag during takeoff!
Dp_takeoff, C_Dp_takeoff = calc_C_Dp(rho_takeoff, V_takeoff, mu_takeoff)
Di_takeoff, C_Di_takeoff = calc_C_Di(C_L_takeoff, rho_takeoff, V_takeoff, e_takeoff)

print('-------')
print('takeoff!')
print('profile drag coefficient =', round(C_Dp_takeoff,3))
print('induced drag coefficient =', round(C_Di_takeoff,3))
print('total drag coefficient =', round(C_Dp_takeoff+C_Di_takeoff,3))


# calculate drag during landing!
Dp_landing, C_Dp_landing = calc_C_Dp(rho_landing, V_landing, mu_landing)
Di_landing, C_Di_landing = calc_C_Di(C_L_landing, rho_landing, V_landing, e_landing)

# print('-------')
# print('landinggg')
# print('profile drag coefficient =', round(C_Dp_landing,3))
# print('induced drag coefficient =', round(C_Di_landing,3))
# print('total drag coefficient =', round(C_Dp_landing+C_Di_landing,3))

# drag polar??
C_D = np.linspace(0,2,1000)
C_L = np.sqrt((C_D - C_Dp_cruise) * np.pi * AR * e_cruise)
C_L2 = np.sqrt((C_D - C_Dp_takeoff) * np.pi * AR * e_takeoff)

plt.figure()
plt.plot(C_D, C_L, label = 'cruise')
plt.plot(C_D, C_L2, label = 'takeoff')
plt.xlabel('$C_D$')
plt.ylabel('$C_L$')
plt.legend()
plt.show()




# print(C_L_cruise/(C_Dp_cruise+C_Di_cruise))



