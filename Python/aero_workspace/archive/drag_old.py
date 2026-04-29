#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:34:04 2026

@author: brendachow
"""

import matplotlib.pyplot as plt
import numpy as np
from ambiance import Atmosphere
from conceptual_design import ureg

# TODO: Check new plane inputs (for each plane) in aero_main and update drag build up.
# I started roughly changing some things

# parameters
b = 17.75  # [m]
c = 2.54   # [m]
S = b*c
AR = b**2/S
# x_to = 45  # [m]

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
S_h = 13.33  # [m^2]
c_h = 2.11 # [m]

S_v = 5.99  # [m^2]
c_v = 2.23 # [m]

l_v = 8.0  # [m]
vt_ar = 1.2
vt_c = np.sqrt(S_v / vt_ar)
l_h = (l_v + 0.25 * vt_c) # [m]

V_h = S_h * l_h / (S * c)
V_v = S_v * l_v / (S * b)


# landing gear (still need to add!!)
tire_width = 0.284 # [m] = 11.2 [in]
tire_radius = 0.818/2 # [m] = 32.2 [in] diameter

# trouser_length = ...


# # TODO: strut
# c_s = ... # chord of strut
# l_s = ... # length of strut
# S_s = c_s * l_s


# wetted area
S_wet_fuse = np.pi * D_f * l_f * (1 - 2 / lambda_f) ** (2 / 3) * (1 + 1 / lambda_f**2)
S_wet_fancowl = (l_n * D_n * (2 + 0.35 * l_1 * l_n + 0.8 * l_1 * D_hl / l_n * D_n + 1.15 * (1 - l_1 / l_n) * D_ef / D_n))
S_wet_wing = 2*S
S_wet_htail = 2*S_h
S_wet_vtail = 2*S_v

S_wet_gear = tire_width * (np.pi*tire_radius*2)
# S_wet_strut = 2*S_s #TODO: Implement
# s_wet_trousers = ...

# S_wet_total = S_wet_fuse + 8*S_wet_fancowl + S_wet_wing + S_wet_htail + S_wet_vtail + S_wet_gear


# cruise
h_cruise = 5486.4 # [m] = 18000 ft
rho_cruise = Atmosphere(h=h_cruise).density[0]  # [kg/m^3]
# print('air density cruise', rho_cruise)
mu_cruise = Atmosphere(h=h_cruise).dynamic_viscosity[0] # [Pa * s]
print('dynamic viscosity cruise', mu_cruise)
V_cruise = 125  # [m/s]
C_L_cruise = W / (1 / 2 * rho_cruise * V_cruise**2 * S)
print('C_L cruise', C_L_cruise)
e_cruise = 0.95 # NOTE: should be using jvl e calculated in aero_main


# takeoff
h_takeoff = 0 # [m]
rho_takeoff = Atmosphere(h=h_takeoff).density[0]
# print('air density takeoff', rho_takeoff)
mu_takeoff = Atmosphere(h=h_takeoff).dynamic_viscosity[0]  # [Ns/m^2] dynamic viscosity at sea level
# print('dynamic viscosity takeoff', mu_takeoff)
# # print('takeoff velocity', V_takeoff)
# C_L_takeoff = cl_max  # from wt data???
# e_takeoff = 0.9

# print(C_L_takeoff)


# landing
h_landing = h_takeoff
rho_landing = rho_takeoff
mu_landing = mu_takeoff
V_landing = 20 # NOTE: should be using jvl velocity in aero_main
C_L_landing = 5 # NOTE: should be using jvl CL in aero_main
e_landing = 0.9 # NOTE: should be using jvl e in aero_main


print("AR =", AR)
print("wing area =", round(S, 2), "m^2")
# print("total wetted area =", round(S_wet_total, 2), "m^2")
# print("takeoff distance =", x_to, "m")
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

    # S_wets = [S_wet_fuse, S_wet_fancowl, S_wet_wing, S_wet_htail, S_wet_vtail, S_wet_gear, S_wet_strut] # TODO: implement S_wet_strut
    # ls = [l_f, l_n, c, c_h, c_v, tire_radius*2, D_s] # TODO: implement S_wet_strut
    S_wets = [S_wet_fuse, S_wet_fancowl, S_wet_wing, S_wet_htail, S_wet_vtail, S_wet_gear]
    ls = [l_f, l_n, c, c_h, c_v, tire_radius*2]
    fineness_ratios = [D_f/l_f, D_n/(2*l_n), 0.12, 0.10, 0.10, 0.10, 0.25] # need to check fineness for landing gear

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

        if i==1: # 8 total fans
            CDA*=8

        CDAs.append(CDA)

    Dps = []
    C_Dps = []
    for CDA in CDAs:
        Dp = CDA * (1/2* rho * V_i**3 / V_inf)
        Dps.append(Dp)
        C_Dp = Dp / (1/2 * rho * V**2 * S)
        C_Dps.append(C_Dp)

    C_Dps = np.array(C_Dps)
    print(np.round(C_Dps,3))

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
    Di = C_Di * 1/2 * rho * V**2 * S

    return Di, C_Di

# print(C_L_cruise)

# calcualte drag during cruise!
Dp_cruise, C_Dp_cruise = calc_C_Dp(rho_cruise, V_cruise, mu_cruise)
Di_cruise, C_Di_cruise = calc_C_Di(C_L_cruise, rho_cruise, V_cruise, e_cruise)

print('-------')
print('cruising :)')
print('profile drag coefficient =', round(C_Dp_cruise,3))
print('induced drag coefficient =', round(C_Di_cruise,3))
print('total drag coefficient =', round(C_Dp_cruise+C_Di_cruise,3))


# # calculate drag during takeoff!
# Dp_takeoff, C_Dp_takeoff = calc_C_Dp(rho_takeoff, V_takeoff, mu_takeoff)
# Di_takeoff, C_Di_takeoff = calc_C_Di(C_L_takeoff, rho_takeoff, V_takeoff, e_takeoff)

# print('-------')
# print('takeoff!')
# print('profile drag coefficient =', round(C_Dp_takeoff,3))
# print('induced drag coefficient =', round(C_Di_takeoff,3))
# print('total drag coefficient =', round(C_Dp_takeoff+C_Di_takeoff,3))


# calculate drag during landing!
Dp_landing, C_Dp_landing = calc_C_Dp(rho_landing, V_landing, mu_landing)
Di_landing, C_Di_landing = calc_C_Di(C_L_landing, rho_landing, V_landing, e_landing)

print('-------')
print('landinggg')
print('profile drag coefficient =', round(C_Dp_landing,3))
print('induced drag coefficient =', round(C_Di_landing,3))
print('total drag coefficient =', round(C_Dp_landing+C_Di_landing,3))



print('cruise L/D', C_L_cruise/(C_Dp_cruise+C_Di_cruise))
# print(C_L_takeoff/(C_Dp_takeoff+C_Di_takeoff))


# # sensitivity plots! (for design decision memo)

# # wing

# # profile drag coefficient vs t/c

# tc_wing = np.linspace(0,1,1000)
# Re_wing = rho_cruise*V_cruise*c/mu_cruise
# cd_wing = 2*0.455/(np.log10(Re_wing))**2.58 * (1+2*(tc_wing)+60*(tc_wing)**4)

# x = 0.12
# y = 2*0.455/(np.log10(Re_wing))**2.58 * (1+2*(x)+60*(x)**4)

# plt.figure()
# plt.plot(tc_wing, cd_wing)
# plt.plot(x, y, 'o')
# plt.text(x+0.01, y+0.01, r'$t/c = 0.12$')
# plt.xlabel(r'Thickness/Chord Ratio, $t/c$')
# plt.ylabel(r'$C_{D_0}$')
# plt.show()


# # # fuselage
# fineness = np.linspace(0,1,1000)
# S_wetted = np.pi*fineness*l_f**2 * (1-2*fineness)**(2/3) * (1+fineness**2)
# form_factor = calc_K_f_axi(fineness)

# reynolds = calc_Re_l(rho_cruise, V_cruise, mu_cruise, l_f)
# skin_friction = calc_C_f(reynolds)

# cd_fuse = S_wetted*form_factor*skin_friction / S

# plt.figure()
# plt.plot(fineness, cd_fuse)
# plt.xlabel('fineness ratio')
# plt.show()


# pie charts

# plt.figure()
# vals_cruise = [0.003, 0.006, 0.007, 0.003]
# components = ['Fuselage', 'Nacelles', 'Wing', 'Tail']
# plt.pie(vals_cruise, labels = components, autopct='%1.1f%%')
# plt.title('Cruise')
# plt.show()

# plt.figure()
# vals_TO = [0.004, 0.008, 0.009, 0.003]
# plt.pie(vals_TO, labels = components)
# plt.title('Takeoff')
# plt.show()


# # nacelles!

# dl_nacelles = np.linspace(0,1,1000)
# Re_nacelles = rho_cruise*V_cruise*l_n/mu_cruise
# cd_nacelles = 1/S * (dl_nacelles*l_n**2 * (2+0.175*l_n**2+0.4*dl_nacelles*D_hl*l_n+0.575*D_hl/l_n*1/dl_nacelles) * 0.455/(np.log10(Re_nacelles))**2.58 * (1+1.5*dl_nacelles**(3/2)+7*dl_nacelles**3))

# x = 0.986
# y = 1/S * (x*l_n**2 * (2+0.175*l_n**2+0.4*x*D_hl*l_n+0.575*D_hl/l_n*1/x) * 0.455/(np.log10(Re_nacelles))**2.58 * (1+1.5*x**(3/2)+7*x**3))

# plt.figure()
# plt.plot(dl_nacelles, cd_nacelles)
# plt.plot(x, y, 'o')
# plt.text(x-0.2, y, r'$d/\ell = 0.986$')
# plt.xlabel(r'Diameter/Length Ratio, $d/\ell$')
# plt.ylabel(r'$C_{D_0}$')
# plt.show()


# # fuselage

# dl_fuse = np.linspace(0,1,1000)
# Re_fuse = rho_cruise*V_cruise*l_f/mu_cruise
# cd_fuse = 1/S * (np.pi*D_f*l_f)*(1-2*dl_fuse)**(2/3)*(1+dl_fuse**2) * (1+1.5*dl_fuse**(3/2)+7*dl_fuse**3) * 0.455/(np.log10(Re_fuse))**2.58
# x = 0.107
# y = 1/S * (np.pi*D_f*l_f)*(1-2*x)**(2/3)*(1+x**2) * (1+1.5*x**(3/2)+7*x**3) * 0.455/(np.log10(Re_fuse))**2.58

# plt.figure()
# plt.plot(dl_fuse, cd_fuse)
# plt.plot(x, y, 'o')
# plt.text(x-0.04, y-0.0003, r'$d/\ell = 0.107$')
# plt.xlabel(r'Diameter/Length Ratio, $d/\ell$')
# plt.ylabel(r'$C_{D_0}$')
# plt.show()


# # tail

# tc_tail = np.linspace(0,1,1000)

# Re_htail = rho_cruise*V_cruise*c_h/mu_cruise
# cd_htail = 1/S * 2*S_h * (1+2*tc_tail+60*tc_tail**4) * 0.455/(np.log10(Re_htail))**2.58

# Re_vtail = rho_cruise*V_cruise*c_v/mu_cruise
# cd_vtail = 1/S * 2*S_v * (1+2*tc_tail+60*tc_tail**4) * 0.455/(np.log10(Re_vtail))**2.58

# x = 0.12
# y_h = 1/S * 2*S_h * (1+2*x+60*x**4) * 0.455/(np.log10(Re_htail))**2.58
# y_v = 1/S * 2*S_h * (1+2*x+60*x**4) * 0.455/(np.log10(Re_vtail))**2.58

# plt.figure()
# plt.plot(tc_tail, cd_htail+cd_vtail)
# plt.plot(x, y_h+y_v, 'o')
# plt.text(x+0.01, y_h+y_v+0.01, r'$t/c = 0.12$')
# plt.xlabel(r'Thickness/Chord Ratio, $t/c$')
# plt.ylabel(r'$C_{D_0}$')
# plt.show()