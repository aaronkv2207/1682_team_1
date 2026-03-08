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
AR = 8  # from twin otter
h_cruise = 10000 * ureg("ft").to("m").magnitude
rho_cruise = Atmosphere(h=h_cruise).density[0]  # [kg/m^2]

b = 19.9  # [m]
c = 1.98  # [m]
S = b**2 / AR
R = 2.39 / 2  # estimated
x_to = 366  # [m]

m = 5670  # [kg]
W = m * 9.81  # [N]


# fuselage (roskam pt 2)
D_f = 1.75  # [m] taken from the twin otter
l_f = 15.77  # [m] taken from the twin otter
lambda_f = l_f / D_f

C_f = 0.002  # skin friction coefficient, estimated from the graph on pg 4 of airmodels.pdf


# nacelles (roskam pt 2)
l_n = 2  # [m] length of cowl (est)
D_n = 1.8  # [m] diameter of cowl (est)
l_1 = l_n / 2
D_hl = 1.295  # [m] diameter of fan (1/2 of twin otter prop)
D_ef = D_hl


# tail
S_h = 12.54  # [m^2] from twin otter
S_v = 7.61  # [m^2] from twin otter

l_h = 9.06  # [m] estimated from twin otter
l_v = 9.06  # [m] estimated from twin otter

V_h = S_h * l_h / (S * c)
V_v = S_v * l_v / (S * b)


# landing gear


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
S_wet_wing = S * 2
S_wet_htail = S_h * 2
S_wet_vtail = S_v * 2

S_wet_total = S_wet_fuse + S_wet_fancowl + S_wet_wing + S_wet_htail + S_wet_vtail

f = 3  # estimated equivalent parasite area (Roskam Airplane Design I, pg 119)


# cruise
mu_cruise = Atmosphere(h=h_cruise).dynamic_viscosity[0]
V_cruise = 93.63  # [m/s]
C_L_cruise = W / (1 / 2 * rho_cruise * V_cruise**2 * S)


# takeoff
rho_sl = Atmosphere(h=0).density[0]
mu_sl = Atmosphere(h=0).dynamic_viscosity[0]  # [Ns/m^2] dynamic viscosity at sea level
V_takeoff = 1.2 * V_STALL  # [m/s]
C_L_takeoff = cl_max  # from wt data???


print("AR =", AR)
print("wing area =", round(S, 2), "m^2")
print("total wetted area =", round(S_wet_total, 2), "m^2")
print("takeoff distance =", x_to, "m")
print("MTOW =", round(W, 2), "N")

print("h tail coefficient =", round(V_h, 3))
print("v tail coefficient =", round(V_v, 3))


# calculate zero lift drag
def C_D0():  # roskam's method
    C_D0 = f / S_wet_total  # zero-lift drag coefficient

    print("C_D0 (roskam) =", round(C_D0, 3))

    return C_D0


C_D_cruise = C_D0()


# calculate profile drag (dissipation summation buildup)
def C_Dp(rho, V, mu, l, x_tr):
    Re_l = rho * V * l / mu
    Re_x_tr = rho * V * x_tr / mu

    C_fl = 1.328 / Re_l ** (1 / 2)
    C_ft = 0.455 / (np.log10(Re_l) ** 2.58)

    Cf = max(C_fl, C_ft - (Re_x_tr / 320 - 39) / Re_l)

    CDA = S_wet * C_f * Kf

    return C_Dp


# = sum(CDA * 1/2 * Vi**3)

# calculate induced drag
def C_Di(C_L):
    C_Di = C_L**2 / (np.pi*AR*e)
    return C_Di
