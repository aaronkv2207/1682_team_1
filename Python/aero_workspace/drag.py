#!/usr/bin/env python3

from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import numpy as np
from ambiance import Atmosphere


@dataclass
class AircraftConstants:
    """Defaults to v2 Plane 1. Update constants via class call if changed."""

    b: float = 17.75
    MAC: float = 2.54

    D_f: float = 1.6
    l_f: float = 15

    D_n: float = 1.1 * 1.116
    D_hl: float = 1.116

    S_h: float = 13.33
    ht_MAC: float = 2.11
    S_v: float = 5.99
    vt_MAC: float = 2.23
    l_h: float = 8.56
    l_v: float = 8.0

    tire_width: float = 0.284
    tire_radius: float = 0.818 / 2

    c_s: float = 0.5
    l_s: float = 1.0

    # initialized values
    l_n: float = 0.0
    S: float = 0.0
    AR: float = 0.0
    lambda_f: float = 0.0
    S_s: float = 0.0
    V_h: float = 0.0
    V_v: float = 0.0
    S_wet_fuse: float = 0.0
    S_wet_fancowl: float = 0.0
    S_wet_wing: float = 0.0
    S_wet_htail: float = 0.0
    S_wet_vtail: float = 0.0
    S_wet_gear: float = 0.0
    S_wet_strut: float = 0.0

    def __init__(self, **kwargs):
        # just deriving inital params based on inputs
        for k, v in kwargs.items():
            setattr(self, k, v)
        self.l_n = self.MAC / 2
        self.S = self.b * self.MAC
        self.AR = self.b**2 / self.S
        self.lambda_f = self.l_f / self.D_f
        self.S_s = self.c_s * self.l_s

        self.V_h = self.S_h * self.l_h / (self.S * self.MAC)
        self.V_v = self.S_v * self.l_v / (self.S * self.b)

        l_1 = self.l_n / 2

        self.S_wet_fuse = (
            np.pi
            * self.D_f
            * self.l_f
            * (1 - 2 / self.lambda_f) ** (2 / 3)
            * (1 + 1 / self.lambda_f**2)
        )

        self.S_wet_fancowl = (
            self.l_n
            * self.D_n
            * (
                2
                + 0.35 * l_1 * self.l_n
                + 0.8 * l_1 * self.D_hl / self.l_n * self.D_n
                + 1.15 * (1 - l_1 / self.l_n) * self.D_hl / self.D_n
            )
        )

        self.S_wet_wing = 2 * self.S
        self.S_wet_htail = 2 * self.S_h
        self.S_wet_vtail = 2 * self.S_v
        self.S_wet_gear = self.tire_width * (np.pi * self.tire_radius * 2)
        self.S_wet_strut = 2 * self.S_s


C = AircraftConstants()  # instantiate constants


# some helper functions
def calc_Re_l(rho, V, mu, l):
    return rho * V * l / mu


def calc_C_f(Re_l):
    return 0.455 / (np.log10(Re_l) ** 2.58)


def calc_CDA(S_wet, C_f, K_f):
    return S_wet * C_f * K_f


def calc_K_f_airfoil(tc):
    return 1 + 2.0 * tc + 60 * tc**4


def calc_K_f_axi(dl):
    return 1 + 1.5 * dl ** (3 / 2) + 7 * dl**3


# Drag buildup
def calc_C_Dp(rho, V, mu, C: AircraftConstants = C):
    S_wets = [
        C.S_wet_fuse,
        C.S_wet_fancowl,
        C.S_wet_wing,
        C.S_wet_htail,
        C.S_wet_vtail,
        C.S_wet_gear,
        C.S_wet_strut,
    ]

    ls = [
        C.l_f,
        C.l_n,
        C.MAC,
        C.ht_MAC,
        C.vt_MAC,
        C.tire_radius * 2,
        C.l_s,
    ]

    fineness_ratios = [
        C.D_f / C.l_f,
        C.D_n / (2 * C.l_n),
        0.12,
        0.10,
        0.10,
        0.10,
        0.25,
    ]

    CDAs = []

    for i, S_wet in enumerate(S_wets):
        Re_l = calc_Re_l(rho, V, mu, ls[i])
        C_f = calc_C_f(Re_l)

        if i <= 1:
            K_f = calc_K_f_axi(fineness_ratios[i])
        else:
            K_f = calc_K_f_airfoil(fineness_ratios[i])

        CDA = calc_CDA(S_wet, C_f, K_f)

        if i == 1:  # 8 nacelles
            CDA *= 8

        CDAs.append(CDA)

    Dps = [CDA * 0.5 * rho * V**2 for CDA in CDAs]
    C_Dps = [Dp / (0.5 * rho * V**2 * C.S) for Dp in Dps]

    return sum(Dps), sum(C_Dps)


# Just for testing
h_cruise = 5486.4
rho_cruise = Atmosphere(h=h_cruise).density[0]
mu_cruise = Atmosphere(h=h_cruise).dynamic_viscosity[0]
V_cruise = 125


Dp_cruise, C_Dp_cruise = calc_C_Dp(rho_cruise, V_cruise, mu_cruise, C)
print("-------")
print("cruise")
print("Cd_p =", round(C_Dp_cruise, 3))
