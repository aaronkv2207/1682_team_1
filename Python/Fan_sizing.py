import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# ATMOSPHERE & AIRCRAFT CONSTANTS (SI UNITS)
# ============================================================

rho = 1.225              # Air density [kg/m^3]
g = 9.81                 # Gravity [m/s^2]

W = 84516.21             # Aircraft weight [N]
m = W / g                # Aircraft mass [kg]

S = 20.0                 # Wing area [m^2]
AR = 20.0                # Aspect ratio [-]
e = 0.7                  # Oswald efficiency factor [-]

CL = 0.2                 # Lift coefficient (TO config) [-]
CD0 = 0.1                # Parasite drag coefficient [-]

mu = 0.02                # Rolling friction coefficient [-]
x_runway = 91.44         # Runway length (300 ft) [m]

V_stall = 20.0           # Stall speed [m/s]


# ============================================================
# DERIVED AERODYNAMIC COEFFICIENTS
# ============================================================

CDi = CL**2 / (np.pi * AR * e)     # Induced drag coefficient [-]
CD = CD0 + CDi                    # Total drag coefficient [-]


# ============================================================
# AERODYNAMIC FORCES
# ============================================================

def L(v):
    """Lift force [N]"""
    return 0.5 * rho * v**2 * S * CL


def D(v):
    """Drag force [N]"""
    return 0.5 * rho * v**2 * S * CD


# ============================================================
# TAKEOFF GROUND ROLL MODEL
# ============================================================

def a(v):
    """Required average acceleration to reach v in runway distance [m/s^2]."""
    return v**2 / (2 * x_runway)


def T(v):
    """Thrust required during ground roll [N]."""
    return (
        D(v)
        + mu * (W - L(v))
        + m * a(v)
    )


# ============================================================
# ACTUATOR DISK MODEL
# ============================================================

def Disk_area(R):
    """Disk area [m^2]"""
    return np.pi * R**2


def T_c(v, R):
    """Velocity-based thrust coefficient Tc [-]"""
    if v == 0:
        return np.inf
    return T(v) / (0.5 * rho * v**2 * Disk_area(R))


def Eta_ideal(v, R):
    """Ideal propulsive efficiency [-]"""
    Tc = T_c(v, R)
    return 2 / (1 + np.sqrt(1 + Tc))


def P_shaft_required(v, R):
    """Ideal shaft power required [W]"""
    if v == 0:
        # Static case
        T = T(0)
        A = Disk_area(R)
        delta_V = np.sqrt(2 * T / (rho * A))
        return T * delta_V / 2
    else:
        return T(v) * v / Eta_ideal(v, R)






    