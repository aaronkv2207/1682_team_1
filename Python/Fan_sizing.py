import numpy as np
import matplotlib.pyplot as plt

# GOAL: Look at tradeoffs between radius for P_required during takeoff
# and eta_idela during cruise

# ============================================================
# ATMOSPHERE & AIRCRAFT CONSTANTS (SI UNITS)
# ============================================================

rho = 1.225              # Air density [kg/m^3]
g = 9.81                 # Gravity [m/s^2]

W = 73618.07             # Aircraft weight [N] --> 16550 lbs
m = W / g                # Aircraft mass [kg]

wing_loading = 14.8456   # Wing loading [N/m^2]
S = W / wing_loading     # Wing area [m^2]
AR = 8                   # Aspect ratio [-]
e = 0.7                  # Oswald efficiency factor [-]

CL = 6.1                 # Lift coefficient (TO config) [-]
CD0 = 0.32               # Parasite drag coefficient [-]

mu = 0.02                # Rolling friction coefficient [-]
x_runway = 52            # Runway length (171 ft) [m]

V_stall = 19.933         # Stall speed [m/s]
V_to = 1.2*V_stall       # Takeoff speed [m/s]
V_cruise = 125           # Cruise Speed [m/s] --> 280 mph

RPM = 2100               # RPM of propeller [rev/min]
omega = RPM*2*np.pi/60   # [rad/s]


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
# CRUISE MODEL (STEADY LEVEL FLIGHT)
# ============================================================

def CL_cruise(v):
    """Lift coefficient required for steady level flight [-]"""
    return W / (0.5 * rho * v**2 * S)


def CD_cruise(v):
    """Total drag coefficient at cruise [-]"""
    CLc = CL_cruise(v)
    CDi_c = CLc**2 / (np.pi * AR * e)
    return CD0 + CDi_c


def D_cruise(v):
    """Drag force at cruise [N]"""
    return 0.5 * rho * v**2 * S * CD_cruise(v)


def T_cruise(v):
    """Thrust required at cruise [N]"""
    return D_cruise(v)

def Tc_cruise(v, R):
    """Velocity-based thrust coefficient at cruise [-]"""
    A = Disk_area(R)
    return T_cruise(v) / (0.5 * rho * v**2 * A)

def Eta_ideal_cruise(v, R):
    Tc = Tc_cruise(v, R)
    return 2 / (1 + np.sqrt(1 + Tc))

def P_cruise(v, R):
    eta = Eta_ideal_cruise(v, R)
    return T_cruise(v) * v / eta


# ============================================================
# ACTUATOR DISK MODEL
# ============================================================

def Disk_area(R):
    """Disk area [m^2]"""
    return np.pi * R**2


def T_c(v, R):
    """Velocity-based thrust coefficient Tc [-]"""
    
    v = np.asarray(v)
    A = Disk_area(R)
    
    Tc = np.zeros_like(v)
    
    # Avoid divide-by-zero
    small_mask = v < 1e-6
    dynamic_mask = ~small_mask
    
    if np.any(dynamic_mask):
        Tc[dynamic_mask] = (
            T(v[dynamic_mask]) /
            (0.5 * rho * v[dynamic_mask]**2 * A)
        )
    
    # For very small velocity → treat as large loading
    Tc[small_mask] = 1e6
    
    return Tc


def Eta_ideal(v, R):
    """Ideal propulsive efficiency [-]"""
    Tc = T_c(v, R)
    return 2 / (1 + np.sqrt(1 + Tc))


# ============================================================
# POWER FUNCTION (FIXED FOR ARRAYS)
# ============================================================

def P_shaft_required(v, R):
    """Ideal shaft power required [W]"""
    
    v = np.asarray(v)
    P = np.zeros_like(v)

    # Static case (v very small)
    static_mask = v < 1e-6
    if np.any(static_mask):
        T_static = T(0)
        A = Disk_area(R)
        delta_V = np.sqrt(2 * T_static / (rho * A))
        P[static_mask] = T_static * delta_V / 2

    # Non-static case
    dynamic_mask = ~static_mask
    if np.any(dynamic_mask):
        eta = Eta_ideal(v[dynamic_mask], R)
        P[dynamic_mask] = T(v[dynamic_mask]) * v[dynamic_mask] / eta

    return P





vel = np.linspace(0, V_to, 200)
R_values = np.linspace(1.524, 3.048, 5) # Radius sweep (5 ft to 10 ft)


# ============================================================
# PLOT 1: THRUST COEFFICIENT VS VELOCITY
# ============================================================

vel_tc = np.linspace(1, V_to, 200)  # avoid zero

plt.figure(figsize=(8,6))

for r in R_values:
    plt.plot(vel_tc, T_c(vel_tc, r), label=f"R = {r:.2f} m")

plt.xlabel("Velocity [m/s]")
plt.ylabel("Velocity-Based Thrust Coefficient $T_c$ [-]")
plt.title("Thrust Coefficient vs Velocity (Ground Roll)")
plt.legend()
plt.grid(True)
plt.show()

# ============================================================
# PLOT 2: SHAFT POWER REQUIRED VS VELOCITY
# ============================================================

plt.figure(figsize=(8,6))

for r in R_values:
    plt.plot(vel, P_shaft_required(vel, r)/1000, label=f"R = {r:.2f} m")

plt.xlabel("Velocity [m/s]")
plt.ylabel("Shaft Power Required [kW]")
plt.title("Shaft Power Required vs Velocity (Ground Roll)")
plt.legend()
plt.grid(True)
plt.show()


# ============================================================
# PLOT 3: Cruise Efficiency vs Advance Ratio
# ============================================================
plt.figure(figsize=(8,6))
for r in R_values:
    J = V_cruise / (omega * r)
    eta = Eta_ideal_cruise(V_cruise, r)
    plt.scatter(J, eta, label=f"R = {r:.2f} m")

plt.xlabel("Advance Ratio J = V/(ΩR)")
plt.ylabel("Ideal Efficiency η")
plt.title("Cruise Efficiency vs Advance Ratio")
plt.legend()
plt.grid(True)
plt.show()



    