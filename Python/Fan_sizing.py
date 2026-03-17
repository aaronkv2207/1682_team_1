import numpy as np
import matplotlib.pyplot as plt

# GOAL: SIZE & DETERMINE THE AMOUNT OF FANS FOR OUR 19 PAX PLANE

# ============================================================
# ATMOSPHERE & AIRCRAFT CONSTANTS (SI UNITS)
# ============================================================

rho = 1.225   # Air density [kg/m^3]
rho_cruise = .549 # Air density at 25k ft [kg/m^3]
g = 9.81  # Gravity [m/s^2]
a = 343  # Speed of sound sea level [m/s^2]
a_cruise = 309 # Speed of sound cruise [m/s^2] --> at 25k ft

W = 73618.07  # Aircraft weight [N] --> 16550 lbs
m = W / g  # Aircraft mass [kg]

wing_loading = 1484.56  # Wing loading [N/m^2]
S = W / wing_loading  # Wing area [m^2]

AR = 8  # Aspect ratio [-]
e = 0.8  # Oswald efficiency factor [-]

CL = 6.1  # Lift coefficient (TO config) [-]
CD0 = 0.032  # Parasite drag coefficient [-]

mu = 0.02  # Rolling friction coefficient [-]
x_runway = 52  # Runway length (171 ft) [m]

V_stall = 19.933  # Stall speed [m/s]
V_to = 1.1*V_stall  # Takeoff speed [m/s]
V_cruise = 125  # Cruise Speed [m/s] --> 280 mph

RPM = 6000 # RPM of propeller [rev/min]
omega = RPM * 2 * np.pi / 60  # [rad/s]
N_fans = 8  # number of fans


# ============================================================
# DERIVED AERODYNAMIC COEFFICIENTS
# ============================================================

CDi = CL**2 / (np.pi * AR * e)  # Induced drag coefficient [-]
CD = CD0 + CDi  # Total drag coefficient [-]


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


def acc(v):
    """Required average acceleration to reach v in runway distance [m/s^2]."""
    return v**2 / (2 * x_runway)


def T(v):
    """Thrust required during ground roll [N]."""
    return D(v) + mu * (W - L(v)) + m * acc(v)


# At takeoff, (W-L(v)) is not included, as there is no friction force when we are off the ground
T_to = D(V_to) + acc(V_to)
print("Velocity at Takeoff:", V_to)
print("BASIC/NAIVE Thrust at takeoff:", T_to)
print("T/W at takeoff: ", T_to/W)

# ============================================================
# CRUISE MODEL (STEADY LEVEL FLIGHT)
# ============================================================


def CL_cruise(v):
    """Lift coefficient required for steady level flight [-]"""
    return W / (0.5 * rho_cruise * v**2 * S)


def CD_cruise(v):
    """Total drag coefficient at cruise [-]"""
    CLc = CL_cruise(v)
    CDi_c = CLc**2 / (np.pi * AR * e)
    return CD0 + CDi_c


def D_cruise(v):
    """Drag force at cruise [N]"""
    return 0.5 * rho_cruise * v**2 * S * CD_cruise(v)


def T_cruise(v):
    """Thrust required at cruise [N]"""
    return D_cruise(v)


def Tc_cruise(v, R):
    """Velocity-based thrust coefficient at cruise [-]"""
    A = Disk_area(R)
    T_per_fan = T_cruise(v) / N_fans
    return T_per_fan / (0.5 * rho * v**2 * A)


def Eta_ideal_cruise(v, R):
    Tc = Tc_cruise(v, R)
    return 2 / (1+np.sqrt(1 + Tc))


def P_cruise(v, R):
    eta_cruise = Eta_ideal_cruise(v, R)
    return T_cruise(v) * v / eta_cruise


def M_tip(R):
    """Tip Mach number at cruise"""
    return np.sqrt((omega * R) ** 2 + V_cruise**2) / a_cruise


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
        T_per_fan = T(v[dynamic_mask]) / N_fans
        Tc[dynamic_mask] = T_per_fan / (0.5 * rho * v[dynamic_mask] ** 2 * A)

    # For very small velocity → treat as large loading
    Tc[small_mask] = 1e6

    return Tc


def Eta_ideal(v, R):
    """Ideal propulsive efficiency [-]"""
    Tc = T_c(v, R)
    return 2 / (1 + np.sqrt(1 + Tc))


print("Effeciency during cruise=", Eta_ideal_cruise(125,3.5))
print("Cl=", CL_cruise(130))
print("Cd=", CD_cruise(130))
print("T_cruise=",T_cruise(130), "N")
print("P_cruise=", P_cruise(130, 1.8)/1000, "kW")


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
        T_static = T(0) / N_fans
        A = Disk_area(R)
        delta_V = np.sqrt(2 * T_static / (rho * A))
        P[static_mask] = T_static * delta_V / 2

    # Non-static case
    dynamic_mask = ~static_mask
    if np.any(dynamic_mask):
        eta = Eta_ideal(v[dynamic_mask], R)
        T_per_fan = T(v[dynamic_mask]) / N_fans
        P[dynamic_mask] = T_per_fan * v[dynamic_mask] / eta

    return P



# # ============================================================
# # PLOT: Tradeoff Plot
# # ============================================================

# R_plot = np.linspace(0.5, 5, 200)  # [m]
# vel = np.linspace(0.1, V_to, 200) # [m/s]

# P_vals = [] 
# for r in R_plot:
#     P_curve = N_fans * P_shaft_required(vel, r)
#     P_vals.append(np.max(P_curve) / 1000)
# eta_vals = [Eta_ideal_cruise(V_cruise, r) for r in R_plot]
# M_vals = [M_tip(r) for r in R_plot]

# fig, ax1 = plt.subplots(figsize=(9, 6))
# # Power curve
# ax1.plot(np.pi * R_plot**2, P_vals, label="Takeoff Power Required (kW)")
# ax1.set_xlabel("Effective Total Propeller Area [m^2]")
# ax1.set_ylabel("Power (kW)")

# motor_power = 2050  # kW 
# # ax1.axhline(motor_power, linestyle=":", label="Motor Power")

# # Propulsor Values from radius chosen using plots above
# R_selected = .583 # [m]
# A_selected = N_fans * np.pi * R_selected**2 # [m^2]
# lam = V_cruise / (omega*R_selected) 
# J = np.pi*lam
# torque_TO = 2.1e6/omega
# print("Ideal Effeciency at Takeoff:", Eta_ideal(V_to, R_selected))
# print("Ideal Effeciency at Cruise:", Eta_ideal_cruise(V_cruise, R_selected))
# print("Effective Total Fan Area:", A_selected)
# print("Advance Ratio: ", J)
# print("Individual Mach Tip number per fan for chosen radius:",M_tip(R_selected))
# print("Torque at takeoff:", torque_TO)

# # # Second axis for tip Mach
# # ax2 = ax1.twinx()
# # ax2.plot(R_plot, M_vals, color="red", label="Tip Mach")
# # ax2.set_ylabel("Tip Mach Number")
# # Mach limits
# # ax2.axhline(1.1, color="red", linestyle="--", label=f"Mach Limit={1.1}")

# # Combine legends
# lines1, labels1 = ax1.get_legend_handles_labels()
# # lines2, labels2 = ax2.get_legend_handles_labels()
# # ax1.legend(lines1 + lines2, labels1 + labels2)
# ax1.legend(lines1, labels1)

# plt.title("Power vs Propeller Effective Area")
# plt.show()





# ============================================================
# THRUST INTERPOLATOR (first used in takeoff model) (0–22 m/s)
# ============================================================

V_data = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0])
T_data = 8 * np.array([3604.0, 3280.0, 2999.0, 2759.0, 2557.0, 2381.0, 2224.0, 2088.0, 1969.0, 1866.0, 1776.0]) # Thrust per fan [N]
degree = 5          # change this to 1,2,3,4,... to test fits
coeffs = np.polyfit(V_data, T_data, degree)
T_poly = np.poly1d(coeffs)
print("Current Thrust at Takeoff:", T_poly(V_to)/1000, " [kN]")

class ThrustInterpolator:
    def __init__(self, V_data=V_data, T_data=T_data, degree=degree):
        self.coeffs = np.polyfit(V_data, T_data, degree)
        self.poly = np.poly1d(self.coeffs)

    def get_T(self, v):
        """Interpolated thrust per fan [N]"""
        return self.poly(v)

def T_fan_interp(v):
        """Interpolated thrust per fan [N]"""
        return T_poly(v)

vel = np.linspace(0,130,500)
# plt.figure(figsize=(8,6))
# plt.scatter(V_data, T_data/1000, label="Original Data")
# plt.plot(vel, T_fan_interp(vel)/1000, label=f"Polynomial Fit (deg={degree})")
# plt.xlabel("Velocity [m/s]")
# plt.ylabel("Thrust per fan [kN]")
# plt.title("Fan Thrust Interpolation")
# plt.legend()
# plt.show()



# Size ducted fan A/2, R/sqrt(2)
# cl_avg about 0.5


# ============================================================
# Fan Blade Distribution Plot
# ============================================================
n_blades = 5 # number of blades
R_selected = .583 # [m]
r = R_selected* np.array([
    0.05944, 0.07832, 0.0972, 0.11608, 0.13496, 0.15384, 0.17272, 0.1916,
    0.21048, 0.22936, 0.24824, 0.26712, 0.286, 0.30488, 0.32376, 0.34264,
    0.36152, 0.3804, 0.39928, 0.41816, 0.43704, 0.45592, 0.4748, 0.49368,
    0.51256, 0.522
])
c = R_selected* np.array([
    0.1856, 0.111, 0.17096, 0.22447, 0.26976, 0.30662, 0.33561, 0.35768,
    0.36382, 0.37497, 0.38193, 0.38534, 0.3857, 0.38336, 0.37849, 0.37115,
    0.36122, 0.34838, 0.33212, 0.3116, 0.29555, 0.28201, 0.26759, 0.25536,
    0.24299, 0.22757
])
beta = np.deg2rad(np.array([
    85.0212, 80.8544, 77.7618, 74.7565, 71.849, 69.047, 66.3561, 63.7793,
    61.3179, 58.9711, 56.7372, 54.6132, 52.5953, 50.6792, 48.8601, 47.1332,
    45.4936, 43.9361, 42.456, 41.0487, 39.7095, 38.4342, 37.2188, 36.0593,
    34.9523, 32.4185
]))

y_half = 1/2*c*np.cos(beta)
y_full = np.concatenate([y_half, -y_half[::-1]])
r_full = np.concatenate([r, r[::-1]])

# PLOTTING 1 BLADE
plt.figure(figsize=(8, 3))
plt.plot(r_full, y_full, '-o', color='blue')
plt.fill(r_full, y_full, alpha=0.3, color='skyblue')
plt.xlabel('Radius (m)')
plt.ylabel('y (m)')
plt.title('Fan Blade Distribution')
plt.grid(True)
plt.axis('equal')
plt.show()

# PLOTTING ALL BLADES
plt.figure(figsize=(6,6))
for i in range(n_blades):
    angle = i * 2*np.pi / n_blades
    # rotate blade coordinates
    x_rot = r_full * np.cos(angle) - y_full * np.sin(angle)
    y_rot = r_full * np.sin(angle) + y_full * np.cos(angle)
    plt.fill(x_rot, y_rot, alpha=0.4, label=f'Blade {i+1}')

plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.title('Full Fan with All Blades')
plt.axis('equal')
plt.grid(True)
plt.show()




















# # ============================================================
# # GRAVEYARD PLOTS/FUNCTIONS
# # ============================================================



# # ============================================================
# # Sanity Check Plot: CRUISE EFFICIENCY VS RADIUS
# # ============================================================

# R_sweep = np.linspace(1.0, 10.0, 200)   # sweep radius from 1m to 4m

# eta_cruise = Eta_ideal_cruise(V_cruise, R_sweep)

# plt.figure(figsize=(8,6))
# plt.plot(R_sweep, eta_cruise)

# plt.xlabel("Propeller Radius [m]")
# plt.ylabel("Ideal Propulsive Efficiency η")
# plt.title("Ideal Cruise Efficiency vs Propeller Radius")

# plt.grid(True)
# plt.show()


# vel = np.linspace(0, V_to, 200)
# R_values = np.linspace(1.524, 3.048, 5) # Radius sweep (5 ft to 10 ft)


# # ============================================================
# # PLOT 1: THRUST COEFFICIENT VS VELOCITY
# # ============================================================

# vel_tc = np.linspace(1, V_to, 200)  # avoid zero

# plt.figure(figsize=(8,6))

# for r in R_values:
#     plt.plot(vel_tc, T_c(vel_tc, r), label=f"R = {r:.2f} m")

# plt.xlabel("Velocity [m/s]")
# plt.ylabel("Velocity-Based Thrust Coefficient $T_c$ [-]")
# plt.title("Thrust Coefficient vs Velocity (Ground Roll)")
# plt.legend()
# plt.grid(True)
# plt.show()

# # ============================================================
# # PLOT 2: SHAFT POWER REQUIRED VS VELOCITY
# # ============================================================

# plt.figure(figsize=(8,6))

# for r in R_values:
#     plt.plot(vel, P_shaft_required(vel, r)/1000, label=f"R = {r:.2f} m")

# plt.xlabel("Velocity [m/s]")
# plt.ylabel("Shaft Power Required [kW]")
# plt.title("Shaft Power Required vs Velocity (Ground Roll)")
# plt.legend()
# plt.grid(True)
# plt.show()


# # ============================================================
# # PLOT 3: Cruise Efficiency vs Advance Ratio
# # ============================================================
# plt.figure(figsize=(8,6))
# for r in R_values:
#     J = V_cruise / (omega * r)
#     eta = Eta_ideal_cruise(V_cruise, r)
#     plt.scatter(J, eta, label=f"R = {r:.2f} m")

# plt.xlabel("Advance Ratio J = V/(ΩR)")
# plt.ylabel("Ideal Efficiency η")
# plt.title("Cruise Efficiency vs Advance Ratio")
# plt.legend()
# plt.grid(True)
# plt.show()


# # ============================================================
# # PLOT 4: Tip Mach # During Cruise
# # ============================================================

# plt.figure(figsize=(8,6))
# R_plot = np.linspace(1.0, 4.0, 200)
# plt.plot(R_plot, M_tip(R_plot), label="Tip Mach")
# plt.axhline(0.65, linestyle="--", label="Efficiency limit")
# plt.axhline(0.85, linestyle="--", label="Compressibility limit")
# plt.xlabel("Propeller Radius [m]")
# plt.ylabel("Tip Mach Number")
# plt.title("Propeller Tip Mach at Cruise")
# plt.legend()
# plt.grid(True)
# plt.show()