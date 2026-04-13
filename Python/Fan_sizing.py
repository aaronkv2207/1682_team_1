import numpy as np
import matplotlib.pyplot as plt
# from aero_workspace.aero_dict import AircraftConfig

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
V_to = 1.1*V_stall # Takeoff speed [m/s]
V_cruise = 125  # Cruise Speed [m/s] --> 280 mph

RPM = 7500 # RPM of propeller [rev/min]
omega = RPM * 2 * np.pi / 60  # [rad/s]
N_fans = 8  # number of fans

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

print("Thrust at takeoff with friction force:", T(V_to))
print("Thrust at takeoff no friction force:", D(V_to)+m*acc(V_to))
print("Friction force:", mu * (W - L(V_to)))
print("Drag:", D(V_to))
print("ma force:", m*(acc(V_to)))
print("Acceleration:", acc(V_to))


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
    A_total = N_fans*Disk_area(R)  # total disk area for all fans

    # Avoid divide-by-zero
    if v < 1e-6:
        return 1e6
    else:
        return T(v) / (0.5 * rho * v ** 2 * A_total)


def Eta_ideal(v, R):
    """Ideal propulsive efficiency [-]"""
    Tc = T_c(v, R)
    return 2 / (1 + np.sqrt(1 + Tc))


# ============================================================
# POWER FUNCTION (TOTAL, SCALAR VERSION)
# ============================================================

def P_shaft_required(v, R):
    """Ideal shaft power required [W]"""
    A = Disk_area(R)

    if v < 1e-6:  # Static case
        T_static = T(0)
        delta_V = np.sqrt(2 * T_static / (rho * A))
        return T_static * delta_V / 2
    else:  # Non-static case
        ##### TRY ASSUMING THE EFFECINECY QFAN GIVES AND SEE IF THAT WORKS WITH 27.5E3
        eta = Eta_ideal(v, R)
        T_total = 62e3
        return T_total * v / eta

# print("P_required: ", P_shaft_required(V_to, 0.583))

# ============================================================
# PLOT: Tradeoff Plot
# ============================================================

R_plot = np.linspace(0.5, 5, 200)   # radius range [m]                       
P_vals = []
A_vals = []

for R in R_plot:
    P_total = P_shaft_required(V_to, R)  # total power
    A_total = np.pi * R**2                         # total area

    P_vals.append(P_total / 1000)  # convert to kW
    A_vals.append(A_total)
P_vals = np.array(P_vals)
A_vals = np.array(A_vals)

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(A_vals, P_vals, label="Takeoff Power Required (kW)")

motor_power = 2100  # kW
A_chosen = 8.542
# ax.axhline(motor_power, linestyle="--", label="Motor Limit (2100 kW)")
ax.axvline(A_chosen, linestyle="--", label="Area chosen")

ax.set_xlabel("Effective Total Propeller Area [m²]")
ax.set_ylabel("Power (kW)")
ax.set_title("Takeoff Power vs Effective Area")
ax.legend()
plt.show()

# Propulsor Values from radius chosen using plots above
R_selected = .583 # [m]
A_selected = N_fans * np.pi * R_selected**2 # [m^2]
lam = V_cruise / (omega*R_selected) 
J = np.pi*lam




# def runway_dist(TW, CL_val, v_eval):
#     """
#     Compute takeoff distance assuming constant average acceleration.
#     """

#     # Induced drag
#     CDi = CL_val**2 / (np.pi * AR * e)
#     CD = CD0 + CDi

#     # Local lift and drag functions
#     def L_local(v):
#         return 0.5 * rho * v**2 * S * CL_val

#     def D_local(v):
#         return 0.5 * rho * v**2 * S * CD

#     # Forces at representative velocity
#     L_eval = L_local(v_eval)
#     D_eval = D_local(v_eval)

#     # Thrust from T/W
#     T_const = TW * W

#     # Acceleration including rolling friction
#     a = (T_const - D_eval - mu*(W - L_eval)) / m

#     if a <= 0:
#         return np.nan  # cannot take off

#     # Runway distance using kinematics
#     x = V_to**2 / (2 * a)

#     return x

# TW_range = np.linspace(0.2, 1.5, 60)
# CL_range = np.linspace(1.5, 10.5, 60)

# X, Y = np.meshgrid(TW_range, CL_range)
# Z = np.zeros_like(X)

# # Representative velocity ~ 60% of takeoff speed
# v_eval = 0.6 * V_to

# for i in range(len(CL_range)):
#     for j in range(len(TW_range)):
#         Z[i, j] = runway_dist(TW_range[j], CL_range[i], v_eval)


# Z_plot = np.ma.masked_where(Z > 200, Z) # Mask extreme values for plotting

# # Plot
# plt.figure(figsize=(10, 7))
# cp = plt.contourf(X, Y, Z_plot, levels=np.linspace(0, 200, 11))
# plt.colorbar(cp, label="Runway Distance [m]")

# plt.xlabel("Thrust-to-Weight Ratio (T/W)")
# plt.ylabel("Lift Coefficient (CL)")
# plt.title("Takeoff Distance Trade Study (Constant Acceleration)")
# plt.grid(alpha=0.3)
# plt.show()


# ============================================================
# Printing Section for important values
# ============================================================

# TAKEOFF
# At takeoff, (W-L(v)) is not included, as there is no friction force when we are off the ground
# T_to = D(V_to) + m*acc(V_to)
T_to = 40e3
print("T/W at takeoff: ", T_to/W)
if False:
    print("--------------------------------")
    print("TAKOFF VALUES")
    print("Velocity at Takeoff:", V_to)
    print("BASIC/NAIVE Thrust at takeoff:", T_to)
    print("T/W at takeoff: ", T_to/W)
    print("Ideal Effeciency at Takeoff:", Eta_ideal(V_to, R_selected))
    print("Tc taekoff = ", T_c(V_to, R_selected))
    print("-------------------------------")

    # CRUISE
    print("CRUISE VALUES")
    print("Effeciency during cruise =", Eta_ideal_cruise(125,3.5))
    print("Cl_cruise =", CL_cruise(V_cruise))
    print("Cd_cruise =", CD_cruise(V_cruise))
    print("T_cruise =",T_cruise(V_cruise), "N")
    print("Tc cruise:, ", Tc_cruise(V_to, R_selected))
    print("P_cruise =", P_cruise(V_cruise, 1.8)/1000, "kW")
    print("Ideal Effeciency at Cruise:", Eta_ideal_cruise(V_cruise, R_selected))
    print("Cruise Advance Ratio: ", J)
    print("-------------------------------")


    # FAN GEOMETRY
    print("FAN GEOMETRY VALUES")
    print("Effective Total Fan Area:", A_selected)
    print("Individual Mach Tip number per fan for chosen radius:",M_tip(R_selected)) # USE QFAN FOR M_TIP
    print("-------------------------------")




# # ============================================================
# # Fan Blade Distribution Plot
# # ============================================================
# n_blades = 5 # number of blades
# R_selected = .583 # [m]
# r = R_selected* np.array([
#     0.05944, 0.07832, 0.0972, 0.11608, 0.13496, 0.15384, 0.17272, 0.1916,
#     0.21048, 0.22936, 0.24824, 0.26712, 0.286, 0.30488, 0.32376, 0.34264,
#     0.36152, 0.3804, 0.39928, 0.41816, 0.43704, 0.45592, 0.4748, 0.49368,
#     0.51256, 0.522
# ])
# c = R_selected* np.array([
#     0.1856, 0.111, 0.17096, 0.22447, 0.26976, 0.30662, 0.33561, 0.35768,
#     0.36382, 0.37497, 0.38193, 0.38534, 0.3857, 0.38336, 0.37849, 0.37115,
#     0.36122, 0.34838, 0.33212, 0.3116, 0.29555, 0.28201, 0.26759, 0.25536,
#     0.24299, 0.22757
# ])
# beta = np.deg2rad(np.array([
#     85.0212, 80.8544, 77.7618, 74.7565, 71.849, 69.047, 66.3561, 63.7793,
#     61.3179, 58.9711, 56.7372, 54.6132, 52.5953, 50.6792, 48.8601, 47.1332,
#     45.4936, 43.9361, 42.456, 41.0487, 39.7095, 38.4342, 37.2188, 36.0593,
#     34.9523, 32.4185
# ]))

# y_half = 1/2*c*np.cos(beta)
# y_full = np.concatenate([y_half, -y_half[::-1]])
# r_full = np.concatenate([r, r[::-1]])

# # PLOTTING 1 BLADE
# plt.figure(figsize=(8, 3))
# plt.plot(r_full, y_full, '-', color='blue')
# plt.fill(r_full, y_full, alpha=0.3, color='skyblue')
# plt.xlabel('Radius (m)')
# plt.ylabel('y (m)')
# plt.title('Fan Blade Distribution')
# # plt.grid(True)
# plt.axis('equal')
# plt.show()

# # PLOTTING ALL BLADES
# plt.figure(figsize=(6,6))
# for i in range(n_blades):
#     angle = i * 2*np.pi / n_blades
#     # rotate blade coordinates
#     x_rot = r_full * np.cos(angle) - y_full * np.sin(angle)
#     y_rot = r_full * np.sin(angle) + y_full * np.cos(angle)
#     plt.fill(x_rot, y_rot, alpha=0.4, color="skyblue", label=f'Blade {i+1}')

# plt.xlabel('X (m)')
# plt.ylabel('Y (m)')
# plt.title('Full Fan with All Blades')
# plt.axis('equal')
# # plt.grid(True)
# plt.show()


# ============================================================
# Mass calculation for ducted fan system
# ============================================================

# def fan_mass_composites(R, chord_distribution, t, N_blades=5, materials=None):
#     """
#     Estimate the mass of a fan for different composite materials.

#     Parameters:
#     R : float
#         Fan radius [m]
#     chord_distribution : array-like
#         Chord at each radial station (m)
#     N_blades : int
#         Number of blades
#     t : float
#         Thickness of the blade (m)
#     materials : dict
#         Dictionary of materials with density [kg/m^3].
#         Example: {'Carbon Fiber': 1600, 'Aluminum': 2700}
    
#     Returns:
#     dict
#         Estimated mass of a single fan for each material [kg]
#     """
#     if materials is None:
#         materials = {'Carbon Fiber': 1600}  # default
    
#     chord_distribution = np.array(chord_distribution)
#     c_avg = np.mean(chord_distribution)
#     volume_blade = c_avg * t * R  # m^3
#     mass_dict = {}
    
#     for mat, rho in materials.items():
#         mass_dict[mat] = N_blades * volume_blade * rho
    
#     return mass_dict

# materials = {
#     'Carbon Fiber': 1600,
#     'Fiberglass': 2000,
#     'Aluminum': 2700,
#     'Titanium': 4500
# }

# print("MASS VALUES FOR EACH MATERIAL")
# t = 0.0121 # blade thickness [cm]
# fan_masses = fan_mass_composites(R_selected, c, t, N_blades=5)
# for mat, mass in fan_masses.items():
#     print(f"{mat}: {8*mass:.2f} kg")






# # ============================================================
# # ============================================================
# # ============================================================
# # GRAVEYARD PLOTS/FUNCTIONS
# # ============================================================
# # ============================================================
# # ============================================================



# # ============================================================
# # THRUST INTERPOLATOR (first used in takeoff model) (0–22 m/s)
# # ============================================================

# V_data = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0])
# T_data = 8 * np.array([3604.0, 3280.0, 2999.0, 2759.0, 2557.0, 2381.0, 2224.0, 2088.0, 1969.0, 1866.0, 1776.0]) # Thrust per fan [N]
# degree = 5          # change this to 1,2,3,4,... to test fits
# coeffs = np.polyfit(V_data, T_data, degree)
# T_poly = np.poly1d(coeffs)
# print("Current Thrust at Takeoff:", T_poly(V_to)/1000, " [kN]")

# class ThrustInterpolator:
#     def __init__(self, V_data=V_data, T_data=T_data, degree=degree):
#         self.coeffs = np.polyfit(V_data, T_data, degree)
#         self.poly = np.poly1d(self.coeffs)

#     def get_T(self, v):
#         """Interpolated thrust per fan [N]"""
#         return self.poly(v)

# def T_fan_interp(v):
#         """Interpolated thrust per fan [N]"""
#         return T_poly(v)

# vel = np.linspace(0,130,500)
# plt.figure(figsize=(8,6))
# plt.scatter(V_data, T_data/1000, label="Original Data")
# plt.plot(vel, T_fan_interp(vel)/1000, label=f"Polynomial Fit (deg={degree})")
# plt.xlabel("Velocity [m/s]")
# plt.ylabel("Thrust per fan [kN]")
# plt.title("Fan Thrust Interpolation")
# plt.legend()
# plt.show()



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