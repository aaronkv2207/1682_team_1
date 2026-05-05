import numpy as np
import matplotlib.pyplot as plt
# from aero_workspace.aero_dict import AircraftConfig

# print(AircraftConfig.v_t0)

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
# print(m)

# wing_loading = 1484.56  # Wing loading [N/m^2]
# S = W / wing_loading  # Wing area [m^2]
S = 45
c = 2.54 # Wing chord [m] --> FAN BLADE RADIUS SHOULD BE 25% OR LESS OF THIS NUMBER
b = 17.75 # Wingspan [m] --> Fans should fit probs on 14 m

AR = 7  # Aspect ratio [-] was 8
e = 0.8  # Oswald efficiency factor [-]

CL = 6.22  # Lift coefficient (TO config) [-]
CD0 = 0.032  # Parasite drag coefficient [-]

mu = 0.02  # Rolling friction coefficient [-]
x_runway = 57.87  # Runway length (171 ft) [m]

# V_stall = 19.933  # Stall speed [m/s]
# V_to = 1.1*V_stall # Takeoff speed [m/s]
V_to = 19.47 #NOTE: IN QMIL/QFAN USE HALF OF V_to BC YOU WANT TO CONSIDER AN AVERGAE OF STATIC AND TAKEOFF CONDITIONS
V_cruise = 125  # Cruise Speed [m/s] --> 280 mph

RPM = 6500 # RPM of propeller [rev/min]
omega = RPM * 2 * np.pi / 60  # [rad/s]


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
# TAKEOFF GROUND ROLL MODEL (W/ FOCUS ON TAKEOFF)
# ============================================================
# At takeoff, (W-L(v)) is not included, as there is no friction force when we are off the ground
# T_to = D(V_to) + m*acc(V_to)

def acc(v):
    """Required average acceleration to reach v in runway distance [m/s^2]."""
    return v**2 / (2 * x_runway)


def T(v):
    """Thrust required during ground roll [N]."""
    return D(v) + mu * (W - L(v)) + m * acc(v)

D_to = 22.05e3
# print("Thrust at takeoff with friction force:", T(V_to))
# print("Thrust at takeoff no friction force:", D(V_to)+m*acc(V_to))
# print("Friction force:", mu * (W - L(V_to)))
# print("Drag:", D(V_to))
print("ma force:", m*(acc(V_to)))
# print("Acceleration:", acc(V_to))

print("Thrust at takeoff: " , D_to + m*(acc(V_to)))
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
    return T_cruise(v) / (0.5 * rho * v**2 * A)


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
    A_total = Disk_area(R)  # total disk area for all fans

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
    A = Disk_area(R) # TOTAL EFFECTIVE AREA

    if v < 1e-6:  # Static case
        T_static = T(0)
        delta_V = np.sqrt(2 * T_static / (rho * A))
        return T_static * delta_V / 2
    else:  # Non-static case
        ##### TRY ASSUMING THE EFFECINECY QFAN GIVES AND SEE IF THAT WORKS WITH 27.5E3
        eta = Eta_ideal(v, R)
        # T_total = T(v)
        T_total = 41e3
        return T_total * v / eta

# def T_required(P_shaft, v, R):
#     return P_shaft / 


print("T_produced: ", 0.5*3.5e6/20 * 0.25 / 1000, "kN")

# ============================================================
# AREA PER FAN FUNCTION
# ============================================================

def area_per_fan(A_eff, N):
    """For a given total effective area and number of fans, 
    this returns the area AND radius per fan"""

    A_per_fan = A_eff / N
    r_per_fan = np.sqrt(A_per_fan / (np.pi))
    return (A_per_fan, r_per_fan)

# # Test case
# A,r = area_per_fan(8.54, 8)
# print("Area:", A) # Should print 1.0675
# print("Radius:", r) # Should print 0.583


# ============================================================
# PLOT: Tradeoff Plot
# ============================================================

R_plot = np.linspace(0.5, 5, 200)   # TOTAL EFFECTIVE radius range [m]                       
P_vals = []
A_vals = []

for R in R_plot:
    P_total = P_shaft_required(V_to, R)  # total power
    A_total = np.pi * R**2               # total area

    P_vals.append(P_total / 1000)  # convert to kW
    A_vals.append(A_total)
P_vals = np.array(P_vals)
A_vals = np.array(A_vals)

fig, ax = plt.subplots(figsize=(9, 6))
ax.plot(A_vals, P_vals, label="Takeoff Power Required (kW)")

motor_power = 3.5e3  # kW
A_chosen = 3.65
# ax.axhline(motor_power, color = "k", linestyle="--", label=f"Chosen Motor Power: {motor_power} kW")
# ax.axvline(A_chosen, linestyle="--", label=f"Area chosen: {A_chosen} m^2")

ax.set_xlabel("Effective Total Propeller Area [m²]", fontsize=15)
ax.set_ylabel("Power (kW)", fontsize=15)
ax.set_title("Takeoff Power vs Effective Area", fontsize=20)
ax.legend()
plt.show()




# ============================================================
# Printing Section for important values
# ============================================================

N_fans = 14
R_selected = .282 # [m]
A_selected = N_fans * np.pi * R_selected**2 # [m^2]
lam = V_cruise / (omega*R_selected) 
J = np.pi*lam

# Test case: NOTE: Max radius to meet 25% diameter requirement = 0.3175 [m]
#            NOTE: b = 17.75 m --> -2 m due to fuselage, -1 m due to spacing for wingtips(0.5 each side) 
#                              --> have 14.75 m to work with 
A_eff = 3.49
A,r = area_per_fan(A_eff, N_fans)
print(f"{N_fans = }")
print(f"Area per fan: {A:.3f} [m^2]")
print(f"Radius per fan : {r:.3f} [m]") 
print(f"Total distance taken up by fans (no spacing): {(2*r*N_fans):.3f} [m]")
T_to = 41e3
if True:
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


# ============================================================
# Performance Curve for Final Fan
# ============================================================
# Plot
# NOTE: add double axes, one for C_Q, C_T, the other for eta_prop
# NOTE: using text, add lines that show C_L at chosen takeoff and cruise conditions OF THE PROP/BLADES (NOT THE AIRCRAFT)
# TODO: ADD CLIMB CONDITION AS A VERTICAL LINE, MAKE SURE IN SLIDES TO ADD A CL, TC (HEAVY OR LIGHTLY LOADED), AND OTHER USEFUL VARIABLES 
#       AT EACH OF THE 3 CONDITIONS
adv = [
    0.00000, 0.06074, 0.12173, 0.18318, 0.24490, 0.30702, 0.36874,
    0.43030, 0.49093, 0.55015, 0.60781, 0.66338, 0.71680, 0.76781
]

Ct = [
    0.8418, 0.7777, 0.7225, 0.6744, 0.6317, 0.5924, 0.5558,
    0.5213, 0.4891, 0.4577, 0.4284, 0.4004, 0.3742, 0.3497
]

Cp = [
    0.3330, 0.3334, 0.3348, 0.3369, 0.3388, 0.3408, 0.3414,
    0.3415, 0.3404, 0.3377, 0.3339, 0.3287, 0.3225, 0.3153
]

eta_fan = [
    0.0000, 0.1417, 0.2627, 0.3666, 0.4566, 0.5337, 0.6004,
    0.6568, 0.7054, 0.7457, 0.7798, 0.8081, 0.8318, 0.8516
]

fig, ax1 = plt.subplots(figsize=(9,6))

# Left axis: Ct and Cp
ax1.plot(adv, Ct, "r-", label="C_T")
ax1.plot(adv, Cp, "m-", label="C_Q)")
ax1.set_xlabel("Advance Ratio (λ)", fontsize=15)
ax1.set_ylabel("C_T , C_Q", fontsize=15)

# Vertical lines
lambda_to = 0.12173
lambda_climb = 0.435
lambda_cruise = 0.74

ax1.axvline(lambda_to, linestyle="--", color='k', label="Takeoff")
ax1.axvline(lambda_climb, linestyle="--", color='k', label="Climb")
ax1.axvline(lambda_cruise, linestyle="--", color='k', label="Cruise")
# Right axis: efficiency
ax2 = ax1.twinx()
ax2.plot(adv, eta_fan, "g", label="η_fan")
ax2.set_ylabel("Efficiency", fontsize=15)

# Combine legends
# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()
# ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right")

# plt.title("Fan Performance")
plt.show()

# ============================================================
# Fan Blade Distribution Plot
# ============================================================

# # FAN FOR PDR
# n_blades = 5 # number of blades
# R_selected = .583 # [m]
# r = np.array([
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

# ACTUAL FAN
n_blades = 9 # number of blades
R_selected = .282 # [m]
r = np.array([
    0.054640, 0.063920, 0.073200, 0.082480, 0.091760, 0.101040,
    0.110320, 0.119600, 0.128880, 0.138160, 0.147440, 0.156720,
    0.166000, 0.175280, 0.184560, 0.193840, 0.203120, 0.212400,
    0.221680, 0.230960, 0.240240, 0.249520, 0.258800, 0.268080,
    0.277360, 0.282000
])
c = 0.775*np.array([
    0.054093, 0.066389, 0.078653, 0.090437, 0.101988, 0.112621,
    0.120846, 0.128034, 0.135261, 0.140406, 0.145471, 0.148895,
    0.151238, 0.152252, 0.152559, 0.151416, 0.148560, 0.144630,
    0.139252, 0.131152, 0.123128, 0.113221, 0.101506, 0.087691,
    0.073276, 0.060063
])
beta = np.deg2rad(np.array([
    82.5032, 80.5276, 78.5657, 76.6199, 74.6928, 72.7873,
    70.9066, 69.0537, 67.2318, 65.4437, 63.6921, 61.9793,
    60.3073, 58.6778, 57.0921, 55.5510, 54.0552, 52.6048,
    51.1997, 49.8398, 48.5243, 47.2526, 46.0238, 44.8368,
    43.6904, 43.1325
]))

# # BACKUP FAN
# n_blades = 7 # number of blades
# R_selected = .282 # [m]
# r = np.array([
#     0.054640, 0.063920, 0.073200, 0.082480, 0.091760, 0.101040,
#     0.110320, 0.119600, 0.128880, 0.138160, 0.147440, 0.156720,
#     0.166000, 0.175280, 0.184560, 0.193840, 0.203120, 0.212400,
#     0.221680, 0.230960, 0.240240, 0.249520, 0.258800, 0.268080,
#     0.277360, 0.282000
# ])
# c = R_selected* np.array([0.22000, 0.25000, 0.28000, 0.30500, 0.32500, 
#               0.33500, 0.34000, 0.34500, 0.35000, 0.35000, 
#               0.35000, 0.35000, 0.34500, 0.34000, 0.33500, 
#               0.32500, 0.31500, 0.30000, 0.28500, 0.26500, 
#               0.24000, 0.21000, 0.18000, 0.15000, 0.12000, 
#               0.10000])

# beta = np.deg2rad(np.array([
#     82.5032, 80.5276, 78.5657, 76.6199, 74.6928, 72.7873,
#     70.9066, 69.0537, 67.2318, 65.4437, 63.6921, 61.9793,
#     60.3073, 58.6778, 57.0921, 55.5510, 54.0552, 52.6048,
#     51.1997, 49.8398, 48.5243, 47.2526, 46.0238, 44.8368,
#     43.6904, 43.1325
# ]))


y_half = 1/2*c
y_full = np.concatenate([y_half, -y_half[::-1]])
r_full = np.concatenate([r, r[::-1]])
# PLOTTING 1 BLADE
plt.figure()
plt.plot(r_full, y_full, '-', color='blue')
plt.fill(r_full, y_full, alpha=0.3, color='skyblue')
plt.xlabel('Radius (m)', fontsize=15)
plt.ylabel('y (m)', fontsize=15)
plt.title('Single Fan Blade Distribution', fontsize=15)
# plt.grid(True)
# plt.axis('equal')
plt.show()


# PLOTTING ALL BLADES
y_half = 1/2*c*np.cos(beta)
y_full = np.concatenate([y_half, -y_half[::-1]])
r_full = np.concatenate([r, r[::-1]])
plt.figure(figsize=(6,6))
for i in range(n_blades):
    angle = i * 2*np.pi / n_blades
    # rotate blade coordinates
    x_rot = r_full * np.cos(angle) - y_full * np.sin(angle)
    y_rot = r_full * np.sin(angle) + y_full * np.cos(angle)
    plt.fill(x_rot, y_rot, alpha=0.4, color="skyblue", label=f'Blade {i+1}')
plt.xlabel('x (m)',fontsize=15)
plt.ylabel('y (m)',fontsize=15)
plt.title(f'Full fan with {n_blades} blades', fontsize=20)
plt.axis('equal')
plt.show()


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



# # ============================================================
# # PLOT: Attempt at sweeping over some paramters to get x_TO curves
# # ============================================================
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
