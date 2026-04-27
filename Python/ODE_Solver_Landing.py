from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from ambiance import Atmosphere
from scipy.integrate import solve_ivp
from ThrustVelocity import ThrustVelocity


# ---------------- DATA STRUCTURES ----------------
@dataclass
class Mass:
    m: float
    Iy: float


@dataclass
class Geometry:
    S: float
    AR: float
    h: float
    lv: float
    Vv: float
    vt_ar: float
    Vh: float
    Sh: float


@dataclass
class Aero:
    CL0: float
    Cm0: float
    CD0: float
    CD_ind: float
    CD_profile: float
    cd_w: float
    cd_t: float
    CL_to: np.ndarray
    CL_td: np.ndarray
    Cm: np.ndarray
    alpha_list: np.ndarray
    CL_alpha_list: np.ndarray
    flap_deg_list: np.ndarray
    CD_tot_list: np.ndarray


@dataclass
class Landing_Gear:
    m1: float
    m2: float
    b: float
    k1: float
    k2: float


# ---------------- AIRCRAFT ----------------
class Aircraft:
    def __init__(self, phase, mass, geom, aero, prop, landing, ic,
                 deltae0=0, deltaf0=0, alpha0=np.radians(2)):

        self.phase = phase
        self.mass = mass
        self.geom = geom
        self.aero = aero
        self.prop = prop
        self.landing = landing
        self.ic = ic

        self.deltae0 = deltae0
        self.deltaf0 = deltaf0
        self.alpha0 = alpha0

    # ---------------- CONTROLS ----------------
    def controls(self, t):
        delta_e = np.radians(5)

        if t < 8:
            delta_f = np.radians(15)
        else:
            delta_f = np.radians(10)

        throttle = 0.15 * np.exp(-t / 10)

        return delta_e, delta_f, throttle

    # ---------------- AERO ----------------
    def aero_coefficents(self, x, controls, t):
        u, w, q, theta, xe, ze = x
        delta_e, delta_f, _ = controls

        # PURE LANDING AOA (relative wind)
        alpha = np.arctan2(w, u)

        CL0 = self.aero.CL0[1]
        CL_vec = self.aero.CL_td

        delta_alpha = alpha - self.alpha0
        delta_dele = delta_e - self.deltae0
        delta_delf = delta_f - self.deltaf0

        state_vec = np.array([delta_alpha, delta_dele, delta_delf])

        CL = CL0 + CL_vec @ state_vec
        CD = self.aero.CD0 + CL**2 / (np.pi * self.geom.AR)

        Cm = 0

        return CL, Cm, CD

    # ---------------- FORCES ----------------
    def forces(self, x, t):

        u, w, q, theta, xe, ze = x
        delta_e, delta_f, throttle = self.controls(t)

        V = np.sqrt(u**2 + w**2)
        alpha = np.arctan2(w, u)

        rho = Atmosphere(h=-ze).density[0]
        Q = 0.5 * rho * V**2

        CL, _, CD = self.aero_coefficents(x, (delta_e, delta_f, throttle), t)

        L = Q * self.geom.S * CL
        D = Q * self.geom.S * CD
        T = self.prop.get_T(V) * throttle
        W = self.mass.m * 9.81

        # ground model
        if ze >= 0:
            mu = 0.02
            N = max(W - L, 0)
            Roll = mu * N
        else:
            N = 0
            Roll = 0

        # force components in body frame
        X = L * np.sin(alpha) + T - D * np.cos(alpha) - Roll
        Z = -L * np.cos(alpha) + D * np.sin(alpha) - N

        return X, Z, L, D, T, N, Roll

    # ---------------- AIRCRAFT DYNAMICS ----------------
    def dynamics(self, t, x):

        u, w, q, theta, xe, ze = x
        g = 9.81

        X, Z, *_ = self.forces(x, t)

        u_dot = X / self.mass.m - g * np.sin(theta) - w * q
        w_dot = Z / self.mass.m + g * np.cos(theta) + u * q
        q_dot = 0  # no pitch dynamics in landing model (controlled/trimmed assumption)
        theta_dot = q

        u_e = u * np.cos(theta) + w * np.sin(theta)
        w_e = -u * np.sin(theta) + w * np.cos(theta)

        x_dot = u_e
        z_dot = w_e

        return [u_dot, w_dot, q_dot, theta_dot, x_dot, z_dot]

    # ---------------- TOUCHDOWN EVENT ----------------
    def touchdown_event(self, t, x):
        return x[5]


Aircraft.touchdown_event.terminal = True
Aircraft.touchdown_event.direction = 1


# ---------------- LANDING GEAR MODEL ----------------
def landing_dynamics(plane, x_touch, z_init, t_span):

    def dyn(t, z):
        z1, z2, z1_dot, z2_dot = z

        u, w = x_touch[0], x_touch[1]
        V = np.sqrt(u**2 + w**2)
        rho = Atmosphere(h=-x_touch[5]).density[0]
        CL, _, _ = plane.aero_coefficents(x_touch,
                                           plane.controls(t),
                                           t)

        L = 0.5 * rho * V**2 * plane.geom.S * CL

        F_shock = plane.landing.b * (z1_dot - z2_dot) + plane.landing.k1 * (z1 - z2)
        F_tire = plane.landing.k2 * z2

        g = 9.81

        z1_dd = (-F_shock + plane.landing.m1 * g - L) / plane.landing.m1
        z2_dd = (F_shock - F_tire) / plane.landing.m2

        return [z1_dot, z2_dot, z1_dd, z2_dd]

    return solve_ivp(dyn, t_span, z_init, max_step=0.02)


# ---------------- PARAMETERS ----------------
geom = Geometry(50.1, 8, 5.5, 10, 0.1, 1.2, 1.05, 11.98)
mass = Mass(7400, 40000)
prop = ThrustVelocity()

landing = Landing_Gear(3300, 120, 3.5e4, 3e5, 5e6)

aero = Aero(
    CL0=[2.96, 0.109],
    Cm0=0.168,
    CD0=0.02,
    cd_w=0.01,
    cd_t=0.01,
    CD_ind=1.645,
    CD_profile=0.02,
    CL_to=np.array([6.5, 0.4, 6.5]),
    CL_td=np.array([6.5, 0.4, 1.2]),
    Cm=np.array([-1.2, -1, -0.1]),
    alpha_list=np.radians([-15, -13,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,13,15]),
    CL_alpha_list=[1.99,2.32,2.64,2.96,3.27,3.58,3.88,4.17,4.46,4.74,5.01,5.27,5.52,5.77,6.00,6.22],
    flap_deg_list=np.radians([40, 50, 60]),
    CD_tot_list=[0.31, 0.38, 0.45]
)

# ---------------- INITIAL CONDITION ----------------
x0 = [30, 1, 0, np.radians(-2), 0, -100]

plane = Aircraft("landing", mass, geom, aero, prop, landing, x0)

# ---------------- SIMULATION ----------------
t_span = (0, 12)
t_eval = np.linspace(0, 12, 1000)

sol = solve_ivp(
    plane.dynamics,
    t_span,
    x0,
    t_eval=t_eval,
    events=plane.touchdown_event,
    max_step=0.01
)

# ---------------- TOUCHDOWN ----------------
x_touch = sol.y_events[0][0]

u_td, w_td, theta_td = x_touch[0], x_touch[1], x_touch[3]
we_td = -u_td*np.sin(theta_td) + w_td*np.cos(theta_td)

z0 = [0, 0, 0.5*we_td, 0]

sol_landing = landing_dynamics(
    plane,
    x_touch,
    z0,
    (0, 5)
)

# ---------------- BASIC OUTPUT ----------------
print("Touchdown completed")
print(f"Impact vertical velocity: {we_td:.3f} m/s")

# ---------------- PLOT ----------------
#Trajectory
plt.plot(sol.y[4], -sol.y[5])
plt.xlabel("Downrange Distance (m)")
plt.ylabel("Altitude (m)")
plt.title("Trajectory")
plt.grid(True)
plt.axis("equal")
plt.legend(loc='lower right')
plt.show()


#Landing Gear Force
F_shock=plane.landing.b*(sol_landing.y[2]-sol_landing.y[3])+plane.landing.k1*(sol_landing.y[0]-sol_landing.y[1])
plt.plot(sol_landing.t, F_shock/10**3, label='Current loading')
plt.axhline(y=151, color='r', linestyle='--', label='2 g')
plt.xlabel("Time (s)")
plt.ylabel("Gear Load (kN)")
plt.title("Landing Gear Force")
plt.legend()
plt.grid()
plt.show()



#Landing Gear Compression
plt.plot(sol_landing.t, sol_landing.y[0] - sol_landing.y[1])
plt.axhline(y=0.4, color='r', linestyle='--', label='0.4 m goal')
plt.title("Landing Gear Compression")
plt.xlabel('Times (s)')
plt.ylabel('Compression (m)')
plt.legend()
plt.grid()
plt.show()
