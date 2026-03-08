import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class Aircraft:
    def __init__(
        self,
        m_total,
        I_y,
        rho,
        S,
        St,
        R,
        omega,
        CL0,
        CLa,
        CLq,
        CL_de,
        CL_df,
        cx,
        cx_prime,
        Cm0,
        Cma,
        Cmq,
        Cm_de,
        Cm_df,
        cd_profile,
        cd_w,
        cd_t,
        CT,
        static_condition,
        alpha_0,
        h,
        b,
        c_bar,
        initial_conditions,
    ):
        self.m = m_total
        self.I_y = I_y
        self.S = S
        self.St = St
        self.R = R
        self.rho = rho
        self.omega = omega
        self.CL0 = CL0
        self.CLa = CLa
        self.CLq = CLq
        self.CL_de = CL_de
        self.CL_df = CL_df
        self.cx = cx
        self.cx_prime = cx_prime
        self.Cm0 = Cm0
        self.Cma = Cma
        self.Cmq = Cmq
        self.Cm_de = Cm_de
        self.Cm_df = Cm_df
        self.cd_profile = cd_profile
        self.cd_w = cd_w
        self.cd_t = cd_t
        self.CT = CT
        self.static = static_condition
        self.alpha_0 = alpha_0
        self.h = h
        self.b = b
        self.AR = b**2 / S
        self.c_bar = c_bar
        self.ic = initial_conditions

    def aircraft_dynamics(self, t, x):
        """
        x = [u, w, q, theta] state vector
        params = dictionary of aircraft parameters
        """
        # State vector
        u, w, q, theta, xe, ze = x

        # Kinematics
        alpha = np.arctan2(w, u)
        V = np.sqrt(u**2 + w**2)
        Q = 0.5 * self.rho * V**2

        # Transform velocities to earth frame
        u_e = u * np.cos(theta) + w * np.sin(theta)
        w_e = w * np.cos(theta) - u * np.sin(theta)
        gamma = np.arctan2(-w_e, u_e)

        # Forces
        # Thrust / drag
        # How do I fix this, de coupling T and D to get rid of cx
        # Are control surfaces inputs into thrust, drag, or both
        if self.static:
            T_D = -Q * self.S * self.cx
        elif ze == 0:  # if at takeoff, consider phi component on induced drag
            phi = (16 * self.h / self.b) ** 2 / (1 + (16 * self.h / self.b) ** 2)
            T_D = -0.5 * self.rho * (self.omega * self.R**2) * (
                self.cx_prime
            ) + phi * self.cl**2 / (np.pi * self.AR)
        else:
            T_D = -0.5 * self.rho * (self.omega * self.R**2) * self.cx_prime

        # Lift
        L = Q * self.S * self.CL(t, x)

        # Rolling resistance at takeoff
        if ze == 0:
            W = self.m * 9.81
            mu_r = 0.02  # coefficent of friction
            R = mu_r * (W - L)  # eventually this will decrease to 0

        # Forces in body frame
        if ze == 0:
            X = T_D * np.cos(alpha) + L * np.sin(alpha) - R
        else:
            X = T_D * np.cos(alpha) + L * np.sin(alpha)
        Z = T_D * np.sin(alpha) - L * np.cos(alpha)

        # Moment
        M = Q * self.S * self.Cm(t, x)  # added damping term, will go away in real model

        # Equations of motion
        u_dot = (1 / self.m) * X - 9.81 * np.sin(theta)
        w_dot = (1 / self.m) * Z + 9.81 * np.cos(theta) + u * q
        q_dot = 1 / self.I_y * M
        theta_dot = q
        x_dot = u_e
        z_dot = w_e

        return [u_dot, w_dot, q_dot, theta_dot, x_dot, z_dot]

    def takeoff_event(self, t, x):
        u, w = x[0], x[1]
        V = np.sqrt(u**2 + w**2)
        Q = 0.5 * self.rho * V**2
        L = Q * self.CL(t, x) * self.S
        W = self.m * 9.81
        return L - W

    def CL(self, t, x):
        # Control surface inputs
        u, w, q, theta, xe, ze = x
        alpha = np.arctan2(w, u)
        (
            ui,
            wi,
            qi,
            thetai,
            xei,
            zei,
        ) = self.ic
        Vi = np.sqrt(ui**2 + wi**2)
        q_bari = (qi * self.c_bar) / (2 * Vi)
        V = np.sqrt(u**2 + w**2)
        q_bar = (q * self.c_bar) / (2 * V)

        del_e = self.elevator(t)
        del_f = self.flaps(t)
        del_alpha = alpha - self.alpha_0
        del_qbar = q_bar - q_bari

        return (
            self.CL0
            + self.CLa * del_alpha
            + self.CLq * del_qbar
            + self.CL_df * del_f
            + self.CL_de * del_e
        )

    def Cm(self, t, x):
        u, w, q, theta, xe, ze = x
        alpha = np.arctan2(w, u)
        (
            ui,
            wi,
            qi,
            thetai,
            xei,
            zei,
        ) = self.ic
        Vi = np.sqrt(ui**2 + wi**2)
        q_bari = (qi * self.c_bar) / (2 * Vi)
        V = np.sqrt(u**2 + w**2)
        q_bar = (q * self.c_bar) / (2 * V)

        # Include nominal conditions
        del_e = self.elevator(t)
        del_f = self.flaps(t)
        del_alpha = alpha - self.alpha_0
        del_qbar = q_bar - q_bari

        return (
            self.Cm0
            + self.Cma * del_alpha
            + self.Cmq * del_qbar
            + self.Cm_df * del_f
            + self.Cm_de * del_e
        )

    def CD(self, t, x):
        return (
            self.cd_profile / self.S
            + self.cd_w
            + self.cd_t * self.St / self.S
            + self.CL(t, x) ** 2 / (np.pi * self.AR)
        )

    """
    Modeling elevator surface delfection over time.
    Given a time input, elevator deflection will change, which will go into Cm, CL, etc caclulations"""

    def elevator(self, t):
        if t < 5:
            return 0
        else:
            return np.radians(-5)  # modellig pull up at takeoff after 5 sec

    def flaps(self, t):
        if t < 20:
            return np.radians(10)
        elif t < 25:
            return np.radians(10) * (1 - (t - 20) / 5)  # linearly decreasing funnction
        else:
            return 0


Aircraft.takeoff_event.terminal = False  # doesn't stop integrating with ivp once takeoff occurs but records when event happens
Aircraft.takeoff_event.direction = (
    1  # ensures that crossing goes from L-W<0 to L-W>0 to find zero
)

# Initial state
# Example initial conditions
u0 = 40  # forward velocity in m/s
w0 = 0  # vertical velocity in m/s
q0 = 0  # pitch rate in rad/s
theta0 = 0.05  # pitch angle in rad
x_e0 = 0  # initial x-position in meters
z_e0 = 0  # initial altitude in meters
x0 = [u0, w0, q0, theta0, x_e0, z_e0]


# Create new versions of aircraft

plane_1 = Aircraft(
    m_total=5670,
    I_y=500,
    rho=1.225,
    S=10,
    St=...,
    R=2,
    omega=100,  # propeller rotation rate
    CL0=0.3,
    CLa=...,
    CLq=...,
    CL_de=...,
    CL_df=...,
    cx=0.3,  # eventually cx either gets split into compoents or becomes a function  or a function of cl, dt, q_bar, df, Re, Ma
    cx_prime=0.2,
    Cm0=0,
    Cma=...,
    Cmq=...,
    Cm_de=...,
    Cm_df=...,
    cd_profile=...,
    cd_w=...,
    cd_t=...,
    CT=...,
    static_condition=False,
    alpha_0=0,
    h=5.5,
    b=19.8,
    c_bar=4,  # mean aerodynamic cord
    initial_conditions=x0,
)

# Intergate for 20 seconds
t_span = (0, 20)
t_eval = np.linspace(0, 20, 1000)

sol = solve_ivp(
    plane_1.aircraft_dynamics,
    t_span,
    plane_1.ic,
    t_eval=t_eval,
    events=plane_1.takeoff_event,
)

# sol.y[i] gives time series of ith state

# RESULTS
# x takeoff and at what time
if sol.t_events[0].size > 0:
    t_to = sol.t_events[0][0]
    x_to = sol.y_events[0][0][4]

    print(f"Takeoff time: {t_to:.3f} s")
    print(f"Takeoff distance: {x_to:.2f} m")

else:
    print("Aircraft did not lift off.")

# --- Plot results ---
# Positions

plt.plot(sol.y[4], sol.y[5])
plt.xlabel("x Position (m)")
plt.ylabel(" z Position (m)")
plt.title("Trajectory")
plt.grid(True)
plt.show()

# Velocities
u = sol.y[0]
w = sol.y[1]
v = np.sqrt(u**2 + w**2)

plt.plot(sol.t, v)
plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")
plt.title("Airspeed vs Time")
plt.grid(True)
plt.show()
