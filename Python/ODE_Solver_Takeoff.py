import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass
from ambiance import Atmosphere
from ThrustVelocity import ThrustVelocity

@dataclass
class Mass:
    m: float
    Iy:float

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
    CL_to: np.ndarray #[CLa, CLq, CLde, CLdf]
    Cm: np.ndarray #[CMa, Cmq, Cmde, Cmddf]
    alpha_list: np.ndarray
    CL_alpha_list: np.ndarray
    flap_deg_list: np.ndarray
    CD_tot_list: np.ndarray

class Aircraft:
    def __init__(self, mass, geom, aero, prop, rho, ic, deltae0, deltaf0, alpha0, qbar0):
        self.mass=mass
        self.geom=geom
        self.aero=aero
        self.prop=prop
        self.rho=rho
        self.ic=ic
        self.b=np.sqrt(geom.S*geom.AR)
        self.c=self.b/geom.AR

        self.deltae0=deltae0
        self.deltaf0=deltaf0
        self.alpha0=alpha0
        self.qbar0=qbar0

    def theta_schedule(self, t):
        # stays 0, then ramps to 15 deg
        return np.radians(15) * (1 / (1 + np.exp(-(t -3.5)/0.8)))

    def controls(self, t):
        throttle = 1/(1+np.exp(-(t-1)/.3))

        ## ELEVATOR LOGIC
        #V1
        if t>4:
            delta_e=np.radians(-20)*1/(1+np.exp(-(t-6)/.6))
        else:
            delta_e=0

        theta = self.theta_schedule(t)

        # V2
        # (super basic controller; don't have to use rn)
        # Kp = 0.5
        # e = gamma_current - gamma_initial
        # if V>15:
        #     delta_e=Kp * e
        # else:
        #     delta_e=0

        ## FLAP LOGIC
        #delta_f=np.radians(20)*1/(1+np.exp(-(t-6)/.8))
        delta_f=np.radians(50)
        '''change this so that it's zero until velocity=15 m/s'''

        return delta_e, delta_f, throttle

    def aero_coefficents(self, x, controls, t):
        u, w, q, theta, xe, ze=x
        delta_e, delta_f, throttle=controls

        V = np.sqrt(u**2 + w**2)
        # alpha=np.clip(np.arctan2(-w, u), self.aero.alpha_list[0], self.aero.alpha_list[-1])
        u_e = u*np.cos(theta) + w*np.sin(theta)
        w_e = -u*np.sin(theta) + w*np.cos(theta)
        gamma = np.arctan2(-w_e, u_e)


        # alpha=np.clip(alpha, np.radians(-15), np.radians(15))
        # if V<1e-6:
        #     q_bar=0
        # else:
        #     q_bar=(q*self.c)/(2*V)
        # if V<15:
        #     alpha=0
        # else:
        #     alpha=np.radians(14)


        # delta_qbar=q_bar-self.qbar0

        theta = self.theta_schedule(t)
        alpha = theta - gamma
        # alpha=np.clip(np.arctan2(-w, u), self.aero.alpha_list[0], np.radians(14))
        CL=np.interp(alpha, self.aero.alpha_list, self.aero.CL_alpha_list)
        CD=np.interp(delta_f, self.aero.flap_deg_list, self.aero.CD_tot_list)

        Cm=0 #self.aero.Cm0+self.aero.Cm@state_vec


        # if ze>=0: #include ground effect factor due to rolling resistance
        #     phi=(16*self.geom.h/self.b)**2/(1+(16*self.geom.h/self.b)**2)
        # else:
        #     phi=1

        # Sv=self.geom.Vv*self.geom.S*self.b/self.geom.lv
        # vt_c=np.sqrt(Sv/self.geom.vt_ar)
        # lh=self.geom.lv+vt_c/3
        # Sh=self.geom.Vh*self.geom.S*self.c/lh

        # CD=self.aero.CD0+self.aero.cd_w+self.aero.cd_t*self.geom.Sh/self.geom.S+phi*CL**2/(np.pi*self.geom.AR)
        # CD=self.aero.CD0+self.aero.CD_ind
        # CD=self.aero.CD_profile+self.aero.CD_ind

        return CL, Cm, CD

    #Forces
    def forces(self, x, t):

        u, w, q, theta, xe, ze = x

        # print(f'altitude: {-ze}')

        controls = self.controls(t)
        delta_e, delta_f, throttle = controls

        V = np.sqrt(u**2 + w**2)

        # alpha=np.clip(np.arctan2(-w, u), self.aero.alpha_list[0], self.aero.alpha_list[-1])
        u_e = u*np.cos(theta) + w*np.sin(theta)
        w_e = -u*np.sin(theta) + w*np.cos(theta)
        gamma = np.arctan2(-w_e, u_e)
        theta = self.theta_schedule(t)
        alpha = theta - gamma

        rho=Atmosphere(h=ze*-1).density[0]
        Q = 0.5 * rho * V**2

        CL, Cm, CD= self.aero_coefficents(x, controls, t)


        L = Q * self.geom.S * CL
        D = Q * self.geom.S * CD
        T = self.prop.get_T(V)* throttle
        W = self.mass.m * 9.81

        '''TAKEOFF: if altitude for zome reason negative or zero
        add normal force and roll resistance to clamp airplane'''
        if ze >= 0:
            mu = 0.02
            Roll_resist = mu * max(W - L, 0)
            N=max(W-L, 0) #if the altitude is for some reason negative, use normal force to clamp on ground
        else:
            Roll_resist = 0

            N=0 #otherwise, no normal force

        # X and Z  are only the components of the aerodynamic force
        X = L*np.sin(alpha)+ T - D*np.cos(alpha) - Roll_resist
        Z = -L*np.cos(alpha) + D*np.sin(alpha) - N
        M = 0 #Q * self.geom.S * self.c * Cm

        return X, Z, M, L, D, T, N, Roll_resist

    def dynamics(self, t, x):
        g=9.81
        u, w, q, theta, xe, ze= x

        X, Z, M, _, _, _, _, _ = self.forces(x, t)

        u_dot = X/self.mass.m - g*np.sin(theta) -w*q
        w_dot = Z/self.mass.m + g*np.cos(theta) + u*q
        q_dot = M/self.mass.Iy
        theta_dot = q

        u_e = u*np.cos(theta) + w*np.sin(theta)
        w_e = -u*np.sin(theta) + w*np.cos(theta)

        x_dot = u_e
        z_dot = w_e

        return [u_dot, w_dot, q_dot, theta_dot, x_dot, z_dot]

    def takeoff_event(self, t, x):
        u, w = x[0], x[1]
        ze=x[5]
        V = np.sqrt(u**2 + w**2)
        rho=Atmosphere(h=ze*-1).density[0]
        Q = 0.5 * rho * V**2
        CL, Cm, CD= self.aero_coefficents(x, self.controls(t), t)

        L=Q*self.geom.S*CL
        W = self.mass.m * 9.81
        return min(L - W, -w)


Aircraft.takeoff_event.terminal=False # doesn't stop integrating with ivp once takeoff occurs but records when event happens
Aircraft.takeoff_event.direction=1 #ensures that crossing goes from L-W<0 to L-W>0 to find zero

'''PARAMETER INPUTS'''
geom = Geometry(
    S=50.1,
    AR=8,
    h=5.5,
    lv=10,
    Vv=.1,
    vt_ar=1.2,
    Vh=1.05,
    Sh=11.98
    # S=AircraftConfig.s_ref
    # AR=AircraftConfig.AR
    )

mass = Mass(
    m=7400,
    Iy=40000
    )

prop = ThrustVelocity()

aero = Aero(
    CL0=[2.96,.109], #[with flaps, w/out flaps]
    # CL0 = 1.2, #without flaps
    # CL0=.109, #without flaps with blown
    Cm0=0.168,
    CD0=0.02,
    cd_w=0.01,
    cd_t=0.01,
    CD_ind=1.645,
    CD_profile=0.02,
    # CL=np.array([6.5, 5, 0.4, 1.2]), use if CL_qbar known
    CL_to=np.array([6.5, 0.4, 6.5]), #change for takeoff to be 3 for CLdeltaf
    # Cm=np.array([-1.2, -12, -1, -0.1]),
    Cm=np.array([-1.2, -1, -0.1]),
    # Cm=np.array([-0.1, -0.2, -0.05]),
    alpha_list=np.radians([-15, -13,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,13,15]),
    CL_alpha_list=[1.99602,2.32163,2.64366,2.96144,3.2743,3.58159,3.88268,4.17698,4.46391,4.74295,5.0136,5.2754,5.52793,5.77082,6.00373,6.22638],
    flap_deg_list=np.radians([40, 50, 60]),
    CD_tot_list=[0.31508, 0.3842, 0.45865] #CD_prof+CD_ind
)

'''SOLVING IVP FOR GENERAL'''
# Initial state General
u0_to = 0       # forward velocity in m/s
w0_to = 0        # vertical velocity in m/s
q0_to = 0        # pitch rate in rad/s
theta0_to = 0   # pitch angle in rad
x_e0_to = 0      # initial x-position in meters
z_e0_to = 0  # initial altitude above runwayin meters


x0 = [u0_to, w0_to, q0_to, theta0_to, x_e0_to, z_e0_to]
plane_1=Aircraft(mass, geom, aero, prop,rho=1.225, ic=x0, deltae0=0, deltaf0=0, alpha0=np.radians(2), qbar0=0) #change alpha_0 to 0

#Intergate over time span
t = 10.5
t_span = (0, t)
t_eval = np.linspace(0, t, 1000)

event=plane_1.takeoff_event


sol = solve_ivp(
    plane_1.dynamics,
    t_span,
    plane_1.ic,
    t_eval=t_eval,
    events=event,
    # rtol=1e-8,
    # atol=1e-8,
    max_step=.01
)


# x takeoff and at what time, only do this is event is takeoff
if sol.t_events[0].size > 0:
    t_to=sol.t_events[0][0]
    x_to=sol.y_events[0][0][4]
    theta_to=plane_1.theta_schedule(t_to)
    print(f"theta at takeoff: {np.degrees(theta_to):.2f} deg")
    # theta_to=sol.y_events[0][0][3]
    u_to=sol.y_events[0][0][0]
    w_to=sol.y_events[0][0][1]

    u_eto = u_to*np.cos(theta_to) + w_to*np.sin(theta_to)
    w_eto = -u_to*np.sin(theta_to) + w_to*np.cos(theta_to)
    gamma_to = np.arctan2(-w_eto, u_eto)
    print(f"gamma at takeoff: {np.degrees(gamma_to):.2f} deg")
    alpha_to=theta_to-gamma_to
    v_to=np.sqrt(sol.y_events[0][0][0]**2+sol.y_events[0][0][1]**2)

    print(f"Takeoff time: {t_to:.3f} s")
    print(f"Takeoff distance: {x_to:.2f} m")
    print(f"Takeoff velocity: {v_to:.2f} m/s")
    print(f"Takeoff AoA: {np.degrees(alpha_to):.2f} deg")

else:
    print("Aircraft did not lift off.")


#alpha = np.clip(alpha, np.radians(-15), np.radians(15))

'''PLOTS'''
# Trajectory
plt.plot(sol.y[4], -sol.y[5])
plt.scatter(x_to, 0, label='Takeoff',c='red', marker='*')
plt.axhline(y=15.24, color='red', linestyle='--', label='50 ft clearing')
plt.xlabel("Downrange Distance (m)")
plt.ylabel("Altitude (m)")
plt.title("Trajectory")
plt.grid(True)
plt.axis("equal")
plt.legend(loc='lower right')
plt.show()

#x and z Positions
fig, ax, =plt.subplots(2)
ax[0].plot(sol.t, sol.y[4], label='x Position')
ax[0].plot(sol.t, -sol.y[5], label='z Position')
ax[0].legend(loc='upper left')
ax[0].grid(True)
# ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('Position (m)')
ax[0].set_title("Position")

# Velocities
u = sol.y[0]
w = sol.y[1]
v = np.sqrt(u**2 + w**2)
ax[1].plot(sol.t, u, label='Forward Velocity (u)')
ax[1].plot(sol.t, -w, label='Vertical Velocity (w)')
ax[1].plot(sol.t, v, label='Total Velocity (v)')
ax[1].set_xlabel("Time (s)")
ax[1].set_ylabel("Velocity (m/s)")
ax[1].set_title("Velocity")
ax[1].legend(loc='upper left')
ax[1].grid(True)
plt.show()


CL_clean = []
L_clean=[]
D_clean=[]
T_clean=[]
N_clean=[]
Roll_resist_clean=[]
alpha_clean=[]
theta_clean=[]

alpha_list=np.radians([-15, -13,-11,-9,-7,-5,-3,-1,1,3,5,7,9,11,13,15])
for i in range(len(sol.t)):
    x = sol.y[:, i] #grab current state vector)
    u, w, ze, theta =x[0], x[1], x[5], x[3]
    # alpha_i=np.clip(np.arctan2(-w, u), alpha_list[0], alpha_list[-1])
    u_e = u*np.cos(theta) + w*np.sin(theta)
    w_e = -u*np.sin(theta) + w*np.cos(theta)
    gamma = np.arctan2(-w_e, u_e)
    theta = plane_1.theta_schedule(t)
    alpha = theta - gamma


    t = sol.t[i]
    controls = plane_1.controls(t)
    CL, _, _ = plane_1.aero_coefficents(x, controls, t)
    CL_clean.append(CL)
    _, _, _, L, D, T, N, Roll_resist=plane_1.forces(x, t)

    alpha_clean.append(alpha)
    theta_clean.append(theta)
    L_clean.append(L/10**3)
    D_clean.append(D/10**3)
    T_clean.append(T/10**3)
    N_clean.append(N/10**3)
    Roll_resist_clean.append(Roll_resist/10**3)

plt.plot(sol.t, np.degrees(alpha_clean), label='AoA')
plt.plot(sol.t, np.degrees(theta_clean), label='Theta')
plt.plot(sol.t, np.degrees(theta_clean)-np.degrees(alpha_clean), label='Climb Angle')
plt.xlabel("Time (s)")
plt.ylabel("Angles [deg]")
plt.ylim(-30, 50)
plt.title("Angles of Interest over time")
plt.legend(loc='lower right')
plt.grid(True)
plt.show()

plt.plot(sol.t, CL_clean)
plt.xlabel("Time (s)")
plt.ylabel("Lift Coefficient (CL)")
plt.title("Lift Coefficient over Time")
plt.grid(True)
plt.show()

# Forces
plt.axvline(x=t_to, label='Takeoff',c='red', linestyle='--')
plt.plot(sol.t, L_clean, label='Lift')
plt.plot(sol.t, D_clean, label='Drag')
plt.plot(sol.t, T_clean, label='Thrust')
plt.xlabel("Time (s)")
plt.ylabel("Force (kN)")
plt.title("Forces over Time")
plt.legend(loc='upper left')
plt.grid(True)
plt.show()

plt.axvline(x=t_to, label='Takeoff',c='red', linestyle='--')
plt.plot(sol.t, N_clean, label='Normal Force')
plt.plot(sol.t, Roll_resist_clean, label='Rolling Resistance')
plt.xlabel("Time (s)")
plt.ylabel("Force (kN)")
plt.title("Forces over Time")
plt.legend(loc='lower right')
plt.grid(True)
plt.show()
