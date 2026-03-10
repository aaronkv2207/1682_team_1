import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass
from ambiance import Atmosphere

# from aero_workspace.aero_main import AircraftConfig, TakeoffCoeff

#eigenmodes can be calculated from JVL!
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

@dataclass
class Aero:
    CL0: float
    Cm0: float
    CD0: float
    cd_w: float
    cd_t: float
    CL: np.ndarray #[CLa, CLq, CLde, CLdf]
    Cm: np.ndarray #[CMa, Cmq, Cmde, Cmddf]
    cl: float
    cm: float
    cx:float

@dataclass
class Propulsion:
    Tmax: float
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

    def controls(self, t):
        if t<3: #chang delta_e to change based off takeoff velcoity, maybe have CL_max known
            delta_e=0
        # elif t<5:
        #     delta_e=np.radians(-10)*(t-3)/2 #ermoving pitch spike
        else:
            delta_e=np.radians(-10)
        if t<20:
            delta_f=np.radians(10)
        elif t<25:
            delta_f=np.radians(10)*(1-(t-20)/5)
        else:
            delta_f=0

        throttle=1-np.exp(-t/2)
        return delta_e, delta_f, throttle

    def aero_coefficents(self, x, controls):
        u, w, q, theta, xe, ze=x
        delta_e, delta_f, throttle=controls

        V = np.sqrt(u**2 + w**2)
        alpha=np.arctan2(w, u)
        alpha=np.clip(alpha, np.radians(-15), np.radians(15))
        # if V<1e-6:
        #     q_bar=0
        # else:
        #     q_bar=(q*self.c)/(2*V)

        # delta_dele=delta_e-self.deltae0
        # delta_delf=delta_f-self.deltaf0
        # delta_alpha=alpha-self.alpha0
        # delta_qbar=q_bar-self.qbar0

        # state_vec=np.array([delta_alpha, delta_qbar, delta_dele, delta_delf])

        # CL=self.aero.CL0+self.aero.CL@state_vec
        # Cm=self.aero.Cm0+self.aero.Cm@state_vec


        # if ze<=0: #include ground effect factor due to rolling resistance
        #     phi=(16*self.geom.h/self.b)**2/(1+(16*self.geom.h/self.b)**2)
        # else:
        #     phi=1

        # Sv=self.geom.Vv*self.geom.S*self.b/self.geom.lv
        # vt_c=np.sqrt(Sv/self.geom.vt_ar)
        # lh=self.geom.lv+vt_c/3
        # Sh=self.geom.Vh*self.geom.S*self.c/lh

        # CD=self.aero.CD0+self.aero.cd_w+self.aero.cd_t*Sh/self.geom.S+phi*CL**2/(np.pi*self.geom.AR)
        CL=2*np.pi*alpha+0.3
        CL=np.clip(CL, -1.5, 1.5)

        Cm=0
        CD=CL**2/(np.pi*self.geom.AR)+.02


        # CL=TakeoffCoeff.CL_alpha(alpha)+TakeoffCoeff.CL_flap(delta_f)+TakeoffCoeff.CL_velocity(V)
        # Cm=TakeoffCoeff.CM_alpha(alpha)+TakeoffCoeff.CM_flap(delta_f)+TakeoffCoeff.CM_velocity(V)
        # CD=TakeoffCoeff.CD_alpha(alpha)+TakeoffCoeff.CD_flap(delta_f)+TakeoffCoeff.CD_velocity(V)
        return CL, Cm, CD

    #Forces
    def forces(self, x, t):

        u, w, q, theta, xe, ze = x

        controls = self.controls(t)
        delta_e, delta_f, throttle = controls

        V = np.sqrt(u**2 + w**2)

        alpha = np.arctan2(w, u)
        alpha=np.clip(alpha, np.radians(-15), np.radians(15))


        rho=Atmosphere(h=ze*-1).density[0]
        Q = 0.5 * rho * V**2

        CL, Cm, CD= self.aero_coefficents(x, controls)


        L = Q * self.geom.S * CL
        D = Q * self.geom.S * CD

        T = self.prop.Tmax * throttle
        # # omega=200
        # # T_W=0.3
        # T=omega*T_W

        W = self.mass.m * 9.81

        # if ze <= 0:
        #     mu = 0.02
        #     R = mu * max(W - L, 0)
        # else:
        #     R = 0

        # X and Z  are only the components of the aerodynamic force
        # X = (T_D)*np.cos(alpha) + L*np.sin(alpha) - R
        X = -L*np.sin(alpha)+ T - D*np.cos(alpha)
        Z = -L*np.cos(alpha) - D*np.sin(alpha)

        M = Q * self.geom.S * self.c * Cm

        return X, Z, M

    def dynamics(self, t, x):
        u, w, q, theta, xe, ze = x
        # print(f'u: {u}, w: {w}')

        X, Z, M = self.forces(x, t)

        u_dot = X/self.mass.m - 9.81*np.sin(theta) -w*q
        w_dot = Z/self.mass.m + 9.81*np.cos(theta) + u*q
        q_dot = M/self.mass.Iy
        theta_dot = q

        u_e = u*np.cos(theta) + w*np.sin(theta)
        w_e = -u*np.sin(theta) + w*np.cos(theta)

        x_dot = u_e
        z_dot = w_e

        # For simplicity, if we are on the ground and not lifting off:
        L_minus_W = self.takeoff_event(t, x)
        if L_minus_W < 0:
            w_dot = 0
            # Force pitch rate to zero if you don't want it to 'nose dive' into the dirt
            if theta <= 0 and q < 0:
                q_dot = 0
                theta_dot = 0


        # print(u_dot, -w_dot, q_dot, theta_dot, x_dot, -z_dot)

        return [u_dot, w_dot, q_dot, theta_dot, x_dot, z_dot]




    def takeoff_event(self, t, x):
        u, w = x[0], x[1]
        ze=x[5]
        V = np.sqrt(u**2 + w**2)
        rho=Atmosphere(h=ze*-1).density[0]
        Q = 0.5 * rho * V**2
        CL, Cm, CD= self.aero_coefficents(x, self.controls(t))

        L=Q*self.geom.S*CL
        W = self.mass.m * 9.81
        return L - W



Aircraft.takeoff_event.terminal=False # doesn't stop integrating with ivp once takeoff occurs but records when event happens
Aircraft.takeoff_event.direction=1 #ensures that crossing goes from L-W<0 to L-W>0 to find zero

# Initial state
# UNKOWN initial conditions
u0 = 0       # forward velocity in m/s
w0 = 0        # vertical velocity in m/s
q0 = 0        # pitch rate in rad/s
theta0 = np.radians(5)   # pitch angle in rad
x_e0 = 0      # initial x-position in meters
z_e0 = 0  # initial altitude in meters
x0 = [u0, w0, q0, theta0, x_e0, z_e0]


geom = Geometry(
    S=49.5,
    AR=8,
    h=5.5,
    lv=10,
    Vv=.1,
    vt_ar=1.2,
    Vh=1.05,
    # S=AircraftConfig.s_ref
    # AR=AircraftConfig.AR

    )

mass = Mass(
    m=55622.7/9.81,
    Iy=40000
    )

prop = Propulsion(Tmax=77000)

aero = Aero(
    CL0=0.3,
    Cm0=0,
    CD0=0.013,
    cd_w=0,
    cd_t=0.01,
    CL=np.array([6.5, 5, 0.4, 1.2]),
    Cm=np.array([-1.2, -12, -1, -0.1]),
    cl=4,
    cm=-0.3,
    cx=-0.6

    # CD0=AircraftConfig.Cd0_takeoff
)

plane_1=Aircraft(mass, geom, aero, prop, rho=1.225, ic=x0, deltae0=0, deltaf0=0, alpha0=np.radians(2), qbar0=0)

#Intergate for 20 seconds
t_span = (0, 10)
t_eval = np.linspace(0, 10, 1000)

sol = solve_ivp(
    plane_1.dynamics,
    t_span,
    plane_1.ic,
    t_eval=t_eval,
    events=plane_1.takeoff_event,
    rtol=1e-8,
    atol=1e-8,
    max_step=.02
)

# sol.y[i] gives time series of ith state

# RESULTS
# x takeoff and at what time
if sol.t_events[0].size > 0:
    t_to=sol.t_events[0][0]
    x_to=sol.y_events[0][0][4]
    v_to=np.sqrt(sol.y_events[0][0][0]**2+sol.y_events[0][0][1]**2)

    print(f"Takeoff time: {t_to:.3f} s")
    print(f"Takeoff distance: {x_to:.2f} m")
    print(f"Takeoff velocity: {v_to:.2f} m/s")

else:
    print("Aircraft did not lift off.")

alpha=np.arctan2(sol.y[1], sol.y[0])
alpha = np.clip(alpha, np.radians(-15), np.radians(15))


# print(f'w_dot: {w_dot_list}')
# print(f'w: {sol.y[1]}') #without a clamp the vertical acceleration of gravity allows the plane to go to negative z
#if w_dot>0 then acceleration will be downwards
# --- Plot results ---
# Positions

plt.plot(sol.y[4], -sol.y[5])
# plt.scatter(x_to, 0, label='Takeoff')
plt.xlabel("Downrange DIstance (m)")
plt.ylabel("Altitude (m)")
plt.title("Takeoff Trajectory")
plt.grid(True)
plt.axis("equal")
plt.show()

# Velocities
u = sol.y[0]
w = sol.y[1]
v = np.sqrt(u**2 + w**2)

# plt.plot(sol.t, sol.y[4], label='x Position')
plt.plot(sol.t, u, label='Forward Velocity (u)')
plt.plot(sol.t, -w, label='Vertical Velocity (w)')
plt.plot(sol.t, v, label='Total Velocity (v)')
plt.xlabel("Time (s)")
plt.ylabel("Velocity (m/s)")
plt.title("Kinematics")
plt.legend()
plt.grid(True)
plt.show()

plt.plot(sol.t, np.degrees(alpha))
plt.xlabel("Time (s)")
plt.ylabel("Angle of Attack (deg)")
plt.title("Angle of Attack")
plt.grid(True)
plt.show()
