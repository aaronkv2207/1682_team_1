import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from dataclasses import dataclass
from ambiance import Atmosphere

# Imports from other subteam dependencies
from ThrustVelocity import ThrustVelocity

# from aero_workspace.aero_main import AircraftConfig, TakeoffCoeff
#eigenmodes can be calculated from JVL!  lowkey might still caculate them in this
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

@dataclass
class Landing_Gear:
    m1: float #mass of aircraft per landing gear
    m2: float #sum of wheel, tires, and axel mass
    b: float
    k1: float
    k2: float
class Aircraft:
    def __init__(self, mass, geom, aero, prop, landing, rho, ic, deltae0, deltaf0, alpha0, qbar0):
        self.mass=mass
        self.geom=geom
        self.aero=aero
        self.prop=prop
        self.landing=landing
        self.rho=rho
        self.ic=ic

        self.b=np.sqrt(geom.S*geom.AR)
        self.c=self.b/geom.AR

        self.deltae0=deltae0
        self.deltaf0=deltaf0
        self.alpha0=alpha0
        self.qbar0=qbar0


    def controls(self, t):
        '''
        ##### TAKEOFF ##########
        max 25 degrees for ailerons
        max 25-35 for rudders/elevator
        '''
        # if t<3: #chang delta_e to change based off takeoff velcoity, maybe have CL_max known
        #     delta_e=0
        # # elif t<5:
        # #     delta_e=np.radians(-10)*(t-3)/2 #ermoving pitch spike
        # elif t<6:
        #     delta_e = np.radians(-10) * (1/(1+np.exp(-(t-4)))) #negative elevator deflection -> pitch up
        # # elif t<8:
        # #     delta_e=np.radians(-10)
        # else:
        #     delta_e=0
        # if t<6:
        #     delta_f=np.radians(10)
        # elif t<9:
        #     delta_f=np.radians(10)*(1-(t-8)/2)
        # else:
        #     delta_f=0

        '''
        using all linear functions to model takeoff controls inputs
        '''
        # if t<4: #chang delta_e to change based off takeoff velcoity, maybe have CL_max known
        #     delta_e=0
        # # elif t<5:
        # #     delta_e=np.radians(-10)*(t-3)/2 #ermoving pitch spike
        # elif t<6:
        #     delta_e = np.radians(-10) * (t-4)/2 #negative elevator deflection -> pitch up
        # else:
        #     delta_e=np.radians(-10)
        # if t<10:
        #     delta_f=np.radians(10)
        # elif t<14:
        #     delta_f=np.radians(10)*(1-(t-10)/4)
        # else:
        #     delta_f=0

        '''
        ##### LANDING #####
        '''
        # small constant nose-down trim,
        '''do i want nose up or down during landing?'''
        delta_e = np.radians(1)

        # flaps deployed for landing
        if t < 8:
            delta_f = np.radians(15)
        else:
            delta_f = np.radians(10)

        # throttle slowly reducing
        throttle = .1
        throttle=0.15*np.exp(-t/10)

        return delta_e, delta_f, throttle

    def aero_coefficents(self, x, controls):
        u, w, q, theta, xe, ze=x
        delta_e, delta_f, throttle=controls

        V = np.sqrt(u**2 + w**2)
        alpha=np.arctan2(w, u)
        alpha=np.clip(alpha, np.radians(-15), np.radians(15))
        if V<1e-6:
            q_bar=0
        else:
            q_bar=(q*self.c)/(2*V)

        delta_dele=delta_e-self.deltae0
        delta_delf=delta_f-self.deltaf0
        delta_alpha=alpha-self.alpha0
        # delta_qbar=q_bar-self.qbar0

        # state_vec=np.array([delta_alpha, delta_qbar, delta_dele, delta_delf])
        state_vec=np.array([delta_alpha, delta_dele, delta_delf]) #in future ingnore CL due to elevator deflection

        CL=self.aero.CL0+self.aero.CL@state_vec
        Cm=self.aero.Cm0+self.aero.Cm@state_vec


        if ze>=0: #include ground effect factor due to rolling resistance
            phi=(16*self.geom.h/self.b)**2/(1+(16*self.geom.h/self.b)**2)
        else:
            phi=1

        Sv=self.geom.Vv*self.geom.S*self.b/self.geom.lv
        vt_c=np.sqrt(Sv/self.geom.vt_ar)
        lh=self.geom.lv+vt_c/3
        Sh=self.geom.Vh*self.geom.S*self.c/lh

        CD=self.aero.CD0+self.aero.cd_w+self.aero.cd_t*Sh/self.geom.S+phi*CL**2/(np.pi*self.geom.AR)


        '''once JVL inputs can be used'''
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
        M = Q * self.geom.S * self.c * Cm

        return X, Z, M

    def dynamics(self, t, x):
        g=9.81
        u, w, q, theta, xe, ze= x

        X, Z, M = self.forces(x, t)

        u_dot = X/self.mass.m - g*np.sin(theta) -w*q
        w_dot = Z/self.mass.m + g*np.cos(theta) + u*q
        q_dot = M/self.mass.Iy
        theta_dot = q

        u_e = u*np.cos(theta) + w*np.sin(theta)
        w_e = -u*np.sin(theta) + w*np.cos(theta)

        x_dot = u_e
        z_dot = w_e

        return [u_dot, w_dot, q_dot, theta_dot, x_dot, z_dot]

    def landing_dynamics(self, controls, t, x, z):
        #landing gear system
        g=9.81
        u, w, q, theta, xe, ze = x
        z1, z2, z1_dot, z2_dot=z

        CL, _, _= self.aero_coefficents(x, controls)
        rho=Atmosphere(h=ze*-1).density[0]
        V=np.sqrt(u**2+w**2)
        Q = 0.5 * rho * V**2
        L = Q * self.geom.S * CL

        F_shock=self.landing.b*(z1_dot-z2_dot)+self.landing.k1*(z1-z2)
        F_tire=self.landing.k2*z2
        z1_dd=(-F_shock+self.landing.m1*g-L)/self.landing.m1
        z2_dd=(F_shock-F_tire)/self.landing.m2

        return [z1_dot, z2_dot, z1_dd, z2_dd]

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

    def touchdown_event(self, t, x):
        ze=x[5]
        return ze #crossover occurs when altitude goes sfrom positive to negative


Aircraft.takeoff_event.terminal=False # doesn't stop integrating with ivp once takeoff occurs but records when event happens
Aircraft.takeoff_event.direction=1 #ensures that crossing goes from L-W<0 to L-W>0 to find zero

Aircraft.touchdown_event.terminal=True #sttop general solver at touchdown so landing dyanmcis can begin
Aircraft.touchdown_event.direction=1

'''PARAMETER INPUTS'''
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
    m=7500,
    Iy=40000
    )

prop = Propulsion(Tmax=77000)

#landing values taken from study
landing = Landing_Gear(
    m1=7500,
    m2=120,
    b=8928.98,
    k1=50002.31,
    k2=223224.59
)

aero = Aero(
    # CL0=2.96, #with flaps
    CL0=.109, #without flaps with blown
    Cm0=0.168,
    CD0=0.02,
    cd_w=0.01,
    cd_t=0.01,
    # CL=np.array([6.5, 5, 0.4, 1.2]), use if CL_qbar known
    CL=np.array([6.5, 0.4, 1.2]),
    # Cm=np.array([-1.2, -12, -1, -0.1]),
    Cm=np.array([-1.2, -1, -0.1]),
    cl=4,
    cm=-0.3,
    cx=-0.6

    # CD0=AircraftConfig.Cd0_takeoff
)

'''SOLVING IVP FOR GENERAL'''
# Initial state General
u0 = 50       # forward velocity in m/s
w0 = -2        # vertical velocity in m/s
q0 = 0        # pitch rate in rad/s
theta0 = np.radians(3)   # pitch angle in rad
x_e0 = 0      # initial x-position in meters
z_e0 = -100  # initial altitude above runwayin meters
x0 = [u0, w0, q0, theta0, x_e0, z_e0]
plane_1=Aircraft(mass, geom, aero, prop, landing, rho=1.225, ic=x0, deltae0=0, deltaf0=0, alpha0=np.radians(2), qbar0=0)

#Intergate over time span
t_span = (0, 100)
t_eval = np.linspace(0, 100, 500)

sol = solve_ivp(
    plane_1.dynamics,
    t_span,
    plane_1.ic,
    t_eval=t_eval,
    events=plane_1.touchdown_event,
    rtol=1e-8,
    atol=1e-8,
    max_step=.02
)

'''SOLVING IVP FOR LANDING'''
#Before touchdown, what is the max and min altitude, min should be about zero
print(f"Min ze: {-np.min(sol.y[5]):.3f} m")
print(f"Final ze: {-sol.y[5][-1]:.3f} m")

if sol.y_events[0].size == 0:
    print("Touchdown not detected")
    exit()

x_touch = sol.y_events[0][0]

u_td = x_touch[0]
w_td = x_touch[1]
theta_td = x_touch[3]
we_td=-u_td*np.sin(theta_td)+w_td*np.cos(theta_td)

#Landing Initial Conditions
z10=0 #no shock absorber compression initally
z20=0 #no tire compression initially
z1_dot0=we_td
z2_dot0=0
z0=[z10, z20, z1_dot0, z2_dot0]

landing_tspan=(0,10)
# landing_teval=np.linspace(0,5,300)
sol_landing=solve_ivp(
    lambda t,z: plane_1.landing_dynamics(plane_1.controls(t), t, x_touch, z),
    landing_tspan,
    z0,
    # landing_teval,
    rtol=1e-8,
    atol=1e-8,
    max_step=.02
)


# RESULTS
# x takeoff and at what time, only do this is event is takeoff
# if sol.t_events[0].size > 0:
#     t_to=sol.t_events[0][0]
#     x_to=sol.y_events[0][0][4]
#     v_to=np.sqrt(sol.y_events[0][0][0]**2+sol.y_events[0][0][1]**2)

#     print(f"Takeoff time: {t_to:.3f} s")
#     print(f"Takeoff distance: {x_to:.2f} m")
#     print(f"Takeoff velocity: {v_to:.2f} m/s")

# else:
#     print("Aircraft did not lift off.")

alpha=np.arctan2(sol.y[1], sol.y[0])
alpha = np.clip(alpha, np.radians(-15), np.radians(15))

'''PLOTS'''
# Trajectory
plt.plot(sol.y[4], -sol.y[5])
# plt.scatter(x_to, 0, label='Takeoff',c='red', marker='*')
plt.xlabel("Downrange Distance (m)")
plt.ylabel("Altitude (m)")
plt.title("Trajectory")
plt.grid(True)
plt.axis("equal")
plt.show()



#x and z Positions
# # plt.plot(sol.t, sol.y[4], label='x Position')
# plt.plot(sol.t, sol.y[4], label='x Position')
# plt.plot(sol.t, -sol.y[5], label='z Position')
# plt.legend()
# plt.grid(True)
# plt.xlabel('Time (s)')
# plt.ylabel('Position (m)')
# plt.title("Position")
# plt.show()

# Velocities
# u = sol.y[0]
# w = sol.y[1]
# v = np.sqrt(u**2 + w**2)
# plt.plot(sol.t, u, label='Forward Velocity (u)')
# plt.plot(sol.t, -w, label='Vertical Velocity (w)')
# plt.plot(sol.t, v, label='Total Velocity (v)')
# plt.xlabel("Time (s)")
# plt.ylabel("Velocity (m/s)")
# plt.title("Velocity")
# plt.legend()
# plt.grid(True)
# plt.show()

#AoA
# plt.plot(sol.t, np.degrees(alpha))
# plt.xlabel("Time (s)")
# plt.ylabel("Angle of Attack (deg)")
# plt.title("Angle of Attack")
# plt.grid(True)
# plt.show()

#Landing
z1=sol_landing.y[0]
z2=sol_landing.y[1]
compression=z1-z2

plt.plot(sol_landing.t, compression)
plt.xlabel("Time (s)")
plt.ylabel("Shock Compression (m)")
plt.title("Landing Gear Compression")
plt.grid()
plt.show()

plt.plot(sol_landing.t, z1, label="Aircraft Mass")
plt.plot(sol_landing.t, z2, label="Wheel")
plt.xlabel("Time (s)")
plt.ylabel("Displacement (m)")
plt.title("Landing Gear Motion")
plt.legend()
plt.grid()
plt.show()


F_tire = plane_1.landing.k2 * sol_landing.y[1]

plt.plot(sol_landing.t, F_tire)
plt.xlabel("Time (s)")
plt.ylabel("Gear Load (N)")
plt.title("Landing Gear Force")
plt.grid()
plt.show()
