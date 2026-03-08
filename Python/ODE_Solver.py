import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

class Aircraft:
    def __init__(self, m_total, I_y, rho, S, St, R, omega, CL0, CLa,CLq, CL_de, CL_df, Cm0, Cma, Cmq, Cm_de, Cm_df, cd_w, cd_t, CD0, TCprime, static_condition, alpha_0, h, AR,deltae0, deltaf0, initial_conditions):
        self.m=m_total
        self.I_y=I_y
        self.S=S
        self.St=St
        self.R=R
        self.rho=rho
        self.omega=omega
        self.CL0=CL0
        self.CLa=CLa
        self.CLq=CLq
        self.CL_de=CL_de
        self.CL_df=CL_df
        # self.cx=cx
        # self.cx_prime=cx_prime
        self.Cm0=Cm0
        self.Cma=Cma
        self.Cmq=Cmq
        self.Cm_de=Cm_de
        self.Cm_df=Cm_df
        self.cd_w=cd_w
        self.cd_t=cd_t
        self.CD0=CD0
        self.TCprime=TCprime
        self.static=static_condition
        self.alpha_0=alpha_0
        self.h=h
        self.AR=AR
        self.b=np.sqrt(AR*S)
        self.c_bar=self.b/AR
        self.deltae0=deltae0
        self.deltaf0=deltaf0
        self.ic=initial_conditions

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


        #Forces
        L=Q*self.S*self.CL(t, x)


        # Thrust-drag
        if self.static:
            T_D=.5*self.rho*(self.omega*self.R)**2*np.pi*self.R**2*(self.TCprime-self.CD) #ask trevor how net cx changes
        else:
            T=15000*self.throttle(t) #need to change
            T_D=T-Q*self.S*self.CD(t, x)

        #Rolling resistance on ground
        W=self.m*9.81
        if ze<=0:
            mu_r=0.02 #coefficent of friction
            roll_resist=mu_r*max(W-L, 0) #eventually this will decrease to 0
        else:
            roll_resist=0


        # Forces in body frame
        X = T_D * np.cos(alpha) + L * np.sin(alpha)-roll_resist
        Z = T_D * np.sin(alpha) - L * np.cos(alpha)

        #Moment
        M=Q*self.S*self.Cm(t, x)*self.c_bar


        # Equations of motion
        u_dot = (1 / self.m) * X - 9.81 * np.sin(theta)
        w_dot = (1 / self.m) * Z + 9.81 * np.cos(theta) + u * q
        q_dot = 1 / self.I_y * M
        theta_dot = q
        x_dot = u_e
        z_dot = -w_e

        if ze<=0 and z_dot<0: #keeps solver from going to negative altitude
            z_dot=0

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
        ui, wi, qi, thetai, xei, zei, =self.ic
        Vi=np.sqrt(ui**2+wi**2)

        if Vi<1e-6:
            q_bari=0
        else:
            q_bari=(qi*self.c_bar)/(2*Vi)

        V = np.sqrt(u**2 + w**2)
        if V<1e-6:
            q_bar=0
        else:
            q_bar=(q*self.c_bar)/(2*V)

        del_deltae=self.elevator(t)-self.deltae0
        del_deltaf=self.flaps(t)-self.deltaf0
        del_alpha=alpha-self.alpha_0
        del_qbar=q_bar-q_bari

        return self.CL0+self.CLa*del_alpha+self.CLq*del_qbar+self.CL_df*del_deltaf+self.CL_de*del_deltae


    def Cm(self, t, x):
        u, w, q, theta, xe, ze = x
        alpha = np.arctan2(w, u)
        ui, wi, qi, thetai, xei, zei, =self.ic
        Vi=np.sqrt(ui**2+wi**2)
        q_bari=(qi*self.c_bar)/(2*Vi)

        if Vi<1e-6:
            q_bari=0
        else:
            q_bari=(qi*self.c_bar)/(2*Vi)

        V = np.sqrt(u**2 + w**2)
        if V<1e-6:
            q_bar=0
        else:
            q_bar=(q*self.c_bar)/(2*V)

        del_deltae=self.elevator(t)-self.deltae0
        del_deltaf=self.flaps(t)-self.deltaf0
        del_alpha=alpha-self.alpha_0
        del_qbar=q_bar-q_bari

        return self.Cm0+self.Cma*del_alpha+self.Cmq*del_qbar+self.Cm_df*del_deltaf+self.Cm_de*del_deltae
    def CD(self, t, x):
        u, w, q, theta, xe, ze = x
        if ze<=0: #include ground effect factor due to rolling resistance
            phi=(16*self.h/self.b)**2/(1+(16*self.h/self.b)**2)
        else:
            phi=1
        return self.CD0+self.cd_w+self.cd_t*self.St/self.S+phi*self.CL(t, x)**2/(np.pi*self.AR)

    """
    Modeling elevator surface delfection over time.
    Given a time input, elevator deflection will change, which will go into Cm, CL, etc caclulations"""

    def elevator(self, t):
        if t<3:
            return 0
        else:
            return np.radians(-10) #modellig pull up at takeoff after 5 sec
    def flaps(self, t):
        if t < 20:
            return np.radians(10)
        elif t < 25:
            return np.radians(10) * (1 - (t - 20) / 5)  # linearly decreasing funnction
        else:
            return 0
    def throttle(self, t):
        tau=2 #engine response time (sec)
        return 1-np.exp(-t/tau) #where does this go in thrust/drag calculation


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


#Create new versions of aircraft

#Area calculations from geom.py have to fix this or get a dict for the whole group
S = 49.5 #twin otter wing area is 39 m^2
AR = 8 #twin otter is 10.05
b = np.sqrt(AR*S)
MAC = b/AR
lv = 10 #distance from wing quarter chord to vertical tail quarter chord
Vv = .1 #Vertical tail volume coefficient
vt_ar = 1.2
tail_hinge = 0.7

Sv = Vv*S*b/lv
vt_c = np.sqrt(Sv/vt_ar) #constant chord
vt_b = vt_ar*vt_c
vt_x = lv + 1/4 * MAC - 1/4 * vt_c #distance from wing LE to tail (avg) LE

# horizontal tail geometry
lh = lv + vt_c/3 #include some sweep in the vertical tail to extend the ht moment arm
Vh = 1.05
ht_ar = 2

Sh = Vh*S*MAC/lh



plane_1=Aircraft(
    m_total=55622.7/9.81,
    I_y=40000, #UNKOWN
    rho=1.225,
    S=49.5,
    St=Sh,
    R=1.5, #UNKOWN
    omega=220, #propeller rotation rate, deltaT have to change this to a function of time, UNKOWN
    CL0=0.3,#UNKOWN
    CLa= 6.5, #U
    CLq=5, #U
    CL_de= 0.4, #U
    CL_df=1.2, #U
    # cx=0.3, #eventually cx either gets split into compoents or becomes a function  or a function of cl, dt, q_bar, df, Re, Ma
    # cx_prime=0.2,
    Cm0=0, #U
    Cma=-1.2,
    Cmq=-12,
    Cm_de=-1,
    Cm_df=-.1,
    cd_w=0,
    cd_t=0.01,
    CD0=0.013,
    TCprime=0.3, #U
    static_condition=False,
    alpha_0=np.radians(2),#U
    h=5.5,#U
    AR=8,
    deltae0=0, #U
    deltaf0=0, #U
    initial_conditions=x0)

#Intergate for 20 seconds
t_span = (0, 15)
t_eval = np.linspace(0, 15, 1000)

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
    t_to=sol.t_events[0][0]
    x_to=sol.y_events[0][0][4]
    v_to=np.sqrt(sol.y_events[0][0][0]**2+sol.y_events[0][0][1]**2)

    print(f"Takeoff time: {t_to:.3f} s")
    print(f"Takeoff distance: {x_to:.2f} m")
    print(f"Takeoff velocity: {v_to:.2f} m/s")

else:
    print("Aircraft did not lift off.")

# --- Plot results ---
# Positions

plt.plot(sol.y[4], sol.y[5])
plt.xlabel("x Position (m)")
plt.ylabel(" z Position (m)")
plt.title("Trajectory")
plt.grid(True)
plt.axis("equal")
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
