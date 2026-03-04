import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class Aircraft:
    def __init__(self, m_total, I_y, rho, S, R, omega, cl, cx, cx_prime, cm, static_condition, h, b, initial_conditions):
        self.m=m_total
        self.I_y=I_y
        self.S=S
        self.R=R
        self.rho=rho
        self.omega=omega
        self.cl=cl
        self.cx=cx
        self.cx_prime=cx_prime
        self.cm=cm
        self.static=static_condition
        self.h=h
        self.b=b
        self.AR=b**2/S
        self.ic=initial_conditions

    def aircraft_dynamics(self, t, x):
        """
        x = [u, w, q, theta] state vector
        params = dictionary of aircraft parameters
        """
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
        # Thrust / drag
        if self.static:
            T_D = -Q * self.S * self.cx
        elif ze==0: #if at takeoff, consider phi component on induced drag
            phi=(16*self.h/self.b)**2/(1+(16*self.h/self.b)**2)
            T_D=-0.5 * self.rho * (self.omega * self.R**2) * (self.cx_prime)+phi*self.cl**2/(np.pi*self.AR)
        else:
            T_D = -0.5 * self.rho * (self.omega * self.R**2) * self.cx_prime

        #Lift
        L=Q*self.S*self.cl

        #Rolling resistance at takeoff
        if ze==0:
            W=self.m*9.81
            mu_r=0.02 #coefficent of friction
            R=mu_r*(W-L) #eventually this will decrease to 0

        # Forces in body frame
        if ze==0:
            X = T_D * np.cos(alpha) + L * np.sin(alpha)-R
        else:
            X = T_D * np.cos(alpha) + L * np.sin(alpha)
        Z = T_D * np.sin(alpha) - L * np.cos(alpha)

        #Moment
        M=Q*self.S*(self.cm-.5*q) #added damping term, will go away in real model


        # Equations of motion
        u_dot = (1/self.m) * X - 9.81 * np.sin(theta)
        w_dot = (1/self.m) * Z + 9.81 * np.cos(theta) + u * q
        q_dot = 1/self.I_y * M
        theta_dot = q
        x_dot = u_e
        z_dot = w_e

        return [u_dot, w_dot, q_dot, theta_dot, x_dot, z_dot]
    def takeoff_event(self, t, x):
        u, w =x[0], x[1]
        V=np.sqrt(u**2+w**2)
        Q=0.5*self.rho*V**2
        L=Q*self.cl*self.S
        W=self.m*9.81
        return L-W

Aircraft.takeoff_event.terminal=False # doesn't stop integrating with ivp once takeoff occurs but records when event happens
Aircraft.takeoff_event.direction=1 #ensures that crossing goes from L-W<0 to L-W>0 to find zero

# Initial state
# Example initial conditions
u0 = 40       # forward velocity in m/s
w0 = 0        # vertical velocity in m/s
q0 = 0        # pitch rate in rad/s
theta0 = 0.05    # pitch angle in rad
x_e0 = 0      # initial x-position in meters
z_e0 = 0  # initial altitude in meters
x0 = [u0, w0, q0, theta0, x_e0, z_e0]


#Create new versions of aircraft

plane_1=Aircraft(
    m_total=5670,
    I_y=500,
    rho=1.225,
    S=10,
    R=2,
    omega=100, #propeller rotation rate
    cl=0.3,
    cx=0.3, #eventually cx either gets split into compoents or becomes a function  or a function of cl, dt, q_bar, df, Re, Ma
    cx_prime=0.2,
    cm=0,
    static_condition=False,
    h=5.5,
    b=19.8,
    initial_conditions=x0)

#Intergate for 20 seconds
t_span = (0, 20)
t_eval = np.linspace(0, 20, 1000)

sol = solve_ivp(plane_1.aircraft_dynamics, t_span, plane_1.ic, t_eval=t_eval, events=plane_1.takeoff_event)

# sol.y[i] gives time series of ith state

#RESULTS
#x takeoff and at what time
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
plt.xlabel('x Position (m)')
plt.ylabel(' z Position (m)')
plt.title('Trajectory')
plt.grid(True)
plt.show()

# Velocities
u=sol.y[0]
w=sol.y[1]
v=np.sqrt(u**2+w**2)

plt.plot(sol.t, v)
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Airspeed vs Time')
plt.grid(True)
plt.show()
