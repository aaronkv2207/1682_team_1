import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


#del_T=omega


def aircraft_dynamics(t, x, params):
    """
    x = [u, w, q, theta] state vector
    params = dictionary of aircraft parameters
    """
    u, w, q, theta, xe, ze = x
    m = params['m']
    I_y = params['I_y']
    rho = params['rho']
    S = params['S']
    R = params['R']
    omega = params['omega']
    cl=params['cl']
    cx = params['cx']
    cx_prime = params['cx_prime']
    M_alpha=params['M_alpha']
    M_q=params['M_q']
    static = params['static']
    h=params['h']
    b=params['b']
    AR=b**2/S

    # Kinematics
    alpha = np.arctan2(w, u)
    V = np.sqrt(u**2 + w**2)
    Q = 0.5 * rho * V**2

    # Transform velocities to earth frame
    u_e = u * np.cos(theta) + w * np.sin(theta)
    w_e = w * np.cos(theta) - u * np.sin(theta)
    gamma = np.arctan2(-w_e, u_e)

    #Forces
    # Thrust / drag
    if static:
        T_D = -Q * S * cx
    elif ze==0: #if at takeoff, consider phi component on induced drag
        phi=(16*h/b)**2/(1+(16*h/b)**2)
        T_D=-0.5 * rho * (omega * R**2) * (cx_prime)+phi*cl**2/(np.pi*AR)
    else:
        T_D = -0.5 * rho * (omega * R**2) * cx_prime

    #Lift
    L=Q*S*cl

    #Rolling resistance at takeoff
    if ze==0:
        W=m*9.81
        mu_r=0.02 #coefficent of friction
        R=mu_r*(W-L) #eventually this will decrease to 0

    # Forces in body frame
    X = T_D * np.cos(alpha) + L * np.sin(alpha)-R
    Z = T_D * np.sin(alpha) - L * np.cos(alpha)

    #Moment
    M=M_alpha*alpha+M_q*q


    # Equations of motion
    u_dot = (1/m) * X - 9.81 * np.sin(theta)
    w_dot = (1/m) * Z + 9.81 * np.cos(theta) + u * q
    q_dot = 1/I_y * M
    theta_dot = q
    x_dot = u_e
    z_dot = w_e

    return [u_dot, w_dot, q_dot, theta_dot, x_dot, z_dot]


# Initial state
# Example initial conditions
u0 = 0       # forward velocity in m/s
w0 = 0        # vertical velocity in m/s
q0 = 0        # pitch rate in rad/s
theta0 = 0    # pitch angle in rad
x_e0 = 0      # initial x-position in meters
z_e0 = 0  # initial altitude in meters
x0 = [u0, w0, q0, theta0, x_e0, z_e0]

# Parameters
params = {
    'm': 5670,
    'I_y': 500,
    'rho': 1.225,
    'S': 10,
    'R': 2,
    'omega': 100,
    'cl': 0.3,
    'cx': 0.3, #eventually cx either gets split into compoents or becomes a function  or a function of cl, dt, q_bar, df, Re, Ma
    'cx_prime': 0.2,
    'M_alpha': -2000, #restoring stiffness
    'M_q':-500, #pitch damping
    'static': False,
    'h':5.5,
    'b': 19.8



}
def takeoff_event(t, x):
    u, w =x[0], x[1]
    rho=params['rho']
    S=params['S']
    cl=params['cl']
    m=params['m']
    V=np.sqrt(u**2+w**2)
    Q=0.5*rho*V**2
    L=Q*cl*S
    W=m*9.81
    return L-W

takeoff_event.terminal=False # doesn't stop integrating with ivp once takeoff occurs but records when event happens
takeoff_event.direction=1 #ensures that crossing goes from L-W<0 to L-W>0 to find zero


# Integrate over 10 seconds
t_span = (0, 20)
t_eval = np.linspace(0, 10, 1000)

sol = solve_ivp(lambda t, x: aircraft_dynamics(t, x, params),t_span, x0, t_eval=t_eval, events=takeoff_event)

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
