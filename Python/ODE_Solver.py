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
    cx = params['cx']
    cx_prime = params['cx_prime']
    L = params['L']  # could be a function of alpha, q, etc.
    static = params['static']

    # Kinematics
    alpha = np.arctan2(w, u)
    V = np.sqrt(u**2 + w**2)
    Q = 0.5 * rho * V**2

    # Transform velocities to earth frame
    u_e = u * np.cos(theta) + w * np.sin(theta)
    w_e = w * np.cos(theta) - u * np.sin(theta)
    gamma = np.arctan2(-w_e, u_e)

    # Thrust / drag
    if static:
        T_D = -Q * S * cx
    else:
        T_D = -0.5 * rho * (omega * R**2) * cx_prime


    # Forces in body frame
    X = T_D * np.cos(alpha) + L * np.sin(alpha)
    Z = T_D * np.sin(alpha) - L * np.cos(alpha)

    # Equations of motion
    u_dot = (1/m) * X - 9.81 * np.sin(theta)
    w_dot = (1/m) * Z + 9.81 * np.cos(theta) + u * q
    q_dot = 1/I_y * params['M']  # M should be moment
    theta_dot = q
    x_dot = u_e
    z_dot = w_e

    return [u_dot, w_dot, q_dot, theta_dot, x_dot, z_dot]


# Initial state
# Example initial conditions
u0 = 50       # forward velocity in m/s
w0 = 0        # vertical velocity in m/s
q0 = 0        # pitch rate in rad/s
theta0 = 0    # pitch angle in rad
x_e0 = 0      # initial x-position in meters
z_e0 = 1000   # initial altitude in meters
x0 = [u0, w0, q0, theta0, x_e0, z_e0]

# Parameters
params = {
    'm': 1000,
    'I_y': 500,
    'rho': 1.225,
    'S': 10,
    'R': 2,
    'omega': 100,
    'cx': 0.3,
    'cx_prime': 0.2,
    'L': 500,      # or a function of alpha
    'M': 50,       # pitching moment
    'static': True
}

# Integrate over 10 seconds
t_span = (0, 10)
t_eval = np.linspace(0, 10, 1000)

sol = solve_ivp(lambda t, x: aircraft_dynamics(t, x, params),
                t_span, x0, t_eval=t_eval)

print(sol.y)  # sol.y[i] gives time series of ith state

# --- Plot results ---
plt.figure(figsize=(12,6))

# Positions
plt.subplot(2,1,1)
plt.plot(sol.t, sol.y[4], label='x position (m)')
plt.plot(sol.t, sol.y[5], label='z position (m)')
plt.xlabel('Time (s)')
plt.ylabel('Position (m)')
plt.title('Aircraft Position vs Time')
plt.legend()
plt.grid(True)

# Velocities
plt.subplot(2,1,2)
plt.plot(sol.t, sol.y[0], label='u (m/s)')
plt.plot(sol.t, sol.y[1], label='w (m/s)')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.title('Aircraft Velocity vs Time')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
