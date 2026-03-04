import numpy as np
import matplotlib.pyplot as plt

# Constants
V_stall = 20 # m/s
rho = 1.225
mu = 0.02 # Rolling friction coefficient
e = 0.7 # Oswald efficiency factor
S = 20 
AR = 20
W = 84516.21
m = 8618.255 # 19,000 lbs in kg
g = 9.81 # gravitational accelearation [m/s^2]


x = 91.44 # 300 ft in m


CL = 0.2
CD0 = 0.1
CDi = CL**2/(np.pi * AR * e)
CD = CD0 + CDi




def D(v):
    return 1/2 * rho * v**2 * S * CDi
def L(v):
    return 1/2 * rho * v**2 * S * CL
def T(v):
    return D(v) + mu*(W - L(v)) + m*v/(2*x)

def T_c(v, R):
    return T(v) / (1/2*rho*(v**2) * np.pi * R**2)
def eta_ideal(v, R):
    return 2 / (1+ np.sqrt(1+T_c(v)))
def P_shaft_required(v,R):
    return (T(v) * v) / eta_ideal(v)






    