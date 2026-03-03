# Team 1 - Structural Sizing
import numpy as np

# Assumptions/variables
L = ... # lift
W_total = ... # total weight
W_wing = ... # wing weight
W_cent = W_total - W_wing # center weight
a_z = ... # acceleration at landing
g = 9.81
w_b_max = ... # maximum deflection to span ratio - not sure if this should be assumed or calculated
rho = ... # density (double check if we need material vs. air density in different places)
V = ... # airspeed
q = 0.5*rho*V**2 # dynamic pressure
twist_max = ... # maximum tolerable tip twist angle

# Sizing Parameters
b = ... # span
h = ... # spar height
c_tip = ... # chord at the tip
c_0 = ... # chord at the root
b_ail = ... # span from root to aileron?
c_ail = ... # chord at the aileron
c_m = ... # local pitching moment coefficient
A = ... # shell enclosed area
y_ail = ... # spanwise coordinate of aileron?
s_tot = ... # not sure exactly what this is
R_outer = ... # tube outer radius
l = ... # length

# Elasticity Values
axial_max = ... # maximum axial stress
shear_max = ... # maximum shear stress
torsion_max = ... # maximum torsional moment
bending_max = ... # maximum bending moment
E = ... # stiffness modulus
G = ... # shear modulus
I = ... # bending moment of inertia
C = 0.25 # depends on how tube ends are supported (this is for fixed at one end only)
P_max = C * (np.pi**2*E*I)/l**2 # not sure what the difference is between P_cr and P_max in the pdf


# Load Factor
N = L/W_total # flight load factor
N_land = a_z/g # landing load factor
N_max = np.max(N, N_land) # maximum load factor



# Spar sizing
taper_ratio = c_tip/c_0
A_cap_0_strength = (N_max*W_cent)/(12*axial_max)*(b/h)*(1+2*taper_ratio)/(1+taper_ratio)
A_cap_0_stiffness = (N*W_cent)/(48*E) * b**2/h**2 * 1/w_b_max * (1+2*taper_ratio)/(1+taper_ratio)
A_cap_0 = np.max(A_cap_0_strength, A_cap_0_stiffness)

A_web_0 = (N_max*W_cent)/(2*shear_max)

# Skin Sizing
t_skin_strength = q*b_ail*c_ail**2*c_m*1/(2*A)*1/(1/shear_max)
t_skin_stiffness = 1*b_ail*c_ail**2*y_ail*c_m*s_tot/(4*A**2*G)*(1/twist_max)

t_skin = np.max(t_skin_strength, t_skin_stiffness)

# Tube Sizing
t_tube_torsion_strength = torsion_max/(2*np.pi*R_outer**2*shear_max)
t_tube_torsion_stiffness = (torsion_max*l)/(2*np.pi*R_outer**3*G)

t_tube_bending_strength = bending_max/(np.pi*R_outer**2*axial_max)

# NOTE: I think more complex equations need to be written out for k, F, and w depending on which we case we choose
t_tube_bending_deflection = (k*F_max*l**3)/(np.pi*R_outer**2*w_max*E)

t_tube_buckling = 1/C * (P_max*l**2)/(np.pi**3*R_outer**3*E)

tube_weight = rho*g*(torsion_max*l**2)/(R_outer**2*G) # not sure if we should do this calc or not

t_tube = np.max(t_tube_torsion_strength, t_tube_torsion_stiffness, t_tube_bending_strength, t_tube_bending_deflection, t_tube_buckling)
