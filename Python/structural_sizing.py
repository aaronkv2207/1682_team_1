# Team 1 - Structural Sizing
import numpy as np

# Assumptions/variables
L =
W_total =
W_wing =
W_cent = W_total - W_wing
a_z =
g = 9.81
w_b_max = ?? # maximum deflection to span ratio - not sure if this should be assumed or calculated
rho =
V =
q = 0.5*rho*V**2
twist_max =

# Sizing Parameters
b =
h =
c_tip =
c_0 =
b_ail =
c_ail =
c_m =
A = # shell enclosed area
y_ail =
s_tot = # not sure exactly what this is
R_outer =
l =

# Elasticity Values
stress_max = ???
shear_max = ??? # this comes up in a few different equations, not sure if it's the same or different for sizing different components
torsion_max = ???
bending_max = ???
E =
G =
I =
C = 0.25 # depends on how tube ends are supported (this is for fixed at one end only)
P_max = C * (np.pi**2*E*I)/l**2 # not sure what the difference is between P_cr and P_max in the pdf


# Load Factor
N = L/W_total
N_land = a_z/g
N_max = np.max(N, N_land)



# Spar sizing
taper_ratio = c_tip/c_0
A_cap_0_strength = (N_max*W_cent)/(12*stress_max)*(b/h)*(1+2*taper_ratio)/(1+taper_ratio)
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

t_tube_bending_strength = bending_max/(np.pi*R_outer**2*stress_max)

# NOTE: I think more complex equations need to be written out for k, F, and w depending on which we case we choose
t_tube_bending_deflection = (k*F_max*l**3)/(np.pi*R_outer**2*w_max*E)

t_tube_buckling = 1/C * (P_max*l**2)/(np.pi**3*R_outer**3*E)

tube_weight = rho*g*(torsion_max*l**2)/(R_outer**2*G) # not sure if we should do this calc or not

t_tube = np.max(t_tube_torsion_strength, t_tube_torsion_stiffness, t_tube_bending_strength, t_tube_bending_deflection, t_tube_buckling)
