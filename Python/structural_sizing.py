# Team 1 - Structural Sizing
import numpy as np
from ... import ... as beam

class Wing():
    """Wing class:
    - Takes in general wing sizing and material
    - Sizes structural components on the wing
    - Outputs weight of the wing"""

    # Set constants
    g = 9.81

    def __init__(self, aero, loading, materials, weight_estimate):
        self.aero = aero
        self.loading = loading
        self.material = materials
        self.weight_estimate = weight_estimate

    def get_N(self, L, a_z):
        """Gets flight load factor, landing load factor, and max load factor"""
        # Get required variable
        W_total = self.aero["W_total"]

        # Calculate max load factor
        N = L/W_total # flight load factor
        N_land = a_z/g # landing load factor
        N_max = np.max(N, N_land) # maximum load factor

        return N, N_land, N_max

    def spar_cap_area(self):
        """Find necessary spar cap area based on loading"""
        # Get required variables
        c_tip = self.aero["c_tip"]
        c_0 = self.aero["c_0"]
        taper_ratio = c_tip/c_0
        N_max = self.get_N()[2]
        W_cent = self.aero["W_total"] - self.weight_estimate
        axial_max = ... # TODO: fill in function from other script
        b = self.aero["b"]
        h = self.aero["h"]
        N = self.get_N()[0]
        E = self.materials["spar_car_E"] # maybe remain this variable if it shows up again just for clarity
        w_b_max = self.aero["w_b_max"] # not sure where we'll actually get this

        # Calculate area based on strength and stiffness requirements
        A_cap_0_strength = (N_max*W_cent)/(12*axial_max)*(b/h)*(1+2*taper_ratio)/(1+taper_ratio)
        A_cap_0_stiffness = (N*W_cent)/(48*E) * b**2/h**2 * 1/w_b_max * (1+2*taper_ratio)/(1+taper_ratio)
        return np.max(A_cap_0_strength, A_cap_0_stiffness)

    def spar_web_area(self):
        """Find necessary spar web area based on loading"""
        # Get required variables
        N_max = self.get_N()[2]
        W_cent = self.aero["W_total"] - self.weight_estimate
        shear_max = ... # TODO: fill in function from other script

        # Return spar web area
        return (N_max*W_cent)/(2*shear_max)

    def skin_thickness(self):
        """Find necessary skin thickness based on loading"""


        t_skin_strength = q*b_ail*c_ail**2*c_m*1/(2*A)*1/(1/shear_max)
        t_skin_stiffness = 1*b_ail*c_ail**2*y_ail*c_m*s_tot/(4*A**2*G)*(1/twist_max)

        return np.max(t_skin_strength, t_skin_stiffness)

    def tube_thickness():
        # TODO: implement this (not sure where this component goes so not doing this yet)
        pass

    # add methods for each flight stage that gets the correct loads and then sizes based off that
    # also add method to run through each flight stage and then get max loading
    # weight method takes max loading and uses that sizing to get weight of components


class Tail():
    pass

class Fuselage():
    pass

class LandingGear():
    pass


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
