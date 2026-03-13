# Wing Sizing

import numpy as np
import math as math
import matplotlib.pyplot as plt

# from aero_workspace.conceptual_design import MTOW, V_CRUISE, V_STALL, W_S, S, ureg
# from aero_workspace.aero_main import AircraftConfig as aircraft
# from aero_workspace.aero_main import TakeoffCoeff as takeoff


# TODO: go through and check for b/2, L/2, etc. stuff (depending on how we define b, L, etc. in the code)
# check signs in stress equations

# using aero_main dataclass here with just the values they gave us
from dataclasses import dataclass

import numpy as np
import pandas as pd
# from ambiance import Atmosphere
from aero_workspace.conceptual_design import MTOW, V_CRUISE, V_STALL, W_S, S, ureg
from scipy.interpolate import interp1d

# get jvl cg --> mset, run file, input file, ...; can shift around masses to shift cg in .mass file


# TODO: update constants in dataclass
@dataclass
class aircraft:
    """All aircraft definitions"""

    #### GLOBAL DEFINITIONS ###
    AR: int = 8
    s_ref: float = S.magnitude   # convert from pint quantity

    #### TAKEOFF DEFINITIONS ###
    v_takeoff: float = (1.2 * V_STALL).magnitude

    #### CRUISE DEFINITIONS ###
    v_cruise: float = V_CRUISE.magnitude
    h_cruise: float = (18000 * ureg("ft")).to("m").magnitude

    #### LANDING DEFINITIONS ###
    v_landing: float = (1.2 * V_STALL).magnitude


class Wing():
    """Wing class:
    - Takes in general wing sizing and material
    - Sizes structural components on the wing
    - Outputs weight of the wing"""

    # Set constants
    g = 9.81
    rho_ground = 1.225
    rho_cruise = 0.699 # at 18,000 ft
    W_total = 7438.9149 * 9.8
    W_duct = 25 * 9.81 # kg - guess
    W_fuel_per_wing = 750 # kg

    def __init__(self, aero, loading, materials, weight_estimate = 0.07*W_total):
        self.aero = aero
        self.loading = loading
        self.materials = materials
        self.weight_estimate = weight_estimate

    # NOTE: commented out functions do better estimations for stresses uses distributed loads
    # but give pretty seem values to the functions currently being used

    def trapz_manual(self, y_vals, x_vals):
        """Manual trapezoidal integration: y_vals over x_vals"""
        integral = 0
        for i in range(len(y_vals) - 1):
            dx = x_vals[i+1] - x_vals[i]
            integral += 0.5 * (y_vals[i] + y_vals[i+1]) * dx
        return integral

    def axial_stress(self, L, T, D):
        # Span and chord
        b = np.sqrt(aircraft.AR*aircraft.s_ref)/2
        c = np.sqrt(aircraft.s_ref/aircraft.AR)
        x_1 = 0.1*c
        x_2 = 0.7*c

        t_x1 = self.aero["t_x1"]
        t_x2 = self.aero["t_x2"]
        h = 0.5*(t_x1 + t_x2)
        w = x_2 - x_1
        t = t_x1 * 0.1  # box thickness estimate

        # Span discretization
        n_span = 50
        y_span = np.linspace(0, b, n_span)

        # Elliptical distributions along half-span
        lift_dist = L * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)
        drag_dist = D * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)
        weight_dist = self.weight_estimate * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)

        # Add point masses for ducted fans
        thrust_dist = np.zeros_like(y_span)
        fan_positions = [b/5, b*2/5, b*3/5, b*4/5]
        thrust_per_fan = T/4
        for pos in fan_positions:
            idx = (np.abs(y_span - pos)).argmin()
            thrust_dist[idx] += thrust_per_fan

        # Bending moments
        M_y = self.trapz_manual(drag_dist * y_span, y_span)  # about y-axis
        M_z = self.trapz_manual(lift_dist * y_span, y_span) + self.trapz_manual(weight_dist * y_span, y_span) - self.trapz_manual(thrust_dist * y_span, y_span)

        # Moments of inertia
        I_y = (h*w**3)/12 - ((h-2*t)*(w-2*t)**3)/12
        I_z = (w*h**3)/12 - ((w-2*t)*(h-2*t)**3)/12

        # Section distances from neutral axis
        z = w/2
        y = h/2

        # Axial stresses
        axial_yy = -(M_z * y)/I_z
        axial_zz = (M_y * z)/I_y
        axial_max = max(abs(axial_yy), abs(axial_zz))

        return 1.5*axial_max

    def shear_stress(self, L, T, D):
        # Span and chord
        b = np.sqrt(aircraft.AR*aircraft.s_ref)/2
        c = np.sqrt(aircraft.s_ref/aircraft.AR)
        x_1 = 0.1*c
        x_2 = 0.7*c

        t_x1 = self.aero["t_x1"]
        t_x2 = self.aero["t_x2"]
        h = 0.5*(t_x1 + t_x2)
        w = x_2 - x_1
        t = t_x1 * 0.1  # box thickness estimate

        # Span discretization
        n_span = 50
        y_span = np.linspace(0, b, n_span)

        # Elliptical distributions along half-span
        lift_dist = L * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)
        drag_dist = D * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)
        weight_dist = self.weight_estimate * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)

        # Add point thrusts
        thrust_dist = np.zeros_like(y_span)
        fan_positions = [b/5, b*2/5, b*3/5, b*4/5]
        thrust_per_fan = T/4
        for pos in fan_positions:
            idx = (np.abs(y_span - pos)).argmin()
            thrust_dist[idx] += thrust_per_fan

        # Shear forces
        F_y = self.trapz_manual(weight_dist, y_span) + self.trapz_manual(thrust_dist, y_span) - self.trapz_manual(lift_dist, y_span)
        F_z = self.trapz_manual(drag_dist, y_span)

        # First moments of area
        Q_y = (h*w**2)/8 - ((h-2*t)*(w-2*t)**2)/8
        Q_z = (w*h**2)/8 - ((w-2*t)*(h-2*t)**2)/8

        # Moments of inertia
        I_y = (h*w**3)/12 - ((h-2*t)*(w-2*t)**3)/12
        I_z = (w*h**3)/12 - ((w-2*t)*(h-2*t)**3)/12

        shear_xy = (F_y * Q_z)/(I_z * t)
        shear_xz = (F_z * Q_y)/(I_y * t)
        shear_max = max(abs(shear_xy), abs(shear_xz))

        return 1.5*shear_max


    def get_N(self, L, a_z):
        """Gets flight load factor, landing load factor, and max load factor"""
        # Get required variable
        # W_total = self.aero["W_total"]

        # Calculate max load factor
        N = L/self.W_total # flight load factor
        N_land = a_z/self.g # landing load factor
        N_max = max(N, N_land) # maximum load factor

        return N, N_land, N_max

    def spar_cap_area(self, L, a_z, axial_stress):
        """Find necessary spar cap area based on loading"""
        # Get required variables
        c = np.sqrt(aircraft.s_ref/aircraft.AR)
        c_tip = c
        c_0 = c
        taper_ratio = c_tip/c_0
        N_max = self.get_N(L, a_z)[2]
        W_cent = self.W_total - self.weight_estimate
        b = np.sqrt(aircraft.AR*aircraft.s_ref)/2
        t_x1 = self.aero["t_x1"]
        t_x2 = self.aero["t_x2"]
        x_1 = 0.1*c
        x_2 = 0.7*c
        h = 0.5*(t_x1+t_x2)
        N = self.get_N(L, a_z)[0]
        E = self.materials["spar_cap_E"] # maybe remain this variable if it shows up again just for clarity
        w_b_max = self.aero["w_b_max"] # not sure where we'll actually get this

        # Calculate area based on strength and stiffness requirements
        A_cap_0_strength = (N_max*W_cent)/(12*axial_stress)*(b/h)*(1+2*taper_ratio)/(1+taper_ratio)
        A_cap_0_stiffness = (N*W_cent)/(48*E) * b**2/h**2 * 1/w_b_max * (1+2*taper_ratio)/(1+taper_ratio)
        return max(A_cap_0_strength, A_cap_0_stiffness)

    def spar_web_area(self, L, a_z, shear_stress):
        """Find necessary spar web area based on loading"""
        # Get required variables
        N_max = self.get_N(L, a_z)[2]
        W_cent = self.W_total - self.weight_estimate

        # Return spar web area
        return (N_max*W_cent)/(2*shear_stress)

    def skin_thickness(self, q, shear_stress):
        """Find necessary skin thickness based on loading"""
        # Get required variables
        b = np.sqrt(aircraft.AR*aircraft.s_ref)/2
        c = np.sqrt(aircraft.s_ref/aircraft.AR)

        # NOTE: all of these are guesses
        b_ail = 0.25*b
        c_ail = 0.75*c
        y_ail = 0.65*b
        c_m = -0.05


        # x1 = self.aero["x1"]
        # x2= self.aero["x2"]
        # tx1 = self.aero["tx1"]
        # tx2 = self.aero["tx2"]
        t_avg = self.aero["t_avg"] # Average thickness to chord of the airfoil
        s_tot = None
        # s_tot = self.aero["airfoil_perimeter"] # we'll estimate this as 2 * c_ail for now
        G = self.materials["skin_G"]
        twist_max = self.aero["twist_max"]
        A = None
        # A = self.aero["airfoil_cross_section_area"] # we'll estimate this as b_ail * c_ail for now
        T = abs(c_m) * 0.5*self.rho_cruise*aircraft.v_cruise**2 * aircraft.s_ref * c
        print("T:", T)

        if A is None:
            A = t_avg * c_ail

        if s_tot is None:
            s_tot = 2 * c_ail

        # Calculate torsional stiffness J per thickness assuming the wing skin is a box beam with width b_ail and height c_ail, and perimeter 2*c_ail
        J_over_thickness = 4 * A ** 2 / s_tot
        # Calculate skin strength requirement
        # t_skin_strength = 0 # making this 0 for rn because it is causing issues
        t_skin_strength = q*b_ail*c_ail**2*abs(c_m)*1/(2*A)*(1/shear_stress)
        # Calculate skin stiffness requirement
        t_skin_stiffness = T/G * y_ail / J_over_thickness * 1 / twist_max

        max_skin_thickness = max(t_skin_strength, t_skin_stiffness)
        if max(t_skin_strength, t_skin_stiffness) < 0.0035:
            max_skin_thickness = 0.0035
        # Return max thickness
        return max_skin_thickness

    def tube_thickness(self):
        # TODO: implement this (not sure where this component goes so not doing this yet)
        pass

    def takeoff(self):
        # Find forces - do we need to incorporate some takeoff angle into this? maybe for all the stuff angle is an option and during cruise it's just zero?
        q_takeoff = 0.5*self.rho_ground*aircraft.v_takeoff**2
        L_takeoff = self.aero["CL_takeoff"]*q_takeoff*aircraft.s_ref # or lift distribution from aero
        T_takeoff = self.loading["T_takeoff"] # not sure where this will come from exactly
        D_takeoff = self.aero["CD_takeoff"]*q_takeoff*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_takeoff = self.axial_stress(L_takeoff, T_takeoff, D_takeoff)
        shear_stress_takeoff = self.shear_stress(L_takeoff, T_takeoff, D_takeoff)
        print("axial stress takeoff:", axial_stress_takeoff)
        print("shear stress takeoff:", shear_stress_takeoff)

        # Find component sizing based on calculated loading
        # NOTE: setting a_z = 0 on all cases but landing so that N_land is not considered
        takeoff_spar_cap_area = self.spar_cap_area(L_takeoff, 0, axial_stress_takeoff)
        takeoff_spar_web_area = self.spar_web_area(L_takeoff, 0, shear_stress_takeoff)
        takeoff_skin_thickness = self.skin_thickness(q_takeoff, shear_stress_takeoff)
        # takeoff_tube_thickness = self.tube_thickness() # still not sure what this is for

        return takeoff_spar_cap_area, takeoff_spar_web_area, takeoff_skin_thickness #, takeoff_tube_thickness


    def climb(self):
        # Find forces - do we need to incorporate climb angle into this? maybe for all the stuff angle is an option and during cruise it's just zero?
        q_climb = 0.5*(0.5*(self.rho_ground+self.rho_cruise))*self.aero["v_climb"]**2
        L_climb = self.aero["CL_climb"]*q_climb*aircraft.s_ref # or lift distribution from aero
        T_climb = self.loading["T_climb"] # not sure where this will come from exactly
        D_climb = self.aero["CD_climb"]*q_climb*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_climb = self.axial_stress(L_climb, T_climb, D_climb)
        shear_stress_climb = self.shear_stress(L_climb, T_climb, D_climb)
        print("axial stress climb:", axial_stress_climb)
        print("shear stress climb:", shear_stress_climb)

        # Find component sizing based on calculated loading
        # NOTE: setting a_z = 0 on all cases but landing so that N_land is not considered
        climb_spar_cap_area = self.spar_cap_area(L_climb, 0, axial_stress_climb)
        climb_spar_web_area = self.spar_web_area(L_climb, 0, shear_stress_climb)
        climb_skin_thickness = self.skin_thickness(q_climb, shear_stress_climb)
         # climb_tube_thickness = self.tube_thickness()

        return climb_spar_cap_area, climb_spar_web_area, climb_skin_thickness #, climb_tube_thickness

    def cruise(self):
        # Find forces
        q_cruise = 0.5*self.rho_cruise*aircraft.v_cruise**2
        L_cruise = self.aero["CL_cruise"]*q_cruise*aircraft.s_ref # or lift distribution from aero
        T_cruise = self.loading["T_cruise"] # not sure where this will come from exactly
        D_cruise = self.aero["CD_cruise"]*q_cruise*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_cruise = self.axial_stress(L_cruise, T_cruise, D_cruise)
        shear_stress_cruise = self.shear_stress(L_cruise, T_cruise, D_cruise)
        print("axial stress cruise:", axial_stress_cruise)
        print("shear stress cruise:", shear_stress_cruise)

        # Find component sizing based on calculated loading
        # NOTE: setting a_z = 0 on all cases but landing so that N_land is not considered
        cruise_spar_cap_area = self.spar_cap_area(L_cruise, 0, axial_stress_cruise)
        cruise_spar_web_area = self.spar_web_area(L_cruise, 0, shear_stress_cruise)
        cruise_skin_thickness = self.skin_thickness(q_cruise, shear_stress_cruise)
         # cruise_tube_thickness = self.tube_thickness()

        return cruise_spar_cap_area, cruise_spar_web_area, cruise_skin_thickness # , cruise_tube_thickness

    def descent(self):
        # Find forces - do we need to incorporate descent angle into this? maybe for all the stuff angle is an option and during cruise it's just zero?
        q_descent = 0.5*(0.5*(self.rho_ground+self.rho_cruise))*self.aero["v_descent"]**2
        L_descent = self.aero["CL_descent"]*q_descent*aircraft.s_ref # or lift distribution from aero
        T_descent = self.loading["T_descent"] # not sure where this will come from exactly
        D_descent = self.aero["CD_descent"]*q_descent*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_descent = self.axial_stress(L_descent, T_descent, D_descent)
        shear_stress_descent = self.shear_stress(L_descent, T_descent, D_descent)
        print("axial stress descent:", axial_stress_descent)
        print("shear stress descent:", shear_stress_descent)

        # Find component sizing based on calculated loading
        # NOTE: setting a_z = 0 on all cases but landing so that N_land is not considered
        descent_spar_cap_area = self.spar_cap_area(L_descent, 0, axial_stress_descent)
        descent_spar_web_area = self.spar_web_area(L_descent, 0, shear_stress_descent)
        descent_skin_thickness = self.skin_thickness(q_descent, shear_stress_descent)
        # descent_tube_thickness = self.tube_thickness()

        return descent_spar_cap_area, descent_spar_web_area, descent_skin_thickness # , descent_tube_thickness

    def landing(self):
        # Find forces - do we need to incorporate landing angle into this? maybe for all the stuff angle is an option and during cruise it's just zero?
        q_landing = 0.5*self.rho_ground*aircraft.v_landing**2
        L_landing = self.aero["CL_landing"]*q_landing*aircraft.s_ref # or lift distribution from aero
        T_landing = self.loading["T_landing"] # not sure where this will come from exactly
        D_landing = self.aero["CD_landing"]*q_landing*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_landing = self.axial_stress(L_landing, T_landing, D_landing)
        shear_stress_landing = self.shear_stress(L_landing, T_landing, D_landing)
        print("axial stress landing:", axial_stress_landing)
        print("shear stress landing:", shear_stress_landing)

        # Find component sizing based on calculated loading
        # NOTE: setting a_z = 0 on all cases but landing so that N_land is not considered
        landing_spar_cap_area = self.spar_cap_area(L_landing, self.aero["a_z"], axial_stress_landing)
        landing_spar_web_area = self.spar_web_area(L_landing, self.aero["a_z"], shear_stress_landing)
        landing_skin_thickness = self.skin_thickness(q_landing, shear_stress_landing)

        return landing_spar_cap_area, landing_spar_web_area, landing_skin_thickness

    def max_load_sizing(self):
        sizing_takeoff = self.takeoff()
        sizing_climb = self.climb()
        sizing_cruise = self.cruise()
        sizing_descent = self.descent()
        # sizing_descent = 0,0,0
        sizing_landing = self.landing()

        print("sizing_takeoff:", sizing_takeoff)
        print("sizing_climb:", sizing_climb)
        print("sizing_cruise:", sizing_cruise)
        print("sizing_descent:", sizing_descent)
        print("sizing_landing:", sizing_landing)

        spar_cap_area = max(sizing_takeoff[0], sizing_climb[0], sizing_cruise[0], sizing_descent[0], sizing_landing[0])
        spar_web_area = max(sizing_takeoff[1], sizing_climb[1], sizing_cruise[1], sizing_descent[1], sizing_landing[1])
        skin_thickness = max(sizing_takeoff[2], sizing_climb[2], sizing_cruise[2], sizing_descent[2], sizing_landing[2])

        return spar_cap_area, spar_web_area, skin_thickness


    def wing_weight(self):
        b = np.sqrt(aircraft.AR*aircraft.s_ref)/2
        # for one of two wings
        sizing = self.max_load_sizing()
        print(sizing)
        spar_cap_weight = sizing[0]*b*self.materials["spar_cap_density"]
        spar_web_weight = sizing[1]*b*self.materials["spar_web_density"]
        skin_weight = sizing[2]*self.aero["airfoil_surface_area"]*self.materials["skin_density"]

        print("spar_cap_weight", spar_cap_weight)
        print("spar_web_weight", spar_web_weight)
        print("skin_weight", skin_weight)

        # should this weight just be structural stuff or also fuel and motors?
        return 2*spar_cap_weight + 2*spar_web_weight + skin_weight


if __name__ == "__main__":
# Create dictionaries for testing
    c = np.sqrt(aircraft.s_ref/aircraft.AR)
    aero = {
        "t_x1": (0.7279684E-01 + 0.2596068E-01)*c,
        "t_x2": (0.2139701E-01 + 0.4671981E-01)*c,
        # all below values are guesses
        "w_b_max": 0.025,
        "twist_max": 0.0524, # radians
        "CL_takeoff": 6,
        "CL_climb": 1,
        "CL_cruise": 0.35,
        "CL_descent": 0.2,
        "CL_landing": 6,
        "CD_takeoff": 0.7,
        "CD_climb": 0.1,
        "CD_cruise": 0.05,
        "CD_descent": 0.045,
        "CD_landing": 1,
        "v_takeoff": 25,
        "v_climb": 60,
        "v_cruise": 125,
        "v_descent": 75,
        "v_landing": 25,
        "drag_area": 2.5,
        "a_z": 3*9.8,
        "airfoil_surface_area": 100, # m^2
        "t_avg": 0.12,
    }

    loading = {
        "T_takeoff": 27 * 10**3,
        "T_climb": 18 * 10**3, # guess
        "T_cruise": 8.5 * 10**3,
        "T_descent": 5 * 10**3, # guess
        "T_landing": 12 * 10**3, # guess
    }

    materials = {
        "spar_cap_E": 6.9 * 10**9,
        "skin_G": 26 * 10**9,
        "spar_cap_density": 2700, #kg/m^3
        # saying spar web is also aluminum for initial testing
        "spar_web_density": 2700, #kg/m^3
        "skin_density": 2700 #kg/m^3
    }

    # wing weight estimate
    # weight_estimate = 0.07*aero["W_total"]

    test_wing = Wing(aero, loading, materials)

    test_wing_weight = test_wing.wing_weight()
    print("Wing weight (kg)", test_wing_weight/9.81)
