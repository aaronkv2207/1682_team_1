# Team 1 - Structural Sizing
import sys
import os

# base_dir = os.path.dirname(__file__)
# sys.path.append(base_dir)
# sys.path.append(os.path.join(base_dir, "aero_workspace"))

import numpy as np
import math as math

# from aero_workspace.conceptual_design import MTOW, V_CRUISE, V_STALL, W_S, S, ureg
# from aero_workspace.aero_main import AircraftConfig as aircraft
# from aero_workspace.aero_main import TakeoffCoeff as takeoff


# TODO: go through and check for b/2, L/2, etc. stuff (depending on how we define b, L, etc. in the code)
# check signs in stress equations

# using aero_main dataclass here with just the values they gave us
from dataclasses import dataclass

import numpy as np
import pandas as pd
from ambiance import Atmosphere
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

    # def trapz_manual(self, y_vals, x_vals):
    #     """Manual trapezoidal integration: y_vals over x_vals"""
    #     integral = 0
    #     for i in range(len(y_vals) - 1):
    #         dx = x_vals[i+1] - x_vals[i]
    #         integral += 0.5 * (y_vals[i] + y_vals[i+1]) * dx
    #     return integral

    # def axial_stress(self, L, T, D):
    #     # Span and chord
    #     b = np.sqrt(aircraft.AR*aircraft.s_ref)/2
    #     c = np.sqrt(aircraft.s_ref/aircraft.AR)
    #     x_1 = 0.1*c
    #     x_2 = 0.7*c

    #     t_x1 = self.aero["t_x1"]
    #     t_x2 = self.aero["t_x2"]
    #     h = 0.5*(t_x1 + t_x2)
    #     w = x_2 - x_1
    #     t = t_x1 * 0.1  # box thickness estimate
    #     print("thickness of box:", t)

    #     # Span discretization
    #     n_span = 50
    #     y_span = np.linspace(0, b, n_span)

    #     # Elliptical distributions along half-span
    #     lift_dist = L * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)
    #     drag_dist = D * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)
    #     weight_dist = self.weight_estimate * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)

    #     # Add point masses for ducted fans
    #     thrust_dist = np.zeros_like(y_span)
    #     fan_positions = [b/6, b*2/6, b*4/6, b*5/6]  # example positions
    #     thrust_per_fan = T
    #     for pos in fan_positions:
    #         idx = (np.abs(y_span - pos)).argmin()
    #         thrust_dist[idx] += thrust_per_fan

    #     # Bending moments
    #     M_y = self.trapz_manual(drag_dist * y_span, y_span)  # about y-axis
    #     M_z = self.trapz_manual(lift_dist * y_span, y_span) + self.trapz_manual(weight_dist * y_span, y_span) - self.trapz_manual(thrust_dist * y_span, y_span)

    #     # Moments of inertia
    #     I_y = (h*w**3)/12 - ((h-2*t)*(w-2*t)**3)/12
    #     I_z = (w*h**3)/12 - ((w-2*t)*(h-2*t)**3)/12

    #     # Section distances from neutral axis
    #     z = w/2
    #     y = h/2

    #     # Axial stresses
    #     axial_yy = -(M_z * y)/I_z
    #     axial_zz = (M_y * z)/I_y
    #     axial_max = max(abs(axial_yy), abs(axial_zz))

    #     return axial_max

    # def shear_stress(self, L, T, D):
    #     # Span and chord
    #     b = np.sqrt(aircraft.AR*aircraft.s_ref)/2
    #     c = np.sqrt(aircraft.s_ref/aircraft.AR)
    #     x_1 = 0.1*c
    #     x_2 = 0.7*c

    #     t_x1 = self.aero["t_x1"]
    #     t_x2 = self.aero["t_x2"]
    #     h = 0.5*(t_x1 + t_x2)
    #     w = x_2 - x_1
    #     t = t_x1 * 0.1  # box thickness estimate

    #     # Span discretization
    #     n_span = 50
    #     y_span = np.linspace(0, b, n_span)

    #     # Elliptical distributions along half-span
    #     lift_dist = L * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)
    #     drag_dist = D * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)
    #     weight_dist = self.weight_estimate * (4 / (np.pi * b**2)) * np.sqrt(b**2 - y_span**2)

    #     # Add point thrusts
    #     thrust_dist = np.zeros_like(y_span)
    #     fan_positions = [b/6, b*2/6, b*4/6, b*5/6]  # example positions
    #     thrust_per_fan = T
    #     for pos in fan_positions:
    #         idx = (np.abs(y_span - pos)).argmin()
    #         thrust_dist[idx] += thrust_per_fan

    #     # Shear forces
    #     F_y = self.trapz_manual(weight_dist, y_span) + self.trapz_manual(thrust_dist, y_span) - self.trapz_manual(lift_dist, y_span)
    #     F_z = self.trapz_manual(drag_dist, y_span)

    #     # First moments of area
    #     Q_y = (h*w**2)/8 - ((h-2*t)*(w-2*t)**2)/8
    #     Q_z = (w*h**2)/8 - ((w-2*t)*(h-2*t)**2)/8

    #     # Moments of inertia
    #     I_y = (h*w**3)/12 - ((h-2*t)*(w-2*t)**3)/12
    #     I_z = (w*h**3)/12 - ((w-2*t)*(h-2*t)**3)/12

    #     shear_xy = (F_y * Q_z)/(I_z * t)
    #     shear_xz = (F_z * Q_y)/(I_y * t)
    #     shear_max = max(abs(shear_xy), abs(shear_xz))

    #     return shear_max

    def axial_stress(self, L, T, D):
        # Get required variables
        b = np.sqrt(aircraft.AR*aircraft.s_ref)/2
        c = np.sqrt(aircraft.s_ref/aircraft.AR)
        x_T1 = b/6
        x_T2 = b*(2/6)
        x_T3 = b*(3/6)
        x_T4 = b*(4/6)
        x_T5 = b*(5/6)
        t_x1 = self.aero["t_x1"]
        t_x2 = self.aero["t_x2"]
        x_1 = 0.1*c
        x_2 = 0.7*c

        h = 0.5*(t_x1+t_x2)
        w = x_2-x_1
        t = t_x1*0.1 # box thickness - estimate

        # NOTE: not sure if this is right z and y are dist. from the neutral axis
        z = w/2
        y = h/2


        # # Moment equations for final model:
        # M_y = T*(x_T1+x_T2+x_T3+x_T4) # - drag integral
        # M_z = self.weight_estimate*(b/4) + W_duct*(x_T1+x_T2+x_T3+x_T4) # - lift integral
        # # TODO: add actual integral stuff for D, L, W

        # Moment equations as of now (3/8) to test model, assume uniform distribution:
        M_y = T*(x_T1+x_T2+x_T3+x_T4+x_T5) - (D*(b)**2)/2
        M_z = self.weight_estimate*(b/2) + self.W_duct*(x_T1+x_T2+x_T3+x_T4+x_T5) + (self.W_fuel_per_wing*(b)**2)- (L*(b)**2)/2

        # Moments of inertia
        I_y = (h*w**3)/12 - ((h-2*t)*(w-2*t)**3)/12
        I_z = (w*h**3)/12 - ((w-2*t)*(h-2*t)**3)/12

        # Stress equations
        axial_yy = -(M_z*y)/I_z
        axial_zz = (M_y*z)/I_y

        axial_max = max(abs(axial_yy), abs(axial_zz))

        # return axial_yy + axial_zz
        return axial_max

    def shear_stress(self, L, T, D):
        # Get required variables
        b = np.sqrt(aircraft.AR*aircraft.s_ref)/2
        c = np.sqrt(aircraft.s_ref/aircraft.AR)
        x_T1 = b/6
        x_T2 = b*(2/6)
        x_T3 = b*(3/6)
        x_T4 = b*(4/6)
        x_T5 = b*(5/6)
        t_x1 = self.aero["t_x1"]
        t_x2 = self.aero["t_x2"]
        x_1 = 0.1*c
        x_2 = 0.7*c

        h = 0.5*(t_x1+t_x2)
        w = x_2-x_1
        t = t_x1*0.1 # box thickness - estimate

        # Force equations
        F_y = self.weight_estimate*(b/2) + 5*self.W_duct + self.W_fuel_per_wing - L # - lift integral
        F_z = -4*T + D # + drag integral
        # TODO: add actual integral stuff for D, L, W

        # 1st moments of area
        Q_y = (h*w**2)/8 - ((h-2*t)*(w-2*t)**2)/8
        Q_z = (w*h**2)/8 - ((w-2*t)*(h-2*t)**2)/8

        # Moments of inertia
        I_y = (h*w**3)/12 - ((h-2*t)*(w-2*t)**3)/12
        I_z = (w*h**3)/12 - ((w-2*t)*(h-2*t)**3)/12

        shear_xy = (F_y*Q_z)/(I_z*t)
        shear_xz = (F_z*Q_y)/(I_y*t)

        shear_max = max(abs(shear_xy), abs(shear_xz))

        # return shear_xy + shear_xz
        return shear_max


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
        s_tot = 2 * c_ail
        # s_tot = self.aero["airfoil_perimeter"] # we'll estimate this as 2 * c_ail for now
        G = self.materials["skin_G"]
        twist_max = self.aero["twist_max"]
        A = b_ail*c_ail
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
        # landing_tube_thickness = self.tube_thickness()

        return landing_spar_cap_area, landing_spar_web_area, landing_skin_thickness # , landing_tube_thickness

    def max_load_sizing(self):
        sizing_takeoff = self.takeoff()
        sizing_climb = self.climb()
        sizing_cruise = self.cruise()
        sizing_descent = self.descent()
        sizing_landing = self.landing()

        print("sizing_takeoff:", sizing_takeoff)
        print("sizing_climb:", sizing_climb)
        print("sizing_cruise:", sizing_cruise)
        print("sizing_descent:", sizing_descent)
        print("sizing_landing:", sizing_landing)

        spar_cap_area = max(sizing_takeoff[0], sizing_climb[0], sizing_cruise[0], sizing_descent[0], sizing_landing[0])
        spar_web_area = max(sizing_takeoff[1], sizing_climb[1], sizing_cruise[1], sizing_descent[1], sizing_landing[1])
        skin_thickness = max(sizing_takeoff[2], sizing_climb[2], sizing_cruise[2], sizing_descent[2], sizing_landing[2])
        # tube_thickness = max(sizing_takeoff[3], sizing_climb[3], sizing_cruise[3], sizing_descent[3], sizing_landing[3])

        return spar_cap_area, spar_web_area, skin_thickness # , tube_thickness


    def wing_weight(self):
        b = np.sqrt(aircraft.AR*aircraft.s_ref)/2
        # for one of two wings
        sizing = self.max_load_sizing()
        print(sizing)
        spar_cap_weight = sizing[0]*b*self.materials["spar_cap_density"]
        spar_web_weight = sizing[1]*b*self.materials["spar_web_density"]
        skin_weight = sizing[2]*self.aero["airfoil_surface_area"]*self.materials["skin_density"]
        # tube_weight = sizing[3]*b*self.materials["tube_density"]

        print("spar_cap_weight", spar_cap_weight)
        print("spar_web_weight", spar_web_weight)
        print("skin_weight", skin_weight)

        # should this weight just be structural stuff or also fuel and motors?
        return spar_cap_weight + spar_web_weight + skin_weight
        # spar_cap_weight + spar_web_weight should = tube_weight right??


class Tail():

    def __init__(self, aero, loading, materials, weight_estimate):
        self.aero = aero
        self.loading = loading
        self.material = materials
        self.weight_estimate = weight_estimate

    def axial_stress_vert_ho (self, Lv, Lh, Dv, Dh, tv, th, hweight):
        # Get required variables
        b_vert = self.aero["b_vert"]
        b_ho = self.aero["b_ho"]

        t_x1_v = self.aero["t_x1_v"]
        t_x2_v = self.aero["t_x2_v"]
        x_1_v = self.aero["x_1_v"]
        x_2_v = self.aero["x_2_v"]

        x_1_h = self.aero["x_1_h"]
        x_2_h = self.aero["x_2_h"]
        t_x1_h = self.aero["t_x1_h"]
        t_x2_h = self.aero["t_x2_h"]

        #vertical tail calcs:
        hv = 0.5*(t_x1_v+t_x2_v) #along y axis
        wv = x_2_v-x_1_v #along x axis
        crossA_v = (hv*wv)-(hv-2*tv)*(wv-2*tv)

        #horizontal tail calcs:
        hh = 0.5*(t_x1_h+t_x2_h) # along z axis
        wh = x_2_h-x_1_h # along x axis
        crossA_h = (hh*wh)-(hh-2*th)*(wh-2*th)

        # NOTE: not sure if this is right z and y are dist. from the neutral axis
        xv = wv/2 #vert
        yv = hv/2 #vert

        xh = wh/2 #ho
        zh = hh/2 #ho

        # Moment for vertical tail (accounting for T-tail, but not including angle):
        M_yv = (Dv*(b_vert/2)**2)/2 #moment from drag of v and h
        M_xv = (Lv*(b_vert/2)**2)/2 #lift moment

        # Moment for horizontal tail (uniform lift and drag for now- 3/9)
        M_xh = (Lh*(b_ho/2)**2)/2
        M_zh = (Dh*(b_ho/2)**2)/2

        #compressive and tensile loads on vertical tail along z (not including angle):
        V_compress = hweight  #h weight
        V_tensile = Lh #h lift

        # Moments of inertia vertical
        I_y_v = (hv*wv**3)/12 - ((hv-2*tv)*(wv-2*tv)**3)/12
        I_x_v = (wv*hv**3)/12 - ((wv-2*tv)*(hv-2*tv)**3)/12

        # Moments of inertia horizontal
        I_x_v = (hh*wh**3)/12 - ((hh-2*th)*(wh-2*th)**3)/12
        I_z_v = (wh*hh**3)/12 - ((wh-2*th)*(hh-2*th)**3)/12

        # Stress equations vertical NOTE check again
        axial_yy_v = -(M_xv*yv)/I_x_v
        axial_xx_v = ((M_yv*xv)/I_y_v)+(V_tensile/crossA_v)-(V_compress/crossA_v)

        # Stress equations horizontal
        axial_xx_h = (M_xh*zh)/I_x_v
        axial_zz_h = ((M_zh*xh)/I_z_v)


        axial_max_v = max(abs(axial_xx_v), abs(axial_yy_v))
        axial_max_h = max(abs(axial_xx_h), abs(axial_zz_h))

        return axial_max_v, axial_max_h

    def shear_stress_vert_ho (self, Lv, Lh, Dv, Dh, tv, th, hweight):
        b_vert = self.aero["b_vert"]
        b_ho = self.aero["b_ho"]

        t_x1_v = self.aero["t_x1_v"]
        t_x2_v = self.aero["t_x2_v"]
        x_1_v = self.aero["x_1_v"]
        x_2_v = self.aero["x_2_v"]

        x_1_h = self.aero["x_1_h"]
        x_2_h = self.aero["x_2_h"]
        t_x1_h = self.aero["t_x1_h"]
        t_x2_h = self.aero["t_x2_h"]

        #vertical tail calcs:
        hv = 0.5*(t_x1_v+t_x2_v) #along y axis
        wv = x_2_v-x_1_v #along x axis
        crossA_v = (hv*wv)-(hv-2*tv)*(wv-2*tv)

        #horizontal tail calcs:
        hh = 0.5*(t_x1_h+t_x2_h) # along z axis
        wh = x_2_h-x_1_h # along x axis
        crossA_h = (hh*wh)-(hh-2*th)*(wh-2*th)

        # Force equations vertical
        F_y_v = Lv
        F_x_v = Dv + Dh
        # TODO: add actual integral stuff for D, L, W

        # Force equations horizontal
        F_z_h = Lh
        F_x_h = Dh

        # 1st moments of area vertical
        Q_x = (hv*wv**2)/8 - ((hv-2*tv)*(wv-2*tv)**2)/8
        Q_y = (wv*hv**2)/8 - ((wv-2*tv)*(hv-2*tv)**2)/8

        # 1st moments of area horizontal
        Q_x = (hh*wh**2)/8 - ((hh-2*th)*(wh-2*th)**2)/8
        Q_z = (wh*hh**2)/8 - ((wh-2*th)*(hh-2*th)**2)/8

        # Moments of inertia vertical
        I_y_v = (hv*wv**3)/12 - ((hv-2*tv)*(wv-2*tv)**3)/12
        I_x_v = (wv*hv**3)/12 - ((wv-2*tv)*(hv-2*tv)**3)/12

        # Moments of inertia horizontal
        I_x_h = (hh*wh**3)/12 - ((hh-2*th)*(wh-2*th)**3)/12
        I_z_h = (wh*hh**3)/12 - ((wh-2*th)*(hh-2*th)**3)/12

        shear_zy_v = (F_y_v*Q_x)/(I_x_v*tv)  #shear zy
        shear_zx_v = (F_x_v*Q_y)/(I_y_v*tv)  #shear zx

        shear_yz_h = (F_z_h*Q_x)/(I_x_h*th)  #shear yz
        shear_zx_h = (F_x_h*Q_z)/(I_z_h*th)  #shear yx


        shear_max_v = max(abs(shear_zy_v), abs(shear_zx_v))
        shear_max_h = max(abs(shear_yz_h), abs(shear_zx_h))

        return shear_max_v, shear_max_h
    

    def get_N_vert_ho(self, L, a_z, cg, fuse_rad, pitch, yaw):

        """
        input: L, a_z, mass, fuselage radius, pitch (deg), yaw (deg)
        Output: flight load factor, landing load factor, and max load factor
        """
        # Get required variable
        W_total = self.aero["W_total"]
        fuse_I = 0.5*W_total*(fuse_rad**2)
        # moment = Lh_req*fuse_I = 
        Lh_req = ...
        Lv_req = ...
        # Calculate max load factor
        Nh = Lh_req/W_total # flight load factor
        # N_land = a_z/self.g # landing load factor
        # N_max = np.max(N, N_land) # maximum load factor
        Nv = Lh_req/W_total

        return Nv, Nh

    def spar_cap_area_vert_ho(self, L, a_z, axial_stress_v, axial_stress_h):
        """Find necessary spar cap area based on loading"""
        # global variables 
        W_cent = self.aero["W_total"] - self.weight_estimate

        # vertical get required variables 
        c_tip_v = self.aero["c_tip_v"]
        c_0_v = self.aero["c_0_v"]
        taper_ratio_v = c_tip_v/c_0_v
        Nv = self.get_N(L, a_z)[0]
        bv = self.aero["b_vert"] # NOTE: might need to do b/2 for half of wing?
        hv = self.aero["h_vert"]
        E = self.materials["spar_car_E"] # maybe remain this variable if it shows up again just for clarity
        w_b_max_v = self.aero["w_b_max_v"] # not sure where we'll actually get this

        # horizontal get required variables 
        c_tip_h = self.aero["c_tip_h"]
        c_0_h = self.aero["c_0_h"]
        taper_ratio_h = c_tip_h/c_0_h
        Nh = self.get_N(L, a_z)[1]
        bh = self.aero["b_ho"] # NOTE: might need to do b/2 for half of wing?
        hh = self.aero["h_ho"]
        E = self.materials["spar_car_E"] # maybe remain this variable if it shows up again just for clarity
        w_b_max_h = self.aero["w_b_max_h"] # not sure where we'll actually get this

        # Vertical Calculate area based on strength and stiffness requirements 
        A_cap_0_strength_v = (Nv*W_cent)/(12*axial_stress_v)*(bv/hv)*(1+2*taper_ratio_v)/(1+taper_ratio_v)
        A_cap_0_stiffness_v = (Nv*W_cent)/(48*E) * bv**2/hv**2 * 1/w_b_max_v * (1+2*taper_ratio_v)/(1+taper_ratio_v)

        # Horizontal Calculate area based on strength and stiffness requirements 
        A_cap_0_strength_h = (Nh*W_cent)/(12*axial_stress_h)*(bh/hh)*(1+2*taper_ratio_h)/(1+taper_ratio_h)
        A_cap_0_stiffness_h = (Nh*W_cent)/(48*E) * bh**2/hh**2 * 1/w_b_max_h * (1+2*taper_ratio_h)/(1+taper_ratio_h)
        return np.max(A_cap_0_strength_v, A_cap_0_stiffness_v), np.max(A_cap_0_strength_h, A_cap_0_stiffness_h)

    def spar_web_area_vert_ho(self, L, a_z, shear_stress_v, shear_stress_h):
        """Find necessary spar web area based on loading"""
        # Get required variables
        Nv = self.get_N(L, a_z)[0]
        Nh = self.get_N(L, a_z)[1]
        W_cent = self.aero["W_total"] - self.weight_estimate

        # Return spar web area
        return (Nv*W_cent)/(2*shear_stress_v), (Nh*W_cent)/(2*shear_stress_h)

    def skin_thickness(self, q, shear_stress_v, shear_stress_h):
        """Find necessary skin thickness based on loading"""
        # Calculated same way as main wing 
        # May need editing since the tail doesn't have ailerons driving torsion

        #global variables 
        y_ail = self.aero["y_ail"] #not sure how to do this with tail 
        c_m = self.aero["c_m"] # I don't think this changes based on flight conditions but not sure
        s_tot = self.aero["s_tot"] # not sure what this is
        G = self.materials["skin_G"]
        twist_max = self.aero["twist_max"]

        # Vertical required variables
        bv = self.aero["b_vert"]
        c_tip_v = self.aero["c_tip_v"]
        c_0_v = self.aero["c_0_v"]
        cv = (c_tip_v+c_0_v)*0.5
        Av = bv*cv

        # Horizontal required variables
        bh = self.aero["b_ho"]
        c_tip_h = self.aero["c_tip_h"]
        c_0_h = self.aero["c_0_h"]
        ch = (c_tip_h+c_0_h)*0.5
        Ah = bh*ch


        # Vertical calculate skin thickness requirements
        t_skin_strength_v = q*bv*cv**2*c_m*1/(2*Av)*1/(1/shear_stress_v)
        t_skin_stiffness_v = 1*bv*cv**2*y_ail*c_m*s_tot/(4*Av**2*G)*(1/twist_max)

        # Horizontal calculate skin thickness requirements
        t_skin_strength_h = q*bh*ch**2*c_m*1/(2*Ah)*1/(1/shear_stress_h)
        t_skin_stiffness_h = 1*bh*ch**2*y_ail*c_m*s_tot/(4*Ah**2*G)*(1/twist_max)

        # Return max thickness
        return np.max(t_skin_strength_v, t_skin_stiffness_v), np.max(t_skin_strength_h, t_skin_stiffness_h)
    
    def ho_max_pitch_velNE(self, velNE, Cl_v_max, Cd_v):
        # Find forces - do we need to incorporate some takeoff angle into this? maybe for all the stuff angle is an option and during cruise it's just zero?
        bv = self.aero["b_vert"]
        c_tip_v = self.aero["c_tip_v"]
        c_0_v = self.aero["c_0_v"]
        cv = (c_tip_v+c_0_v)*0.5
        Av = bv*cv
        L_req = 0.5*self.rho_cruise*(velNE**2)*Cl_v_max*Av
        D = 0.5*self.rho_cruise*(velNE**2)*Cd_v*Av


        # Calculate stresses and torsion
        def axial_stress_vert_ho (self, Lv, Lh, Dv, Dh, tv, th, hweight):
        axial_stress_takeoff = self.axial_stress_vert_ho(L_takeoff, T_takeoff, D_takeoff)
        shear_stress_takeoff = self.shear_stress(L_takeoff, T_takeoff, D_takeoff)

        # Find component sizing based on calculated loading
        # NOTE: setting a_z = 0 on all cases but landing so that N_land is not considered
        takeoff_spar_cap_area = self.spar_cap_area(L_takeoff, 0, axial_stress_takeoff)
        takeoff_spar_web_area = self.spar_web_area(L_takeoff, 0, shear_stress_takeoff)
        takeoff_skin_thickness = self.skin_thickness(q_takeoff, shear_stress_takeoff)
        takeoff_tube_thickness = self.tube_thickness() # still not sure what this is for

        return takeoff_spar_cap_area, takeoff_spar_web_area, takeoff_skin_thickness, takeoff_tube_thickness


    def climb(self):
        # Find forces - do we need to incorporate climb angle into this? maybe for all the stuff angle is an option and during cruise it's just zero?
        q_climb = 0.5*(0.5*(self.rho_ground+self.rho_cruise))*self.aero["v_climb"]**2
        L_climb = self.aero["CL_climb"]*q_climb*self.aero["S"] # or lift distribution from aero
        T_climb = self.loading["T_climb"] # not sure where this will come from exactly
        D_climb = self.aero["CD_climb"]*q_climb*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_climb = self.axial_stress(L_climb, T_climb, D_climb)
        shear_stress_climb = self.shear_stress(L_climb, T_climb, D_climb)

        # Find component sizing based on calculated loading
        # NOTE: setting a_z = 0 on all cases but landing so that N_land is not considered
        climb_spar_cap_area = self.spar_cap_area(L_climb, 0, axial_stress_climb)
        climb_spar_web_area = self.spar_web_area(L_climb, 0, shear_stress_climb)
        climb_skin_thickness = self.skin_thickness(q_climb, shear_stress_climb)
        climb_tube_thickness = self.tube_thickness()

        return climb_spar_cap_area, climb_spar_web_area, climb_skin_thickness, climb_tube_thickness

    def cruise(self):
        # Find forces
        q_cruise = 0.5*self.rho_cruise*self.aero["v_cruise"]**2
        L_cruise = self.aero["CL_cruise"]*q_cruise*self.aero["S"] # or lift distribution from aero
        T_cruise = self.loading["T_cruise"] # not sure where this will come from exactly
        D_cruise = self.aero["CD_cruise"]*q_cruise*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_cruise = self.axial_stress(L_cruise, T_cruise, D_cruise)
        shear_stress_cruise = self.shear_stress(L_cruise, T_cruise, D_cruise)

        # Find component sizing based on calculated loading
        # NOTE: setting a_z = 0 on all cases but landing so that N_land is not considered
        cruise_spar_cap_area = self.spar_cap_area(L_cruise, 0, axial_stress_cruise)
        cruise_spar_web_area = self.spar_web_area(L_cruise, 0, shear_stress_cruise)
        cruise_skin_thickness = self.skin_thickness(q_cruise, shear_stress_cruise)
        cruise_tube_thickness = self.tube_thickness()

        return cruise_spar_cap_area, cruise_spar_web_area, cruise_skin_thickness, cruise_tube_thickness

    def descent(self):
        # Find forces - do we need to incorporate descent angle into this? maybe for all the stuff angle is an option and during cruise it's just zero?
        q_descent = 0.5*(0.5*(self.rho_ground+self.rho_cruise))*self.aero["v_descent"]**2
        L_descent = self.aero["CL_descent"]*q_descent*self.aero["S"] # or lift distribution from aero
        T_descent = self.loading["T_descent"] # not sure where this will come from exactly
        D_descent = self.aero["CD_descent"]*q_descent*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_descent = self.axial_stress(L_descent, T_descent, D_descent)
        shear_stress_descent = self.shear_stress(L_descent, T_descent, D_descent)

        # Find component sizing based on calculated loading
        # NOTE: setting a_z = 0 on all cases but landing so that N_land is not considered
        descent_spar_cap_area = self.spar_cap_area(L_descent, 0, axial_stress_descent)
        descent_spar_web_area = self.spar_web_area(L_descent, 0, shear_stress_descent)
        descent_skin_thickness = self.skin_thickness(q_descent, shear_stress_descent)
        descent_tube_thickness = self.tube_thickness()

        return descent_spar_cap_area, descent_spar_web_area, descent_skin_thickness, descent_tube_thickness

    def landing(self):
        # Find forces - do we need to incorporate landing angle into this? maybe for all the stuff angle is an option and during cruise it's just zero?
        q_landing = 0.5*self.rho_ground*self.aero["v_landing"]**2
        L_landing = self.aero["CL_landing"]*q_landing*self.aero["S"] # or lift distribution from aero
        T_landing = self.loading["T_landing"] # not sure where this will come from exactly
        D_landing = self.aero["CD_landing"]*q_landing*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_landing = self.axial_stress(L_landing, T_landing, D_landing)
        shear_stress_landing = self.shear_stress(L_landing, T_landing, D_landing)

        # Find component sizing based on calculated loading
        # NOTE: setting a_z = 0 on all cases but landing so that N_land is not considered
        landing_spar_cap_area = self.spar_cap_area(L_landing, self.aero["a_z"], axial_stress_landing)
        landing_spar_web_area = self.spar_web_area(L_landing, self.aero["a_z"], shear_stress_landing)
        landing_skin_thickness = self.skin_thickness(q_landing, shear_stress_landing)
        landing_tube_thickness = self.tube_thickness()

        return landing_spar_cap_area, landing_spar_web_area, landing_skin_thickness, landing_tube_thickness

    def max_load_sizing(self):
        sizing_takeoff = self.takeoff()
        sizing_climb = self.climb()
        sizing_cruise = self.cruise()
        sizing_descent = self.descent()
        sizing_landing = self.landing()

        spar_cap_area = max(sizing_takeoff[0], sizing_climb[0], sizing_cruise[0], sizing_descent[0], sizing_landing[0])
        spar_web_area = max(sizing_takeoff[1], sizing_climb[1], sizing_cruise[1], sizing_descent[1], sizing_landing[1])
        skin_thickness = max(sizing_takeoff[2], sizing_climb[2], sizing_cruise[2], sizing_descent[2], sizing_landing[2])
        tube_thickness = max(sizing_takeoff[3], sizing_climb[3], sizing_cruise[3], sizing_descent[3], sizing_landing[3])

        return spar_cap_area, spar_web_area, skin_thickness, tube_thickness


    def wing_weight(self):
        b = self.aero["b"]
        # for one of two wings
        sizing = self.max_load_sizing()
        spar_cap_weight = sizing[0]*b*self.materials["spar_cap_density"]
        spar_web_weight = sizing[1]*b*self.materials["spar_web_density"]
        skin_weight = sizing[2]*self.aero["airfoil_surface_area"]*self.materials["skin_density"]
        tube_weight = sizing[3]*b*self.materials["tube_density"]

        spar_weight = max(spar_cap_weight + spar_web_weight, tube_weight)

        # should this weight just be structural stuff or also fuel and motors?
        return spar_weight + skin_weight
        # spar_cap_weight + spar_web_weight should = tube_weight right??






class Fuselage:

    def __init__(self, length, radius, n_people, cabin_width):
        self.length = length
        self.R = radius
        self.n = n_people
        self.cabin_width = cabin_width
        self.weight = 0.0 # to be derived



    def pressure_at_altitude(self, h):
        """
        Returns atmospheric pressure (Pa) at altitude h (meters)
        Valid up to 11 km (troposphere)
        """
        P0 = 101325      # Sea level pressure (Pa)
        T0 = 288.15      # Sea level temperature (K)
        L  = 0.0065      # Lapse rate (K/m)
        g  = 9.80665     # Gravity (m/s^2)
        R  = 287.05      # Gas constant for air (J/kg*K)

        return P0 * (1 - (L * h) / T0)**(g / (R * L))




    # -------------------------
    # Geometry functions
    # -------------------------

    def shell_area(self):
        return self.circumference() * self.length

    def circumference(self):
        return 2.0 * np.pi * self.R

    def count_members(self, length, spacing):
        return int(np.floor(length / spacing)) + 1

    def required_thickness_hoop(self,
                        yield_strength,
                        safety_factor=2):
        """
        Computes required wall thickness required (meters)

        radius: cabin radius (m)
        altitude: flight altitude (m)
        cabin_pressure: desired internal pressure (Pa)
        yield_strength: material yield strength (Pa)
        safety_factor: structural safety factor
        """

        max_altitude = 5500 # m (~18,000 ft)
        cabin_pressure = 75150   # Pa (~8,000 ft)

        # Outside pressure
        P_out = self.pressure_at_altitude(max_altitude)

        # Pressure differential
        delta_P = cabin_pressure - P_out

        # Allowable stress
        sigma_allow = yield_strength / safety_factor

        # Thin wall hoop stress formula
        self.hoop_t = (delta_P * self.R) / sigma_allow

        return self.hoop_t

    def required_thickness_moments(self, yield_stress, safety_factor=2):
        pass



    
    # -------------------------
    # Component mass functions
    # -------------------------


    def get_structural_mass(self):
        # assume some materials
        # aluminum 7075-T6
        yield_strength_al70 = 490e6    # Pa
        rho_al70 = 2810    # kg/m^3

        # honey comb fiberglass/carbonfiber
        rho_hc = 3 # kg/m^3

        # skin mass
        # skin mass
        # self.thickness = self.required_thickness_hoop(yield_strengths)+self.required_thickness_moments(yield_strengths) # solve for later
        thickness = 0.002 # meter
        skin_volume = self.cylinder_area(self.R)*thickness
        self.skin_mass = skin_volume*rho_al70

        # frame mass
        # frame specs
        frame_spacing = .5 # assumed
        frame_disk_depth = .1 # assumed
        frame_thickness = 0.002 # assumed
        frame_area = (np.pi*self.R**2-np.pi*(self.R-frame_disk_depth)**2)
        n_frames = self.count_members(self.length, frame_spacing)
        total_volume = n_frames * frame_thickness * frame_area
        self.frame_mass = total_volume * rho_al70


        # stringer mass
        # stringer specs
        stringer_area = 1e-4 # assumed
        stringer_spacing = 0.25 # assumed
        # stringer mass
        n_stringers = self.count_members(self.circumference(), stringer_spacing)
        total_volume = n_stringers * self.length * stringer_area
        self.stringer_mass = total_volume * rho_al70

        # floor panel mass
        # floor specs
        floor_thickness = 0.01 # assumed
        area = self.length * self.cabin_width
        self.fpanel_mass = area * floor_thickness * rho_hc
    

        # floor fframe mass
        # floor fframe specs
        n_fframes = n_frames # assume beams at every frame
        fframe_depth = 0.07 # assumed
        fframe_thickness = 0.002 # assumed
        total_volume = n_fframes * self.cabin_width * fframe_thickness * fframe_depth
        self.fframe_mass = total_volume * rho_al70


        
        # floor beam mass
        # floor beam specs
        beam_spacing = .5
        n_fbeams = self.count_members(self.cabin_width, beam_spacing)
        fbeam_depth = 0.07 # assumed
        fbeam_thickness= 0.002 # assumed
        total_volume = n_fbeams * self.length * fbeam_thickness * fbeam_depth
        self.fbeam_mass = total_volume * rho_al70


        return self.skin_mass+self.frame_mass+self.stringer_mass+self.fpanel_mass+self.fframe_mass+self.fbeam_mass
    
    def get_dead_mass(self):
        seat_weight = self.n*13.0 # kg (average modern aircraft seat weight)
        person_weight = self.n*100.0 # kg (person + luggage)
        return seat_weight + person_weight


    def get_total_weight(self):
        g  = 9.80665     # Gravity (m/s^2)
        return self.get_dead_mass()*g + self.get_structural_mass()*g



class LandingGear():
    pass


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
        "CL_climb": 1.5,
        "CL_cruise": 0.35,
        "CL_descent": 0.2,
        "CL_landing": 6,
        "CD_takeoff": 0.7,
        "CD_climb": 0.1,
        "CD_cruise": 0.05,
        "CD_descent": 0.045,
        "CD_landing": 1,
        "v_takeoff": 25,
        "v_climb": 75,
        "v_cruise": 125,
        "v_descent": 75,
        "v_landing": 25,
        "drag_area": 2.5,
        "a_z": 2*9.8,
        "airfoil_surface_area": 100, # m^2
        # "b_vert":,
        # "b_ho":,
        # "t_x1_v":,
        # "t_x2_v":,
        # "x_1_v":,
        # "x_2_v":,
        # "t_x1_h":,
        # "t_x2_h":,
        # "x_1_h":,
        # "x_2_h":,
    }

    loading = {
        "T_takeoff": 60 * 10**3,
        "T_climb": 40 * 10**3, # guess
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



# # Assumptions/variables
# L = ... # lift
# W_total = ... # total weight
# W_wing = ... # wing weight
# W_cent = W_total - W_wing # center weight
# a_z = ... # acceleration at landing
# g = 9.81
# w_b_max = ... # maximum deflection to span ratio - not sure if this should be assumed or calculated
# rho = ... # density (double check if we need material vs. air density in different places)
# V = ... # airspeed
# q = 0.5*rho*V**2 # dynamic pressure
# twist_max = ... # maximum tolerable tip twist angle

# # Sizing Parameters - done
# b = ... # span
# h = ... # spar height
# c_tip = ... # chord at the tip
# c_0 = ... # chord at the root
# b_ail = ... # span from root to aileron?
# c_ail = ... # chord at the aileron
# c_m = ... # local pitching moment coefficient
# A = ... # shell enclosed area
# y_ail = ... # spanwise coordinate of aileron?
# s_tot = ... # not sure exactly what this is
# R_outer = ... # tube outer radius
# l = ... # length

# # Elasticity Values
# axial_max = ... # maximum axial stress
# shear_max = ... # maximum shear stress
# torsion_max = ... # maximum torsional moment
# bending_max = ... # maximum bending moment
# E = ... # stiffness modulus
# G = ... # shear modulus
# I = ... # bending moment of inertia
# C = 0.25 # depends on how tube ends are supported (this is for fixed at one end only)
# P_max = C * (np.pi**2*E*I)/l**2 # not sure what the difference is between P_cr and P_max in the pdf


# # Load Factor
# N = L/W_total # flight load factor
# N_land = a_z/g # landing load factor
# N_max = np.max(N, N_land) # maximum load factor



# # Spar sizing
# taper_ratio = c_tip/c_0
# A_cap_0_strength = (N_max*W_cent)/(12*axial_max)*(b/h)*(1+2*taper_ratio)/(1+taper_ratio)
# A_cap_0_stiffness = (N*W_cent)/(48*E) * b**2/h**2 * 1/w_b_max * (1+2*taper_ratio)/(1+taper_ratio)
# A_cap_0 = np.max(A_cap_0_strength, A_cap_0_stiffness)

# A_web_0 = (N_max*W_cent)/(2*shear_max)

# # Skin Sizing
# t_skin_strength = q*b_ail*c_ail**2*c_m*1/(2*A)*1/(1/shear_max)
# t_skin_stiffness = 1*b_ail*c_ail**2*y_ail*c_m*s_tot/(4*A**2*G)*(1/twist_max)

# t_skin = np.max(t_skin_strength, t_skin_stiffness)

# # Tube Sizing
# t_tube_torsion_strength = torsion_max/(2*np.pi*R_outer**2*shear_max)
# t_tube_torsion_stiffness = (torsion_max*l)/(2*np.pi*R_outer**3*G)

# t_tube_bending_strength = bending_max/(np.pi*R_outer**2*axial_max)

# # NOTE: I think more complex equations need to be written out for k, F, and w depending on which we case we choose
# t_tube_bending_deflection = (k*F_max*l**3)/(np.pi*R_outer**2*w_max*E)

# t_tube_buckling = 1/C * (P_max*l**2)/(np.pi**3*R_outer**3*E)

# tube_weight = rho*g*(torsion_max*l**2)/(R_outer**2*G) # not sure if we should do this calc or not

# t_tube = np.max(t_tube_torsion_strength, t_tube_torsion_stiffness, t_tube_bending_strength, t_tube_bending_deflection, t_tube_buckling)
