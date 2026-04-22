# Wing Sizing

import numpy as np
import math as math
import matplotlib.pyplot as plt

# using aero_main dataclass here with just the values they gave us
from dataclasses import dataclass

import numpy as np
import pandas as pd
# from ambiance import Atmosphere
from aero_workspace.conceptual_design import MTOW, V_CRUISE, V_STALL, W_S, S, ureg
from scipy.interpolate import interp1d


# TODO: update constants in dataclass
@dataclass
class aircraft:
    """All aircraft definitions"""

    #### GLOBAL DEFINITIONS ###
    AR: int = 8
    s_ref: float = 49.6
    b: float = 19.92

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
    W_total = 7500 * 9.8
    W_duct = 25 * 9.81 # kg - guess
    W_fuel_per_wing = 750 * 9.8 # kg

    def __init__(self, aero, loading, materials, weight_estimate = 0.07*W_total):
        self.aero = aero
        self.loading = loading
        self.materials = materials
        self.weight_estimate = weight_estimate

    def axial_stress(self):
        return 276e6


    def shear_stress(self):
        return 207e6

    def max_moment(self, L, N_landing, strut_loc, strut_area):

        # Geometry
        b = np.sqrt(aircraft.AR * aircraft.s_ref) / 2

        # Fan locations
        x_T1 = b/5
        x_T2 = b*(2/5)
        x_T3 = b*(3/5)
        x_T4 = b*(4/5)

        # weights
        W_struct = N_landing * self.weight_estimate
        W_fuel = N_landing * self.W_fuel_per_wing
        W_duct = N_landing * self.W_duct

        # Elliptical integrals
        lift_moment = (2 * N_landing * L * b) / (3 * np.pi)
        struct_moment = (2 * W_struct * b) / (3 * np.pi)
        fuel_moment = (2 * W_fuel * b) / (3 * np.pi)

        # Base cantilever moment
        M_cantilever = (
            struct_moment
            + fuel_moment
            + W_duct * (x_T1 + x_T2 + x_T3 + x_T4)
            - lift_moment
        )

        # --- STRUT FORCE ---
        E_strut = 6.9e10
        L_strut = np.sqrt(1.1**2 + (b * strut_loc)**2)
        print("L_strut", L_strut)
        A_strut = strut_area
        d_strut = 2*np.sqrt(A_strut/np.pi)

        k = (A_strut * E_strut) / L_strut
        compression = 0.001 * N_landing
        R_s_spring = k * compression
        print("Spring R", R_s_spring)

        # Buckling
        I_strut = (np.pi*d_strut**4)/64
        P_cr = (np.pi**2 * E_strut * I_strut) / (L_strut**2)
        print("Buckling", P_cr)

        R_s = min(R_s_spring, P_cr)
        # R_s = 0  # for validation

        # Strut location
        x_s = b * strut_loc

        print("M_cantilever:", M_cantilever)
        # print("R_spring:", R_s_spring)
        # Final moment
        M_total = abs(M_cantilever) - R_s * x_s

        return abs(M_total), 0

    def max_shear_load(self, L, N_landing, strut_loc, strut_area):

        b = np.sqrt(aircraft.AR * aircraft.s_ref) / 2

        # weights
        W_struct = N_landing * self.weight_estimate
        W_fuel = N_landing * self.W_fuel_per_wing
        W_duct = N_landing * self.W_duct

        # --- SAME cantilever shear ---
        V_cantilever = W_struct + 4 * W_duct + W_fuel - L * N_landing

        # --- STRUT FORCE ---
        E_strut = 6.9e10
        L_strut = np.sqrt(1.1**2 + (b * strut_loc)**2)
        A_strut = strut_area
        d_strut = 2*np.sqrt(A_strut/np.pi)

        k = (A_strut * E_strut) / L_strut
        compression = 0.001 * N_landing
        #print(k, compression)
        R_s_spring = k * compression
        #print("R_s_spring:", R_s_spring)

        I_strut = (np.pi*d_strut**4)/64
        P_cr = (np.pi**2 * E_strut * I_strut) / (L_strut**2)
        print("Buckling", P_cr)

        R_s = min(R_s_spring, P_cr)
        # R_s = 0

        # Final shear
        print("V_cantilever:", V_cantilever)
        V_total = abs(V_cantilever) - R_s

        return abs(V_total), 0


    def get_N(self, L, a_z):
        """Gets flight load factor, landing load factor, and max load factor"""
        # Get required variable
        # W_total = self.aero["W_total"]

        # Calculate max load factor
        N = L/self.W_total # flight load factor
        N_land = a_z/self.g # landing load factor
        N_max = max(N, N_land) # maximum load factor

        return N, N_land, N_max

    def spar_cap_area(self, L, a_z, axial_stress, max_moment):
        """Find necessary spar cap area based on loading"""
        # Get required variables
        c = 2.49
        c_tip = c
        c_0 = c
        taper_ratio = c_tip/c_0
        N_max = self.get_N(L, a_z)[2]
        W_cent = self.W_total - self.weight_estimate
        b = aircraft.b
        t_x1 = self.aero["t_x1"]
        t_x2 = self.aero["t_x2"]
        x_1 = 0.1*c
        x_2 = 0.7*c
        # h = 2.64*.15
        h = 0.5*(t_x1+t_x2)
        N = self.get_N(L, a_z)[0]
        E = self.materials["spar_cap_E"] # maybe remain this variable if it shows up again just for clarity
        w_b_max = self.aero["w_b_max"] # not sure where we'll actually get this

        # Calculate area based on strength and stiffness requirements
        # A_cap_0_strength = (N_max*W_cent)/(12*axial_stress)*(b/h)*(1+2*taper_ratio)/(1+taper_ratio)
        # A_cap_0_stiffness = (N*W_cent)/(48*E) * b**2/h**2 * 1/w_b_max * (1+2*taper_ratio)/(1+taper_ratio)
        A_cap_0_strength = max_moment/(h*axial_stress)
        A_cap_0_stiffness = (max_moment*b**2)/(4*h**2*E*w_b_max)
        # return A_cap_0_strength
        return max(A_cap_0_strength, A_cap_0_stiffness)

    def spar_web_area(self, L, a_z, shear_stress, max_shear):
        """Find necessary spar web area based on loading"""
        # Get required variables
        N_max = self.get_N(L, a_z)[2]
        W_cent = self.W_total - self.weight_estimate

        # Return spar web area
        # return (N_max*W_cent)/(2*shear_stress)
        return max_shear/shear_stress

    def skin_thickness(self, q, shear_stress):
        """Find necessary skin thickness based on loading"""
        # Get required variables
        b = aircraft.b/2
        c = np.sqrt(aircraft.s_ref/aircraft.AR)

        # NOTE: all of these are guesses
        b_ail = self.aero["b_ail"]
        c_ail = self.aero["c_ail"]
        y_ail = self.aero["y_ail"]
        c_m = self.aero["c_m"]
        t_avg = self.aero["t_avg"] # Average thickness to chord of the airfoil
        s_tot = self.aero["s_tot"]
        # s_tot = self.aero["airfoil_perimeter"] # we'll estimate this as 2 * c_ail for now
        G = self.materials["skin_G"]
        twist_max = self.aero["twist_max"]
        A = self.aero["A"]
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

    def takeoff(self, A_strut):
        print("Takeoff Values:")
        # Find forces - do we need to incorporate some takeoff angle into this? maybe for all the stuff angle is an option and during cruise it's just zero?
        q_takeoff = 0.5*self.rho_ground*aircraft.v_takeoff**2
        L_takeoff = self.aero["CL_takeoff"]*q_takeoff*aircraft.s_ref # or lift distribution from aero
        T_takeoff = self.loading["T_takeoff"] # not sure where this will come from exactly
        D_takeoff = self.aero["CD_takeoff"]*q_takeoff*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_takeoff = self.axial_stress()
        shear_stress_takeoff = self.shear_stress()
        moment_takeoff_pre_strut = [self.max_moment(L_takeoff, 1, (1/5), A_strut)[0], self.max_moment(L_takeoff, 1, (1/4), A_strut)[0], self.max_moment(L_takeoff, 1, (1/3), A_strut)[0]]
        moment_takeoff_post_strut = [self.max_moment(L_takeoff, 1, (1/5), A_strut)[1], self.max_moment(L_takeoff, 1, (1/4), A_strut)[1], self.max_moment(L_takeoff, 1, (1/3), A_strut)[1]]
        shear_load_takeoff_pre_strut = [self.max_shear_load(L_takeoff, 1, (1/5), A_strut)[0], self.max_shear_load(L_takeoff, 1, (1/4), A_strut)[0], self.max_shear_load(L_takeoff, 1, (1/3), A_strut)[0]]
        shear_load_takeoff_post_strut = [self.max_shear_load(L_takeoff, 1, (1/5), A_strut)[1], self.max_shear_load(L_takeoff, 1, (1/4), A_strut)[1], self.max_shear_load(L_takeoff, 1, (1/3), A_strut)[1]]

        # Find component sizing based on calculated loading
        takeoff_spar_cap_area_pre_strut = [self.spar_cap_area(L_takeoff, 0, axial_stress_takeoff, moment_takeoff_pre_strut[0]), self.spar_cap_area(L_takeoff, 0, axial_stress_takeoff, moment_takeoff_pre_strut[1]), self.spar_cap_area(L_takeoff, 0, axial_stress_takeoff, moment_takeoff_pre_strut[2])]
        takeoff_spar_cap_area_post_strut = [self.spar_cap_area(L_takeoff, 0, axial_stress_takeoff, moment_takeoff_post_strut[0]), self.spar_cap_area(L_takeoff, 0, axial_stress_takeoff, moment_takeoff_post_strut[1]), self.spar_cap_area(L_takeoff, 0, axial_stress_takeoff, moment_takeoff_post_strut[2])]
        takeoff_spar_web_area_pre_strut = [self.spar_web_area(L_takeoff, 0, shear_stress_takeoff, shear_load_takeoff_pre_strut[0]), self.spar_web_area(L_takeoff, 0, shear_stress_takeoff, shear_load_takeoff_pre_strut[1]), self.spar_web_area(L_takeoff, 0, shear_stress_takeoff, shear_load_takeoff_pre_strut[2])]
        takeoff_spar_web_area_post_strut = [self.spar_web_area(L_takeoff, 0, shear_stress_takeoff, shear_load_takeoff_post_strut[0]), self.spar_web_area(L_takeoff, 0, shear_stress_takeoff, shear_load_takeoff_post_strut[1]), self.spar_web_area(L_takeoff, 0, shear_stress_takeoff, shear_load_takeoff_post_strut[2])]
        takeoff_skin_thickness = self.skin_thickness(q_takeoff, shear_stress_takeoff)

        return takeoff_spar_cap_area_pre_strut, takeoff_spar_cap_area_post_strut, takeoff_spar_web_area_pre_strut, takeoff_spar_web_area_post_strut, takeoff_skin_thickness


    def landing(self, A_strut):
        print("Landing Values:")
        # Find forces - do we need to incorporate landing angle into this? maybe for all the stuff angle is an option and during cruise it's just zero?
        q_landing = 0.5*self.rho_ground*aircraft.v_landing**2
        L_landing = self.aero["CL_landing"]*q_landing*aircraft.s_ref # or lift distribution from aero
        T_landing = self.loading["T_landing"] # not sure where this will come from exactly
        D_landing = self.aero["CD_landing"]*q_landing*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_landing = self.axial_stress()
        shear_stress_landing = self.shear_stress()
        moment_landing_pre_strut = [self.max_moment(L_landing, 1.8, (1/5), A_strut)[0], self.max_moment(L_landing, 1.8, (1/4), A_strut)[0], self.max_moment(L_landing, 1.8, (1/3), A_strut)[0]]
        moment_landing_post_strut = [self.max_moment(L_landing, 1.8, (1/5), A_strut)[1], self.max_moment(L_landing, 1.8, (1/4), A_strut)[1], self.max_moment(L_landing, 1.8, (1/3), A_strut)[1]]
        shear_load_landing_pre_strut = [self.max_shear_load(L_landing, 1.8, (1/5), A_strut)[0], self.max_shear_load(L_landing, 1.8, (1/4), A_strut)[0], self.max_shear_load(L_landing, 1.8, (1/3), A_strut)[0]]
        shear_load_landing_post_strut = [self.max_shear_load(L_landing, 1.8, (1/5), A_strut)[1], self.max_shear_load(L_landing, 1.8, (1/4), A_strut)[1], self.max_shear_load(L_landing, 1.8, (1/3), A_strut)[1]]

        # Find component sizing based on calculated loading
        landing_spar_cap_area_pre_strut = [self.spar_cap_area(L_landing, 0, axial_stress_landing, moment_landing_pre_strut[0]), self.spar_cap_area(L_landing, 0, axial_stress_landing, moment_landing_pre_strut[1]), self.spar_cap_area(L_landing, 0, axial_stress_landing, moment_landing_pre_strut[2])]
        landing_spar_cap_area_post_strut = [self.spar_cap_area(L_landing, 0, axial_stress_landing, moment_landing_post_strut[0]), self.spar_cap_area(L_landing, 0, axial_stress_landing, moment_landing_post_strut[1]), self.spar_cap_area(L_landing, 0, axial_stress_landing, moment_landing_post_strut[2])]
        landing_spar_web_area_pre_strut = [self.spar_web_area(L_landing, 0, shear_stress_landing, shear_load_landing_pre_strut[0]), self.spar_web_area(L_landing, 0, shear_stress_landing, shear_load_landing_pre_strut[1]), self.spar_web_area(L_landing, 0, shear_stress_landing, shear_load_landing_pre_strut[2])]
        landing_spar_web_area_post_strut = [self.spar_web_area(L_landing, 0, shear_stress_landing, shear_load_landing_post_strut[0]), self.spar_web_area(L_landing, 0, shear_stress_landing, shear_load_landing_post_strut[1]), self.spar_web_area(L_landing, 0, shear_stress_landing, shear_load_landing_post_strut[2])]
        landing_skin_thickness = self.skin_thickness(q_landing, shear_stress_landing)

        return landing_spar_cap_area_pre_strut, landing_spar_cap_area_post_strut, landing_spar_web_area_pre_strut, landing_spar_web_area_post_strut, landing_skin_thickness

    def max_load_sizing(self, A_strut):
        sizing_takeoff = self.takeoff(A_strut)
        sizing_landing = self.landing(A_strut)

        print("sizing_takeoff:", sizing_takeoff)
        print("sizing_landing:", sizing_landing)


        # TODO: make sure this works correctly
        spar_cap_area_pre_strut = [max(a, b) for a, b in zip(sizing_takeoff[0], sizing_landing[0])]
        spar_cap_area_post_strut = [max(a, b) for a, b in zip(sizing_takeoff[1], sizing_landing[1])]
        spar_web_area_pre_strut = [max(a, b) for a, b in zip(sizing_takeoff[2], sizing_landing[2])]
        spar_web_area_post_strut = [max(a, b) for a, b in zip(sizing_takeoff[3], sizing_landing[3])]
        # spar_cap_area_pre_strut = max(sizing_takeoff[0], sizing_landing[0])
        # spar_cap_area_post_strut = max(sizing_takeoff[1], sizing_landing[1])
        # spar_web_area_pre_strut = max(sizing_takeoff[2], sizing_landing[2])
        # spar_web_area_post_strut = max(sizing_takeoff[3], sizing_landing[3])
        skin_thickness = max(sizing_takeoff[4], sizing_landing[4])

        return spar_cap_area_pre_strut, spar_cap_area_post_strut, spar_web_area_pre_strut, spar_web_area_post_strut, skin_thickness


    def wing_weight(self, A_strut):
        b = aircraft.b
        sizing = self.max_load_sizing(A_strut)
        print("Sizing:", sizing)

        wing_weights = []

        skin_weight = sizing[4] * self.aero["airfoil_surface_area"] * 0.4 * self.materials["skin_density"] * 1.5

        for i in range(3):
            # individual components
            spar_cap_pre = sizing[0][i] * (b/2) * self.materials["spar_cap_density"] * 0.5
            # spar_cap_post = sizing[1][i] * (b/4) * self.materials["spar_cap_density"]

            spar_web_pre = sizing[2][i] * (b/2) * self.materials["spar_web_density"] * 0.5
            # spar_web_post = sizing[3][i] * (b/4) * self.materials["spar_web_density"]

            # total per configuration
            spar_and_skin = (
                2 * 1.5 * spar_cap_pre
                + 2 * 1.5 * spar_web_pre
                + skin_weight
            )

            wing_weight = spar_and_skin / 0.6
            wing_weights.append(wing_weight)

            print(f"\nOption {i+1}:")
            print("  spar caps:", spar_cap_pre)
            print("  spar cap pre:", spar_cap_pre)
            # print("  spar cap post:", spar_cap_post)
            print("  spar webs:", spar_web_pre)
            print("  skin weight:", skin_weight)
            print("  total wing weight:", wing_weight)

        return wing_weights


if __name__ == "__main__":
# Create dictionaries for testing
    c = np.sqrt(aircraft.s_ref/aircraft.AR)
    aero = {
        "t_x1": (0.7279684E-01 + 0.2596068E-01)*c,
        "t_x2": (0.2139701E-01 + 0.4671981E-01)*c,
        # all below values are guesses
        "w_b_max": 19.92/2*1.05,
        "twist_max": 0.0523599, # radians
        "b_ail":2.5, #should this for one aileron or both? this is for one aileron, so total aileron span is 5 m
        "c_ail":0.75,
        "y_ail":7.47 + 2.5/2,
        "c_m":-0.2,
        "A":.275*1.6,
        "s_tot":2.0700619*2.49,
        "x_1":.06*2.49,
        "x_2":(.06+.6)*2.49,
        "t_x1":.275,
        "t_x2":.275,
        "CL_takeoff":6.1,
        "CL_climb":.426, #but depends where in climb we are
        "CL_cruise":.237,
        "CL_descent":.4, #tbh we don't really know yet
        "CL_landing":6.1, #also subject to significant change
        "CD_takeoff":1.665,#including induced drag
        "CD_climb":0.027,#including induced drag
        "CD_cruise":.019,
        "CD_descent":.04,
        "CD_landing":1.665,
        "v_takeoff":20,
        "v_climb":80,
        "v_cruise":125,
        "v_descent":80,
        "v_landing":20,
        "drag_area": 2.5,
        "a_z": 2*9.8,
        "airfoil_surface_area":2.0700619*2.49,
        "t_avg": 0.12,
    }

    loading = {
        "T_takeoff": 22.3 * 10**3,
        "T_climb": 18645,
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

    test_wing = Wing(aero, loading, materials, 1192.42*9.8)

    test_wing_weight = test_wing.wing_weight(1e-3)
    print("Wing weight (kg)", test_wing_weight)

    # # =========================
    # # PLOTTING SECTION
    # # =========================

    # import matplotlib.pyplot as plt
    # import numpy as np

    # # --- Sweep over strut area ---
    # strut_areas = np.linspace(1e-4, 5e-3, 25)  # m^2

    # # Store results for 3 configurations (1/5, 1/4, 1/3)
    # weights_vs_area = [[], [], []]

    # # Run sweep
    # for A in strut_areas:
    #     weights = test_wing.wing_weight(A)  # returns [w1, w2, w3]
    #     for i in range(3):
    #         weights_vs_area[i].append(weights[i])


    # # =========================
    # # Plot 1: Wing Weight vs Strut Area
    # # =========================
    # plt.figure(figsize=(8,5))

    # labels = ['Strut @ 1/5 span', 'Strut @ 1/4 span', 'Strut @ 1/3 span']

    # for i in range(3):
    #     plt.plot(strut_areas, weights_vs_area[i], marker='o', label=labels[i])

    # plt.xlabel("Strut Area (m²)")
    # plt.ylabel("Wing Weight (kg)")
    # plt.title("Wing Weight vs Strut Area")
    # plt.legend()
    # plt.grid(True)
    # plt.show()


    # # =========================
    # # Plot 2: Wing Weight vs Strut Location
    # # =========================
    # A_ref = 1e-3  # choose a representative area
    # wing_weights_locations = test_wing.wing_weight(A_ref)

    # plt.figure(figsize=(8,5))
    # plt.plot(['1/5', '1/4', '1/3'], wing_weights_locations, marker='o')

    # plt.xlabel("Strut Location (fraction of span)")
    # plt.ylabel("Wing Weight (kg)")
    # plt.title(f"Wing Weight vs Strut Location (A_strut = {A_ref:.1e} m²)")
    # plt.grid(True)
    # plt.show()


    # # =========================
    # # Plot 3: With vs Without Strut
    # # =========================

    # # Best case WITH strut (minimum across all sweeps)
    # wing_weight_with_strut = min([min(w) for w in weights_vs_area])

    # # Approximate WITHOUT strut → use near-zero area
    # wing_weight_no_strut = min(test_wing.wing_weight(1e-6))

    # plt.figure(figsize=(8,5))
    # plt.bar(['With Strut (optimized)', 'Without Strut'],
    #         [wing_weight_with_strut, wing_weight_no_strut])

    # plt.ylabel("Wing Weight (kg)")
    # plt.title("Wing Weight: With vs Without Strut")
    # plt.show()


    # # def max_moment(self, L, N_landing, strut_loc, strut_area):

    # #     # Geometry
    # #     b = np.sqrt(aircraft.AR*aircraft.s_ref)
    # #     c = np.sqrt(aircraft.s_ref/aircraft.AR)

    # #     x_T1 = 2.78
    # #     x_T2 = 3.46
    # #     x_T3 = 6.14
    # #     x_T4 = 7.82

    # #     # weights
    # #     W_struct = N_landing*self.weight_estimate
    # #     W_fuel = N_landing*self.W_fuel_per_wing
    # #     W_struct_tot = (W_struct+W_fuel)*2
    # #     W_duct = N_landing*self.W_duct

    # #     # Strut spring force
    # #     E_strut = 6.9e10  # <-- FIX THIS (you had 6.9e9, aluminum is 69 GPa)
    # #     L_strut = np.sqrt(1.1**2 + (b*strut_loc)**2)
    # #     A_strut = strut_area

    # #     # Spring stiffness
    # #     k = (A_strut * E_strut) / L_strut
    # #     compression = 0.01 * N_landing
    # #     R_s_spring = k * compression

    # #     # --- Euler buckling limit ---
    # #     K = 1.0  # pinned-pinned

    # #     # second moment of area (circular assumption)
    # #     I_strut = (A_strut**2) / (4 * np.pi)

    # #     P_cr = (np.pi**2 * E_strut * I_strut) / ((K * L_strut)**2)

    # #     # --- FINAL STRUT FORCE ---
    # #     R_s = min(R_s_spring, P_cr)
    # #     # R_s = 0

    # #     # print("k", k)
    # #     # print("Spring force", R_s_spring)
    # #     # print("Buckling limit", P_cr)
    # #     # print("Final strut force", R_s)

    # #     if strut_loc == (1/5):
    # #         M_z_1 = abs(-0.125*b*W_struct_tot-2*W_duct*(x_T1+x_T2) + 0.125*L) - (b/5)*R_s
    # #         # print("M_1/5 breakdown", -0.125*b*W_struct_tot, -2*W_duct*(x_T1+x_T2), 0.125*L, - (b/5)*R_s)
    # #         M_z_2 = -0.1795*b*W_struct_tot - 2*W_duct*(x_T3 + x_T4) + 0.1795*L

    # #     elif strut_loc == (1/4):
    # #         M_z_1 = abs(-0.1582*b*W_struct_tot-2*W_duct*(x_T1+x_T2) + 0.1582*L) - (b/4)*R_s
    # #         M_z_2 = -0.1463*b*W_struct_tot - 2*W_duct*(x_T3 + x_T4) + 0.1463*L

    # #     elif strut_loc == (1/3):
    # #         M_z_1 = abs(-0.2113*b*W_struct_tot-3*W_duct*(x_T1+x_T2+x_T3) + 0.2113*L) - (b/3)*R_s
    # #         M_z_2 = -0.0931*b*W_struct_tot - W_duct*(x_T4) + 0.0931*L

    # #     # print("Moment Pre Strut:", abs(M_z_1))
    # #     # print("Moment Post Strut:", abs(M_z_2))
    # #     return abs(M_z_1), abs(M_z_2)

    # # def max_shear_load(self, L, N_landing, strut_loc, strut_area):

    # #     b = np.sqrt(aircraft.AR*aircraft.s_ref)
    # #     c = np.sqrt(aircraft.s_ref/aircraft.AR)

    # #     # weights
    # #     W_struct = N_landing*self.weight_estimate
    # #     W_fuel = N_landing*self.W_fuel_per_wing
    # #     W_struct_tot = (W_struct+W_fuel)*2
    # #     W_duct = N_landing*self.W_duct

    # #     # Strut spring force
    # #     E_strut = 6.9e10  # <-- FIX THIS (you had 6.9e9, aluminum is 69 GPa)
    # #     L_strut = np.sqrt(1.1**2 + (b*strut_loc)**2)
    # #     A_strut = strut_area

    # #     # Spring stiffness
    # #     k = (A_strut * E_strut) / L_strut
    # #     compression = 0.01 * N_landing
    # #     R_s_spring = k * compression

    # #     # --- Euler buckling limit ---
    # #     K = 1.0  # pinned-pinned

    # #     # second moment of area (circular assumption)
    # #     I_strut = (A_strut**2) / (4 * np.pi)

    # #     P_cr = (np.pi**2 * E_strut * I_strut) / ((K * L_strut)**2)

    # #     # --- FINAL STRUT FORCE ---
    # #     R_s = min(R_s_spring, P_cr)
    # #     # R_s = 0

    # #     # print("k", k)
    # #     # print("Spring force", R_s_spring)
    # #     # print("Buckling limit", P_cr)
    # #     # print("Final strut force", R_s)

    # #     if strut_loc == (1/5):
    # #         F_y_1 = abs(-0.086*W_struct_tot-2*W_duct + 0.086*L) - R_s
    # #         F_y_2 = -0.1*W_struct_tot - 2*W_duct + 0.1*L

    # #     elif strut_loc == (1/4):
    # #         F_y_1 = abs(-0.1071*W_struct_tot-2*W_duct + 0.1071*L) - R_s
    # #         F_y_2 = -0.079*W_struct_tot - 2*W_duct + 0.079*L

    # #     elif strut_loc == (1/3):
    # #         F_y_1 = abs(-0.1391*W_struct_tot-3*W_duct + 0.1391*L) - R_s
    # #         F_y_2 = -0.047*W_struct_tot - W_duct + 0.047*L

    # #     # print("Force Pre Strut:", abs(F_y_1))
    # #     # print("Force Post Strut:", abs(F_y_2))

    # #     return abs(F_y_1), abs(F_y_2)
