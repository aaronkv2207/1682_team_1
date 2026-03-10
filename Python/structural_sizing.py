# Team 1 - Structural Sizing
import numpy as np
import math as math

from aero_workspace.aero_main import AircraftConfig as aircraft
from aero_workspace.aero_main import TakeoffCoeff as takeoff

# TODO: go through and check for b/2, L/2, etc. stuff (depending on how we define b, L, etc. in the code)
# check signs in stress equations

class Wing():
    """Wing class:
    - Takes in general wing sizing and material
    - Sizes structural components on the wing
    - Outputs weight of the wing"""

    # Set constants
    g = 9.81
    rho_ground = 1.225
    rho_cruise = 0.699 # at 18,000 ft
    W_total = 7438.9149
    W_duct = 30 # kg - guess

    def __init__(self, aero, loading, materials, weight_estimate = 0.07*W_total):
        self.aero = aero
        self.loading = loading
        self.materials = materials
        self.weight_estimate = weight_estimate

    def axial_stress(self, L, T, D):
        # Get required variables
        b = np.sqrt(aircraft.AR*aircraft.s_ref)
        c = np.sqrt(aircraft.S/aircraft.AR)
        x_T1 = b/5
        x_T2 = b*(2/5)
        x_T3 = b*(3/5)
        x_T4 = b*(4/5)
        t_x1 = self.aero["t_x1"]
        t_x2 = self.aero["t_x2"]
        x_1 = 0.1*c
        x_2 = 0.7*c

        h = 0.5*(t_x1+t_x2)
        w = x_2-x_1
        t = ... # box thickness - estimate

        # NOTE: not sure if this is right z and y are dist. from the neutral axis
        z = w/2
        y = h/2


        # # Moment equations for final model:
        # M_y = T*(x_T1+x_T2+x_T3+x_T4) # - drag integral
        # M_z = self.weight_estimate*(b/4) + W_duct*(x_T1+x_T2+x_T3+x_T4) # - lift integral
        # # TODO: add actual integral stuff for D, L, W

        # Moment equations as of now (3/8) to test model, assume uniform distribution:
        M_y = T*(x_T1+x_T2+x_T3+x_T4) - (D*(b/2)**2)/2
        M_z = self.weight_estimate*(b/4) + self.W_duct*(x_T1+x_T2+x_T3+x_T4) - (L*(b/2)**2)/2

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
        b = np.sqrt(aircraft.AR*aircraft.s_ref)
        c = np.sqrt(aircraft.S/aircraft.AR)
        x_T1 = b/5
        x_T2 = b*(2/5)
        x_T3 = b*(3/5)
        x_T4 = b*(4/5)
        t_x1 = self.aero["t_x1"]
        t_x2 = self.aero["t_x2"]
        x_1 = 0.1*c
        x_2 = 0.7*c

        h = 0.5*(t_x1+t_x2)
        w = x_2-x_1
        t = ... # box thickness - estimate

        # Force equations
        F_y = self.weight_estimate*(b/2) + 4*self.W_duct - L # - lift integral
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
        N_max = np.max(N, N_land) # maximum load factor

        return N, N_land, N_max

    def spar_cap_area(self, L, a_z, axial_stress):
        """Find necessary spar cap area based on loading"""
        # Get required variables
        c = np.sqrt(aircraft.S/aircraft.AR)
        c_tip = c
        c_0 = c
        taper_ratio = c_tip/c_0
        N_max = self.get_N(L, a_z)[2]
        W_cent = self.W_total - self.weight_estimate
        b = np.sqrt(aircraft.AR*aircraft.s_ref)
        t_x1 = self.aero["t_x1"]
        t_x2 = self.aero["t_x2"]
        x_1 = 0.1*c
        x_2 = 0.7*c
        h = 0.5*(t_x1+t_x2)
        N = self.get_N(L, a_z)[0]
        E = self.materials["spar_car_E"] # maybe remain this variable if it shows up again just for clarity
        w_b_max = self.aero["w_b_max"] # not sure where we'll actually get this

        # Calculate area based on strength and stiffness requirements
        A_cap_0_strength = (N_max*W_cent)/(12*axial_stress)*(b/h)*(1+2*taper_ratio)/(1+taper_ratio)
        A_cap_0_stiffness = (N*W_cent)/(48*E) * b**2/h**2 * 1/w_b_max * (1+2*taper_ratio)/(1+taper_ratio)
        return np.max(A_cap_0_strength, A_cap_0_stiffness)

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
        b_ail = self.aero["b_ail"]
        c_ail = self.aero["c_ail"]
        y_ail = self.aero["y_ail"]
        c_m = self.aero["c_m"] # I don't think this changes based on flight conditions but not sure
        x1 = self.aero["x1"]
        x2= self.aero["x2"]
        tx1 = self.aero["tx1"]
        tx2 = self.aero["tx2"]
        s_tot = self.aero["s_tot"] # not sure what this is
        G = self.materials["skin_G"]
        twist_max = self.aero["twist_max"]
        A = self.aero["airfoil_surface_area"]


        # Calculate skin thickness requirements
        t_skin_strength = q*b_ail*c_ail**2*c_m*1/(2*A)*1/(1/shear_stress)
        t_skin_stiffness = 1*b_ail*c_ail**2*y_ail*c_m*s_tot/(4*A**2*G)*(1/twist_max)

        # Return max thickness
        return np.max(t_skin_strength, t_skin_stiffness)

    def tube_thickness(self):
        # TODO: implement this (not sure where this component goes so not doing this yet)
        pass

    def takeoff(self):
        # Find forces - do we need to incorporate some takeoff angle into this? maybe for all the stuff angle is an option and during cruise it's just zero?
        q_takeoff = 0.5*self.rho_ground*aircraft.v_takeoff**2
        L_takeoff = self.aero["CL_takeoff"]*q_takeoff*aircraft.S # or lift distribution from aero
        T_takeoff = self.loading["T_takeoff"] # not sure where this will come from exactly
        D_takeoff = self.aero["CD_takeoff"]*q_takeoff*self.aero["drag_area"]


        # Calculate stresses and torsion
        axial_stress_takeoff = self.axial_stress(L_takeoff, T_takeoff, D_takeoff)
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
        L_climb = self.aero["CL_climb"]*q_climb*aircraft.S # or lift distribution from aero
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
        q_cruise = 0.5*self.rho_cruise*aircraft.v_cruise**2
        L_cruise = self.aero["CL_cruise"]*q_cruise*aircraft.S # or lift distribution from aero
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
        L_descent = self.aero["CL_descent"]*q_descent*aircraft.S # or lift distribution from aero
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
        L_landing = self.aero["CL_landing"]*q_landing*aircraft.S # or lift distribution from aero
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






class Fuselage:

    def __init__(self, length, radius, n):
        self.length = length
        self.R = radius
        self.n = n
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

    def required_thickness_moment(self, yield_stress, safety_factor=2):
        pass

    def get_dead_weight(self):
        seat_weight = self.n*13.0 # kg (average modern aircraft seat weight)
        person_weight = self.n*100.0 # kg
        return seat_weight + person_weight


    def get_structural_weight(self, safety_factor):
        # assume some materials

        # skin material aluminum 7075-T6
        yield_strength = 490e6    # Pa
        skin_rho = 2810    # kg/m^3
        # skin weight
        self.thickness = self.required_thickness_hoop(yield_strength)+self.required_thickness_moments(yield_stress)
        self.skin_volume = self.length*2*math.pi*self.R*self.thickness
        self.skin_weight = self.skin_volume*self.skin_rho

        # frame weight
        self.frame_weight

        # stringer weight
        self.stringer_weight

        return self.skin_weight+self.frame_weight+self.stringer_weight


    def get_total_weight(self, n, safety_factor):
        g  = 9.80665     # Gravity (m/s^2)
        self.dead_weight = self.get_dead_weight()
        self.structural_weight = self.get_structural_weight(safety_factor)
        return self.dead_weight*g + self.structural_weight*g



class LandingGear():
    pass


if __name__ == "__main__":
    # Create dictionaries for testing
    # NOTE: h here is spar height which is maybe something we calculate ourselves
    aero = {
        "b_ail":,
        "c_ail":,
        "c_m":,
        "A":,
        "y_ail":,
        "s_tot":,
        "t_x1": 0.7279684E-01 + 0.2596068E-01,
        "t_x2": 0.2139701E-01 + 0.4671981E-01,
        "W_total":,
        "w_b_max":,
        "twist_max":,
        "airfoil_surface_area":,
        "CL_takeoff":,
        "CL_climb":,
        "CL_cruise":,
        "CL_descent":,
        "CL_landing":,
        "CD_takeoff":,
        "CD_climb":,
        "CD_cruise":,
        "CD_descent":,
        "CD_landing":,
        "v_takeoff":,
        "v_climb":,
        "v_cruise":,
        "v_descent":,
        "v_landing":,
        "drag_area":,
        "a_z":,
        "b_vert":,
        "b_ho":,
        "t_x1_v":,
        "t_x2_v":,
        "x_1_v":,
        "x_2_v":,
        "t_x1_h":,
        "t_x2_h":,
        "x_1_h":,
        "x_2_h":,
    }

    loading = {
        "T_takeoff":,
        "T_climb":,
        "T_cruise":,
        "T_descent":,
        "T_landing":,
    }

    materials = {
        "spar_cap_E": 6.9 * 10**9,
        "skin_G": 26 * 10**9
    }

    # wing weight estimate
    # weight_estimate = 0.07*aero["W_total"]

    test_wing = Wing(aero, loading, materials)


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
