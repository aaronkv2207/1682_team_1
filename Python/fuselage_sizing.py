import numpy as np
import math
import matplotlib.pyplot as plt


class Fuselage:

    def __init__(self, length, radius, n_people, cabin_width):
        self.length = length
        self.R = radius
        self.n = n_people
        self.cabin_width = cabin_width
        self.weight = 0.0

        # sized values (filled in later)
        self.skin_t = None
        self.stringer_area = None
        self.n_stringers = None

    # -------------------------
    # Atmosphere
    # -------------------------
    def pressure_at_altitude(self, h):
        """
        Returns atmospheric pressure (Pa) at altitude h (meters)
        Valid up to 11 km (troposphere)
        """
        P0 = 101325
        T0 = 288.15
        L = 0.0065
        g = 9.80665
        R_air = 287.05

        return P0 * (1 - (L * h) / T0) ** (g / (R_air * L))

    # -------------------------
    # Geometry functions
    # -------------------------
    def circumference(self):
        return 2.0 * np.pi * self.R

    def shell_area(self):
        return self.circumference() * self.length

    def enclosed_area(self):
        """Area enclosed by circular skin midline approximation"""
        return np.pi * self.R**2

    def count_members(self, length, spacing):
        return int(np.floor(length / spacing)) + 1

    def stringer_angles(self, n_stringers):
        """Evenly spaced angles around circle"""
        return np.linspace(0.0, 2.0 * np.pi, n_stringers, endpoint=False)

    def stringer_y_coords(self, n_stringers):
        """
        y-coordinates of longerons for vertical bending.
        Top/bottom longerons have largest |y|.
        """
        theta = self.stringer_angles(n_stringers)
        return self.R * np.sin(theta)

    # -------------------------
    # Section inertias
    # -------------------------
    def skin_bending_inertia(self, skin_thickness):
        """
        Thin-walled circular shell bending inertia about horizontal centroidal axis:
            I_skin = pi * R^3 * t
        """
        return np.pi * self.R**3 * skin_thickness

    def stringer_bending_inertia(self, stringer_area, n_stringers):
        """
        Boom/stringer contribution:
            I_stringers = sum(A_s * y_i^2)
        """
        y = self.stringer_y_coords(n_stringers)
        return stringer_area * np.sum(y**2)

    # -------------------------
    # Allowables
    # -------------------------
    def allowable_normal_stress(self, yield_strength, safety_factor=1.5):
        return yield_strength / safety_factor

    def allowable_shear_stress(self, yield_strength, safety_factor=1.5, shear_ratio=0.577):
        """
        Approx shear yield from von Mises: tau_y ~ 0.577*sigma_y
        """
        return shear_ratio * yield_strength / safety_factor

    # -------------------------
    # Pressure sizing
    # -------------------------
    def required_thickness_hoop(self, yield_strength, safety_factor=1.5,
                                max_altitude=5500.0, cabin_pressure=75150.0):
        """
        Thin-wall hoop-stress sizing of skin.
        """
        P_out = self.pressure_at_altitude(max_altitude)
        delta_P = max(cabin_pressure - P_out, 0.0)

        sigma_allow = self.allowable_normal_stress(yield_strength, safety_factor)

        t_hoop = (delta_P * self.R) / sigma_allow
        return t_hoop
    
    # -------------------------
    # Longeron sizing from bending
    # -------------------------
    def required_stringer_area_bending(self, M_max, skin_thickness, yield_strength,
                                       n_stringers=10, safety_factor=1.5,
                                       min_stringer_area=1e-5):
        """
        Size equal-area stringers with skin included in global bending stiffness.

        Total section inertia:
            I_total = I_skin + I_stringers
                    = pi*R^3*t + A_s * sum(y_i^2)

        Require:
            sigma_max = M*R / I_total <= sigma_allow

        Solve for A_s:
            A_s = (M*R/sigma_allow - I_skin) / sum(y_i^2)

        If skin alone is sufficient, return minimum practical stringer area.
        """
        sigma_allow = self.allowable_normal_stress(yield_strength, safety_factor)

        y = self.stringer_y_coords(n_stringers)
        y_max = np.max(np.abs(y))
        sum_y2 = np.sum(y**2)

        if sum_y2 <= 0.0:
            raise ValueError("Invalid stringer geometry for bending.")

        I_skin = self.skin_bending_inertia(skin_thickness)
        I_required = abs(M_max) * y_max / sigma_allow

        A_s = (I_required - I_skin) / sum_y2

        # If skin alone satisfies global bending, keep a minimum practical area
        A_s = max(A_s, min_stringer_area)

        return A_s

    # -------------------------
    # Skin sizing from shear
    # -------------------------
    def required_skin_thickness_shear(self, V_max, yield_strength, safety_factor=1.5):
        """
        First-pass max shear-flow approximation for circular closed section:
            q_max ~ V / (pi*R)
            tau = q/t  =>  t = q/tau_allow
        """
        tau_allow = self.allowable_shear_stress(yield_strength, safety_factor)

        q_max = abs(V_max) / (np.pi * self.R)
        t_shear = q_max / tau_allow

        return t_shear

    # -------------------------
    # Skin sizing from torsion
    # -------------------------
    def required_skin_thickness_torsion(self, T_max, yield_strength,
                                        safety_factor=1.5):
        """
        Closed single-cell thin-wall torsion:
            q = T / (2*A_m)
            tau = q/t  =>  t = q/tau_allow
        """
        tau_allow = self.allowable_shear_stress(yield_strength, safety_factor)
        A_m = self.enclosed_area()

        q_t = abs(T_max) / (2.0 * A_m)
        t_torsion = q_t / tau_allow

        return t_torsion

    # -------------------------
    # Overall structural sizing
    # -------------------------
    def size_primary_structure(self,
                               M_max,
                               V_max,
                               T_max=0.0,
                               yield_strength=400e6,
                               density=2700,
                               n_stringers=10,
                               safety_factor=1.5,
                               min_skin_gauge=0.0010):
        """
        Primary first-pass sizing:
        - longerons sized by bending
        - skin sized by max of hoop, shear, torsion, minimum gauge
        """
        self.n_stringers = n_stringers
        self.rho = density
        self.yield_strength = yield_strength
        self.safety_factor = safety_factor

        # skin
        t_hoop = self.required_thickness_hoop(
            yield_strength=yield_strength,
            safety_factor=safety_factor
        )
        t_shear = self.required_skin_thickness_shear(
            V_max=V_max,
            yield_strength=yield_strength,
            safety_factor=safety_factor
        )
        t_torsion = self.required_skin_thickness_torsion(
            T_max=T_max,
            yield_strength=yield_strength,
            safety_factor=safety_factor
        )

        self.skin_t = max(t_hoop, t_shear, t_torsion, min_skin_gauge)

        # keep individual drivers for debugging
        self.t_hoop = t_hoop
        self.t_shear = t_shear
        self.t_torsion = t_torsion
        self.t_min = min_skin_gauge

        # longerons
        self.stringer_area = self.required_stringer_area_bending(
            M_max=M_max, 
            skin_thickness=self.skin_t,
            yield_strength=yield_strength,
            n_stringers=n_stringers,
            safety_factor=safety_factor
        )

        return {
            "stringer_area_each_m2": self.stringer_area,
            "skin_thickness_m": self.skin_t,
            "t_hoop_m": self.t_hoop,
            "t_shear_m": self.t_shear,
            "t_torsion_m": self.t_torsion,
            "t_min_gauge_m": self.t_min
        }

    # -------------------------
    # Component masses
    # -------------------------
    def get_structural_mass(self,
                            frame_spacing=0.5,
                            frame_depth=0.10,
                            frame_thickness=0.002,
                            floor_thickness=0.01,
                            rho_floor=3.0,
                            floor_beam_depth=0.07,
                            floor_beam_thickness=0.002,
                            floor_beam_spacing=0.5):
        """
        Uses previously sized skin_t and stringer_area.
        """
        if self.skin_t is None or self.stringer_area is None or self.n_stringers is None:
            raise ValueError("Run size_primary_structure() before get_structural_mass().")

        rho = self.rho

        # skin mass
        skin_volume = self.shell_area() * self.skin_t
        self.skin_mass = skin_volume * rho

        # fuselage frames
        n_frames = self.count_members(self.length, frame_spacing)
        frame_ring_area = np.pi * self.R**2 - np.pi * (self.R - frame_depth)**2
        frame_volume = n_frames * frame_thickness * frame_ring_area
        self.frame_mass = frame_volume * rho

        # longeron mass
        stringer_volume = self.n_stringers * self.length * self.stringer_area
        self.stringer_mass = stringer_volume * rho

        # floor panel mass
        floor_area = self.length * self.cabin_width
        self.fpanel_mass = floor_area * floor_thickness * rho_floor

        # floor transverse members at fuselage frames
        n_fframes = n_frames
        self.fframe_mass = n_fframes * self.cabin_width * frame_thickness * floor_beam_depth * rho

        # floor longitudinal beams
        n_fbeams = self.count_members(self.cabin_width, floor_beam_spacing)
        self.fbeam_mass = n_fbeams * self.length * floor_beam_thickness * floor_beam_depth * rho

        self.structural_mass = (
            self.skin_mass
            + self.frame_mass
            + self.stringer_mass
            + self.fpanel_mass
            + self.fframe_mass
            + self.fbeam_mass
        )

        return self.structural_mass

    def get_dead_mass(self):
        seat_weight = self.n * 13.0
        person_weight = self.n * 100.0
        return seat_weight + person_weight

    def get_total_weight(self):
        g = 9.80665
        return (self.get_dead_mass() + self.structural_mass) * g

    # -------------------------
    # Nice summary print
    # -------------------------
    def summary(self):
        print("---- Fuselage Summary ----")
        if self.stringer_area is not None:
            print(f"Stringer area each      = {self.stringer_area:.6e} m^2")
        if self.skin_t is not None:
            print(f"Skin thickness          = {self.skin_t:.6e} m")
            print(f"  Hoop-driven t         = {self.t_hoop:.6e} m")
            print(f"  Shear-driven t        = {self.t_shear:.6e} m")
            print(f"  Torsion-driven t      = {self.t_torsion:.6e} m")
            print(f"  Minimum gauge         = {self.t_min:.6e} m")
        if hasattr(self, "structural_mass"):
            print(f"Structural mass         = {self.structural_mass:.2f} kg")
        if hasattr(self, "skin_mass"):
            print(f"  Skin mass             = {self.skin_mass:.2f} kg")
            print(f"  Frame mass            = {self.frame_mass:.2f} kg")
            print(f"  Stringer mass         = {self.stringer_mass:.2f} kg")
            print(f"  Floor panel mass      = {self.fpanel_mass:.2f} kg")
            print(f"  Floor frame mass      = {self.fframe_mass:.2f} kg")
            print(f"  Floor beam mass       = {self.fbeam_mass:.2f} kg")



# ============================================================
# Test script
# ============================================================

if __name__ == "__main__":

    # -------------------------
    # Aircraft geometry
    # -------------------------

    length = 9          # m  
    radius = 1.1        # m  
    n_people = 19        # 
    cabin_width = 1.25    # m usable cabin width

    fuse = Fuselage(length, radius, n_people, cabin_width)

    # -------------------------
    # Structural loads
    # (placeholder estimates — replace with your team's values)
    # -------------------------

    M_max = 1e6   # N*m  max fuselage bending moment
    V_max = 1e5     # N    max shear force
    h_Tforce = 1.6   # m  height from neutral axis where rudder force applied
    T_max = V_max*h_Tforce     # N*m  torsion

    # -------------------------
    # Size structure
    # -------------------------

    sizing_results = fuse.size_primary_structure(
        M_max=M_max,
        V_max=V_max,
        T_max=T_max
    )

    print("\n--- Structural sizing ---")
    for key, value in sizing_results.items():
        print(f"{key}: {value}")

    # -------------------------
    # Compute structural mass
    # -------------------------

    mass = fuse.get_structural_mass()

    print("\nStructural mass:", mass, "kg")

    # -------------------------
    # Total loaded fuselage mass
    # -------------------------

    total_weight = fuse.get_total_weight()
    total_mass = total_weight / 9.80665

    # -------------------------
    # Add landing gear
    # -------------------------

    landing_gear_mass = 102.43   # kg
    total_mass_with_lg = total_mass + landing_gear_mass

    print("Passenger + seat mass:", n_people*113, "kg")
    print("Landing gear mass:", landing_gear_mass, "kg")
    print("Total fuselage mass (including passengers and landing gear):", total_mass_with_lg, "kg")

    # -------------------------
    # Detailed summary
    # -------------------------

    print()
    fuse.summary()