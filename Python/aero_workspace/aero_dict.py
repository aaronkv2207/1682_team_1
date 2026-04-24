"""First plane uses a slightly smaller wing area. Performance for each plane may be evaluated, but start with plane 1 (AircraftConfig)"""

from dataclasses import dataclass

from ambiance import Atmosphere
from conceptual_design import ureg
from drag import (
    calc_C_Dp,
)  # TODO: Brenda's drag model needs to be updated based on revised plane geometry


# TODO: update constants in dataclass
@dataclass
class AircraftConfig:
    """All V2 aircraft definitons --> all in SI units"""

    # NOTE: Required inputs --> update when calling class!!!
    v_cruise: float = 0.0
    v_t0: float = 0.0

    # V2 Plane Geometry
    AR: float = 7.0
    s_ref: float = 45.0
    V_h: float = 1.0
    V_v: float = 0.06

    ail_hinge: float = 0.7
    flap_hinge: float = 0.75
    wing_incidence: float = 0.0
    aileron_fraction: float = 0.75

    lv: float = 8.0
    vt_ar: float = 1.2
    tail_hinge: float = 0.7
    ht_ar: float = 3.0

    fan_radius: float = 0.583
    n_fans: int = 8

    b: float = 17.75

    vt_area: float = 5.99
    vt_chord: float = 2.23
    mac: float = 2.54
    blown_span: float = 11.71
    fan_length: float = 9.33
    hdisk: float = 0.73
    blowing_dy: float = 1.46
    ht_area: float = 13.33
    ht_mac: float = 2.11

    #### TAKEOFF DEFINITIONS ###
    # NOTE: See post-processor
    h_t0: float = 0.0
    rho_t0 = Atmosphere(h=h_t0).density[0]  # [kg/m^3]
    mu_t0: float = Atmosphere(h=h_t0).dynamic_viscosity[0]

    # #### CLIMB DEFINITIONS (at a single point) ###
    # NOTE: See post-processor

    # #### CRUISE DEFINITIONS ###
    # NOTE: See post-processor
    h_cruise: int = 18000 * ureg("ft").to("m").magnitude

    def __post_init__(self):
        # Takeoff
        atm_t0 = Atmosphere(h=self.h_t0)
        self.rho_t0 = atm_t0.density[0]
        self.mu_t0 = atm_t0.dynamic_viscosity[0]

        self.Dp_t0, self.C_Dp_t0 = calc_C_Dp(
            self.rho_t0,
            self.v_t0,
            self.mu_t0,
        )

        # Cruise
        atm_cruise = Atmosphere(h=self.h_cruise)
        self.rho_cruise = atm_cruise.density[0]
        self.mu_cruise = atm_cruise.dynamic_viscosity[0]

        self.Dp_cruise, self.C_Dp_cruise = calc_C_Dp(
            self.rho_cruise,
            self.v_cruise,
            self.mu_cruise,
        )

    # #### LANDING DEFINITIONS ###
    # NOTE: See post-processor


@dataclass
class AircraftConfig2:
    """All V2 aircraft definitons w/ smaller wing area --> all in SI units"""

    # NOTE: Required inputs --> update when calling class!!!
    v_cruise: float = 0.0
    v_t0: float = 0.0

    # V2 Plane Geometry
    AR: float = 7.0
    s_ref: float = 42.0
    V_h: float = 1.0
    V_v: float = 0.06

    nose_x: float = -5.0
    fuse_width: float = 1.6
    ail_hinge: float = 0.7
    flap_hinge: float = 0.75
    wing_incidence: float = 0.0
    aileron_fraction: float = 0.75

    lv: float = 8.0
    vt_ar: float = 1.2
    tail_hinge: float = 0.7
    ht_ar: float = 3.0

    fan_radius: float = 0.583
    n_fans: int = 8

    c: float = 1.98
    b: float = 17.15

    vt_area: float = 5.40
    vt_chord: float = 2.12
    mac: float = 2.45
    blown_span: float = 11.26
    fan_length: float = 9.33
    hdisk: float = 0.76
    blowing_dy: float = 1.41
    ht_area: float = 12.06
    ht_mac: float = 2.01

    #### TAKEOFF DEFINITIONS ###
    # NOTE: See post-processor
    h_t0: float = 0.0
    rho_t0 = Atmosphere(h=h_t0).density[0]  # [kg/m^3]
    mu_t0: float = Atmosphere(h=h_t0).dynamic_viscosity[0]

    # #### CLIMB DEFINITIONS (at a single point) ###
    # NOTE: See post-processor

    # #### CRUISE DEFINITIONS ###
    # NOTE: See post-processor
    h_cruise: int = 18000 * ureg("ft").to("m").magnitude

    def __post_init__(self):
        # Takeoff
        atm_t0 = Atmosphere(h=self.h_t0)
        self.rho_t0 = atm_t0.density[0]
        self.mu_t0 = atm_t0.dynamic_viscosity[0]

        self.Dp_t0, self.C_Dp_t0 = calc_C_Dp(
            self.rho_t0,
            self.v_t0,
            self.mu_t0,
        )

        # Cruise
        atm_cruise = Atmosphere(h=self.h_cruise)
        self.rho_cruise = atm_cruise.density[0]
        self.mu_cruise = atm_cruise.dynamic_viscosity[0]

        self.Dp_cruise, self.C_Dp_cruise = calc_C_Dp(
            self.rho_cruise,
            self.v_cruise,
            self.mu_cruise,
        )

    # #### LANDING DEFINITIONS ###
    # NOTE: See post-processor
