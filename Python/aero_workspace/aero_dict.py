"""First plane uses a slightly smaller wing area. Performance for each plane may be evaluated, but start with plane 1 (AircraftConfig)"""

from dataclasses import dataclass

from ambiance import Atmosphere
from aero_workspace.conceptual_design import ureg
from aero_workspace.drag import (
    calc_C_Dp,
)  # TODO: Brenda's drag model needs to be updated based on revised plane geometry


# TODO: update constants in dataclass
@dataclass
class AircraftConfig:
    """All V2 aircraft definitons --> all in SI units"""

    # NOTE: Required inputs --> obtained from aero_main.py dataframe
    # All potential operating points can be seen under "NOTE(s)" below
    v_cruise: float = 0.0
    v_t0: float = 0.0

    # Main Wing
    AR: float = 7.0
    S: float = 45.0
    MAC: float = 2.54  # mean chord
    b: float = 17.75
    wing_incidence: float = 0.0

    # Control surfaces
    ail_hinge: float = 0.7
    flap_hinge: float = 0.75
    aileron_fraction: float = 0.75
    tail_hinge: float = 0.7

    # Horizontal Tail
    ht_AR: float = 3.0
    S_h: float = 13.33
    ht_MAC: float = 2.11
    V_h: float = 1.0  # S_h * l_h / (S * c)
    lh: float = 8.56

    # Vertical Tail
    vt_AR: float = 1.2
    S_v: float = 5.99
    vt_MAC: float = 2.23
    V_v: float = 0.06  # S_v * l_v / (S * b)
    lv: float = 8.0

    # Fuselage
    nose_x: float = -5.0
    fuse_width: float = 1.6

    # Prop
    tc_prime: float = 2.0  # prop value
    fan_radius: float = 0.583
    n_fans: int = 8
    blown_span: float = 11.71
    fan_length: float = 9.33
    hdisk: float = 0.73
    blowing_dy: float = 1.46

    # #### TAKEOFF DEFINITIONS ###
    # NOTE: Any operating point can be called from aero_main.py.
    #       Current sweep:
    #           "alphas": np.array([5, 10, 15, 20, 25])
    #           "velocities": np.array([20, 24, 28])
    #           "flap_deflections": np.array([50, 60, 65])
    h_t0: float = 0.0
    rho_t0 = Atmosphere(h=h_t0).density[0]  # [kg/m^3]
    mu_t0: float = Atmosphere(h=h_t0).dynamic_viscosity[0]

    #### CLIMB DEFINITIONS ###
    # NOTE: Any operating point can be called from aero_main.py.
    #       Current sweep:
    #           "alphas": np.array([10, 15, 20])
    #           "velocities": np.array([20, 30, 40, 50])
    #           "flap_deflections": np.array([0, 10, 20, 30, 40, 50])

    # #### CRUISE DEFINITIONS ###
    # NOTE: Any operating point can be called from aero_main.py.
    #       Current sweep:
    #            "alphas": np.array([0])
    #            "velocities": np.array([80.0, 100.0, 120.0, 125.0, 130.0, 140.0, 150.0])
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
    # NOTE: Any operating point can be called from aero_main.py.
    #       Current sweep:
    #           "alphas": np.array([1, 5, 10, 15, 20, 25])
    #           "velocities": np.linspace(1, 80, 9)
    #           "flap_deflections": np.array([50, 60, 65])


@dataclass
class AircraftConfig2:
    """All V2 aircraft definitons w/ smaller wing area --> all in SI units"""

    # NOTE: Required inputs --> obtained from aero_main.py dataframe
    # All potential operating points can be seen under "NOTE(s)" below
    v_cruise: float = 0.0
    v_t0: float = 0.0

    # Main Wing
    AR: float = 7.0
    S: float = 42.0
    MAC: float = 2.45
    b: float = 17.15
    wing_incidence: float = 0.0

    # Control surfaces
    ail_hinge: float = 0.7
    flap_hinge: float = 0.75
    aileron_fraction: float = 0.75
    tail_hinge: float = 0.7

    # Horizontal Tail
    ht_AR: float = 3.0
    S_h: float = 12.06
    ht_MAC: float = 2.01
    V_h: float = 1.0
    lh: float = 8.53

    # Vertical Tail
    vt_AR: float = 1.2
    S_v: float = 5.40
    vt_MAC: float = 2.12
    V_v: float = 0.06
    lv: float = 8.0

    # Fuselage
    nose_x: float = -5.0
    fuse_width: float = 1.6
    # fuselage (roskam pt 2) - (For drag modeling)
    D_f: float = 1.6  # [m]
    l_f: float = 15  # [m] f
    lambda_f: float = l_f / D_f

    # Prop
    tc_prime: float = 2.0
    fan_radius: float = 0.583
    n_fans: int = 8
    blown_span: float = 11.26
    fan_length: float = 9.33
    hdisk: float = 0.76
    blowing_dy: float = 1.41

    # #### TAKEOFF DEFINITIONS ###
    # NOTE: Any operating point can be called from aero_main.py.
    #       Current sweep:
    #           "alphas": np.array([5, 10, 15, 20, 25])
    #           "velocities": np.array([20, 24, 28])
    #           "flap_deflections": np.array([50, 60, 65])
    h_t0: float = 0.0
    rho_t0 = Atmosphere(h=h_t0).density[0]  # [kg/m^3]
    mu_t0: float = Atmosphere(h=h_t0).dynamic_viscosity[0]

    #### CLIMB DEFINITIONS ###
    # NOTE: Any operating point can be called from aero_main.py.
    #       Current sweep:
    #           "alphas": np.array([10, 15, 20])
    #           "velocities": np.array([20, 30, 40, 50])
    #           "flap_deflections": np.array([0, 10, 20, 30, 40, 50])

    # #### CRUISE DEFINITIONS ###
    # NOTE: Any operating point can be called from aero_main.py.
    #       Current sweep:
    #            "alphas": np.array([0])
    #            "velocities": np.array([80.0, 100.0, 120.0, 125.0, 130.0, 140.0, 150.0])
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
    # NOTE: Any operating point can be called from aero_main.py.
    #       Current sweep:
    #           "alphas": np.array([1, 5, 10, 15, 20, 25])
    #           "velocities": np.linspace(1, 80, 9)
    #           "flap_deflections": np.array([50, 60, 65])
