from dataclasses import dataclass

from ambiance import Atmosphere
from conceptual_design import MTOW, V_CRUISE, V_STALL, W_S, X_TAKEOFF, S, cl_max, ureg
from drag import V_h, V_v, calc_C_Dp


# TODO: update constants in dataclass
@dataclass
class AircraftConfig:
    """All aircraft definitons --> all in SI units"""

    #### GLOBAL DEFINITIONS ###
    AR: float = 8.0
    s_ref: float = S.magnitude  # NOTE: S_ref may change after Brenda's drag update
    V_h: float = V_h  # TODO: will need to update based on stability analyses
    V_v: float = V_v  # TODO: will need to update based on stability analyses
    # parameters

    b = (s_ref * AR) ** 0.5  # [m]
    W_S: float = W_S.magnitude
    MTOW: float = MTOW.magnitude * 9.81
    c = 1.98  # NOTE: Brenda got this from Aaron --> Look into this
    # R = 2.39 / 2  # estimated
    x_to: float = X_TAKEOFF  # [m]

    #### TAKEOFF DEFINITIONS ###
    h_t0: float = 0.0
    rho_t0 = Atmosphere(h=h_t0).density[0]  # [kg/m^2]
    v_t0: float = 1.1 * V_STALL.magnitude
    mu_t0: float = Atmosphere(h=h_t0).dynamic_viscosity[0]

    Dp_t0, C_Dp_t0 = calc_C_Dp(rho_t0, v_t0, mu_t0)
    # Cd0_takeoff: float = ...  # TODO: couple with Brenda's drag model output
    # Cdv_takeoff: float = ...  # TODO: couple with drag model outputs
    # CDi_takeoff: float = ...  # TODO: couple with drag model outputs
    # e_takeoff: float = ...  # NOTE: specify e or CDi; make sure to pass in relevant parameter in function call
    # CL_takeoff: float = ...  # TODO: fill-in based on JVL outputs
    # CM_takeoff: float = ...  # TODO: fill-in based on JVL outputs

    # #### CLIMB DEFINITIONS ###
    # v_climb: int = ...
    # weight_climb: ...  # TODO: update to varied model
    # h_dot: ...
    # Cd0_climb: ...
    # Cdv_climb: ...
    # CDi_climb: ...
    # CM_climb: float = ...  # TODO: fill-in based on JVL outputs
    # # TODO: will need to integrate lift from control surfaces

    # #### CRUISE DEFINITIONS ###
    # weight_cruise: ...  # TODO: update to varied model
    h_cruise: int = 18000 * ureg("ft").to("m").magnitude
    rho_cruise = Atmosphere(h=h_cruise).density[0]  # [kg/m^2]
    v_cruise: float = 125  # [m/s]
    mu_cruise: float = Atmosphere(h=h_cruise).dynamic_viscosity[0]
    Dp_cruise, C_Dp_cruise = calc_C_Dp(rho_cruise, v_cruise, mu_cruise)

    # TODO: see if Brenda's model can define all stage parameters at the top, so function call only takes in a stage
    # Her script will also add up the drags for all components

    # Cd0_cruise: float = C_Dp(stage="cruise")
    # Cdv_cruise: float = 0.0
    # CDi_cruise: ...

    # CM_cruise: float = ...  # TODO: fill-in based on JVL outputs
    # # TODO: will need to integrate lift from control surfaces

    # #### LANDING DEFINITIONS ###
    # v_landing: float = 1.2 * V_STALL
    # Cd0_landing: ...  # TODO: couple with Brenda's drag model output
    # Cdv_landing: ...  # TODO: couple with drag model outputs
    # CDi_landing: ...  # TODO: couple with drag model outputs
    # CL_landing: float = ...  # TODO: fill-in based on JVL outputs
    # CM_landing: float = ...  # TODO: fill-in based on JVL outputs


# v, alpha, delta_flap angle --> takeoff
