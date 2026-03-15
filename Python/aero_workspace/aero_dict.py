from dataclasses import dataclass

from ambiance import Atmosphere
from conceptual_design import V_CRUISE, V_STALL, W_S, S, cl_max, ureg, MTOW, X_TAKEOFF


# TODO: update constants in dataclass
@dataclass
class AircraftConfig:
    """All aircraft definitons --> all in SI units"""

    #### GLOBAL DEFINITIONS ###
    AR: int = 8.0
    s_ref: float = S  # NOTE: S_ref may change after Brenda's drag update
    V_h: float = ...  # TODO: will need to update Brenda's script
    V_v: float = ...
    # parameters

    b = (S * AR) ** 0.5  # [m]
    S: float = b**2 / AR
    W_S: float = W_S
    MTOW: float = MTOW.magnitude * 9.81
    c = 1.98  # NOTE: Brenda got this from Aaron --> Look into this
    # R = 2.39 / 2  # estimated
    x_to: float = X_TAKEOFF  # [m]

    #### TAKEOFF DEFINITIONS ###
    v_takeoff: float = 1.2 * V_STALL
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
    # v_cruise: int = V_CRUISE
    # weight_cruise: ...  # TODO: update to varied model
    h_cruise: int = 18000 * ureg("ft").to("m").magnitude
    rho_cruise = Atmosphere(h=h_cruise).density[0]  # [kg/m^2]

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
