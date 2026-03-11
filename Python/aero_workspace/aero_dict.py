from dataclasses import dataclass


# TODO: update constants in dataclass
@dataclass
class AircraftConfig:
    """All aircraft definitons"""

    #### GLOBAL DEFINITIONS ###
    AR: int = 8
    s_ref: float = S  # NOTE: S_ref may change after Brenda's drag update
    V_h: float = ...  # TODO: will need to update Brenda's script
    V_v: float = ...
    # parameters
    AR = AircraftConfig.AR # [m^2]
    h_cruise = AircraftConfig.h_cruise # [m]
    rho_cruise = Atmosphere(h=h_cruise).density[0]  # [kg/m^2]

    b = 19.9  # [m]
    c = 1.98  # [m]
    S = b**2 / AR
    # R = 2.39 / 2  # estimated
    x_to = 45  # [m]

    m = 5670  # [kg]
    W = m * 9.81  # [N]



    #### TAKEOFF DEFINITIONS ###
    v_takeoff: float = 1.2 * V_STALL
    Cd0_takeoff: float = ...  # TODO: couple with Brenda's drag model output
    Cdv_takeoff: float = ...  # TODO: couple with drag model outputs
    CDi_takeoff: float = ...  # TODO: couple with drag model outputs
    e_takeoff: float = ...  # NOTE: specify e or CDi; make sure to pass in relevant parameter in function call
    CL_takeoff: float = ...  # TODO: fill-in based on JVL outputs
    CM_takeoff: float = ...  # TODO: fill-in based on JVL outputs

    #### CLIMB DEFINITIONS ###
    v_climb: int = ...
    weight_climb: ...  # TODO: update to varied model
    h_dot: ...
    Cd0_climb: ...
    Cdv_climb: ...
    CDi_climb: ...
    CM_climb: float = ...  # TODO: fill-in based on JVL outputs
    # TODO: will need to integrate lift from control surfaces

    #### CRUISE DEFINITIONS ###
    v_cruise: int = V_CRUISE
    weight_cruise: ...  # TODO: update to varied model
    h_cruise: int = 18000 * ureg("ft").to("m")

    # TODO: see if Brenda's model can define all stage parameters at the top, so function call only takes in a stage
    # Her script will also add up the drags for all components

    Cd0_cruise: float = C_Dp(stage="cruise")
    Cdv_cruise: float = 0.0
    CDi_cruise: ...

    CM_cruise: float = ...  # TODO: fill-in based on JVL outputs
    # TODO: will need to integrate lift from control surfaces

    #### LANDING DEFINITIONS ###
    v_landing: float = 1.2 * V_STALL
    Cd0_landing: ...  # TODO: couple with Brenda's drag model output
    Cdv_landing: ...  # TODO: couple with drag model outputs
    CDi_landing: ...  # TODO: couple with drag model outputs
    CL_landing: float = ...  # TODO: fill-in based on JVL outputs
    CM_landing: float = ...  # TODO: fill-in based on JVL outputs


# v, alpha, delta_flap angle --> takeoff