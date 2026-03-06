from dataclasses import dataclass

import numpy as np
from ambiance import Atmosphere
from conceptual_design import MTOW, V_STALL, V_CRUISE, W_S, S, ureg


# TODO: update constants in dataclass
@dataclass
class AircraftConfig:
    """All aircraft definitons"""

    #### GLOBAL DEFINITIONS ###
    AR: int = 8
    s_ref: float = S

    #### TAKEOFF DEFINITIONS ###
    v_takeoff: float = 1.2 * V_STALL
    Cd0_takeoff: ...  # TODO: couple with Brenda's drag model output
    Cdv_takeoff: ...  # TODO: couple with drag model outputs
    CL_takeoff: float = ...  # TODO: fill-in based on JVL outputs
    CM_takeoff: float = ...  # TODO: fill-in based on JVL outputs

    #### CRUISE DEFINITIONS ###
    v_cruise: int = V_CRUISE
    weight_cruise: ...  # TODO: update to varied model
    h_cruise: ...
    Cd0_cruise: ...
    Cdv_cruise: ...

    #### LANDING DEFINITIONS ###
    v_landing: float = 1.2 * V_STALL
    Cd0_landing: ...  # TODO: couple with Brenda's drag model output
    Cdv_landing: ...  # TODO: couple with drag model outputs
    CL_landing: float = ...  # TODO: fill-in based on JVL outputs
    CM_landing: float = ...  # TODO: fill-in based on JVL outputs


class TakeOff:
    """Will read a summary of JVL output data as a .txt and interpolate results from JVL outputs at various operating points.
    Return polynomial fit of operating points.

    Returns:
        array: Functions of CL, CD, CM
    """


class CruiseModel:
    def __init__(self, s_ref, weight, v_cruise, h_cruise, AR, Cd0, Cdv, e=None) -> None:
        self.s_ref = s_ref * ureg("m^2")
        self.weight = weight * ureg("N")  # newtons
        self.v_cruise = v_cruise.magnitude * ureg("m/s")
        self.h_cruise = h_cruise
        self.density = Atmosphere(h=h_cruise.magnitude).density[0] * ureg("kg/m^3")
        self.AR = AR
        self.e = e
        self.q = 0.5 * self.density * (self.v_cruise**2)
        self.Cd0 = Cd0
        self.Cdv = Cdv

    def cl(self):
        L = self.weight  # assumes I have weight as a function of time
        return (L / (self.q * self.s_ref)).magnitude

    def cd_induced(self):
        if self.e is None:
            return "Query JVL output"
        else:
            return (self.cl() ** 2) / (np.pi * self.AR * self.e)

    def cd_parasitic(self):
        return self.Cd0

    def cd_viscous(self):
        return self.Cdv

    def cd_total(self):
        return self.cd_induced() + self.cd_parasitic() + self.cd_viscous()

    def drag_total(self):
        return self.cd_total() * self.q * self.s_ref


class Landing:
    """Will read a summary of JVL output data as a .txt and interpolate results from JVL outputs at various operating points.
    Return polynomial fit of operating points.

    Returns:
        array: Functions of CL, CD, CM
    """


# Runner script
if __name__ == "__main__":
    cruise_cls = CruiseModel(
        s_ref=AircraftConfig.s_ref,
        weight=AircraftConfig.weight_cruise,  # TODO: update to varied model
        v_cruise=AircraftConfig.v_cruise,  # # TODO: update to varied model
        h_cruise=AircraftConfig.h_cruise,
        AR=AircraftConfig.AR,
        Cd0=AircraftConfig.Cd0_cruise,
        Cdv=AircraftConfig.Cdv_cruise,
    )
    CD_total_cruise = cruise_cls.cd_total()
    L_over_D_cruise = cruise_cls.cl() / CD_total_cruise
