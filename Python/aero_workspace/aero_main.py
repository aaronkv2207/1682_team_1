from dataclasses import dataclass

import numpy as np
import pandas as pd
from ambiance import Atmosphere
from conceptual_design import MTOW, V_CRUISE, V_STALL, W_S, S, ureg
from drag import C_Dp, V_v, V_h
from scipy.interpolate import interp1d

# get jvl cg --> mset, run file, input file, ...; can shift around masses to shift cg in .mass file


# TODO: update constants in dataclass
@dataclass
class AircraftConfig:
    """All aircraft definitons"""

    #### GLOBAL DEFINITIONS ###
    AR: int = 8
    s_ref: float = S  # NOTE: S_ref may change after Brenda's drag update
    V_h: float = V_h  # TODO: will need to update Brenda's script
    V_v: float = V_v

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


class TakeoffCoeff:
    """Will read a summary of JVL output data as a .txt and interpolate results from JVL outputs at various operating points.
    Defines functions based on fit of operating points --> CL, CD, CM."""

    df = pd.read_csv("./takeoff_coefficients")
    alphas = df["ALPHA"]
    flap_angle = df["BETAS"]
    velocities = df["VELOCITY"]

    # TODO: Will need definitions for every control surface
    CL_alpha = np.vstack((alphas, df["CL"]))  # convert from degrees --> radians
    CD_alpha = np.vstack((alphas, df["CD"]))
    CM_alpha = np.vstack((alphas, df["CM"]))

    CL_flap = np.vstack((flap_angle, df["CL"]))
    CD_flap = np.vstack((flap_angle, df["CD"]))
    CM_flap = np.vstack((flap_angle, df["CM"]))

    CL_velocity = np.vstack((velocities, df["CL"]))
    CD_velocity = np.vstack((velocities, df["CD"]))
    CM_velocity = np.vstack((velocities, df["CM"]))


class ClimbCoeff:  # TODO: NOT IMPLEMENTED; NEEDS UPDATE
    """Will read a summary of JVL output data as a .txt and interpolate results from JVL outputs at various operating points.
    Defines functions based on fit of operating points --> CL, CD, CM."""


class CruiseCoeff:  # TODO: NOT IMPLEMENTED; NEEDS UPDATE
    """Will read a summary of JVL output data as a .txt and interpolate results from JVL outputs at various operating points.
    Defines functions based on fit of operating points --> CL, CD, CM."""

    df = pd.read_csv("./cruise_coefficients")
    alphas = df["ALPHA"]
    flap_angle = df["BETAS"]
    velocities = df["VELOCITY"]

    CL_alpha = np.vstack((alphas, df["CL"]))
    CD_alpha = np.vstack((alphas, df["CD"]))
    CM_alpha = np.vstack((alphas, df["CM"]))

    CL_flap = np.vstack((flap_angle, df["CL"]))
    CD_flap = np.vstack((flap_angle, df["CD"]))
    CM_flap = np.vstack((flap_angle, df["CM"]))

    CL_velocity = np.vstack((velocities, df["CL"]))
    CD_velocity = np.vstack((velocities, df["CD"]))
    CM_velocity = np.vstack((velocities, df["CM"]))


class CruiseModel:
    def __init__(
        self, s_ref, weight, v_cruise, h_cruise, AR, Cd0, Cdv, CDi, e=0.7
    ) -> None:
        self.s_ref = s_ref * ureg("m^2")
        self.weight = weight * ureg("N")  # newtons
        self.v_cruise = v_cruise.magnitude * ureg("m/s")
        self.h_cruise = h_cruise
        self.density = Atmosphere(h=h_cruise.magnitude).density[0] * ureg("kg/m^3")
        self.AR = AR
        self.q = 0.5 * self.density * (self.v_cruise**2)
        self.Cd0 = Cd0
        self.Cdv = Cdv
        self.CDi = CDi
        self.e = e

    def cl(self):
        L = self.weight  # assumes I have weight as a function of time
        return (L / (self.q * self.s_ref)).magnitude

    def cd_induced(self):
        if self.CDi is None:
            return (self.cl() ** 2) / (np.pi * self.AR * self.e)
        else:
            return self.CDi

    def cd_total(self):
        return self.cd_induced() + self.Cd0 + self.Cdv

    def drag_total(self):
        return self.cd_total() * self.q * self.s_ref


class LandingCoeff:  # TODO: NOT IMPLEMENTED; NEEDS UPDATE
    """Will read a summary of JVL output data as a .txt and interpolate results from JVL outputs at various operating points.
    Defines functions based on fit of operating points --> CL, CD, CM."""

    df = pd.read_csv("./landing_coefficients")
    alphas = df["ALPHA"]
    flap_angle = df["BETAS"]
    velocities = df["VELOCITY"]

    CL_alpha = np.vstack((alphas, df["CL"]))
    CD_alpha = np.vstack((alphas, df["CD"]))
    CM_alpha = np.vstack((alphas, df["CM"]))

    CL_flap = np.vstack((flap_angle, df["CL"]))
    CD_flap = np.vstack((flap_angle, df["CD"]))
    CM_flap = np.vstack((flap_angle, df["CM"]))

    CL_velocity = np.vstack((velocities, df["CL"]))
    CD_velocity = np.vstack((velocities, df["CD"]))
    CM_velocity = np.vstack((velocities, df["CM"]))


# Runner script
if __name__ == "__main__":
    cruise_cls = CruiseModel(
        s_ref=AircraftConfig.s_ref,
        weight=AircraftConfig.weight_cruise,
        v_cruise=AircraftConfig.v_cruise,
        h_cruise=AircraftConfig.h_cruise,
        AR=AircraftConfig.AR,
        Cd0=AircraftConfig.Cd0_cruise,
        Cdv=AircraftConfig.Cdv_cruise,
        CDi=AircraftConfig.CDi_cruise,
    )  #  default assumes e=0.7; change this parameter if different
    CD_total_cruise = cruise_cls.cd_total()
    L_over_D_cruise = cruise_cls.cl() / CD_total_cruise
