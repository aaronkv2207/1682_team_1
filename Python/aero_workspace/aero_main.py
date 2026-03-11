import numpy as np
import pandas as pd
from ambiance import Atmosphere
from conceptual_design import MTOW, V_CRUISE, V_STALL, W_S, S, ureg
from drag import C_Dp, V_v, V_h
from scipy.interpolate import interp1d
from aero_dict import AircraftConfig
# get jvl cg --> mset, run file, input file, ...; can shift around masses to shift cg in .mass file


class TakeoffCoeff:
    """Will read a summary of JVL output data as a .txt and interpolate results from JVL outputs at various operating points.
    Defines functions based on fit of operating points --> CL, CD, CM."""

    df = pd.read_csv("./takeoff_coefficients")
    alphas = df["ALPHA"]
    flap_angle = df["BETAS"]
    velocities = df["VELOCITY"]

    # TODO: Will need definitions for every control surface
    CL_alpha = np.hstack((alphas, df["CL"]))  # convert from degrees --> radians
    CD_alpha = np.hstack((alphas, df["CD"]))
    CM_alpha = np.hstack((alphas, df["CM"]))

    CL_flap = np.hstack((flap_angle, df["CL"]))
    CD_flap = np.hstack((flap_angle, df["CD"]))
    CM_flap = np.hstack((flap_angle, df["CM"]))

    CL_velocity = np.hstack((velocities, df["CL"]))
    CD_velocity = np.hstack((velocities, df["CD"]))
    CM_velocity = np.hstack((velocities, df["CM"]))


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

    CL_alpha = np.hstack((alphas, df["CL"]))
    CD_alpha = np.hstack((alphas, df["CD"]))
    CM_alpha = np.hstack((alphas, df["CM"]))

    CL_flap = np.hstack((flap_angle, df["CL"]))
    CD_flap = np.hstack((flap_angle, df["CD"]))
    CM_flap = np.hstack((flap_angle, df["CM"]))

    CL_velocity = np.hstack((velocities, df["CL"]))
    CD_velocity = np.hstack((velocities, df["CD"]))
    CM_velocity = np.hstack((velocities, df["CM"]))


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

    CL_alpha = np.hstack((alphas, df["CL"]))
    CD_alpha = np.hstack((alphas, df["CD"]))
    CM_alpha = np.hstack((alphas, df["CM"]))

    CL_flap = np.hstack((flap_angle, df["CL"]))
    CD_flap = np.hstack((flap_angle, df["CD"]))
    CM_flap = np.hstack((flap_angle, df["CM"]))

    CL_velocity = np.hstack((velocities, df["CL"]))
    CD_velocity = np.hstack((velocities, df["CD"]))
    CM_velocity = np.hstack((velocities, df["CM"]))


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
