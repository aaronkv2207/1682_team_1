import pickle
from dataclasses import dataclass

import jvl_run_outputs
import numpy as np
import pandas as pd
from aero_dict import AircraftConfig
from ambiance import Atmosphere
from conceptual_design import MTOW, V_CRUISE, V_STALL, W_S, S, ureg

# get jvl cg --> mset, run file, input file, ...; can shift around masses to shift cg in .mass file

DEG2RAD_CONV = ureg("deg").to("rad").magnitude


@dataclass
class TakeoffCoeff:
    """Reads a summary of JVL output dataframe at various operating points. Defines functions
    based on operating points --> CL, CD, CM. If other parameters are desired, see data dictionary."""

    FILE_NAME = "Python/aero_workspace/jvl_run_outputs/takeoff.pkl"
    with open(FILE_NAME, "rb") as f:
        data = pickle.load(f)

    alphas, velocities, betas = [], [], []
    CL, CD, Cm, CDind, e = [], [], [], [], []

    for _, run in enumerate(data):
        alphas.append(run["alpha"] * DEG2RAD_CONV)
        velocities.append(run["velocity"])
        betas.append(run["beta"] * DEG2RAD_CONV)
        CL.append(run["CL"])  # NOTE: double check definitions: Cl vs CL
        # CD.append(run["CD"])
        Cm.append(run["Cm"])

        CDind.append(run["CDind"])
        e.append(run["e"])

        # NOTE: For drag buildup, CD_ind is obtained from JVL;
        #       CD_form and CD_visc is obtained from Brenda's drag build up
    CD_DP = AircraftConfig.C_Dp_t0  # profile drag from Brenda's model
    CD_tot = CDind + CD_DP

@dataclass
class ClimbCoeff:  # TODO: NO DATA IMPLEMENTED; NEEDS UPDATE
    """Will read a summary of JVL output dataframe at various operating points. Defines functions
    based on operating points --> CL, CD, CM. If other parameters are desired, see data dictionary."""

    FILE_NAME = "Python/aero_workspace/jvl_run_outputs/climb.pkl"
    with open(FILE_NAME, "rb") as f:
        data = pickle.load(f)

    alphas, velocities, betas = [], [], []
    CL, CD, Cm, CDind, e = [], [], [], [], []

    for _, run in enumerate(data):
        alphas.append(run["alpha"] * DEG2RAD_CONV)
        velocities.append(run["velocity"])
        betas.append(run["beta"] * DEG2RAD_CONV)
        CL.append(run["CL"])  # NOTE: double check definitions: Cl vs CL
        # CD.append(run["CD"])
        Cm.append(run["Cm"])

        CDind.append(run["CDind"])
        e.append(run["e"])

        # NOTE: For drag buildup, CD_ind is obtained from JVL;
        #       CD_form and CD_visc is obtained from Brenda's drag build up
    # CD_DP = AircraftConfig.C_Dp_climb  # profile drag from Brenda's model
    # CD_tot = (CDind + CD_DP)

    # TODO: CD_DP = AircraftConfig.C_Dp_climb # profile drag from Brenda's model


@dataclass
class CruiseCoeff:  # TODO: NO DATA IMPLEMENTED; NEEDS UPDATE
    """Will read a summary of JVL output dataframe at various operating points. Defines functions
    based on operating points --> CL, CD, CM. If other parameters are desired, see data dictionary."""

    FILE_NAME = "Python/aero_workspace/jvl_run_outputs/cruise.pkl"
    with open(FILE_NAME, "rb") as f:
        data = pickle.load(f)

    alphas, velocities, betas = [], [], []
    CL, CD, Cm, CDind, e = [], [], [], [], []

    for _, run in enumerate(data):
        alphas.append(run["alpha"] * DEG2RAD_CONV)
        velocities.append(run["velocity"])
        betas.append(run["beta"] * DEG2RAD_CONV)
        CL.append(run["CL"])  # NOTE: double check definitions: Cl vs CL
        # CD.append(run["CD"])
        Cm.append(run["Cm"])

        CDind.append(run["CDind"])
        e.append(run["e"])

        # NOTE: For drag buildup, CD_ind is obtained from JVL;
        #       CD_form and CD_visc is obtained from Brenda's drag build up
    CD_DP = AircraftConfig.C_Dp_cruise  # profile drag from Brenda's model
    CD_tot = CDind + CD_DP


class CruiseModel:
    def __init__(
        self, s_ref, weight, v_cruise, h_cruise, AR, Cd0, Cdv, CDi, e=0.7
    ) -> None:
        self.s_ref = s_ref * ureg("m^2")
        self.weight = weight * ureg("N")  # newtons
        self.v_cruise = v_cruise
        self.h_cruise = h_cruise
        self.density = Atmosphere(h=h_cruise).density[0]
        self.AR = AR
        self.q = 0.5 * self.density * (self.v_cruise**2)
        self.Cd0 = Cd0
        self.Cdv = Cdv
        self.CDi = CDi
        self.e = e

    def cl(self):
        L = self.weight  # assumes I have weight as a function of time
        return L / (self.q * self.s_ref)

    def cd_induced(self):
        if self.CDi is None:
            return (self.cl() ** 2) / (np.pi * self.AR * self.e)
        else:
            return self.CDi

    def cd_total(self):
        return self.cd_induced() + self.Cd0 + self.Cdv

    def drag_total(self):
        return self.cd_total() * self.q * self.s_ref


class LandingCoeff:  # TODO: NO DATA IMPLEMENTED; NEEDS UPDATE
    """Will read a summary of JVL output dataframe at various operating points. Defines functions
    based on operating points --> CL, CD, CM. If other parameters are desired, see data dictionary."""

    FILE_NAME = "Python/aero_workspace/jvl_run_outputs/landing.pkl"
    with open(FILE_NAME, "rb") as f:
        data = pickle.load(f)

    alphas, velocities, betas = [], [], []
    CL, CD, Cm, CDind, e = [], [], [], [], []

    for _, run in enumerate(data):
        alphas.append(run["alpha"] * DEG2RAD_CONV)
        velocities.append(run["velocity"])
        betas.append(run["beta"] * DEG2RAD_CONV)
        CL.append(run["CL"])  # NOTE: double check definitions: Cl vs CL
        # CD.append(run["CD"])
        Cm.append(run["Cm"])

        CDind.append(run["CDind"])
        e.append(run["e"])
    # # TODO: CD_DP = AircraftConfig.C_Dp_landing # profile drag from Brenda's model
    #     # NOTE: For drag buildup, CD_ind is obtained from JVL;
    #     #       CD_form and CD_visc is obtained from Brenda's drag build up
    # CD_DP = AircraftConfig.C_Dp_landing  # profile drag from Brenda's model
    # CD_tot = (CDind + CD_DP)


# Runner script
if __name__ == "__main__":
    cruise_cls = CruiseModel(
        s_ref=AircraftConfig.s_ref,
        weight=0.85
        * MTOW.magnitude,  # TODO: needed from prop; UPDATE to AircraftConfig.weight_cruise
        v_cruise=CruiseCoeff.velocities[
            0
        ],  # NOTE: change for a particular run of interest
        h_cruise=AircraftConfig.h_cruise,
        AR=AircraftConfig.AR,
        Cd0=0.0,
        Cdv=TakeoffCoeff.CD_DP,
        CDi=CruiseCoeff.CDind[0],  # NOTE: change for a particular run of interest
        e=CruiseCoeff.e[0],  # NOTE: change for a particular run of interest
    )  #  default assumes e=0.7; change this parameter if different
    CD_total_cruise = cruise_cls.cd_total()
    L_over_D_cruise = cruise_cls.cl() / CD_total_cruise
