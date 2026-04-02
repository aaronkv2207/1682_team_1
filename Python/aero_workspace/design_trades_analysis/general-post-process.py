import os
import pickle
import sys
import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from ambiance import Atmosphere
from pint import UnitRegistry

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from aero_dict import AircraftConfig
from conceptual_design import MTOW, V_CRUISE, V_STALL, W_S, S, ureg

DEG2RAD_CONV = ureg("deg").to("rad").magnitude


@dataclass
class DesiredCoeff:
    """Reads a summary of JVL output dataframe at various operating points. Defines functions
    based on operating points --> CL, CD, CM. If other parameters are desired, see data dictionary."""

    try:
        FILE_NAME = "..."  # TODO: User-specified based on the case you're interested in analyzing
        with open(FILE_NAME, "rb") as f:
            data = pickle.load(f)

    except FileNotFoundError:
        raise NotImplementedError(
            "User must specify FILE_NAME based on the case you're interested in analyzing"
        )

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
    CD_tot = (CDind + CD_DP) * 1.2
