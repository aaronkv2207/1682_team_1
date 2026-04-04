import os
import pickle
import sys
from dataclasses import dataclass

import numpy as np
import pandas as pd
from ambiance import Atmosphere
from pint import UnitRegistry

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from aero_dict import AircraftConfig
from conceptual_design import ureg

DEG2RAD_CONV = ureg("deg").to("rad").magnitude
S = 52  # NOTE: User-specified value based on the case you're interested in analyzing. See available options in runner script.
AR = 8


@dataclass
class TakeoffCoeff_config:
    """Reads a summary of JVL output dataframe at various operating points. Defines functions
    based on operating points --> CL, CD, CM. If other parameters are desired, see data dictionary."""

    try:
        FILE_NAME = (
            f"JVL_writer/sref_trades/run_outputs/coeff_results/Sref-{S}/takeoff.pkl"
        )
        with open(FILE_NAME, "rb") as f:
            data = pickle.load(f)

    except FileNotFoundError:
        raise FileNotFoundError(
            "User must correctly specify FILE_NAME based on the case you're interested in analyzing"
        )

    alphas, velocities, betas = [], [], []
    CL, CD, Cm, CDind, e = [], [], [], [], []

    for _, run in enumerate(data):
        alphas.append(run["alpha"] * DEG2RAD_CONV)
        velocities.append(run["velocity"])
        CL.append(run["CL"])  # NOTE: double check definitions: Cl vs CL
        # CD.append(run["CD"])
        Cm.append(run["Cm"])

        _CDind = run["CL"] ** 2 / (AR * np.pi * run["e"])
        CDind.append(_CDind)

        # NOTE: For drag buildup, CD_ind is obtained from JVL;
        #       CD_form and CD_visc is obtained from Brenda's drag build up
    CD_DP = AircraftConfig.C_Dp_t0  # profile drag from Brenda's model
    CD_tot = (CDind + CD_DP) * 1.2


@dataclass
class CruiseCoeff_config:
    """Reads a summary of JVL output dataframe at various operating points. Defines functions
    based on operating points --> CL, CD, CM. If other parameters are desired, see data dictionary."""

    try:
        FILE_NAME = (
            f"JVL_writer/sref_trades/run_outputs/coeff_results/Sref-{S}/cruise.pkl"
        )
        with open(FILE_NAME, "rb") as f:
            data = pickle.load(f)

    except FileNotFoundError:
        raise FileNotFoundError(
            "User must correctly specify FILE_NAME based on the case you're interested in analyzing"
        )

    alphas, velocities, betas = [], [], []
    CL, CD, Cm, CDind, e = [], [], [], [], []

    for _, run in enumerate(data):
        alphas.append(run["alpha"] * DEG2RAD_CONV)
        velocities.append(run["velocity"])
        CL.append(run["CL"])  # NOTE: double check definitions: Cl vs CL
        # CD.append(run["CD"])
        Cm.append(run["Cm"])

        _CDind = run["CL"] ** 2 / (AR * np.pi * run["e"])
        CDind.append(_CDind)

        # NOTE: For drag buildup, CD_ind is obtained from JVL;
        #       CD_form and CD_visc is obtained from Brenda's drag build up
    CD_DP = AircraftConfig.C_Dp_cruise  # profile drag from Brenda's model
    CD_tot = (CDind + CD_DP) * 1.2


@dataclass  # NOTE: Validate drag build-up
class LandingCoeff_config:
    """Reads a summary of JVL output dataframe at various operating points. Defines functions
    based on operating points --> CL, CD, CM. If other parameters are desired, see data dictionary."""

    # raise NotImplementedError("Not implemented in V1 sizing. Will need to update configuration in runner.py")
    try:
        FILE_NAME = (
            f"JVL_writer/sref_trades/run_outputs/coeff_results/Sref-{S}/landing.pkl"
        )
        with open(FILE_NAME, "rb") as f:
            data = pickle.load(f)

    except FileNotFoundError:
        raise FileNotFoundError(
            "User must correctly specify FILE_NAME based on the case you're interested in analyzing"
        )

    alphas, velocities, betas = [], [], []
    CL, CD, Cm, CDind, e = [], [], [], [], []

    for _, run in enumerate(data):
        alphas.append(run["alpha"] * DEG2RAD_CONV)
        velocities.append(run["velocity"])
        CL.append(run["CL"])  # NOTE: double check definitions: Cl vs CL
        # CD.append(run["CD"])
        Cm.append(run["Cm"])

        _CDind = run["CL"] ** 2 / (AR * np.pi * run["e"])
        CDind.append(_CDind)

        # NOTE: For drag buildup, CD_ind is obtained from JVL;
        #       CD_form and CD_visc is obtained from Brenda's drag build up
    CD_DP = AircraftConfig.C_Dp_t0  # profile drag from Brenda's model
    CD_tot = (CDind + CD_DP) * 1.2


if __name__ == "__main__":
    print("")
    for config in [TakeoffCoeff_config, CruiseCoeff_config, LandingCoeff_config]:
        name = config.__name__.replace("Coeff_config", "")
        print(f"\n*** {name} Configuration ***")

        print(f"Velocities: {np.round(config.velocities, 4)}")
        print(f"AOAs: {np.round(config.alphas, 4)}")
        print(f"CL: {np.round(config.CL, 4)}")
        print(f"CD_tot: {np.round(config.CD_tot, 4)}")
        print(f"Cm: {np.round(config.Cm, 4)}")
