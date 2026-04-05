import os
import pickle
import sys
from dataclasses import dataclass

import numpy as np
from ambiance import Atmosphere
from pint import UnitRegistry

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from aero_dict import AircraftConfig
from conceptual_design import ureg

DEG2RAD_CONV = ureg("deg").to("rad").magnitude
AR = 8


@dataclass
class AeroCoeffConfig:
    """Reads a summary of JVL output dataframe at various operating points. Defines functions
    based on operating points --> CL, CD, CM. If other parameters are desired, see data dictionary."""

    S: int
    phase: str  # 'takeoff', 'cruise', or 'landing'

    def __post_init__(self):
        file_path = f"JVL_writer/sref_trades/run_outputs/coeff_results/Sref-{self.S}/{self.phase}.pkl"

        with open(file_path, "rb") as f:
            data = pickle.load(f)

        self.alphas, self.velocities, self.d_elevator = [], [], []
        self.CL, self.CD, self.Cm, self.CDind = [], [], [], []

        for run in data:
            self.alphas.append(run["alpha"] * DEG2RAD_CONV)
            self.velocities.append(run["velocity"])
            self.d_elevator.append(run['d6'])
            self.CL.append(run["CL"])
            self.Cm.append(run["Cm"])

            _CDind = run["CL"] ** 2 / (AR * np.pi * run["e"])
            self.CDind.append(_CDind)

        if self.phase == "cruise":
            CD_DP = AircraftConfig.C_Dp_cruise
        else:
            CD_DP = (
                AircraftConfig.C_Dp_t0
            )  # Default for takeoff/landing; may need to modify in future iterations

        self.CD_tot = (self.CDind + CD_DP) * 1.2

        self.alphas = np.array(self.alphas)
        self.velocities = np.array(self.velocities)
        self.d_elevator = np.array(self.d_elevator)
        self.CL = np.array(self.CL)
        self.Cm = np.array(self.Cm)
        self.CDind = np.array(self.CDind)



if __name__ == "__main__":
    S_list = [44, 48, 52]  # NOTE: See available options in runner script.

    for S in S_list:
        print("")
        for phase in ["takeoff", "cruise", "landing"]:
            config = AeroCoeffConfig(S=S, phase=phase)
            print(f"\n=== {phase.capitalize()} (S = {S} [m^2]) ===")

            print(f"Velocities: {np.round(config.velocities, 4)}")
            print(f"AOAs: {np.round(config.alphas * (1 / DEG2RAD_CONV), 4)} [degrees]")
            print(f"d_elevator: {np.round(config.d_elevator, 4)}")
            print(f"CL: {np.round(config.CL, 4)}")
            print(f"CD_tot: {np.round(config.CD_tot, 4)}")
            print(f"Cm: {np.round(config.Cm, 4)}")
