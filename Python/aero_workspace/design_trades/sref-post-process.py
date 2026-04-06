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
    max_elevator: float = 20.0  # degrees

    def __post_init__(self):
        file_path = f"JVL_writer/sref_trades/run_outputs/coeff_results/Sref-{self.S}/{self.phase}.pkl"

        with open(file_path, "rb") as f:
            data = pickle.load(f)

        self.alphas, self.velocities, self.d_elevator = [], [], []
        self.CL, self.CD, self.Cm, self.CDind, self.e = [], [], [], [], []

        for run in data:
            self.alphas.append(run["alpha"] * DEG2RAD_CONV)
            self.velocities.append(run["velocity"])
            self.d_elevator.append(run["d6"])
            self.CL.append(run["CL"])
            self.Cm.append(run["Cm"])
            self.e.append(run["e"])

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
        self.e = np.array(self.e)

        # Filter based on reasonable elevator deflections
        elevator_mask = np.abs(self.d_elevator) < self.max_elevator
        cl_best_real_mask = np.argmax(self.CL[elevator_mask])

        # produces best
        if self.phase == "takeoff":  # for cl_max
            self.alphas = self.alphas[cl_best_real_mask]
            self.velocities = self.velocities[cl_best_real_mask]
            self.d_elevator = self.d_elevator[cl_best_real_mask]
            self.CL = self.CL[cl_best_real_mask]
            self.Cm = self.Cm[cl_best_real_mask]
            self.CDind = self.CDind[cl_best_real_mask]
            self.CD_tot = self.CD_tot[cl_best_real_mask]
            self.e = self.e[1]

        if self.phase == "cruise":  # for cl_max; assuming vcruise of 125 [m/s]
            self.alphas = self.alphas[1]
            self.velocities = self.velocities[1]
            self.d_elevator = self.d_elevator[1]
            self.CL = self.CL[1]
            self.Cm = self.Cm[1]
            self.CDind = self.CDind[1]
            self.CD_tot = self.CD_tot[1]
            self.e = self.e[1]


if __name__ == "__main__":
    S_list = [44, 48, 52]  # NOTE: See available options in runner script.

    for idx, S in enumerate(S_list):
        print(f"\n{'=' * 40}\nS = {S} [m^2] Performance")
        for phase in ["takeoff", "cruise", "landing"]:
            config = AeroCoeffConfig(S=S, phase=phase)
            print(f"\n=== {phase.upper()} (S = {S} [m^2]) ===")

            print(f"Velocities: {np.round(config.velocities, 4)}")
            print(f"AOAs: {np.round(config.alphas * (1 / DEG2RAD_CONV), 4)} [degrees]")
            print(f"d_elevator: {np.round(config.d_elevator, 4)}")
            print(f"CL: {np.round(config.CL, 4)}")
            print(f"CD_tot: {np.round(config.CD_tot, 4)}")
            print(f"Cm: {np.round(config.Cm, 4)}")
            print(f"e: {np.round(config.e, 4)}")
