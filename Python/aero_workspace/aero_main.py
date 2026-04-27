import pickle
from dataclasses import dataclass

import numpy as np
from aero_dict import AircraftConfig, AircraftConfig2
from pint import UnitRegistry

ureg = UnitRegistry()
DEG2RAD_CONV = ureg("deg").to("rad").magnitude


@dataclass
class AeroCoeffConfig:
    """Reads a summary of JVL output dataframe at various operating points. Defines functions
    based on operating points --> CL, CD, CM. If other parameters are desired, see data dictionary."""

    phase: str  # 'takeoff', 'cruise', or 'landing'
    # NOTE: CHANGE TO either AircraftConfig or AircraftConfig2
    aircraft: object = AircraftConfig
    max_elevator: float = 20.0  # degrees
    max_flap_unblown: float = 40.0  # degrees # TODO: Needs more rigorous analysis

    S = int(aircraft.S)
    AR = aircraft.AR

    def __post_init__(self):
        file_path = f"JVL_writer/v2_plane/run_outputs/coeff_results/Sref-{self.S}/{self.phase}.pkl"

        with open(file_path, "rb") as f:
            data = pickle.load(f)

        self.alphas, self.velocities, self.d_elevator = [], [], []
        self.CL, self.CD, self.Cm = [], [], []
        self.CDind, self.CD_tot, self.e = [], [], []
        self.flap_1, self.flap_2 = [], []

        for run in data:
            self.alphas.append(run["alpha"] * DEG2RAD_CONV)
            self.velocities.append(run["velocity"])
            self.d_elevator.append(run["d6"])
            self.CL.append(run["CL"])
            self.Cm.append(run["Cm"])
            self.e.append(run["e"])
            self.flap_1.append(run["d1"])
            self.flap_2.append(run["d2"])

            _CDind = run["CL"] ** 2 / (self.AR * np.pi * run["e"])
            self.CDind.append(_CDind)

            if self.phase == "cruise":
                CD_DP = AircraftConfig(v_cruise=run["velocity"]).C_Dp_cruise
            else:
                CD_DP = AircraftConfig(v_t0=run["velocity"]).C_Dp_t0

            self.CD_tot.append(
                [(_CDind + CD_DP) * 1.2]
            )  # TODO: Double check drag treatment

        self.alphas = np.array(self.alphas)
        self.velocities = np.array(self.velocities)
        self.d_elevator = np.array(self.d_elevator)
        self.CL = np.array(self.CL)
        self.Cm = np.array(self.Cm)
        self.CDind = np.array(self.CDind)
        self.CD_tot = np.array(self.CD_tot)
        self.e = np.array(self.e)
        self.flap_1 = np.array(self.flap_1)
        self.flap_2 = np.array(self.flap_2)

        # Filter based on reasonable elevator deflections
        elevator_mask = np.abs(self.d_elevator) < self.max_elevator
        flap_mask = (np.abs(self.flap_1) < self.max_flap_unblown) & (
            np.abs(self.flap_2) < self.max_flap_unblown
        )
        # cl_best_real_mask = np.argmax(self.CL[elevator_mask])

        if self.phase == "takeoff":
            self.alphas = self.alphas[elevator_mask]
            self.velocities = self.velocities[elevator_mask]
            self.d_elevator = self.d_elevator[elevator_mask]
            self.CL = self.CL[elevator_mask]
            self.Cm = self.Cm[elevator_mask]
            self.CDind = self.CDind[elevator_mask]
            self.CD_tot = self.CD_tot[elevator_mask]
            self.e = self.e[elevator_mask]
            self.flap_1 = self.flap_1[elevator_mask]
            self.flap_2 = self.flap_2[elevator_mask]

        if self.phase == "cruise":
            self.alphas = self.alphas[flap_mask]
            self.velocities = self.velocities[flap_mask]
            self.d_elevator = self.d_elevator[flap_mask]
            self.CL = self.CL[flap_mask]
            self.Cm = self.Cm[flap_mask]
            self.CDind = self.CDind[flap_mask]
            self.CD_tot = self.CD_tot[flap_mask]
            self.e = self.e[flap_mask]
            self.flap_1 = self.flap_1[flap_mask]
            self.flap_2 = self.flap_2[flap_mask]


if __name__ == "__main__":
    # NOTE: See all available operating conditions in JVL_writer/v2_plane/main_runner.py
    # for plane_idx, plane in enumerate([AircraftConfig, AircraftConfig2]):
    for plane_idx, plane in enumerate([AircraftConfig]):
        print(f"\n\n{'*' * 50}\n{'*' * 50}\nPlane #{plane_idx}")
        S = plane.S
        print(f"\n{'=' * 40}\nS = {S} [m^2] Performance")
        for phase in ["takeoff", "cruise", "climb", "landing"]:
            config = AeroCoeffConfig(phase=phase, aircraft=plane)
            print(f"\n=== {phase.upper()} (S = {S} [m^2]) ===")

            print(f"Velocities: {np.round(config.velocities, 4)}")
            print(f"AOAs: {np.round(config.alphas * (1 / DEG2RAD_CONV), 4)} [°]")
            print(f"d_elevator: {np.round(config.d_elevator, 4)}")
            print(f"d_flap: {np.round(config.flap_1, 4)}")
            print(f"CL: {np.round(config.CL, 4)}")
            print(f"CD_tot: {np.round(config.CD_tot, 4)}")
            print(f"Cm: {np.round(config.Cm, 4)}")
            print(f"e: {np.round(config.e, 4)}")
