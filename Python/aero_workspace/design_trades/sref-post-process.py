import os
import pickle
import sys
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from pint import UnitRegistry

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from aero_dict import AircraftConfig
from conceptual_design import ureg

DEG2RAD_CONV = ureg("deg").to("rad").magnitude
METERS2FT_CONV = ureg("m").to("ft").magnitude
AR = 8


def get_runway_length(
    v,
    S,
    CL,
    CD,
    mu=AircraftConfig.mu_t0,
    rho=AircraftConfig.rho_t0,
    T=27496,
    m=7504,
    g=9.81,
):
    W = m * g
    L = 0.5 * rho * v**2 * S * CL
    D = 0.5 * rho * v**2 * S * CD

    net_force = T - D - mu * (W - L)

    if net_force <= 0:
        raise ValueError("Insufficient thrust to accelerate (denominator <= 0).")

    x_runway = (m * v**2) / (2 * net_force)
    return x_runway


def get_climb_gradient(
    v,
    S,
    CD,
    rho=AircraftConfig.rho_climb,
    T=27496,
    m=7504,
    g=9.81,
):
    W = m * g
    q = 0.5 * rho * v**2
    D = q * S * CD
    gamma = (T - D) / W

    return gamma


@dataclass
class AeroCoeffConfig:
    """Reads a summary of JVL output dataframe at various operating points. Defines functions
    based on operating points --> CL, CD, CM. If other parameters are desired, see data dictionary."""

    S: int
    phase: str  # 'takeoff', 'cruise', or 'landing'
    name: str  # specify sub-folder
    max_elevator: float = 20.0  # degrees

    def __post_init__(self):
        file_path = f"JVL_writer/sref_trades/run_outputs/{self.name}/coeff_results/Sref-{self.S}/{self.phase}.pkl"

        with open(file_path, "rb") as f:
            data = pickle.load(f)

        self.alphas, self.velocities, self.d_elevator = [], [], []
        self.CL, self.CD, self.Cm, self.CDind, self.e, self.xto, self.gamma = (
            [],
            [],
            [],
            [],
            [],
            [],
            [],
        )

        for run in data:
            self.alphas.append(run["alpha"] * DEG2RAD_CONV)
            self.velocities.append(run["velocity"])
            self.d_elevator.append(run["d6"])
            self.CL.append(run["CL"])
            self.Cm.append(run["Cm"])
            self.e.append(run["e"])

            _CDind = run["CL"] ** 2 / (AR * np.pi * run["e"])
            self.CDind.append(_CDind)

            if self.phase in ("takeoff", "landing"):
                self.xto.append(
                    get_runway_length(
                        v=run["velocity"], S=self.S, CL=run["CL"], CD=run["CD"]
                    )
                    * METERS2FT_CONV
                )

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

        # produces best - for cl_max
        if self.phase == "takeoff":
            self.alphas = self.alphas[cl_best_real_mask]
            self.velocities = self.velocities[cl_best_real_mask]
            self.d_elevator = self.d_elevator[cl_best_real_mask]
            self.CL = self.CL[cl_best_real_mask]
            self.Cm = self.Cm[cl_best_real_mask]
            self.CDind = self.CDind[cl_best_real_mask]
            self.CD_tot = self.CD_tot[cl_best_real_mask]
            self.e = self.e[cl_best_real_mask]
            self.xto = np.array(self.xto)[cl_best_real_mask]

        if self.phase == "cruise":  # assuming vcruise of 125 [m/s]
            self.alphas = self.alphas[1]
            self.velocities = self.velocities[1]
            self.d_elevator = self.d_elevator[1]
            self.CL = self.CL[1]
            self.Cm = self.Cm[1]
            self.CDind = self.CDind[1]
            self.CD_tot = self.CD_tot[1]
            self.e = self.e[1]

        if self.phase == "climb":
            # picks best climb condition (max gamma)
            gamma_all = get_climb_gradient(
                self.velocities,
                self.S,
                self.CD_tot,
                rho=AircraftConfig.rho_climb,
                T=27496,  # TODO: update with correct climb thrust
                m=7504,  # TODO: update with correct climb thrust
                g=9.81,
            )
            self.gamma = gamma_all

            # idx = np.argmax(gamma_all)

            # self.alphas = self.alphas[idx]
            # self.velocities = self.velocities[idx]
            # self.d_elevator = self.d_elevator[idx]
            # self.CL = self.CL[idx]
            # self.Cm = self.Cm[idx]
            # self.CDind = self.CDind[idx]
            # self.CD_tot = self.CD_tot[idx]
            # self.e = self.e[idx]


if __name__ == "__main__":
    S_list = np.linspace(40, 52, 3)  # NOTE: See available options in runner script.
    # S_list = np.linspace(40, 52, 13)  # NOTE: See available options in runner script.
    blow_configs = ["full_blow", "half_blow"]

    for name in blow_configs:
        print(f"\n\n{'*' * 50}\n{'*' * 50}\n{name.upper()}")
        for idx, S in enumerate(S_list):
            print(f"\n{'=' * 40}\nS = {S} [m^2] Performance")
            for phase in ["takeoff", "cruise", "climb"]:
                config = AeroCoeffConfig(S=S, phase=phase, name=name)
                print(f"\n=== {phase.upper()} (S = {S} [m^2]) ===")

                print(f"Velocities: {np.round(config.velocities, 4)}")
                print(
                    f"AOAs: {np.round(config.alphas * (1 / DEG2RAD_CONV), 4)} [degrees]"
                )
                print(f"d_elevator: {np.round(config.d_elevator, 4)}")
                print(f"CL: {np.round(config.CL, 4)}")
                print(f"CD_tot: {np.round(config.CD_tot, 4)}")
                print(f"Cm: {np.round(config.Cm, 4)}")
                print(f"e: {np.round(config.e, 4)}")
                if phase in ("takeoff", "landing"):
                    print(f"xto: {np.round(config.xto, 4)} [ft]")

    # plt.figure()

    # for name in blow_configs:
    #     xto_vals = []

    #     for S in S_list:
    #         config = AeroCoeffConfig(S=S, phase="takeoff", name=name)
    #         xto_vals.append(config.xto)

    #     plt.plot(S_list, xto_vals, label=name, marker='o')

    # plt.xlabel("Wing Area S [m²]")
    # plt.ylabel("Runway Length [ft]")
    # plt.title("Takeoff Distance vs Wing Area")
    # plt.grid(True)
    # plt.legend()

    plt.figure()
    for name in blow_configs:
        ld_vals = []

        for S in S_list:
            config = AeroCoeffConfig(S=S, phase="takeoff", name=name)
            ld_vals.append(config.CL / config.CD_tot)

        plt.plot(S_list, ld_vals, label=name, marker="o")

    plt.xlabel("Wing Area S [m²]")
    plt.ylabel("L/D")
    plt.title("L/D vs Wing Area (Takeoff)")
    # plt.grid(True)
    plt.legend()

    REQUIRED_GAMMA = 0.286  # TODO: based on clearing 50 [ft] tree
    gamma_vals = {}
    cmap = plt.get_cmap("cividis")
    plt.figure()
    for S in S_list:
        config = AeroCoeffConfig(S=S, phase="climb", name="full_blow")
        alphas_deg = config.alphas / DEG2RAD_CONV

        unique_alphas = np.unique(alphas_deg)

        for i, alpha in enumerate(unique_alphas):
            mask = alphas_deg == alpha

            plt.scatter(
                [S] * np.sum(mask),
                config.gamma[mask],
                color=cmap(i / len(unique_alphas)),
                s=70,
                label=f"AoA = {alpha:.1f}°" if S == S_list[0] else None,
            )

    plt.axhline(
        REQUIRED_GAMMA, linestyle="--", color="red", label="Requirement (minimum)"
    )

    plt.xlabel("Wing Area S [m^2]")
    plt.ylabel("Climb Gradient")
    plt.title("Climb Gradient vs Wing Area")
    plt.legend()
    plt.show()

    # REQUIRED_XTO = 300  # ft
    # for name in blow_configs:
    #     plt.figure()
    #     for S in S_list:
    #         config = AeroCoeffConfig(S=S, phase="takeoff", name=name)
    #         plt.scatter(
    #             config.CL,
    #             config.xto,
    #             c=[S],
    #             cmap="magma",
    #             label=f"{name}, S={S}",
    #             vmin=40,
    #             vmax=60,
    #             s=70,
    #         )
    #     plt.xlabel("CL")
    #     plt.ylabel("Runway Length [ft]")
    #     plt.title("Runway Length vs CL")
    # plt.show()
