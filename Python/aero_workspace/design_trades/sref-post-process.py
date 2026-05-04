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

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from range_model import get_range

DEG2RAD_CONV = ureg("deg").to("rad").magnitude
METERS2FT_CONV = ureg("m").to("ft").magnitude
AR = 7
Thrust_takeoff = 50000  # 35000*0.8


def get_runway_length(
    v=22,
    S=50,
    CL=6.1,
    CD=1.88,
    mu=0.01,
    rho=AircraftConfig.rho_t0,
    T=Thrust_takeoff,  # 27496,
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


def T_climb(v=20, P_available_kW=1525, eta=0.8):
    P_W = P_available_kW * 1000.0
    return eta * P_W / v


def get_climb_gradient(
    vref,
    S,
    CL,
    CD,
    rho=AircraftConfig.rho_t0,
    m=7504,
    g=9.81,
):
    W = m * g
    q = 0.5 * rho * vref**2
    D = q * S * CD
    L = q * S * CL

    # Lift balance for physically consistent aircrafts only NOTE: using a loose tolerance atm
    mask = np.abs(L / W - 1) < 1e-2
    # mask = np.argmin(np.abs(L / W - 1))
    gamma = (T_climb(vref[mask]) - D[mask]) / W

    return mask, gamma


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
        if self.phase == "cruise":
            CD_DP = AircraftConfig.C_Dp_cruise
        else:
            CD_DP = (
                AircraftConfig.C_Dp_t0
            )  # Default for takeoff/landing; may need to modify in future iterations

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
                        v=run["velocity"], S=self.S, CL=run["CL"], CD=(CD_DP + _CDind)
                    )
                    * METERS2FT_CONV
                )

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
        self.alphas = self.alphas[elevator_mask]
        self.velocities = self.velocities[elevator_mask]
        self.d_elevator = self.d_elevator[elevator_mask]
        self.CL = self.CL[elevator_mask]
        self.Cm = self.Cm[elevator_mask]
        self.CDind = self.CDind[elevator_mask]
        self.e = self.e[elevator_mask]

        if self.phase == "cruise" and S in (49, 50):
            print(8)

        cl_best_real_mask = np.argmax(self.CL)
        cl_oper_real_mask = np.argwhere(
            (self.alphas / DEG2RAD_CONV == 10) & (self.velocities == 22)
        ).flatten()

        # produces best - for cl_max
        if self.phase == "takeoff":
            self.alphas_best = self.alphas[cl_best_real_mask]
            self.velocities_best = self.velocities[cl_best_real_mask]
            self.d_elevator_best = self.d_elevator[cl_best_real_mask]
            self.CL_best = self.CL[cl_best_real_mask]
            self.Cm_best = self.Cm[cl_best_real_mask]
            self.CDind_best = self.CDind[cl_best_real_mask]
            self.CD_tot_best = self.CD_tot[cl_best_real_mask]
            self.e_best = self.e[cl_best_real_mask]
            self.xto_best = np.array(self.xto)[cl_best_real_mask]

            self.alphas_oper = self.alphas[cl_oper_real_mask]
            self.velocities_oper = self.velocities[cl_oper_real_mask]
            self.d_elevator_oper = self.d_elevator[cl_oper_real_mask]
            self.CL_oper = self.CL[cl_oper_real_mask]
            self.Cm_oper = self.Cm[cl_oper_real_mask]
            self.CDind_oper = self.CDind[cl_oper_real_mask]
            self.CD_tot_oper = self.CD_tot[cl_oper_real_mask]
            self.e_oper = self.e[cl_oper_real_mask]
            self.xto_oper = np.array(self.xto)[cl_oper_real_mask]

        if self.phase == "climb":
            # picks best climb condition (max gamma)
            mask, gamma_all = get_climb_gradient(
                self.velocities,
                self.S,
                self.CL,
                self.CD_tot,
                rho=AircraftConfig.rho_climb,
                m=7504,  # TODO: update with correct mass
                g=9.81,
            )
            self.gamma = gamma_all

            # for valid climb angle
            self.alphas = self.alphas[mask]
            self.velocities = self.velocities[mask]
            self.d_elevator = self.d_elevator[mask]
            self.CL = self.CL[mask]
            self.Cm = self.Cm[mask]
            self.CDind = self.CDind[mask]
            self.CD_tot = self.CD_tot[mask]
            self.e = self.e[mask]

            # for steepest climb angle
            idx = np.argmax(gamma_all)
            self.alphas_best = self.alphas[idx]
            self.velocities_best = self.velocities[idx]
            self.d_elevator_best = self.d_elevator[idx]
            self.CL_best = self.CL[idx]
            self.Cm_best = self.Cm[idx]
            self.CDind_best = self.CDind[idx]
            self.CD_tot_best = self.CD_tot[idx]
            self.e_best = self.e[idx]
            self.gamma_best = self.gamma[idx]

        self.xto = np.array(self.xto)


if __name__ == "__main__":
    # S_list = np.linspace(40, 48, 9)  # NOTE: See available options in runner script.
    S_list = np.linspace(40, 52, 13)  # NOTE: See available options in runner script.
    # blow_configs = ["full_blow"]
    blow_configs = ["full_blow", "half_blow"]

    # for name in blow_configs:
    #     print(f"\n\n{'*' * 50}\n{'*' * 50}\n{name.upper()}")
    #     for idx, S in enumerate(S_list):
    #         print(f"\n{'=' * 40}\nS = {S} [m^2] Performance")
    #         # for phase in ["climb"]:
    #         for phase in ["takeoff", "cruise", "climb"]:
    #             config = AeroCoeffConfig(S=S, phase=phase, name=name)
    #             print(f"\n=== {phase.upper()} (S = {S} [m^2]) ===")

    #             print(f"Velocities: {np.round(config.velocities, 4)}")
    #             print(f"AOAs: {np.round(config.alphas * (1 / DEG2RAD_CONV), 4)} [°]")
    #             print(f"d_elevator: {np.round(config.d_elevator, 4)}")
    #             print(f"CL: {np.round(config.CL, 4)}")
    #             print(f"CD_tot: {np.round(config.CD_tot, 4)}")
    #             print(f"Cm: {np.round(config.Cm, 4)}")
    #             print(f"e: {np.round(config.e, 4)}")
    #             if phase in ("takeoff", "landing"):
    #                 print(f"xto: {np.round(config.xto, 4)} [ft]")

    # # # # ####################################################################
    # # L/D @Takeoff Trade
    # plt.figure(figsize=(8, 6))
    # for name in blow_configs:
    #     ld_vals = []

    #     for S in S_list:
    #         config = AeroCoeffConfig(S=S, phase="takeoff", name=name)
    #         ld_vals.append(config.CL_best / config.CD_tot_best)

    #     plt.plot(
    #         S_list, ld_vals, label=({name.replace("_", " ").capitalize()}), marker="o"
    #     )

    # # plt.xlabel("Wing Area, S [m^2]")
    # plt.xlabel("Wing Area, S [m^2]")
    # plt.ylabel("L/D")
    # plt.title("L/D vs Wing Area (Takeoff)")
    # # plt.grid(True)
    # plt.legend()
    # plt.show()
    # # # ####################################################################

    # # # ####################################################################
    # # Climb Trade
    # REQUIRED_GAMMA = 10  # TODO: based on clearing 50 [ft] tree
    # gamma_vals = {}
    # cmap = plt.get_cmap("cividis")
    # plt.figure(figsize=(8, 6))
    # for S in S_list:
    #     config = AeroCoeffConfig(S=S, phase="climb", name="full_blow")
    #     alphas_deg = config.alphas / DEG2RAD_CONV

    #     unique_alphas = np.unique(alphas_deg)

    #     for i, alpha in enumerate(unique_alphas):
    #         mask = alphas_deg == alpha

    #         plt.scatter(
    #             [S] * np.sum(mask),
    #             config.gamma[mask] / DEG2RAD_CONV,
    #             color=cmap(i / len(unique_alphas)),
    #             s=70,
    #             label=f"AoA = {alpha:.1f}°, V = {config.velocities[mask][0]:.1f}"
    #             if S == S_list[0]
    #             else None,
    #         )
    # plt.axhline(REQUIRED_GAMMA, linestyle="--", color="red", label="Requirement")
    # plt.xlabel("Wing Area, S [m^2]")
    # plt.ylabel("Climb Gradient [°]")
    # plt.title("Climb Gradient vs Wing Area")
    # plt.legend()
    # plt.show()
    # # # ####################################################################

    # # ####################################################################
    # # Range Trade
    # cmap = plt.get_cmap("cividis")
    # plt.figure(figsize=(8, 6))

    # S_vals = []
    # range_vals_plot = []
    # colors = []

    # for i, S in enumerate(S_list):
    #     color = cmap(i / len(S_list))
    #     config = AeroCoeffConfig(S=S, phase="cruise", name="full_blow")
    #     config_clmax = AeroCoeffConfig(S=S, phase="takeoff", name="full_blow")
    #     config_climb = AeroCoeffConfig(S=S, phase="climb", name="full_blow")

    #     range_vals = get_range(
    #         mass=7600,
    #         Cd0=AircraftConfig.C_Dp_cruise,
    #         Cdi_cruise=config.CDind,
    #         CLmax=config_clmax.CL_best,
    #         AR=AR,
    #         wing_area=S,
    #         v_cruise=config.velocities,
    #         takeoff_power=2500.0,
    #         climb_angle=10,#config_climb.gamma[0],
    #         # climb_angle=10*DEG2RAD_CONV,#max(config_climb.gamma),
    #         cruise_altitude=AircraftConfig.h_cruise,
    #     )

    #     cruise_range = range_vals["cruise_range_km"]
    #     max_range = np.max(cruise_range)
    #     S_vals.append(S)
    #     range_vals_plot.append(max_range)
    #     colors.append(color)

    # plt.bar([str(S) for S in S_vals], range_vals_plot, color=colors)
    # plt.xlabel("Wing Area, S [m^2]")
    # plt.ylabel("Max Cruise Range (km)")
    # plt.title("Max Range vs Wing Area")
    # plt.show()

    # # ####################################################################
    # # x_TO & CL
    # REQUIRED_XTO = 300  # ft
    # # for name in blow_configs:
    # for name in ["full_blow"]:
    #     plt.figure(figsize=(8, 6))
    #     for S in S_list:
    #         config = AeroCoeffConfig(S=S, phase="takeoff", name=name)
    #         sc = plt.scatter(
    #             config.alphas / DEG2RAD_CONV,
    #             config.xto,
    #             c=[S] * len(config.alphas),
    #             cmap="plasma",
    #             vmin=min(S_list),
    #             vmax=max(S_list),
    #             s=70,
    #         )

    #     plt.xlabel("CL")
    #     plt.ylabel("Runway Length [ft]")
    #     plt.title(f"Runway Length vs CL - ({name.replace('_', ' ').capitalize()})")
    #     cbar = plt.colorbar(sc)
    #     cbar.set_label("Wing Area, S [m^2]")

    # plt.show()

    # # ####################################################################
    # # x_TO Trade w/o CL
    # plt.figure(figsize=(8, 6))

    # for name in ["full_blow"]:
    #     xto_vals = []
    #     CL_vals = []
    #     alphas = []

    #     for S in S_list:
    #         config = AeroCoeffConfig(S=S, phase="takeoff", name=name)

    #         if len(config.xto_oper) == 0:
    #             xto_vals.append(np.nan)
    #             CL_vals.append(np.nan)
    #         else:
    #             xto_vals.append(config.xto_oper[0])# * 0.33)
    #             CL_vals.append(config.CL_oper[0])  # <-- key addition
    #             alphas.append(config.alphas_oper[0])  # <-- key addition

    #     xto_vals = np.array(xto_vals)
    #     CL_vals = np.array(CL_vals)
    #     alphas = np.array(alphas)

    #     # Base line
    #     plt.plot(S_list, xto_vals, color="black", linewidth=1.0)#, alpha=0.7)

    #     # Color band via scatter
    #     sc = plt.scatter(
    #         S_list,
    #         xto_vals,
    #         c=CL_vals,
    #         cmap="viridis",
    #         s=90, #vmin=5.4,vmax=6.8,
    #         # edgecolor="k",
    #     )

    # plt.xlabel("Wing Area, S [m^2]")
    # plt.ylabel("Takeoff Length [ft]")
    # plt.title("Takeoff Distance vs Wing Area (Full blow)")

    # cbar = plt.colorbar(sc)
    # cbar.set_label("CL")
    # plt.show()
    # # # Multiple oper points
    # # # x_TO Trade w/o CL
    # # plt.figure(figsize=(8, 6))

    # # name = "full_blow"
    # # config_ref = AeroCoeffConfig(S=S_list[0], phase="takeoff", name=name)

    # # alphas = config_ref.alphas / DEG2RAD_CONV  # deg
    # # velocities = config_ref.velocities

    # # for k in range(len(alphas)):
    # #     xto_curve = []

    # #     for S in S_list:
    # #         config = AeroCoeffConfig(S=S, phase="takeoff", name=name)

    # #         # assume SAME ordering of runs for each S
    # #         xto_curve.append(config.xto[k])

    # #     plt.plot(
    # #         S_list,
    # #         xto_curve,
    # #         marker="o",
    # #         label=f"α={alphas[k]:.1f}°, V={velocities[k]:.1f} m/s",
    # #     )

    # # plt.xlabel("Wing Area, S [m^2]")
    # # plt.ylabel("Runway Length [ft]")
    # # plt.title("Takeoff Distance vs Wing Area")
    # # # plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")  # put legend on side
    # # # plt.tight_layout()

    # # plt.show()

    # -------------------------------------------------
    # Pareto front: Range vs Takeoff Distance
    # -------------------------------------------------
    S_vals = []
    range_vals = []
    xto_vals = []

    for S in S_list:
        # configs
        config_cruise = AeroCoeffConfig(S=S, phase="cruise", name="full_blow")
        config_to = AeroCoeffConfig(S=S, phase="takeoff", name="full_blow")
        config_climb = AeroCoeffConfig(S=S, phase="climb", name="full_blow")

        # range
        range_out = get_range(
            mass=7500,
            Cd0=AircraftConfig.C_Dp_cruise,
            Cdi_cruise=config_cruise.CDind,
            CLmax=config_to.CL_best,
            AR=AR,
            wing_area=S,
            v_cruise=config_cruise.velocities,
            takeoff_power=2500.0,
            climb_angle=15,
            cruise_altitude=AircraftConfig.h_cruise,
        )

        max_range = np.max(range_out["cruise_range_km"])*0.75

        # takeoff distance (operating point)
        if len(config_to.xto_oper) == 0:
            continue

        xto = config_to.xto_oper[0]

        S_vals.append(S)
        range_vals.append(max_range)
        xto_vals.append(xto)

    range_vals = np.array(range_vals)
    xto_vals = np.array(xto_vals)*0.92
    S_vals = np.array(S_vals)

    # -------------------------------------------------
    # Pareto filter (maximize range, minimize xto)
    # -------------------------------------------------
    pareto_mask = np.ones(len(range_vals), dtype=bool)

    for i in range(len(range_vals)):
        for j in range(len(range_vals)):
            if (range_vals[j] >= range_vals[i] and xto_vals[j] <= xto_vals[i]) and (
                range_vals[j] > range_vals[i] or xto_vals[j] < xto_vals[i]
            ):
                pareto_mask[i] = False
                break

    # -------------------------------------------------
    # Plot
    # -------------------------------------------------
    plt.figure(figsize=(8, 6))

    # all designs (faded)
    plt.scatter(xto_vals, range_vals, c=S_vals, cmap="viridis", alpha=0.4, s=60)

    # pareto front (bold)
    pareto_sorted = np.argsort(xto_vals[pareto_mask])
    plt.plot(
        xto_vals[pareto_mask][pareto_sorted],
        range_vals[pareto_mask][pareto_sorted],
        color="red",
        linewidth=2.5,
        label="Pareto front",
    )

    plt.scatter(
        xto_vals[pareto_mask],
        range_vals[pareto_mask],
        c=S_vals[pareto_mask],
        cmap="viridis",
        edgecolor="k",
        s=100,
    )

    cbar = plt.colorbar()
    cbar.set_label("Wing Area S [m²]")

    plt.xlabel("Takeoff Distance [ft]")
    plt.ylabel("Cruise Range [km]")
    plt.title("Design Trade: Range vs Takeoff Distance")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.show()