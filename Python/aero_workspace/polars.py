import matplotlib.pyplot as plt
import numpy as np
from aero_dict import AircraftConfig, AircraftConfig2
from aero_main import DEG2RAD_CONV, AeroCoeffConfig
from matplotlib import cm, colors


def plot_drag_polars(config: AeroCoeffConfig, title: str):
    # NOTE: some simulated points are only for plotting purposes

    unique_vels = np.unique(config.velocities)
    chosen_flap_deflection = 50.0

    fig, ax = plt.subplots()

    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    all_alphas = config.alphas / DEG2RAD_CONV
    norm = colors.Normalize(vmin=np.min(all_alphas), vmax=np.max(all_alphas))
    cmap = cm.viridis

    for i, V in enumerate(unique_vels):
        mask = np.isclose(config.velocities, V)

        if config.phase == "takeoff":
            mask &= np.isclose(config.flap_1, chosen_flap_deflection)

        CL = config.CL[mask]
        CD = config.CD_tot[mask].flatten()
        AoA = config.alphas[mask] / DEG2RAD_CONV

        if len(CL) < 3:
            continue

        sort_idx = np.argsort(CL)
        CL = CL[sort_idx]
        CD = CD[sort_idx]
        AoA = AoA[sort_idx]

        color = color_cycle[i % len(color_cycle)]

        # quadratic drag polar fit
        CL_smooth = np.linspace(CL.min(), CL.max(), 400)
        coeffs = np.polyfit(CL, CD, 2)
        CD_smooth = np.polyval(coeffs, CL_smooth)
        ax.plot(
            CD_smooth,
            CL_smooth,
            color=color,
            linewidth=2,
            label=f"V = {V:.1f} m/s",
            zorder=2,
        )
        sc = ax.scatter(
            CD, CL, c=AoA, cmap=cmap, norm=norm, s=30, edgecolors="none", zorder=3
        )

    ax.set_xlabel("CD")
    ax.set_ylabel("CL")
    ax.set_title(title)
    ax.legend()
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.02)
    cbar.set_label("Angle of Attack (deg)")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    for plane_idx, plane in enumerate([AircraftConfig, AircraftConfig2]):
        S = plane.S

        for phase in ["takeoff", "cruise"]:
            config = AeroCoeffConfig(phase=phase, aircraft=plane)

            plot_drag_polars(
                config,
                title=f"{phase.capitalize()} Drag Polar (Plane {plane_idx + 1}, S={S})",
            )
