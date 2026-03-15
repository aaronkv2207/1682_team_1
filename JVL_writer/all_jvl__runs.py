"""
Main Run File. File is structured as follows:
    1. Geometries & control-surfaces of interest
    2. Operating conditions
    3. Parameter sweep
    4. Run jvl and save as dataframes
    NOTE:   Only use this script to produce new dataframes; otherwise use current run data in aero_main,
            or run output folder. All post-processing should be done outside of this script.
"""

import pickle
from dataclasses import dataclass
from itertools import product

import aerosandbox as asb
import aerosandbox.numpy as np
from geom import planeB
from J import JVL, JWing, WingJSec

# planes/control-surface definitions
DEFAULT_ANALYSIS_OPTIONS_PLANE = {
    asb.Airplane: dict(profile_drag_coefficient=0),
    JWing: dict(
        wing_level_spanwise_spacing=True,
        spanwise_resolution=25,
        spanwise_spacing="cosine",
        chordwise_resolution=16,
        chordwise_spacing="cosine",
        component=None,  # This is an int
        no_wake=False,
        no_alpha_beta=False,
        no_load=False,
        drag_polar=None,
    ),
    asb.Wing: dict(
        wing_level_spanwise_spacing=True,
        spanwise_resolution=10,
        spanwise_spacing="cosine",
        chordwise_resolution=8,
        chordwise_spacing="cosine",
        component=None,  # This is an int
        no_wake=False,
        no_alpha_beta=False,
        no_load=False,
        drag_polar=None,
    ),
    WingJSec: dict(
        spanwise_resolution=12,
        spanwise_spacing="cosine",
        cl_alpha_factor=None,  # This is a float
        drag_polar=None,
    ),
    # asb.Fuselage: dict(panel_resolution=24, panel_spacing="cosine"), # NOTE: fuselage creation doesn't work on Mac silicon jvl version
}


def create_plane(cond):
    # assumed for entire plane at the moment
    jvl = JVL(
        airplane=planeB,  # NOTE: assumes plane geometry made in geom.py file
        op_point=asb.OperatingPoint(
            velocity=cond["velocity"],
            alpha=cond["alpha"],
            beta=cond["beta"],
            p=0,
            q=0,
            r=0,
        ),
        avl_command="jvl",
    )

    jvl.default_analysis_specific_options = DEFAULT_ANALYSIS_OPTIONS_PLANE
    return jvl


def run_case(cond):
    surface = create_plane(cond)

    out = surface.run()

    return {
        "velocity": cond["velocity"],
        "alpha": cond["alpha"],
        "beta": cond["beta"],
        **out,
    }


def save_results(data, filename):
    """Save data as dataframe."""
    with open(filename, "wb") as f:
        pickle.dump(data, f)


def takeoff():
    alphas = np.linspace(8, 30, 1)
    betas = np.linspace(9, 30, 1)
    velocities = np.linspace(120, 150, 1)

    cases = [
        {"velocity": v, "alpha": a, "beta": b}
        for v, a, b in product(velocities, alphas, betas)
    ]

    results = [run_case(case) for case in cases]

    save_results(results, "Python/aero_workspace/jvl_run_outputs/takeoff.pkl")


if __name__ == "__main__":
    takeoff()
