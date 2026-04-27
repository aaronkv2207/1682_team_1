"""
Main Run File. File is structured as follows:
    1. Plane config of interest. (Can edit control surface deflections in geom.py)
    2. Operating conditions
    3. Parameter sweep
    4. Run jvl and save as dataframes
    NOTE:   Only use this script to produce new dataframes; otherwise use current run data in aero_main,
            or run output folder. All post-processing should be done outside of this script.
"""

import os
import pickle
import sys
from dataclasses import dataclass
from itertools import product

import aerosandbox as asb
import aerosandbox.numpy as np
from geom import analysis_options, planeB
from J import JVL, JWing, WingJSec

DEPRECATION_MESSAGE = "DEPRECATED! You will need to update script logic to sweep over deflections, creating different airplane configs instead of importing planeB. "


# planes/control-surface definitions
def create_plane(cond):
    # assumed for entire plane at the moment
    jvl = JVL(
        airplane=planeB,  # TODO: Update with desired plane config. Assumes plane geometry made in geom.py file
        op_point=asb.OperatingPoint(
            velocity=cond["velocity"],
            alpha=cond["alpha"],
            p=0,
            q=0,
            r=0,
        ),
        avl_command="jvl",
    )

    jvl.default_analysis_specific_options = analysis_options
    return jvl


def run_case(cond):
    surface = create_plane(cond)

    out = surface.run()

    return {
        "velocity": cond["velocity"],
        "alpha": cond["alpha"],
        **out,
    }


def save_results(data, filename):
    """Save data as dataframe."""
    with open(filename, "wb") as f:
        pickle.dump(data, f)


def takeoff():
    alphas = np.linspace(8, 20, 1)
    velocities = np.linspace(18, 25, 3)
    cases = [{"velocity": v, "alpha": a} for v, a in product(velocities, alphas)]
    results = []

    for case_idx, case in enumerate(cases):
        status = run_case(case)
        print(f"Successfully ran case {case_idx + 1} of {len(cases)} in takeoff")
        results.append(status)
    save_results(results, "JVL_writer/sref_design-trades/run_file_ouputs/coefficient_results/takeoff.pkl")


def cruise():
    # TODO: Define operating points
    alphas = np.array([0])
    velocities = np.array([80, 125, 150])
    cases = [{"velocity": v, "alpha": a} for v, a in product(velocities, alphas)]
    results = []

    for case_idx, case in enumerate(cases):
        status = run_case(case)
        print(f"Successfully ran case {case_idx + 1} of {len(cases)} in cruise")
        results.append(status)
    save_results(results, "JVL_writer/sref_design-trades/run_file_ouputs/coefficient_results/cruise.pkl")


def landing():
    # TODO: Define operating points
    alphas = -np.linspace(8, 20, 1)
    velocities = np.linspace(18, 25, 3)
    cases = [{"velocity": v, "alpha": a} for v, a in product(velocities, alphas)]
    results = []

    for case_idx, case in enumerate(cases):
        status = run_case(case)
        print(f"Successfully ran case {case_idx + 1} of {len(cases)} in landing")
        results.append(status)
    save_results(results, "JVL_writer/sref_design-trades/run_file_ouputs/coefficient_results/landing.pkl")


if __name__ == "__main__":
    # sys.exit(DEPRECATION_MESSAGE)
    takeoff()
    cruise()
    landing()
