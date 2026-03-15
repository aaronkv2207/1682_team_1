"""
Main Run File. File is structured as follows:
    1. Plane config of interest. (Can edit control surface deflections in geom.py)
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
from geom import analysis_options, planeB
from J import JVL, JWing, WingJSec


# planes/control-surface definitions
def create_plane(cond):
    # assumed for entire plane at the moment
    jvl = JVL(
        airplane=planeB,  # TODO: Update with desired plane config. Assumes plane geometry made in geom.py file
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

    jvl.default_analysis_specific_options = analysis_options
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

    results = []

    for case_idx, case in enumerate(cases):
        status = run_case(case)
        print(f"Successfully ran case {case_idx + 1} of {len(cases)} in takeoff")
        results.append(status)
        save_results(results, "Python/aero_workspace/jvl_run_outputs/takeoff.pkl")


def climb():
    # TODO: Define operating points
    alphas = np.linspace(8, 30, 1)
    betas = np.linspace(9, 30, 1)
    velocities = np.linspace(120, 150, 1)

    cases = [
        {"velocity": v, "alpha": a, "beta": b}
        for v, a, b in product(velocities, alphas, betas)
    ]

    results = []

    for case_idx, case in enumerate(cases):
        status = run_case(case)
        print(f"Successfully ran case {case_idx + 1} of {len(cases)} in climb")
        results.append(status)
        save_results(results, "Python/aero_workspace/jvl_run_outputs/climb.pkl")


def cruise():
    # TODO: Define operating points
    alphas = np.linspace(8, 30, 1)
    betas = np.linspace(9, 30, 1)
    velocities = np.linspace(120, 150, 1)

    cases = [
        {"velocity": v, "alpha": a, "beta": b}
        for v, a, b in product(velocities, alphas, betas)
    ]

    results = []

    for case_idx, case in enumerate(cases):
        status = run_case(case)
        print(f"Successfully ran case {case_idx + 1} of {len(cases)} in cruise")
        results.append(status)
        save_results(results, "Python/aero_workspace/jvl_run_outputs/cruise.pkl")


def landing():
    # TODO: Define operating points
    alphas = np.linspace(8, 30, 1)
    betas = np.linspace(9, 30, 1)
    velocities = np.linspace(120, 150, 1)

    cases = [
        {"velocity": v, "alpha": a, "beta": b}
        for v, a, b in product(velocities, alphas, betas)
    ]

    results = []

    for case_idx, case in enumerate(cases):
        status = run_case(case)
        print(f"Successfully ran case {case_idx + 1} of {len(cases)} in landing")
        results.append(status)
        save_results(results, "Python/aero_workspace/jvl_run_outputs/landing.pkl")


if __name__ == "__main__":
    takeoff()
    climb()
    cruise()
    landing()
