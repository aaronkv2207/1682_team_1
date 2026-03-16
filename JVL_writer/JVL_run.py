from typing import List

import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.tools.pretty_plots as p
from aerosandbox.tools import units
from geom import planeB as plane  # , jvl_plane
from J import JVL, JetControl, JetParam, JWing, WingJSec
from scipy.interpolate import Akima1DInterpolator

# jvl_plane.run()

jvl_plane = JVL(
    airplane=plane,
    op_point=asb.OperatingPoint(
        velocity=100,
        alpha=5,
        beta=0,
        p=0,
        q=0,
        r=0,
    ),
    avl_command="jvl",  # NOTE: may change depending on your executable
)
jvl_plane.default_analysis_specific_options = {
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

result = jvl_plane.run()  # dictionary of results
print(result)
