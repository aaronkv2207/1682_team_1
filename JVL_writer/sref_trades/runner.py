import os
import pickle
import sys
from itertools import product

import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.tools.pretty_plots as p
from aerosandbox import Atmosphere
from aerosandbox.tools import units
from scipy.interpolate import Akima1DInterpolator

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from J import JVL, JetControl, JetParam, JWing, WingJSec

# assumed cruise height
h_cruise = 5486.4  # m
h_climb = 3135  # TODO: completely arbitrary --> needs correction

takeoff_deflection = 50  # degrees
climb_deflection = 30
landing_deflection = 65
trim_variable = "d6"  # elevator

# TODO: Fix inner and outer flap spans

# fuselage (not very well parametrized at the moment)
nose_x = -5
fuse_width = 1.6  # width at the wing: 2 rows of seats
AR = 8  # twin otter is 10.05
ail_hinge = 0.7
flap_hinge = 0.75
wing_incidence = 0
aileron_fraction = 0.75

# vertical tail geometry - fixed parameters
lv = 8  # distance from wing quarter chord to vertical tail quarter chord
Vv = 0.06  # Vertical tail volume coefficient
vt_ar = 1.2
tail_hinge = 0.7
Vh = 0.75

# h_tail - fixed parameters
ht_ar = 2

# fan - fixed parameters
fan_radius = 1.166 / 2
n_fans = 8

# Wing geometry
main_foil = asb.Airfoil(coordinates="./JVL_writer/jw05.dat")
# S = 49.6  # twin otter wing area is 39 m^2


def plane_operating_point(cond, plane, analysis_options):
    jvl = JVL(
        airplane=plane,
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


def generate_fuselage_xsecs(
    N: int, vt_x, vt_sweep_x: float, vt_c: float
) -> list[asb.FuselageXSec]:
    """
    Generates a fuselage with N sections, transitioning from a small circular nose,
    to a large square-like midsection, and tapering into a smaller square tail.

    Args:
        N (int): Number of sections defining the fuselage.

    Returns:
        List[FuselageXSec]: List of fuselage cross-sections.
    """
    fuse_end = (
        vt_x - vt_sweep_x + vt_c
    )  # end of the fuselage, at the end of the vertical tail
    x_positions = np.linspace(
        nose_x, fuse_end, N
    )  # Generate N sections along the fuselage
    # print(nose_x, fuse_end)
    r = np.array(
        [
            0.3,
            2 / 3 * fuse_width,
            fuse_width,
            fuse_width,
            fuse_width,
            3 / 4 * fuse_width,
            0.5,
        ]
    )
    z = -r / 2
    zs = np.interp(
        x_positions,
        [nose_x, 7.8 * nose_x, nose_x / 2, 0, vt_x / 2, vt_x, fuse_end],
        z,
    )  # Interpolate z position
    radius = np.interp(
        x_positions,
        [nose_x, 7 / 8 * nose_x, 1 / 2 * nose_x, 0, vt_x / 2, vt_x, fuse_end],
        r,
    )  # Width transition
    shapes = np.interp(
        x_positions,
        [nose_x, nose_x / 2, vt_x / 2, vt_x, vt_x + vt_c],
        [1, 2.5, 2.5, 2, 1],
    )  # Shape transition

    xsecs = []
    for x, z, r, shape in zip(x_positions, zs, radius, shapes):
        xsecs.append(
            asb.FuselageXSec(
                xyz_c=[x, 0, z],
                xyz_normal=[1, 0, 0],  # Assume fuselage is aligned along x-axis
                width=r,
                height=1.1 * r,
                shape=2.3,
            )
        )

    return xsecs


def save_results(data, filename):
    """Save data as dataframe."""
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    with open(filename, "wb") as f:
        pickle.dump(data, f)


def run_case(phase, surface, velocity, alpha, blowing_perturb):
    if phase == "takeoff":
        out = surface.run(
            flap_deflections={"d1": takeoff_deflection, "d2": takeoff_deflection},
            trim_Cm_to_zero=True,
            trim_variable=trim_variable,
            blowing={
                "J1": 2.0 * blowing_perturb,
                "J2": 2.0 * blowing_perturb,
                "J3": 2.0 * blowing_perturb,
                "J4": 2.0 * blowing_perturb,
            },
        )
    elif phase == "climb":
        out = surface.run(
            flap_deflections={"d1": climb_deflection, "d2": climb_deflection},
            trim_Cm_to_zero=True,
            trim_variable=trim_variable,
            blowing={
                "J1": 1.0 * blowing_perturb,
                "J2": 1.0 * blowing_perturb,
                "J3": 1.0 * blowing_perturb,
                "J4": 1.0 * blowing_perturb,
            },
        )
    elif phase == "cruise":
        out = surface.run(
            run_command=None,
            trim_Cm_to_zero=True,
            trim_variable=trim_variable,
        )
    else:  # landing
        out = surface.run(
            run_command=None,
            trim_Cm_to_zero=True,
            trim_variable=trim_variable,
            flap_deflections={"d1": landing_deflection, "d2": landing_deflection},
            blowing={
                "J1": 1.5 * blowing_perturb,
                "J2": 1.5 * blowing_perturb,
                "J3": 1.5 * blowing_perturb,
                "J4": 1.5 * blowing_perturb,
            },
        )

    return {
        "velocity": velocity,
        "alpha": alpha,
        "trim-var (CHECK-FEASIBILITY)": trim_variable,
        **out,
    }


# NOTE: I am extending on planeB, since I am not looking at differential blowing
def run_sref_cases(S_list, oper_dict, folder_name, blowing_perturb):  # noqa: PLR0915
    for S in S_list:
        for phase in oper_dict:
            velocities, alphas = (
                oper_dict[phase]["velocities"],
                oper_dict[phase]["alphas"],
            )
            cases = [
                {"velocity": v, "alpha": a} for v, a in product(velocities, alphas)
            ]
            results = []
            for case_idx, case in enumerate(cases):
                velocity, alpha = case.values()
                config = f"{phase}_Sref_{S}_v_{int(velocity)}_aoa_{int(alpha)}"
                b = np.sqrt(AR * S)
                aileron_y = (
                    aileron_fraction * b / 2
                )  # start aileron at 2/3 of halfspan. No blowing at this point.
                MAC = b / AR
                tip_chord = 9 / 10 * MAC  # set taper ratio if desired
                root_chord = MAC  # (2*MAC - tip_chord) #root chord needs to be adjusted for lost area from taper at the wingtips #TODO: fix

                blown_span = (aileron_fraction) * b - fuse_width

                fan_length = n_fans * 2 * fan_radius
                disc_area = n_fans * np.pi * fan_radius**2
                hdisk = disc_area / blown_span
                blowing_dy = blown_span / n_fans

                Sv = Vv * S * b / lv
                vt_c = np.sqrt(Sv / vt_ar)  # constant chord
                vt_b = vt_ar * vt_c
                vt_x = (
                    lv + 1 / 4 * MAC - 1 / 4 * vt_c
                )  # distance from wing LE to tail (avg) LE

                # horizontal tail geometry
                lh = (
                    lv + 0.25 * vt_c
                )  # include some sweep in the vertical tail to extend the ht moment arm

                Sh = Vh * S * MAC / lh
                ht_mac = np.sqrt(
                    Sh / ht_ar
                )  # technically not mac but geometry mean chord, but whatever
                ht_rc = vt_c
                ht_tc = 2 * ht_mac - ht_rc
                ht_b = ht_ar * ht_mac

                ht_x = lh + 1 / 4 * MAC - 1 / 4 * ht_mac
                vt_sweep_x = ht_x - vt_x

                if False:
                    print(f"Vertical tail area: {Sv:.2f} m^2")
                    print(f"Vertical tail chord: {vt_c:.2f} m")
                    print(f"Wing span: {b:.2f} m")
                    print(f"Wing aspect ratio: {AR:.2f}")
                    print(f"Wing mean aerodynamic chord: {MAC:.2f} m")
                    print(f"blown_span: {blown_span: .2f} m")
                    print(f"fan_length: {fan_length: .2f} m")
                    print(f"hdisk: {hdisk: .2f} m")
                    print(f"blowing_dy: {blowing_dy: .2f} m")
                    print(f"Horizontal tail area: {Sh:.2f} m^2")
                    print(f"Horizontal tail MA chord: {ht_mac:.2f} m")

                wingB_unblown = JWing(
                    name="Main Wing",
                    symmetric=True,
                    JetParam=JetParam(
                        hdisk=hdisk,
                        fh=0.0,
                        djet0=0.0,
                        djet1=0.0,
                        djet3=0.0,
                        dxdsk=0.07,
                        dndsk=-0.51 * hdisk,
                    ),
                    xsecs=[
                        WingJSec(
                            xyz_le=[0, 0, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1", hinge_point=flap_hinge, deflection=0
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet1", gain=0, sgn_dup=1)
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, fuse_width / 2 + blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1", hinge_point=flap_hinge, deflection=0
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet1", gain=0, sgn_dup=1),
                                JetControl(jet_name="FlapJet2", gain=0, sgn_dup=1),
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, fuse_width / 2 + 2 * blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1", hinge_point=flap_hinge, deflection=0
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet2", gain=0, sgn_dup=1),
                                JetControl(jet_name="FlapJet3", gain=0, sgn_dup=1),
                            ],
                        ),
                        WingJSec(  # flap1 goes to end of 3rd duct
                            xyz_le=[0, fuse_width / 2 + 3 * blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1", hinge_point=flap_hinge, deflection=0
                                ),
                                asb.ControlSurface(
                                    name="Flap2", hinge_point=flap_hinge, deflection=0
                                ),
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet3", gain=0, sgn_dup=1),
                                JetControl(jet_name="FlapJet4", gain=0, sgn_dup=1),
                            ],
                        ),
                        # WingJSec(
                        #     xyz_le=[0, fuse_width/2 + 4*blowing_dy, 0],
                        #     chord=root_chord,
                        #     twist = wing_incidence,
                        #     airfoil=main_foil,
                        #     control_surfaces = [asb.ControlSurface(name="Flap2", hinge_point=flap_hinge, deflection=0)],
                        #     JetControls= [JetControl(jet_name="FlapJet4", gain=1, sgn_dup=1), JetControl(jet_name="FlapJet5", gain=1, sgn_dup=1)],
                        # ),
                        WingJSec(
                            xyz_le=[0, aileron_y, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap2", hinge_point=0.8, deflection=0
                                ),
                                asb.ControlSurface(
                                    name="Aileron",
                                    symmetric=False,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                                asb.ControlSurface(
                                    name="Flaperon",
                                    symmetric=True,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet4", gain=0, sgn_dup=1)
                            ],
                        ),
                        WingJSec(
                            xyz_le=[
                                ail_hinge * root_chord - ail_hinge * tip_chord,
                                1 / 2 * b,
                                0,
                            ],
                            chord=tip_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Aileron",
                                    symmetric=False,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                                asb.ControlSurface(
                                    name="Flaperon",
                                    symmetric=True,
                                    hinge_point=0.8,
                                    deflection=0,
                                ),
                            ],
                        ),
                    ],
                )

                wingB_takeoff = JWing(
                    name="Main Wing",
                    symmetric=True,
                    JetParam=JetParam(
                        hdisk=hdisk,
                        fh=0.0,
                        djet0=0.0,
                        djet1=0.0,
                        djet3=0.0,
                        dxdsk=0.07,
                        dndsk=-0.51 * hdisk,
                    ),
                    xsecs=[
                        WingJSec(
                            xyz_le=[0, 0, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=takeoff_deflection,
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1)
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, fuse_width / 2 + blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=takeoff_deflection,
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1),
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, fuse_width / 2 + 2 * blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=takeoff_deflection,
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet3", gain=1, sgn_dup=1),
                            ],
                        ),
                        WingJSec(  # flap1 goes to end of 3rd duct
                            xyz_le=[0, fuse_width / 2 + 3 * blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=takeoff_deflection,
                                ),
                                asb.ControlSurface(
                                    name="Flap2",
                                    hinge_point=flap_hinge,
                                    deflection=takeoff_deflection,
                                ),
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet3", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet4", gain=1, sgn_dup=1),
                            ],
                        ),
                        # WingJSec(
                        #     xyz_le=[0, fuse_width/2 + 4*blowing_dy, 0],
                        #     chord=root_chord,
                        #     twist = wing_incidence,
                        #     airfoil=main_foil,
                        #     control_surfaces = [asb.ControlSurface(name="Flap2", hinge_point=flap_hinge, deflection=0)],
                        #     JetControls= [JetControl(jet_name="FlapJet4", gain=1, sgn_dup=1), JetControl(jet_name="FlapJet5", gain=1, sgn_dup=1)],
                        # ),
                        WingJSec(
                            xyz_le=[0, aileron_y, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap2",
                                    hinge_point=0.8,
                                    deflection=takeoff_deflection,
                                ),
                                asb.ControlSurface(
                                    name="Aileron",
                                    symmetric=False,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                                asb.ControlSurface(
                                    name="Flaperon",
                                    symmetric=True,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet4", gain=1, sgn_dup=1)
                            ],
                        ),
                        WingJSec(
                            xyz_le=[
                                ail_hinge * root_chord - ail_hinge * tip_chord,
                                1 / 2 * b,
                                0,
                            ],
                            chord=tip_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Aileron",
                                    symmetric=False,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                                asb.ControlSurface(
                                    name="Flaperon",
                                    symmetric=True,
                                    hinge_point=0.8,
                                    deflection=0,
                                ),
                            ],
                        ),
                    ],
                )

                wingB_landing = JWing(
                    name="Main Wing",
                    symmetric=True,
                    JetParam=JetParam(
                        hdisk=hdisk,
                        fh=0.0,
                        djet0=0.0,
                        djet1=0.0,
                        djet3=0.0,
                        dxdsk=0.07,
                        dndsk=-0.51 * hdisk,
                    ),
                    xsecs=[
                        WingJSec(
                            xyz_le=[0, 0, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=landing_deflection,
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1)
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, fuse_width / 2 + blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=landing_deflection,
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1),
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, fuse_width / 2 + 2 * blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=landing_deflection,
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet3", gain=1, sgn_dup=1),
                            ],
                        ),
                        WingJSec(  # flap1 goes to end of 3rd duct
                            xyz_le=[0, fuse_width / 2 + 3 * blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=landing_deflection,
                                ),
                                asb.ControlSurface(
                                    name="Flap2",
                                    hinge_point=flap_hinge,
                                    deflection=landing_deflection,
                                ),
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet3", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet4", gain=1, sgn_dup=1),
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, aileron_y, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap2",
                                    hinge_point=0.8,
                                    deflection=landing_deflection,
                                ),
                                asb.ControlSurface(
                                    name="Aileron",
                                    symmetric=False,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                                asb.ControlSurface(
                                    name="Flaperon",
                                    symmetric=True,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet4", gain=1, sgn_dup=1)
                            ],
                        ),
                        WingJSec(
                            xyz_le=[
                                ail_hinge * root_chord - ail_hinge * tip_chord,
                                1 / 2 * b,
                                0,
                            ],
                            chord=tip_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Aileron",
                                    symmetric=False,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                                asb.ControlSurface(
                                    name="Flaperon",
                                    symmetric=True,
                                    hinge_point=0.8,
                                    deflection=0,
                                ),
                            ],
                        ),
                    ],
                )

                wingB_climb = JWing(
                    name="Main Wing",
                    symmetric=True,
                    JetParam=JetParam(
                        hdisk=hdisk,
                        fh=0.0,
                        djet0=0.0,
                        djet1=0.0,
                        djet3=0.0,
                        dxdsk=0.07,
                        dndsk=-0.51 * hdisk,
                    ),
                    xsecs=[
                        WingJSec(
                            xyz_le=[0, 0, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=climb_deflection,
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1)
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, fuse_width / 2 + blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=climb_deflection,
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1),
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, fuse_width / 2 + 2 * blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=climb_deflection,
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet3", gain=1, sgn_dup=1),
                            ],
                        ),
                        WingJSec(  # flap1 goes to end of 3rd duct
                            xyz_le=[0, fuse_width / 2 + 3 * blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1",
                                    hinge_point=flap_hinge,
                                    deflection=climb_deflection,
                                ),
                                asb.ControlSurface(
                                    name="Flap2",
                                    hinge_point=flap_hinge,
                                    deflection=climb_deflection,
                                ),
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet3", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet4", gain=1, sgn_dup=1),
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, aileron_y, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap2",
                                    hinge_point=0.8,
                                    deflection=climb_deflection,
                                ),
                                asb.ControlSurface(
                                    name="Aileron",
                                    symmetric=False,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                                asb.ControlSurface(
                                    name="Flaperon",
                                    symmetric=True,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet4", gain=1, sgn_dup=1)
                            ],
                        ),
                        WingJSec(
                            xyz_le=[
                                ail_hinge * root_chord - ail_hinge * tip_chord,
                                1 / 2 * b,
                                0,
                            ],
                            chord=tip_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Aileron",
                                    symmetric=False,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                                asb.ControlSurface(
                                    name="Flaperon",
                                    symmetric=True,
                                    hinge_point=0.8,
                                    deflection=0,
                                ),
                            ],
                        ),
                    ],
                )

                # standard config
                wingB = JWing(
                    name="Main Wing",
                    symmetric=True,
                    JetParam=JetParam(
                        hdisk=hdisk,
                        fh=0.0,
                        djet0=0.0,
                        djet1=0.0,
                        djet3=0.0,
                        dxdsk=0.07,
                        dndsk=-0.51 * hdisk,
                    ),
                    xsecs=[
                        WingJSec(
                            xyz_le=[0, 0, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1", hinge_point=flap_hinge, deflection=0
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1)
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, fuse_width / 2 + blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1", hinge_point=flap_hinge, deflection=0
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1),
                            ],
                        ),
                        WingJSec(
                            xyz_le=[0, fuse_width / 2 + 2 * blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1", hinge_point=flap_hinge, deflection=0
                                )
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet3", gain=1, sgn_dup=1),
                            ],
                        ),
                        WingJSec(  # flap1 goes to end of 3rd duct
                            xyz_le=[0, fuse_width / 2 + 3 * blowing_dy, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap1", hinge_point=flap_hinge, deflection=0
                                ),
                                asb.ControlSurface(
                                    name="Flap2", hinge_point=flap_hinge, deflection=0
                                ),
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet3", gain=1, sgn_dup=1),
                                JetControl(jet_name="FlapJet4", gain=1, sgn_dup=1),
                            ],
                        ),
                        # WingJSec(
                        #     xyz_le=[0, fuse_width/2 + 4*blowing_dy, 0],
                        #     chord=root_chord,
                        #     twist = wing_incidence,
                        #     airfoil=main_foil,
                        #     control_surfaces = [asb.ControlSurface(name="Flap2", hinge_point=flap_hinge, deflection=0)],
                        #     JetControls= [JetControl(jet_name="FlapJet4", gain=1, sgn_dup=1), JetControl(jet_name="FlapJet5", gain=1, sgn_dup=1)],
                        # ),
                        WingJSec(
                            xyz_le=[0, aileron_y, 0],
                            chord=root_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Flap2", hinge_point=0.8, deflection=0
                                ),
                                asb.ControlSurface(
                                    name="Aileron",
                                    symmetric=False,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                                asb.ControlSurface(
                                    name="Flaperon",
                                    symmetric=True,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                            ],
                            JetControls=[
                                JetControl(jet_name="FlapJet4", gain=1, sgn_dup=1)
                            ],
                        ),
                        WingJSec(
                            xyz_le=[
                                ail_hinge * root_chord - ail_hinge * tip_chord,
                                1 / 2 * b,
                                0,
                            ],
                            chord=tip_chord,
                            twist=wing_incidence,
                            airfoil=main_foil,
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Aileron",
                                    symmetric=False,
                                    hinge_point=ail_hinge,
                                    deflection=0,
                                ),
                                asb.ControlSurface(
                                    name="Flaperon",
                                    symmetric=True,
                                    hinge_point=0.8,
                                    deflection=0,
                                ),
                            ],
                        ),
                    ],
                )

                vertical_tail = asb.Wing(
                    name="Vertical Tail",
                    symmetric=False,
                    xsecs=[
                        WingJSec(
                            xyz_le=[-vt_sweep_x, 0, 0],
                            chord=vt_c,
                            twist=0,
                            airfoil=asb.Airfoil(name="NACA0012"),
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Rudder", hinge_point=tail_hinge, deflection=0
                                )
                            ],
                        ),
                        WingJSec(
                            xyz_le=[vt_sweep_x, 0, vt_b],
                            chord=vt_c,
                            twist=0,
                            airfoil=asb.Airfoil(name="NACA0012"),
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Rudder", hinge_point=tail_hinge, deflection=0
                                )
                            ],
                        ),
                    ],
                ).translate([vt_x, 0, 0])  # guess
                horizontal_tail = asb.Wing(
                    name="Horizontal Tail",
                    symmetric=True,
                    xsecs=[
                        WingJSec(
                            xyz_le=[0, 0, 0],
                            chord=ht_rc,
                            twist=0,
                            airfoil=asb.Airfoil(name="NACA0012"),
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Elevator",
                                    hinge_point=tail_hinge,
                                    deflection=0,
                                )
                            ],
                        ),
                        WingJSec(
                            xyz_le=[
                                ht_rc - ht_tc,
                                ht_b,
                                0,
                            ],  # constant TE, okay for now but should update to be constant hinge line
                            chord=ht_tc,
                            twist=0,
                            airfoil=asb.Airfoil(name="NACA0012"),
                            control_surfaces=[
                                asb.ControlSurface(
                                    name="Elevator",
                                    hinge_point=tail_hinge,
                                    deflection=0,
                                )
                            ],
                        ),
                    ],
                ).translate([ht_x, 0, vt_b])
                fuselage_xsecs = generate_fuselage_xsecs(
                    N=10, vt_x=vt_x, vt_sweep_x=vt_sweep_x, vt_c=vt_c
                )
                fuselage = asb.Fuselage(name="Fuselage", xsecs=fuselage_xsecs)
                analysis_options = {
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
                    asb.Fuselage: dict(panel_resolution=24, panel_spacing="cosine"),
                }

                # NOTE: Uncomment fuselage in asb.Airplane(...) below if you are a non-Mac user
                planeB = asb.Airplane(
                    name="Initial Aircraft",
                    xyz_ref=[MAC / 4, 0, 0],
                    wings=[wingB, vertical_tail, horizontal_tail],
                    # fuselages=[fuselage]
                )
                planeB_takeoff = asb.Airplane(
                    name="Initial Aircraft",
                    xyz_ref=[MAC / 4, 0, 0],
                    wings=[wingB_takeoff, vertical_tail, horizontal_tail],
                    # fuselages=[fuselage]
                )
                planeB_climb = asb.Airplane(
                    name="Initial Aircraft",
                    xyz_ref=[MAC / 4, 0, 0],
                    wings=[wingB_climb, vertical_tail, horizontal_tail],
                    # fuselages=[fuselage]
                )
                planeB_landing = asb.Airplane(
                    name="Initial Aircraft",
                    xyz_ref=[MAC / 4, 0, 0],
                    wings=[wingB_landing, vertical_tail, horizontal_tail],
                    # fuselages=[fuselage]
                )
                planeB_unblown = asb.Airplane(
                    name="Initial Aircraft",
                    xyz_ref=[MAC / 4, 0, 0],
                    wings=[wingB_unblown, vertical_tail, horizontal_tail],
                    # fuselages=[fuselage]
                )

                if phase.lower() == "cruise":
                    jvl_plane = JVL(
                        airplane=planeB_unblown,
                        op_point=asb.OperatingPoint(
                            velocity=velocity,
                            alpha=alpha,
                            atmosphere=Atmosphere(altitude=h_cruise),
                            beta=0,
                            p=0,
                            q=0,
                            r=0,
                        ),  # no ground effect
                        avl_command="jvl",
                    )
                elif phase.lower() == "climb":
                    jvl_plane = JVL(
                        airplane=planeB_climb,
                        op_point=asb.OperatingPoint(
                            velocity=velocity,
                            alpha=alpha,
                            atmosphere=Atmosphere(altitude=h_climb),
                            beta=0,
                            p=0,
                            q=0,
                            r=0,
                        ),  # no ground effect
                        avl_command="jvl",
                    )
                else:
                    plane = [planeB_landing if phase == "landing" else planeB_takeoff][
                        0
                    ]
                    jvl_plane = JVL(
                        airplane=plane,
                        op_point=asb.OperatingPoint(
                            velocity=velocity,
                            alpha=alpha,
                            beta=0,
                            p=0,
                            q=0,
                            r=0,
                        ),
                        ground_effect=True,
                        ground_effect_height=-(1 + fuse_width),  # rough estimate
                        avl_command="jvl",
                    )

                jvl_plane.default_analysis_specific_options = analysis_options
                jvl_plane.write_jvl(
                    f"./JVL_writer/sref_trades/run_outputs/{folder_name}/geoms/{config}",
                    # f"JVL_writer/sref_design-trades/run_outputs/geoms/{config}",
                    CLAF=False,
                    j=True,
                )

                results.append(
                    run_case(phase, jvl_plane, velocity, alpha, blowing_perturb)
                )  # run JVL
                print(
                    f"Successfully ran case {case_idx + 1} of {len(cases)} in {phase}"
                )
            save_results(
                results,
                f"./JVL_writer/sref_trades/run_outputs/{folder_name}/coeff_results/Sref-{S}/{phase}.pkl",
            )


if __name__ == "__main__":
    S_list = np.linspace(40, 52, 13)
    # S_list = np.linspace(40, 52, 13)

    # NOTE: Technically stall velocity slightly changes, but so does mass. At this stage, we can't
    # accurately predict the relative magnitude of W/S & which effect dominates due to lack of correct mass models.
    # Assumed that v_stall is largely constant. (W / (0.5 * rho * S))^0.5
    v_stall = 20  # [m/s]; (7600 / (0.5 * 1.225 * S_list)) ** 0.5 --> exact but rounded to 20 for simplicitly

    oper_dict = {
        "takeoff": {
            "alphas": np.linspace(14, 20, 3),
            "velocities": np.array(
                [v_stall * 1.1]
            ),  # largely unaffected by V_inf; dominated by blowing
        },
        "climb": {
            "alphas": np.array([20, 25, 30]),
            "velocities": np.array([20]),
        },
        "cruise": {
            "alphas": np.array([0]),
            "velocities": np.array([80, 125, 150]),
        },
        "landing": {
            "alphas": np.linspace(12, 18, 3),
            "velocities": np.array(
                [v_stall * 1.1]  # TODO: Look into precise value
            ),
        },
    }

    for folder_name, val in {"full_blow": 1.0, "half_blow": 0.5}.items():
        run_sref_cases(
            S_list, oper_dict=oper_dict, folder_name=folder_name, blowing_perturb=val
        )  # NOTE: blowing_perturb ~ [0, 1]
