import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox.tools import units
from typing import List
from J import JetParam, JetControl, WingJSec, JVL, JWing


jw01 = asb.Airfoil(coordinates='.\\jw01.dat')
wing = JWing(
    name="Main Wing",
    symmetric=True,
    JetParam=JetParam(hdisk=0.188, fh=0.0, djet0=0.0, djet1=0.0, djet3=0.0),
    xsecs=[
        # WingJSec(
        #     xyz_le=[0, 0, 0],
        #     chord=15*units.inch,
        #     twist=0,
        #     airfoil=jw01
        # ),
        WingJSec(
            xyz_le=[0, 0, 0],
            chord=15*units.inch,
            twist=0,
            airfoil=jw01,
            control_surfaces = [asb.ControlSurface(name="Flap1", hinge_point=0.66, deflection=0)],
            JetControls= [JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1)],
        ),
        WingJSec(
            xyz_le=[0, 23*units.inch, 0],
            chord=15*units.inch,
            twist=0,
            airfoil=jw01,
            control_surfaces = [asb.ControlSurface(name="Flap1", hinge_point=0.66, deflection=0), asb.ControlSurface(name="Flap2", hinge_point=0.66, deflection=0)],
            JetControls= [JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1), JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1)],
        ),
        WingJSec(
            xyz_le=[0, 41*units.inch, 0],
            chord=15*units.inch,
            twist=0,
            airfoil=jw01,
            control_surfaces = [asb.ControlSurface(name="Flap2", hinge_point=0.66, deflection=0), asb.ControlSurface(name="Aileron", symmetric=False, hinge_point=0.66, deflection=0)],
            JetControls= [JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1), JetControl(jet_name="AilJet", gain=1, sgn_dup=-1)],
        ),
        WingJSec(
            xyz_le=[6*units.inch, 60*units.inch, 0],
            chord=10*units.inch,
            twist=0,
            airfoil=jw01,
            control_surfaces = [asb.ControlSurface(name="Aileron", symmetric=False, hinge_point=0.66, deflection=0)],
            JetControls= [JetControl(jet_name="AilJet", gain=1, sgn_dup=-1)],
        )
    ]
)
vertical_tail = JWing(
    name="Vertical Tail",
    symmetric=False,
    xsecs=[
        WingJSec(
            xyz_le=[0, 0, 0],
            chord=15*units.inch,
            twist=0,
            airfoil=asb.Airfoil(name="NACA0012"),
            control_surfaces = [asb.ControlSurface(name="Rudder", hinge_point=0.66, deflection=0)]
        ),
        WingJSec(
            xyz_le=[10.5*units.inch, 0, 15*units.inch],
            chord=11*units.inch,
            twist=0,
            airfoil=asb.Airfoil(name="NACA0012"),
            control_surfaces = [asb.ControlSurface(name="Rudder", hinge_point=0.66, deflection=0)]
        )
    ]
).translate([52.5*units.inch, 0, 0])
    
horizinatal_tail = JWing(
    name="Horizontal Tail",
    symmetric=True,
    xsecs=[
        WingJSec(
            xyz_le=[0, 0, 0],
            chord=12*units.inch,
            twist=0,
            airfoil=asb.Airfoil(name="NACA0012"),
            control_surfaces = [asb.ControlSurface(name="Elevator", hinge_point=0.5, deflection=0)]
        ),
        WingJSec(
            xyz_le=[3*units.inch, 25*units.inch, 0],
            chord=7*units.inch,
            twist=0,
            airfoil=asb.Airfoil(name="NACA0012"),
            control_surfaces = [asb.ControlSurface(name="Elevator", hinge_point=0.5, deflection=0)]
        )
    ]
).translate([65*units.inch, 0, 15*units.inch])

def generate_fuselage_xsecs(N: int) -> List[asb.FuselageXSec]:
    """
    Generates a fuselage with N sections, transitioning from a small circular nose,
    to a large square-like midsection, and tapering into a smaller square tail.

    Args:
        N (int): Number of sections defining the fuselage.

    Returns:
        List[FuselageXSec]: List of fuselage cross-sections.
    """
    x_positions = np.linspace(-43, 60, N)  # Generate N sections along the fuselage
    zs = np.interp(x_positions, [-43, -20, 25, 52], [-12, -12, -5, -5])  # Interpolate height
    widths = np.interp(x_positions, [-43, -20, 25, 52, 60], [3.0, 10.0, 10.0, 3.0, 1.0])  # Width transition
    heights = np.interp(x_positions, [-43, -20, 25, 52, 60], [5.0, 15, 15.0, 5.0, 2.0])  # Height transition
    shapes = np.interp(x_positions, [-43, -20, 25, 52], [2, 50, 50, 50])  # Shape transition

    xsecs = []
    for x, z, width, height, shape in zip(x_positions, zs, widths, heights, shapes):
        xsecs.append(asb.FuselageXSec(
            xyz_c=[x*units.inch, 0, z*units.inch],
            xyz_normal=[1, 0, 0],  # Assume fuselage is aligned along x-axis
            width=width*units.inch,
            height=height*units.inch,  # Rectangular cross-section
            shape=shape
        ))

    return xsecs

fuselage_xsecs = generate_fuselage_xsecs(10)  # Generates 10 sections
fuselage = asb.Fuselage(
    name='Fuselage',
    xsecs=fuselage_xsecs
)

plane = asb.Airplane(
    name="Initial Aircraft",
    xyz_ref=[0, 0, 0],
    wings=[wing, vertical_tail, horizinatal_tail],
    # fuselages=[fuselage]
)

avl_plane = JVL(
    airplane=plane,
    op_point=asb.OperatingPoint(
        velocity=100,
        alpha=5,
        beta=0,
        p=0,
        q=0,
        r=0,
    ),
    avl_command='.\\jvl2.20.exe')
avl_plane.default_analysis_specific_options = {
        asb.Airplane: dict(profile_drag_coefficient=0),
        JWing: dict(
            wing_level_spanwise_spacing=True,
            spanwise_resolution=25,
            spanwise_spacing="cosine",
            chordwise_resolution=25,
            chordwise_spacing="cosine",
            component=None,  # This is an int
            no_wake=False,
            no_alpha_beta=False,
            no_load=False,
            drag_polar=None,
        ),
        asb.Wing: dict(
            wing_level_spanwise_spacing=True,
            spanwise_resolution=12,
            spanwise_spacing="cosine",
            chordwise_resolution=12,
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

avl_plane.write_jvl('821p1', CLAF=False, j=True)

# """ Defining Plane Based on Sections"""

# Wing

## Section 1:

# jw01 = asb.Airfoil(name='.\\jw01.dat')
# jw01.draw()