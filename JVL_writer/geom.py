import aerosandbox as asb
import aerosandbox.numpy as np
from aerosandbox.tools import units
import aerosandbox.tools.pretty_plots as p
from typing import List
from J import JetParam, JetControl, WingJSec, JVL, JWing
from scipy.interpolate import Akima1DInterpolator

# Wing geometry
main_foil = asb.Airfoil(coordinates='./JVL_writer/jw05.dat')
S = 49.6 #twin otter wing area is 39 m^2
AR = 8 #twin otter is 10.05
ail_hinge = 0.7

b = np.sqrt(AR*S)
aileron_y = 2/3 * b/2 #start aileron at 2/3 of halfspan. No blowing at this point.
MAC = b/AR
tip_chord = 4/4*MAC #set taper ratio if desired
root_chord = 2*MAC - tip_chord #root chord needs to be adjusted for lost area from taper at the wingtips

# vertical tail geometry
lv = 10 #distance from wing quarter chord to vertical tail quarter chord
#TODO: update VT moment arm to account for the sweep
Vv = .1 #Vertical tail volume coefficient 
vt_ar = 1.2
tail_hinge = 0.7

Sv = Vv*S*b/lv
vt_c = np.sqrt(Sv/vt_ar) #constant chord
vt_b = vt_ar*vt_c
vt_x = lv + 1/4 * MAC - 1/4 * vt_c #distance from wing LE to tail (avg) LE

# horizontal tail geometry
lh = lv + vt_c/3 #include some sweep in the vertical tail to extend the ht moment arm
Vh = 1.05
ht_ar = 2

Sh = Vh*S*MAC/lh
ht_mac = np.sqrt(Sh/ht_ar) #technically not mac but geometry mean chord, but whatever
ht_rc = vt_c
ht_tc = 2*ht_mac - ht_rc
ht_b = ht_ar*ht_mac

ht_x = lh + 1/4 * MAC - 1/4 * ht_mac
vt_sweep_x = ht_x - vt_x

# fuselage (not very well parametrized at the moment)
nose_x = -5
fuse_width = 2.5 #width at the wing

wing = JWing(
    name="Main Wing",
    symmetric=True,
    JetParam=JetParam(hdisk=0.188, fh=0.0, djet0=0.0, djet1=0.0, djet3=0.0),
    xsecs=[
        WingJSec(
            xyz_le=[0, 0, 0],
            chord=root_chord,
            twist=0,
            airfoil=main_foil,
        ),
        WingJSec(
            xyz_le=[0, 1/2*fuse_width, 0],
            chord=root_chord,
            twist=0,
            airfoil=main_foil,
            control_surfaces = [asb.ControlSurface(name="Flap1", hinge_point=0.8, deflection=0)],
            JetControls= [JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1)],
        ),
        WingJSec(
            xyz_le=[0, (aileron_y + 1/2*fuse_width)/2, 0],
            chord=root_chord,
            twist=0,
            airfoil=main_foil,
            control_surfaces = [asb.ControlSurface(name="Flap1", hinge_point=0.8, deflection=0), asb.ControlSurface(name="Flap2", hinge_point=0.8, deflection=0)],
            JetControls= [JetControl(jet_name="FlapJet1", gain=1, sgn_dup=1), JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1)],
        ),
        WingJSec(
            xyz_le=[0, aileron_y, 0],
            chord=root_chord,
            twist=0,
            airfoil=main_foil,
            control_surfaces = [asb.ControlSurface(name="Flap2", hinge_point=0.8, deflection=0), asb.ControlSurface(name="Aileron", symmetric=False, hinge_point=0.8, deflection=0)],
            JetControls= [JetControl(jet_name="FlapJet2", gain=1, sgn_dup=1), JetControl(jet_name="AilJet", gain=1, sgn_dup=-1)],
        ),
        WingJSec(
            xyz_le=[ail_hinge*root_chord-ail_hinge*tip_chord, 1/2*b, 0],
            chord=tip_chord,
            twist=0,
            airfoil=main_foil,
            control_surfaces = [asb.ControlSurface(name="Aileron", symmetric=False, hinge_point=ail_hinge, deflection=0)],
            JetControls= [JetControl(jet_name="AilJet", gain=1, sgn_dup=-1)],
        )
    ]
)
vertical_tail = JWing(
    name="Vertical Tail",
    symmetric=False,
    xsecs=[
        WingJSec(
            xyz_le=[-vt_sweep_x, 0, 0],
            chord=vt_c,
            twist=0,
            airfoil=asb.Airfoil(name="NACA0012"),
            control_surfaces = [asb.ControlSurface(name="Rudder", hinge_point=tail_hinge, deflection=0)]
        ),
        WingJSec(
            xyz_le=[vt_sweep_x, 0, vt_b],
            chord=vt_c,
            twist=0,
            airfoil=asb.Airfoil(name="NACA0012"),
            control_surfaces = [asb.ControlSurface(name="Rudder", hinge_point=tail_hinge, deflection=0)]
        )
    ]
).translate([vt_x, 0, 0]) #guess
    
horizinatal_tail = JWing(
    name="Horizontal Tail",
    symmetric=True,
    xsecs=[
        WingJSec(
            xyz_le=[0, 0, 0],
            chord=ht_rc,
            twist=0,
            airfoil=asb.Airfoil(name="NACA0012"),
            control_surfaces = [asb.ControlSurface(name="Elevator", hinge_point=tail_hinge, deflection=0)]
        ),
        WingJSec(
            xyz_le=[ht_rc - ht_tc, ht_b, 0], #constant TE, okay for now but should update to be constant hinge line
            chord=ht_tc,
            twist=0,
            airfoil=asb.Airfoil(name="NACA0012"),
            control_surfaces = [asb.ControlSurface(name="Elevator", hinge_point=tail_hinge, deflection=0)]
        )
    ]
).translate([ht_x, 0, vt_b])

def generate_fuselage_xsecs(N: int) -> List[asb.FuselageXSec]:
    """
    Generates a fuselage with N sections, transitioning from a small circular nose,
    to a large square-like midsection, and tapering into a smaller square tail.

    Args:
        N (int): Number of sections defining the fuselage.

    Returns:
        List[FuselageXSec]: List of fuselage cross-sections.
    """
    fuse_end = vt_x - vt_sweep_x + vt_c #end of the fuselage, at the end of the vertical tail
    x_positions = np.linspace(nose_x, fuse_end, N)  # Generate N sections along the fuselage
    r = np.array([0.3, 2/3*fuse_width, fuse_width, fuse_width, fuse_width, 3/4*fuse_width, 1.0])
    z = -r/2
    zs = np.interp(x_positions, [nose_x, 7.8*nose_x, nose_x/2, 0, vt_x/2, vt_x, fuse_end], z)  # Interpolate z position
    widths = np.interp(x_positions, [nose_x, 7/8*nose_x, 1/2*nose_x, 0, vt_x/2, vt_x, fuse_end], r)  # Width transition
    heights = Akima1DInterpolator([nose_x, 7/8*nose_x, 1/2*nose_x, 0, vt_x/2, vt_x, fuse_end], r, method="akima")(x_positions)
    # heights = np.interp(x_positions, [nose_x, 7/8*nose_x, nose_x/2, vt_x/2, vt_x, vt_x + vt_c], [0, 2, 3, 3, 2.5, 2.0])  # Height transition
    shapes = np.interp(x_positions, [nose_x, nose_x/2, vt_x/2, vt_x, vt_x + vt_c], [1, 2.5, 2.5, 2, 1])  # Shape transition

    xsecs = []
    for x, z, width, height, shape in zip(x_positions, zs, widths, heights, shapes):
        xsecs.append(asb.FuselageXSec(
            xyz_c=[x, 0, z],
            xyz_normal=[1, 0, 0],  # Assume fuselage is aligned along x-axis
            width=width,
            height=height,  # Rectangular cross-section
            shape=shape
        ))

    return xsecs

fuselage_xsecs = generate_fuselage_xsecs(15)
fuselage = asb.Fuselage(
    name='Fuselage',
    xsecs=fuselage_xsecs
)

plane = asb.Airplane(
    name="Initial Aircraft",
    xyz_ref=[0, 0, 0],
    wings=[wing, vertical_tail, horizinatal_tail],
    fuselages=[fuselage]
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

avl_plane.write_jvl('./JVL_files/1682_v0', CLAF=False, j=True)

plane.draw_three_view(show=False)
p.show_plot(tight_layout=False, savefig="3view.png")

print(f'Wing span: {b:.2f} m')
print(f'Wing aspect ratio: {AR:.2f}')
print(f'Wing mean aerodynamic chord: {MAC:.2f} m')