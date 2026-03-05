import aerosandbox as asb
from aerosandbox.tools import units
from geom import plane

total_mass_properties = asb.MassProperties(mass=0)
num_motors = 12
spacing_between_motors = 5*units.inch
motor_mass_props = asb.MassProperties(mass=0)
fan_mass_properties = asb.MassProperties(mass=0)
motor_mass = 0.02 #from propulsion
wing = plane.wings[0]
ht = plane.wings[2]
vt = plane.wings[1]
S = plane.s_ref
b = plane.b_ref

for i in range(int(num_motors/2)):

    motor_mass_props = motor_mass_props + asb.MassProperties(
        mass=motor_mass,
        x_cg = 0, y_cg=10*units.inch+spacing_between_motors*i, z_cg=0
    ) + asb.MassProperties(mass=motor_mass,
        x_cg = 0, y_cg=-10*units.inch+spacing_between_motors*i, z_cg=0
    )


    fan_mass_properties = fan_mass_properties + asb.MassProperties(
        mass = 0.36/num_motors,
        x_cg = 0, y_cg=10*units.inch+spacing_between_motors*i, z_cg=0
    ) + asb.MassProperties(mass = 0.36/num_motors, 
        x_cg = 0, y_cg=-10*units.inch+spacing_between_motors*i, z_cg=0
    )
total_mass_properties = total_mass_properties + motor_mass_props + fan_mass_properties

fuselage_mass_props = asb.MassProperties(
    mass = 0.2,
    x_cg = wing.xsecs[0].chord*0.5

)

servo_mass_props = asb.MassProperties(
    mass=0.3,
    x_cg = 0.25
)

motor_wiring_props = asb.MassProperties(
    mass = 0.18,
    x_cg=0.1
)

avionics_mass_props = asb.MassProperties(
    mass = 0.18,
    x_cg=0.25)

lipos_mass_props = asb.MassProperties(
    mass = 1.1,
    x_cg=0.25)

boom_mass_props = asb.MassProperties(
    mass = 0.7,
    x_cg = (ht.xsecs[0].xyz_le[0]-wing.xsecs[0].xyz_le[0])
)

total_mass_properties = total_mass_properties + fuselage_mass_props + servo_mass_props + motor_wiring_props + avionics_mass_props + lipos_mass_props + boom_mass_props

foam_density = 48/10 #kg/m^3
fiberglass_density = 2.6*1000 #kg/m^3
fiberglass_youngs_modulus = 70e9 #Pa
fiberglass_shear_modulus = 30e9 #Pa
fiberglass_layup_epoxy_factor = 2.2

carbon_fiber_density = 1.75*1000 #kg/m^3
carbon_fiber_yield_strength = 600e6 #Pa
carbon_fiber_youngs_modulus = 230e9 #Pa
carbon_fiber_shear_modulus = 50e9 #Pa
carbon_fiber_layup_epoxy_factor = 2.2

plywood_density = 550 #kg/m^3

# main wing mass (foam estimate)
ht_surface_area = plane.wings[2].area(type='wetted')
vt_surface_area = plane.wings[1].area(type='wetted')
layup_thickness = 0.00004 # 0.04mm

balsa_density = 160 #kg/m^3
balsa_youngs_modulus = 4e9 #Pa
balsa_shear_modulus = 1.5e9 #Pa
balsa_layup_epoxy_factor = 1.2

# spar dimensions (estimate)
spar_width = 0.005
cap_height = 0.001
spar_height = wing.xsecs[0].airfoil.max_thickness()*wing.xsecs[0].chord
spar_core_volume = spar_height*spar_width*b
spar_surface_area = 2*(spar_height+spar_width)*b
core_mass = spar_core_volume*balsa_density*balsa_layup_epoxy_factor
print(f'core mass: {core_mass}')
spar_fiberglass_mass = spar_surface_area*fiberglass_density*layup_thickness*fiberglass_layup_epoxy_factor
print(f'spar fiberglass mass: {spar_fiberglass_mass}')
spar_cap_mass = 2*(spar_width*cap_height)*b*carbon_fiber_density*carbon_fiber_layup_epoxy_factor
print(f'spar cap mass: {spar_cap_mass}')
spar_mass = core_mass + spar_fiberglass_mass + spar_cap_mass

S_fiberglass_mass = S*fiberglass_density*layup_thickness*fiberglass_layup_epoxy_factor/10
S_foam_mass = wing.volume()*foam_density
foam_wing_mass = S_fiberglass_mass + S_foam_mass

ht_fiberglass_mass = ht_surface_area*fiberglass_density*layup_thickness*fiberglass_layup_epoxy_factor
ht_foam_mass = ht.volume()*foam_density

vt_fiberglass_mass = vt_surface_area*fiberglass_density*layup_thickness*fiberglass_layup_epoxy_factor
vt_foam_mass = vt.volume()*foam_density

foam_wing_mass_props = asb.MassProperties(mass=foam_wing_mass, x_cg=wing.xsecs[0].chord*0.5)
ht_mass_props = asb.MassProperties(mass=ht_fiberglass_mass + ht_foam_mass, x_cg=ht.xsecs[0].chord*0.5+ht.xsecs[0].xyz_le[0])
vt_mass_props = asb.MassProperties(mass=vt_fiberglass_mass + vt_foam_mass, x_cg=vt.xsecs[0].chord*0.5+vt.xsecs[0].xyz_le[0])

total_mass_properties = (total_mass_properties  + 
     servo_mass_props + motor_wiring_props + avionics_mass_props + lipos_mass_props + ht_mass_props + vt_mass_props + foam_wing_mass_props)

total_mass_properties.export_AVL_mass_file('jvl_test.mass')

# airframe_mass_props = airframe_mass_props - foam_wing_mass_props - ht_mass_props - vt_mass_props

# main wing mass (center built up)
# rib_spacing = 0.15
# num_ribs = round(b/rib_spacing)
# rib_thickness = 0.004
# rib_density = plywood_density
# # approximating that  60% of wing is built up, the rest is foam + spar
# rib_area = wing.xsecs[0].chord*wing.xsecs[0].airfoil.area()*0.6 * 1/3 # 1/3 factor from experience, if structure is efficient (ie cut out plywood where not needed)
# rib_mass = rib_area*rib_thickness*rib_density*num_ribs

# part_foam_mass = wing.xsecs[0].chord*wing.xsecs[0].airfoil.area()*0.4*foam_density*b


# print(part_foam_mass)
# print(f'rib mass: {rib_mass}')
# print(f'part foam mass: {part_foam_mass}')
# print(f'spar mass: {spar_mass}')
# print(f'S fiberglass mass: {S_fiberglass_mass}')

# part_built_wing_mass = rib_mass + part_foam_mass + spar_mass + S_fiberglass_mass # monokote mass is comparableto foam, so just use the same area
# part_built_mass_props = asb.MassProperties(mass=part_built_wing_mass, x_cg=wing.xsecs[0].chord*0.5)

# print(f'htht_mass_props: {ht_mass_props}')
# print(f'vt_mass_props: {vt_mass_props}')
# print(f'part_built_mass_props: {part_built_mass_props}')

# print(foam_wing_mass_props)
# print(part_built_mass_props)
# print(ht_mass_props)
# print(vt_mass_props)

# print(total_mass_properties)