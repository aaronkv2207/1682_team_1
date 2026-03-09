import matplotlib.pyplot as plt
from ambiance import Atmosphere
from pint import UnitRegistry

ureg = UnitRegistry()

# mission requirements
g = 9.81 * ureg("m/s^2")
PAYLOAD = 19 * 100 * ureg("kg")
RANGE = 1500 * ureg("miles").to("m")
V_CRUISE = 125 * ureg("m/s")
X_TAKEOFF = 300 * ureg("ft").to("m")  # will need to be larger

# twin otter 300 model: https://www.719skvadron.no/dhc6/dhc6-spec.htm
wing_span_ot = 65 * ureg("feet").to_base_units()
length_ot = 52 * ureg("feet").to_base_units()
height_ot = 19 * ureg("feet").to_base_units()
tailplane_span_ot = 21 * ureg("feet").to_base_units()
wheel_track_ot = 12 * ureg("feet").to_base_units()
wheel_base_ot = 15 * ureg("feet").to_base_units()

wing_area_ot = 420.0 * ureg("ft^2").to_base_units()
aileron_area_ot = 33.2 * ureg("ft^2").to_base_units()
flap_area_ot = 112.2 * ureg("ft^2").to_base_units()
fin_area_ot = 48.0 * ureg("ft^2").to_base_units()
rudder_area_incl_tab_ot = 34.0 * ureg("ft^2").to_base_units()
tailplane_area_ot = 100.0 * ureg("ft^2").to_base_units()
elevator_area_incl_tabs_ot = 35.0 * ureg("ft^2").to_base_units()
n_props_ot = 2
p_shaft_ot = 620 * ureg("hp").to("W")
w_s_ot = 145.50 * ureg("kg / m^2")
w_p_ot = 10.1 * ureg("lbs/hp").to("kg/W")  # power loading

weight_pgen = (
    50 * ureg("lbs").to_base_units()
) * 2  # https://www.startergenerator.com/inventory/250SG114Q

mtow_ot = 12500 * ureg("lbs").to_base_units()
# mlw_ot = 12300 * ureg("lbs").to_base_units()
weight_pay_fuel = (
    4400 * ureg("lb").to_base_units() * 2.0
)  # https://www.nohrsc.noaa.gov/snowsurvey/twin_otter.html; added fuel weight for improved range


# NOTE: preliminary assumption; will inform mass fractions for fuel and structures
weight_fixed = (
    PAYLOAD * 1.10
)  # assumed slight mass increase for other equipment (avionics, e.g.)
weight_pow = weight_pgen + (weight_pay_fuel - weight_fixed)
weight_struct = mtow_ot * 0.6

# xcg = ...  # TODO
MTOW = weight_pow + weight_fixed + weight_struct
pow_f = (weight_pow / MTOW).magnitude
struct_f = (weight_struct / MTOW).magnitude
fix_f = (weight_fixed / MTOW).magnitude

#############
v_stall_ot = 28.8 * ureg("m/s") * (MTOW / mtow_ot)
v_approach_ot = 1.3 * v_stall_ot
v_cruise_ot = 182 * ureg("knots").to_base_units()

x_to = 366 * ureg("m")  # src: https://dehavilland.com/twin-otter-classic-300-g/
x_land = 320 * ureg("m")

# modified parameters to size our plane
cl_max = 6.1  # TODO: Fill in based on wind tunnel data
atm = Atmosphere(h=0)
rho = atm.density[0] * ureg("kg/m^3")

# NOTE: x_to proportional to v^3 (1st order approx.)
beta = (1 / 7) ** (1 / 3)  # will enable for reduction in takeoff length
V_STALL = v_stall_ot * beta

W_S = ((0.5 * rho * V_STALL**2) * cl_max).to("N/m^2")
S = (MTOW * g).to("N") / W_S

print_bool = True
if print_bool:
    print(
        f"{35 * '='}\nMass Fractions from correlation\n"
        f"MTOW: {round(MTOW.to('lbs'), 2)}\n"
        f"pow_f: {round(pow_f, 2)}\n"
        f"struct_f: {round(struct_f, 2)}\n"
        f"fix_f: {round(fix_f, 2)}\n{35 * '-'}\n"
        f"V_STALL: {round(V_STALL, 2)}\n"
        f"W/S: {round(W_S, 2)}\n"
        f"S: {round(S, 2)}\n{35 * '='}"
    )
