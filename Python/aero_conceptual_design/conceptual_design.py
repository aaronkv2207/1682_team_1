import numpy as np
import matplotlib.pyplot as plt
from pint import UnitRegistry
from ambiance import Atmosphere

ureg = UnitRegistry()

# mission requirements
g = 9.81 * ureg("m/s^2")
N = 19
PAY = 100 * ureg("kg")
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

w_s_ot = 145.50 * ureg("kg / m^2")
w_p_ot = 10.1 * ureg("lbs/shp").to_base_units()  # power loading
mtow_ot = 5670 * ureg("kg")
v_stall_ot = 28.8 * ureg("m/s")
v_approach_ot = 1.3 * v_stall_ot
v_cruise_ot = 182 * ureg("knots").to_base_units()

n_props_ot = 2
p_shaft_ot = 620 * ureg("shp").to_base_units()

# modified parameters to size our plane
cl_max = ...  # TODO: Fill in based on wind tunnel data
atm = Atmosphere(h=0)
rho = atm.density[0]

w_s = (0.5 * rho * v_stall_ot**2) * cl_max
mtow = mtow_ot  # NOTE: preliminary assumption; will inform mass fractions for fuel and structures
xcg = ...
fuel_f = ...
struct_f = ...
