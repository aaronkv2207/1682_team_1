import numpy as np
import math

'''
Inputs: Fan mass (biruk), Cd0 and Cdi_cruise (osahon), CLmax (osahon), AR (osahon), wing area (osahon), v_cruise (osahon), takeoff power (biruk)

'''

def get_range(
    mass,
    Cd0,
    Cdi_cruise,
    CLmax,
    AR,
    wing_area,
    v_cruise,
    takeoff_power,
    climb_angle,
    cruise_altitude,
    mass_budget=2700,
    motor_mass=35,
    fan_mass=10,
    duct_mass=35,
    num_fans=8,
):
    """
    Unified mission-range sizing function.

    Inputs
    ------
    mass : float
        Aircraft mass [kg]
    Cd0 : float
        Zero-lift drag coefficient [-]
    Cdi_cruise : float
        Induced drag coefficient at the specified cruise condition [-]
    CLmax : float
        Maximum lift coefficient used for stall / climb-speed floor [-]
    AR : float
        Aspect ratio [-] (used only to report implied Oswald efficiency)
    wing_area : float
        Wing reference area [m^2]
    v_cruise : float
        Cruise true airspeed [m/s]
    takeoff_power : float
        Required takeoff electrical power, used directly as generator nameplate [kW]
    climb_angle : float
        Steady climb flight-path angle [deg]
    cruise_altitude : float
        Cruise altitude [m]
    mass_budget : float, optional
        Total propulsion-system mass budget [kg]
    motor_mass : float, optional
        Per-fan motor mass [kg]
    fan_mass : float, optional
        Per-fan fan mass [kg]
    duct_mass : float, optional
        Per-fan duct mass [kg]
    num_fans : int, optional
        Number of fans [-]

    Returns
    -------
    results : dict
        Dictionary containing range and all major intermediate quantities.
    """

    # -------------------------
    # Fixed assumptions
    # -------------------------
    g = 9.80665
    eta_prop = 0.80
    eta_motor = 0.95
    eta_inverter = 0.98
    eta_gearbox = 1.00
    eta_chain = eta_prop * eta_motor * eta_inverter * eta_gearbox

    thermal_efficiency = 0.30          # fuel -> electrical
    LHV_MJ_per_kg = 42.0               # jet fuel lower heating value
    generator_specific_power = 6.0     # kW/kg

    # -------------------------
    # Helpers
    # -------------------------
    def isa_density(alt_m):
        """
        ISA density up to 11 km. If altitude exceeds 11 km,
        it is clipped to 11 km to avoid breaking.
        """
        R = 287.05287
        T0 = 288.15
        p0 = 101325.0
        L = 0.0065

        alt_m = max(0.0, float(alt_m))
        alt_m = min(alt_m, 11000.0)

        T = T0 - L * alt_m
        p = p0 * (T / T0) ** (g / (R * L))
        return p / (R * T)

    def fuel_mass_from_electrical_energy_kwh(E_elec_kwh):
        """
        Converts electrical energy delivered by the turbogenerator [kWh]
        into fuel mass [kg] using thermal efficiency and LHV.
        """
        E_MJ = E_elec_kwh * 3.6
        return E_MJ / max(thermal_efficiency * LHV_MJ_per_kg, 1e-12)

    # -------------------------
    # Basic quantities
    # -------------------------
    W = mass * g
    gamma = math.radians(climb_angle)

    # -------------------------
    # Cruise condition
    # -------------------------
    rho_cruise = isa_density(cruise_altitude)
    q_cruise = 0.5 * rho_cruise * v_cruise**2

    CL_cruise = W / max(q_cruise * wing_area, 1e-12)

    # Calibrate parabolic drag model from the given cruise induced drag:
    # Cdi = k * CL^2  =>  k = Cdi_cruise / CL_cruise^2
    k_induced = Cdi_cruise / max(CL_cruise**2, 1e-12)

    CD_cruise = Cd0 + k_induced * CL_cruise**2
    D_cruise = q_cruise * wing_area * CD_cruise

    # Electrical cruise power required at generator output
    P_thrust_cruise_W = D_cruise * v_cruise
    P_elec_cruise_kW = (P_thrust_cruise_W / max(eta_chain, 1e-12)) / 1000.0

    # Cruise fuel flow
    fuel_flow_cruise_kg_per_hr = fuel_mass_from_electrical_energy_kwh(P_elec_cruise_kW)

    # -------------------------
    # Generator sizing from TAKEOFF power
    # -------------------------
    generator_mass_kg = takeoff_power / max(generator_specific_power, 1e-12)

    # -------------------------
    # Fixed propulsion hardware mass
    # -------------------------
    dry_propulsion_mass_kg = (
        generator_mass_kg
        + num_fans * (motor_mass + fan_mass + duct_mass)
    )

    total_fuel_available_kg = mass_budget - dry_propulsion_mass_kg

    if total_fuel_available_kg <= 0.0:
        return {
            "feasible": False,
            "reason": "No mass budget left for fuel after propulsion hardware.",
            "generator_mass_kg": generator_mass_kg,
            "dry_propulsion_mass_kg": dry_propulsion_mass_kg,
            "total_fuel_available_kg": total_fuel_available_kg,
        }

    # -------------------------
    # Climb model
    # Use midpoint density as a simple approximation for the full climb.
    # -------------------------
    rho_climb = isa_density(0.5 * cruise_altitude)

    # Stall-based lower bound for climb-speed sweep
    V_stall_climb = math.sqrt(
        2.0 * W * math.cos(gamma) / max(rho_climb * wing_area * CLmax, 1e-12)
    )

    V_climb_min = 1.2 * V_stall_climb
    V_climb_max = max(1.6 * v_cruise, 2.5 * V_stall_climb)
    V_climb_candidates = np.linspace(V_climb_min, V_climb_max, 400)

    best = None

    for V_climb in V_climb_candidates:
        q_climb = 0.5 * rho_climb * V_climb**2
        CL_climb = (W * math.cos(gamma)) / max(q_climb * wing_area, 1e-12)
        CD_climb = Cd0 + k_induced * CL_climb**2
        D_climb = q_climb * wing_area * CD_climb

        T_req_climb = D_climb + W * math.sin(gamma)

        P_thrust_climb_W = T_req_climb * V_climb
        P_elec_climb_kW = (P_thrust_climb_W / max(eta_chain, 1e-12)) / 1000.0

        # Generator is limited by takeoff-sized nameplate power
        if P_elec_climb_kW > takeoff_power:
            continue

        ROC = V_climb * math.sin(gamma)
        climb_time_hr = (cruise_altitude / max(ROC, 1e-12)) / 3600.0
        climb_energy_kWh = P_elec_climb_kW * climb_time_hr
        climb_fuel_kg = fuel_mass_from_electrical_energy_kwh(climb_energy_kWh)

        if (best is None) or (climb_energy_kWh < best["climb_energy_kWh"]):
            best = {
                "V_climb_opt_mps": V_climb,
                "CL_climb_opt": CL_climb,
                "CD_climb_opt": CD_climb,
                "thrust_required_climb_N": T_req_climb,
                "P_elec_climb_opt_kW": P_elec_climb_kW,
                "ROC_opt_mps": ROC,
                "climb_time_hr": climb_time_hr,
                "climb_energy_kWh": climb_energy_kWh,
                "fuel_mass_climb_kg": climb_fuel_kg,
            }

    if best is None:
        return {
            "feasible": False,
            "reason": "No feasible climb speed found within the takeoff-power generator limit.",
            "generator_mass_kg": generator_mass_kg,
            "dry_propulsion_mass_kg": dry_propulsion_mass_kg,
            "total_fuel_available_kg": total_fuel_available_kg,
            "P_elec_cruise_kW": P_elec_cruise_kW,
            "fuel_flow_cruise_kg_per_hr": fuel_flow_cruise_kg_per_hr,
            "V_stall_climb_mps": V_stall_climb,
        }

    # -------------------------
    # Fuel bookkeeping
    # -------------------------
    fuel_mass_cruise_kg = total_fuel_available_kg - best["fuel_mass_climb_kg"]

    if fuel_mass_cruise_kg < 0.0:
        return {
            "feasible": False,
            "reason": "All available fuel is consumed during climb.",
            "generator_mass_kg": generator_mass_kg,
            "dry_propulsion_mass_kg": dry_propulsion_mass_kg,
            "total_fuel_available_kg": total_fuel_available_kg,
            **best,
        }

    # -------------------------
    # Cruise range
    # -------------------------
    cruise_time_hr = fuel_mass_cruise_kg / max(fuel_flow_cruise_kg_per_hr, 1e-12)
    cruise_range_km = cruise_time_hr * 3600.0 * v_cruise / 1000.0

    # Optional: climb horizontal distance, if you want total mission distance too
    climb_distance_km = cruise_altitude / max(math.tan(gamma), 1e-12) / 1000.0
    total_range_km = climb_distance_km + cruise_range_km

    # Implied Oswald factor from the calibrated k, for reference only
    e_implied = np.nan
    if AR > 0.0 and k_induced > 0.0:
        e_implied = 1.0 / (math.pi * AR * k_induced)

    return {
        "feasible": True,

        # Polar / aero
        "rho_cruise_kgm3": rho_cruise,
        "rho_climb_kgm3": rho_climb,
        "CL_cruise": CL_cruise,
        "CD_cruise": CD_cruise,
        "k_induced": k_induced,
        "e_implied": e_implied,
        "thrust_required_cruise_N": D_cruise,

        # Cruise power / fuel
        "P_elec_cruise_kW": P_elec_cruise_kW,
        "fuel_flow_cruise_kg_per_hr": fuel_flow_cruise_kg_per_hr,

        # Generator and mass bookkeeping
        "generator_mass_kg": generator_mass_kg,
        "dry_propulsion_mass_kg": dry_propulsion_mass_kg,
        "total_fuel_available_kg": total_fuel_available_kg,

        # Climb optimum
        "V_stall_climb_mps": V_stall_climb,
        **best,

        # Remaining cruise fuel
        "fuel_mass_cruise_kg": fuel_mass_cruise_kg,
        "cruise_time_hr": cruise_time_hr,

        # Range
        "cruise_range_km": cruise_range_km,
        "climb_distance_km": climb_distance_km,
        "total_range_km": total_range_km,
    }

if __name__ == "__main__":
    result = get_range(
    mass=16550 / 2.2,
    Cd0=0.019,
    Cdi_cruise=0.003,
    CLmax=6.1,
    AR=8.0,
    wing_area=(16550 / 2.2) / (1484.56 / 9.81),
    v_cruise=125.0,
    takeoff_power=2500.0,      # kW
    climb_angle=10.0,          # deg
    cruise_altitude=18000.0 * 0.3048,
    mass_budget=2700,
    motor_mass=35,
    fan_mass=10,
    duct_mass=35,
    num_fans=8,
)

    print(result["total_range_km"])
    print(result["V_climb_opt_mps"])
    print(result["fuel_mass_climb_kg"])
    print(result["P_elec_cruise_kW"])