import numpy as np
import math


def get_range(
    mass,
    Cd0,
    Cdi_cruise,  # kept for compatibility, not used directly
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
    e_oswald=0.80,
    n_climb_steps=300,
    descent_reserve_fraction_of_climb_fuel=0.50,
):
    """
    Turbogenerator-only mission range sizing function.

    Architecture
    ------------
    - The turbogenerator is sized by sea-level nameplate power, takeoff_power [kW].
    - No battery is included.
    - Generator available power lapses with altitude using:

          P_available(h) = P_SL * delta(h) / sqrt(theta(h))

      where:
          delta = p / p0
          theta = T / T0

    - During climb, the code integrates through altitude and checks whether
      required electrical power exceeds available generator power at each point.
    - Climb speed is selected by minimizing total climb energy among feasible speeds.
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
    def isa_state(alt_m):
        """
        ISA state up to 11 km.

        Returns:
            T      [K]
            p      [Pa]
            rho    [kg/m^3]
            theta  [-]
            delta  [-]
            lapse  [-] = delta / sqrt(theta)
        """

        R = 287.05287
        T0 = 288.15
        p0 = 101325.0
        L = 0.0065

        alt_m = max(0.0, float(alt_m))
        alt_m = min(alt_m, 11000.0)

        T = T0 - L * alt_m
        p = p0 * (T / T0) ** (g / (R * L))
        rho = p / (R * T)

        theta = T / T0
        delta = p / p0
        lapse = delta / max(math.sqrt(theta), 1e-12)

        return {
            "T": T,
            "p": p,
            "rho": rho,
            "theta": theta,
            "delta": delta,
            "lapse": lapse,
        }

    def fuel_mass_from_electrical_energy_kwh(E_elec_kwh):
        """
        Converts electrical energy delivered by the turbogenerator [kWh]
        into fuel mass [kg].
        """
        E_MJ = E_elec_kwh * 3.6
        return E_MJ / max(thermal_efficiency * LHV_MJ_per_kg, 1e-12)

    def aero_at_condition(rho, V, W_normal):
        """
        Computes CL, CDi, CD, and drag for a given density, speed, and normal force.
        """

        q = 0.5 * rho * V**2
        CL = W_normal / max(q * wing_area, 1e-12)
        CDi = k_induced * CL**2
        CD = Cd0 + CDi
        D = q * wing_area * CD

        return CL, CDi, CD, D

    # -------------------------
    # Basic quantities
    # -------------------------
    W = mass * g
    gamma = math.radians(climb_angle)

    # Fixed parabolic drag polar
    k_induced = 1.0 / max(math.pi * e_oswald * AR, 1e-12)

    # -------------------------
    # Generator sizing
    # -------------------------
    # Turbogenerator-only architecture:
    # takeoff_power is the sea-level nameplate rating.
    P_generator_SL_rating_kW = takeoff_power
    generator_mass_kg = P_generator_SL_rating_kW / max(generator_specific_power, 1e-12)

    fixed_propulsion_hardware_mass_kg = (
        generator_mass_kg
        + num_fans * (motor_mass + fan_mass + duct_mass)
    )

    total_fuel_available_kg = mass_budget - fixed_propulsion_hardware_mass_kg

    if total_fuel_available_kg <= 0.0:
        return {
            "feasible": False,
            "reason": "No mass budget left for fuel after generator and propulsion hardware.",
            "generator_mass_kg": generator_mass_kg,
            "fixed_propulsion_hardware_mass_kg": fixed_propulsion_hardware_mass_kg,
            "total_fuel_available_kg": total_fuel_available_kg,
        }

    # -------------------------
    # Cruise condition
    # -------------------------
    cruise_state = isa_state(cruise_altitude)
    rho_cruise = cruise_state["rho"]
    theta_cruise = cruise_state["theta"]
    delta_cruise = cruise_state["delta"]
    lapse_cruise = cruise_state["lapse"]

    CL_cruise, CDi_cruise, CD_cruise, D_cruise = aero_at_condition(
        rho=rho_cruise,
        V=v_cruise,
        W_normal=W,
    )

    P_thrust_cruise_W = D_cruise * v_cruise
    P_elec_cruise_kW = (P_thrust_cruise_W / max(eta_chain, 1e-12)) / 1000.0

    P_generator_cruise_available_kW = P_generator_SL_rating_kW * lapse_cruise
    cruise_power_margin_kW = P_generator_cruise_available_kW - P_elec_cruise_kW

    if P_elec_cruise_kW > P_generator_cruise_available_kW:
        return {
            "feasible": False,
            "reason": "Cruise power required exceeds lapsed generator power available at cruise altitude.",
            "P_elec_cruise_kW": P_elec_cruise_kW,
            "P_generator_cruise_available_kW": P_generator_cruise_available_kW,
            "cruise_power_margin_kW": cruise_power_margin_kW,
            "lapse_cruise": lapse_cruise,
            "generator_mass_kg": generator_mass_kg,
            "fixed_propulsion_hardware_mass_kg": fixed_propulsion_hardware_mass_kg,
            "total_fuel_available_kg": total_fuel_available_kg,
        }

    fuel_flow_cruise_kg_per_hr = fuel_mass_from_electrical_energy_kwh(P_elec_cruise_kW)

    # Minimum-drag cruise speed
    V_min_drag_mps = math.sqrt(
        (2.0 * W / max(rho_cruise * wing_area, 1e-12))
        * math.sqrt(k_induced / max(Cd0, 1e-12))
    )

    # -------------------------
    # Climb speed sweep with altitude-varying lapse
    # -------------------------
    # Use sea-level density for a conservative stall lower bound.
    sea_level_state = isa_state(0.0)
    rho_sl = sea_level_state["rho"]

    V_stall_climb_SL = math.sqrt(
        2.0 * W * math.cos(gamma) / max(rho_sl * wing_area * CLmax, 1e-12)
    )

    V_climb_min = 1.2 * V_stall_climb_SL
    V_climb_max = max(1.6 * v_cruise, 2.5 * V_stall_climb_SL)
    V_climb_candidates = np.linspace(V_climb_min, V_climb_max, 500)

    altitude_grid = np.linspace(0.0, cruise_altitude, n_climb_steps + 1)

    best = None

    for V_climb in V_climb_candidates:
        ROC = V_climb * math.sin(gamma)

        if ROC <= 0.0:
            continue

        total_climb_energy_kWh = 0.0
        total_climb_fuel_kg = 0.0
        total_climb_time_hr = 0.0

        max_climb_power_required_kW = 0.0
        min_power_margin_kW = float("inf")

        feasible_climb = True

        # Store representative values near worst power margin
        worst_margin_alt_m = None
        worst_margin_lapse = None
        worst_margin_required_power_kW = None
        worst_margin_available_power_kW = None

        for i in range(n_climb_steps):
            h0 = altitude_grid[i]
            h1 = altitude_grid[i + 1]
            h_mid = 0.5 * (h0 + h1)
            dh = h1 - h0

            state = isa_state(h_mid)
            rho = state["rho"]
            lapse = state["lapse"]

            # Along a steady climb:
            # Lift balances W cos(gamma)
            CL_climb, CDi_climb, CD_climb, D_climb = aero_at_condition(
                rho=rho,
                V=V_climb,
                W_normal=W * math.cos(gamma),
            )

            T_req_climb = D_climb + W * math.sin(gamma)

            P_thrust_climb_W = T_req_climb * V_climb
            P_elec_climb_kW = (P_thrust_climb_W / max(eta_chain, 1e-12)) / 1000.0

            P_generator_available_kW = P_generator_SL_rating_kW * lapse
            power_margin_kW = P_generator_available_kW - P_elec_climb_kW

            if P_elec_climb_kW > P_generator_available_kW:
                feasible_climb = False
                break

            dt_hr = (dh / ROC) / 3600.0
            dE_kWh = P_elec_climb_kW * dt_hr
            dfuel_kg = fuel_mass_from_electrical_energy_kwh(dE_kWh)

            total_climb_energy_kWh += dE_kWh
            total_climb_fuel_kg += dfuel_kg
            total_climb_time_hr += dt_hr

            max_climb_power_required_kW = max(
                max_climb_power_required_kW,
                P_elec_climb_kW,
            )

            if power_margin_kW < min_power_margin_kW:
                min_power_margin_kW = power_margin_kW
                worst_margin_alt_m = h_mid
                worst_margin_lapse = lapse
                worst_margin_required_power_kW = P_elec_climb_kW
                worst_margin_available_power_kW = P_generator_available_kW

        if not feasible_climb:
            continue

        fuel_mass_cruise_kg = (
            total_fuel_available_kg
            - total_climb_fuel_kg
            - descent_reserve_fraction_of_climb_fuel * total_climb_fuel_kg
        )

        if fuel_mass_cruise_kg <= 0.0:
            continue

        cruise_time_hr = fuel_mass_cruise_kg / max(fuel_flow_cruise_kg_per_hr, 1e-12)
        cruise_range_km = cruise_time_hr * 3600.0 * v_cruise / 1000.0

        climb_distance_km = cruise_altitude / max(math.tan(gamma), 1e-12) / 1000.0
        total_range_km = climb_distance_km + cruise_range_km

        candidate = {
            "V_climb_opt_mps": V_climb,
            "ROC_opt_mps": ROC,
            "climb_time_hr": total_climb_time_hr,
            "climb_energy_kWh": total_climb_energy_kWh,
            "fuel_mass_climb_kg": total_climb_fuel_kg,
            "max_climb_power_required_kW": max_climb_power_required_kW,
            "min_climb_power_margin_kW": min_power_margin_kW,
            "worst_margin_alt_m": worst_margin_alt_m,
            "worst_margin_lapse": worst_margin_lapse,
            "worst_margin_required_power_kW": worst_margin_required_power_kW,
            "worst_margin_available_power_kW": worst_margin_available_power_kW,
            "fuel_mass_cruise_kg": fuel_mass_cruise_kg,
            "cruise_time_hr": cruise_time_hr,
            "cruise_range_km": cruise_range_km,
            "climb_distance_km": climb_distance_km,
            "total_range_km": total_range_km,
        }

        # Minimum-energy climb speed among feasible speeds.
        if (best is None) or (candidate["climb_energy_kWh"] < best["climb_energy_kWh"]):
            best = candidate

    if best is None:
        return {
            "feasible": False,
            "reason": "No feasible climb speed found after accounting for altitude lapse during climb.",
            "P_generator_SL_rating_kW": P_generator_SL_rating_kW,
            "generator_mass_kg": generator_mass_kg,
            "fixed_propulsion_hardware_mass_kg": fixed_propulsion_hardware_mass_kg,
            "total_fuel_available_kg": total_fuel_available_kg,
            "P_elec_cruise_kW": P_elec_cruise_kW,
            "P_generator_cruise_available_kW": P_generator_cruise_available_kW,
            "cruise_power_margin_kW": cruise_power_margin_kW,
            "V_stall_climb_SL_mps": V_stall_climb_SL,
            "V_climb_min_mps": V_climb_min,
            "V_climb_max_mps": V_climb_max,
        }

    return {
        "feasible": True,

        # Atmosphere / lapse
        "rho_cruise_kgm3": rho_cruise,
        "theta_cruise": theta_cruise,
        "delta_cruise": delta_cruise,
        "lapse_cruise": lapse_cruise,

        # Aero
        "CL_cruise": CL_cruise,
        "CDi_cruise": CDi_cruise,
        "CD_cruise": CD_cruise,
        "k_induced": k_induced,
        "e_oswald": e_oswald,
        "thrust_required_cruise_N": D_cruise,
        "V_min_drag_mps": V_min_drag_mps,
        "wing_area": wing_area,

        # Generator
        "P_generator_SL_rating_kW": P_generator_SL_rating_kW,
        "P_generator_cruise_available_kW": P_generator_cruise_available_kW,
        "cruise_power_margin_kW": cruise_power_margin_kW,
        "generator_specific_power_kW_per_kg": generator_specific_power,
        "generator_mass_kg": generator_mass_kg,

        # Cruise power / fuel
        "P_elec_cruise_kW": P_elec_cruise_kW,
        "fuel_flow_cruise_kg_per_hr": fuel_flow_cruise_kg_per_hr,

        # Mass bookkeeping
        "fixed_propulsion_hardware_mass_kg": fixed_propulsion_hardware_mass_kg,
        "total_fuel_available_kg": total_fuel_available_kg,

        # Climb sweep bounds
        "V_stall_climb_SL_mps": V_stall_climb_SL,
        "V_climb_min_mps": V_climb_min,
        "V_climb_max_mps": V_climb_max,

        # Best climb and range
        **best,
    }


if __name__ == "__main__":
    result = get_range(
        mass=16550 / 2.2,
        Cd0=0.019,
        Cdi_cruise=0.003,          # legacy input, not used directly
        CLmax=6.1,
        AR=7.0,
        wing_area=45,
        v_cruise=125,
        takeoff_power=3250.0,      # sea-level turbogenerator rating [kW]
        climb_angle=15.0,          # deg
        cruise_altitude=18000.0 * 0.3048,
        mass_budget=2700,
        motor_mass=35,
        fan_mass=10,
        duct_mass=10,
        num_fans=12,
        e_oswald=0.80,
        n_climb_steps=300,
    )

    if not result["feasible"]:
        print("INFEASIBLE CASE")
        print(f"Reason: {result['reason']}")
        print()
        print("---- Generator / Mass ----")
        print(f"Generator sea-level rating: {result['P_generator_SL_rating_kW']:.1f} kW")
        print(f"Generator mass: {result['generator_mass_kg']:.1f} kg")
        print(f"Fixed propulsion hardware mass: {result['fixed_propulsion_hardware_mass_kg']:.1f} kg")
        print(f"Total fuel available: {result['total_fuel_available_kg']:.1f} kg")
        print()
        print("---- Cruise ----")
        print(f"Cruise electrical power required: {result['P_elec_cruise_kW']:.1f} kW")
        print(f"Generator power available at cruise: {result['P_generator_cruise_available_kW']:.1f} kW")
        print(f"Cruise power margin: {result['cruise_power_margin_kW']:.1f} kW")
        print()
        print("---- Climb Sweep ----")
        print(f"Sea-level climb stall speed: {result['V_stall_climb_SL_mps']:.1f} m/s")
        print(f"Climb speed sweep min: {result['V_climb_min_mps']:.1f} m/s")
        print(f"Climb speed sweep max: {result['V_climb_max_mps']:.1f} m/s")

    else:
        print("FEASIBLE CASE")
        print()

        print("---- Range ----")
        print(f"Total range: {result['total_range_km']:.1f} km")
        print(f"Cruise range: {result['cruise_range_km']:.1f} km")
        print(f"Climb distance: {result['climb_distance_km']:.1f} km")
        print(f"Cruise time: {result['cruise_time_hr']:.2f} hr")
        print()

        print("---- Cruise Aero / Power ----")
        print(f"Cruise electrical power required: {result['P_elec_cruise_kW']:.1f} kW")
        print(f"Generator power available at cruise: {result['P_generator_cruise_available_kW']:.1f} kW")
        print(f"Cruise power margin: {result['cruise_power_margin_kW']:.1f} kW")
        print(f"Cruise lapse factor delta/sqrt(theta): {result['lapse_cruise']:.3f}")
        print(f"Cruise fuel flow: {result['fuel_flow_cruise_kg_per_hr']:.1f} kg/hr")
        print(f"Minimum-drag cruise speed: {result['V_min_drag_mps']:.1f} m/s")
        print(f"Cruise CL: {result['CL_cruise']:.3f}")
        print(f"Cruise CDi: {result['CDi_cruise']:.4f}")
        print(f"Cruise CD: {result['CD_cruise']:.4f}")
        print()

        print("---- Generator / Mass ----")
        print(f"Generator sea-level rating: {result['P_generator_SL_rating_kW']:.1f} kW")
        print(f"Generator mass: {result['generator_mass_kg']:.1f} kg")
        print(f"Fixed propulsion hardware mass: {result['fixed_propulsion_hardware_mass_kg']:.1f} kg")
        print(f"Total fuel available before mission burn: {result['total_fuel_available_kg']:.1f} kg")
        print()

        print("---- Climb with Lapse ----")
        print(f"Sea-level climb stall speed: {result['V_stall_climb_SL_mps']:.1f} m/s")
        print(f"Climb speed sweep min: {result['V_climb_min_mps']:.1f} m/s")
        print(f"Climb speed sweep max: {result['V_climb_max_mps']:.1f} m/s")
        print(f"Minimum-energy climb speed: {result['V_climb_opt_mps']:.1f} m/s")
        print(f"Climb time: {result['climb_time_hr']:.3f} hr")
        print(f"Climb time: {result['climb_time_hr'] * 60:.1f} min")
        print(f"Climb energy: {result['climb_energy_kWh']:.1f} kWh")
        print(f"Climb fuel mass: {result['fuel_mass_climb_kg']:.1f} kg")
        print(f"Max climb electrical power required: {result['max_climb_power_required_kW']:.1f} kW")
        print(f"Minimum climb power margin: {result['min_climb_power_margin_kW']:.1f} kW")
        print(f"Worst margin altitude: {result['worst_margin_alt_m']:.0f} m")
        print(f"Worst margin lapse factor: {result['worst_margin_lapse']:.3f}")
        print(f"Power required at worst margin: {result['worst_margin_required_power_kW']:.1f} kW")
        print(f"Power available at worst margin: {result['worst_margin_available_power_kW']:.1f} kW")
        print()

        print("---- Fuel Bookkeeping ----")
        print(f"Fuel available before mission burn: {result['total_fuel_available_kg']:.1f} kg")
        print(f"Fuel burned in climb: {result['fuel_mass_climb_kg']:.1f} kg")
        print(f"Fuel remaining for cruise: {result['fuel_mass_cruise_kg']:.1f} kg")