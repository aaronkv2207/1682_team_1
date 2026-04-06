import numpy as np
import matplotlib.pyplot as plt

from range_model import get_range


# ============================================================
# Assumptions
# ============================================================

GEN_SPECIFIC_POWER_KW_PER_KG = 6.0

BATTERY_SPECIFIC_ENERGY_WH_PER_KG = 272.0
BATTERY_DOD = 0.80

THERMAL_EFFICIENCY = 0.30
LHV_MJ_PER_KG = 42.0


# ============================================================
# Helper functions
# ============================================================

def fuel_mass_from_electrical_energy_kwh(
    E_elec_kwh,
    eta_thermal=THERMAL_EFFICIENCY,
    LHV_MJ_per_kg=LHV_MJ_PER_KG,
):
    """
    Convert electrical energy [kWh] into fuel mass [kg].
    """
    E_MJ = E_elec_kwh * 3.6
    return E_MJ / max(eta_thermal * LHV_MJ_per_kg, 1e-12)


def battery_mass_from_energy_kwh(
    E_batt_kwh,
    battery_specific_energy_Wh_per_kg=BATTERY_SPECIFIC_ENERGY_WH_PER_KG,
    battery_dod=BATTERY_DOD,
):
    """
    Convert required battery energy [kWh] into battery mass [kg].

    Uses usable specific energy:
        specific_energy * depth_of_discharge
    """
    usable_specific_energy_kWh_per_kg = (
        battery_specific_energy_Wh_per_kg / 1000.0
    ) * battery_dod

    return E_batt_kwh / max(usable_specific_energy_kWh_per_kg, 1e-12)


# ============================================================
# Hybrid architecture model
# ============================================================

def get_range_battery_assist(
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
    Hybrid architecture:

    - generator is sized only to cruise electrical power
    - takeoff power and climb power are assumed equal
    - during climb:
        generator provides P_cruise
        battery provides P_takeoff - P_cruise
    - battery shortfall is sustained for the full climb duration
    - propulsion mass budget is fixed:
        mass_budget = fixed_hardware + generator + battery + fuel

    This function uses the baseline get_range() call as the mission/aero solver
    so the hybrid and baseline stay consistent on:
    - cruise power
    - climb time
    - climb distance
    """

    baseline_like = get_range(
        mass=mass,
        Cd0=Cd0,
        Cdi_cruise=Cdi_cruise,
        CLmax=CLmax,
        AR=AR,
        wing_area=wing_area,
        v_cruise=v_cruise,
        takeoff_power=takeoff_power,
        climb_angle=climb_angle,
        cruise_altitude=cruise_altitude,
        mass_budget=mass_budget,
        motor_mass=motor_mass,
        fan_mass=fan_mass,
        duct_mass=duct_mass,
        num_fans=num_fans,
    )

    if not baseline_like.get("feasible", False):
        return {
            "feasible": False,
            "reason": f"Reference mission solver infeasible: {baseline_like.get('reason', 'unknown')}",
        }

    # Cruise generator sizing
    P_cruise_kW = baseline_like["P_elec_cruise_kW"]
    generator_mass_kg = P_cruise_kW / max(GEN_SPECIFIC_POWER_KW_PER_KG, 1e-12)

    # Climb duration from same aero/mission logic as baseline
    climb_time_hr = baseline_like["climb_time_hr"]

    # Battery only covers the shortfall during climb
    P_batt_climb_kW = max(0.0, takeoff_power - P_cruise_kW)
    E_batt_climb_kWh = P_batt_climb_kW * climb_time_hr
    battery_mass_kg = battery_mass_from_energy_kwh(E_batt_climb_kWh)

    # Fixed propulsion hardware that exists in both architectures
    fixed_propulsion_mass_kg = num_fans * (motor_mass + fan_mass + duct_mass)

    # Dry hybrid propulsion mass
    dry_mass_hybrid_kg = (
        fixed_propulsion_mass_kg
        + generator_mass_kg
        + battery_mass_kg
    )

    # Fuel loaded inside the fixed mass budget
    total_fuel_available_kg = mass_budget - dry_mass_hybrid_kg

    if total_fuel_available_kg <= 0.0:
        return {
            "feasible": False,
            "reason": "No mass budget left for fuel after adding cruise generator and climb battery.",
            "generator_mass_kg": generator_mass_kg,
            "battery_mass_kg": battery_mass_kg,
            "fixed_propulsion_mass_kg": fixed_propulsion_mass_kg,
            "dry_mass_hybrid_kg": dry_mass_hybrid_kg,
            "total_fuel_available_kg": total_fuel_available_kg,
        }

    # During climb the generator is still operating at cruise power
    E_gen_climb_kWh = P_cruise_kW * climb_time_hr
    fuel_mass_climb_kg = fuel_mass_from_electrical_energy_kwh(E_gen_climb_kWh)

    # Remaining fuel after climb is available for cruise
    fuel_mass_cruise_kg = total_fuel_available_kg - fuel_mass_climb_kg

    if fuel_mass_cruise_kg <= 0.0:
        return {
            "feasible": False,
            "reason": "All available fuel is consumed during climb generator contribution.",
            "generator_mass_kg": generator_mass_kg,
            "battery_mass_kg": battery_mass_kg,
            "fixed_propulsion_mass_kg": fixed_propulsion_mass_kg,
            "dry_mass_hybrid_kg": dry_mass_hybrid_kg,
            "total_fuel_available_kg": total_fuel_available_kg,
            "fuel_mass_climb_kg": fuel_mass_climb_kg,
        }

    # Cruise fuel flow from baseline cruise condition
    fuel_flow_cruise_kg_per_hr = baseline_like["fuel_flow_cruise_kg_per_hr"]
    cruise_time_hr = fuel_mass_cruise_kg / max(fuel_flow_cruise_kg_per_hr, 1e-12)
    cruise_range_km = cruise_time_hr * 3600.0 * v_cruise / 1000.0

    total_range_km = cruise_range_km + baseline_like["climb_distance_km"]

    return {
        "feasible": True,
        "architecture": "cruise_generator_plus_climb_battery",

        # Reused performance quantities
        "P_elec_cruise_kW": P_cruise_kW,
        "climb_time_hr": climb_time_hr,
        "V_climb_opt_mps": baseline_like["V_climb_opt_mps"],
        "climb_distance_km": baseline_like["climb_distance_km"],

        # Generator / battery split during climb
        "P_gen_climb_kW": P_cruise_kW,
        "P_batt_climb_kW": P_batt_climb_kW,
        "E_gen_climb_kWh": E_gen_climb_kWh,
        "E_batt_climb_kWh": E_batt_climb_kWh,

        # Masses
        "generator_mass_kg": generator_mass_kg,
        "battery_mass_kg": battery_mass_kg,
        "fixed_propulsion_mass_kg": fixed_propulsion_mass_kg,
        "dry_mass_hybrid_kg": dry_mass_hybrid_kg,
        "total_fuel_available_kg": total_fuel_available_kg,

        # Fuel bookkeeping
        "fuel_mass_climb_kg": fuel_mass_climb_kg,
        "fuel_mass_cruise_kg": fuel_mass_cruise_kg,
        "fuel_flow_cruise_kg_per_hr": fuel_flow_cruise_kg_per_hr,

        # Range
        "cruise_time_hr": cruise_time_hr,
        "cruise_range_km": cruise_range_km,
        "total_range_km": total_range_km,
    }


# ============================================================
# Comparison wrapper
# ============================================================

def compare_architectures(**kwargs):
    """
    Compare:
    1) baseline: generator sized to takeoff/climb power
    2) hybrid: generator sized to cruise power, battery covers climb shortfall
    """
    baseline = get_range(**kwargs)
    hybrid = get_range_battery_assist(**kwargs)

    out = {
        "baseline": baseline,
        "hybrid": hybrid,
    }

    if baseline.get("feasible", False) and hybrid.get("feasible", False):
        out["delta_range_km"] = hybrid["total_range_km"] - baseline["total_range_km"]
        out["range_ratio"] = hybrid["total_range_km"] / max(baseline["total_range_km"], 1e-12)
    else:
        out["delta_range_km"] = np.nan
        out["range_ratio"] = np.nan

    return out


# ============================================================
# Sweep utility
# ============================================================

def sweep_parameter(base_kwargs, param_name, values):
    """
    Sweep one input parameter and compare baseline vs hybrid architectures.
    """
    results = {
        "param_name": param_name,
        "values": np.asarray(values, dtype=float),

        "range_baseline_km": [],
        "range_hybrid_km": [],
        "delta_range_km": [],

        "dry_mass_baseline_kg": [],
        "dry_mass_hybrid_kg": [],

        "fuel_mass_baseline_kg": [],
        "fuel_mass_hybrid_kg": [],

        "generator_mass_baseline_kg": [],
        "generator_mass_hybrid_kg": [],
        "battery_mass_hybrid_kg": [],

        "P_cruise_hybrid_kW": [],
        "P_batt_climb_hybrid_kW": [],
        "E_batt_climb_hybrid_kWh": [],
        "climb_time_hybrid_hr": [],
    }

    for value in values:
        case = dict(base_kwargs)
        case[param_name] = float(value)

        comp = compare_architectures(**case)
        baseline = comp["baseline"]
        hybrid = comp["hybrid"]

        if baseline.get("feasible", False):
            results["range_baseline_km"].append(baseline["total_range_km"])
            results["dry_mass_baseline_kg"].append(baseline["dry_propulsion_mass_kg"])
            results["fuel_mass_baseline_kg"].append(baseline["total_fuel_available_kg"])
            results["generator_mass_baseline_kg"].append(baseline["generator_mass_kg"])
        else:
            results["range_baseline_km"].append(np.nan)
            results["dry_mass_baseline_kg"].append(np.nan)
            results["fuel_mass_baseline_kg"].append(np.nan)
            results["generator_mass_baseline_kg"].append(np.nan)

        if hybrid.get("feasible", False):
            results["range_hybrid_km"].append(hybrid["total_range_km"])
            results["dry_mass_hybrid_kg"].append(hybrid["dry_mass_hybrid_kg"])
            results["fuel_mass_hybrid_kg"].append(hybrid["total_fuel_available_kg"])
            results["generator_mass_hybrid_kg"].append(hybrid["generator_mass_kg"])
            results["battery_mass_hybrid_kg"].append(hybrid["battery_mass_kg"])
            results["P_cruise_hybrid_kW"].append(hybrid["P_elec_cruise_kW"])
            results["P_batt_climb_hybrid_kW"].append(hybrid["P_batt_climb_kW"])
            results["E_batt_climb_hybrid_kWh"].append(hybrid["E_batt_climb_kWh"])
            results["climb_time_hybrid_hr"].append(hybrid["climb_time_hr"])
        else:
            results["range_hybrid_km"].append(np.nan)
            results["dry_mass_hybrid_kg"].append(np.nan)
            results["fuel_mass_hybrid_kg"].append(np.nan)
            results["generator_mass_hybrid_kg"].append(np.nan)
            results["battery_mass_hybrid_kg"].append(np.nan)
            results["P_cruise_hybrid_kW"].append(np.nan)
            results["P_batt_climb_hybrid_kW"].append(np.nan)
            results["E_batt_climb_hybrid_kWh"].append(np.nan)
            results["climb_time_hybrid_hr"].append(np.nan)

        results["delta_range_km"].append(comp["delta_range_km"])

    for k, v in results.items():
        if isinstance(v, list):
            results[k] = np.asarray(v, dtype=float)

    return results


# ============================================================
# Plotting
# ============================================================

def plot_delta_range_vs_parameter(results, xlabel=None):
    x = results["values"]
    name = results["param_name"]
    xlabel = xlabel or name

    plt.figure()
    plt.plot(x, results["delta_range_km"], linewidth=2)
    plt.axhline(0.0, linestyle="--")
    plt.xlabel(xlabel)
    plt.ylabel("Hybrid - baseline range (km)")
    plt.title(f"Delta Range vs {name}")
    plt.tight_layout()


def plot_range_comparison(results, xlabel=None):
    x = results["values"]
    name = results["param_name"]
    xlabel = xlabel or name

    plt.figure()
    plt.plot(x, results["range_baseline_km"], label="Baseline range")
    plt.plot(x, results["range_hybrid_km"], label="Hybrid range")
    plt.xlabel(xlabel)
    plt.ylabel("Total range (km)")
    plt.title(f"Range Comparison vs {name}")
    plt.legend()
    plt.tight_layout()


def plot_dry_mass_and_fuel_trade(results, xlabel=None):
    x = results["values"]
    name = results["param_name"]
    xlabel = xlabel or name

    plt.figure()
    plt.plot(x, results["dry_mass_baseline_kg"], label="Baseline dry propulsion mass")
    plt.plot(x, results["dry_mass_hybrid_kg"], label="Hybrid dry propulsion mass")
    plt.xlabel(xlabel)
    plt.ylabel("Dry mass (kg)")
    plt.title(f"Dry Propulsion Mass vs {name}")
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(x, results["fuel_mass_baseline_kg"], label="Baseline fuel mass")
    plt.plot(x, results["fuel_mass_hybrid_kg"], label="Hybrid fuel mass")
    plt.xlabel(xlabel)
    plt.ylabel("Fuel mass within fixed budget (kg)")
    plt.title(f"Fuel Mass Available vs {name}")
    plt.legend()
    plt.tight_layout()


def plot_component_masses(results, xlabel=None):
    x = results["values"]
    name = results["param_name"]
    xlabel = xlabel or name

    plt.figure()
    plt.plot(x, results["generator_mass_baseline_kg"], label="Baseline generator mass")
    plt.plot(x, results["generator_mass_hybrid_kg"], label="Hybrid generator mass")
    plt.plot(x, results["battery_mass_hybrid_kg"], label="Hybrid battery mass")
    plt.xlabel(xlabel)
    plt.ylabel("Mass (kg)")
    plt.title(f"Component Mass Breakdown vs {name}")
    plt.legend()
    plt.tight_layout()


def plot_hybrid_battery_burden(results, xlabel=None):
    x = results["values"]
    name = results["param_name"]
    xlabel = xlabel or name

    plt.figure()
    plt.plot(x, results["P_batt_climb_hybrid_kW"])
    plt.xlabel(xlabel)
    plt.ylabel("Battery climb power (kW)")
    plt.title(f"Hybrid Battery Power Requirement vs {name}")
    plt.tight_layout()

    plt.figure()
    plt.plot(x, results["E_batt_climb_hybrid_kWh"])
    plt.xlabel(xlabel)
    plt.ylabel("Battery climb energy (kWh)")
    plt.title(f"Hybrid Battery Energy Requirement vs {name}")
    plt.tight_layout()

    plt.figure()
    plt.plot(x, results["climb_time_hybrid_hr"])
    plt.xlabel(xlabel)
    plt.ylabel("Climb time (hr)")
    plt.title(f"Hybrid Climb Time vs {name}")
    plt.tight_layout()


def plot_all(results, xlabel=None):
    plot_range_comparison(results, xlabel=xlabel)
    plot_delta_range_vs_parameter(results, xlabel=xlabel)
    plot_dry_mass_and_fuel_trade(results, xlabel=xlabel)
    plot_component_masses(results, xlabel=xlabel)
    plot_hybrid_battery_burden(results, xlabel=xlabel)
    
# ============================================================
# Requirement-centered study tools
# ============================================================

TARGET_RANGE_KM = 2400.0
FIXED_PROPULSION_BUDGET_KG = 2700.0
FIXED_V_CRUISE_MPS = 125.0


def _evaluate_architecture(architecture_name, case):
    """
    architecture_name: "baseline" or "hybrid"
    """
    if architecture_name == "baseline":
        return get_range(**case)
    elif architecture_name == "hybrid":
        return get_range_battery_assist(**case)
    else:
        raise ValueError(f"Unknown architecture_name: {architecture_name}")


def sweep_requirement_margin_vs_power(
    base_case,
    power_vals,
    target_range_km=TARGET_RANGE_KM,
):
    """
    Sweep takeoff/climb power while holding v_cruise fixed at the design value.
    Returns range and range margin for baseline and hybrid.
    """
    results = {
        "power_kW": np.asarray(power_vals, dtype=float),

        "range_baseline_km": [],
        "range_hybrid_km": [],

        "margin_baseline_km": [],
        "margin_hybrid_km": [],

        "dry_mass_baseline_kg": [],
        "dry_mass_hybrid_kg": [],

        "fuel_mass_baseline_kg": [],
        "fuel_mass_hybrid_kg": [],

        "battery_mass_hybrid_kg": [],
    }

    for power in power_vals:
        case = dict(base_case)
        case["takeoff_power"] = float(power)
        case["v_cruise"] = FIXED_V_CRUISE_MPS
        case["mass_budget"] = FIXED_PROPULSION_BUDGET_KG

        baseline = get_range(**case)
        hybrid = get_range_battery_assist(**case)

        if baseline.get("feasible", False):
            Rb = baseline["total_range_km"]
            results["range_baseline_km"].append(Rb)
            results["margin_baseline_km"].append(Rb - target_range_km)
            results["dry_mass_baseline_kg"].append(baseline["dry_propulsion_mass_kg"])
            results["fuel_mass_baseline_kg"].append(baseline["total_fuel_available_kg"])
        else:
            results["range_baseline_km"].append(np.nan)
            results["margin_baseline_km"].append(np.nan)
            results["dry_mass_baseline_kg"].append(np.nan)
            results["fuel_mass_baseline_kg"].append(np.nan)

        if hybrid.get("feasible", False):
            Rh = hybrid["total_range_km"]
            results["range_hybrid_km"].append(Rh)
            results["margin_hybrid_km"].append(Rh - target_range_km)
            results["dry_mass_hybrid_kg"].append(hybrid["dry_mass_hybrid_kg"])
            results["fuel_mass_hybrid_kg"].append(hybrid["total_fuel_available_kg"])
            results["battery_mass_hybrid_kg"].append(hybrid["battery_mass_kg"])
        else:
            results["range_hybrid_km"].append(np.nan)
            results["margin_hybrid_km"].append(np.nan)
            results["dry_mass_hybrid_kg"].append(np.nan)
            results["fuel_mass_hybrid_kg"].append(np.nan)
            results["battery_mass_hybrid_kg"].append(np.nan)

    for k, v in results.items():
        if isinstance(v, list):
            results[k] = np.asarray(v, dtype=float)

    return results


def find_required_mass_budget(
    architecture_name,
    case,
    target_range_km=TARGET_RANGE_KM,
    lower_budget_kg=500.0,
    upper_budget_kg=4000.0,
    max_budget_kg=20000.0,
    tol_kg=1.0,
    max_iter=60,
):
    """
    Finds the minimum propulsion mass budget needed to hit target_range_km
    for the requested architecture.

    Returns np.nan if the target cannot be met even after expanding the budget.
    """
    test_case = dict(case)
    test_case["v_cruise"] = FIXED_V_CRUISE_MPS

    def meets_target(budget_kg):
        test_case["mass_budget"] = float(budget_kg)
        out = _evaluate_architecture(architecture_name, test_case)
        return out.get("feasible", False) and (out["total_range_km"] >= target_range_km)

    lo = lower_budget_kg
    hi = upper_budget_kg

    while (not meets_target(hi)) and (hi < max_budget_kg):
        hi *= 1.5

    if not meets_target(hi):
        return np.nan

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        if meets_target(mid):
            hi = mid
        else:
            lo = mid

        if abs(hi - lo) <= tol_kg:
            break

    return hi


def sweep_required_budget_vs_power(
    base_case,
    power_vals,
    target_range_km=TARGET_RANGE_KM,
):
    """
    Sweep takeoff/climb power and compute the minimum mass budget required
    for each architecture to meet the 2400 km requirement at 125 m/s.
    """
    results = {
        "power_kW": np.asarray(power_vals, dtype=float),
        "required_budget_baseline_kg": [],
        "required_budget_hybrid_kg": [],
        "budget_margin_baseline_kg": [],
        "budget_margin_hybrid_kg": [],
    }

    for power in power_vals:
        case = dict(base_case)
        case["takeoff_power"] = float(power)
        case["v_cruise"] = FIXED_V_CRUISE_MPS

        req_base = find_required_mass_budget(
            "baseline",
            case,
            target_range_km=target_range_km,
        )
        req_hyb = find_required_mass_budget(
            "hybrid",
            case,
            target_range_km=target_range_km,
        )

        results["required_budget_baseline_kg"].append(req_base)
        results["required_budget_hybrid_kg"].append(req_hyb)
        results["budget_margin_baseline_kg"].append(FIXED_PROPULSION_BUDGET_KG - req_base if np.isfinite(req_base) else np.nan)
        results["budget_margin_hybrid_kg"].append(FIXED_PROPULSION_BUDGET_KG - req_hyb if np.isfinite(req_hyb) else np.nan)

    for k, v in results.items():
        if isinstance(v, list):
            results[k] = np.asarray(v, dtype=float)

    return results


def best_case_envelope_vs_power(
    base_case,
    power_vals,
    CLmax_vals,
    wing_area_vals,
    altitude_vals,
    target_range_km=TARGET_RANGE_KM,
):
    """
    For each takeoff/climb power, compute the best achievable range for
    baseline and hybrid over a plausible aerodynamic envelope.

    This is the strongest version of the argument:
    if the best hybrid range never reaches 2400 km, batteries are ruled out
    over the entire chosen sweep space.
    """
    results = {
        "power_kW": np.asarray(power_vals, dtype=float),

        "best_range_baseline_km": [],
        "best_range_hybrid_km": [],

        "best_margin_baseline_km": [],
        "best_margin_hybrid_km": [],

        "best_CLmax_baseline": [],
        "best_CLmax_hybrid": [],
        "best_wing_area_baseline": [],
        "best_wing_area_hybrid": [],
        "best_altitude_baseline_m": [],
        "best_altitude_hybrid_m": [],
    }

    for power in power_vals:
        best_R_base = -np.inf
        best_R_hyb = -np.inf

        best_cfg_base = (np.nan, np.nan, np.nan)
        best_cfg_hyb = (np.nan, np.nan, np.nan)

        for CLmax in CLmax_vals:
            for S in wing_area_vals:
                for alt in altitude_vals:
                    case = dict(base_case)
                    case["takeoff_power"] = float(power)
                    case["CLmax"] = float(CLmax)
                    case["wing_area"] = float(S)
                    case["cruise_altitude"] = float(alt)
                    case["v_cruise"] = FIXED_V_CRUISE_MPS
                    case["mass_budget"] = FIXED_PROPULSION_BUDGET_KG

                    base_out = get_range(**case)
                    if base_out.get("feasible", False):
                        Rb = base_out["total_range_km"]
                        if Rb > best_R_base:
                            best_R_base = Rb
                            best_cfg_base = (CLmax, S, alt)

                    hyb_out = get_range_battery_assist(**case)
                    if hyb_out.get("feasible", False):
                        Rh = hyb_out["total_range_km"]
                        if Rh > best_R_hyb:
                            best_R_hyb = Rh
                            best_cfg_hyb = (CLmax, S, alt)

        results["best_range_baseline_km"].append(best_R_base if np.isfinite(best_R_base) else np.nan)
        results["best_range_hybrid_km"].append(best_R_hyb if np.isfinite(best_R_hyb) else np.nan)

        results["best_margin_baseline_km"].append(best_R_base - target_range_km if np.isfinite(best_R_base) else np.nan)
        results["best_margin_hybrid_km"].append(best_R_hyb - target_range_km if np.isfinite(best_R_hyb) else np.nan)

        results["best_CLmax_baseline"].append(best_cfg_base[0])
        results["best_CLmax_hybrid"].append(best_cfg_hyb[0])
        results["best_wing_area_baseline"].append(best_cfg_base[1])
        results["best_wing_area_hybrid"].append(best_cfg_hyb[1])
        results["best_altitude_baseline_m"].append(best_cfg_base[2])
        results["best_altitude_hybrid_m"].append(best_cfg_hyb[2])

    for k, v in results.items():
        if isinstance(v, list):
            results[k] = np.asarray(v, dtype=float)

    return results


# ============================================================
# Plotting for requirement study
# ============================================================

def plot_requirement_margin_vs_power(results):
    x = results["power_kW"]

    plt.figure()
    plt.plot(x, results["margin_baseline_km"], label="Baseline margin")
    plt.plot(x, results["margin_hybrid_km"], label="Hybrid margin")
    plt.axhline(0.0, linestyle="--")
    plt.xlabel("Takeoff / climb power (kW)")
    plt.ylabel("Range margin to 2400 km (km)")
    plt.title("Range Margin to Requirement at 125 m/s")
    plt.legend()
    plt.tight_layout()


def plot_required_budget_vs_power(results):
    x = results["power_kW"]

    plt.figure()
    plt.plot(x, results["required_budget_baseline_kg"], label="Baseline required budget")
    plt.plot(x, results["required_budget_hybrid_kg"], label="Hybrid required budget")
    plt.axhline(FIXED_PROPULSION_BUDGET_KG, linestyle="--", label="Available budget = 2700 kg")
    plt.xlabel("Takeoff / climb power (kW)")
    plt.ylabel("Required propulsion mass budget (kg)")
    plt.title("Minimum Budget Required to Reach 2400 km at 125 m/s")
    plt.legend()
    plt.tight_layout()


def plot_best_case_envelope(results):
    x = results["power_kW"]

    plt.figure()
    plt.plot(x, results["best_range_baseline_km"], label="Baseline best-case range")
    plt.plot(x, results["best_range_hybrid_km"], label="Hybrid best-case range")
    plt.axhline(TARGET_RANGE_KM, linestyle="--", label="Requirement = 2400 km")
    plt.xlabel("Takeoff / climb power (kW)")
    plt.ylabel("Best achievable range (km)")
    plt.title("Best-Case Range Envelope at 125 m/s and 2700 kg Budget")
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(x, results["best_margin_baseline_km"], label="Baseline best-case margin")
    plt.plot(x, results["best_margin_hybrid_km"], label="Hybrid best-case margin")
    plt.axhline(0.0, linestyle="--")
    plt.xlabel("Takeoff / climb power (kW)")
    plt.ylabel("Best-case margin to 2400 km (km)")
    plt.title("Best-Case Requirement Margin")
    plt.legend()
    plt.tight_layout()



# ============================================================
# Example run
# ============================================================


if __name__ == "__main__":
    base_case = {
        "mass": 16550 / 2.2,
        "Cd0": 0.019,
        "Cdi_cruise": 0.003,
        "CLmax": 6.1,
        "AR": 8.0,
        "wing_area": (16550 / 2.2) / (1484.56 / 9.81),
        "v_cruise": FIXED_V_CRUISE_MPS,
        "takeoff_power": 2500.0,
        "climb_angle": 10.0,
        "cruise_altitude": 18000.0 * 0.3048,
        "mass_budget": FIXED_PROPULSION_BUDGET_KG,
        "motor_mass": 35,
        "fan_mass": 10,
        "duct_mass": 35,
        "num_fans": 8,
    }

    # 1) Direct requirement margin sweep
    power_vals = np.linspace(500.0, 4000.0, 80)
    margin_sweep = sweep_requirement_margin_vs_power(base_case, power_vals)
    plot_requirement_margin_vs_power(margin_sweep)

    # 2) Required budget sweep
    budget_sweep = sweep_required_budget_vs_power(base_case, power_vals)
    plot_required_budget_vs_power(budget_sweep)

    # 3) Best-case envelope over plausible aero ranges
    #    Start coarse, then refine if needed.
    CLmax_vals = np.linspace(5.0, 8.0, 7)
    wing_area_vals = np.linspace(base_case["wing_area"] * 0.85, base_case["wing_area"] * 1.20, 7)
    altitude_vals = np.linspace(4000.0, 7000.0, 7)

    envelope = best_case_envelope_vs_power(
        base_case,
        power_vals=np.linspace(1000.0, 3500.0, 40),
        CLmax_vals=CLmax_vals,
        wing_area_vals=wing_area_vals,
        altitude_vals=altitude_vals,
    )
    plot_best_case_envelope(envelope)

    plt.show()

# if __name__ == "__main__":
#     base_case = {
#         "mass": 16550 / 2.2,
#         "Cd0": 0.019,
#         "Cdi_cruise": 0.003,
#         "CLmax": 6.1,
#         "AR": 8.0,
#         "wing_area": (16550 / 2.2) / (1484.56 / 9.81),
#         "v_cruise": 125.0,
#         "takeoff_power": 200.0,             # kW, also used as climb power
#         "climb_angle": 10.0,                 # deg
#         "cruise_altitude": 18000.0 * 0.3048, # m
#         "mass_budget": 2700,
#         "motor_mass": 35,
#         "fan_mass": 10,
#         "duct_mass": 35,
#         "num_fans": 8,
#     }

#     # --------------------------------
#     # Single-point comparison
#     # --------------------------------
#     comp = compare_architectures(**base_case)

#     print("\n=== Single Point Comparison ===")

#     if comp["baseline"].get("feasible", False):
#         print(f"Baseline range:                 {comp['baseline']['total_range_km']:.1f} km")
#         print(f"Baseline dry propulsion mass:   {comp['baseline']['dry_propulsion_mass_kg']:.1f} kg")
#         print(f"Baseline fuel mass:             {comp['baseline']['total_fuel_available_kg']:.1f} kg")
#         print(f"Baseline generator mass:        {comp['baseline']['generator_mass_kg']:.1f} kg")
#     else:
#         print("Baseline infeasible:", comp["baseline"].get("reason", "unknown"))

#     if comp["hybrid"].get("feasible", False):
#         print(f"Hybrid range:                   {comp['hybrid']['total_range_km']:.1f} km")
#         print(f"Hybrid dry propulsion mass:     {comp['hybrid']['dry_mass_hybrid_kg']:.1f} kg")
#         print(f"Hybrid fuel mass:               {comp['hybrid']['total_fuel_available_kg']:.1f} kg")
#         print(f"Hybrid generator mass:          {comp['hybrid']['generator_mass_kg']:.1f} kg")
#         print(f"Hybrid battery mass:            {comp['hybrid']['battery_mass_kg']:.1f} kg")
#         print(f"Hybrid cruise power:            {comp['hybrid']['P_elec_cruise_kW']:.1f} kW")
#         print(f"Hybrid battery climb power:     {comp['hybrid']['P_batt_climb_kW']:.1f} kW")
#         print(f"Hybrid battery climb energy:    {comp['hybrid']['E_batt_climb_kWh']:.1f} kWh")
#         print(f"Hybrid climb time:              {comp['hybrid']['climb_time_hr']:.3f} hr")
#     else:
#         print("Hybrid infeasible:", comp["hybrid"].get("reason", "unknown"))

#     print(f"Delta range (hybrid - baseline): {comp['delta_range_km']:.1f} km")

#     # --------------------------------
#     # Main sweep: delta range vs takeoff/climb power
#     # --------------------------------
#     sweep_takeoff_power = sweep_parameter(
#         base_case,
#         "takeoff_power",
#         np.linspace(500.0, 4000.0, 80),
#     )
#     plot_all(sweep_takeoff_power, xlabel="Takeoff / climb power (kW)")

#     # --------------------------------
#     # Optional extra sweeps
#     # Uncomment any of these as needed
#     # --------------------------------

#     # sweep_clmax = sweep_parameter(base_case, "CLmax", np.linspace(4.0, 8.0, 50))
#     # plot_all(sweep_clmax, xlabel="CLmax")

#     # sweep_area = sweep_parameter(base_case, "wing_area", np.linspace(35.0, 70.0, 50))
#     # plot_all(sweep_area, xlabel="Wing area (m^2)")

#     # sweep_v = sweep_parameter(base_case, "v_cruise", np.linspace(90.0, 170.0, 50))
#     # plot_all(sweep_v, xlabel="Cruise speed (m/s)")

#     # sweep_alt = sweep_parameter(base_case, "cruise_altitude", np.linspace(2000.0, 8000.0, 50))
#     # plot_all(sweep_alt, xlabel="Cruise altitude (m)")

#     plt.show()