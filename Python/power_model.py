"""
DEP sizing + takeoff metric + cruise fuel burn + parameter sweeps w/ contour plots.

Updates:
- Parabolic drag model: CD = CD0 + k*CL^2
  (k can be provided directly or computed from AR, e)
- Steady climb model with climb angle gamma (parameter, default 15 deg)
- Generator sized by *cruise* electrical power requirement (with margin)
- Batteries provide any *excess* electrical power above generator nameplate
  (e.g., in climb)
- Added climb-speed sweep to find minimum total climb energy
- Added plots of CL and CD vs climb speed

Notes:
- LHV should be ~42 MJ/kg (not kJ/kg).
- eta_apu_overall is "fuel -> electrical bus" overall efficiency.
- battery_specific_mass_kg_per_kWh is in kg/kWh
"""

import math
import numpy as np
import matplotlib.pyplot as plt

'''
cruise Cd0 = 0.016
induced = .003
total = .019

takeoff Cd0 = .02
induced = 1.645
total = 1.665
'''

# -------------------------
# Atmosphere (ISA, up to 11 km)
# -------------------------
def isa_density(alt_m: float) -> float:
    g0 = 9.80665
    R = 287.05287
    T0 = 288.15
    p0 = 101325.0
    L = 0.0065

    alt_m = max(0.0, alt_m)
    if alt_m > 11000.0:
        raise ValueError("isa_density only implemented up to 11,000 m")

    T = T0 - L * alt_m
    p = p0 * (T / T0) ** (g0 / (R * L))
    return p / (R * T)


# -------------------------
# Drag polar helper
# -------------------------
def induced_factor_k(params: dict) -> float:
    """
    Returns k in CD = CD0 + k*CL^2.

    You can either provide:
      - params["k_induced"] directly, OR
      - params["AR"] and params["e"] and we compute k = 1/(pi*e*AR)
    """
    if "k_induced" in params and params["k_induced"] is not None:
        return float(params["k_induced"])

    AR = float(params.get("AR", 0.0))
    e = float(params.get("e", 0.0))
    if AR > 0.0 and e > 0.0:
        return 1.0 / (math.pi * e * AR)

    raise KeyError("Provide either k_induced or (AR and e) for the parabolic drag model.")


def CD_parabolic(CL: float, CD0: float, k: float) -> float:
    return CD0 + k * CL * CL


# -------------------------
# Core calculations
# -------------------------
def generator_kW_and_xto(params: dict) -> dict:
    # Mission
    R_m = float(params["range_m"])
    V = float(params["V_cruise"])
    V_climb = float(params["V_climb"])
    alt_m = float(params["alt_cruise_m"])

    # Climb
    gamma_deg = float(params.get("climb_angle_deg", 15.0))
    gamma = math.radians(gamma_deg)

    # Aero
    CD0 = float(params["CD0"])
    k = induced_factor_k(params)

    # Aircraft/design variables
    m_kg = float(params["mass_kg"])
    g = float(params["g"])
    m_over_S = float(params["wing_loading_kgm2"])   # kg/m^2
    T_over_W = float(params["thrust_loading_TW"])   # -

    # Takeoff
    CLmax = float(params["CLmax"])
    rho_to = params.get("rho_to_kgm3", None)
    if rho_to is None:
        rho_to = isa_density(0.0)

    # Efficiencies + margin
    eta_prop = float(params["eta_prop"])
    eta_motor = float(params["eta_motor"])
    eta_inv = float(params["eta_inverter"])
    eta_gb = float(params["eta_gearbox"])
    eta_gen = float(params["eta_generator"])
    margin = float(params["power_margin"])

    # Battery specific mass (kg/kWh)
    battery_specific_mass = float(params["battery_specific_mass_kg_per_kWh"])

    # 1) Wing area from mass wing loading
    S = m_kg / m_over_S   # m^2

    # Weights and dynamic pressure
    W = m_kg * g
    rho_cr = isa_density(alt_m)
    q = 0.5 * rho_cr * V**2

    # -------------------------
    # 2) CRUISE (level flight)
    # -------------------------
    CL_cruise = W / max(q * S, 1e-12)
    CD_cruise = CD_parabolic(CL_cruise, CD0=CD0, k=k)
    D_cruise = q * S * CD_cruise
    T_req_cruise = D_cruise

    # -------------------------
    # 3) Cruise electrical power
    # Generator sized by cruise load with margin
    # -------------------------
    P_thrust_cruise = T_req_cruise * V
    P_shaft_cruise = P_thrust_cruise / max(eta_prop, 1e-12)
    eta_down = max(eta_motor * eta_inv * eta_gb, 1e-12)
    P_bus_cruise = P_shaft_cruise / eta_down
    P_gen_elec_cruise = P_bus_cruise / max(eta_gen, 1e-12)
    P_gen_nameplate = 2500000 #P_gen_elec_cruise * margin

    # -------------------------
    # 4) CLIMB (steady climb at fixed gamma)
    # Perpendicular: L = W cos(gamma)
    # Along path:    T = D + W sin(gamma)
    # -------------------------
    q_climb = 0.5 * rho_cr * V_climb**2
    CL_climb = (W * math.cos(gamma)) / max(q_climb * S, 1e-12)
    CD_climb = CD_parabolic(CL_climb, CD0=CD0, k=k)
    D_climb = q_climb * S * CD_climb
    T_req_climb = D_climb + W * math.sin(gamma)

    P_thrust_climb = T_req_climb * V_climb
    P_shaft_climb = P_thrust_climb / max(eta_prop, 1e-12)
    P_bus_climb = P_shaft_climb / eta_down
    P_elec_req_climb = P_bus_climb / max(eta_gen, 1e-12)

    # Batteries cover only excess beyond generator nameplate
    P_batt_climb = max(0.0, P_elec_req_climb - P_gen_nameplate)
    P_gen_used_climb = min(P_gen_nameplate, P_elec_req_climb)

    # Climb kinematics
    ROC_mps = V_climb * math.sin(gamma)
    climb_time_hr = (alt_m / max(ROC_mps, 1e-12)) / 3600.0

    # Battery sizing for climb excess only
    battery_kwh = P_batt_climb / 1000.0 * climb_time_hr
    battery_mass = battery_specific_mass * battery_kwh

    # -------------------------
    # 5) Takeoff roll metric
    # x_to = (W/S)/(T/W) * (1/(rho*g))*(1/CLmax)
    # -------------------------
    W_over_S = m_over_S * g
    x_to = (
        (W_over_S / max(T_over_W, 1e-12))
        * (1.0 / (rho_to * g))
        * (1.0 / max(CLmax, 1e-12))
    )

    # -------------------------
    # 6) Mission cruise energy
    # -------------------------
    t_cruise = R_m / max(V, 1e-12)
    E_cruise_J = P_gen_elec_cruise * t_cruise
    E_cruise_MWh = E_cruise_J / 3.6e9

    return {
        # Atmosphere / geometry
        "rho_cruise_kgm3": rho_cr,
        "wing_area_m2": S,

        # Drag polar
        "CD0": CD0,
        "k_induced": k,

        # Cruise aero + power
        "CL_cruise": CL_cruise,
        "CD_cruise": CD_cruise,
        "thrust_required_cruise_N": T_req_cruise,
        "P_gen_elec_kW": P_gen_elec_cruise / 1000.0,
        "P_gen_sized_kW": P_gen_nameplate / 1000.0,

        # Climb aero + power
        "climb_angle_deg": gamma_deg,
        "ROC_mps": ROC_mps,
        "climb_time_hr": climb_time_hr,
        "CL_climb": CL_climb,
        "CD_climb": CD_climb,
        "thrust_required_climb_N": T_req_climb,
        "P_elec_required_climb_kW": P_elec_req_climb / 1000.0,
        "P_gen_used_climb_kW": P_gen_used_climb / 1000.0,
        "P_batt_climb_kW": P_batt_climb / 1000.0,
        "battery_energy_climb_kWh": battery_kwh,
        "battery_mass_kg": battery_mass,

        # Takeoff + mission
        "x_to_metric": x_to,
        "cruise_time_hr": t_cruise / 3600.0,
        "cruise_energy_MWh": E_cruise_MWh,
    }


# -------------------------
# Fuel burn from cruise electrical energy
# -------------------------
def cruise_fuel_from_energy(
    E_elec_MWh: float,
    eta_apu: float,
    LHV_MJ_per_kg: float = 42.0,
) -> float:
    E_MJ = E_elec_MWh * 3600.0
    return E_MJ / max(eta_apu * LHV_MJ_per_kg, 1e-12)


def add_cruise_fuel_metrics(out: dict, params: dict) -> dict:
    eta_apu = float(params.get("eta_apu_overall", 0.30))
    LHV = float(params.get("LHV_MJ_per_kg", 42.0))

    fuel_kg = cruise_fuel_from_energy(
        out["cruise_energy_MWh"],
        eta_apu=eta_apu,
        LHV_MJ_per_kg=LHV,
    )
    t_s = out["cruise_time_hr"] * 3600.0

    out2 = dict(out)
    out2["fuel_mass_cruise_kg"] = fuel_kg
    out2["fuel_flow_cruise_kg_s"] = fuel_kg / max(t_s, 1e-12)
    out2["fuel_flow_cruise_kg_hr"] = out2["fuel_flow_cruise_kg_s"] * 3600.0
    return out2


# -------------------------
# Objective function
# -------------------------
def objective_function(
    out: dict,
    params: dict,
    x_to_ref: float,
    V_ref: float = 125.0,
    w_to: float = 1.0,
    w_V: float = 1.0,
    wing_loading_limit: float = 250.0,
    penalty_weight: float = 50.0,
) -> float:
    x_to_hat = out["x_to_metric"] / max(x_to_ref, 1e-12)
    V_hat = float(params["V_cruise"]) / max(V_ref, 1e-12)

    m_over_S = float(params["wing_loading_kgm2"])
    if m_over_S <= wing_loading_limit:
        penalty = 0.0
    else:
        penalty = penalty_weight * ((m_over_S - wing_loading_limit) / wing_loading_limit) ** 2

    return w_to * x_to_hat - w_V * V_hat + penalty


def add_objective_metric(
    out: dict,
    params: dict,
    x_to_ref: float,
    V_ref: float = 125.0,
    w_to: float = 1.0,
    w_V: float = 1.0,
    wing_loading_limit: float = 250.0,
    penalty_weight: float = 50.0,
) -> dict:
    out2 = dict(out)
    out2["objective_J"] = objective_function(
        out,
        params,
        x_to_ref=x_to_ref,
        V_ref=V_ref,
        w_to=w_to,
        w_V=w_V,
        wing_loading_limit=wing_loading_limit,
        penalty_weight=penalty_weight,
    )
    return out2


# -------------------------
# Generic 2D sweep + contour plot
# -------------------------
def sweep_2d(
    params_base: dict,
    x_key: str,
    x_vals: np.ndarray,
    y_key: str,
    y_vals: np.ndarray,
    z_key: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    X, Y = np.meshgrid(x_vals, y_vals)
    Z = np.full_like(X, np.nan, dtype=float)

    base_out = add_cruise_fuel_metrics(generator_kW_and_xto(params_base), params_base)
    x_to_ref = base_out["x_to_metric"]

    for i in range(Y.shape[0]):
        for j in range(X.shape[1]):
            p = dict(params_base)
            p[x_key] = float(X[i, j])
            p[y_key] = float(Y[i, j])

            out = generator_kW_and_xto(p)
            out = add_cruise_fuel_metrics(out, p)
            out = add_objective_metric(
                out,
                p,
                x_to_ref=x_to_ref,
                V_ref=params_base.get("V_ref", 125.0),
                w_to=params_base.get("w_to", 1.0),
                w_V=params_base.get("w_V", 1.0),
                wing_loading_limit=params_base.get("wing_loading_limit", 250.0),
                penalty_weight=params_base.get("penalty_weight", 50.0),
            )

            if z_key not in out:
                raise KeyError(f"z_key='{z_key}' not found. Available: {list(out.keys())}")

            Z[i, j] = float(out[z_key])

    return X, Y, Z


def contour_plot(
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    x_label: str,
    y_label: str,
    title: str,
    levels: int = 20,
    add_limit_line: dict | None = None,
):
    plt.figure()
    cs = plt.contourf(X, Y, Z, levels=levels)
    plt.colorbar(cs, label=title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    if add_limit_line is not None:
        if add_limit_line["orientation"] == "v":
            plt.axvline(add_limit_line["value"], linestyle="--")
            if add_limit_line.get("label"):
                plt.text(
                    add_limit_line["value"],
                    np.nanmin(Y),
                    add_limit_line["label"],
                    rotation=90,
                    va="bottom",
                )
        elif add_limit_line["orientation"] == "h":
            plt.axhline(add_limit_line["value"], linestyle="--")
            if add_limit_line.get("label"):
                plt.text(
                    np.nanmin(X),
                    add_limit_line["value"],
                    add_limit_line["label"],
                    va="bottom",
                )

    plt.tight_layout()


# -------------------------
# Climb speed sweep
# -------------------------
def climb_energy_at_speed(params: dict, V_climb: float) -> dict:
    """
    Compute climb performance/energy for a single climb speed while holding
    climb angle gamma fixed.

    Hybrid climb power split:
    - generator and battery both contribute during climb
    - battery supplies a fixed fraction of climb electrical power
    """

    p = dict(params)
    p["V_climb"] = float(V_climb)

    # Mission / aircraft
    alt_m = float(p["alt_cruise_m"])
    gamma_deg = float(p.get("climb_angle_deg", 15.0))
    gamma = math.radians(gamma_deg)

    m_kg = float(p["mass_kg"])
    g = float(p["g"])
    W = m_kg * g

    m_over_S = float(p["wing_loading_kgm2"])
    S = m_kg / m_over_S

    # Atmosphere
    rho = isa_density(alt_m)

    # Aero
    CD0 = float(p["CD0"])
    k = induced_factor_k(p)

    # Efficiencies
    eta_prop = float(p["eta_prop"])
    eta_motor = float(p["eta_motor"])
    eta_inv = float(p["eta_inverter"])
    eta_gb = float(p["eta_gearbox"])
    eta_gen = float(p["eta_generator"])

    # Hybrid split
    battery_fraction = float(p.get("battery_power_fraction_climb", 0.30))
    battery_fraction = min(max(battery_fraction, 0.0), 1.0)

    # Generator sized off cruise
    base_out = generator_kW_and_xto(p)
    P_gen_nameplate = base_out["P_gen_sized_kW"] * 1000.0  # W

    # Climb state
    q = 0.5 * rho * V_climb**2
    CL = (W * math.cos(gamma)) / max(q * S, 1e-12)
    CD = CD_parabolic(CL, CD0=CD0, k=k)
    D = q * S * CD
    T_req = D + W * math.sin(gamma)

    # Power chain
    P_thrust = T_req * V_climb
    P_shaft = P_thrust / max(eta_prop, 1e-12)
    eta_down = max(eta_motor * eta_inv * eta_gb, 1e-12)
    P_bus = P_shaft / eta_down
    P_elec_req = P_bus / max(eta_gen, 1e-12)

    # Fixed hybrid split
    P_batt = battery_fraction * P_elec_req
    P_gen_used = (1.0 - battery_fraction) * P_elec_req

    # Optional overload check
    generator_overload = max(0.0, P_gen_used - P_gen_nameplate)

    # Climb duration
    ROC = V_climb * math.sin(gamma)
    climb_time_s = alt_m / max(ROC, 1e-12)

    # Energies
    E_total_kWh = P_elec_req * climb_time_s / 3.6e6
    E_batt_kWh = P_batt * climb_time_s / 3.6e6
    E_gen_kWh = P_gen_used * climb_time_s / 3.6e6

    return {
        "V_climb_mps": V_climb,
        "q_Pa": q,
        "CL_climb": CL,
        "CD_climb": CD,
        "D_climb_N": D,
        "T_req_climb_N": T_req,
        "P_elec_req_climb_kW": P_elec_req / 1000.0,
        "P_batt_climb_kW": P_batt / 1000.0,
        "P_gen_used_climb_kW": P_gen_used / 1000.0,
        "P_gen_nameplate_kW": P_gen_nameplate / 1000.0,
        "generator_overload_kW": generator_overload / 1000.0,
        "ROC_mps": ROC,
        "climb_time_s": climb_time_s,
        "climb_time_hr": climb_time_s / 3600.0,
        "E_total_climb_kWh": E_total_kWh,
        "E_batt_climb_kWh": E_batt_kWh,
        "E_gen_climb_kWh": E_gen_kWh,
    }
    
def battery_mass_for_climb_and_generator(
    params: dict,
    V_climb: float,
    P_gen_nameplate_kW: float,
) -> dict:
    """
    For a chosen climb speed and generator nameplate power, compute:
    - climb electrical power required
    - climb time
    - battery power required
    - battery energy required
    - battery mass required
    - generator mass required
    - climb fuel mass required for the generator contribution

    Assumes generator can contribute up to P_gen_nameplate_kW continuously in climb.
    Battery supplies any remaining shortfall.
    """

    p = dict(params)
    p["V_climb"] = float(V_climb)

    alt_m = float(p["alt_cruise_m"])
    gamma_deg = float(p.get("climb_angle_deg", 15.0))
    gamma = math.radians(gamma_deg)

    m_kg = float(p["mass_kg"])
    g = float(p["g"])
    W = m_kg * g

    m_over_S = float(p["wing_loading_kgm2"])
    S = m_kg / m_over_S

    rho = isa_density(alt_m)

    CD0 = float(p["CD0"])
    k = induced_factor_k(p)

    eta_prop = float(p["eta_prop"])
    eta_motor = float(p["eta_motor"])
    eta_inv = float(p["eta_inverter"])
    eta_gb = float(p["eta_gearbox"])
    eta_gen = float(p["eta_generator"])

    battery_specific_mass = float(p["battery_specific_mass_kg_per_kWh"])

    # Fuel-to-electric conversion assumptions
    eta_apu = float(p.get("eta_apu_overall", 0.30))
    LHV_MJ_per_kg = float(p.get("LHV_MJ_per_kg", 42.0))

    q = 0.5 * rho * V_climb**2
    CL = (W * math.cos(gamma)) / max(q * S, 1e-12)
    CD = CD_parabolic(CL, CD0=CD0, k=k)
    D = q * S * CD
    T_req = D + W * math.sin(gamma)

    P_thrust = T_req * V_climb
    P_shaft = P_thrust / max(eta_prop, 1e-12)
    eta_down = max(eta_motor * eta_inv * eta_gb, 1e-12)
    P_bus = P_shaft / eta_down
    P_elec_req_kW = (P_bus / max(eta_gen, 1e-12)) / 1000.0

    ROC_mps = V_climb * math.sin(gamma)
    climb_time_hr = (alt_m / max(ROC_mps, 1e-12)) / 3600.0

    P_gen_used_kW = min(P_gen_nameplate_kW, P_elec_req_kW)
    P_batt_kW = max(0.0, P_elec_req_kW - P_gen_nameplate_kW)

    E_batt_kWh = P_batt_kW * climb_time_hr
    E_gen_kWh = P_gen_used_kW * climb_time_hr
    E_total_kWh = P_elec_req_kW * climb_time_hr

    battery_mass_kg = battery_specific_mass * E_batt_kWh

    # Fuel mass required for generator-provided climb energy
    fuel_mass_climb_kg = fuel_mass_from_electrical_energy_kWh(
        E_gen_kWh,
        eta_apu=eta_apu,
        LHV_MJ_per_kg=LHV_MJ_per_kg,
    )

    out = {
        "V_climb_mps": V_climb,
        "P_gen_nameplate_kW": P_gen_nameplate_kW,
        "CL_climb": CL,
        "CD_climb": CD,
        "thrust_required_climb_N": T_req,
        "P_elec_required_climb_kW": P_elec_req_kW,
        "P_gen_used_climb_kW": P_gen_used_kW,
        "P_batt_climb_kW": P_batt_kW,
        "ROC_mps": ROC_mps,
        "climb_time_hr": climb_time_hr,
        "E_total_climb_kWh": E_total_kWh,
        "E_gen_climb_kWh": E_gen_kWh,
        "E_batt_climb_kWh": E_batt_kWh,
        "battery_mass_kg": battery_mass_kg,
        "fuel_mass_climb_kg": fuel_mass_climb_kg,
    }

    # Generator mass from specific power (kW/kg)
    if "generator_specific_power_kW_per_kg" in p:
        gen_specific_power = float(p["generator_specific_power_kW_per_kg"])
        generator_mass_kg = P_gen_nameplate_kW / max(gen_specific_power, 1e-12)
        out["generator_mass_kg"] = generator_mass_kg

        # Old metric
        out["total_gen_plus_batt_mass_kg"] = generator_mass_kg + battery_mass_kg

        # New mission-aware metric
        out["total_gen_plus_batt_plus_fuel_mass_kg"] = (
            generator_mass_kg + battery_mass_kg + fuel_mass_climb_kg
        )

    return out  
def sweep_climb_speed_and_generator(
    params: dict,
    V_climb_vals: np.ndarray,
    P_gen_vals_kW: np.ndarray,
    z_key: str = "battery_mass_kg",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    2D sweep over climb speed and generator nameplate power.

    Returns:
    X = climb speed mesh
    Y = generator power mesh
    Z = selected output metric

    Useful z_key values:
    - "battery_mass_kg"
    - "generator_mass_kg"
    - "total_gen_plus_batt_mass_kg"
    - "E_batt_climb_kWh"
    - "P_batt_climb_kW"
    - "P_elec_required_climb_kW"
    """

    X, Y = np.meshgrid(V_climb_vals, P_gen_vals_kW)
    Z = np.full_like(X, np.nan, dtype=float)

    for i in range(Y.shape[0]):
        for j in range(X.shape[1]):
            Vc = float(X[i, j])
            Pg = float(Y[i, j])

            out = battery_mass_for_climb_and_generator(params, Vc, Pg)

            if z_key not in out:
                raise KeyError(f"z_key='{z_key}' not found. Available: {list(out.keys())}")

            Z[i, j] = float(out[z_key])

    return X, Y, Z

def plot_climb_generator_contour(
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    title: str,
    z_label: str,
):
    plt.figure()
    cs = plt.contourf(X, Y, Z, levels=25)
    plt.colorbar(cs, label=z_label)
    plt.xlabel("Climb speed (m/s)")
    plt.ylabel("Generator nameplate power (kW)")
    plt.title(title)
    plt.tight_layout()

def plot_power_breakdown_vs_climb_speed(
    params: dict,
    V_climb_vals: np.ndarray,
    P_gen_nameplate_kW: float,
):
    """
    Plot instantaneous climb power breakdown vs climb speed for a fixed generator size.
    """

    V_arr = np.asarray(V_climb_vals, dtype=float)

    P_total = np.full_like(V_arr, np.nan, dtype=float)
    P_gen = np.full_like(V_arr, np.nan, dtype=float)
    P_batt = np.full_like(V_arr, np.nan, dtype=float)
    CL_arr = np.full_like(V_arr, np.nan, dtype=float)
    CD_arr = np.full_like(V_arr, np.nan, dtype=float)
    E_total = np.full_like(V_arr, np.nan, dtype=float)

    for i, Vc in enumerate(V_arr):
        out = battery_mass_for_climb_and_generator(params, Vc, P_gen_nameplate_kW)

        P_total[i] = out["P_elec_required_climb_kW"]
        P_gen[i] = out["P_gen_used_climb_kW"]
        P_batt[i] = out["P_batt_climb_kW"]
        CL_arr[i] = out["CL_climb"]
        CD_arr[i] = out["CD_climb"]
        E_total[i] = out["E_total_climb_kWh"]

    i_opt = int(np.nanargmin(E_total))
    V_opt = V_arr[i_opt]

    plt.figure()
    plt.plot(V_arr, P_total, label="Total climb electrical power")
    plt.plot(V_arr, P_gen, label="Generator power")
    plt.plot(V_arr, P_batt, label="Battery power")
    plt.axvline(V_opt, linestyle="--", label=f"Min-energy climb speed = {V_opt:.1f} m/s")
    plt.xlabel("Climb speed (m/s)")
    plt.ylabel("Power (kW)")
    plt.title(f"Instantaneous Climb Power vs Climb Speed\nGenerator = {P_gen_nameplate_kW:.1f} kW")
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(V_arr, CL_arr)
    plt.axvline(V_opt, linestyle="--")
    plt.xlabel("Climb speed (m/s)")
    plt.ylabel("CL")
    plt.title("CL vs Climb Speed")
    plt.tight_layout()

    plt.figure()
    plt.plot(V_arr, CD_arr)
    plt.axvline(V_opt, linestyle="--")
    plt.xlabel("Climb speed (m/s)")
    plt.ylabel("CD")
    plt.title("CD vs Climb Speed")
    plt.tight_layout()
    
    
def plot_mass_trade_vs_generator(
    params: dict,
    V_climb: float,
    P_gen_vals_kW: np.ndarray,
) -> dict:
    """
    Plot battery mass, generator mass, climb fuel mass, and total mission mass
    vs generator size at a fixed climb speed.

    Returns a dictionary containing the full sweep data and the optimal point.
    """

    batt_mass = []
    gen_mass = []
    fuel_mass = []
    total_mass_old = []
    total_mass_new = []
    batt_power = []
    gen_energy = []

    for Pg in P_gen_vals_kW:
        out = battery_mass_for_climb_and_generator(params, V_climb, Pg)
        batt_mass.append(out["battery_mass_kg"])
        gen_mass.append(out["generator_mass_kg"])
        fuel_mass.append(out["fuel_mass_climb_kg"])
        total_mass_old.append(out["total_gen_plus_batt_mass_kg"])
        total_mass_new.append(out["total_gen_plus_batt_plus_fuel_mass_kg"])
        batt_power.append(out["P_batt_climb_kW"])
        gen_energy.append(out["E_gen_climb_kWh"])

    batt_mass = np.asarray(batt_mass)
    gen_mass = np.asarray(gen_mass)
    fuel_mass = np.asarray(fuel_mass)
    total_mass_old = np.asarray(total_mass_old)
    total_mass_new = np.asarray(total_mass_new)
    batt_power = np.asarray(batt_power)
    gen_energy = np.asarray(gen_energy)

    i_best = int(np.nanargmin(total_mass_new))
    Pg_best = float(P_gen_vals_kW[i_best])

    gen_mass_best = float(gen_mass[i_best])
    fuel_mass_best = float(fuel_mass[i_best])
    batt_mass_best = float(batt_mass[i_best])
    total_mass_best = float(total_mass_new[i_best])

    plt.figure()
    plt.plot(P_gen_vals_kW, batt_mass, label="Battery mass")
    plt.plot(P_gen_vals_kW, gen_mass, label="Generator mass")
    plt.plot(P_gen_vals_kW, fuel_mass, label="Climb fuel mass")
    plt.plot(P_gen_vals_kW, total_mass_old, label="Generator + battery mass")
    plt.plot(P_gen_vals_kW, total_mass_new, label="Generator + battery + climb fuel mass")
    plt.axvline(Pg_best, linestyle="--", label=f"Best mission mass = {Pg_best:.1f} kW")
    plt.xlabel("Generator nameplate power (kW)")
    plt.ylabel("Mass (kg)")
    plt.title(f"Mass Trade vs Generator Size at V_climb = {V_climb:.1f} m/s")
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(P_gen_vals_kW, batt_power)
    plt.axvline(Pg_best, linestyle="--")
    plt.xlabel("Generator nameplate power (kW)")
    plt.ylabel("Battery climb power required (kW)")
    plt.title(f"Battery Power Shortfall vs Generator Size at V_climb = {V_climb:.1f} m/s")
    plt.tight_layout()

    plt.figure()
    plt.plot(P_gen_vals_kW, gen_energy)
    plt.axvline(Pg_best, linestyle="--")
    plt.xlabel("Generator nameplate power (kW)")
    plt.ylabel("Generator-supplied climb energy (kWh)")
    plt.title(f"Generator Climb Energy vs Generator Size at V_climb = {V_climb:.1f} m/s")
    plt.tight_layout()

    return {
        "P_gen_vals_kW": np.asarray(P_gen_vals_kW),
        "battery_mass_kg": batt_mass,
        "generator_mass_kg": gen_mass,
        "fuel_mass_climb_kg": fuel_mass,
        "total_mass_old_kg": total_mass_old,
        "total_mass_new_kg": total_mass_new,
        "battery_power_kW": batt_power,
        "generator_energy_kWh": gen_energy,
        "i_best": i_best,
        "P_gen_best_kW": Pg_best,
        "battery_mass_best_kg": batt_mass_best,
        "generator_mass_best_kg": gen_mass_best,
        "fuel_mass_best_kg": fuel_mass_best,
        "total_mass_best_kg": total_mass_best,
    }
    
    
def get_cruise_baseline_generator_kW(params: dict) -> float:
    out = generator_kW_and_xto(params)
    return float(out["P_gen_sized_kW"])

def fuel_mass_from_electrical_energy_kWh(
    E_elec_kWh: float,
    eta_apu: float,
    LHV_MJ_per_kg: float = 42.0,
) -> float:
    """
    Convert required electrical energy (kWh on bus/generator side)
    into fuel mass, using overall fuel->electrical efficiency.
    """
    E_MJ = E_elec_kWh * 3.6
    return E_MJ / max(eta_apu * LHV_MJ_per_kg, 1e-12)

def sweep_climb_speed_energy(
    params: dict,
    V_climb_vals: np.ndarray,
    plot: bool = True,
) -> dict:
    """
    Sweep climb speed and identify the minimum total electrical climb energy.

    Plots:
      1) Total climb energy vs climb speed
      2) CL vs climb speed
      3) CD vs climb speed
      4) Instantaneous electrical power vs climb speed
         (total, generator, battery)
    """

    V_arr = np.asarray(V_climb_vals, dtype=float)

    E_total = np.full_like(V_arr, np.nan, dtype=float)
    E_batt = np.full_like(V_arr, np.nan, dtype=float)
    E_gen = np.full_like(V_arr, np.nan, dtype=float)

    P_total = np.full_like(V_arr, np.nan, dtype=float)
    P_batt = np.full_like(V_arr, np.nan, dtype=float)
    P_gen = np.full_like(V_arr, np.nan, dtype=float)

    ROC = np.full_like(V_arr, np.nan, dtype=float)
    t_hr = np.full_like(V_arr, np.nan, dtype=float)
    CL_arr = np.full_like(V_arr, np.nan, dtype=float)
    CD_arr = np.full_like(V_arr, np.nan, dtype=float)

    for i, Vc in enumerate(V_arr):

        r = climb_energy_at_speed(params, Vc)

        E_total[i] = r["E_total_climb_kWh"]
        E_batt[i] = r["E_batt_climb_kWh"]
        E_gen[i] = r["E_gen_climb_kWh"]

        P_total[i] = r["P_elec_req_climb_kW"]
        P_batt[i] = r["P_batt_climb_kW"]
        P_gen[i] = r["P_gen_used_climb_kW"]

        ROC[i] = r["ROC_mps"]
        t_hr[i] = r["climb_time_hr"]

        CL_arr[i] = r["CL_climb"]
        CD_arr[i] = r["CD_climb"]

    # Find minimum climb energy speed
    i_opt = int(np.nanargmin(E_total))
    V_opt = V_arr[i_opt]

    results = {
        "V_climb_mps": V_arr,
        "E_total_climb_kWh": E_total,
        "E_batt_climb_kWh": E_batt,
        "E_gen_climb_kWh": E_gen,
        "P_total_kW": P_total,
        "P_batt_kW": P_batt,
        "P_gen_kW": P_gen,
        "ROC_mps": ROC,
        "climb_time_hr": t_hr,
        "CL_climb": CL_arr,
        "CD_climb": CD_arr,
        "V_opt_mps": V_opt,
        "V_opt_kts": V_opt * 1.94384,
        "E_total_opt_kWh": E_total[i_opt],
        "E_batt_opt_kWh": E_batt[i_opt],
        "E_gen_opt_kWh": E_gen[i_opt],
        "CL_opt": CL_arr[i_opt],
        "CD_opt": CD_arr[i_opt],
        "ROC_opt_mps": ROC[i_opt],
        "climb_time_opt_hr": t_hr[i_opt],
    }

    if plot:

        # --------------------------------
        # Climb Energy Plot
        # --------------------------------
        plt.figure()
        plt.plot(V_arr, E_total, label="Total climb energy")
        plt.plot(V_arr, E_gen, label="Generator energy")
        plt.plot(V_arr, E_batt, label="Battery energy")
        plt.axvline(V_opt, linestyle="--", label=f"Optimal = {V_opt:.1f} m/s")
        plt.xlabel("Climb speed (m/s)")
        plt.ylabel("Energy (kWh)")
        plt.title("Climb Energy vs Climb Speed")
        plt.legend()
        plt.tight_layout()

        # --------------------------------
        # CL vs climb speed
        # --------------------------------
        plt.figure()
        plt.plot(V_arr, CL_arr)
        plt.axvline(V_opt, linestyle="--")
        plt.xlabel("Climb speed (m/s)")
        plt.ylabel("CL")
        plt.title("CL vs Climb Speed")
        plt.tight_layout()

        # --------------------------------
        # CD vs climb speed
        # --------------------------------
        plt.figure()
        plt.plot(V_arr, CD_arr)
        plt.axvline(V_opt, linestyle="--")
        plt.xlabel("Climb speed (m/s)")
        plt.ylabel("CD")
        plt.title("CD vs Climb Speed")
        plt.tight_layout()

        # --------------------------------
        # ⭐ NEW: Instantaneous Power Plot
        # --------------------------------
        plt.figure()
        plt.plot(V_arr, P_total, label="Total climb electrical power")
        plt.plot(V_arr, P_gen, label="Generator power")
        plt.plot(V_arr, P_batt, label="Battery power")
        plt.axvline(V_opt, linestyle="--", label=f"Optimal = {V_opt:.1f} m/s")
        plt.xlabel("Climb speed (m/s)")
        plt.ylabel("Power (kW)")
        plt.title("Instantaneous Climb Power vs Climb Speed")
        plt.legend()
        plt.tight_layout()

    return results


# -------------------------
# Example run
# -------------------------
if __name__ == "__main__":
    params = {
        # Mission
        "range_m": 2000e3,
        "V_cruise": 125.0,
        "alt_cruise_m": 18000.0 * 0.3048,

        # Drag polar
        "CD0": 0.019,
        "AR": 8.0,
        "e": 0.80,
        # "k_induced": 0.04,

        # Climb model
        "climb_angle_deg": 10,
        "V_climb": 82,

        # Aircraft/design variables
        "mass_kg": 16550 / 2.2,
        "wing_loading_kgm2": 1484.56 / 9.81,
        "thrust_loading_TW": 0.5,

        # Takeoff
        "CLmax": 6.1,

        # Constants
        "g": 9.80665,

        # Efficiencies + margin
        "eta_prop": 0.8,
        "eta_motor": 0.95,
        "eta_inverter": 0.98,
        "eta_gearbox": 1.0,
        "eta_generator": 0.95,
        "power_margin": 1.15,

        # Fuel/APU conversion
        "eta_apu_overall": 0.30,
        "LHV_MJ_per_kg": 42.0,

        # Battery
        "battery_specific_mass_kg_per_kWh": 8.0,
        "generator_specific_power_kW_per_kg": 6,

        # Objective tuning
        "V_ref": 125.0,
        "w_to": 1.0,
        "w_V": 1.0,
        "wing_loading_limit": 250.0,
        "penalty_weight": 50.0,
    }

    base_out = add_cruise_fuel_metrics(generator_kW_and_xto(params), params)
    x_to_ref = base_out["x_to_metric"]
    out = add_objective_metric(
        base_out,
        params,
        x_to_ref=x_to_ref,
        V_ref=params["V_ref"],
        w_to=params["w_to"],
        w_V=params["w_V"],
        wing_loading_limit=params["wing_loading_limit"],
        penalty_weight=params["penalty_weight"],
    )

    print("=== DEP Generator Sizing + Parabolic Drag + Climb + Cruise Fuel + Objective ===")
    print(f"Cruise density rho:              {out['rho_cruise_kgm3']:.3f} kg/m^3")
    print(f"Wing area S:                    {out['wing_area_m2']:.2f} m^2")
    print(f"Cruise CL, CD:                  {out['CL_cruise']:.3f}, {out['CD_cruise']:.4f}")
    print(f"Cruise thrust required:         {out['thrust_required_cruise_N']:.0f} N")
    print(f"Generator req (cruise):         {out['P_gen_elec_kW']:.1f} kW")
    print(f"Generator nameplate (sized):    {out['P_gen_sized_kW']:.1f} kW")

    print(f"Climb angle gamma:              {out['climb_angle_deg']:.1f} deg")
    print(f"Climb ROC (from gamma, V):      {out['ROC_mps']:.2f} m/s")
    print(f"Climb time:                     {out['climb_time_hr']:.2f} hr")
    print(f"Climb CL, CD:                   {out['CL_climb']:.3f}, {out['CD_climb']:.4f}")
    print(f"Climb thrust required:          {out['thrust_required_climb_N']:.0f} N")
    print(f"Elec required (climb):          {out['P_elec_required_climb_kW']:.1f} kW")
    print(f"Gen used (climb):               {out['P_gen_used_climb_kW']:.1f} kW")
    print(f"Battery make-up (climb):        {out['P_batt_climb_kW']:.1f} kW")
    print(f"Battery climb energy:           {out['battery_energy_climb_kWh']:.1f} kWh")
    print(f"Battery mass:                   {out['battery_mass_kg']:.1f} kg")

    print(f"x_to metric:                    {out['x_to_metric']:.4f}")
    print(f"Cruise time:                    {out['cruise_time_hr']:.2f} hr")
    print(f"Cruise energy:                  {out['cruise_energy_MWh']:.3f} MWh")
    print(f"Cruise fuel mass:               {out['fuel_mass_cruise_kg']:.1f} kg")
    print(f"Avg cruise fuel flow:           {out['fuel_flow_cruise_kg_hr']:.1f} kg/hr")
    print(f"Objective J (lower=better):     {out['objective_J']:.3f}")
    
    

    # Climb speed sweep
    V_climb_sweep = np.linspace(40.0, 140.0, 120)
    climb_sweep = sweep_climb_speed_energy(params, V_climb_sweep, plot=True)

    print("\n=== Climb Speed Sweep ===")
    print(f"Optimal climb speed (min total energy): {climb_sweep['V_opt_mps']:.2f} m/s ({climb_sweep['V_opt_kts']:.1f} kt)")
    print(f"Minimum total climb energy:             {climb_sweep['E_total_opt_kWh']:.2f} kWh")
    print(f"Battery climb energy at optimum:        {climb_sweep['E_batt_opt_kWh']:.2f} kWh")
    print(f"CL at optimum:                          {climb_sweep['CL_opt']:.3f}")
    print(f"CD at optimum:                          {climb_sweep['CD_opt']:.4f}")
    print(f"ROC at optimum:                         {climb_sweep['ROC_opt_mps']:.2f} m/s")
    print(f"Climb time at optimum:                  {climb_sweep['climb_time_opt_hr']:.2f} hr")

    # Assumptions
    params["battery_specific_mass_kg_per_kWh"] = 5.0
    params["generator_specific_power_kW_per_kg"] = 6

    # Cruise baseline generator size
    cruise_baseline_kW = get_cruise_baseline_generator_kW(params)
    print(f"Cruise baseline generator size = {cruise_baseline_kW:.1f} kW")

    # Sweep ranges
    V_climb_vals = np.linspace(40.0, 140.0, 80)
    P_gen_vals_kW = np.linspace(0.6 * cruise_baseline_kW, 2.0 * cruise_baseline_kW, 80)

    # Contour: battery mass
    X, Y, Z_batt = sweep_climb_speed_and_generator(
        params,
        V_climb_vals,
        P_gen_vals_kW,
        z_key="battery_mass_kg",
    )
    # plot_climb_generator_contour(
    #     X, Y, Z_batt,
    #     title="Battery Mass vs Climb Speed and Generator Power",
    #     z_label="Battery mass (kg)",
    # )

    # Contour: total generator + battery + climb fuel mass
    X, Y, Z_total = sweep_climb_speed_and_generator(
        params,
        V_climb_vals,
        P_gen_vals_kW,
        z_key="total_gen_plus_batt_plus_fuel_mass_kg",
    )
    # plot_climb_generator_contour(
    #     X, Y, Z_total,
    #     title="Generator + Battery + Climb Fuel Mass vs Climb Speed and Generator Power",
    #     z_label="Total mission mass (kg)",
    # )
    plot_climb_generator_contour(
        X, Y, Z_total,
        title="Generator + Battery Mass vs Climb Speed and Generator Power",
        z_label="Total mass (kg)",
    )

    # Power breakdown vs climb speed at cruise baseline generator size
    # plot_power_breakdown_vs_climb_speed(
    #     params,
    #     V_climb_vals,
    #     P_gen_nameplate_kW=cruise_baseline_kW,
    # )

    # Mass trade vs generator size at a chosen climb speed
    trade = plot_mass_trade_vs_generator(
        params,
        V_climb=82,
        P_gen_vals_kW=P_gen_vals_kW,
    )
    plt.show()
    
    print("Optimal generator size (kW):", trade["P_gen_best_kW"])
    print("Optimal generator mass (kg):", trade["generator_mass_best_kg"])
    print("Optimal fuel mass (kg):", trade["fuel_mass_best_kg"])
    print("Optimal battery mass (kg):", trade["battery_mass_best_kg"])
    print("Optimal total mass (kg):", trade["total_mass_best_kg"])
    
    mass_budget = 2700 #kg
    motor_mass = 35 #kg
    number_fans = 8
    duct_mass_factor = 1.5
    duct_mass = motor_mass * duct_mass_factor
    generator_mass = trade["generator_mass_best_kg"]
    dry_mass = generator_mass + number_fans * (duct_mass + motor_mass)
    wet_mass = mass_budget - dry_mass
    cruise_fuel = wet_mass - trade["fuel_mass_best_kg"]
    cruise_time_hours = cruise_fuel / out["fuel_flow_cruise_kg_hr"]
    final_range = cruise_time_hours * 3600 * params["V_cruise"] / 1000
    print(f"Dry Propulsion Mass: {dry_mass} kg")
    print(f"Total Range = {final_range} km")
    plt.show()