from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np

# ============================================================
# Data containers
# ============================================================

@dataclass
class Aircraft:
    mass: float          # kg
    Ixx: float           # kg m^2
    Iyy: float           # kg m^2
    Izz: float           # kg m^2
    Ixz: float           # kg m^2

    S: float             # m^2
    b: float             # m
    cbar: float          # m

    rho: float           # kg/m^3
    g: float = 9.80665   # m/s^2


@dataclass
class TrimPoint:
    U0: float            # m/s
    theta0: float        # rad


# ============================================================
# Helpers
# ============================================================

def qbar(rho: float, U0: float) -> float:
    return 0.5 * rho * U0**2


def _cget(d: Dict[str, float], key: str) -> float:
    """Safe getter with default 0.0."""
    return float(d.get(key, 0.0))


# ============================================================
# Longitudinal A matrix
#
# States: x_lon = [u, w, q, theta]^T
#
# Expected derivative dictionary keys:
#   CX_u,  CX_alpha, CX_q
#   CZ_u,  CZ_alpha, CZ_q
#   Cm_u,  Cm_alpha, Cm_q
#
# Assumptions:
#   - alpha ~= w / U0
#   - q-derivatives are with respect to qhat = q*cbar/(2U0)
#   - *_u derivatives are wrt dimensional u [m/s]
# ============================================================

def build_A_longitudinal(
    ac: Aircraft,
    trim: TrimPoint,
    derivs: Dict[str, float],
) -> np.ndarray:
    U0 = trim.U0
    theta0 = trim.theta0
    Q = qbar(ac.rho, U0)

    # coefficient -> translational acceleration scaling
    Fscale = Q * ac.S / ac.mass

    # coefficient -> pitch angular acceleration scaling
    Mscale = Q * ac.S * ac.cbar / ac.Iyy

    Xu = Fscale * _cget(derivs, "CX_u")
    Xw = Fscale * _cget(derivs, "CX_alpha") / U0
    Xq = Fscale * _cget(derivs, "CX_q") * (ac.cbar / (2.0 * U0))

    Zu = Fscale * _cget(derivs, "CZ_u")
    Zw = Fscale * _cget(derivs, "CZ_alpha") / U0
    Zq = Fscale * _cget(derivs, "CZ_q") * (ac.cbar / (2.0 * U0))

    Mu = Mscale * _cget(derivs, "Cm_u")
    Mw = Mscale * _cget(derivs, "Cm_alpha") / U0
    Mq = Mscale * _cget(derivs, "Cm_q") * (ac.cbar / (2.0 * U0))

    A_lon = np.array([
        [Xu, Xw, Xq, -ac.g * np.cos(theta0)],
        [Zu, Zw, U0 + Zq, -ac.g * np.sin(theta0)],
        [Mu, Mw, Mq, 0.0],
        [0.0, 0.0, 1.0, 0.0],
    ], dtype=float)

    return A_lon


# ============================================================
# Lateral-directional A matrix
#
# States: x_lat = [v, p, r, phi]^T
#
# Expected derivative dictionary keys:
#   CY_beta, CY_p, CY_r
#   Cl_beta, Cl_p, Cl_r
#   Cn_beta, Cn_p, Cn_r
#
# Assumptions:
#   - beta ~= v / U0
#   - p/r derivatives are wrt phat = p*b/(2U0), rhat = r*b/(2U0)
#   - starred derivatives optionally include Ixz coupling
# ============================================================

def build_A_lateral(
    ac: Aircraft,
    trim: TrimPoint,
    derivs: Dict[str, float],
    include_Ixz_coupling: bool = True,
) -> np.ndarray:
    U0 = trim.U0
    theta0 = trim.theta0
    Q = qbar(ac.rho, U0)

    Fscale = Q * ac.S / ac.mass
    Lscale = Q * ac.S * ac.b / ac.Ixx
    Nscale = Q * ac.S * ac.b / ac.Izz

    Yv = Fscale * _cget(derivs, "CY_beta") / U0
    Yp = Fscale * _cget(derivs, "CY_p") * (ac.b / (2.0 * U0))
    Yr = Fscale * _cget(derivs, "CY_r") * (ac.b / (2.0 * U0))

    Lv = Lscale * _cget(derivs, "Cl_beta") / U0
    Lp = Lscale * _cget(derivs, "Cl_p") * (ac.b / (2.0 * U0))
    Lr = Lscale * _cget(derivs, "Cl_r") * (ac.b / (2.0 * U0))

    Nv = Nscale * _cget(derivs, "Cn_beta") / U0
    Np = Nscale * _cget(derivs, "Cn_p") * (ac.b / (2.0 * U0))
    Nr = Nscale * _cget(derivs, "Cn_r") * (ac.b / (2.0 * U0))

    if include_Ixz_coupling and abs(ac.Ixz) > 1e-12:
        den = ac.Ixx * ac.Izz - ac.Ixz**2

        Lv_star = (ac.Izz * Lv + ac.Ixz * Nv) / den
        Lp_star = (ac.Izz * Lp + ac.Ixz * Np) / den
        Lr_star = (ac.Izz * Lr + ac.Ixz * Nr) / den

        Nv_star = (ac.Ixz * Lv + ac.Ixx * Nv) / den
        Np_star = (ac.Ixz * Lp + ac.Ixx * Np) / den
        Nr_star = (ac.Ixz * Lr + ac.Ixx * Nr) / den

        Lv, Lp, Lr = Lv_star, Lp_star, Lr_star
        Nv, Np, Nr = Nv_star, Np_star, Nr_star

    A_lat = np.array([
        [Yv, Yp, Yr - U0, ac.g * np.cos(theta0)],
        [Lv, Lp, Lr, 0.0],
        [Nv, Np, Nr, 0.0],
        [0.0, 1.0, np.tan(theta0), 0.0],
    ], dtype=float)

    return A_lat

####################

# from __future__ import annotations

# from typing import Dict, List
# import numpy as np


# ============================================================
# Eigen-analysis and labeling
# ============================================================

def eig_summary(A: np.ndarray, state_names: List[str]) -> List[Dict]:
    eigvals, eigvecs = np.linalg.eig(A)

    modes = []
    for i, lam in enumerate(eigvals):
        sigma = float(np.real(lam))
        omega_d = float(np.imag(lam))
        omega_n = float(np.sqrt(sigma**2 + omega_d**2))

        if abs(omega_d) > 1e-12:
            zeta = -sigma / omega_n
            period_s = 2.0 * np.pi / abs(omega_d)
        else:
            zeta = np.nan
            period_s = np.inf

        if sigma < 0.0:
            half_or_double_time_s = np.log(2.0) / (-sigma)
        elif sigma > 0.0:
            half_or_double_time_s = np.log(2.0) / sigma
        else:
            half_or_double_time_s = np.inf

        vec = eigvecs[:, i]
        vec_norm = vec / np.max(np.abs(vec))

        mode = {
            "eigval": lam,
            "real": sigma,
            "imag": omega_d,
            "omega_n": omega_n,
            "zeta": zeta,
            "period_s": period_s,
            "time_const_like_s": half_or_double_time_s,
            "vector": vec,
            "vector_norm": vec_norm,
            "state_names": state_names,
        }
        modes.append(mode)

    modes.sort(key=lambda m: (abs(m["imag"]) < 1e-8, abs(m["imag"]), m["real"]))
    return modes


def classify_longitudinal_modes(modes: List[Dict]) -> Dict[str, Dict]:
    complex_positive_imag = [m for m in modes if m["imag"] > 1e-8]
    if len(complex_positive_imag) < 2:
        return {}

    complex_positive_imag.sort(key=lambda m: abs(m["imag"]))
    return {
        "phugoid": complex_positive_imag[0],
        "short_period": complex_positive_imag[-1],
    }


def classify_lateral_modes(modes: List[Dict]) -> Dict[str, Dict]:
    out: Dict[str, Dict] = {}

    complex_positive_imag = [m for m in modes if m["imag"] > 1e-8]
    real_modes = [m for m in modes if abs(m["imag"]) <= 1e-8]

    if complex_positive_imag:
        out["dutch_roll"] = max(complex_positive_imag, key=lambda m: abs(m["imag"]))

    if real_modes:
        out["roll"] = min(real_modes, key=lambda m: m["real"])
        out["spiral"] = min(real_modes, key=lambda m: abs(m["real"]))

    return out


def print_mode(mode_name: str, mode: Dict) -> None:
    lam = mode["eigval"]
    print(f"{mode_name}:")
    print(f"  eigenvalue      = {lam.real:+.6f} {lam.imag:+.6f}j")
    print(f"  omega_n [rad/s] = {mode['omega_n']:.6f}")
    print(f"  zeta            = {mode['zeta']:.6f}")

    if np.isfinite(mode["period_s"]):
        print(f"  period [s]      = {mode['period_s']:.6f}")
    else:
        print("  period [s]      = inf (real mode)")

    print(f"  half/double [s] = {mode['time_const_like_s']:.6f}")
    print("  eigenvector (normalized to max abs = 1):")

    for name, value in zip(mode["state_names"], mode["vector_norm"]):
        print(f"    {name:>6s} : {value.real:+.6f} {value.imag:+.6f}j")
    print()


def print_labeled_modes(title: str, labeled_modes: Dict[str, Dict]) -> None:
    print(f"\n=== {title} ===\n")
    for mode_name, mode in labeled_modes.items():
        print_mode(mode_name, mode)

##########################################

# from __future__ import annotations

# import numpy as np


# ============================================================
# Sample run with example derivative values
#
# Replace `get_current_stability_derivatives()` with your
# operating-point database lookup.
# ============================================================

def get_current_stability_derivatives() -> dict:
    # Example values only.
    # These are assumed to already be the derivatives at the
    # current operating point from your external lookup/database.
    return {
        # Longitudinal
        "CX_u":      -0.030,
        "CX_alpha":  +0.180,
        "CX_q":      +0.000,

        "CZ_u":      -0.280,
        "CZ_alpha":  -5.200,
        "CZ_q":      -7.500,

        "Cm_u":      +0.012,
        "Cm_alpha":  -1.050,
        "Cm_q":      -16.000,

        # Lateral-directional
        "CY_beta":   -0.850,
        "CY_p":      +0.000,
        "CY_r":      +0.220,

        "Cl_beta":   -0.135,
        "Cl_p":      -0.620,
        "Cl_r":      +0.095,

        "Cn_beta":   +0.240,
        "Cn_p":      -0.028,
        "Cn_r":      -0.210,
    }


if __name__ == "__main__":
    # Example aircraft / trim
    ac = Aircraft(
        mass=1450.0,
        Ixx=980.0,
        Iyy=1520.0,
        Izz=2300.0,
        Ixz=45.0,
        S=16.8,
        b=10.9,
        cbar=1.62,
        rho=1.225,
    )

    trim = TrimPoint(
        U0=58.0,                  # m/s
        theta0=np.deg2rad(1.5),   # rad
    )

    derivs = get_current_stability_derivatives()

    A_lon = build_A_longitudinal(ac, trim, derivs)
    A_lat = build_A_lateral(ac, trim, derivs, include_Ixz_coupling=True)

    lon_modes = eig_summary(A_lon, ["u", "w", "q", "theta"])
    lat_modes = eig_summary(A_lat, ["v", "p", "r", "phi"])

    lon_labeled = classify_longitudinal_modes(lon_modes)
    lat_labeled = classify_lateral_modes(lat_modes)

    np.set_printoptions(precision=6, suppress=True)

    print("A_longitudinal =")
    print(A_lon)
    print()

    print("A_lateral =")
    print(A_lat)
    print()

    print_labeled_modes("LONGITUDINAL MODES", lon_labeled)
    print_labeled_modes("LATERAL MODES", lat_labeled)