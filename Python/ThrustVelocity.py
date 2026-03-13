import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# THRUST INTERPOLATOR (first used in takeoff model) (0–22 m/s)
# ============================================================

V_data = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130])
T_data = 8 * 10 * np.array([709.4, 414.2, 303.2, 257.9, 225.0, 200.0, 182.0, 170.0, 160.0, 152.0, 145.0, 139.0, 135.0, 131.0]) # Thrust per fan [N]
degree = 5          # change this to 1,2,3,4,... to test fits
coeffs = np.polyfit(V_data, T_data, degree)
T_poly = np.poly1d(coeffs)


class ThrustVelocity:
    def __init__(self, V_data=V_data, T_data=T_data, degree=degree):
        self.coeffs = np.polyfit(V_data, T_data, degree)
        self.poly = np.poly1d(self.coeffs)

    def get_T(self, v):
        """Interpolated thrust per fan [N]"""
        return self.poly(v)

def T_fan_interp(v):
        """Interpolated thrust per fan [N]"""
        return T_poly(v)

vel = np.linspace(0,130,500)
plt.figure(figsize=(8,6))
plt.scatter(V_data, T_data, label="Original Data")
plt.plot(vel, T_fan_interp(vel), label=f"Polynomial Fit (deg={degree})")
plt.xlabel("Velocity [m/s]")
plt.ylabel("Thrust per fan [N]")
plt.title("Fan Thrust Interpolation")
plt.legend()
plt.grid(True)
plt.show()
