import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# THRUST INTERPOLATOR (first used in takeoff model) (0–22 m/s)
# ============================================================
N_fans = 14
V_data = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 120.0])
T_data = N_fans * np.array([3503.0, 3231.0, 2989.0, 2777.0, 2582.0, 2411.0, 2256.0, 2115.0, 1990.0, 1877.0, 1775.0, 1686.0, 1605.0, 1533.0]) # Thrust per fan [N]
degree = 5          # change this to 1,2,3,4,... to test fits
coeffs = np.polyfit(V_data, T_data, degree)
T_poly = np.poly1d(coeffs)


# # Values assuming takeoff thrust of about 41kN
# V_data = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
# T_data = [5936.0, 5511.0, 5135.0, 4791.0, 4485.0, 4209.0, 3962.0, 3744.0, 3549.0, 3372.0, 3209.0]

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
# plt.figure(figsize=(8,6))
# plt.scatter(V_data, T_data/1000, label="Original Data")
# plt.plot(vel, T_fan_interp(vel)/1000, label=f"Polynomial Fit (deg={degree})")
# plt.xlabel("Velocity [m/s]")
# plt.ylabel("Thrust per fan [kN]")
# plt.title("Fan Thrust Interpolation")
# plt.legend()
# plt.grid(True)
# plt.show()
