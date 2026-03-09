# Climb
import numpy as np


def climb_rate(params):
    T_W = params["T_W"]
    rho_inf = params["rho_inf"]
    v_inf = params["v_inf"]
    W_S = params["W_S"]
    AR = params["AR"]
    cd0 = params["cd0"]
    return v_inf * (
        T_W
        - 0.5 * rho_inf * v_inf**2 * cd0 / W_S
        - W_S * (2 / (rho_inf * v_inf**2 * AR))
    )


params = {
    "T_W": ...,
    "rho_inf": ...,
    "v_inf": ...,
    "W_S": ...,
    "cd0": ...,
    "S": ...,
    "gamma": ...,
    "AR": ...,
    "L_D": ...,
}


def max_climb_angle(gamma):
    T_W = params["T_W"]
    L_D = params["L_D"]
    return np.sin(gamma) - T_W + 1 / L_D


# should i be solving for gamma(t)
# theroretically know everything except for v_C, could that be something we also want to assume
# are gamma and pitch attitude different from each other
# assuming non accelreaing clmb
# alpha=pitch-gamma
