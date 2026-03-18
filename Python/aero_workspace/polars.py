from aero_main import TakeoffCoeff, CruiseCoeff
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline


def takeoff_polars():
    CD_ind = [
        0.004911222748,
        0.07646259794,
        0.1189664504,
        0.1704539104,
        0.2306020724,
        0.3752174848,
        0.458684567,
        0.5488187362,
        0.6449712989,
        0.7464521674,
        0.852526822,
        0.9624413741,
        1.075416152,
        1.190660597,
        1.307381279,
        1.42479456,
        1.542124715,
        1.658625877,
    ]
    C_L = [
        0.33881,
        1.33686,
        1.66753,
        1.99602,
        2.32163,
        2.96144,
        3.2743,
        3.58159,
        3.88268,
        4.17698,
        4.46391,
        4.74295,
        5.0136,
        5.2754,
        5.52793,
        5.77082,
        6.00373,
        6.22638,
    ]

    C_D = CruiseCoeff.CD_DP + CD_ind
    cd_smooth = np.linspace(min(C_D), max(C_D), 500)
    spl = make_interp_spline(C_D, C_L, k=4)
    cl_smooth = spl(cd_smooth)

    # add angle of attack points
    y = [4.46391, 5.0136, 5.77082, 6.22638]
    x = [0.852526822, 1.075416152, 1.42479456, 1.658625877]
    labels = [r'$\alpha$=0$^{\circ}$', r'$\alpha$=5$^{\circ}$', r'$\alpha$=10$^{\circ}$', r'$\alpha$=15$^{\circ}$']

    plt.figure(figsize=(8,7))
    plt.rcParams["font.size"] = 14
    plt.title("Takeoff")
    plt.plot(cd_smooth, cl_smooth, color="blue")
    plt.xlabel("$C_D$")
    plt.ylabel("$C_L$")

    plt.scatter(x,y)
    for i, label in enumerate(labels):
        plt.text(x[i], y[i], labels[i])

    plt.show()

def cruise_polars():
    CD_tot = [
        0.01971088318,
        0.01672137502,
        0.01605249978,
        0.01770689687,
        0.02166429446,
        0.02788089001,
        0.03629081639,
        0.04680722209,
        0.05932263602,
        0.07371162078,
        0.08983174688,
    ]
    C_L = [
        -0.29451,
        -0.12985,
        0.03503,
        0.19974,
        0.36386,
        0.52697,
        0.68867,
        0.84857,
        1.00628,
        1.16143,
        1.31366,
    ]

    # Use a 2nd-degree polynomial to get the smooth U-shape
    cl_smooth = np.linspace(min(C_L), max(C_L), 500)
    coeffs = np.polyfit(C_L, CD_tot, 2)
    poly_func = np.poly1d(coeffs)
    cd_smooth = poly_func(cl_smooth)

    # add angle of attack points
    y = [-0.29451, 0.19974, 0.52697, 1.00628]
    x = [0.01971088318, 0.01770689687, 0.02788089001, 0.05932263602]
    labels = [r'$\alpha$=-5$^{\circ}$', r'$\alpha$=0$^{\circ}$', r'$\alpha$=5$^{\circ}$', r'$\alpha$=10$^{\circ}$']

    plt.figure(figsize=(8,7))
    plt.rcParams["font.size"] = 14
    plt.plot(cd_smooth, cl_smooth, color="blue")
    plt.title("Cruise")
    plt.xlabel("$C_D$")
    plt.ylabel("$C_L$")

    plt.scatter(x,y)
    for i, label in enumerate(labels):
        plt.text(x[i], y[i], labels[i])

    plt.show()


if __name__ == "__main__":
    cruise_polars()
    takeoff_polars()
