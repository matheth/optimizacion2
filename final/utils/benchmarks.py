import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the functions


def ackley(x, y, a=20.0, b=0.2, c=2 * np.pi):
    """
    minimiser = (0,0)
    """
    d = len(x) if isinstance(x, list) or isinstance(x, np.ndarray) else 1
    sum_1 = np.sqrt(1 / d * np.sum(np.array(x) ** 2))
    sum_2 = np.sum(np.cos(c * np.array(x)))
    return -a * np.exp(-b * sum_1) - np.exp(1 / d * sum_2) + a + np.exp(1)


def bukin(x, y):
    """
    minimiser: x = (-10,1)
    """
    return np.sqrt(np.abs(y - 0.01 * x**2)) + 0.01 * np.abs(x + 10)


def cross_in_tray(x, y):
    """
    minimisers: x = (1.3491, -1.3491), (1.3491, 1.3491), (-1.3491, 1.3491), (-1.3491, -1.3491)
    """
    return -1e-4 * (
        (
            np.abs(
                np.sin(x)
                * np.sin(y)
                * np.exp(np.abs(100 - np.sqrt(x**2 + y**2) / np.pi))
            )
            + 1
        )
        ** 0.1
    )


def drop_wave(x, y):
    """
    minimiser: x = (0,0)
    """
    return -(1 + np.cos(12 * np.sqrt(x**2 + y**2))) / (0.5 * (x**2 + y**2) + 2)


def eggholder(x, y):
    """
    minimiser: x = (512, 404.2319)
    """
    return -(x + 47) * np.sin(np.sqrt(np.abs(y + x / 2 + 47))) - x * np.sin(
        np.sqrt(np.abs(x - (y + 47)))
    )


def griewank(x, y):
    """
    minimiser : x = 0
    """
    return (
        np.sum(np.array(x) ** 2) / 4000
        - np.prod(np.cos(np.array(x) / np.sqrt(np.arange(1, len(x) + 1))))
        + 1
    )


def holder_table(x, y):
    """
    minimisers : x = (8.05502, 9.66459), (8.05502, -9.66459), (-8.05502, 9.66459), (-8.05502, -9.66459)
    """
    return -np.abs(
        np.sin(x) * np.cos(y) * np.exp(np.abs(1 - np.sqrt(x**2 + y**2) / np.pi))
    )


def langermann(x, y, c=[1, 2, 5, 2, 3], a=np.array([[3, 5, 2, 1, 7], [5, 2, 1, 4, 9]])):
    """
    minimiser: x= (2.00299219, 1.006096)
    """
    m, d = a.shape
    return np.sum(
        c[i]
        * np.exp(-1 / np.pi * np.sum((x[j] - a[i, j]) ** 2 for j in range(d)))
        * np.cos(np.pi * np.sum((x[j] - a[i, j]) ** 2 for j in range(d)))
        for i in range(m)
    )


def levy(x, y):
    """
    minimiser: x = (1,...,1)
    """

    def w(x):
        return 1 + (x - 1) / 4 if x > 1 else 1

    if len(x) > 2:
        s = np.sum(
            (w(x[i]) - 1) ** 2 * (1 + 10 * np.sin(np.pi * w(x[i]) + 1) ** 2)
            for i in range(1, len(x) - 1)
        )
    else:
        s = 0
    return (
        np.sin(np.pi * w(x[0])) ** 2
        + s
        + (w(y) - 1) ** 2 * (1 + np.sin(2 * np.pi * w(y)) ** 2)
    )

# -------------------------------------BOWL-SHAPED-------------------------------------------------------
# -------------------------------------------------------------------------------------------------------

"""
    bochavesky functions
    minimiser: (x,y) = (0,0)
    n is the number of the bochavesky function
"""


def bochavesky(x, y, n=1):
    if n == 1:
        return (
            x**2
            + 2 * (y**2)
            - 0.3 * np.cos(3 * np.pi * x)
            - 0.4 * np.cos(4 * np.pi * y)
            + 0.7
        )
    elif n == 2:
        return (
            x**2
            + 2 * (y**2)
            - 0.3 * np.cos(3 * np.pi * x) * np.cos(4 * np.pi * y)
            + 0.3
        )
    elif n == 3:
        return x**2 + 2 * (y**2) - 0.3 * np.cos(3 * np.pi * x + 4 * np.pi * y) + 0.3


"""
    perm_zero_d_beta function
    minimiser: (1, 1/2, ..., 1/d)
"""


def perm_zero_d_β(x, β=1):
    return np.sum(
        np.sum((j + β) * (x[j] ** i - 1 / (j**i)) for j in range(1, len(x) + 1)) ** 2
        for i in range(1, len(x) + 1)
    )


"""
    rot_hyper_ellipsoid function
    minimiser: x = 0
"""


def rot_hyper_ellipsoid(x):
    return np.sum(np.sum(x[:i] ** 2 for i in range(1, len(x) + 1)))


"""
    sphere
    minimiser: x = 0
"""


def sphere(x, y):
    return (
        np.sum(np.array(x) ** 2)
        if isinstance(x, list) or isinstance(x, np.ndarray)
        else x**2 + y**2
    )


"""
    su_powers function
    minimiser: x = 0
"""


def sum_powers(x, y):
    return x**2 + np.abs(y) ** 3


"""
    sum_squares
    minimiser: x = 0
"""


def sum_squares(x, y):
    return x**2 + 2 * y**2


"""
    trid function
    minimiser: xᵢ   = i(d+1-i) ∀ i=1,2,...,d
"""


def trid(x, y):
    return (
        np.sum((np.array(x) - 1) ** 2) - np.sum(np.array(x[1:]) * np.array(x[:-1]))
        if isinstance(x, list) or isinstance(x, np.ndarray)
        else (x - 1) ** 2 - x * y
    )


# -------------------------------------PLATE-SHAPED------------------------------------------------------
# -------------------------------------------------------------------------------------------------------

"""
    booth function
    minimiser: x= (1,3)
"""


def booth(x, y):
    return (x + 2 * y - 7) ** 2 + (2 * x + y - 5) ** 2


""" 
    matyas function
    minimiser: x= (0,0)
"""


def matyas(x, y):
    return 0.26 * (x**2 + y**2) - 0.48 * x * y


"""
    mccormick function
    minimiser: x= (-0.54719, -1.54719)
"""


def mccormick(x, y):
    return np.sin(x + y) + (x - y) ** 2 - 1.5 * x + 2.5 * y + 1


"""
    power_sum function
"""


def power_sum(x, y, b=[8, 18]):
    return np.sum(
        (np.sum(x[j] ** (i + 1) for j in range(len(x))) - b[i]) ** 2
        for i in range(len(x))
    )


"""
    zakharov function
    minimiser: x = 0
"""


def zakharov(x, y):
    return (
        np.sum(x**2)
        + np.sum(0.5 * (np.arange(1, len(x) + 1)) * x) ** 2
        + np.sum(0.5 * (np.arange(1, len(x) + 1)) * x) ** 4
    )


# -------------------------------------VALLEY-SHAPED-----------------------------------------------------
# -------------------------------------------------------------------------------------------------------

"""
    three_hump_cael
    global minimiser: x= (0,0)
    Has three local minima
"""


def three_hump_camel(x, y):
    return 2 * x**2 - 1.05 * x**4 + (x**6) / 6 + x * y + y**2


"""
    six_hump_camel
    global minimisers : x = (0.0898, -0.7126)  and  x = (-0.0898, 0.7126)
    Has six local minima
"""


def six_hump_camel(x, y):
    return (4 - 2.1 * x**2 + (x**4) / 3) * (x**2) + x * y + (-4 + 4 * y**2) * (y**2)


"""
    dixon_price function
    minimiser (for d=2): x = (1,1/sqrt(2))
"""


def dixon_price(x, y):
    return (
        sum(i * (2 * x[i] ** 2 - x[i - 1]) ** 2 for i in range(1, len(x)))
        + (x[0] - 1) ** 2
    )


"""
    rosenbrock function
    minimiser : x = (1,..., 1)
"""


def rosenbrock(x, y):
    return sum(
        100 * (x[i + 1] - x[i] ** 2) ** 2 + (x[i] - 1) ** 2 for i in range(len(x) - 1)
    )


# -------------------------------------STEEP RIDGES/DROPS-----------------------------------------
# ------------------------------------------------------------------------------------------------


def de_jong_5(x):
    r1 = np.repeat([-32, -16, 0, 16, 32], 5)
    r2 = np.tile([-32, -16, 0, 16, 32], 5)
    a = np.vstack((r1, r2)).T
    return 1 / (
        0.002
        + np.sum(
            1 / (i + (x[0] - a[i, 0]) ** 6 + (x[1] - a[i, 1]) ** 6) for i in range(25)
        )
    )


"""
    easom function
    minimiser : x = (pi,pi)
"""


def easom(x, y):
    return -np.cos(x) * np.cos(y) * np.exp(-((x - np.pi) ** 2) - (y - np.pi) ** 2)


"""
    michalewicz function
    m defines steepness of valleys and ridges
    Has d! local minima
"""


def michalewicz(x, y, m=10):
    return -np.sum(
        np.sin(x[i]) * (np.sin((i + 1) * x[i] ** 2 / np.pi)) ** (2 * m)
        for i in range(len(x))
    )


# -------------------------------------OTHERS------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------

"""
    beale function
    minimiser : x = (3,0.5)
"""


def beale(x, y):
    return (
        (1.5 - x + x * y) ** 2
        + (2.25 - x + x * y**2) ** 2
        + (2.625 - x + x * y**3) ** 2
    )


"""
    branin function
    minimisers : x = (-pi,12.275) , (pi, 2.275), (9.42478, 2.475)
"""
pi = math.pi


def branin(x, y, a=1, b=5.1 / (4 * pi**2), c=5 / pi, r=6, s=10, t=1 / (8 * pi)):
    return a * (y - b * x + c * x - r) ** 2 + s * (1 - t) * np.cos(x) + s


"""
    colville function
    minimiser : x = (1,1,1,1)
"""


def colville(x):
    return (
        100 * (x[0] ** 2 - x[1]) ** 2
        + (x[0] - 1) ** 2
        + (x[1] - 1) ** 2
        + 90 * (x[2] ** 2 - x[3]) ** 2
        + 10.1 * ((x[1] - 1) ** 2 + (x[3] - 1) ** 2)
        + 19.8 * (x[1] - 1) * (x[3] - 1)
    )


"""
    goldstein_price
    minimiser : x = (0,-1)
"""


def goldstein_price(x, y):
    return (
        1
        + (x + y + 1) ** 2 * (19 - 14 * x + 3 * x**2 - 14 * y + 6 * x * y + 3 * y**2)
        + 30
        + (2 * x - 3 * y) ** 2
        * (18 - 32 * x + 12 * y**2 + 48 * y - 36 * x * y + 27 * y**2)
    )


"""
    perm_d_beta
    minimiser : x = (1,2,...,d)
"""


def perm_d_beta(x, β=1):
    return np.sum(
        np.sum((j**i + β) * ((x[j - 1] / j) ** i - 1) for j in range(1, len(x) + 1))
        ** 2
        for i in range(1, len(x) + 1)
    )


"""
    minimiser : x = (-2.903534,...,-2.903534)
"""


def styblinski_tang(x):
    return 0.5 * np.sum(x**4 - 16 * x**2 + 5 * x)


"""
    parsopoulos
    ∞ global minimisers: x = (k*(pi/2), l*pi) for k=1,2,3,... and l= 0, 1, 2, ...
"""


def parsopoulos(x, y):
    return np.cos(x) ** 2 + np.sin(y) ** 2


limits_dict = {
    "ackley": (-5, 5, -5, 5),
    "bukin": (-15, -5, -3, 3),
    "cross_in_tray": (-5, 5, -5, 5),
    "drop_wave": (-5.2, 5.2, -5.2, 5.2),
    "eggholder": (-520, 520, -520, 520),
    "griewank": (-5, 5, -5, 5),
    "holder_table": (-10, 10, -10, 10),
    "langermann": (0, 10, 0, 10),
    "levy": (-10, 10, -10, 10),
    "bochavesky": (-100, 100, -100, 100),
    "perm_zero_d_beta": (-2, 2, -2, 2),
    "rot_hyper_ellipsoid": (-100, 100, -100, 100),
    "sphere": (-100, 100, -100, 100),
    "sum_powers": (-1, 1, -1, 1),
    "sum_squares": (-10, 10, -10, 10),
    "trid": (-4, 4, -4, 4),
    "booth": (-10, 10, -10, 10),
    "matyas": (-10, 10, -10, 10),
    "mccormick": (-1.5, 4, -3, 4),
    "power_sum": (-2, 2, -2, 2),
    "zakharov": (-5, 10, -5, 10),
    "three_hump_camel": (-5, 5, -5, 5),
    "six_hump_camel": (-1.5, 1.5, -1.5, 1.5),
    "dixon_price": (-10, 10, -10, 10),
    "rosenbrock": (-2.048, 2.048, -2.048, 2.048),
    "de_jong_5": (-50, 50, -50, 50),
    "easom": (-20, 20, -20, 20),
    "michalewicz": (0, 4, 0, 4),
    "beale": (-4.5, 4.5, -4.5, 4.5),
    "branin": (-5, 15, -5, 15),
    "goldstein_price": (-2, 2, -2, 2),
    "perm_d_beta": (-3, 3, -3, 3),
    "styblinski_tang": (-5, 5, -5, 5),
    "parsopoulos": (-5, 5, -5, 5),
}

# Define the plotting function
def plot_fun(f, points, limits, params):
    if limits is None:
        limits = limits_dict[f.__name__]

    x = np.linspace(limits[0], limits[1], 100)
    y = np.linspace(limits[2], limits[3], 100)
    X, Y = np.meshgrid(x, y)
    Z = f(X, Y)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z, cmap="viridis")

    if points is not None:
        ax.scatter(points[:, 0], points[:, 1], points[:, 2], color="red")

    plt.show()
