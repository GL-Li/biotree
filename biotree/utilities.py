"""
define utility function for use in the package
"""

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats


def read_1567(unit="mmHg"):
    """read BT1567_SF perfusion data into a pandas series"""
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, 'BT1567_SF.xlsx')
    df = pd.read_excel(filename)
    df.index = df.time

    if unit == "mmHg":
        df["pressure"] = df["BT1567_SF"] / (3.1416 * 0.0074**2) / 133.3
        return(df.pressure)

    # return pressure as a series
    return(df.BT1567_SF)


def model_1567():
    # match to BT1567 mannually
    BT1567 = read_1567()
    sect_0 = BT1567[85:110]
    sect_1 = BT1567[116.8:130.4]
    sect_2 = BT1567[136.0:360.0]

    fit_1 = linear_fit(sect_1)
    fit_2 = linear_fit(sect_2)

    # section 0 fit, use round() to avoid uncertain floating numbers for index
    x0 = np.arange(40, 136, 0.1).round(1)
    y0 = (sect_0.mean() * np.ones(1360 - 400)).round(1)

    # section 1 fit
    x1 = np.arange(100, 140, 0.1).round(1)
    y1 = (fit_1[0] + fit_1[1] * x1).round(1)

    # section 2 fit
    x2 = np.arange(112, 580, 0.1).round(1)
    y2 = (fit_2[0] + fit_2[1] * x2).round(1)

    x = np.concatenate((x0, x1, x2))
    y = np.concatenate((y0, y1, y2))

    return(pd.Series(data=y, index=x))


def linear_fit(series):
    """
    Get slope and intercept of linear fit to a series

    Parameters
    ----------
    series: numerical series with numeric index

    Returns
    -------
    a list of [intercept, slope]
    """

    slope, intercept, r_value, p_value, stderr = stats.linregress(
        series.index.values, series
    )

    return((intercept, slope))


def add_legend(loc=(1, 0), linewidth=2):
    """
    Add legend to a plot with controlled position and linewidth
    """
    leg = plt.legend(loc=loc)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(linewidth)


def exponential_func(x, a, b, c):
    return a * np.exp(-b * x) + c
