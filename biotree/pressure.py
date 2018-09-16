"""
Process xlsx file of perfusion data.
"""

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from biotree.utilities import model_1567, linear_fit, add_legend
plt.style.use("ggplot")


class Pressure:

    def __init__(self, description="perfusion data"):
        self.description = description
        self.path = ""

    def read(self, xlsx_file="perfusion_pressure_sensor.xlsx"):
        xlsx_file = self.path + xlsx_file
        data = pd.read_excel(xlsx_file)
        data.index = data.time
        data = data.drop("time", 1)
        self.data = data
        self.columns = data.columns
        self.index = data.index

    def to_pressure(self, base=270):
        # use mmHg
        self.data = (self.data - base) / 30.7 * 43.6
        return(self)

    def to_force(self, base=270):
        # so we can compare force measured with force gauge
        self.data = (self.data - base) / 30.7
        return(self)

    def plot_sample(self, sample, base=270, with_model=True, x_shift=0,
                    color="blue", figsize=(8, 5), fit_from=30, fit_to=90):
        # plot sample
        spl = self.data[sample].dropna()
        spl = (spl - base) / 30.7 * 43.6  # to mmHg
        plt.figure(figsize=figsize)
        plt.plot(spl.index, spl, c=color)

        # plot model
        if with_model:
            mdl = model_1567()
            plt.scatter(mdl.index + x_shift, mdl, c="red", s=0.5)

        # plot fit line to mercox perfusion
        mercox = spl.loc[(x_shift + fit_from + 120):(x_shift + fit_to + 120)]
        fit = linear_fit(mercox)
        x = np.arange(x_shift + 100, x_shift + 580, 1)
        y = fit[0] + x * fit[1]
        plt.plot(x,  y, c="black", lw=1)

        # difference between the flat stage of sample and model
        flat_median = spl.loc[(x_shift + 90):(x_shift + 106)].median()
        y_diff = round(flat_median - 41.7, 2)

        plt.annotate(sample, xy=(10, 1110), fontsize=20)
        plt.annotate("x_shift: " + str(x_shift) + " s", xy=(10, 990),
                     fontsize=16)
        plt.annotate("y_diff: " + str(y_diff) + " mmHg", xy=(10, 870),
                     fontsize=16)
        plt.annotate("slope: " + str(round(fit[1], 2)) + " mmHg/s",
                     xy=(10, 750), fontsize=16, color="black")

        plt.xlabel("Time (sec)", fontsize=16)
        plt.ylabel("Pressure (mmHg)", fontsize=16)
        plt.ylim(-100, 1200)
        plt.xlim(-20, spl.index.max() + 20)

        # save x_shift for future use in plot
        dirname = os.path.dirname(__file__)
        filename = os.path.join(dirname, 'x_shift.csv')
        sample_xshift = pd.read_csv(filename, index_col="sample")
        sample_xshift.loc[sample] = [x_shift, base]
        sample_xshift.to_csv(filename)

    def select(self, cols):
        return(self.data[cols])

    def compare_to_1567(self, column):
        pass

    def plot_together(self, samples, x_shifts=None, colors=None, alphas=None,
                      figsize=(8, 5)):
        """
        sample_xshift_color : a list of tuples like
            [("BT1902_HL", 355, "red"), ("BT1903_HL", 467, "blue")]
        """
        # read sample, x_shift, and base
        dirname = os.path.dirname(__file__)
        filename = os.path.join(dirname, 'x_shift.csv')
        sxb = pd.read_csv(filename, index_col="sample")

        if alphas is None:
            alphas = [1] * len(samples)

        plt.figure(figsize=figsize)

        if x_shifts is None:
            if colors is None:
                for spl, alpha in zip(samples, alphas):
                    plt.plot(self.index - sxb.loc[spl, "x_shift"],
                             (self.data[spl] - sxb.loc[spl, "base"]) * 1.42,
                             alpha=alpha)
            else:
                for spl, color, alpha in zip(samples, colors, alphas):
                    plt.plot(self.index - sxb.loc[spl, "x_shift"],
                             (self.data[spl] - sxb.loc[spl, "base"]) * 1.42,
                             c=color, alpha=alpha)
        else:
            if colors is None:
                for spl, x_shift, alpha in zip(samples, x_shifts, alphas):
                    plt.plot(self.index - x_shift,
                             (self.data[spl] - sxb.loc[spl, "base"]) * 1.42,
                             alpha=alpha)
            else:
                for spl, x_shift, color, alpha in zip(samples, x_shifts,
                                                      colors, alphas):
                    plt.plot(self.index - x_shift,
                             (self.data[spl] - sxb.loc[spl, "base"]) * 1.42,
                             c=color, alpha=alpha)

        plt.xlabel("Time (sec)", fontsize=16)
        plt.ylabel("Pressure (mmHg)", fontsize=16)
        plt.ylim(-100, 1200)
        add_legend(loc="upper left")
