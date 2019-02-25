"""
Process xlsx file of perfusion data.
"""

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from biotree.utilities import model_1567, linear_fit, add_legend
from biotree.sensor import Sensor
plt.style.use("ggplot")


class Perfusion(Sensor):
    # path to the directory holding the perfusion data.
    path = ""

    def __getitem__(self, samples):
        """
        select samples and return data as a pandas dataframe
        """
        return(self.data[samples])

    def plot_sample(self, sample, base=270, with_model=True,
                    x_shift=0, y_shift=0,
                    color="blue", figsize=(8, 5), fit_from=30, fit_to=90):
        """
        Plot a perfusion curve with comparison to model curve and fitting to
        mercox perfusion stage.

        Parameters
        ----------
        sample : sample name such as "BT1880_HL".
        base : base measurement of pressure sensor such as 270.45.
        with_model : whether to plot the model curve for comparison.
        x_shift : shift the model curve along x-axis, for example, by 538.
        color : color of the curve.
        figsize : size of the plot
        fit_from,fit_to : fit range after mercox perfusion started.
        """

        if self.data is None:
            self.read(sample)
        else:
            if sample not in self.data.columns:
                self.read(sample)

        spl = self.data[sample].dropna()
        spl = (spl - base) / 30.7 * 43.6  # to mmHg
        plt.figure(figsize=figsize)
        plt.plot(spl.index, spl, c=color)

        # plot model
        if with_model:
            mdl = model_1567()
            plt.scatter(mdl.index + x_shift, mdl + y_shift, c="red", s=0.5)

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

    def plot_together(self, samples, x_shifts=None, colors=None, alphas=None,
                      figsize=(8, 5)):
        """
        Plot perfusion pressure curves together aligned.

        Parameters
        ----------
        samples : list of columns such as ["BT1880_HL", "BT1881_HL"].
        x_shifts : list of shifts along x-axis of each curve in samples.
            If None, curves will be aligned automatically at the beginning
            of Mercox perfusion.
        colors : list of colors for each curve. Use default color if None.
        alphas : list of alphas for each curve. Use 1 if None.

        Returns
        -------
        No return but make a matplotlib.pyplot plot.
        """
        # read sample, x_shift, and base
        dirname = os.path.dirname(__file__)
        filename = os.path.join(dirname, 'x_shift.csv')
        sxb = pd.read_csv(filename, index_col="sample")

        for spl in samples:
            if self.data is None:
                self.read(spl)
            else:
                if spl not in self.data.columns:
                    self.read(spl)
            if spl not in sxb.index:
                sxb.loc[spl, "x_shift"] = 0
                sxb.loc[spl, "base"] = 270
                print("Please use plot_sample() to get x_shift for " + spl)

        if alphas is None:
            alphas = [1] * len(samples)

        plt.figure(figsize=figsize)

        if x_shifts is None:
            if colors is None:
                for spl, alpha in zip(samples, alphas):
                    # 1.42 = 43.7 / 30.7
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
