"""
Process csv file of perfusion data.
"""

import os
import glob
from colorama import Fore
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from biotree.utilities import model_1567, linear_fit
from biotree.utilities import exponential_func
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

    def get_fullname(self, sample):
        """
        Get the full sample name. For example, if the argument sample is
        "BT2157" and the full sample is "BT2157_HL", it will return the
        full name.
        """
        file_name = self.path + sample
        file_name = glob.glob(file_name + "*")
        if len(file_name) > 1:
            print(Fore.RED, sample + " has multiple matches. "
                  + "Please select one from below.\n")
            match = [fn.split("/")[-1][:-4] for fn in file_name]
            match.sort()
            print(match)
            return(None)

        file_name = file_name[0].split("/")[-1]
        fullname = file_name[:-4]

        return(fullname)

    def read_shift_base(self):
        """
        Read the file containing x_shift and base of each perfusion.

        Return
        ------
        A dataframe having two columns "x_shift" and "base". Its rows are
        indexed with sample name such "BT2240_AL".
        """
        dirname = os.path.dirname(__file__)
        filename = os.path.join(dirname, 'x_shift.csv')
        sxb = pd.read_csv(filename, index_col="sample")
        return sxb

    def get_highflashing(self, sample, smooth=True):
        """
        Get perfusion data at high flashing stage
        """
        sample = self.get_fullname(sample)
        self.read(sample)
        sxb = self.read_shift_base()
        base = sxb.loc[sample, "base"]
        dat = self.data[sample].dropna()
        dat = (dat - base) * 1.42  # to mmHg

        hf_first = dat[dat > 50].index[0]
        if smooth:
            dat = dat.rolling(40).mean()
        hf = dat[hf_first:(hf_first + 250)]
        hf.index = hf.index - hf_first
        hf.iloc[-1] = 0

        return hf

    def get_mercox(self, sample, smooth=True):
        sample = self.get_fullname(sample)
        self.read(sample)
        sxb = self.read_shift_base()
        base = sxb.loc[sample, "base"]
        x_shift = sxb.loc[sample, "x_shift"]
        dat = self.data[sample].dropna()
        dat = (dat - base) * 1.42  # to mmHg
        if smooth:
            dat = dat.rolling(150).mean()
        mercox = dat[x_shift:]
        mercox.index = mercox.index - x_shift + 1050
        mercox.iloc[0] = 0

        return mercox

    def get_pfa(self, sample, cut_index=0, smooth=True):
        """
        ignore section whose index smaller than cut_index
        """
        sample = self.get_fullname(sample)
        self.read(sample)
        sxb = self.read_shift_base()
        base = sxb.loc[sample, "base"]
        dat = self.data[sample].dropna()
        dat = (dat - base) * 1.42  # to mmHg
        dat = dat[cut_index:]

        dat_diff = dat[dat.diff(30).abs() > 20]
        idx_diff = dat_diff.index.to_series().diff()
        idx_diff_max = idx_diff.max()
        idx_diff_max_idx = idx_diff.idxmax()

        pfa_first = idx_diff_max_idx - idx_diff_max
        if smooth:
            dat = dat.rolling(150).mean()
        pfa_dat = dat.loc[pfa_first:(pfa_first + 700)]

        pfa_dat.index = pfa_dat.index - pfa_first + 300
        pfa_dat.iloc[0] = 0
        pfa_dat.iloc[-1] = 0
        return pfa_dat

    def plot_hf_pfa_mercox(self, sample, cut_index=0, smooth=True):
        hf = self.get_highflashing(sample, smooth)
        pfa = self.get_pfa(sample, cut_index, smooth)
        mrx = self.get_mercox(sample, smooth)
        dat = pd.concat([hf, pfa, mrx])
        plt.plot(dat.index, dat, label=sample)
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (mmHg)")
        plt.legend()

    def fit_highflashing(self, sample, t0=0, tn=50, p0=(10, 0.1, 80)):
        """
        Fit the exponential decay of the high flashing stage.

        Parameters
        ----------
        sample: BT nnumber such as "BT1234".
        base: base pressure.
        t0: time where the peak is placed.
        tn: time in second after peak used for fitting.
        p0: initial guess of fitting paraters of a, b, c in exponential decay.
        """
        sample = self.get_fullname(sample)

        # read sample, x_shift, and base
        # dirname = os.path.dirname(__file__)
        # filename = os.path.join(dirname, 'x_shift.csv')
        # sxb = pd.read_csv(filename, index_col="sample")
        sxb = self.read_shift_base()
        base = sxb.loc[sample, "base"]

        self.read(sample)
        dat = self.data[sample].dropna()
        dat = (dat - base) * 1.42  # to mmHg

        x0 = np.asarray(dat.index)[0:2500].squeeze()
        y0 = np.asarray(dat)[0:2500].squeeze()

        # set the peak at 1 second for all perfusion curves
        n_peak = np.argmax(y0)
        x0 = x0 - n_peak / 10 + t0

        x = x0[n_peak:(n_peak + 10 * tn)]
        y = y0[n_peak:(n_peak + 10 * tn)]

        abc, _ = curve_fit(exponential_func, x, y, p0=p0)
        xx = np.linspace(0, tn, 1000)
        yy = exponential_func(xx, *abc)

        # force all plot in the same color with c=p1[0].get_color()
        p1 = plt.plot(x0, y0, linewidth=1, label="_nolegend_", alpha=0.5)
        plt.plot(xx, yy, c=p1[0].get_color(), linewidth=1,
                 label="{}: a = {:.1f}, b = {:.3f}, c = {:.1f}".
                       format(sample, abc[0], abc[1], abc[2]))
        plt.scatter(x, y, s=10, alpha=0.5)
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (mmHg)")
        plt.legend()

    def plot_sample(self, sample, base=270, n_smooth=1, with_model="BT1567",
                    x_shift=0, y_shift=0, diff_method="mean",
                    color="blue", figsize=(8, 5), fit_from=30, fit_to=300,
                    xlim=None, ylim=None):
        """
        Plot a perfusion curve with comparison to model curve and fitting to
        mercox perfusion stage.

        Parameters
        ----------
        sample : sample name such as "BT1880_HL" or "BT1880".
        base : base measurement of pressure sensor such as 270.45.
        n_smooth : number of point used to smooth the plot.
        with_model : whether to plot the model curve for comparison.
        x_shift : shift the model curve along x-axis, for example, by 538.
        color : color of the curve.
        figsize : size of the plot
        fit_from,fit_to : fit range after mercox perfusion started.
        xlim,ylim : plot range of x and y axis.
        """
        # get full name
        sample = self.get_fullname(sample)

        # save x_shift for future use in plot
        dirname = os.path.dirname(__file__)
        filename = os.path.join(dirname, 'x_shift.csv')
        sample_xshift = pd.read_csv(filename, index_col="sample")
        sample_xshift.loc[sample] = [x_shift, base]
        sample_xshift.to_csv(filename)

        if sample is None:
            return(None)

        if self.data is None:
            self.read(sample)
        else:
            if sample not in self.data.columns:
                self.read(sample)

        spl = self.data[sample].dropna()
        spl = (spl - base) * 1.42  # to mmHg
        spl = spl.rolling(n_smooth).mean()
        plt.figure(figsize=figsize)
        # convert to np.array which works for all matplotlib version
        plt.plot(np.asarray(spl.index), np.asarray(spl), c=color)

        # plot model
        if with_model == "BT1567":
            mdl = model_1567()
            plt.scatter(np.asarray(mdl.index) + x_shift,
                        np.asarray(mdl) + y_shift, c="red", s=0.5)
        elif with_model == "BT2239":
            print("to be completed")
            pass

        # plot fit line to mercox perfusion
        mercox = spl.loc[(x_shift + fit_from + 120):(x_shift + fit_to + 120)]
        fit = linear_fit(mercox)
        x = np.arange(x_shift + 100, x_shift + 580, 1)
        y = fit[0] + x * fit[1]
        plt.plot(x,  y, c="black", lw=1)

        # difference between the flat stage of sample and model
        flat_median = spl.loc[(x_shift + 90):(x_shift + 106)].median()
        flat_mean = spl.loc[(x_shift + 90):(x_shift + 106)].mean()

        if diff_method == "median":
            y_diff = round(flat_median - 41.7, 2)
        elif diff_method == "mean":
            y_diff = round(flat_mean - 41.7, 2)

        # annotate position relative to axis
        plt.annotate(sample, xy=(0.05, 0.9), xycoords="axes fraction",
                     fontsize=20)
        plt.annotate("x_shift: " + str(x_shift) + " s", xy=(0.05, 0.8),
                     xycoords="axes fraction", fontsize=16)
        plt.annotate("y_diff: " + str(y_diff) + " mmHg", xy=(0.05, 0.7),
                     xycoords="axes fraction", fontsize=16)
        plt.annotate("slope: " + str(round(fit[1], 2)) + " mmHg/s",
                     xy=(0.05, 0.6), xycoords="axes fraction", fontsize=16,
                     color="black")

        plt.xlabel("Time (sec)", fontsize=16)
        plt.ylabel("Pressure (mmHg)", fontsize=16)

        if ylim is None:
            plt.ylim(-20, 1100)
        else:
            plt.ylim(ylim[0], ylim[1])
        if xlim is None:
            plt.xlim(-20, spl.index.max() + 20)
        else:
            plt.xlim(xlim[0], xlim[1])

    def plot_together(self, samples, align=None, n_smooth=1, colors=None,
                      alphas=None, figsize=(8, 5)):
        """
        Plot perfusion pressure curves together aligned.

        Parameters
        ----------
        samples : list of columns such as ["BT1880_HL", "BT1881_HL"].
        align : align curves by shifting along x-axis of each curve in samples.
            If None, aligned at time 0. If "hf", aligned at peak of high
            flashing. If "mrx", aligned at the rising of Merocox perfusion.
            If list, shift each curve according to the elements.
        n_smooth : number of point used to smooth the plot.
        colors : list of colors for each curve. Use default color if None.
        alphas : list of alphas for each curve. Use 1 if None.

        Returns
        -------
        No return but make a matplotlib.pyplot plot.
        """
        samples = [self.get_fullname(sample) for sample in samples]

        # read sample, x_shift, and base
        # dirname = os.path.dirname(__file__)
        # filename = os.path.join(dirname, 'x_shift.csv')
        # sxb = pd.read_csv(filename, index_col="sample")
        sxb = self.read_shift_base()

        if None in samples:
            return(None)

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

        df = self.data.rolling(n_smooth).mean()
        x = np.asarray(df.index)

        plt.figure(figsize=figsize)

        if align is None:  # start from time zero
            if colors is None:
                for spl, alpha in zip(samples, alphas):
                    plt.plot(x,
                             (np.asarray(df[spl])
                              - sxb.loc[spl, "base"]) * 1.42,
                             alpha=alpha, label=spl)
            else:
                for spl, color, alpha in zip(samples, colors, alphas):
                    plt.plot(x,
                             (np.asarray(df[spl])
                              - sxb.loc[spl, "base"]) * 1.42,
                             c=color, alpha=alpha, label=spl)
        elif align == "mrx":  # mercox
            if colors is None:
                for spl, alpha in zip(samples, alphas):
                    # 1.42 = 43.7 / 30.7
                    plt.plot(x - sxb.loc[spl, "x_shift"],
                             (np.asarray(df[spl])
                             - sxb.loc[spl, "base"]) * 1.42,
                             alpha=alpha, label=spl)
            else:
                for spl, color, alpha in zip(samples, colors, alphas):
                    plt.plot(x - sxb.loc[spl, "x_shift"],
                             (np.asarray(df[spl])
                              - sxb.loc[spl, "base"]) * 1.42,
                             c=color, alpha=alpha, label=spl)
        elif align == "hf":  # for high flashing
            if colors is None:
                for spl, alpha in zip(samples, alphas):
                    x0 = x - np.nanargmax(np.asarray(df[spl])[0:1000]) / 10
                    print(np.argmax(np.asarray(df[spl])[0:1000]) / 10)
                    plt.plot(x0,
                             (np.asarray(df[spl])
                              - sxb.loc[spl, "base"]) * 1.42,
                             alpha=alpha, label=spl)
            else:
                for spl, color, alpha in zip(samples, colors, alphas):
                    x0 = x - np.nanargmax(np.asarray(df[spl])[0:1000]) / 10
                    plt.plot(x0,
                             (np.asarray(df[spl])
                              - sxb.loc[spl, "base"]) * 1.42,
                             c=color, alpha=alpha, label=spl)
        elif len(align) > 1:  # for pfa
            if colors is None:
                for spl, x_shift, alpha in zip(samples, align, alphas):
                    plt.plot(x - x_shift,
                             (np.asarray(df[spl])
                              - sxb.loc[spl, "base"]) * 1.42,
                             alpha=alpha, label=spl)
            else:
                for spl, x_shift, color, alpha in zip(samples, align,
                                                      colors, alphas):
                    plt.plot(x - x_shift,
                             (np.asarray(df[spl])
                              - sxb.loc[spl, "base"]) * 1.42,
                             c=color, alpha=alpha, label=spl)

        plt.xlabel("Time (sec)", fontsize=16)
        plt.ylabel("Pressure (mmHg)", fontsize=16)
        plt.ylim(-20, 1100)
        plt.legend()
