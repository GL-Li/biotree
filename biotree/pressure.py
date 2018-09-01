"""
Process xlsx file of perfusion data.
"""

import pandas as pd
from matplotlib import pyplot as plt
from biotree.utilities import model_1567
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
                    color="blue", figsize=(9, 6)):
        spl = self.data[sample]

        plt.figure(figsize=figsize)
        plt.plot(spl.index, (spl - base) / 30.7 * 43.6, c=color)

        if with_model:
            mdl = model_1567()
            plt.scatter(mdl.index + x_shift, mdl, c="red", s=1)

    def select(self, cols):
        return(self.data[cols])

    def compare_to_1567(self, column):
        pass

    def plot_together(self, sample_xshift_color):
        """
        sample_xshift_color : a list of tuples like
            [("BT1902_HL", 355, "red"), ("BT1903_HL", 467, "blue")]
        """
        for ele in sample_xshift_color:
            plt.plot(self.index + ele[1], self.data[ele[0]], c=ele[2])
