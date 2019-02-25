"""
Read and process individual csv data file generated by the pressure sensor.
"""
import pandas as pd
from matplotlib import pyplot as plt


class Sensor:

    def __init__(self, description="pressure sensor data"):
        self.description = description
        self.data = None
        self.path = ""

    def __repr__(self):
        return(self.description)

    def read(self, csv_file, sample=None):
        """
        Read a csv file generated by pressure sensor and add it to existing
        data.

        Parameters
        ----------
        csv_file: string, csv file name such as "xxx.csv" or "xxx".
        sample: string, column name for this data. Use file name if sample is
            not given.
        """

        if not csv_file.endswith(".csv"):
            csv_file = csv_file + ".csv"

        if sample is None:
            sample = csv_file[:-4]

        csv_file = self.path + csv_file
        data = pd.read_csv(csv_file,
                           skiprows=2,
                           header=None,
                           names=["time", sample])
        data.index = data.time
        data = data.drop("time", 1)

        if self.data is None:
            data_combined = data
        else:
            if sample not in self.data.columns:
                data_combined = pd.concat([self.data, data], 1)
            else:
                data_combined = self.data

        self.data = data_combined
        self.index = data_combined.index

    def trim(self, position="first"):
        # remove data point after first or last maximun
        for col in self.data.columns:
            value = self.data[col]
            if position == "first":
                index_of_max = value[value == max(value)].index[0]
            elif position == "last":
                index_of_max = value[value == max(value)].index[-1]
            self.data[col] = value[:index_of_max]
        return(self)

    def to_pressure(self, base=272):
        # use mmHg
        self.data = (self.data - base) / 30.7 * 43.6
        return(self.data)

    def to_force(self, base=272):
        # so we can compare force measured with force gauge
        self.data = (self.data - base) / 30.7
        return(self.data)

    def plot(self, samples=None):
        if samples is None:
            samples = self.data.columns
        for spl in samples:
            plt.plot(self.data.index, self.data[spl])
        plt.legend()
