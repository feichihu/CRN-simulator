import matplotlib.pyplot as plt
import pandas as pd

from crn import Species

class Simulation:
    """
    A completed simulation of a CRN. Contains the time-series information
    of every species in the CRN throughout the simulation. Allows for quick
    plotting and extraction of simulation data.

    This class probably won't be constructed by a user, thus it's
    implementation is more internal.

    args:
        sim: Dict[crn.Species, np.ndarray]
            A dictionary of species name to concentration time series.
            This dictionary also contains other fields such as "time" and
            "nothing" which are used for plotting data. This dictionary
            makes no guarantees that these are the only things it will
            contain, it may have new fields that are added if they are
            needed.
    """
    def __init__(self, sim, stochastic=False):
        self.sim = sim
        self.stochastic = stochastic
        #print(sim)
        #self.time = sim["time"]
        #self.reactions = sim.get("reactions", None)

        #del sim["time"]
        #del sim["reactions"]

    def __getitem__(self, s):
        if type(s) is not Species and s!='time':
            raise ValueError(
                "Simulation.__getitem__: tried to get item of non-species. "
                "Type of key must be Species. The type of the key "
                f"passed was {type(s)}")

        return self.sim[s]

    def plot(self, filename=None, title=None):
        """
        Plots the concentration of all of the species over time.

        args:
            filename: Optional[str]
                if present, save the plot to a file `filename`. Otherwise,
                the plot will show up as a new window.

            title: Optional[str]
                if present, the plot will have a title `title`.
        """
        res = pd.DataFrame()
        if filename:
            backend = plt.get_backend()
            plt.switch_backend("Svg")
        #for species, series in sorted(self.sim.items()):
        #print(self.sim)
        res['time'] = self.sim['time']
        i = 1
        for species, series in self.sim.items():
            series = self.sim[species]
            #print(species)
            #print(series)
            #if species.name != "nothing":
            if species != "nothing" and species != "time":
                plt.plot(self.sim['time'], series, label=f"[{species}]")
                name = str(species)
                i = i+1
                res[name] = series
                #print(series)
                #plt.savefig('results.png')
        res.to_csv('res.csv')

        plt.xlabel("time (seconds)")
        if self.stochastic:
            plt.ylabel("molecule counts")
        else:
            plt.ylabel("concentration (M)")

        plt.legend(loc="best")
        if title:
            plt.title(title)

        if filename:
            plt.savefig(filename)
        else:
            plt.show(block=True)

        if filename:
            plt.switch_backend(backend)


