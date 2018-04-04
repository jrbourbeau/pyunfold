
from __future__ import division, print_function
import numpy as np

from .utils import none_to_empty_list


class DataDist(object):
    """Class to define generic distribution w/assoc labels & axes
    """
    def __init__(self, name="", data=None, error=None, axis=None, edges=None,
                 xlabel="", ylabel="", units="", **kwargs):
        data, error, axis, edges = none_to_empty_list(data, error, axis, edges)
        self.name = name
        self.data = data
        self.nbins = len(data)
        self.error = error
        self.axis = axis
        self.edges = edges
        self.width = np.diff(edges)
        self.xlab = xlabel
        self.ylab = ylabel
        self.units = units
        # Systematic and Statistical Error
        self.sterr = np.zeros(self.nbins)
        self.syerr = np.zeros(self.nbins)
        self.CheckArrays()

    # Ensure that arrays are of proper lenghts
    def CheckArrays(self):
        if self.nbins != len(self.error):
            raise ValueError("{} data and error arrays unequal length!".format(self.name))
        if self.nbins != len(self.axis):
            raise ValueError("{} data and axis arrays unequal length!".format(self.name))
        if self.nbins != len(self.edges) - 1:
            raise ValueError("{} data and edges arrays improper length!".format(self.name))

    # Setter functions systematic and statistical error arrays
    def setStatErr(self, staterr):
        self.sterr = staterr.copy()

    def setSysErr(self, syserr):
        self.syerr = syserr.copy()
