
from __future__ import division, print_function
import numpy as np


def jeffreys_prior(norm, xarray):
    """Jeffreys Prior
    """
    # All cause bins are given equal probability.
    # Best prior for x-ranges spanning decades.
    ilen = len(xarray)
    ln_factor = np.log(xarray[ilen - 1] / xarray[0])
    jprior = norm / (ln_factor * xarray)

    return jprior


def user_prior(func, xarray, norm):
    """User provided prior function
    """
    if func.lower() == "jeffreys":
        prior = jeffreys_prior(norm, xarray)
    else:
        exec("def func(x): return {}".format(func))
        prior = func(xarray)

    return prior
