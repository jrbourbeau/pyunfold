
from __future__ import division, print_function
import numpy as np


def JeffreysPrior(norm, xarray):
    """Jeffreys Prior
    """
    # All cause bins are given equal probability.
    # Best prior for x-ranges spanning decades.
    ilen = len(xarray)
    ln_factor = np.log(xarray[ilen - 1] / xarray[0])
    jprior = norm / (ln_factor * xarray)

    return jprior


def UserPrior(FuncList, xarray, norm):
    """User provided prior function
    """
    nAnalysisBins = len(FuncList)
    prior = []
    for ibin in range(nAnalysisBins):
        FuncString = FuncList[ibin]
        ixarray = xarray[ibin].copy()
        if (FuncString.lower() == "jeffreys"):
            iprior = JeffreysPrior(norm, ixarray)
        else:
            exec("def func(x): return {}".format(FuncString))
            iprior = func(ixarray)
        if prior == []:
            prior = iprior
        else:
            prior = np.append(prior, iprior)

    return prior
