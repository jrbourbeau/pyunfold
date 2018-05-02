
from __future__ import division, print_function
import numpy as np
from six import string_types
import pandas as pd


def jeffreys_prior(norm, xarray):
    """Jeffreys Prior
    """
    # All cause bins are given equal probability.
    # Best prior for x-ranges spanning decades.
    ln_factor = np.log(xarray[-1] / xarray[0])
    prior = norm / (ln_factor * xarray)
    # Want to make sure prior is normalized (i.e. prior.sum() == 1)
    prior = prior / prior.sum()

    return prior


def setup_prior(priors='Jeffreys', num_causes=None, num_observations=None):
    """Setup for prior array

    Parameters
    ----------
    priors : str or array_like, optional
        Prior distribution to use (default is 'Jeffreys', will use Jeffreys
        prior).
    num_causes : int, optional
        Number of cause bins. Only needed if priors='Jeffreys' (default is None).
    num_observations : int, optional
        Number of total effect observations. Only needed if priors='Jeffreys'
        (default is None).

    Returns
    -------
    prior : numpy.ndarray
        Normalized prior distribution.
    """
    if isinstance(priors, string_types) and priors.lower() == 'jeffreys':
        assert num_causes is not None, 'num_causes must be specified for Jeffreys prior'
        assert num_observations is not None, 'num_observations must be specified for Jeffreys prior'
        cause_bin_edges = np.arange(num_causes + 1, dtype=float)
        cause_axis = (cause_bin_edges[1:] + cause_bin_edges[:-1]) / 2
        prior = jeffreys_prior(norm=num_observations, xarray=cause_axis)
    elif isinstance(priors, (list, tuple, np.ndarray, pd.Series)):
        prior = np.asarray(priors)
    else:
        raise TypeError('priors must be either "Jeffreys" or array_like, '
                        'but got {}'.format(type(priors)))

    if not np.allclose(np.sum(prior), 1):
        raise ValueError('Prior (which is an array of probabilities) does '
                         'not add to 1. sum(prior) = {}'.format(np.sum(prior)))

    return prior
