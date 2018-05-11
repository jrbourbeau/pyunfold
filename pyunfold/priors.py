
from __future__ import division, print_function
import numpy as np
from six import string_types
import pandas as pd


def uniform_prior(num_causes):
    """Uniform Prior
    """
    # All bins are given equal probability.
    # Most generic normalized prior.
    prior = np.ones(num_causes) / num_causes

    return prior


def jeffreys_prior(xarray):
    """Jeffreys Prior
    """
    # All cause bins are given equal probability mass.
    # Best prior for x-ranges spanning decades.
    ln_factor = np.log(xarray[-1] / xarray[0])
    prior = 1 / (ln_factor * xarray)
    # Want to make sure prior is normalized (i.e. prior.sum() == 1)
    prior = prior / prior.sum()

    return prior


def setup_prior(priors='Jeffreys', num_causes=None, cause_axis=None):
    """Setup for prior array

    Parameters
    ----------
    priors : str or array_like, optional
        Prior distribution to use (default is 'Uniform', will use Uniform
        prior).
    num_causes : int, optional
        Number of cause bins. Only needed if priors='Jeffreys' (default is None).

    Returns
    -------
    prior : numpy.ndarray
        Normalized prior distribution.
    """
    if isinstance(priors, string_types) and priors.lower() == 'uniform':
        assert num_causes is not None, 'num_causes must be specified for uniform prior'
        prior = uniform_prior(num_causes=num_causes)
    elif isinstance(priors, string_types) and priors.lower() == 'jeffreys':
        assert cause_axis is not None, 'cause_axis must be specified for Jeffreys prior'
        prior = jeffreys_prior(xarray=cause_axis)
    elif isinstance(priors, (list, tuple, np.ndarray, pd.Series)):
        prior = np.asarray(priors)
    else:
        raise TypeError('priors must be either "Jeffreys" or array_like, '
                        'but got {}'.format(type(priors)))

    if not np.allclose(np.sum(prior), 1):
        raise ValueError('Prior (which is an array of probabilities) does '
                         'not add to 1. sum(prior) = {}'.format(np.sum(prior)))

    return prior
