
from __future__ import division, print_function
import numpy as np
import pandas as pd


def uniform_prior(num_causes):
    """Uniform Prior
    """
    # All bins are given equal probability.
    # Most generic normalized prior.
    prior = np.full(num_causes, 1/num_causes)

    return prior


def jeffreys_prior(causes):
    """Jeffreys Prior
    """
    # All cause bins are given equal probability mass.
    # Best prior for x-ranges spanning decades.
    ln_factor = np.log(causes[-1] / causes[0])
    prior = 1 / (ln_factor * causes)
    # Want to make sure prior is normalized (i.e. prior.sum() == 1)
    prior = prior / prior.sum()

    return prior


def setup_prior(prior=None, num_causes=None):
    """Setup for prior array

    Parameters
    ----------
    prior : array_like, optional
        Prior distribution to use. If not specified, each cause bin will be
        assigned an equal probability (default is None, and uniform prior will
        be used).
    num_causes : int, optional
        Number of cause bins. Only needed for uniform prior (default is None).

    Returns
    -------
    prior : numpy.ndarray
        Normalized prior distribution.
    """
    if prior is None:
        assert num_causes is not None, 'num_causes must be specified for uniform prior'
        prior = uniform_prior(num_causes=num_causes)
    elif isinstance(prior, (list, tuple, np.ndarray, pd.Series)):
        prior = np.asarray(prior)
    else:
        raise TypeError('prior must be either None or array_like, '
                        'but got {}'.format(type(prior)))

    if not np.allclose(np.sum(prior), 1):
        raise ValueError('Prior (which is an array of probabilities) does '
                         'not add to 1. sum(prior) = {}'.format(np.sum(prior)))

    return prior
