
from __future__ import division, print_function
import numpy as np
import pandas as pd


def uniform_prior(num_causes):
    """Convenience function to calculate uniform prior distribution

    Parameters
    ----------
    num_causes : int
        Number of cause bins.

    Returns
    -------
    prior : numpy.ndarray
        Normalized uniform prior distribution.

    Examples
    --------
    >>> from pyunfold.priors import uniform_prior
    >>> uniform_prior(num_causes=4)
    array([0.25, 0.25, 0.25, 0.25])
    """
    # All bins are given equal probability.
    # Most generic normalized prior.
    prior = np.full(num_causes, 1/num_causes)

    return prior


def jeffreys_prior(causes):
    """Convenience function to calculate Jeffreys prior distribution

    Parameters
    ----------
    causes : array_like
        Midpoint value of cause bins. For instance if cause bin edges are
        given by [0, 2, 4], then ``causes`` is [1, 3].

    Returns
    -------
    prior : numpy.ndarray
        Normalized Jeffreys prior distribution.

    Notes
    -----
    The Jeffreys prior is defined as

    .. math::

        P(C_{\\mu})^{\\text{Jeffreys}} = \\frac{1}{\\log(C_{\\text{max}}/C_\\text{min})C_{\\mu}}

    for cause bin values :math:`C_{\\mu}` and maximum/minimum cause values
    :math:`C_{\\text{max}}`/:math:`C_{\\text{min}}`. For more details regarding
    Jeffreys prior see [1]_.

    References
    ----------
    .. [1] Jeffreys, H. "An Invariant Form for the Prior Probability in Estimation Problems".
        *Proc. of the Royal Society of London A: Mathematical, Physical and Engineering Sciences*
        186 (1007). London, England:453-61. `<https://doi.org/10.1098/rspa.1946.0056>`_.

    Examples
    --------
    >>> from pyunfold.priors import jeffreys_prior
    >>> causes = [1, 2, 3, 4]
    >>> jeffreys_prior(causes=causes)
    array([0.48, 0.24, 0.16, 0.12])

    """
    causes = np.asarray(causes)
    # All cause bins are given equal probability mass.
    # Best prior for x-ranges spanning decades.
    ln_factor = np.log(causes.max() / causes.min())
    prior = 1 / (ln_factor * causes)
    # Want to make sure prior is normalized (i.e. prior.sum() == 1)
    prior = prior / np.sum(prior)

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
    if np.amin(prior) < 0:
        raise ValueError('Input prior has negative values. Since the values '
                         'of prior are interpreted as probabilities, they '
                         'cannot be negative.')

    return prior
