
from __future__ import division, print_function
import numpy as np
import pandas as pd

from .LoadStats import make_mctables
from .Mix import Mixer
from .teststat import get_ts
from .datadist import DataDist
from .priors import UserPrior
from .utils import assert_same_shape, cast_to_array


def iterative_unfold(data, data_err, response, response_err, efficiencies,
                     efficiencies_err, priors='Jeffreys', ts='ks',
                     ts_stopping=0.01, max_iter=100, return_iterations=False):
    """Performs iterative Bayesian unfolding

    Parameters
    ----------
    data : array_like
        Input observed data distribution.
    data_err : array_like
        Uncertainties associated with the input observed data distribution.
        Must be the same shape as data.
    response : array_like
        Response matrix.
    response_err : array_like
        Response matrix errors.
    efficiencies : array_like
        Detection efficiencies for the observed data distribution.
    efficiencies_err : array_like
        Uncertainty in detection efficiencies.
    outfile : str
        Path where output ROOT file will be saved.
    priors : str or array_like, optional
        Prior distribution to use in unfolding. If 'Jeffreys', then the
        Jeffreys (flat) prior is used. Otherwise, must be array_like with
        same shape as data (default is 'Jeffreys').
    ts : {'ks', 'chi2', 'pf', 'rmd'}
        Name of test statistic to use for stopping condition (default is 'ks').
    ts_stopping : float, optional
        Test statistic stopping condition. At each unfolding iteration, the
        test statistic is computed between the current and previous iteration.
        Once the test statistic drops below ts_stopping, the unfolding
        procedure is stopped (default is 0.01).
    max_iter : int, optional
        Maximum number of iterations to allow (default is 100).
    return_iterations : bool, optional
        Whether to return unfolded distributions for each iteration
        (default is False).

    Returns
    -------
    unfolded_result : dict
        Returned if return_iterations is False (default). Final unfolded
        distribution and associated uncertainties.
    unfolding_iters : pandas.DataFrame
        Returned if return_iterations is True. DataFrame containing the
        unfolded distribution and associated uncertainties at each iteration.

    Examples
    --------
    >>> from pyunfold import iterative_unfold
    >>> data = [100, 150]
    >>> data_err = [10, 12.2]
    >>> response = [[0.9, 0.1],
    ...             [0.1, 0.9]]
    >>> response_err = [[0.01, 0.01],
    ...                 [0.01, 0.01]]
    >>> efficiencies = [1, 1]
    >>> efficiencies_err = [0.01, 0.01]
    >>> unfolded = iterative_unfold(data, data_err,
    ...                             response, response_err,
    ...                             efficiencies, efficiencies_err)
    >>> unfolded
    {'unfolded': array([ 94.48002622, 155.51997378]),
    'sys_err': array([0.66204237, 0.6620424 ]),
    'stat_err': array([11.2351567 , 13.75617997])}
    """
    # Validate user input
    data, data_err = cast_to_array(data, data_err)
    response, response_err = cast_to_array(response, response_err)
    efficiencies, efficiencies_err = cast_to_array(efficiencies,
                                                   efficiencies_err)

    assert_same_shape(data, data_err)
    assert_same_shape(efficiencies, efficiencies_err)
    assert_same_shape(response, response_err)
    assert len(data) == response.shape[0]
    assert len(efficiencies) == response.shape[1]
    if priors != 'Jeffreys':
        assert_same_shape(efficiencies, priors)

    mixer, ts_func, n_c = setup_mixer_ts_prior(data=data,
                                               data_err=data_err,
                                               priors=priors,
                                               efficiencies=efficiencies,
                                               efficiencies_err=efficiencies_err,
                                               response=response,
                                               response_err=response_err,
                                               ts=ts,
                                               ts_stopping=ts_stopping,
                                               max_iter=max_iter)
    unfolding_iters = perform_unfolding(n_c=n_c,
                                        mixer=mixer,
                                        ts_func=ts_func,
                                        max_iter=max_iter)

    if return_iterations:
        return unfolding_iters
    else:
        unfolded_result = dict(unfolding_iters.iloc[-1])
        return unfolded_result


def setup_mixer_ts_prior(data=None, data_err=None, priors='Jeffreys',
                         efficiencies=None, efficiencies_err=None,
                         response=None, response_err=None,
                         ts='ks', ts_stopping=0.01,
                         max_iter=100, cov_error='ACM'):

    if cov_error not in ['ACM', 'DCM']:
        raise ValueError('Invalid cov_error entered ({}). Must be in '
                         'either "ACM" or "DCM"'.format(cov_error))

    # Setup the Observed and MC Data Arrays
    # Load MC Stats (NCmc), Cause Efficiency (Eff) and Migration Matrix P(E|C)
    MCStats = make_mctables(efficiencies=efficiencies,
                            efficiencies_err=efficiencies_err,
                            response=response,
                            response_err=response_err)

    Cedges = np.arange(len(efficiencies) + 1, dtype=float)
    # Get bin midpoints
    Caxis = (Cedges[1:] + Cedges[:-1]) / 2

    # Load the Observed Data (data), define total observed events (n_obs)
    Exlab = 'Effects'
    Eylab = 'Counts'
    Etitle = 'effects histogram'
    Eedges = np.arange(len(data) + 1, dtype=float)
    # Get bin midpoints
    Eaxis = (Eedges[1:] + Eedges[:-1]) / 2
    EffDist = DataDist(Etitle, data=data, error=data_err,
                       axis=Eaxis, edges=Eedges, xlabel=Exlab,
                       ylabel=Eylab, units='')
    n_obs = np.sum(data)

    # Setup prior
    if priors == 'Jeffreys':
        n_c = UserPrior(['Jeffreys'], [Caxis], n_obs)
        n_c = n_c / np.sum(n_c)
    elif isinstance(priors, (list, tuple, np.ndarray, pd.Series)):
        n_c = np.asarray(priors)
    else:
        raise TypeError('priors must be either "Jeffreys" or array_like, '
                        'but got {}'.format(type(priors)))

    if not np.allclose(np.sum(n_c), 1):
        raise ValueError('Prior (which is an array of probabilities) does '
                         'not add to 1. sum(priors) = {}'.format(np.sum(n_c)))

    # Prepare Test Statistic-er
    ts_obj = get_ts(ts)
    ts_func = ts_obj(ts,
                     tol=ts_stopping,
                     Xaxis=Caxis,
                     TestRange=[0, 1e2],
                     verbose=False)

    # Prepare Mixer
    mixer = Mixer('SrMixALot',
                  ErrorType=cov_error,
                  MCTables=MCStats,
                  EffectsDist=EffDist)

    return mixer, ts_func, n_c


def perform_unfolding(n_c=None, mixer=None, ts_func=None, max_iter=100):
    """Perform iterative unfolding

    Parameters
    ----------
    n_c : array_like
        Cause distribution array.
    mixer : pyunfold.Mix.Mixer
        Mixer to perform the unfolding.
    ts_func : pyunfold.Utils.TestStat
        Test statistic object.
    max_iter : int, optional
        Maximum allowed number of iterations to perform.

    Returns
    -------
    unfolding_iters : pandas.DataFrame
        DataFrame containing the unfolded result for each iteration.
        Each row in unfolding_result corresponds to an iteration.
    """
    current_n_c = n_c.copy()
    counter = 0

    unfolding_iters = []

    while (not ts_func.pass_tol() and counter < max_iter):
        # Perform unfolding for this iteration
        unfolded_n_c = mixer.smear(current_n_c)

        # Add mixing result to unfolding_result
        status = {'unfolded': unfolded_n_c,
                  'stat_err': mixer.get_stat_err(),
                  'sys_err': mixer.get_MC_err()}
        unfolding_iters.append(status)

        ts_cur, ts_del, ts_prob = ts_func.GetStats(unfolded_n_c,
                                                   current_n_c)
        # Updated current distribution for next iteration of unfolding
        current_n_c = unfolded_n_c.copy()
        counter += 1

    # Convert unfolding_result dictionary to a pandas DataFrame
    columns = ['sys_err', 'unfolded', 'stat_err']
    unfolding_iters = pd.DataFrame.from_records(unfolding_iters,
                                                columns=columns)

    return unfolding_iters
