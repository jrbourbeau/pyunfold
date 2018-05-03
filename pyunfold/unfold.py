
from __future__ import division, print_function
import numpy as np
import pandas as pd

from .mix import Mixer
from .teststat import get_ts
from .priors import setup_prior
from .utils import cast_to_array
from .callbacks import validate_callbacks, extract_regularizer


def iterative_unfold(data, data_err, response, response_err, efficiencies,
                     efficiencies_err, priors='Jeffreys', ts='ks',
                     ts_stopping=0.01, max_iter=100, return_iterations=False,
                     callbacks=None):
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
    callbacks : list, optional
        List of ``pyunfold.callbacks.Callback`` instances to be applied during
        unfolding (default is None, which means no Callbacks are applied).

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

    num_causes = len(efficiencies)

    # Setup prior
    n_c = setup_prior(priors=priors,
                      num_causes=num_causes,
                      num_observations=np.sum(data))

    # Setup Mixer
    mixer = Mixer(error_type='ACM',
                  data=data,
                  data_err=data_err,
                  efficiencies=efficiencies,
                  efficiencies_err=efficiencies_err,
                  response=response,
                  response_err=response_err)

    # Setup test statistic
    ts_obj = get_ts(ts)
    ts_func = ts_obj(ts,
                     tol=ts_stopping,
                     num_causes=num_causes,
                     TestRange=[0, 1e2],
                     verbose=False)

    unfolding_iters = _unfold(prior=n_c,
                              mixer=mixer,
                              ts_func=ts_func,
                              max_iter=max_iter,
                              callbacks=callbacks)

    if return_iterations:
        return unfolding_iters
    else:
        unfolded_result = dict(unfolding_iters.iloc[-1])
        return unfolded_result


def _unfold(prior=None, mixer=None, ts_func=None, max_iter=100,
            callbacks=None):
    """Perform iterative unfolding

    Parameters
    ----------
    prior : array_like
        Initial cause distribution.
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
    callbacks = validate_callbacks(callbacks)
    regularizer = extract_regularizer(callbacks)
    # Will treat regularizer Callback separately
    callbacks = [c for c in callbacks if c is not regularizer]

    current_n_c = prior.copy()
    iteration = 0
    unfolding_iters = []
    while (not ts_func.pass_tol() and iteration < max_iter):

        # Perform unfolding for this iteration
        unfolded_n_c = mixer.smear(current_n_c)
        iteration += 1
        status = {'unfolded': unfolded_n_c,
                  'stat_err': mixer.get_stat_err(),
                  'sys_err': mixer.get_MC_err()}

        if regularizer is not None:
            # Will want the nonregularized distribution for the final iteration
            unfolded_nonregularized = unfolded_n_c.copy()
            unfolded_n_c = regularizer.on_iteration_end(iteration=iteration,
                                                        params=status)
            status['unfolded'] = unfolded_n_c

        ts_cur, ts_del, ts_prob = ts_func.GetStats(unfolded_n_c, current_n_c)
        status['ts_iter'] = ts_cur
        status['ts_stopping'] = ts_func.tol

        for callback in callbacks:
            callback.on_iteration_end(iteration=iteration,
                                      params=status)

        unfolding_iters.append(status)
        # Updated current distribution for next iteration of unfolding
        current_n_c = unfolded_n_c.copy()

    # Convert unfolding_iters list of dictionaries to a pandas DataFrame
    unfolding_iters = pd.DataFrame.from_records(unfolding_iters)

    # Replace final folded iteration with un-regularized distribution
    if regularizer:
        last_iteration_index = unfolding_iters.index[-1]
        unfolding_iters.at[last_iteration_index, 'unfolded'] = unfolded_nonregularized

    return unfolding_iters
