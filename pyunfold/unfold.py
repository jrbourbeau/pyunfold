
from __future__ import division, print_function
import warnings
import numpy as np
import pandas as pd

from .mix import Mixer
from .priors import setup_prior
from .utils import cast_to_array
from .callbacks import setup_callbacks_regularizer, CallbackList, TSStopping


def iterative_unfold(data=None, data_err=None, response=None,
                     response_err=None, efficiencies=None,
                     efficiencies_err=None, prior=None, ts=None,
                     ts_stopping=None, max_iter=100, cov_type='multinomial',
                     return_iterations=False, callbacks=None):
    """Performs iterative unfolding

    Parameters
    ----------
    data : array_like
        Input observed data distribution.
    data_err : array_like
        Uncertainties of the input observed data distribution. Must be the
        same shape as ``data``.
    response : array_like
        Response matrix.
    response_err : array_like
        Uncertainties of response matrix. Must be the same shape as
        ``response``.
    efficiencies : array_like
        Detection efficiencies for the cause distribution.
    efficiencies_err : array_like
        Uncertainties of detection efficiencies. Must be the same shape as
        ``efficiencies``.
    prior : array_like, optional
        Prior distribution to use in unfolding. If None, then a uniform
        (or flat) prior will be used. If array_like, then must have the same
        shape as ``efficiencies`` (default is None).
    ts : {'ks', 'chi2', 'bf', 'rmd'}
        Test statistic to use for stopping condition (default is 'ks').
        For more information about the available test statistics, see the
        `Test Statistics API documentation <api.rst#test-statistics>`__.
    ts_stopping : float, optional
        Test statistic stopping condition. At each unfolding iteration, the
        test statistic is computed between the current and previous iteration.
        Once the test statistic drops below ts_stopping, the unfolding
        procedure is stopped (default is 0.01).
    max_iter : int, optional
        Maximum number of iterations to allow (default is 100).
    cov_type : {'multinomial', 'poisson'}
        Whether to use the Multinomial or Poisson form for the covariance
        matrix (default is 'multinomial').
    return_iterations : bool, optional
        Whether to return unfolded distributions for each iteration
        (default is False).
    callbacks : list, optional
        List of ``pyunfold.callbacks.Callback`` instances to be applied during
        unfolding (default is None, which means no Callbacks are applied).

    Returns
    -------
    unfolded_result : dict
        Returned if ``return_iterations`` is False (default). Dictionary
        containing the final unfolded distribution, associated uncertainties,
        and test statistic information.

        The returned ``dict`` has the following keys:

            unfolded
                Final unfolded cause distribution
            stat_err
                Statistical uncertainties on the unfolded cause distribution
            sys_err
                Systematic uncertainties on the unfolded cause distribution
                associated with limited statistics in the response matrix
            ts_iter
                Final test statistic value
            ts_stopping
                Test statistic stopping criterion
            num_iterations
                Number of unfolding iterations
            unfolding_matrix
                Unfolding matrix
            prior
                Prior distribution used for unfolding

    unfolding_iters : pandas.DataFrame
        Returned if ``return_iterations`` is True. DataFrame containing the
        unfolded distribution, associated uncertainties, test statistic
        information, etc. at each iteration.

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
    >>> unfolded = iterative_unfold(data=data,
    ...                             data_err=data_err,
    ...                             response=response,
    ...                             response_err=response_err,
    ...                             efficiencies=efficiencies,
    ...                             efficiencies_err=efficiencies_err)
    >>> unfolded
    {'num_iterations': 4,
     'stat_err': array([11.16853268, 13.65488168]),
     'sys_err': array([0.65570621, 0.65570621]),
     'ts_iter': 0.0038300087456445975,
     'ts_stopping': 0.01,
     'unfolded': array([ 94.32086967, 155.67913033]),
     'unfolding_matrix': array([[0.8471473 , 0.1528527 ],
                                [0.06404093, 0.93595907]])}
    """
    # Validate user input
    inputs = {'data': data,
              'data_err': data_err,
              'response': response,
              'response_err': response_err,
              'efficiencies': efficiencies,
              'efficiencies_err': efficiencies_err
              }
    for name in inputs:
        if inputs[name] is None:
            raise ValueError('The input for {} must not be None.'.format(name))
        elif np.amin(inputs[name]) < 0:
            raise ValueError('The items in {} must be non-negative.'.format(name))

    data, data_err = cast_to_array(data, data_err)
    response, response_err = cast_to_array(response, response_err)
    efficiencies, efficiencies_err = cast_to_array(efficiencies,
                                                   efficiencies_err)

    num_causes = len(efficiencies)

    # Setup prior
    prior = setup_prior(prior=prior, num_causes=num_causes)

    # Define first prior counts distribution
    n_c = np.sum(data) * prior

    # Setup Mixer
    mixer = Mixer(data=data,
                  data_err=data_err,
                  efficiencies=efficiencies,
                  efficiencies_err=efficiencies_err,
                  response=response,
                  response_err=response_err,
                  cov_type=cov_type)

    # Set up callbacks, regularizer Callbacks are treated separately
    callbacks, regularizer = setup_callbacks_regularizer(callbacks)

    if ts is None and ts_stopping is None:
        # Setup legacy stopping condition behavior
        stopping = TSStopping(ts='ks',
                              ts_stopping=0.01,
                              num_causes=num_causes)
        callbacks = CallbackList([stopping] + callbacks.callbacks)
    elif any([ts, ts_stopping]):
        msg = ('Using ts and ts_stopping parameters in iterative_unfold has '
               'been deprecated. Please use the TSStopping callback instead.')
        warnings.warn(msg, DeprecationWarning)

        if ts is None:
            ts = 'ks'
        if ts_stopping is None:
            ts_stopping = 0.01
        stopping = TSStopping(ts=ts,
                              ts_stopping=ts_stopping,
                              num_causes=num_causes)
        callbacks = CallbackList([stopping] + callbacks.callbacks)

    unfolding_iters = _unfold(prior=n_c,
                              mixer=mixer,
                              max_iter=max_iter,
                              callbacks=callbacks,
                              regularizer=regularizer)

    if return_iterations:
        return unfolding_iters
    else:
        unfolded_result = dict(unfolding_iters.iloc[-1])
        return unfolded_result


def _unfold(prior=None, mixer=None, max_iter=100, callbacks=None, regularizer=None):
    """Perform iterative unfolding

    Parameters
    ----------
    prior : array_like
        Initial cause distribution.
    mixer : pyunfold.Mix.Mixer
        Mixer to perform the unfolding.
    max_iter : int, optional
        Maximum allowed number of iterations to perform.
    callbacks : list, optional
        List of ``pyunfold.callbacks.Callback`` instances to be applied during
        unfolding (default is None, which means no Callbacks are applied).
    regularizer : pyunfold.callbacks.Regularizer, optional
        Regularizer to smooth unfolded distribution at each iteration (default
        is None, which means no regularization is applied).

    Returns
    -------
    unfolding_iters : pandas.DataFrame
        DataFrame containing the unfolded result for each iteration.
        Each row in unfolding_result corresponds to an iteration.
    """
    callbacks.on_unfolding_begin()

    current_n_c = prior.copy()
    iteration = 0
    status = []
    stop_iterations = False
    while not stop_iterations and iteration < max_iter:
        callbacks.on_iteration_begin(iteration=iteration)

        # Perform unfolding for this iteration
        unfolded_n_c = mixer.smear(current_n_c)
        iteration += 1
        iteration_idx = iteration - 1
        status_current = {'unfolded': unfolded_n_c,
                          'stat_err': mixer.get_stat_err(),
                          'sys_err': mixer.get_MC_err(),
                          'num_iterations': iteration,
                          'unfolding_matrix': mixer.Mij,
                          'prior': current_n_c}
        status.append(status_current)

        if regularizer:
            # Will want the nonregularized distribution for the final iteration
            unfolded_nonregularized = status[iteration_idx]['unfolded'].copy()
            regularizer.on_iteration_end(iteration=iteration, status=status)

        callbacks.on_iteration_end(iteration=iteration, status=status)

        # Updated current distribution for next iteration of unfolding
        current_n_c = status[iteration_idx]['unfolded'].copy()

        # Check unfolding stopping condition
        stop_iterations = any([getattr(c, 'stop_iterations', None) for c in callbacks])

    # Convert status (list of dictionaries) to a pandas DataFrame
    unfolding_iters = pd.DataFrame.from_records(status)

    # Replace final folded iteration with un-regularized distribution
    if regularizer:
        last_iteration_index = unfolding_iters.index[-1]
        unfolding_iters.at[last_iteration_index, 'unfolded'] = unfolded_nonregularized

    callbacks.on_unfolding_end(status=status)

    return unfolding_iters
