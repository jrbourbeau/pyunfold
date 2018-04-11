
import numpy as np
import pytest

from pyunfold.unfold import iterative_unfold
from pyunfold.callbacks import Callback, Logger, SplineRegularizer, Regularizer

# Example case
num_samples = int(1e4)
true_samples = np.random.normal(loc=0.0, scale=1.0, size=num_samples)
random_noise = np.random.normal(loc=0.3, scale=0.5, size=num_samples)
observed_samples = true_samples + random_noise
bins = np.linspace(-3, 3, 21)
data_true, _ = np.histogram(true_samples, bins=bins)
data_true = np.array(data_true, dtype=float)
data, _ = np.histogram(observed_samples, bins=bins)
data_err = np.sqrt(data)
response, _, _ = np.histogram2d(observed_samples, true_samples, bins=bins)
response_err = np.sqrt(response)
response = response / response.sum(axis=0)
response_err = response_err / response.sum(axis=0)
efficiencies = response.sum(axis=0)
efficiencies_err = np.full_like(efficiencies, 0.1, dtype=float)


@pytest.mark.parametrize('callbacks', [[Logger()], Logger()])
def test_logger(capsys, callbacks):

    # Perform iterative unfolding
    unfolded_results = iterative_unfold(data, data_err,
                                        response, response_err,
                                        efficiencies, efficiencies_err,
                                        return_iterations=True,
                                        callbacks=callbacks)

    # Get stdout and std err from iterative_unfold
    out, err = capsys.readouterr()

    # Build expected output
    expected_output = ''
    for row_index, row in unfolded_results.iterrows():
        row_output = ('Iteration {}: ts = {:0.4f}, ts_stopping ='
                      ' {}\n'.format(row_index + 1,
                                     row['ts_iter'],
                                     row['ts_stopping']))
        expected_output += row_output

    assert expected_output == out


def test_Logger_is_Callback():

    logger = Logger()
    assert isinstance(logger, Callback)


def test_SplineRegularizer_is_Regularizer():

    spline_reg = SplineRegularizer()
    assert isinstance(spline_reg, Regularizer)


def test_SplineRegularizer_passes():

    spline_reg = SplineRegularizer(degree=3)
    # Perform iterative unfolding
    iterative_unfold(data=data,
                     data_err=data_err,
                     response=response,
                     response_err=response_err,
                     efficiencies=efficiencies,
                     efficiencies_err=efficiencies_err,
                     return_iterations=True,
                     callbacks=[spline_reg])
