
from __future__ import division, print_function, unicode_literals
import os
import sys
import numpy as np
import pandas as pd
import pytest

from .testing_utils import diagonal_response, triangular_response

from pyunfold.unfold import iterative_unfold
from pyunfold.priors import jeffreys_prior

PY_VERSION = sys.version_info.major

columns = ['sys_err', 'unfolded', 'stat_err']


@pytest.mark.parametrize('response_type', ['diagonal', 'triangular'])
def test_iterative_unfold(response_type):
    # Load test counts distribution and diagonal response matrix
    np.random.seed(2)

    samples = np.random.normal(loc=0, scale=1, size=int(1e5))

    bins = np.linspace(-1, 1, 10)
    counts, _ = np.histogram(samples, bins=bins)
    counts_err = np.sqrt(counts)
    if response_type == 'diagonal':
        response, response_err = diagonal_response(len(counts))
    else:
        response, response_err = triangular_response(len(counts))
    efficiencies = np.ones_like(counts, dtype=float)
    efficiencies_err = np.full_like(efficiencies, 0.001)

    unfolded_result = iterative_unfold(counts, counts_err,
                                       response, response_err,
                                       efficiencies, efficiencies_err,
                                       ts='ks',
                                       ts_stopping=0.01,
                                       return_iterations=False)

    unfolded_counts = unfolded_result['unfolded']

    # Given diagonal response matrix, unfolded counts should be same as measured counts
    if response_type == 'diagonal':
        assert np.allclose(counts, unfolded_counts)
    else:
        assert not np.allclose(counts, unfolded_counts)


def test_iterative_unfold_max_iter():
    # Load test counts distribution and diagonal response matrix
    np.random.seed(2)

    samples = np.random.normal(loc=10, scale=1, size=int(1e5))

    bins = np.linspace(8, 12, 10)
    counts, _ = np.histogram(samples, bins=bins)
    counts_err = np.sqrt(counts)
    response, response_err = triangular_response(len(counts))
    efficiencies = np.ones_like(counts, dtype=float)
    efficiencies_err = np.full_like(efficiencies, 0.001)

    causes = (bins[1:] + bins[:-1]) / 2
    prior = jeffreys_prior(causes=causes)

    max_iter = 5
    unfolded_result = iterative_unfold(counts, counts_err,
                                       response, response_err,
                                       efficiencies, efficiencies_err,
                                       prior=prior,
                                       ts='ks',
                                       ts_stopping=0.00001,
                                       max_iter=max_iter,
                                       return_iterations=True)

    assert unfolded_result.shape[0] == max_iter


def test_example():
    pytest.importorskip('tables')
    here = os.path.abspath(os.path.dirname(__file__))
    test_file = os.path.join(here,
                             'test_data',
                             'example1_python{}.hdf'.format(PY_VERSION))
    expected = pd.read_hdf(test_file)

    # Run example case
    data = [100, 150]
    data_err = [10, 12.2]
    response = [[0.9, 0.1],
                [0.1, 0.9]]
    response_err = [[0.01, 0.01],
                    [0.01, 0.01]]
    efficiencies = [0.4, 0.67]
    efficiencies_err = [0.01, 0.01]
    causes = np.arange(len(data)) + 0.5
    prior = jeffreys_prior(causes=causes)

    # Perform iterative unfolding
    unfolded = iterative_unfold(data=data,
                                data_err=data_err,
                                response=response,
                                response_err=response_err,
                                efficiencies=efficiencies,
                                efficiencies_err=efficiencies_err,
                                prior=prior,
                                return_iterations=True)
    unfolded = unfolded[columns]

    pd.testing.assert_frame_equal(unfolded, expected)


def test_example_2():
    pytest.importorskip('tables')
    here = os.path.abspath(os.path.dirname(__file__))
    test_file = os.path.join(here,
                             'test_data',
                             'example2_python{}.hdf'.format(PY_VERSION))
    expected = pd.read_hdf(test_file)

    # Run example case
    data = [100, 150]
    data_err = [10, 12.2]
    response = [[0.8, 0.1],
                [0.2, 0.9]]
    response_err = [[0.01, 0.01],
                    [0.01, 0.01]]
    efficiencies = [0.4, 0.67]
    efficiencies_err = [0.01, 0.01]
    prior = [0.34, 1-0.34]
    # Perform iterative unfolding
    unfolded = iterative_unfold(data, data_err,
                                response, response_err,
                                efficiencies, efficiencies_err,
                                prior=prior,
                                return_iterations=True)
    unfolded = unfolded[columns]

    pd.testing.assert_frame_equal(unfolded, expected)


def test_example_non_square_response():
    pytest.importorskip('tables')
    here = os.path.abspath(os.path.dirname(__file__))
    test_file = os.path.join(here,
                             'test_data',
                             'example3_python{}.hdf'.format(PY_VERSION))
    expected = pd.read_hdf(test_file)

    # Run example case
    data = [100, 150]
    data_err = [10, 12.2]
    response = [[0.8, 0.1, 0.6],
                [0.2, 0.9, 0.4]]
    response_err = [[0.01, 0.01, 0.01],
                    [0.01, 0.01, 0.01]]
    efficiencies = [0.4, 0.67, 0.8]
    efficiencies_err = [0.01, 0.01, 0.01]
    prior = [0.34, 0.21, 1 - (0.34 + 0.21)]
    # Perform iterative unfolding
    unfolded = iterative_unfold(data, data_err,
                                response, response_err,
                                efficiencies, efficiencies_err,
                                prior=prior,
                                return_iterations=True)
    unfolded = unfolded[columns]

    pd.testing.assert_frame_equal(unfolded, expected)


inputs = ['data', 'data_err', 'response', 'response_err', 'efficiencies', 'efficiencies_err']


@pytest.mark.parametrize('none_input', inputs)
def test_iterative_unfold_none_input_raises(none_input, example_dataset):
    inputs = {'data': example_dataset.data,
              'data_err': example_dataset.data_err,
              'response': example_dataset.response,
              'response_err': example_dataset.response_err,
              'efficiencies': example_dataset.efficiencies,
              'efficiencies_err': example_dataset.efficiencies_err
              }
    inputs[none_input] = None
    with pytest.raises(ValueError) as excinfo:
        iterative_unfold(**inputs)
    expected_msg = 'The input for {} must not be None.'.format(none_input)
    assert expected_msg == str(excinfo.value)


@pytest.mark.parametrize('negative_input', inputs)
def test_iterative_unfold_negative_input_raises(negative_input, example_dataset):
    inputs = {'data': example_dataset.data,
              'data_err': example_dataset.data_err,
              'response': example_dataset.response,
              'response_err': example_dataset.response_err,
              'efficiencies': example_dataset.efficiencies,
              'efficiencies_err': example_dataset.efficiencies_err
              }
    inputs[negative_input][0] *= -1
    with pytest.raises(ValueError) as excinfo:
        iterative_unfold(**inputs)
    expected_msg = 'The items in {} must be non-negative.'.format(negative_input)
    assert expected_msg == str(excinfo.value)


@pytest.mark.parametrize('cov_type_1,cov_type_2', [('multinomial', 'Multinomial'),
                                                   ('poisson', 'Poisson')])
def test_iterative_unfold_cov_type_case_insensitive(cov_type_1, cov_type_2, example_dataset):
    # Test that cov_type is case insensitive
    inputs = {'data': example_dataset.data,
              'data_err': example_dataset.data_err,
              'response': example_dataset.response,
              'response_err': example_dataset.response_err,
              'efficiencies': example_dataset.efficiencies,
              'efficiencies_err': example_dataset.efficiencies_err
              }
    result_1 = iterative_unfold(cov_type=cov_type_1, **inputs)
    result_2 = iterative_unfold(cov_type=cov_type_2, **inputs)

    assert result_1.keys() == result_2.keys()
    for key in result_1.keys():
        if isinstance(result_1[key], np.ndarray):
            np.testing.assert_array_equal(result_1[key], result_2[key])
        else:
            assert result_1[key] == result_2[key]


def test_iterative_unfold_cov_type_raises(example_dataset):
    inputs = {'data': example_dataset.data,
              'data_err': example_dataset.data_err,
              'response': example_dataset.response,
              'response_err': example_dataset.response_err,
              'efficiencies': example_dataset.efficiencies,
              'efficiencies_err': example_dataset.efficiencies_err
              }
    cov_type = 'Not a valid cov_type'
    with pytest.raises(ValueError) as excinfo:
        iterative_unfold(cov_type=cov_type, **inputs)
    expected_msg = ('Invalid pec_cov_type entered: {}. Must be '
                    'either "multinomial" or "poisson".'.format(cov_type))
    assert expected_msg == str(excinfo.value)


@pytest.mark.parametrize('return_iterations', [True, False])
def test_iterative_unfold_keys(example_dataset, return_iterations):
    inputs = {'data': example_dataset.data,
              'data_err': example_dataset.data_err,
              'response': example_dataset.response,
              'response_err': example_dataset.response_err,
              'efficiencies': example_dataset.efficiencies,
              'efficiencies_err': example_dataset.efficiencies_err
              }
    unfolded_result = iterative_unfold(return_iterations=return_iterations,
                                       **inputs)
    keys = ['unfolded',
            'stat_err',
            'sys_err',
            'num_iterations',
            'ts_stopping',
            'ts_iter',
            'unfolding_matrix',
            ]
    for key in keys:
        assert key in unfolded_result


def test_unfolding_matrix(example_dataset):
    inputs = {'data': example_dataset.data,
              'data_err': example_dataset.data_err,
              'response': example_dataset.response,
              'response_err': example_dataset.response_err,
              'efficiencies': example_dataset.efficiencies,
              'efficiencies_err': example_dataset.efficiencies_err
              }
    unfolded_result = iterative_unfold(return_iterations=True, **inputs)
    for _, row in unfolded_result.iterrows():
        np.testing.assert_array_equal(row['unfolded'],
                                      np.dot(inputs['data'], row['unfolding_matrix']))
