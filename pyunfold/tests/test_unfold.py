
from __future__ import division, print_function, unicode_literals
import os
import sys
import numpy as np
import pandas as pd
import pytest

from .testing_utils import diagonal_response, triangular_response

from pyunfold.unfold import iterative_unfold

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
                                       priors='Jeffreys',
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

    samples = np.random.normal(loc=0, scale=1, size=int(1e5))

    bins = np.linspace(-1, 1, 10)
    counts, _ = np.histogram(samples, bins=bins)
    counts_err = np.sqrt(counts)
    response, response_err = triangular_response(len(counts))
    efficiencies = np.ones_like(counts, dtype=float)
    efficiencies_err = np.full_like(efficiencies, 0.001)

    max_iter = 5
    unfolded_result = iterative_unfold(counts, counts_err,
                                       response, response_err,
                                       efficiencies, efficiencies_err,
                                       priors='Jeffreys',
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
    # Perform iterative unfolding
    unfolded = iterative_unfold(data, data_err,
                                response, response_err,
                                efficiencies, efficiencies_err,
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
    priors = [0.34, 1-0.34]
    # Perform iterative unfolding
    unfolded = iterative_unfold(data, data_err,
                                response, response_err,
                                efficiencies, efficiencies_err,
                                priors=priors,
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
    priors = [0.34, 0.21, 1 - (0.34 + 0.21)]
    # Perform iterative unfolding
    unfolded = iterative_unfold(data, data_err,
                                response, response_err,
                                efficiencies, efficiencies_err,
                                priors=priors,
                                return_iterations=True)
    unfolded = unfolded[columns]

    pd.testing.assert_frame_equal(unfolded, expected)
