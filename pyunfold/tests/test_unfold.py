
import os
import numpy as np
from ROOT import TFile
from root_numpy import hist2array
import pytest

from .testing_utils import (save_test_root_file, diagonal_response,
                                    triangular_response)

from pyunfold.unfold import iterative_unfold, unfold


@pytest.mark.parametrize('response_type', ['diagonal', 'triangular'])
def test_iterative_unfold(response_type):
    # Load test counts distribution and diagonal response matrix
    np.random.seed(2)

    samples = np.random.normal(loc=0, scale=1, size=int(1e5))

    bins = np.linspace(-5, 5, 100)
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

    bins = np.linspace(-5, 5, 100)
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
