
from __future__ import division, print_function
import pytest
import numpy as np

from pyunfold.priors import jeffreys_prior, setup_prior, uniform_prior


@pytest.mark.parametrize('prior', ['uniform',
                                   123,
                                   123.0])
def test_setup_prior_invalid_prior(prior):
    with pytest.raises(TypeError) as excinfo:
        setup_prior(prior)
    expected_msg = ('prior must be either None or array_like, '
                    'but got {}'.format(type(prior)))
    assert expected_msg == str(excinfo.value)


def test_setup_prior_non_normalized_raises():
    prior = [1, 2, 3, 4]
    with pytest.raises(ValueError) as excinfo:
        setup_prior(prior)
    expected_msg = ('Prior (which is an array of probabilities) does '
                    'not add to 1. sum(prior) = {}'.format(np.sum(prior)))
    assert expected_msg == str(excinfo.value)


def test_setup_prior_negative_raises():
    prior = [2, 0, -1]
    with pytest.raises(ValueError) as excinfo:
        setup_prior(prior)
    expected_msg = ('Input prior has negative values. Since the values '
                    'of prior are interpreted as probabilities, they '
                    'cannot be negative.')
    assert expected_msg == str(excinfo.value)


def test_jeffreys_prior():

    causes = np.linspace(5, 10, 4)

    # Calculate expected prior
    ln_factor = np.log(causes.max() / causes.min())
    prior = 1 / (ln_factor * causes)
    prior = prior / np.sum(prior)

    np.testing.assert_allclose(prior, jeffreys_prior(causes=causes))


def test_jeffreys_prior_normalized():
    causes = np.array([0.5, 1.5])
    prior = jeffreys_prior(causes=causes)

    np.testing.assert_allclose(prior.sum(), 1)


@pytest.mark.parametrize('type_', [list, tuple, np.array])
def test_jeffreys_prior_array_like(type_):
    causes = type_([1, 2, 3, 4, 5])
    jeffreys_prior(causes=causes)


@pytest.mark.parametrize('num_causes', [1, 7, 100])
def test_uniform_prior(num_causes):
    prior = uniform_prior(num_causes)
    # Correct number of cause bins
    assert len(prior) == num_causes
    # Every bin has same probability
    assert len(np.unique(prior))
    # Sum of probabilities add to one
    np.testing.assert_allclose(np.sum(prior), 1)
