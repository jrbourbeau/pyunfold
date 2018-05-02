
from __future__ import division, print_function
import pytest
import numpy as np

from pyunfold.priors import jeffreys_prior, setup_prior


@pytest.mark.parametrize('prior', ['not jeffreys',
                                   123,
                                   123.0])
def test_setup_prior_invalid_prior(prior):
    with pytest.raises(TypeError) as excinfo:
        setup_prior(prior)
    expected_msg = ('priors must be either "Jeffreys" or array_like, '
                    'but got {}'.format(type(prior)))
    assert expected_msg == str(excinfo.value)


def test_setup_prior_non_normalized_raises():
    prior = [1, 2, 3, 4]
    with pytest.raises(ValueError) as excinfo:
        setup_prior(prior)
    expected_msg = ('Prior (which is an array of probabilities) does '
                    'not add to 1. sum(prior) = {}'.format(np.sum(prior)))
    assert expected_msg == str(excinfo.value)


@pytest.mark.parametrize('prior', ['jeffreys',
                                   'Jeffreys',
                                   'JEFFREYS'])
def test_setup_prior_jeffreys_prior_spelling(prior):
    setup_prior(prior, num_causes=10, num_observations=10)


@pytest.mark.parametrize('num_causes,num_observations', [(10, None),
                                                         (None, 10)])
def test_setup_prior_jeffreys_prior_extra_raises(num_causes, num_observations):
    with pytest.raises(AssertionError):
        setup_prior('jeffreys',
                    num_causes=num_causes,
                    num_observations=num_observations)


def test_jeffreys_prior_normalized():
    xarray = np.array([0.5, 1.5])
    prior = jeffreys_prior(norm=10, xarray=xarray)

    np.testing.assert_allclose(prior.sum(), 1)
