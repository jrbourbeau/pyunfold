
from __future__ import division, print_function
import pytest
import numpy as np

from pyunfold.mix import Mixer


num_causes = 5
num_effects = 7
data = np.arange(num_effects)
data_err = np.sqrt(data)
efficiencies = np.ones(num_causes, dtype=float)
efficiencies_err = np.full_like(efficiencies, 0.1)

np.random.seed(2)
response = 1 + np.random.rand(num_effects, num_causes)
response_err = np.sqrt(response)

mixer = Mixer(data=data,
              data_err=data_err,
              efficiencies=efficiencies,
              efficiencies_err=efficiencies_err,
              response=response,
              response_err=response_err,
              cov_type='multinomial')


def test_mixer_invalid_effects_bins_raises():
    # Add extra effect bin
    data = np.arange(num_effects + 1)
    data_err = np.sqrt(data)

    with pytest.raises(ValueError) as excinfo:
        Mixer(data=data,
              data_err=data_err,
              efficiencies=efficiencies,
              efficiencies_err=efficiencies_err,
              response=response,
              response_err=response_err)

    err_msg = ('Inconsistent number of effect bins. Observed data '
               'has {} effect bins, while response matrix has {} '
               'effect bins.'.format(len(data), num_effects))
    assert err_msg == str(excinfo.value)


def test_mixer_smear_invalid_cause_bins_raises():
    # Add extra cause bin
    prior = np.arange(num_causes + 1)

    mixer = Mixer(data=data,
                  data_err=data_err,
                  efficiencies=efficiencies,
                  efficiencies_err=efficiencies_err,
                  response=response,
                  response_err=response_err)
    with pytest.raises(ValueError) as excinfo:
        mixer.smear(prior)

    err_msg = ('Trying to unfold with the wrong number of causes. '
               'Response matrix has {} cause bins, while prior '
               'has {} cause bins.'.format(num_causes, len(prior)))
    assert err_msg == str(excinfo.value)


def test_mixer_zeros_prior():
    np.testing.assert_array_equal(mixer.smear(np.zeros(num_causes)),
                                  np.zeros(num_causes))


def test_mixer_2d_response_check():
    # Tests that a non-2D response matrix raises an error
    # Construct a response matrix that is 3-dimensional
    response_bad = 1 + np.random.rand(num_effects, num_causes, 1)
    response_bad_err = np.sqrt(response_bad)

    with pytest.raises(ValueError) as excinfo:
        Mixer(data=data,
              data_err=data_err,
              efficiencies=efficiencies,
              efficiencies_err=efficiencies_err,
              response=response_bad,
              response_err=response_bad_err)

    raised_msg = str(excinfo.value)
    assert 'Response matrix must be 2-dimensional' in raised_msg
    assert '{}-dimensional response'.format(response_bad.ndim) in raised_msg
