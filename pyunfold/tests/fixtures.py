
from collections import namedtuple
import numpy as np
import pytest


@pytest.fixture()
def example_dataset():
    data_attributes = ['data',
                       'data_err',
                       'response',
                       'response_err',
                       'efficiencies',
                       'efficiencies_err',
                       ]
    ExampleInput = namedtuple('ExampleInput', data_attributes)

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

    example = ExampleInput(data=data,
                           data_err=data_err,
                           response=response,
                           response_err=response_err,
                           efficiencies=efficiencies,
                           efficiencies_err=efficiencies_err)

    return example
