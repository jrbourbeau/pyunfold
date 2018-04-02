
import os
import numpy as np
from ROOT import TH1F, TH2F, TFile


def diagonal_response(size):
    """Returns diagonal response matrix
    """
    # Diagonal response matrix
    response = np.eye(size)
    response_err = np.sqrt(response)

    return response, response_err


def triangular_response(size):
    """Returns upper-triangular response matrix
    """
    # Triangular response matrix
    response = np.zeros((size, size))
    inds = np.triu_indices(size)
    response[inds] = 1
    response_err = np.sqrt(response)

    return response, response_err
