
import os
import numpy as np
from ROOT import TH1F, TH2F, TFile


def generate_diagonal_response(size):
    """Returns diagonal response matrix
    """
    # Diagonal response matrix
    response = np.eye(size)
    response_err = np.sqrt(response)

    return response, response_err


def generate_triangular_response(size):
    """Returns upper-triangular response matrix
    """
    # Triangular response matrix
    response = np.zeros((size, size))
    inds = np.triu_indices(size)
    response[inds] = 1
    response_err = np.sqrt(response)

    return response, response_err


def save_test_root_file(outfile, response='diagonal'):
    """Saves ROOT file with generated test data

    Parameters
    ----------
    outfile : str
        Path where output ROOT file will be saved.
    response : {'diagonal', 'triangular'}
        Whether to use a diagonal or upper-triangular response matrix.
    """
    if response not in ['diagonal', 'triangular']:
        raise ValueError('Invalid response, {}, entered. Must be either '
                         '"diagonal" or "triangular".'.format(response))
    if os.path.exists(outfile):
        os.remove(outfile)

    fout = TFile(outfile , 'UPDATE')
    # Bin Definitions
    binname = 'bin0'
    pdir = fout.mkdir(binname, 'Bin number 0')
    # Go to home of ROOT file
    fout.cd(binname)

    # Load test counts distribution and diagonal response matrix
    np.random.seed(2)

    samples = np.random.normal(loc=0, scale=1, size=int(1e5))

    bins = np.linspace(-5, 5, 100)
    counts, _ = np.histogram(samples, bins=bins)
    counts_err = np.sqrt(counts)

    if response == 'diagonal':
        response_array, response_err_array = generate_diagonal_response(len(counts))
    else:
        response_array, response_err_array = generate_triangular_response(len(counts))

    efficiencies = np.ones_like(counts, dtype=float)
    efficiencies_err = np.full_like(efficiencies, 0.001)

    cbins = len(counts)+1
    carray = np.arange(cbins, dtype=float)

    ebins = len(counts)+1
    earray = np.arange(ebins, dtype=float)
    cbins -= 1
    ebins -= 1

    # Measured effects distribution
    ne_meas = TH1F('ne_meas', 'effects histogram', ebins, earray)
    ne_meas.GetXaxis().SetTitle('Effects')
    ne_meas.GetYaxis().SetTitle('Counts')
    ne_meas.SetStats(0)
    ne_meas.Sumw2()

    # Prepare Combined Weighted Histograms - To be Normalized by Model After Filling
    # Isotropic Weights of Causes - For Calculating Combined Species Efficiency
    eff = TH1F('Eff', 'Non-Normed Combined Efficiency', cbins, carray)
    eff.GetXaxis().SetTitle('Causes')
    eff.GetYaxis().SetTitle('Efficiency')
    eff.SetStats(0)
    eff.Sumw2()

    # Isotropic Weighted Mixing Matrix - For Calculating Combined Species MM
    response = TH2F('MM', 'Weighted Combined Mixing Matrix',
                    cbins, carray, ebins, earray)
    response.GetXaxis().SetTitle('Causes')
    response.GetYaxis().SetTitle('Effects')
    response.SetStats(0)
    response.Sumw2()

    for ci in range(0, cbins):

        # Fill measured effects histogram
        ne_meas.SetBinContent(ci+1, counts[ci])
        if counts_err is None:
            ne_meas.SetBinError(ci+1, np.sqrt(counts[ci]))
        else:
            ne_meas.SetBinError(ci+1, counts_err[ci])

        # Fill response matrix entries
        for ek in range(0, ebins):
            response.SetBinContent(ci+1, ek+1, response_array[ek][ci])
            response.SetBinError(ci+1, ek+1, response_err_array[ek][ci])

        # Fill efficiency histogram from response matrix
        eff.SetBinContent(ci+1, efficiencies[ci])
        eff.SetBinError(ci+1, efficiencies_err[ci])

    # Write measured effects histogram to file
    ne_meas.Write()
    # Write the cause and effect arrays to file
    CARRAY = TH1F('CARRAY','Cause Array', cbins, carray)
    CARRAY.GetXaxis().SetTitle('Causes')
    EARRAY = TH1F('EARRAY','Effect Array', ebins, earray)
    EARRAY.GetXaxis().SetTitle('Effects')
    CARRAY.Write()
    EARRAY.Write()
    # Write efficiencies histogram to file
    eff.Write()
    # Write response matrix to file
    response.Write()

    fout.Write()
    fout.Close()
