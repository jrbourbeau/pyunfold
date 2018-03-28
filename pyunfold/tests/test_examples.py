
import os
import numpy as np
from ROOT import TFile
from root_numpy import hist2array
import pytest

from .generate_testing_data import save_test_root_file

from pyunfold.Unfold import unfold

here = os.path.abspath(os.path.dirname(__file__))


@pytest.mark.parametrize('response', ['diagonal', 'triangular'])
def test_response_example(tmpdir, response):
    # Generate testing input ROOT file
    root_file = os.path.join(str(tmpdir), 'test_file.root')
    save_test_root_file(outfile=root_file, response=response)

    # Run unfolding
    config_file = os.path.abspath(os.path.join(here, '../', 'config.cfg'))
    df_unfolding_iter = unfold(config_name=config_file,
                               EffDist=None,
                               priors='Jeffreys',
                               input_file=root_file,
                               ts='ks',
                               ts_stopping=0.1)
    unfolded_counts = df_unfolding_iter.iloc[-1]['n_c']

    # Extract input counts distribution from input ROOT file
    f = TFile.Open(root_file)
    hist = f.Get('bin0/ne_meas')
    measured_counts = hist2array(hist)

    # Given diagonal response matrix, unfolded counts should be same as measured counts
    if response == 'diagonal':
        assert np.allclose(measured_counts, unfolded_counts)
    else:
        assert not np.allclose(measured_counts, unfolded_counts)
