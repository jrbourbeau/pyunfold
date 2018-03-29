
import os
import tempfile
import shutil
import numpy as np
import pandas as pd

from .LoadStats import MCTables
from .IterUnfold import IterativeUnfolder
from .Mix import Mixer
from .RootReader import get_labels, get1d
from .Utils import (save_input_to_root_file, ConfigFM, get_ts, DataDist,
                    UserPrior, Regularizer)


def iterative_unfold(counts, counts_err, response, response_err, efficiencies,
                     efficiencies_err, priors='Jeffreys', ts='ks',
                     ts_stopping=0.01, max_iter=100, return_iterations=False):
    """Performs iterative Bayesian unfolding

    Parameters
    ----------
    counts : array_like
        Input observed data distribution.
    counts_err : array_like
        Uncertainties associated with the input observed data distribution.
        Must be the same shape as counts.
    response : array_like
        Response matrix.
    response_err : array_like
        Response matrix errors.
    efficiencies : array_like
        Detection efficiencies for the observed data distribution.
    efficiencies_err : array_like
        Uncertainty in detection efficiencies.
    outfile : str
        Path where output ROOT file will be saved.
    priors : str or array_like, optional
        Prior distribution to use in unfolding. If 'Jeffreys', then the
        Jeffreys (flat) prior is used. Otherwise, must be array_like with
        same shape as counts (default is 'Jeffreys').
    ts : {'ks', 'chi2', 'pf', 'rmd'}
        Name of test statistic to use for stopping condition (default is 'ks').
    ts_stopping : float, optional
        Test statistic stopping condition. At each unfolding iteration, the
        test statistic is computed between the current and previous iteration.
        Once the test statistic drops below ts_stopping, the unfolding
        procedure is stopped (default is 0.01).
    max_iter : int, optional
        Maximum number of iterations to allow (default is 100).
    return_iterations : bool, optional
        Whether to return unfolded distributions for each iteration
        (default is False).

    Returns
    -------
    unfolded_result : dict
        Returned if return_iterations is False (default). Final unfolded
        distribution and associated uncertainties.
    unfolding_results : pandas.DataFrame
        Returned if return_iterations is True. DataFrame containing the
        unfolded distribution and associated uncertainties at each iteration.
    """
    # Save counts, response, efficiencies, etc. to a ROOT file
    temp_dir_path = tempfile.mkdtemp()
    root_file = os.path.join(temp_dir_path, 'pyunfold_root_file.root')
    save_input_to_root_file(counts,
                            counts_err,
                            response,
                            response_err,
                            efficiencies,
                            efficiencies_err,
                            outfile=root_file)

    here = os.path.abspath(os.path.dirname(__file__))
    config_file = os.path.join(here, 'config.cfg')

    unfolding_results = unfold(config_name=config_file,
                               EffDist=None,
                               priors=priors,
                               input_file=root_file,
                               ts=ts,
                               ts_stopping=ts_stopping,
                               max_iter=max_iter)
    shutil.rmtree(temp_dir_path)

    if return_iterations:
        return unfolding_results
    else:
        unfolded_result = dict(unfolding_results.iloc[-1])
        return unfolded_result


def unfold(config_name=None, EffDist=None, priors='Jeffreys', input_file=None,
           ts='ks', ts_stopping=0.01, max_iter=100, **kwargs):

    if config_name is None:
        raise ValueError('config_name must be provided')

    assert input_file is not None

    # Get the Configuration Parameters from the Config File
    config = ConfigFM(config_name)

    # Input Data ROOT File Name
    dataHeader = 'data'
    InputFile = input_file
    NE_meas_name = config.get(dataHeader, 'ne_meas', default='', cast=str)
    isMC = config.get_boolean(dataHeader, 'ismc', default=False)

    # Analysis Bin
    binHeader = 'analysisbins'
    binnumberStr = config.get(binHeader, 'bin', default=0, cast=str)
    stackFlag = config.get_boolean(binHeader, 'stack', default=False)
    bin_list = ['bin'+val.replace(' ', '') for val in binnumberStr.split(',')]
    nStack = len(bin_list)

    unfbinname = 'bin0'

    # Unfolder Options
    unfoldHeader = 'unfolder'
    UnfolderName = config.get(unfoldHeader, 'unfolder_name',
                              default='Unfolder', cast=str)
    # UnfMaxIter = config.get(unfoldHeader, 'max_iter', default=100, cast=int)
    UnfMaxIter = max_iter
    UnfSmoothIter = config.get_boolean(unfoldHeader, 'smooth_with_reg',
                                       default=False)
    UnfVerbFlag = config.get_boolean(unfoldHeader, 'verbose', default=False)

    # Mixer Name and Error Propagation Type
    # Options: ACM, DCM
    mixHeader = 'mixer'
    MixName = config.get(mixHeader, 'mix_name', default='', cast=str)
    CovError = config.get(mixHeader, 'error_type', default='', cast=str)

    # Test Statistic - Stat Function & Options
    # Options: chi2, rmd, pf, ks
    tsname = ts
    tsTol = ts_stopping
    tsRange = [0, 1e2]
    tsVerbFlag = False

    # Regularization Function, Initial Parameters, & Options
    regHeader = 'regularization'
    RegFunc = config.get(regHeader, 'reg_func', default='', cast=str)
    #  Param Names
    ConfigParamNames = config.get(regHeader, 'param_names', default='',
                                  cast=str)
    ParamNames = [x.strip() for x in ConfigParamNames.split(',')]
    #  Initial Parameters
    IPars = config.get(regHeader, 'param_init', default='', cast=str)
    InitParams = [float(val) for val in IPars.split(',')]
    #  Limits
    PLow = config.get(regHeader, 'param_lim_lo', default='', cast=str)
    PLimLo = [float(val) for val in PLow.split(',')]
    PHigh = config.get(regHeader, 'param_lim_hi', default='', cast=str)
    PLimHi = [float(val) for val in PHigh.split(',')]
    RegRangeStr = config.get(regHeader, 'reg_range', cast=str)
    RegRange = [float(val) for val in RegRangeStr.split(',')]
    #  Options
    RegPFlag = config.get_boolean(regHeader, 'plot', default=False)
    RegVFlag = config.get_boolean(regHeader, 'verbose', default=False)

    # Get MCInput
    mcHeader = 'mcinput'
    # StatsFile = config.get(mcHeader, 'stats_file', default='', cast=str)
    StatsFile = input_file
    Eff_hist_name = config.get(mcHeader, 'eff_hist', default='', cast=str)
    MM_hist_name = config.get(mcHeader, 'mm_hist', default='', cast=str)

    # Setup the Observed and MC Data Arrays
    # Load MC Stats (NCmc), Cause Efficiency (Eff) and Migration Matrix P(E|C)
    MCStats = MCTables(StatsFile,
                       BinName=bin_list,
                       RespMatrixName=MM_hist_name,
                       EffName=Eff_hist_name,
                       Stack=stackFlag)
    Caxis = []
    Cedges = []
    cutList = []
    for index in range(nStack):
        axis, edge = MCStats.GetCauseAxis(index)
        Caxis.append(axis)
        Cedges.append(edge)
    Eaxis, Eedges = MCStats.GetEffectAxis()
    # Effect and Cause X and Y Labels from Respective Histograms
    Cxlab, Cylab, Ctitle = get_labels(StatsFile, Eff_hist_name, bin_list[0],
                                      verbose=False)

    # Load the Observed Data (n_eff), define total observed events (n_obs)
    # Get from ROOT input file if requested
    if EffDist is None:
        Exlab, Eylab, Etitle = get_labels(InputFile, NE_meas_name, unfbinname,
                                          verbose=False)
        Eaxis, Eedges, n_eff, n_eff_err = get1d(InputFile, NE_meas_name,
                                                unfbinname)
        EffDist = DataDist(Etitle, data=n_eff, error=n_eff_err,
                           axis=Eaxis, edges=Eedges, xlabel=Exlab,
                           ylabel=Eylab, units='')
    Exlab = EffDist.xlab
    Eylab = EffDist.ylab
    Etitle = EffDist.name
    n_eff = EffDist.getData()
    n_eff_err = EffDist.getError()
    n_obs = np.sum(n_eff)

    # Initial best guess (0th prior) expected prob dist
    # default: Jeffrey's Prior
    if isinstance(priors, (list, tuple, np.ndarray, pd.Series)):
        n_c = np.asarray(priors)
    elif priors == 'Jeffreys':
        n_c = UserPrior(['Jeffreys'], Caxis, n_obs)
        n_c = n_c / np.sum(n_c)
    else:
        raise TypeError('priors must be a array_like, '
                        'but got {}'.format(type(priors)))
    np.testing.assert_allclose(np.sum(n_c), 1)

    # Setup the Tools Used in Unfolding
    # Prepare Regularizer
    Rglzr = [Regularizer('REG', FitFunc=[RegFunc], Range=RegRange,
                         InitialParams=InitParams, ParamLo=PLimLo,
                         ParamHi=PLimHi, ParamNames=ParamNames,
                         xarray=Caxis[i], xedges=Cedges[i],
                         verbose=RegVFlag, plot=RegPFlag)
             for i in range(nStack)]
    # Prepare Test Statistic-er
    tsMeth = get_ts(tsname)
    tsFunc = [tsMeth(tsname,
                     tol=tsTol,
                     Xaxis=Caxis[i],
                     TestRange=tsRange,
                     verbose=tsVerbFlag)
              for i in range(nStack)]

    # Prepare Mixer
    mixer = Mixer(MixName,
                  ErrorType=CovError,
                  MCTables=MCStats,
                  EffectsDist=EffDist)

    # Unfolder!!!
    if stackFlag:
        UnfolderName += '_'+''.join(bin_list)

    Unfolder = IterativeUnfolder(UnfolderName,
                                 max_iter=UnfMaxIter,
                                 smooth_iter=UnfSmoothIter,
                                 n_c=n_c,
                                 mix_func=mixer,
                                 reg_func=Rglzr,
                                 ts_func=tsFunc,
                                 stack=stackFlag,
                                 verbose=UnfVerbFlag)
    # Iterate the Unfolder
    unfolding_result = Unfolder.IUnfold()

    return unfolding_result
