
import warnings
import numpy as np
import pandas as pd
from ROOT import gROOT

from .Utils import none_to_empty_list

# Clear ROOT global variables
gROOT.Reset()

# Ignore ROOT fit RuntimeWarning
warnings.filterwarnings("ignore", category=RuntimeWarning, module="ROOT")


class IterativeUnfolder(object):
    """Common base class for iterative unfolder
    """
    def __init__(self, name, max_iter=100, smooth_iter=False, n_c=None,
                 mix_func=None, reg_func=None, ts_func=None, stack=False,
                 verbose=False, **kwargs):

        n_c, mix_func, reg_func, ts_func = none_to_empty_list(n_c,
                                                              mix_func,
                                                              reg_func,
                                                              ts_func)

        self.name = name
        self.max_iter = max_iter
        self.smooth_iter = smooth_iter
        # Stacking info
        self.nstack = len(reg_func)
        # stacking Bin Indices Initializer
        self.stackInd = [[] for i in range(self.nstack)]
        # Mixing, Regularizing, Test Statistic Functions
        self.Mix = mix_func
        self.Rglzr = reg_func
        self.ts_func = ts_func
        # If user doesn't provide prior...
        self.n_c = n_c.copy()
        # Count Iterations
        self.counter = 0

        # Lists to Append Iteration Results
        #  Unfolded Results
        self.n_c_iters = [[] for i in range(self.nstack)]
        self.n_c_final = []
        self.covM = []
        self.statErr = []
        self.sysErr = []
        #  Regularization Results
        self.ParamOut = [[] for i in range(self.nstack)]
        for j in range(self.nstack):
            for i in xrange(0, reg_func[j].nParams):
                self.ParamOut[j].append([])
        self.N_measured_out = []
        self.N_measured_out.append(np.sum(n_c))

        #  Test Statistic Results
        self.TS_out = [[] for i in range(self.nstack)]

        # Set stacking Bin Indices
        self.SetstackInd()

    def SetstackInd(self):
        iind = 0
        for iS in range(self.nstack):
            iaxis = self.Rglzr[iS].xarray.copy()
            self.stackInd[iS] = [range(iind, iind + len(iaxis))]
            iind += len(iaxis)

    def pass_tol(self):
        is_pass_tol = any([self.ts_func[iS].pass_tol()
                           for iS in range(self.nstack)])
        return is_pass_tol

    def IUnfold(self):
        """Perform iterative unfolding

        Parameters
        ----------
        self : self

        Returns
        -------
        unfolding_result : pandas.DataFrame
            DataFrame containing the unfolded result for each iteration.
            Each row in unfolding_result corresponds to an iteration.
        """
        # Create dictionary to store unfolding result at each iteration
        unfolding_result = {}

        # First Mixing
        n_c_update = self.Mix.Smear(self.n_c)
        n_update = np.sum(n_c_update)

        # Add first mixing result to unfolding_result
        unfolding_result[self.counter] = {'unfolded': n_c_update,
                                          'stat_err': self.Mix.getStatErr(),
                                          'sys_err': self.Mix.getMCErr()}

        self.counter += 1

        # Calculate Chi2 of first mix with prior
        for iS in range(self.nstack):
            iind = self.stackInd[iS]
            inc = self.n_c[iind]
            self.n_c_iters[iS].append(inc)
            incu = n_c_update[iind]
            TS_cur, TS_del, TS_prob = self.ts_func[iS].GetStats(incu, inc)
            self.TS_out[iS].append(TS_cur)

        reg_fit = np.zeros(self.n_c.shape)
        # Perform Iterative Unfolding
        while (not self.pass_tol() and self.counter < self.max_iter):
            # Updated unfolded distribution
            n_c = n_c_update.copy()
            # For test statistic purposes.
            # Don't want to compare reg_fit to n_c_update
            n_c_prev = n_c_update.copy()
            # Save n_update for root output
            self.N_measured_out.append(n_update)

            # Regularize unfolded dist in favor of 'smooth physics' prior
            for iS in range(self.nstack):
                iind = self.stackInd[iS][0]
                incu = n_c_update[iind]
                self.n_c_iters[iS].append(incu)
                reg_fit_iS, FitParams = self.Rglzr[iS].Regularize(incu)
                reg_fit[iind] = reg_fit_iS.copy()

                # Save fit parameters for root output
                for j in xrange(0, self.Rglzr[iS].nParams):
                    self.ParamOut[iS][j].append(FitParams[j])

                # Smoothed n_c for next iteration
                if self.smooth_iter:
                    n_c[iind] = reg_fit_iS.copy()

            # Mix w/n_c from previous iter
            n_c_update = self.Mix.Smear(n_c)
            n_update = np.sum(n_c_update)

            # Add mixing result to unfolding_result
            unfolding_result[self.counter] = {
                                        'unfolded': n_c_update,
                                        'stat_err': self.Mix.getStatErr(),
                                        'sys_err': self.Mix.getMCErr()}

            self.counter += 1

            # Calculate TS wrt previous regularization
            for iS in range(self.nstack):
                iind = self.stackInd[iS][0]
                inc = self.n_c[iind]
                incu = n_c_update[iind]
                incp = n_c_prev[iind]
                TS_cur, TS_del, TS_prob = self.ts_func[iS].GetStats(incu, incp)
                self.TS_out[iS].append(TS_cur)

        # Final Iteration
        n_c_final = n_c_update.copy()
        self.n_c_final = n_c_final.copy()
        # Calculate error estimates
        self.covM = self.Mix.getCov()
        self.statErr = self.Mix.getStatErr()
        self.sysErr = self.Mix.getMCErr()
        # Save n_update for root output
        self.N_measured_out.append(n_update)
        # Regularize unfolded dist in favor of 'smooth physics' prior
        for iS in range(self.nstack):
            iind = self.stackInd[iS]
            incu = n_c_final[iind]
            self.n_c_iters[iS].append(incu)
            final_fit, FitParams = self.Rglzr[iS].Regularize(incu)

            for j in xrange(0, self.Rglzr[iS].nParams):
                self.ParamOut[iS][j].append(FitParams[j])

        # Convert unfolding_result dictionary to a pandas DataFrame
        unfolding_result = pd.DataFrame.from_dict(unfolding_result,
                                                  orient='index')

        return unfolding_result
