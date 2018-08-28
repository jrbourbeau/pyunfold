
from __future__ import division, print_function
import numpy as np
from scipy.special import gammaln as lgamma

from .utils import none_to_empty_list


class TestStat(object):
    """Common base class for test statistic methods
    """
    def __init__(self, tol=None, num_causes=None, test_range=None, **kwargs):
        test_range = none_to_empty_list(test_range)
        self.tol = tol
        if num_causes is None:
            raise ValueError('Number of causes (num_causes) must be provided.')
        cause_bin_edges = np.arange(num_causes + 1, dtype=float)
        # Get bin midpoints
        cause_axis = (cause_bin_edges[1:] + cause_bin_edges[:-1]) / 2
        self.cause_axis = cause_axis
        self.ts_range = test_range
        self.has_ts_range = False
        self.ts_bins = self.set_test_range_bins()
        # Initialize Unnatural TS data
        self.stat = np.inf
        self.dof = -1
        self.dofSet = False

    def set_test_range_bins(self):
        bins = [0, -1]
        if self.ts_range != []:
            lTSR = len(self.ts_range)
            err_mess = ("***\n Test stat range can only have two elements. "
                        "This has {}. Exiting...***\n".format(lTSR))
            assert lTSR == 2, err_mess
            xlo = self.ts_range[0]
            xhi = self.ts_range[1]
            err_mess = "***\n Test stat limits reversed. xlo must be < xhi. Exiting...***\n"
            assert xlo < xhi, err_mess

            # Find the bins corresponding to the test range requested
            lobin = np.searchsorted(self.cause_axis, xlo)
            hibin = np.searchsorted(self.cause_axis, xhi)
            bins = [lobin, hibin]
            self.has_ts_range = True
        return bins

    def get_array_range(self, dist1, dist2):
        if self.has_ts_range:
            NR1 = dist1[self.ts_bins[0]:self.ts_bins[1]]
            NR2 = dist2[self.ts_bins[0]:self.ts_bins[1]]
            return NR1, NR2
        else:
            return dist1.copy(), dist2.copy()

    def pass_tol(self):
        """Function testing whether TS < tol
        """
        pass_tol = self.stat < self.tol
        return pass_tol

    def set_dof(self, dof):
        """Set degrees of freedom
        """
        if self.dof == -1:
            self.dof = dof

    def check_lengths(self, dist1, dist2):
        """Test for equal length distributions
        """
        ln1 = len(dist1)
        ln2 = len(dist2)
        err_mess = ("Test Statistic arrays are not equal length. "
                    "{} != {}. Exiting...\n".format(ln1, ln2))
        assert ln1 == ln2, err_mess
        if not self.dofSet:
            self.set_dof(ln1)
            self.dofSet = True

    def calc(self, dist1, dist2):
        """Undefined test statistics calculator
        """
        raise NotImplementedError()


class Chi2(TestStat):
    """Reduced chi-squared test statistic
    """
    def calc(self, dist1, dist2):
        """Calculate the test statistic between two input distributions

        Parameters
        ----------
        dist1 : array_like
            Input distribution.
        dist2 : array_like
            Input distribution.

        Returns
        -------
        stat : float
            Test statistic
        """
        dist1, dist2 = self.get_array_range(dist1, dist2)
        self.check_lengths(dist1, dist2)
        n1 = np.sum(dist1)
        n2 = np.sum(dist2)

        h_sum = dist1 + dist2
        # Don't divide by 0...
        h_sum[(h_sum < 1)] = 1.
        h_dif = n2 * dist1 - n1 * dist2
        h_quot = h_dif * h_dif / h_sum

        stat = np.sum(h_quot)/(n1*n2)/self.dof
        self.stat = stat

        return stat


class BF(TestStat):
    """Bayes factor test statistic

    Notes
    -----
    For details related to the Bayes fator see [1]_.

    References
    ----------
    .. [1] S. Y. BenZvi and B. M. Connolly and C. G. Pfendner and S. Westerhoff.
        "A Bayesian Approach to Comparing Cosmic Ray Energy Spectra".
        *The Astrophysical Journal* 738 (1):82.
        `<https://doi.org/10.1088/0004-637X/738/1/82>`_.
    """
    def calc(self, dist1, dist2):
        """Calculate the test statistic between two input distributions

        Parameters
        ----------
        dist1 : array_like
            Input distribution.
        dist2 : array_like
            Input distribution.

        Returns
        -------
        stat : float
            Test statistic
        """
        dist1, dist2 = self.get_array_range(dist1, dist2)
        self.check_lengths(dist1, dist2)
        lnB = 0
        n1 = np.sum(dist1)
        n2 = np.sum(dist2)
        nFactor = lgamma(n1+n2+2) - lgamma(n1+1) - lgamma(n2+1)

        lnB += nFactor
        for i in range(0, len(dist1)):
            lnB += lgamma(dist1[i]+1) + lgamma(dist2[i]+1) - lgamma(dist1[i]+dist2[i]+2)

        self.stat = lnB

        return lnB


class RMD(TestStat):
    """Maximum relative difference test statistic
    """
    def calc(self, dist1, dist2):
        """Calculate the test statistic between two input distributions

        Parameters
        ----------
        dist1 : array_like
            Input distribution.
        dist2 : array_like
            Input distribution.

        Returns
        -------
        stat : float
            Test statistic
        """
        dist1, dist2 = self.get_array_range(dist1, dist2)
        self.check_lengths(dist1, dist2)

        h_sum = dist1+dist2
        h_sum[(h_sum < 1)] = 1.
        h_dif = np.abs(dist1 - dist2)
        h_quot = h_dif / h_sum
        stat = np.max(h_quot)
        self.stat = stat

        return stat


class KS(TestStat):
    """Kolmogorov-Smirnov (KS) two-sided test statistic
    """
    def calc(self, dist1, dist2):
        """Calculate the test statistic between two input distributions

        Parameters
        ----------
        dist1 : array_like
            Input distribution.
        dist2 : array_like
            Input distribution.

        Returns
        -------
        stat : float
            Test statistic
        """
        dist1, dist2 = self.get_array_range(dist1, dist2)
        self.check_lengths(dist1, dist2)

        n1 = np.sum(dist1)
        n2 = np.sum(dist2)
        cs1 = np.cumsum(dist1)/n1
        cs2 = np.cumsum(dist2)/n2

        len1 = len(dist1)
        self.en = np.sqrt(len1/2)
        stat = np.max(np.abs(cs1-cs2))
        self.stat = stat

        return stat


TEST_STATISTICS = {"chi2": Chi2,
                   "bf": BF,
                   "rmd": RMD,
                   "ks": KS,
                   }


def get_ts(name='ks'):
    """Convenience function for retrieving test statisitc calculators

    Parameters
    ----------
    name : {'ks', 'chi2', 'bf', 'rmd'}
        Name of test statistic.

    Returns
    -------
    ts : TestStat
        Test statistics calculator
    """
    if name in TEST_STATISTICS:
        ts = TEST_STATISTICS[name]
        return ts
    else:
        raise ValueError('Invalid test statistic, {}, entered. Must be '
                         'in {}'.format(name, TEST_STATISTICS.keys()))
