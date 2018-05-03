
from __future__ import division, print_function
import numpy as np
from scipy.special import gammainc as gammaq
from scipy.special import gammaln as lgamma
from scipy.stats import kstwobign

from .utils import none_to_empty_list


class TestStat(object):
    """Common base class for test statistic methods
    """
    def __init__(self, name="TestStat", tol=None, num_causes=None,
                 TestRange=None, verbose=False, **kwargs):
        TestRange = none_to_empty_list(TestRange)
        # TS Name
        self.name = name
        # User Defined TS Tolerance
        self.tol = tol
        if num_causes is None:
            raise ValueError('Number of causes (num_causes) must be provided.')
        cause_bin_edges = np.arange(num_causes + 1, dtype=float)
        # Get bin midpoints
        cause_axis = (cause_bin_edges[1:] + cause_bin_edges[:-1]) / 2
        self.Xaxis = cause_axis
        self.TSRange = TestRange
        self.IsRangeTS = False
        self.TSbins = self.SetTestRangeBins()
        # Verbose Functionality
        self.verbose = verbose
        self.printProbMessage = True
        self.printStatsHeader = True
        if verbose:
            self.PrintName()
        # Initialize Unnatural TS data
        self.stat = np.inf
        self.delstat = 0
        self.prob = -1
        self.dof = -1
        self.dofSet = False

    def SetTestRangeBins(self):
        bins = [0, -1]
        if (not self.TSRange == []):
            lTSR = len(self.TSRange)
            err_mess = ("***\n Test stat range can only have two elements. "
                        "This has {}. Exiting...***\n".format(lTSR))
            assert lTSR == 2, err_mess
            xlo = self.TSRange[0]
            xhi = self.TSRange[1]
            err_mess = "***\n Test stat limits reversed. xlo must be < xhi. Exiting...***\n"
            assert xlo < xhi, err_mess

            # Find the bins corresponding to the test range requested
            lobin = np.searchsorted(self.Xaxis, xlo)
            hibin = np.searchsorted(self.Xaxis, xhi)
            bins = [lobin, hibin]
            self.IsRangeTS = True
        return bins

    def GetArrayRange(self, N1, N2):
        if self.IsRangeTS:
            NR1 = N1[self.TSbins[0]:self.TSbins[1]]
            NR2 = N2[self.TSbins[0]:self.TSbins[1]]
            return NR1, NR2
        else:
            return N1.copy(), N2.copy()

    def pass_tol(self):
        """Function testing whether TS < tol
        """
        pass_tol = self.stat < self.tol
        return pass_tol

    # Set Degrees of Freedom
    def SetDOF(self, dof):
        if (self.dof == -1):
            self.dof = dof

    # Set TS and delTS
    def SetStat(self, stat):
        self.delstat = stat - self.stat
        self.stat = stat

    # Test for Equal Length Distributions
    def TestLengths(self, N1, N2):
        ln1 = len(N1)
        ln2 = len(N2)
        err_mess = ("Test Statistic arrays are not equal length. "
                    "{} != {}. Exiting...\n".format(ln1, ln2))
        assert ln1 == ln2, err_mess
        if not self.dofSet:
            self.SetDOF(ln1)
            self.dofSet = True

    # Calculate the TS
    def TSCalc(self, N1, N2):
        """Undefined Test Statistics Calculator
        """
        raise NotImplementedError()

    # Calculate the TS Probability Function
    def Prob(self):
        if self.verbose and self.printProbMessage:
            print("***No probability function defined for this method, {}. "
                  "It's ok...***".format(self.name))
            self.printProbMessage = False

    def PrintName(self):
        print("\nTest Statistic Method: ", self.__doc__)
        print("Test statistic only valid in range: {} {}\n".format(self.TSRange[0],
                                                                   self.TSRange[1]))

    def PrintStats(self):
        output = "\t\t{:.04f}\t{:.04f}\t{:.03f}\n".format(self.stat, self.delstat, self.prob)
        if self.printStatsHeader:
            print("\t\tStat\tdelStat\tProb")
            self.printStatsHeader = False
        print(output)

    # Calculate and return TS Data
    def GetStats(self, N1, N2):
        # Calculate the Test Statistic
        self.TSCalc(N1, N2)
        # Calculate the Probability of TS
        self.Prob()
        # Print if default
        if self.verbose:
            self.PrintStats()
        # Return TS Data
        return self.stat, self.delstat, self.prob


class Chi2(TestStat):
    """Reduced Chi2 Test Statistic
    """
    def TSCalc(self, N1, N2):
        N1, N2 = self.GetArrayRange(N1, N2)
        self.TestLengths(N1, N2)
        n1 = np.sum(N1)
        n2 = np.sum(N2)

        h_sum = N1 + N2
        # Don't divide by 0...
        h_sum[(h_sum < 1)] = 1.
        h_dif = n2 * N1 - n1 * N2
        h_quot = h_dif * h_dif / h_sum

        stat = np.sum(h_quot)/(n1*n2)/self.dof

        self.SetStat(stat)

    def Prob(self):
        # Chi2 Probability Function
        self.prob = gammaq(0.5*self.dof, 0.5*self.stat)


# Bayes Factor Test - comparing two binned distributions
# Method B - Pfendner et al.
# Recall that lgamma(x) = log ( gamma(x) )
class PF(TestStat):
    """Bayes Factor Test Statistic
    """
    def TSCalc(self, N1, N2):
        N1, N2 = self.GetArrayRange(N1, N2)
        self.TestLengths(N1, N2)
        lnB = 0
        n1 = np.sum(N1)
        n2 = np.sum(N2)
        nFactor = lgamma(n1+n2+2) - lgamma(n1+1) - lgamma(n2+1)

        lnB += nFactor
        for i in range(0, len(N1)):
            lnB += lgamma(N1[i]+1) + lgamma(N2[i]+1) - lgamma(N1[i]+N2[i]+2)

        self.SetStat(lnB)


# Relative Difference Test - comparing two binned distributions
# Taking the maximal relative difference as TS
class RMD(TestStat):
    'Max Relative Difference Test Statistic'
    def TSCalc(self, N1, N2):
        N1, N2 = self.GetArrayRange(N1, N2)
        self.TestLengths(N1, N2)

        h_sum = N1+N2
        h_sum[(h_sum < 1)] = 1.
        h_dif = np.abs(N1 - N2)
        h_quot = h_dif / h_sum
        stat = np.max(h_quot)

        self.SetStat(stat)


# Kuiper's Statistic - comparing two binned distributions
# Taken from https://github.com/scipy/scipy/blob/v0.14.0/scipy/stats/stats.py#L3809
# and
# http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.stats.kstwobign.html
class KS(TestStat):
    """KS Test Statistic
    """
    def TSCalc(self, N1, N2):
        N1, N2 = self.GetArrayRange(N1, N2)
        self.TestLengths(N1, N2)

        n1 = np.sum(N1)
        n2 = np.sum(N2)
        cs1 = np.cumsum(N1)/n1
        cs2 = np.cumsum(N2)/n2

        len1 = len(N1)
        self.en = np.sqrt(len1/2)
        d = np.max(np.abs(cs1-cs2))
        self.d = d
        self.SetStat(d)

    def Prob(self):
        try:
            prob = kstwobign.sf((self.en + .12 + .11 / self.en) * self.d)
        except Exception:
            prob = 1.0

        # KS Probability Function
        self.prob = prob


TEST_STATISTICS = {"chi2": Chi2,
                   "pf": PF,
                   "rmd": RMD,
                   "ks": KS,
                   }


def get_ts(name='ks'):
    """Convenience function for retrieving TestStat object

    Parameters
    ----------
    name : {'ks', 'chi2', 'pf', 'rmd'}
        Name of test statistic.

    Returns
    -------
    ts : TestStat
        Test statistics object
    """
    if name in TEST_STATISTICS:
        ts = TEST_STATISTICS[name]
        return ts
    else:
        raise ValueError('Invalid test statisitc, {}, entered. Must be '
                         'in {}'.format(name, TEST_STATISTICS.keys()))
