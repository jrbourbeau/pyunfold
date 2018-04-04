
from __future__ import division, print_function
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import ROOT
from ROOT import TF1, TH1F
from ROOT import gROOT

gROOT.Reset()
# Turn Off TCanvas Warning at DATA.Fit
ROOT.gErrorIgnoreLevel = ROOT.kWarning


class Regularizer(object):
    """Regularizer Class
    """
    def __init__(self, name, FitFunc, Range, InitialParams, ParamLo, ParamHi,
                 ParamNames, xarray, xedges, verbose=False, plot=False,
                 **kwargs):

        # Prepare string fnc definition for ROOT TF1 declaration
        funcstring = str(FitFunc)
        funcstring = funcstring.strip('[]')
        funcstring = funcstring.strip("''")

        '''From Initialization Inputs'''
        self.name = name
        self.verbose = verbose
        self.plot = plot
        # Number of data bins and xarray and edges of bins
        self.nbins = len(xarray)
        self.xarray = xarray.copy()
        self.xedges = xedges.copy()
        self.Range = self.SetFitRange(Range)
        # Number of fit params, names, init vals, & limits
        self.nParams = len(InitialParams)
        self.ParamNames = ParamNames
        self.Params = np.asarray(InitialParams, dtype=float)
        self.ParamLimitsLo = np.asarray(ParamLo, dtype=float)
        self.ParamLimitsHi = np.asarray(ParamHi, dtype=float)
        self.ParLimProvided = False
        # Chi2 of fit and degrees of freedom
        self.chi2 = -1
        self.dof = 0
        # Fit function string
        self.FitFunc = funcstring
        self.FIT = TF1("FIT", self.FitFunc, self.Range[0], self.Range[1])

        # Ensure that Fit Function is Initialized Properly
        self.TestLengths(ParamNames, InitialParams, "Parameter Names & Initial Value")
        self.TestFuncString(funcstring)
        # Set the Parameter limits if provided
        self.SetParLimits()

        # self.PrintInitMessage()

    # Initialization Message Print Out
    def PrintInitMessage(self):
        stringFunc = self.FitFunc
        for i in xrange(0, self.nParams):
            stringFunc = stringFunc.replace(str(i), self.ParamNames[i])
        stringFunc = stringFunc.replace('[', "")
        stringFunc = stringFunc.replace(']', "")
        print("\nRegularizing {}-parameter function initialized to "
              "form: {}".format(self.nParams, stringFunc))
        print("Can only support fit functions with up to 10 parameters.\n")

    def SetFitRange(self, Range):

        if len(Range) == 0:
            Range = np.zeros(2)
            Range[0] = self.xarray[0]
            Range[1] = self.xarray[-1]
        else:
            lR = len(Range)
            err_mess = ("\n*** Fit range can only have two elements. "
                        "This has {}. Exiting...***\n".format(lR))
            assert lR == 2, err_mess
            xlo = Range[0]
            xhi = Range[1]
            err_mess = "\n*** Fit range limits reversed. xlo must be < xhi. Exiting...***\n"
            assert xlo < xhi, err_mess

        return Range

    def SetParLimits(self):
        # Set Par Limits if provided
        if not len(self.ParamLimitsLo) == 0 and not len(self.ParamLimitsHi) == 0:
            PLo = self.ParamLimitsLo.copy()
            PHi = self.ParamLimitsHi.copy()
            self.TestLengths(PLo, PHi, "Parameter Limit")
            self.TestLengths(PLo, self.Params, "Parameter Limit & Initial Value")
            for i in xrange(0, self.nParams):
                err_mess = ("\n*** Regularizer param {} limits reversed. "
                            "xlo must be < xhi. Exiting...***\n".format(i))
                assert PLo[i] < PHi[i], err_mess
                self.ParLimProvided = True

    # Test ROOT String Func Number of Params
    def TestFuncString(self, string):
        num_nums = sum(c.isdigit() for c in string)

        err_mess = ("\nRegularizer func has {} parameters != {} requested "
                    "in InitialParams keyword. Exiting...\n".format(num_nums, self.nParams))
        assert num_nums == self.nParams, err_mess

        func_pars = np.asarray([float(c) for c in string if c.isdigit()])

        err_mess = ("\nRegularizer func definition is funky. "
                    "Perhaps you counted twice or skipped an integer: "
                    "{}. Exiting...\n".format(string))
        assert np.array_equal(func_pars, np.linspace(0, num_nums-1, num_nums)), err_mess

    # Test for Equal Length Arrays
    def TestLengths(self, N1, N2, message):
        ln1 = len(N1)
        ln2 = len(N2)
        err_mess = ("\nRegularizer {} arrays are not equal length. {} != {}. "
                    "Exiting...\n".format(message, ln1, ln2))
        assert ln1 == ln2, err_mess

    # Get Reduced Chi2 of Fit
    def GetRedChi2(self):
        err_mess = "\nFit not yet performed, dof = {}".format(self.dof)
        assert not self.dof == 0, err_mess

        return self.chi2 / self.dof

    # Evaulate fit at points given by xarray
    def FitEval(self, xarray):
        f_eval = np.zeros(len(xarray))
        for j in xrange(0, len(xarray)):
            f_eval[j] = self.FIT.Eval(xarray[j])
        return f_eval

    # Regularization procedure
    def Regularize(self, ydata, yerr=None):

        if yerr is None:
            yerr = []

        # Local variable for cleanliness
        nPar = self.nParams

        self.TestLengths(self.xarray, ydata, "Cause X-axis & Unfolded Data")
        # ROOT object to fit
        DATA = TH1F("data", "data", self.nbins, self.xedges)
        for i in xrange(0, self.nbins):
            DATA.SetBinContent(i+1, ydata[i])
            if len(yerr) > 0:
                DATA.SetBinError(i+1, yerr[i])
            else:
                DATA.SetBinError(i+1, np.sqrt(ydata[i]))

        # Prepare fit parameters based on user InitialParams
        # or previous fit if used iteratively
        Pars = np.array(self.Params, dtype=float)
        self.FIT.SetParameters(Pars)

        # Set Par Limits if provided
        if self.ParLimProvided:
            PLo = self.ParamLimitsLo.copy()
            PHi = self.ParamLimitsHi.copy()
            for i in xrange(0, self.nParams):
                self.FIT.SetParLimits(i, PLo[i], PHi[i])

        Fopts = "R"
        if not self.verbose:
            # Quiet mode
            Fopts += "Q"
        DATA.Fit(self.FIT, Fopts)

        # Eval fit at xarray
        f_eval = self.FitEval(self.xarray)

        # Store fit parameters
        for i in xrange(0, nPar):
            self.Params[i] = self.FIT.GetParameter(i)

        # Store chi2 and dof of fit
        self.chi2 = self.FIT.GetChisquare()
        self.dof = self.FIT.GetNDF()

        # Print Results Nicely :)
        if self.verbose:
            print("=====================")
            print("Fit Parameter Results")
            for i in xrange(0, nPar):
                print(self.ParamNames[i], " = ", self.Params[i])
            print("Red X2: = ", self.GetRedChi2())
            print("=====================")

        # Plot Comparison Nicely :)
        if self.plot:
            self.Plot(ydata, f_eval)

        # Return fit eval array and fit parameters
        return f_eval, self.Params

    # Plot data and fit
    def Plot(self, data, fitdata):
            fig = plt.figure(figsize=(8, 7))
            ax = fig.add_subplot(111)
            mpl.rc("font", family="serif", size=14)
            ax.plot(self.xarray, data, label="data", color='k', ls='--')
            ax.plot(self.xarray, fitdata, label="fit", color='r', ls=':')
            ax.set_xscale("log")
            ax.set_xlabel("X")
            ax.set_xlim([self.xarray[0], self.xarray[-1]])
            ax.set_yscale("log")
            ax.set_ylim([np.min(fitdata), np.max(fitdata)])
            ax.set_ylabel("Y")
            ax.set_title("Compare Data to Fit")
            handles, labels = ax.get_legend_handles_labels()
            leg = ax.legend(handles, labels, fontsize=10, loc="best")
            leg.get_frame().set_linewidth(0)
            plt.show()
