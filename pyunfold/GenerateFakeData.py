#!/usr/bin/env python
"""
   Routine to generate fake data from
   a user provided response matrix.
   By default looks for config.cfg
   file for response file, zenith
   bin and outfile definitions.
"""

import numpy as np
import argparse
import sys
import os

from .Mix import Mixer
from .Utils import *
from .LoadStats import *
from .IterUnfold import *

from .Plotter import *
from .RootReader import *


# Global program options
args = None

def main():
    global args
    p = argparse.ArgumentParser(description="Script to Generate Fake Reconstructed Energy Distribution")
    p.add_argument("-c", "--config", dest="config_name", default="config.cfg", help="Configuration file")
    p.add_argument("-o", "--outfile", dest="outfile", default="FakeData.root", required=False, help="Output File")
    p.add_argument("-z", "--zbin", dest="zbin", default=0, type=int, required=False, help="Zenith Bin")
    p.add_argument("-p", "--plot", dest="plot_local", action="store_true", help="Friendly plotter")
    args = p.parse_args()
    usage = "%prog -c config_file -o outfile -z zbin -p"

    if (args.config_name is None): raise Exception, "Need a config file!"

    config = ConfigFM(args.config_name)

    # Get Output file name and zbin
    OutFileName = args.outfile
    zbin = args.zbin
    # Define Out hist names (input for unfolding, thus the 'data' header)
    dataHeader = "data"
    NC_outname =  config.get(dataHeader,"nc_thrown",default="",cast=str)
    NE_outname =  config.get(dataHeader,"ne_meas",default="",cast=str)

    # Get MCInput
    mcHeader = "mcinput"
    StatsFile = config.get(mcHeader,"stats_file",default="",cast=str)
    Eff_hist = config.get(mcHeader,"eff_hist",default="",cast=str)
    MM_hist = config.get(mcHeader,"mm_hist",default="",cast=str)

    # If zbin is btw 0 and 2, then we're using the likelihood energy estimator!
    if (zbin>-1 and zbin < 3):
        Eff_hist +="_%s"%zbin
        MM_hist +="_%s"%zbin
        NC_outname += "_%i"%zbin
        NE_outname += "_%i"%zbin

    # Load Effective Area (Aeff) and Migration Matrix ( P(E|C) )
    MCStats = LoadMCTables(StatsFile,RespMatrixName=MM_hist,EffName=Eff_hist)
    Caxis, Cedges = MCStats.GetCauseAxis()
    Eaxis, Eedges = MCStats.GetEffectAxis()
    pec = MCStats.GetPEC()

    # Plotting Labels
    #  Effects
    Exlab = r'%s'%MCStats.GetEffectLabel()
    Eylab = r'$\mathscr{N}$ (%s)'%Exlab
    #  Causes
    Cxlab = r'%s'%MCStats.GetCauseLabel()
    Cylab = r'$\mathscr{N}$ (%s)'%Cxlab

    # Number of effects and causes
    pec_shape = pec.shape
    e_bins = pec_shape[0]
    c_bins = pec_shape[1]

    # Observed effect frequency dist
    Nevents1 = 1e15
    idx1 = 2.6
    p1 = PowerLaw("PL1",Nevents1,idx1,[Cedges[0],Cedges[len(Cedges)-1]])
    N1 = p1.Fill(Cedges,"analytic")
    Nevents2 = 1e15
    idx2 = 2.6
    p2 = PowerLaw("PL2",Nevents2,idx2,[Cedges[0],Cedges[len(Cedges)-1]])
    N2 = p2.Fill(Cedges,"analytic")
    NSim = N1+N2
    #NSim *= np.exp(-Caxis/1e4)
    NSim = np.round(NSim)

    n_thrown = np.sum(NSim)

    ''' No Touchy! Saves fake data appropriately '''
    # Apply Detector Response to Thrown Spectrum
    n_ef = np.dot(pec,NSim)
    # Apply Noise
    n_ef = np.sqrt(n_ef)*np.random.randn(e_bins)+n_ef
    # Make Count Data
    n_ef = np.round(n_ef)
    n_obs = np.sum(n_ef)

    print "\n\n====================================================\n"
    print "Total number of simulated events:\t%e"%n_thrown
    print "Total number of observed events:\t%e"%n_obs
    print "\n====================================================\n\n"

    # Write to ROOT Output
    #  Thrown Causes
    RFile = ROOT.TFile(OutFileName, "UPDATE")
    RFile.cd()
    NTHROWN = TH1F("NC","Thrown Fake Cause Distribution",len(Caxis),Cedges)
    NTHROWN.GetXaxis().SetTitle(Cxlab.strip("$"))
    NTHROWN.GetYaxis().SetTitle(Cylab.strip("$"))
    NTHROWN.SetStats(0)
    for j in xrange(0,len(Cedges)-1):
        binval = NSim[j]
        NTHROWN.SetBinContent(j+1,binval)
        NTHROWN.SetBinError(j+1,np.sqrt(binval))
    #  Observed Effects
    NOBS = TH1F("NE","Observed Fake Effects Distribution",len(Eaxis),Eedges)
    NOBS.GetXaxis().SetTitle(Exlab.strip("$"))
    NOBS.GetYaxis().SetTitle(Eylab.strip("$"))
    NOBS.SetStats(0)
    for j in xrange(0,len(Eedges)-1):
        binval = n_ef[j]
        NOBS.SetBinContent(j+1,binval)
        NOBS.SetBinError(j+1,np.sqrt(binval))

    # Write to file
    NTHROWN.Write(NC_outname)
    NOBS.Write(NE_outname)
    RFile.Close()

    print "\n\n========== Writing ROOT File ==============="
    print "\nFake Data ROOT file saved as %s\n"%OutFileName
    print "============= Done Writing =================\n\n"

    if (args.plot_local):
        print("*********\n\nIf you can't see your distributions, "\
              "check to see if log plot flags are set in script.\n\n*********")
        #  Plot Cause and Effect Distributions
        fig = plt.figure(figsize=(17,6))
        ax = fig.add_subplot(121)
        mpl.rc("font", family="serif", size=14)

        title = "Integral Cause Distribution"
        xlim = [Caxis[0],Caxis[-1]]
        ylim = [np.min(NSim),np.max(NSim)*1.1]
        log = 'x,y'
        ax.step(Caxis, NSim, color='b',ls='-')
        if ('x' in log):
            ax.set_xscale('log', nonposx='clip')
        if ('y' in log):
            ax.set_yscale("log", nonposx='clip')
        ax.set_xlabel(r'%s'%Cxlab)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_ylabel(r'%s'%Cylab)
        ax.set_title(title)
        ax.grid(True)
        fig.subplots_adjust(left=0.125, right=0.85)

        ax = fig.add_subplot(122)
        title = "Integral Effects Distribution"
        xlim = [Eaxis[0],Eaxis[-1]]
        ylim = [np.min(n_ef[(n_ef>0)]),np.max(n_ef)*1.1]
        log = 'x,y'
        ax.step(Eaxis, n_ef, color='b',ls='-')
        if ('x' in log):
            ax.set_xscale('log', nonposx='clip')
        if ('y' in log):
            ax.set_yscale("log", nonposx='clip')
        ax.set_xlabel(Exlab)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_ylabel(Eylab)
        ax.set_title(title)
        ax.grid(True)
        fig.subplots_adjust(left=0.125, right=0.85)
        plt.show()



if __name__ == "__main__":
    main()
