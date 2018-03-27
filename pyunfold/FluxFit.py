#!/usr/bin/env python
"""
   Script to transform fit from counts
   to flux, including a fit to config.cfg
   defined function.
   Requires config script via command line
   option -c.
"""

import ROOT
from ROOT import TCanvas, TFile, TH1D, TH1F, TF1, TGraph, TGraphErrors
from ROOT import gROOT, gPad, gRandom, gSystem, TDirectory
gROOT.Reset()

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

import Utils
import RootReader as rr

import numpy as np
import argparse
import sys
import os
import re


# Global program options
args = None

colors = ["red", "blue", "green"]
styles = ["-", "--", "-."]


def getZenithLimits(cuts):

    # Grab zenith angle cut
    if "theta" not in cuts:
        raise ValueError("\nNeed to define some zenith angle (theta) cut.\nPlease, fix your mistake. Exiting...\n")

    zcutstring = re.finditer("theta([<>=]*)([0-9.e-]*)", cuts)
    #zcutstring = re.finditer("rec.zenithAngle([<>=]*)([0-9.e-]*)", cuts)
    zcuts = []
    for izcut in zcutstring:
        comp, val = izcut.groups()
        zcuts.append(val)
    zcuts = np.asarray(zcuts,dtype=float)

    if len(zcuts)>2:
        print("\n=========================\n")
        raise ValueError("Too many zenith cuts! There can only be one (or two)!\nPlease, fix your mistake. Exiting...\n")

    zlo = np.min(zcuts)
    zhi = np.max(zcuts)

    if len(zcuts) == 1:
        zlo = 0

    return zlo, zhi


# Main Flux Fit
def main():
    global args
    p = argparse.ArgumentParser(description="Script to Calculate a Differential Flux from Unfolded Distribution")
    p.add_argument("-c", "--config", dest="config_name", required=True, help="Configuration file")
    p.add_argument("-p", "--plot", dest="plot_local", action="store_true", help="Friendly plotter")
    p.add_argument("-w", "--write", dest="write_fit", action="store_true", help="Write to data ROOT file")
    p.add_argument("-s", "--scale", dest="scale", type=float, default=0., help="Energy scaling applied to flux plot")

    try:
        args = p.parse_args()
    except:
        p.print_help()
        sys.exit(0)

    if (args.config_name is None): raise Exception, "Need a config file!"

    #### Get the Configuration Parameters from the Config File ####
    config = Utils.ConfigFM(args.config_name)

    # Raw (reco) Data ROOT File Name
    dataHeader = "data"
    RawFile = config.get(dataHeader,"inputfile",default="",cast=str)
    TimeName = config.get(dataHeader,"obstime",default="",cast=str)

    # Input Data ROOT File Name (Is output from Unfolding...)
    inHeader = "output"
    InputFile = config.get(inHeader,"outfile",default="",cast=str)
    NC_meas = config.get(inHeader,"nc_meas",default="",cast=str)

    # Get MCInput - Needed for Zenith Cuts
    mcHeader = "mcinput"
    StatsFile = config.get(mcHeader,"stats_file",default="",cast=str)

    # Regularization Function, Initial Parameters, & Options
    regHeader = "flux"
    RegFunc = config.get(regHeader,"reg_func",default="",cast=str)
    #  Param Names
    ConfigParamNames = config.get(regHeader,"param_names",default="",cast=str)
    ParamNames = [x.strip() for x in ConfigParamNames.split(",")]
    #  Initial Parameters
    IPars = config.get(regHeader,"param_init",default="",cast=str)
    InitParams = [float(val) for val in IPars.split(",")]
    #  Limits
    PLow = config.get(regHeader,"param_lim_lo",default="",cast=str)
    PLimLo = [float(val) for val in PLow.split(",")]
    PHigh = config.get(regHeader,"param_lim_hi",default="",cast=str)
    PLimHi = [float(val) for val in PHigh.split(",")]
    RegRangeStr = config.get(regHeader,"reg_range",cast=str)
    RegRange = [float(val) for val in RegRangeStr.split(",")]


    # Unfolder Options
    unfoldHeader = "unfolder"
    UnfolderName = config.get(unfoldHeader,"unfolder_name",default="Unfolder",cast=str)
    # Analysis Bins
    binHeader = "analysisbins"
    binnumberStr = config.get(binHeader,"bin",default=0,cast=str)
    binnumbers = [val for val in binnumberStr.split(",")]
    stackFlag = config.get_boolean(binHeader,"stack",default=False)
    binList = ["bin"+val.replace(" ","") for val in binnumberStr.split(",")]
    nStack = len(binList)

    binnumber = 0
    if stackFlag is False:
        try: binnumber = np.int(binnumberStr)
        except: raise ValueError("\n\n*** You have requested more than 1 analysis bin without stacking. \n\tPlease fix your mistake. Exiting... ***\n")
        unfbinname = "bin%i"%binnumber
    else:
        if nStack<=1:
            mess = "\n**** You have request to stack analysis bins, but have only requested %i bin, %s."%(nStack,binList[0])
            mess += " You need at least 2 bins to stack.\n"
            raise ValueError(mess+"\n\tPlease correct your mistake. Exiting... ***\n")
        unfbinname = "bin0"

    # Unfolder!!!
    if stackFlag:
        UnfolderName += "_"+"".join(binList)
    # Test Statistic - Stat Function & Options
    tsHeader = "teststatistic"
    tsname = config.get(tsHeader,"ts_name",default="rmd",cast=str)
    # Mixer Name and Error Propagation Type
    mixHeader = "mixer"
    CovError = config.get(mixHeader,"error_type",default="",cast=str)

    Caxis_tab = []
    Cedges_tab = []
    flux_tab = []
    flux_stat_err_tab = []
    flux_sys_err_tab = []
    reg_fit_tab = []

    for iS in range(nStack):
        # Get Cuts from MCInput Stats File
        zlo = 0.
        zhi = 0.2928428
        if rr.objExists(StatsFile,"cuts",binList[iS],verbose=False):
            cuts = rr.get_tnamed(StatsFile,"cuts",binList[iS],verbose=False)
            zlo, zhi = getZenithLimits(cuts)

        # Output Subdirectory Name
        subDirName ="%s_%s_%s/%s"%(UnfolderName,tsname,CovError,binList[iS])

        # Get Measured (Unfolded) Cause Distribution
        Caxis, Cedges, n_meas, n_meas_err = rr.get1d(InputFile, NC_meas, subDirName)
        Cxlab, Cylab, Ctitle = rr.get_labels(InputFile, NC_meas, subDirName)
        Caxis_tab.append(Caxis)
        Cedges_tab.append(Cedges)
        # Get Statistical and MC (Systematic) Errors
        Cax, Ced, stat_err, dstat_err = rr.get1d(InputFile, NC_meas+"_errStat", subDirName)
        Cax, Ced, sys_err, dsys_err = rr.get1d(InputFile, NC_meas+"_errSys", subDirName)
        #  Energy Contribution to Bin Width
        Cdiff = np.diff(Cedges)
        # Get Observation Time and Zenith Bins
        Taxis, Tedges, obstime, obstime_err = rr.get1d(RawFile, TimeName, unfbinname)
        #  Solid Angle Contribution to Bin Width
        del_sr = 2*np.pi*(np.cos(zlo)-np.cos(zhi))
        # Geometrical Area Contribution
        geo_factor = (np.cos(zlo)+np.cos(zhi))/2
        # Total Bin Width
        bin_widths = Cdiff*obstime*del_sr*1e3*1e3*np.pi*geo_factor

        # Convert Unfolded Cause Distribution to Flux
        flux = n_meas/bin_widths
        flux_err = n_meas_err/bin_widths

        flux_stat_err = stat_err/bin_widths
        flux_sys_err = sys_err/bin_widths

        flux_tab.append(flux)
        flux_stat_err_tab.append(flux_stat_err)
        flux_sys_err_tab.append(flux_sys_err)

        NCFLUX = TH1F("%s_flux"%(NC_meas),"Cause Flux",len(Caxis),Cedges)
        NCFLUX.GetXaxis().SetTitle(Cxlab)
        NCFLUX.GetYaxis().SetTitle("Flux")
        NCFLUX.SetStats(0)
        for i in xrange(0,len(Caxis)):
            NCFLUX.SetBinContent(i+1,flux[i])
            NCFLUX.SetBinError(i+1,flux_err[i])

        # Prepare Regularizer
        Rglzr = Utils.Regularizer("REG",FitFunc=[RegFunc],Range=RegRange,InitialParams=InitParams,ParamLo=PLimLo,ParamHi=PLimHi,\
                            ParamNames=ParamNames,xarray=Caxis,xedges=Cedges,verbose=True,plot=False)

        # Fit Flux
        reg_fit, FitParams = Rglzr.Regularize(flux, flux_err)
        FLUXFIT = Rglzr.FIT

        reg_fit_tab.append(reg_fit)

        if (args.write_fit):
            RFile = ROOT.TFile(InputFile, "UPDATE")
            RFile.cd()
            # Save to Subdirectory
            if ( RFile.GetDirectory(subDirName) ):
                pdir = RFile.GetDirectory(subDirName)
                odir = pdir.mkdir("FitResults","Fit Results Directory")
                RFile.cd(subDirName+"/"+"FitResults")
                FLUXFIT.Write("fit")
                NCFLUX.Write()
                print "\n\n========== Writing Fit to ROOT File ==============="
                print "\nFit saved to %s\n"%InputFile
                print "\tin Sub-Directory %s\n"%(subDirName+"/FitResults")
                print "============= Done Writing ================="
            else:
                print "\n\nDirectory %s doesn't exists in file %s!"%(subDirName,InputFile)
                print "Consider running with a different configuration."
                print "Fit will not be written out to %s.\n"%InputFile

            RFile.Close()

    if (args.plot_local):

        scale = args.scale
        fig = plt.figure(figsize=(13.5,9))
        #fig = plt.figure(figsize=(9,6))
        mpl.rc("font", family="serif")
        ax = fig.add_subplot(111)
        ax.set_xscale('log', nonposx='clip')
        ax.set_yscale('log', nonposy='clip')
        ax.grid(True)
        #ax.set_xlabel(r'%s'%Cxlab)
        ax.set_xlabel(r'Energy (GeV)',fontsize=20)
        if (scale != 0.):
            ax.set_ylabel(r'Flux $\times$ E$^{%.2f}$ (s m$^{2}$ sr)$^{-1}$ (GeV)$^{%.02f}$'%(scale,scale-1),fontsize=20)
        else:
            ax.set_ylabel(r'Flux (GeV s m$^{2}$ sr)$^{-1}$',fontsize=20)
        ax.set_xlim(3e4,1.2e6)
        #ax.set_xlim(1e2,1.2e6)
        #ax.set_xlim(RegRange[0],RegRange[1])
        #ax.set_ylim(1e3,1e4)#np.max(top_err+flux))
        #ax.set_ylim(np.min(bot_err),9e-7)#np.max(top_err+flux))
        ax.set_title("All Species Cosmic Ray Energy Spectrum",fontsize=20)
        for iS in range(nStack):
            if (stackFlag):
                 lab_sys_err = r'$\sigma^{%s}_{sys}$'%binnumbers[iS]
                 lab_fit = r'$N_{%s} \times E^{-\gamma}$'%binnumbers[iS]
                 lab_stat_err = r'Flux$_{%s}$ + $\sigma_{stat}$'%binnumbers[iS]
            else:
                 lab_sys_err = r'$\sigma_{sys}$'
                 lab_fit = r'$N \times E^{-\gamma}$'
                 lab_stat_err = r'Flux + $\sigma_{stat}$'

            Caxis = Caxis_tab[iS]
            Cedges = Cedges_tab[iS]
            flux = flux_tab[iS]
            flux_stat_err = flux_stat_err_tab[iS]
            flux_sys_err = 2*flux_sys_err_tab[iS]
            reg_fit = reg_fit_tab[iS]

            reg_fit *= Caxis**scale
            flux *= Caxis**scale
            xerr = np.diff(Cedges)/2
            yerr = flux_stat_err*Caxis**scale
            top_err = flux_sys_err*Caxis**scale
            bot_err = flux_sys_err*Caxis**scale

            #ax.set_ylim(np.min(flux)*.05,1e-6)
            ax.set_ylim(np.min(flux)*.05,np.max(flux)*20)

            eb = ax.errorbar(Caxis, flux, xerr=0, yerr=yerr, fmt='o', \
                        label=lab_stat_err, linewidth=5, capsize=0, \
                        elinewidth=2, markersize=5, color=colors[iS%3])
            ax.bar(Caxis*1e5,np.ones(Caxis.shape), \
                   bottom=flux-bot_err,linewidth=0, \
                   #bottom=flux-bot_err,label=lab_sys_err,linewidth=0, \
                   width=Cdiff,alpha=.3,color='r',zorder=0,align='center')
            ax.fill_between(Caxis, flux-bot_err, flux+top_err, label=lab_sys_err,\
            #ax.fill_between(Caxis, flux-bot_err, flux+top_err,\
                            alpha=.3, edgecolor=colors[iS%3], facecolor=colors[iS%3], interpolate=True)
                            #alpha=.3, edgecolor='b', facecolor='b', interpolate=True)
            fpl = ax.plot(Caxis, reg_fit, linewidth=2, color='k', ls='-', label=lab_fit)
            pp = [eb[0],fpl[0]]
        #ax.grid(True)
        handles, labels = ax.get_legend_handles_labels()
        ncols = 1
        if stackFlag:
            ncols = 3
        ax.legend(handles, labels,loc='upper right', numpoints=1,ncol=ncols)
        #ax.legend(handles[::-1], labels[::-1],loc='upper right', numpoints=1)
        #plt.legend(loc='upper right',numpoints=1)
        #plt.legend(pp,[lab_stat_err],loc='upper right',numpoints=1)
        #plt.legend(pp,[lab_stat_err,lab_fit],loc='upper right',numpoints=1)
        ax.text(0.2, 0.1,'Preliminary', color='r', ha='center', va='center', transform=ax.transAxes, fontsize=28)

        #fig.savefig('CRFlux_wFit.png')
        plt.show()


if __name__ == "__main__":
    main()
