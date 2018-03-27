#!/usr/bin/env python
"""
   Run script to perform unfolding.
   Requires config script via command line
   option -c if standalone, or first argument
   if calling Unfold in a script.
"""

import numpy as np
import argparse
import sys
import os

import Mix
import Utils
import LoadStats
import IterUnfold

import Plotter as pltr
from .Plotter import plt, mpl
import RootReader as rr


# Global program options
global args
args = None


# Allow for function functionality
def Unfold(config_name="",return_dists=False,EffDist=None,plot_local=False, **kwargs):

    if (not config_name): raise Exception, "Need a config file!"

    #### Get the Configuration Parameters from the Config File ####
    config = Utils.ConfigFM(config_name)

    # Input Data ROOT File Name
    dataHeader = "data"
    InputFile = config.get(dataHeader,"inputfile",default="",cast=str)
    NE_meas_name = config.get(dataHeader,"ne_meas",default="",cast=str)
    isMC = config.get_boolean(dataHeader,"ismc",default=False)

    # Output ROOT File Name
    outHeader = "output"
    WriteOutput = config.get_boolean(outHeader,"write_output",default=False)
    OutFile = config.get(outHeader,"outfile",default="",cast=str)
    PlotDir = config.get(outHeader,"plotdir",default=os.path.dirname(OutFile),cast=str)
    NC_meas = config.get(outHeader,"nc_meas",default="",cast=str)

    if (OutFile == ""):
        InDir = os.path.dirname(InputFile)
        InBase = os.path.basename(InputFile)
        OutFile = "unfolded_"+InBase

    # Analysis Bin
    binHeader = "analysisbins"
    binnumberStr = config.get(binHeader,"bin",default=0,cast=str)
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

    # Unfolder Options
    unfoldHeader = "unfolder"
    UnfolderName = config.get(unfoldHeader,"unfolder_name",default="Unfolder",cast=str)
    UnfMaxIter = config.get(unfoldHeader,"max_iter",default=100,cast=int)
    UnfSmoothIter = config.get_boolean(unfoldHeader,"smooth_with_reg",default=False)
    UnfVerbFlag = config.get_boolean(unfoldHeader,"verbose",default=False)

    # Get Prior Definition
    priorHeader = "prior"
    priorString = config.get(priorHeader,"func",default="Jeffreys",cast=str)
    priorList = [val for val in priorString.split(",")]
    nPriorFuncs = len(priorList)
    if nPriorFuncs != nStack:
        if nPriorFuncs == 1:
            for i in range(nStack-1):
                priorList.append(priorList[0])
        else:
            mess = "\n**** You have requested an incorrect number of prior functions (%i), but there are %i analysis bins."%(nPriorFuncs,nStack)
            raise ValueError(mess+"\n\tPlease correct your mistake. Exiting... ***\n")

    # Mixer Name and Error Propagation Type
    '''Options: ACM, DCM'''
    mixHeader = "mixer"
    MixName = config.get(mixHeader,"mix_name",default="",cast=str)
    CovError = config.get(mixHeader,"error_type",default="",cast=str)

    # Test Statistic - Stat Function & Options
    '''Options: chi2, rmd, pf, ks'''
    tsHeader = "teststatistic"
    tsname = config.get(tsHeader,"ts_name",default="rmd",cast=str)
    tsTol = config.get(tsHeader,"ts_tolerance",cast=float)
    tsRangeStr = config.get(tsHeader,"ts_range",cast=str)
    tsRange = [float(val) for val in tsRangeStr.split(",")]
    tsVerbFlag = config.get_boolean(tsHeader,"verbose",default=False)

    # Regularization Function, Initial Parameters, & Options
    regHeader = "regularization"
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
    #  Options
    RegPFlag = config.get_boolean(regHeader,"plot",default=False)
    RegVFlag = config.get_boolean(regHeader,"verbose",default=False)

    # Get MCInput
    mcHeader = "mcinput"
    StatsFile = config.get(mcHeader,"stats_file",default="",cast=str)
    Eff_hist_name = config.get(mcHeader,"eff_hist",default="",cast=str)
    MM_hist_name = config.get(mcHeader,"mm_hist",default="",cast=str)

    #### Setup the Observed and MC Data Arrays ####
    # Load MC Stats (NCmc), Cause Efficiency (Eff) and Migration Matrix ( P(E|C) )
    MCStats = LoadStats.MCTables(StatsFile,BinName=binList,RespMatrixName=MM_hist_name,EffName=Eff_hist_name,Stack=stackFlag)
    Caxis = []
    Cedges = []
    modelList = []
    cutList = []
    for index in range(nStack):
        axis, edge = MCStats.GetCauseAxis(index)
        Caxis.append(axis)
        Cedges.append(edge)
        if rr.objExists(StatsFile,"ModelName",binList[index],verbose=False):
            modelList.append(rr.get_tnamed(StatsFile,"ModelName",binList[index],verbose=False))
    Eaxis, Eedges = MCStats.GetEffectAxis()
    # Effect and Cause X and Y Labels from Respective Histograms
    Cxlab, Cylab, Ctitle = rr.get_labels(StatsFile,Eff_hist_name,binList[0],verbose=False)

    # Load the Observed Data (n_eff), define total observed events (n_obs)
    #  Get from ROOT input file if requested
    if (EffDist is None):
        Exlab, Eylab, Etitle = rr.get_labels(InputFile,NE_meas_name,unfbinname,verbose=False)
        Eaxis, Eedges, n_eff, n_eff_err = rr.get1d(InputFile,NE_meas_name,unfbinname)
        EffDist = Utils.DataDist(Etitle,data=n_eff,error=n_eff_err,axis=Eaxis,edges=Eedges,xlabel=Exlab,ylabel=Eylab,units="")
    #  Otherwise get from input EffDist object
    Exlab = EffDist.xlab
    Eylab = EffDist.ylab
    Etitle = EffDist.name
    n_eff = EffDist.getData()
    n_eff_err = EffDist.getError()
    n_obs = np.sum(n_eff)

    # Initial best guess (0th prior) expected prob dist (Jeffrey's Prior)
    n_c = Utils.UserPrior(priorList, Caxis, n_obs)

    print "\n\n====================================================\n"
    print "Total number of observed events:\t%e"%n_obs
    print "\n====================================================\n\n"
    if (plot_local):
        xlim = [Eedges[0],Eedges[-1]]
        ylim = [np.min(n_eff[(n_eff>0)]-np.sqrt(n_eff[(n_eff>0)])),np.max(n_eff)*1.4]
        figure = pltr.ebar(Eaxis, [n_eff], np.diff(Eedges)/2, [n_eff_err], r'Reconstructed Energy E$_{reco}$ (GeV)', Eylab, \
                  "Observed Counts", [""], xlim, ylim, 'x,y')
        figure.savefig("%s/ObservedCounts_%s.png"%(PlotDir,unfbinname))

    #### Setup the Tools Used in Unfolding ####
    # Prepare Regularizer
    Rglzr = [Utils.Regularizer("REG",FitFunc=[RegFunc],Range=RegRange,InitialParams=InitParams,ParamLo=PLimLo,ParamHi=PLimHi,\
                               ParamNames=ParamNames,xarray=Caxis[i],xedges=Cedges[i],verbose=RegVFlag,plot=RegPFlag) for i in range(nStack)]
    # Prepare Test Statistic-er
    tsMeth = Utils.get_ts(tsname)
    tsFunc = [tsMeth(tsname,tol=tsTol,Xaxis=Caxis[i],TestRange=tsRange,verbose=tsVerbFlag) for i in range(nStack)]

    # Prepare Mixer
    Mixer = Mix.Mixer(MixName,ErrorType=CovError,MCTables=MCStats,EffectsDist=EffDist)

    # Unfolder!!!
    if stackFlag:
        UnfolderName += "_"+"".join(binList)
    Unfolder = IterUnfold.IterativeUnfolder(UnfolderName,maxIter=UnfMaxIter,smoothIter=UnfSmoothIter,n_c=n_c,\
                                            MixFunc=Mixer,RegFunc=Rglzr,TSFunc=tsFunc,Stack=stackFlag,verbose=UnfVerbFlag)

    # Iterate the Unfolder
    Unfolder.IUnfold()

    if (WriteOutput):
        # Output Subdirectory Name
        subDirName = "%s_%s_%s"%(UnfolderName,tsname,CovError)

        # Write the Results to ROOT File
        Unfolder.SetAxesLabels(Exlab,Eylab,Cxlab,Cylab)
        Unfolder.WriteResults(OutFileName=OutFile,NCName=NC_meas,NEName=NE_meas_name,Eaxis=Eedges,BinList=binList,subdirName=subDirName)

        # If MC, Write Thrown Cause Distribution
        if (isMC):
            NC_Thrown = config.get("data","nc_thrown",subDirName,cast=str)
            Thrown_hist = rr.getter(InputFile,NC_Thrown)
            RFile = ROOT.TFile(Unfolder.RootName, "UPDATE")
            pdir = RFile.GetDirectory(Unfolder.SubDirName)
            pdir.cd()
            Thrown_hist.Write("NCThrown")
            RFile.Close()


    #### Plot Some Results of the Unfolding for the User
    # Get the results
    print "\n\n================ Unfolding Stats ================="
    for iS in range(Unfolder.nStack):
        n_c_iters, n_c_labels = Unfolder.GetINC(iS)
        nIter = Unfolder.counter
        Caxis = Rglzr[iS].xarray
        Cedges = Rglzr[iS].xedges

        print "\n"
        print "\tModel Name: \t\t%s"%modelList[iS]
        print "\tBin in Resp File: \t%s"%binList[iS]
        print "\tEvents Observed: \t%e"%n_obs
        print "\tEvents Measured: \t%e"%np.sum(n_c_iters[-1])
        print "\tUnfolding Efficiency:\t%e"%(n_obs/np.sum(n_c_iters[-1]))

        if (plot_local):
            # Plot Some Results
            title = "Unfolded Integral Cause Distribution - Iterations"
            if (stackFlag):
                title = "Unfolded Iterations - Model: %s"%modelList[iS]
            xlim = [Caxis[0],Caxis[-1]]
            ylim = [max(np.min(n_c_iters[:]),1),np.max(n_c_iters[:])*1.1]
            #  Plot 1st 5 iters, last iter, and mult of 10
            pln_c = []
            pln_l = []
            for i in xrange(0,nIter+1):
                if (i<6 or i==nIter or i%10==0):
                    pln_c.append(n_c_iters[i])
                    pln_l.append(n_c_labels[i])
            figure = pltr.oplot(pln_c, Caxis, r'%s'%Cxlab, 'Frequency', r'%s'%title, pln_l, xlim=xlim, ylim=ylim, log='x,y')
            figure.savefig("%s/Iterations_%s.png"%(PlotDir,modelList[iS].replace(" ", "")))

    if (plot_local and stackFlag):
        colors = ["red", "blue", "green"]
        styles = ["-", "--", "-."]
        fig = plt.figure(figsize=(9,6))
        mpl.rc("font", family="serif", size=14)
        ax = fig.add_subplot(111)
        ax.set_xlabel(Cxlab)
        ax.set_ylabel("Frequency")
        ax.set_title("Stack Comparison")
        #ax.set_xlim(1e2,2e6)
        #ax.set_ylim(1e4,1e14)
        ax.set_xscale('log', nonposx='clip')
        ax.set_yscale('log', nonposy='clip')
        for iS in range(Unfolder.nStack):
            n_c_iters, n_c_labels = Unfolder.GetINC(iS)
            Caxis = Rglzr[iS].xarray
            Cedges = Rglzr[iS].xedges
            ax.plot(Caxis,n_c_iters[-1], label="%s"%modelList[iS], color=colors[iS%3],ls=styles[iS%3])
        plt.legend(loc='best')
        plt.show()
        fig.savefig("%s/StackResults_%s.png"%(PlotDir,'_'.join(modelList).replace(" ", "")))

    print "\n\n================ Exiting Unfolder ================\n\n"

    # Return N_C, Error Estimates
    if (return_dists):
        n_c = Unfolder.n_c_final
        n_c_err = np.sqrt(np.diagonal(Unfolder.covM))
        CauseDist = Utils.DataDist(Ctitle,data=n_c,error=n_c_err,axis=Caxis,edges=Cedges,xlabel=Cxlab,ylabel=Cylab,units="GeV")
        CauseDist.setStatErr(Unfolder.statErr)
        CauseDist.setSysErr(Unfolder.sysErr)
        return CauseDist


# Allow for standalone execution
if __name__ == '__main__':
    p = argparse.ArgumentParser(description="Script to Run an Unfolding")
    p.add_argument("-c", "--config", dest="config_name", required=True, help="Configuration file")
    p.add_argument("-p", "--plot", dest="plot", action="store_true", help="Friendly plotter")

    try:
        args = p.parse_args()
    except:
        p.print_help()
        sys.exit(0)

    if (not args.config_name): raise Exception, "Need a config file!"

    Unfold(config_name=args.config_name,return_dists=False,EffDist=None,plot_local=args.plot)
