#!/usr/bin/env python
"""
   Class to perform an iterative unfolding,
   stopping when either max_iter or test
   statistic tolerance is reached.

.. codeauthor: Zig Hampel
"""

__version__ = "$Id"

try:
    import numpy as np
    import sys
    import os

    import ROOT
    from ROOT import TCanvas, TFile, TH2F, TH1F, TF1, TGraph, TGraphErrors, TGraphAsymmErrors
    from ROOT import gROOT, gPad, gRandom, gSystem, TDirectory
    gROOT.Reset()
    
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.ticker import ScalarFormatter
    
    from Mix import Mixer
    from Utils import *
    from Plotter import *
    from RootReader import mkDir
except ImportError as e:
    print e
    raise ImportError

class IterativeUnfolder():
    'Common base class for Unfolder'
    def __init__(self,name,maxIter=100,smoothIter=False,n_c=[],MixFunc=[],RegFunc=[],TSFunc=[],Stack=False,verbose=False,**kwargs):
        self.name = name
        self.maxIter = maxIter
        self.smoothIter = smoothIter
        # Stacking Info
        self.StackFlag = Stack
        self.nStack = len(RegFunc)
        # Stacking Bin Indices Initializer
        self.StackInd = [[] for i in range(self.nStack)]
        # Mixing, Regularizing, Test Statistic Functions
        self.Mix = MixFunc
        self.Rglzr = RegFunc
        self.tsFunc = TSFunc
        # If user doesn't provide prior...
        self.n_c = n_c.copy()
        # Count Iterations
        self.counter = 0

        # Lists to Append Iteration Results
        #  Unfolded Results
        self.n_c_iters = [[] for i in range(self.nStack)]
        self.n_c_labels = []
        self.n_c_labels.append("Prior")
        self.n_c_final = []
        self.covM = []
        self.statErr = []
        self.sysErr = []
        #  Regularization Results
        self.ParamOut = [[] for i in range(self.nStack)]
        for j in range(self.nStack):
            [self.ParamOut[j].append([]) for i in xrange(0,RegFunc[j].nParams)]
        self.N_measured_out = []
        self.N_measured_out.append(np.sum(n_c))

        #  Test Statistic Results
        self.TS_out = [[] for i in range(self.nStack)]

        # Set Stacking Bin Indices
        self.SetStackInd()

        # ROOT Output Info
        self.RootName = ""
        self.SubDirName = ""

        # X & Y labels for causes and effects
        self.Cxlab = r'$Causes$'
        self.Cylab = r'$N_{C}$'
        self.Exlab = r'$Effect$'
        self.Eylab = r'$N_{E}$'

    def GetIParams(self,ind):
        return self.ParamOut[ind]
    def GetINMeasured(self):
        return self.N_measured_out
    def GetITS(self,ind):
        return self.TS_out[ind]
    def GetINC(self,ind):
        return self.n_c_iters[ind], self.n_c_labels

    def SetStackInd(self):
        iind = 0
        for iS in range(self.nStack):
            iaxis = self.Rglzr[iS].xarray.copy()
            self.StackInd[iS] = [range(iind,iind+len(iaxis))]
            iind += len(iaxis)

    def passTol(self):
        passFlag = True
        for iS in range(self.nStack):
            passFlag *= self.tsFunc[iS].PassTol()
        return passFlag

    def IUnfold(self):
        print("\n====Starting Iterative Unfolding! Exciting!!!===\n")
        print("Iter: %i"%self.counter)
        # First Mixing
        n_c_update = self.Mix.Smear(self.n_c)
        self.counter += 1
        n_update = np.sum(n_c_update)

        # Calculate Chi2 of first mix with prior
        for iS in range(self.nStack):
            iind = self.StackInd[iS]
            inc = self.n_c[iind]
            self.n_c_iters[iS].append(inc)
            incu = n_c_update[iind]
            TS_cur, TS_del, TS_prob = self.tsFunc[iS].GetStats(incu,inc) 
            self.TS_out[iS].append(TS_cur)
        prev_fit = self.n_c.copy()

        reg_fit = np.zeros(self.n_c.shape)
        # Perform Iterative Unfolding
        while (not self.passTol() and self.counter < self.maxIter):
            self.n_c_labels.append("Iter %s"%self.counter)

            # Updated unfolded distribution
            n_c = n_c_update.copy()
            # For test statistic purposes. 
            # Don't want to compare reg_fit to n_c_update
            n_c_prev = n_c_update.copy()
            # Save n_update for root output
            self.N_measured_out.append(n_update)
            
            # Regularize unfolded dist in favor of 'smooth physics' prior
            for iS in range(self.nStack):
                iind = self.StackInd[iS][0]
                incu = n_c_update[iind]
                self.n_c_iters[iS].append(incu)
                reg_fit_iS, FitParams = self.Rglzr[iS].Regularize(incu)
                reg_fit[iind] = reg_fit_iS.copy()

                # Save fit parameters for root output
                for j in xrange(0,self.Rglzr[iS].nParams):
                    self.ParamOut[iS][j].append(FitParams[j])

                # Smoothed n_c for next iteration
                if (self.smoothIter):
                    n_c[iind] = reg_fit_iS.copy()
            #n_reg_events = np.sum(reg_fit)

            # Mix w/n_c from previous iter
            n_c_update = self.Mix.Smear(n_c)
            self.counter += 1
            n_update = np.sum(n_c_update)
            
            print("Iter: %i"%self.counter)
            # Calculate Chi2 wrt previous regularization
            for iS in range(self.nStack):
                iind = self.StackInd[iS][0]
                inc = self.n_c[iind]
                incu = n_c_update[iind]
                incp = n_c_prev[iind]
                TS_cur, TS_del, TS_prob = self.tsFunc[iS].GetStats(incu,incp)
                #TS_cur, TS_del, TS_prob = self.tsFunc[iS].GetStats(n_c_update[iind],n_c_prev[iind])
                self.TS_out[iS].append(TS_cur)
            prev_fit = reg_fit.copy()
            

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
        for iS in range(self.nStack):
            iind = self.StackInd[iS]
            incu = n_c_final[iind]
            self.n_c_iters[iS].append(incu)
            final_fit, FitParams = self.Rglzr[iS].Regularize(incu)

            for j in xrange(0,self.Rglzr[iS].nParams):
                self.ParamOut[iS][j].append(FitParams[j])

        ## Plot results of iterations
        self.n_c_labels.append("Final")

    # Set the Labels for ROOT Plots
    def SetAxesLabels(self,Exlab,Eylab,Cxlab,Cylab):
        self.Cxlab = Cxlab
        self.Cylab = Cylab
        self.Exlab = Exlab
        self.Eylab = Eylab

    # Save ROOT Output
    def WriteResults(self,OutFileName="",NCName="NC",NEName="NE",Eaxis=[],BinList=[],subdirName="",**kwargs):
        Cxlab= self.Cxlab
        Cylab= self.Cylab
        Exlab= self.Exlab
        Eylab= self.Eylab

        n_c_final_tot = self.n_c_final.copy()
        Vc = self.covM.copy()
        n_c_err_tot = np.sqrt(Vc.diagonal())

        for iS in range(self.nStack):
            sind = self.StackInd[iS]
            n_c_final = n_c_final_tot[sind]
            n_c_err = n_c_err_tot[sind]
            # Get the results
            ParamOut = self.GetIParams(iS)
            TS_out = self.GetITS(iS)
            N_measured_out = self.GetINMeasured()
            n_c_iters, n_c_labels = self.GetINC(iS)

            nIter = self.counter

            tsFunc = self.tsFunc[iS]

            Rglzr = self.Rglzr[iS]
            Caxis = Rglzr.xarray.copy()
            Cedges = Rglzr.xedges.copy()

            # Final Unfolded Cause Dist
            NCFINAL = TH1F(NCName,NCName,len(Caxis),Cedges)
            for i in xrange(0,len(Caxis)):
                NCFINAL.SetBinContent(i+1,n_c_final[i])
                NCFINAL.SetBinError(i+1,n_c_err[i])
            NCFINAL.SetTitle("Final Unfolded Cause Distribution")
            NCFINAL.GetXaxis().SetTitle(Cxlab)
            NCFINAL.GetYaxis().SetTitle(Cylab)
            NCFINAL.SetStats(0)
            NCFINAL.GetXaxis().SetLimits(Cedges[0],Cedges[-1])
            NCFINAL.GetYaxis().SetRangeUser(np.min(n_c_final[(n_c_final>0)])/2,np.max(n_c_final)*2)
            
            # Save Statistical and MC (Systematic) Errors
            errStat = self.statErr.copy()
            errSys = self.sysErr.copy()
            ERRSTAT = TH1F(NCName+"_errStat",NCName+"_errStat",len(Caxis),Cedges)
            ERRSYS = TH1F(NCName+"_errSys",NCName+"_errSys",len(Caxis),Cedges)
            for i in xrange(0,len(Caxis)):
                ERRSTAT.SetBinContent(i+1,errStat[i])
                ERRSYS.SetBinContent(i+1,errSys[i])

            # Linear space covering iterations
            ArrIter = np.linspace(0,nIter-1,nIter)

            # Observed Effects Dist
            NOBSERVED = TH1F(NEName,"Observed Effects Distribution",len(Eaxis)-1,Eaxis)
            NOBSERVED.GetXaxis().SetTitle(Exlab)
            NOBSERVED.GetYaxis().SetTitle(Eylab)
            NOBSERVED.SetStats(0)
            for j in xrange(0,len(Eaxis)-1):
                binval = self.Mix.NEobs[j]
                NOBSERVED.SetBinContent(j+1,binval)
                NOBSERVED.SetBinError(j+1,np.sqrt(binval))

            # Number of Events Measured vs Iteration
            NMEASURED = TGraph(nIter+1,np.append(ArrIter,nIter),np.asarray(N_measured_out))
            NMEASURED.SetTitle("Number of Events Measured vs Iteration")
            NMEASURED.GetXaxis().SetTitle('Iteration')
            NMEASURED.GetYaxis().SetTitle('N_{measured}')
   
            # Regularization Parameters's Values vs Iteration
            PARAMS = []
            for j in xrange(0,Rglzr.nParams):
                PARAMS.append(TGraph(nIter,ArrIter,np.asarray(ParamOut[j])))
                PARAMS[j].SetTitle(Rglzr.ParamNames[j]+' of Regularization Fit vs Iteration')
                PARAMS[j].GetXaxis().SetTitle('Iteration')
                PARAMS[j].GetYaxis().SetTitle(Rglzr.ParamNames[j])
            
            # Test Statistic vs Iteration
            TS_TGRAPH = TGraph(nIter,ArrIter,np.asarray(TS_out))
            TS_TGRAPH.SetTitle('Test Statistic vs Iteration')
            TS_TGRAPH.GetXaxis().SetTitle('Iteration')
            TS_TGRAPH.GetYaxis().SetTitle(tsFunc.name)

            # Final Covariance and Coefficient of Correlation Matrices
            COV = TH2F("COVM","Covariance Matrix - %s"%self.Mix.ErrorType,len(Caxis),Cedges,len(Caxis),Cedges)
            COV.GetXaxis().SetTitle(Cxlab)
            COV.GetYaxis().SetTitle(Cxlab)
            COV.SetStats(0)
            CORR = TH2F("CORR","Coefficient of Correlation Matrix - %s"%self.Mix.ErrorType,len(Caxis),Cedges,len(Caxis),Cedges)
            CORR.GetXaxis().SetTitle(Cxlab)
            CORR.GetYaxis().SetTitle(Cxlab)
            CORR.SetStats(0)
            for j in xrange(0,len(Caxis)):
                diag = Vc[j,j]
                for k in xrange(0,len(Caxis)):
                    val = Vc[j,k]
                    COV.SetBinContent(j+1,k+1,val)
                    if (diag > 0):
                        CORR.SetBinContent(j+1,k+1,val/diag)

            # Save Each Iteration of True Distribution
            NCUNFOLDED = []
            for j in xrange(0,nIter+1):
                NCUNFOLDED.append(TH1F("%s_%i"%(NCName,j),"Cause Distribution - %s"%n_c_labels[j],len(Caxis),Cedges))
                NCUNFOLDED[j].GetXaxis().SetTitle(Cxlab)
                NCUNFOLDED[j].GetYaxis().SetTitle(Cylab)
                NCUNFOLDED[j].SetStats(0)
                n_c_j = np.asarray(n_c_iters[j]).copy()
                for i in xrange(0,len(Caxis)):
                    NCUNFOLDED[j].SetBinContent(i+1,n_c_j[i])

            # Write to ROOT File
            RFile = ROOT.TFile(OutFileName, "UPDATE")
            RFile.cd()

            self.RootName = OutFileName
            self.SubDirName = subdirName

            iDir = subdirName + "/" + BinList[iS]

            tdir = RFile.GetDirectory(subdirName)
            if ( not tdir ):
                tdir = RFile.mkdir(subdirName,"%s Directory"%self.name)
            RFile.cd(subdirName)

            if ( not RFile.GetDirectory(iDir) ):
                pdir = RFile.mkdir(iDir,"%s Directory"%BinList[iS])
                RFile.cd(iDir)

                NCFINAL.Write(NCName)
                ERRSTAT.Write()
                ERRSYS.Write()
                COV.Write("CovM")
                CORR.Write("CorrM")
                NOBSERVED.Write(NEName)
                NMEASURED.Write("%sMeasured"%NCName.lower())
                TS_TGRAPH.Write("%s"%(tsFunc.name))
                for j in xrange(0,Rglzr.nParams):
                    PARAMS[j].Write("%s"%(Rglzr.ParamNames[j]))

                # Write Unfolded Dists to Nested Directory
                ncDirName = "%s/%s"%(iDir,NCName.lower())
                ncdir = RFile.mkdir(ncDirName,"Unfolded Distributions for All Iterations")
                RFile.cd(ncDirName)
                for j in xrange(0,nIter+1):
                    NCUNFOLDED[j].Write("%s_%i"%(NCName,j))

                print "\n\n========== Writing ROOT File ==============="
                print "\nResults ROOT file saved as %s\n"%OutFileName
                print "\tSub-Directory is %s\n"%(iDir)
                print "============= Done Writing ================="
            else:
                print "\n\nDirectory %s already exists in file %s!"%(iDir,OutFileName)
                print "Consider running with a different configuration."
                print "Results will not be written out to %s.\n"%OutFileName

            RFile.Close()
