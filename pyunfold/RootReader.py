"""
   Root file reading routines for getting
   various plotting objects from specified
   root file.
"""

import numpy as np

import ROOT
from ROOT import gDirectory, TChain, TCanvas, TFile, TProfile, TNtuple
from ROOT import TH1, TH1D, TH2F, TF1, TH1F, TGraph, TGraphErrors, TLine
from ROOT import gROOT, gStyle, gBenchmark, gRandom, gSystem, gDirectory
from ROOT import gPad, TText, TLatex, TMarker, TColor, TNamed

ROOT.gROOT.SetBatch(True)
TH1.AddDirectory(False)


def objExists(infile, objname, subdir='', verbose=False):
    """
    Test if object exists in input file
    objExists(input file, object name, sub-directory='', verbose=False)
    """
    # Try to Open File
    fin = ROOT.TFile(infile, "READ")
    file_open = fin.IsOpen()
    err_mess = "\n***Can not open ROOT file %s***\n" %(infile)
    assert (file_open),err_mess

    if (verbose):
        print "\n============= Opening %s\n" %infile

    direxists = False
    if ( subdir == ''):
        direxists = True
    if ( fin.GetDirectory(subdir) ):
        fin.cd(subdir)
        if (verbose):
            print "\n============= Accessing directory %s\n" %subdir
        direxists = True
    err_mess = "\n***Directory %s not in ROOT file %s. Exiting...***\n" %(subdir,infile)
    assert (direxists),err_mess

    # Check if Object is in File
    objexists = gDirectory.GetListOfKeys().Contains(objname)
    if (verbose):
        if objexists:
            print "\n*Object %s exists and is in %s***\n"%(objname,infile)
        else:
            print "\n***Object %s isn't in ROOT file %s***\n" %(objname,infile)

    return objexists


def getter(infile, objname, subdir='', verbose=False):
    """
    ROOT object getter
    getter(input file, object name, sub-directory='', verbose=False)
    """

    # Try to Open File
    fin = ROOT.TFile(infile, "READ")
    file_open = fin.IsOpen()
    err_mess = "\n***Can not open ROOT file %s***\n" %(infile)
    assert (file_open),err_mess

    if (verbose):
        print "\n============= Opening %s\n" %infile

    direxists = False
    if ( subdir == ''):
        direxists = True
    if ( fin.GetDirectory(subdir) ):
        fin.cd(subdir)
        if (verbose):
            print "\n============= Accessing directory %s\n" %subdir
        direxists = True
    err_mess = "\n***Directory %s not in ROOT file %s. Exiting...***\n" %(subdir,infile)
    assert (direxists),err_mess


    # Check if Object is in File
    objexists = gDirectory.GetListOfKeys().Contains(objname)
    err_mess = "\n***Object %s isn't in ROOT file %s***\n" %(objname,infile)
    assert (objexists),err_mess

    if (verbose):
        print "\n============= Got %s\n" %objname

    # Get the Object
    obj = gDirectory.Get(objname)
    objname = str(obj.IsA())

    if ("TH" in objname):
        obj.SetDirectory(0)

    fin.Close()

    return obj

def mkDir(infile, subdir, verbose=False):
    """
    ROOT directory make-r or cd-er
    mkDir(infile, subdir, verbose=False)
    """

    # Write to ROOT File
    RFile = ROOT.TFile(infile, "UPDATE")
    RFile.cd()

    sdir = []

    if ( not RFile.GetDirectory(subdir) ):
        sdir = RFile.GetDirectory(subdir)
    else:
        sdir = getDir(infile, subdir, verbose=False)

    return sdir

def getDir(infile, subdir, verbose=False):
    """
    ROOT directory cd-er
    getDir(infile, subdir, verbose=False)
    """

    # Write to ROOT File
    RFile = ROOT.TFile(infile, "UPDATE")
    RFile.cd()

    sdir = []

    if ( RFile.GetDirectory(subdir) ):
        sdir = RFile.GetDirectory(subdir)

    err_mess = "\n*** Directory %s doesn't exist in %s. Exiting...***\n"%(subdir,infile)
    assert (sdir),err_mess

    return sdir


def getDirs(infile, verbose=False):
    """
    ROOT directory lister
    getDirs(infile, verbose=False)
    """

    # Try to Open File
    fin = ROOT.TFile(infile, "READ")
    file_open = fin.IsOpen()
    err_mess = "\n***Can not open ROOT file %s***\n" %(infile)
    assert (file_open),err_mess

    # List object for directories
    dir_list = []

    # Loop through keys and get directories
    for key in fin.GetListOfKeys():
        kname = key.GetName()
        if ( fin.GetDirectory(kname) and kname != ''):
            dir_list.append(kname)

    fin.Close()

    err_mess = "\n*** There are no directories in this %s ***\n" %(infile)
    assert (dir_list),err_mess

    return dir_list


def getObjs(infile, verbose=False):
    """
    ROOT object lister
    getObjs(infile, verbose=False)
    """

    # Try to Open File
    fin = ROOT.TFile(infile, "READ")
    file_open = fin.IsOpen()
    err_mess = "\n***Can not open ROOT file %s***\n" %(infile)
    assert (file_open),err_mess

    # List object for directories
    obj_list = []

    # Loop through keys and get directories
    for key in fin.GetListOfKeys():
        kname = key.GetName()
        if key.IsFolder():
            continue
        if ( kname != ''):
            obj_list.append(kname)

    fin.Close()

    err_mess = "\n*** There are no objects in this %s ***\n" %(infile)
    assert (obj_list),err_mess

    return obj_list

def get1d(infile, histname, subdir='',verbose=False):
    """
    1d root Hist getter
    get1d(infile, histname, subdir='',verbose=False)
    """

    ### 1d Histogram
    Hist = getter(infile,histname,subdir,verbose)

    nbinsHist = Hist.GetSize()-2
    axesX = np.zeros(nbinsHist)
    edgesX = np.zeros(nbinsHist+1)
    Arr = np.zeros(nbinsHist)
    dArr = np.zeros(nbinsHist)
    for j in xrange(0,nbinsHist):
        axesX[j] = Hist.GetBinCenter(j+1)
        edgesX[j] = Hist.GetBinLowEdge(j+1)
        Arr[j] = Hist.GetBinContent(j+1)
        dArr[j] = Hist.GetBinError(j+1)
    edgesX[nbinsHist] = Hist.GetBinLowEdge(nbinsHist+1)

    return axesX, edgesX, Arr, dArr

def get2d(infile, histname, subdir='',verbose=False):
    """
    2d Hist getter
    get2d(infile, histname, subdir='',verbose=False)
    """

    ## 2d Histogram
    Hist = getter(infile,histname,subdir,verbose)

    nbinsX, nbinsY = Hist.GetNbinsX(), Hist.GetNbinsY()
    Arr = np.zeros((nbinsY,nbinsX))
    dArr = np.zeros((nbinsY,nbinsX))
    axesX = np.zeros(nbinsX)
    axesY = np.zeros(nbinsY)
    edgesX = np.zeros(nbinsX+1)
    edgesY = np.zeros(nbinsY+1)
    for j in xrange(0,nbinsX):
        axesX[j] = Hist.GetXaxis().GetBinCenter(j+1)
        edgesX[j] = Hist.GetXaxis().GetBinLowEdge(j+1)
    edgesX[nbinsX] = Hist.GetXaxis().GetBinLowEdge(nbinsX+1)

    for j in xrange(0,nbinsY):
        axesY[j] = Hist.GetYaxis().GetBinCenter(j+1)
        edgesY[j] = Hist.GetYaxis().GetBinLowEdge(j+1)
    edgesY[nbinsY] = Hist.GetYaxis().GetBinLowEdge(nbinsY+1)

    axes = [axesX, axesY]
    edges = [edgesX, edgesY]

    for j in xrange(0,nbinsX):
        for k in xrange(0,nbinsY):
            Arr[k,j] = Hist.GetBinContent(j+1,k+1)
            dArr[k,j] = Hist.GetBinError(j+1,k+1)

    return axes, edges, Arr, dArr

def get3d(infile, histname, subdir='',verbose=False):
    """
    3d Hist getter
    get3d(infile, histname, subdir='',verbose=False)
    """

    ## 2d Histogram
    Hist = getter(infile,histname,subdir,verbose)

    nbinsX, nbinsY, nbinsZ = Hist.GetNbinsX(), Hist.GetNbinsY(), Hist.GetNbinsZ()
    Arr = np.zeros((nbinsZ,nbinsY,nbinsX))
    dArr = np.zeros((nbinsZ,nbinsY,nbinsX))
    axesX = np.zeros(nbinsX)
    axesY = np.zeros(nbinsY)
    axesZ = np.zeros(nbinsZ)
    edgesX = np.zeros(nbinsX+1)
    edgesY = np.zeros(nbinsY+1)
    edgesZ = np.zeros(nbinsZ+1)
    for j in xrange(0,nbinsX):
        axesX[j] = Hist.GetXaxis().GetBinCenter(j+1)
        edgesX[j] = Hist.GetXaxis().GetBinLowEdge(j+1)
    edgesX[nbinsX] = Hist.GetXaxis().GetBinLowEdge(nbinsX+1)

    for j in xrange(0,nbinsY):
        axesY[j] = Hist.GetYaxis().GetBinCenter(j+1)
        edgesY[j] = Hist.GetYaxis().GetBinLowEdge(j+1)
    edgesY[nbinsY] = Hist.GetYaxis().GetBinLowEdge(nbinsY+1)

    for j in xrange(0,nbinsZ):
        axesZ[j] = Hist.GetZaxis().GetBinCenter(j+1)
        edgesZ[j] = Hist.GetZaxis().GetBinLowEdge(j+1)
    edgesZ[nbinsZ] = Hist.GetZaxis().GetBinLowEdge(nbinsZ+1)

    axes = [axesX, axesY, axesZ]
    edges = [edgesX, edgesY, edgesZ]

    for j in xrange(0,nbinsX):
        for k in xrange(0,nbinsY):
            for l in xrange(0,nbinsZ):
                Arr[l,k,j] = Hist.GetBinContent(j+1,k+1,l+1)
                dArr[l,k,j] = Hist.GetBinError(j+1,k+1,l+1)

    return axes, edges, Arr, dArr


def getTG(infile, tgname, subdir='', verbose=False):
    """
    1d root TGraphErrors getter
    getTG(infile, tgname, subdir='', verbose=False)
    """

    ## TGraphErrors
    TG = getter(infile,tgname,subdir,verbose)

    nbins = TG.GetN()
    axesX = np.zeros(nbins)
    Arr = np.zeros(nbins)
    edgesX = np.zeros(nbins)
    dArr = np.zeros(nbins)
    for j in xrange(0,nbins):
        x = np.zeros(1)
        y = np.zeros(1)
        TG.GetPoint(j,x,y)
        axesX[j] = x[0]
        Arr[j] = y[0]
        edgesX[j] = TG.GetErrorX(j)
        dArr[j] = TG.GetErrorY(j)

    return axesX, edgesX, Arr, dArr

def get_labels(infile, histname, subdir='',verbose=False):
    """
    ROOT Hist label getter
    get_labels(infile, histname, subdir='',verbose=False)
    """

    ## 1d Histogram
    Hist = getter(infile,histname,subdir,verbose)

    xlab = Hist.GetXaxis().GetTitle()
    ylab = Hist.GetYaxis().GetTitle()
    title = Hist.GetTitle()

    return xlab, ylab, title


def get_tnamed(infile, tnamed, subdir='',verbose=False):
    """
    ROOT Hist label getter
    get_labels(infile, histname, subdir='',verbose=False)
    """

    ## TNamed Object
    tnamed = getter(infile,tnamed,subdir,verbose)
    return tnamed.GetTitle()
