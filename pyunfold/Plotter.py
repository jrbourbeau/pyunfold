#!/usr/bin/env python
"""Compendium of matplolib plotting routines
"""

import numpy as np
import bisect as bs
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import LogNorm
from matplotlib import rcParams

from .Utils import none_to_empty_list


# rcParam adjustments
mpl.rc("font", family="serif", size=18)
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.size': 20})
rcParams.update({'axes.titlesize':18})
rcParams.update({'axes.labelsize':16})
rcParams.update({'xtick.labelsize':15})
rcParams.update({'ytick.labelsize':15})
rcParams.update({'legend.fontsize':14})
rcParams.update({'legend.numpoints':1})
rcParams.update({'grid.linewidth':0.2})

# Colorbar options
cbarfontsize = 14
cbarticksize = 0.5

# Marker/Line Options
colors = ["blue", "red", "green", "black"]
colorsmall = ["b", "r", "g", "k"]
styles = ["-", "--", "-."]

# Get number of columns for legend
def nLegCols(nLines):
    groupsize = np.int(4)
    nLines = np.int(nLines)
    ncols = nLines / groupsize + 1
    return ncols

# Set INF for easier limits check
INF = np.Inf

# Make a pcolormesh plot with options for log scaling
# on all axes based on final input parameter
def cmesh(data,x,y,title,xlabel,ylabel,zlabel,lim,log):
    '''Single colormesh plot.
       cmesh(data,x,y,title,xlabel,ylabel,zlabel,lim,log)'''
    xx, yy = np.meshgrid(x,y)
    fig = plt.figure(figsize=(10,5.75))
    ax = fig.add_subplot(111)
    if ('z' in log):
        pl = ax.pcolormesh(xx,yy,data,norm=LogNorm(vmin=lim[0], vmax=lim[1]))
    else:
        pl = ax.pcolormesh(xx,yy,data,vmin=lim[0], vmax=lim[1])
    # Colorbar and label
    cb = fig.colorbar(pl)
    cb.set_label(zlabel)
    cb.ax.tick_params(labelsize=cbarfontsize, width=cbarticksize)
    if ('x' in log):
        ax.set_xscale('log', nonposx='clip')
    if ('y' in log):
        ax.set_yscale('log', nonposy='clip')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(x[0],x[-1])
    ax.set_ylim(y[0],y[-1])
    plt.show()
    return fig


def cplot(origdata,bins,x,y,title,xlabel,ylabel,zlabel,lim,log):
    '''Side by side comparison colormesh plots.
       cplot(origdata,bins,x,y,name,lim,log)'''
    xx, yy = np.meshgrid(x,y)
    fig = plt.figure(figsize=(17,6))
    ax = fig.add_subplot(121)
    if ('z' in log):
        pl = ax.pcolormesh(xx,yy,origdata,norm=LogNorm(vmin=lim[0], vmax=lim[1]))
    else:
        pl = ax.pcolormesh(xx,yy,origdata,vmin=lim[0], vmax=lim[1])
    #cb = fig.colorbar(pl)
    #cb.set_label(zlabel)
    if ('x' in log):
        ax.set_xscale('log', nonposx='clip')
    if ('y' in log):
        ax.set_yscale('log', nonposy='clip')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(x[0],x[-1])
    ax.set_ylim(y[0],y[-1])
    fig.subplots_adjust(left=0.125, right=0.85)

    ax = fig.add_subplot(122)
    if ('z' in log):
        pl = ax.pcolormesh(xx,yy,bins,norm=LogNorm(vmin=lim[0], vmax=lim[1]))
    else:
        pl = ax.pcolormesh(xx,yy,bins,vmin=lim[0], vmax=lim[1])
    cb = fig.colorbar(pl)
    cb.set_label(zlabel)
    cb.ax.tick_params(labelsize=cbarfontsize, width=cbarticksize)
    if ('x' in log):
        ax.set_xscale('log', nonposx='clip')
    if ('y' in log):
        ax.set_yscale('log', nonposy='clip')
    ax.set_xlabel(xlabel)
    #ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(x[0],x[-1])
    ax.set_ylim(y[0],y[-1])
    fig.subplots_adjust(left=0.125, right=0.85)
    plt.show()
    return fig

# Make errorbar plot
def ebar(x, y, xerr, yerr, xlab, ylab, title, labels, xlim, ylim, log):
    """
    Small errorbar plotting function
    ebar(x, y, xerr, yerr, xlab, ylab, title, labels, xlim, ylim, log)
    """
    ntabs = len(y)
    fig = plt.figure(figsize=(9,6))
    mpl.rc("font", family="serif")
    ax = fig.add_subplot(111)
    if ('x' in log):
        ax.set_xscale('log', nonposx='clip')
    if ('y' in log):
        ax.set_yscale('log', nonposy='clip')
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.set_xlim(xlim[0],xlim[1])
    ax.set_ylim(ylim[0],ylim[1])
    #ax.text(0.8, 0.9,'Preliminary', color='r', ha='center', va='center', transform=ax.transAxes)
    # Loop over provided data
    for i in range(len(y)):
        pl = ax.errorbar(x, y[i], xerr=xerr, yerr=yerr[i], fmt='.', label=labels[i], linewidth=5, capsize=0, elinewidth=1, markersize=0)
    nlegcol = nLegCols(ntabs)
    plt.legend(loc='best',numpoints=1, ncol=nlegcol)
    plt.show()
    return fig

def oplot(tabs, xarray, xlab, ylab, title, labels, xlim, ylim, log):
    """
    Small plotting function
    oplot(tabs, xarray, xlab, ylab, title, labels, xlim, ylim, log)
    """
    ntabs = len(tabs)
    fig = plt.figure(figsize=(9,6))
    mpl.rc("font", family="serif")
    ax = fig.add_subplot(111)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.set_xlim(xlim[0],xlim[1])
    ax.set_ylim(ylim[0],ylim[1])
    if ('x' in log):
        ax.set_xscale('log', nonposx='clip')
    if ('y' in log):
        ax.set_yscale('log', nonposy='clip')
    legflag = True
    if (len(labels) == 0):
        labels = ['x']*len(tabs)
        legflag = False
    # Loop over provided data
    for i in range(ntabs-1):
        pl = ax.plot(xarray, tabs[i], label=labels[i], color=colors[i%4],ls=styles[(i-1)%3])
    pl = ax.plot(xarray, tabs[-1], label=labels[-1], color='k',ls='-')

    if (legflag):
        nlegcol = nLegCols(ntabs)
        plt.legend(loc='best',numpoints=1, ncol=nlegcol)
    plt.show()
    return fig


def Clean(xarray,yarray,data):
    fig = plt.figure(figsize=(11,8))
    ax = fig.add_subplot(111)
    ax.set_title('click to remove line segments')
    rf, qf = np.meshgrid(xarray,yarray)
    pl = ax.pcolormesh(rf,qf,data,vmin=-7,vmax=0,picker=True)
    xdiff = np.diff(xarray)
    ydiff = np.diff(yarray)
    ax.set_xlim(xarray[0]-xdiff[0],xarray[-1]+xdiff[-1])
    ax.set_ylim(yarray[0]-ydiff[0],yarray[-1]+ydiff[-1])
    cleaner = DataCleaner(pl,xarray,yarray)
    plt.show()
    return cleaner.Data()


class DataCleaner:
    def __init__(self, data, xarray, yarray):
        self.data = data
        self.xarray = xarray
        self.yarray = yarray
        self.ind = np.ones((len(yarray),len(xarray)))
        self.cid = data.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        if event.inaxes!=self.data.axes: return
        x = event.xdata
        y = event.ydata
        x_ind = bs.bisect(self.xarray,x)-1
        y_ind = bs.bisect(self.yarray,y)-1
        print x_ind, y_ind
        self.ind[y_ind,x_ind] = 0
        self.data.figure.canvas.draw()

    def Data(self):
        return self.ind

def step(ydata, xarray, xlab, ylab, title, xlim, ylim, log):
    """
    Small step plotter function
    step(ydata, xarray, xlab, ylab, title, xlim, ylim, log)
    """
    fig = plt.figure(figsize=(9,6))
    mpl.rc("font", family="serif")
    ax = fig.add_subplot(111)
    ax.set_xlim(xlim[0],xlim[1])
    ax.set_ylim(ylim[0],ylim[1])
    ax.step(xarray,ydata,'k-', linewidth=1.5)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    if ('x' in log):
        ax.set_xscale('log', nonposx='clip')
    if ('y' in log):
        ax.set_yscale('log', nonposy='clip')
    plt.show()
    return fig

def fbwt(x, y, ybtw, title='X vs Y', xlab='X', ylab='Y', labels=None,
         xlim=None, ylim=None, log=''):
    """
    Fill between plotter function
    fbwt(x,y,err,title,xlab,ylab,labels,log)
    """
    labels, xlim, ylim = none_to_empty_list(labels, xlim, ylim)
    # Limits
    xlo = INF; xhi = -INF; ylo = INF; yhi = -INF
    xlo = np.min([x[0][0],xlo])
    xhi = np.max([x[0][-1],xhi])
    xax = x[0].copy()
    # Important Numbers
    nY = len(y)
    # Plot
    fig = plt.figure(figsize=(9,6))
    mpl.rc("font", family="serif")
    ax = fig.add_subplot(111)
    for i in range(nY):
        if (len(x) > 1):
            xax = x[i].copy()
            xlo = np.min([xax[0],xlo])
            xhi = np.max([xax[-1],xhi])

        if (nY == 1):
            ax.plot(xax,y[i],'%s%s'%(colorsmall[i%4],styles[i%3]))
        else:
            ax.plot(xax,y[i],'%s%s'%(colorsmall[i%4],styles[i%3]),label=labels[i])
        error = ybtw[i]
        if (len(error)==1):
            ax.fill_between(xax, error, error, alpha=.5, edgecolor=colorsmall[i%4], facecolor=colorsmall[i%4], interpolate=True)
            ylo = np.min([np.min(error),ylo])
            yhi = np.max([np.max(error),yhi])
        elif (len(error)==2):
            ax.fill_between(xax, error[0], error[1], alpha=.5, edgecolor=colorsmall[i%4], facecolor=colorsmall[i%4], interpolate=True)
            ylo = np.min([np.min(error[0]),ylo])
            yhi = np.max([np.max(error[1]),yhi])
        else:
            print "Wrong number of uncertainties for error bars! Exiting..."
            import sys
            sys.exit(0)

    if (not xlim):
        xlim.append(xlo)
        xlim.append(xhi)
    if (not ylim):
        ylim.append(ylo)
        ylim.append(yhi)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    if ('x' in log):
        ax.set_xscale('log', nonposx='clip')
    if ('y' in log):
        ax.set_yscale('log', nonposy='clip')
    ax.grid(True)
    plt.legend(loc='best',numpoints=1)
    plt.show()
    return fig


def FILLBTW(X,Y,ERR,title,xlab,ylab,log):
    """
    Small fill between plotter function
    FILLBTW(X,Y,ERR,title,xlab,ylab,log)
    """
    # Plot
    fig = plt.figure(figsize=(9,6))
    mpl.rc("font", family="serif")
    ax = fig.add_subplot(111)
    ax.plot(X,Y,'k.')
    ax.set_xlim([10**2,10**6])
    ax.fill_between(X, Y-ERR, Y+ERR, alpha=.5, edgecolor='#3F7F4C', facecolor='#7EFF99', interpolate=True)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    if ('x' in log):
        ax.set_xscale('log', nonposx='clip')
    if ('y' in log):
        ax.set_yscale('log', nonposy='clip')
    ax.grid(True)
    plt.show()
    return fig
