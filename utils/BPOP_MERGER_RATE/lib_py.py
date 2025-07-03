import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kde
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
import scipy.ndimage
import matplotlib
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy.ndimage
from scipy import interpolate
from scipy import stats
from scipy.stats import kde
from numpy import linspace
from matplotlib.colors import LogNorm
from matplotlib import pyplot as plt
import math
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
import re
import random 
import csv
from scipy import optimize
from matplotlib.collections import PatchCollection
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.ticker as ticker
import os.path
from os import path
import matplotlib.colors as mcolors
#from scipy.interpolate import spline
#from scipy.interpolate import make_interp_spline, BSpline
from numpy import mean
import random
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp



plt.rc('text', usetex=True)
cmap = "viridis"
fsize = 41
lsize = 40
family = 'Times New Roman'
matplotlib.rcParams.update({'font.size': fsize})
matplotlib.rcParams.update({'font.family': family})




def biaser(m1, m2, rgw):

    for i in range(len(m1)):
        if(m1[i] < m2[i]):
            mh = m2[i]
            m2[i] = m1[i]
            m1[i] = mh




    #DO THE FOLLOWING TO GET BIASED SAMPLES
    #Primary Bias

    bias_m = "yes"
    bias_q = "yes"

    minbh = min(m1)
    maxbh = max(m2)
    mslope = 2.0
    if(bias_m == "no"):
        mslope = 0.0
    obsbhx = np.arange(minbh,maxbh,(maxbh-minbh)/100)
    obsbhy = pow(obsbhx,mslope)
    obsmax = max(obsbhy)
    
    redshift_horizon = 2.0


    id_bias=[]
    for i in range(len(m1)):
        if(rgw[i] < redshift_horizon):
            qslope = getq(m1[i],m2[i])
            if(bias_q == "no"):
                qslope = 0.0
            
            mbh_tmp = m1[i]
            fmbh = pow(mbh_tmp,mslope)/obsmax
            fmrej= random.random()
        

        
            if(fmrej <= fmbh):
                mns_tmp = m2[i]
                rat_tmp = mns_tmp / mbh_tmp
                if(rat_tmp > 1.):
                    print ("PLEASE PROVIDE M1 > M2")
                    exit()                
                
                fqns = pow(rat_tmp,qslope)
                fqrn = random.random()
                if(fqrn <= fqns):
                    id_bias.append(i)
                


    

	#print ("Selected binaries", len(id_bias))


    return id_bias



def test():
    zredmax   = 2.036
    psisfrmax = 0.01 * pow(1+zredmax,2.6) / (1. + pow((1+zredmax)/3.2,6.2))
    z = []
    
    for i in range(1000):
        while True:
            zred = 10.*np.random.random();
            psisfr = 0.01 * pow(1+zred,2.6) / (1. + pow((1+zred)/3.2,6.2));    
            psirnd = psisfrmax * np.random.random();
            
            if(psisfr > psirnd):
                break

        z.append(zred)
        
    return z

def tvsred(z):

    tH = 27.0164323345212 # Gyr
    alpha = 1.82262661377844
    scatter =  0.3 # Gyr;
    
    toff = tH / (1.+pow(1.+z, alpha)) + (-scatter + 2.*scatter*random.random());
    
    return toff

def redvst(t):
    
    tH = 27.0164323345212 # Gyr
    alpha = 1.82262661377844
    scatter =  0.03 # Gyr;

    if(t > tH):
        t = tH;
    
    zoff = pow(tH / t - 1, 1./alpha) - 1 + (-scatter + 2.*scatter*random.random());
        
    
    return zoff

def mapping(model, lims, x_iso, y_iso, x_clu, y_clu, xlab, ylab, colmap, filename):

    x = []
    y = []
    for ii in range(len(x_iso)):
        x.append(x_iso[ii])
        y.append(y_iso[ii])
    for ii in range(len(x_clu)):
        x.append(x_clu[ii])
        y.append(y_clu[ii])

    valx = [(1./len(x)) for i in range(len(x)) ]
    valy = [(1./len(y)) for i in range(len(y)) ]


    valxiso = [1./len(x_iso) for iii in range(len(x_iso)) ]
    valyiso = [1./len(y_iso) for iii in range(len(y_iso)) ]     
    valxclu = [1./len(x_clu) for iii in range(len(x_clu)) ]
    valyclu = [1./len(y_clu) for iii in range(len(y_clu)) ]

    lsize = 49
    fsize = 49
    family = 'Times New Roman'
    matplotlib.rcParams.update({'font.family': family})

    xlabel = xlab
    ylabel = ylab
    cmap   = cm.get_cmap(colmap)
    boxcol1 = blu3
    boxcol2 = blu3
    edgecol = dpur

    # Make the plot
    xmin = lims[0]
    xmax = 150 #lims[1]
    ymin = lims[2]
    ymax = lims[3]
 
    xlims = [xmin,xmax]
    ylims = [ymin,ymax]
    
    # Define the locations for the axes
    left0,right = 0.15,0.56
    left, width = 0.15, 0.56
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02
    
    # Set up the geometry of the three plot
    rect_temperature = [left0, bottom, right, height]   # dimensions of temp plot
    rect_histx       = [left, bottom_h, width, 0.25]    # dimensions of x-histogram
    rect_histy       = [left_h, bottom, 0.25, height]   # dimensions of y-histogram
    

    fig = plt.figure(1, figsize=(20,18))
    

    axTemperature = plt.axes(rect_temperature)
    axHistx       = plt.axes(rect_histx) # x histogram
    axHisty       = plt.axes(rect_histy) # y histogram
    
    #axTemperature.set_xscale('log')
    #axHistx.set_xscale('log')
    #axTemperature.set_yscale('log')
    #axHisty.set_yscale('log')
    
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    # Make the 'main' temperature plot
    
    # Define the number of bins
    nxbins = 25
    nybins = 25
    nbins  = 25
    
    xbins   = linspace(start = xmin, stop = xmax, num = nxbins)
    ybins   = linspace(start = ymin, stop = ymax, num = nybins)
    xcenter = (xbins[0:-1]+xbins[1:])/2.0
    ycenter = (ybins[0:-1]+ybins[1:])/2.0
    
    # Find the min/max of the data
    xmin = min(xlims) #xmin #min(xlims)
    xmax = max(xlims) #xmax #max(xlims)
    ymin = min(ylims) #ymin #min(ylims)
    ymax = max(ylims) #ymax #max(ylims)
    
    nbins = 25
    k = kde.gaussian_kde([x_iso,y_iso])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    #mesh = axTemperature.pcolormesh(xi, yi, zi.reshape(xi.shape),cmap=cmap)

    n = 500
    zi = zi/zi.sum()
    t = np.linspace(0, zi.max(), n)
    integral = ((zi >= t[:, None, None]) * zi).sum(axis=(1,2))

    t = t[::-1]
    integral = integral[::-1]

    f = interpolate.interp1d(integral, t)
    fname = 'Times New Roman'
    t_contours = f(np.array([0.99,0.7]))
    fmt = {}
    strs = ["99%","70%"]
    for l, s in zip(t_contours, strs):
        fmt[l] = s

    ccol = blu3
    origin = 'lower'
    CS = axTemperature.contour(xi, yi, zi.reshape(xi.shape), t_contours, colors=(ccol), linewidths=2.3, origin=origin)        
    axTemperature.clabel(CS, inline=1, fmt=fmt, fontsize=fsize)#, fontname=fname)
    
        
    nbins = 25
    k = kde.gaussian_kde([x_clu,y_clu])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    #mesh = axTemperature.pcolormesh(xi, yi, zi.reshape(xi.shape),cmap=cmap)

    n = 500
    zi = zi/zi.sum()
    t = np.linspace(0, zi.max(), n)
    integral = ((zi >= t[:, None, None]) * zi).sum(axis=(1,2))

    t = t[::-1]
    integral = integral[::-1]

    f = interpolate.interp1d(integral, t)
    fname = 'Times New Roman'
    t_contours = f(np.array([0.99,0.7]))
    fmt = {}
    strs = ["99%","70%"]
    for l, s in zip(t_contours, strs):
        fmt[l] = s

    ccol = "orange"
    origin = 'lower'
    CS = axTemperature.contour(xi, yi, zi.reshape(xi.shape), t_contours, colors=(ccol), linewidths=2.3, origin=origin)        
    axTemperature.clabel(CS, inline=1, fmt=fmt, fontsize=fsize)#, fontname=fname)


    #Plot the axes labels
    axTemperature.set_xlabel(xlabel,fontsize=lsize)#,fontname='Times New Roman')    
    axTemperature.set_ylabel(ylabel,fontsize=lsize)#,fontname='Times New Roman')
    
    #Make the tickmarks pretty
    ticklabels = axTemperature.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    
    ticklabels = axTemperature.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')

    #Make the tickmarks pretty
    
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')

    #Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    
    
    #Set up the plot limits
    axTemperature.set_xlim(xlims)
    axTemperature.set_ylim(ylims)
            
    
    #Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax-xmin)*1./nbins)
    ybins = np.arange(ymin, ymax, (ymax-ymin)*1./nbins)
    
    #Plot the histograms
    #axHistx.hist(x, bins=xbins,                           color = boxcol1, edgecolor=edgecol, weights=valx, histtype='stepfilled', linewidth=2,alpha=0.7)

    if(xlab == "$M$ [M$_\odot$]"):
        l1 = "iso"
        l2 = "dyn"
    else:
        l1 = "1st-Gen"
        l2 = "$n$th-Gen"
        
        
    axHistx.hist(x_iso, bins=xbins, range=[xmin,xmax], color=blu3, edgecolor="black", weights=valxiso, histtype='stepfilled', linewidth=4, linestyle="-",label=l1)
    axHistx.hist(x_clu, bins=xbins, range=[xmin,xmax], edgecolor="orange", weights=valxclu, histtype='step', linewidth=4, linestyle="-",label=l2)
    
    #axHisty.hist(y, bins=ybins, orientation='horizontal', color = boxcol2, weights=valy,edgecolor=edgecol, histtype='stepfilled', linewidth=2)

    axHisty.hist(y_iso, bins=ybins, range=[ymin,ymax], color=blu3, orientation='horizontal', weights=valyiso,edgecolor="black", histtype='stepfilled', linewidth=4, linestyle="-")
    axHisty.hist(y_clu, bins=ybins, range=[ymin,ymax], orientation='horizontal', edgecolor = "orange", weights=valyclu, histtype='step', linewidth=4, linestyle="-")

    axHistx.legend(loc="upper right")
    
    axHistx.set_xlim( xmin, xmax )
    axHisty.set_ylim( ymin, ymax )

    axTemperature.set_xlim( xmin, xmax)

    axHisty.xaxis.set_major_locator(MaxNLocator(3))
    axHistx.yaxis.set_major_locator(MaxNLocator(3))

    axHistx.set_yscale("log")
    axHisty.set_xscale("log")

    axHisty.xaxis.set_major_formatter(FormatStrFormatter('%1.0e'))
    axHistx.yaxis.set_major_formatter(FormatStrFormatter('%1.0e'))
    #ltit= " --- "
    #textstr1 = ltit
    #props = dict(boxstyle='round', facecolor='white', alpha=0.8)
    #axHistx.text(9.5, 1.0, textstr1, fontsize=fsize, color='black', fontname='Times New Roman', bbox=props)

    #axTemperature.text(x=0.8,y=0.1, s=model, fontsize=44, 
    #                   bbox = props, horizontalalignment='center',
    #                   verticalalignment='center',
    #                   transform = axTemperature.transAxes)

    axTemperature.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='off',direction="in", pad=10)
    axHistx.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='off',direction="in", pad=10)
    axHisty.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='off',direction="in", pad=10)
    
    if(xlab == "$M_\mathrm{rem}$ [M$_\odot$]"):
        fname = "/home/manuel/Scrivania/ACTIVE_PROJECTS/LIGO_DATA/LIGO_CAT"
        Mfi_O2  = reading(fname,13,"yes")
        MfiL_O2 = reading(fname,14,"yes")
        MfiU_O2 = reading(fname,15,"yes")
        
        Sfi_O2  = reading(fname,16,"yes")
        SfiL_O2 = reading(fname,17,"yes")
        SfiU_O2 = reading(fname,18,"yes")
        
        axTemperature.errorbar(Mfi_O2,Sfi_O2,xerr=[MfiL_O2,MfiU_O2], yerr=[SfiL_O2,SfiU_O2], color="black", fmt='.', capsize=10, elinewidth=2.8,markeredgewidth=1.5, alpha=0.6)
        axTemperature.scatter(Mfi_O2,Sfi_O2,marker="o",s=500,edgecolor="black",color="white",label="GWTC-1")


        xP1 = [142.]
        yP1 = [0.72]

        exPs=[16.]
        exPd=[28.]

        eyPs=[0.12]
        eyPd=[0.06]

        axTemperature.errorbar(xP1,yP1,xerr=[exPs,exPd], yerr=[eyPs,eyPd], color="black", fmt='.', capsize=10, elinewidth=2.8,markeredgewidth=1.5, alpha=0.6)
        axTemperature.scatter(xP1,yP1,marker="o",s=500,edgecolor="black",color="white",label="GW190521")


        
    plt.savefig(filename+'.jpeg')
    plt.savefig(filename+'.pdf')
    plt.close()



def creating_plot(model, lims, x, y, x1, y1, exd, exs, eyd, eys, xlab, ylab, colmap, filename,env):
    
    valx = [(1./len(x)) for i in range(len(x)) ]
    valy = [(1./len(y)) for i in range(len(y)) ]
    valx1 = [(.2/len(x1)) for i in range(len(x1)) ]
    valy1 = [(.2/len(y1)) for i in range(len(y1)) ]


    x_clu = []
    x_iso = []
    y_clu = []
    y_iso = []

    for i in range(len(x)):
        if(env[i] == "none"):
            x_iso.append(x[i])
            y_iso.append(y[i])
        elif(env[i] == "young" or env[i] == "globular" or env[i] == "nuclear"):
            x_clu.append(x[i])
            y_clu.append(y[i])
        else:
            print ("No environs found ..retry")
            exit()

    valxiso = [1./len(x_iso) for iii in range(len(x_iso)) ]
    valyiso = [1./len(y_iso) for iii in range(len(y_iso)) ]     
    valxclu = [1./len(x_clu) for iii in range(len(x_clu)) ]
    valyclu = [1./len(y_clu) for iii in range(len(y_clu)) ]
    
    lsize = 49
    fsize = 49
    family = 'Times New Roman'
    matplotlib.rcParams.update({'font.family': family})

    xlabel = xlab
    ylabel = ylab
    cmap   = cm.get_cmap(colmap)
    boxcol1 = blu3
    boxcol2 = blu3
    edgecol = dpur

    # Make the plot
    xmin = lims[0]
    xmax = lims[1]
    ymin = lims[2]
    ymax = lims[3]
 
    xlims = [xmin,xmax]
    ylims = [ymin,ymax]
    
    
    nbins = 100
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    
    # Define the locations for the axes
    left0,right = 0.12,0.55
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02
    
    # Set up the geometry of the three plot
    rect_temperature = [left0, bottom, right, height]   # dimensions of temp plot
    rect_histx       = [left, bottom_h, width, 0.25]    # dimensions of x-histogram
    rect_histy       = [left_h, bottom, 0.25, height]   # dimensions of y-histogram
    

    fig = plt.figure(1, figsize=(20,18))
    

    axTemperature = plt.axes(rect_temperature)
    axHistx       = plt.axes(rect_histx) # x histogram
    axHisty       = plt.axes(rect_histy) # y histogram
    
    #axTemperature.set_xscale('log')
    #axHistx.set_xscale('log')
    #axTemperature.set_yscale('log')
    #axHisty.set_yscale('log')
    
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    # Make the 'main' temperature plot
    
    # Define the number of bins
    nxbins = 20
    nybins = 20
    nbins  = 50
    
    xbins   = linspace(start = xmin, stop = xmax, num = nxbins)
    ybins   = linspace(start = ymin, stop = ymax, num = nybins)
    xcenter = (xbins[0:-1]+xbins[1:])/2.0
    ycenter = (ybins[0:-1]+ybins[1:])/2.0

    xbins1   = linspace(start = xmin, stop = xmax, num = 15) #xbins #linspace(start = min(x1), stop = max(x1), num = nxbins)
    ybins1   = linspace(start = ymin, stop = ymax, num = 15) #ybins #linspace(start = min(y1), stop = max(y1), num = nybins)

    
    # Find the min/max of the data
    xmin = min(xlims) #xmin #min(xlims)
    xmax = max(xlims) #xmax #max(xlims)
    ymin = min(ylims) #ymin #min(ylims)
    ymax = max(ylims) #ymax #max(ylims)
    
    
    mesh = axTemperature.pcolormesh(xi, yi, zi.reshape(xi.shape),cmap=cmap)
    n = 1000
    zi = zi/zi.sum()
    t = np.linspace(0, zi.max(), n)
    integral = ((zi >= t[:, None, None]) * zi).sum(axis=(1,2))


    t = t[::-1]
    integral = integral[::-1]

    f = interpolate.interp1d(integral, t)
    fname = 'Times New Roman'
    t_contours = f(np.array([0.99999999,0.99,0.9,0.7]))
    fmt = {}
    strs = ["","99%","90%","70%"]
    for l, s in zip(t_contours, strs):
        fmt[l] = s


    ccol = 'white'
    origin = 'lower'
    CS = axTemperature.contour(xi, yi, zi.reshape(xi.shape), t_contours, colors=(ccol), linewidths=2.3, origin=origin)        
    axTemperature.clabel(CS, inline=1, fmt=fmt, fontsize=fsize)#, fontname=fname)


    pcol = []
    
    f = []
    for io in range(len(x1)):
        if(io < len(x1)-8):
            f.append("o")
            pcol.append("white")
        else:
            f.append("s")
            pcol.append("white")

    sss = 25
    s = [sss for io in range(len(x1))]
    z = np.arange(len(x1))
    
    #axTemperature.errorbar(x1, y1, xerr=[exs,exd], yerr=[eys,eyd], ecolor=pcol,fmt='.', capsize=10, elinewidth=2.8,markeredgewidth=1.5, alpha=0.6)
    for i in z:
        axTemperature.errorbar(x1[i], y1[i], xerr=[[exs[i]],[exd[i]]], yerr=[[eys[i]],[eyd[i]]], ecolor=pcol[i], capsize=10, elinewidth=2.8, markeredgewidth=1.5, alpha=0.6, marker=f[i],markersize=s[i], markerfacecolor=pcol[i], markeredgecolor="black")        
        #axTemperature.scatter(x1[i], y1[i], marker=f[i],s=s[i],color=pcol[i],linewidths=1.5, edgecolors='black', alpha=0.6)    

    axTemperature.clabel(CS, inline=1, fmt=fmt, fontsize=fsize)#, fontname=fname)
    
    #LVC CONTOUR PLOTS 
    '''
    nbins_lvc= 8
    xmin_lvc = xmin
    xmax_lvc = xmax
    ymin_lvc = ymin
    ymax_lvc = ymax
    k_lvc = kde.gaussian_kde([x1,y1])
    xi_lvc, yi_lvc = np.mgrid[xmin_lvc:xmax_lvc:nbins_lvc*1j, ymin_lvc:ymax_lvc:nbins_lvc*1j]
    zi_lvc = k_lvc(np.vstack([xi_lvc.flatten(), yi_lvc.flatten()]))    
    n = 1000
    zi_lvc = zi_lvc/zi_lvc.sum()
    t = np.linspace(0, zi_lvc.max(), n)
    integral = ((zi_lvc >= t[:, None, None]) * zi_lvc).sum(axis=(1,2))
    t = t[::-1]
    integral = integral[::-1]
    f = interpolate.interp1d(integral, t)
    fname = 'Times New Roman'
    t_contours = f(np.array([0.99,0.9,0.7]))
    fmt = {}
    strs = ["99%","90%","70%"]
    for l, s in zip(t_contours, strs):
        fmt[l] = s
    CS = axTemperature.contour(xi_lvc, yi_lvc, zi_lvc.reshape(xi_lvc.shape), t_contours, colors="red", linewidths=2.3, origin=origin)        
    #axTemperature.clabel(CS, inline=1, fmt=fmt, fontsize=fsize, fontname=fname)
    '''
    
        
    #Plot the axes labels
    axTemperature.set_xlabel(xlabel,fontsize=lsize,fontname='Times New Roman')  
    axTemperature.set_ylabel(ylabel,fontsize=lsize,fontname='Times New Roman')
    
    #Make the tickmarks pretty
    ticklabels = axTemperature.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    
    ticklabels = axTemperature.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')

    #Make the tickmarks pretty
    
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')

    #Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    
    
    #Set up the plot limits
    axTemperature.set_xlim(xlims)
    axTemperature.set_ylim(ylims)
            
    
    #Set up the histogram bins
    #print ymin, ymax, (ymax-ymin)*1./nbins
    
    xbins = np.arange(xmin, xmax, (xmax-xmin)*1./nbins)
    ybins = np.arange(ymin, ymax, (ymax-ymin)*1./nbins)
    
    #Plot the histograms
    axHistx.hist(x, bins=xbins,                           color = boxcol1, edgecolor=edgecol, weights=valx, histtype='stepfilled', linewidth=2,alpha=0.7)
    axHistx.hist(x1,bins=xbins1,                           color = boxcol1, edgecolor="red", weights=valx1, histtype='step', linewidth=4, linestyle=":", label="LVC")



    axHistx.hist(x_iso, bins=xbins, range=[xmin,xmax], edgecolor=dpur, weights=valxiso, histtype='step', linewidth=4, linestyle="-",label="isolated")
    axHistx.hist(x_clu, bins=xbins, range=[xmin,xmax], edgecolor=dpur, weights=valxclu, histtype='step', linewidth=4, linestyle="--",label="dynamical")

    
    
    axHisty.hist(y, bins=ybins, orientation='horizontal', color = boxcol2, weights=valy,edgecolor=edgecol, histtype='stepfilled', linewidth=2)
    axHisty.hist(y1,bins=ybins1, orientation='horizontal', color = boxcol2, weights=valy1,edgecolor="red", histtype='step', linewidth=4, linestyle=":", label="LVC")

    axHisty.hist(y_iso, bins=ybins, range=[ymin,ymax], orientation='horizontal', color = dpur, weights=valyiso,edgecolor=dpur, histtype='step', linewidth=4, linestyle="-")
    axHisty.hist(y_clu, bins=ybins, range=[ymin,ymax], orientation='horizontal', color = dpur, weights=valyclu,edgecolor=dpur, histtype='step', linewidth=4, linestyle="--")

    axHistx.legend(loc="upper right",bbox_to_anchor=(1.5,1.))
    
 
    #Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(3))
    axHistx.yaxis.set_major_locator(MaxNLocator(3))    
    axHisty.xaxis.set_major_formatter(FormatStrFormatter('%1.0e'))
    axHistx.yaxis.set_major_formatter(FormatStrFormatter('%1.0e'))

    ltit= " --- "
    textstr1 = ltit
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', alpha=0.8)
    # place a text box in upper left in axes coords
    #axHistx.text(9.5, 1.0, textstr1, fontsize=fsize, color='black', fontname='Times New Roman', bbox=props)
    #    axHistx.set_ylim()
    #    axHisty.set_xlim(0,2.5)

    #axHistx.yaxis.set_ticks(np.arange(0.01, 0.05, 0.03))
    #axHisty.xaxis.set_ticks(np.arange(0.03, 0.09, 0.05))

    
    #Show the plot
    
    #plt.draw()
    
    
    #plt.colorbar()

    #axTemperature.legend(bbox_to_anchor=(0.8, 0.95),bbox_transform=plt.gcf().transFigure, title = model, prop={"size":52})
    #axTemperature.legend(loc="lower right", title=str(model), framealpha=0.7)


    axHistx.set_yscale("log")
    axHisty.set_xscale("log")
    
    axHistx.set_ylim(1.e-3,.7e0)
    axHisty.set_xlim(1.e-3,.7e0)


    axTemperature.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='on',direction="in", pad=10)
    axHistx.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='on',direction="in", pad=10)
    axHisty.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='on',direction="in", pad=10)

    xpos = 0.1
    ypos = 0.9
    if(ylab == "$\chi_\mathrm{eff}$"):
        xpos = 0.9
        ypos = 0.1
    
    axTemperature.text(x=xpos,y=ypos, s=model, fontsize=44, 
                       bbox = props, horizontalalignment='center',
                       verticalalignment='center',
                       transform = axTemperature.transAxes)

    #Set up the histogram limits
    axHistx.set_xlim( lims[0], lims[1] )
    axHisty.set_ylim( lims[2], lims[3] )
    axTemperature.set_xlim( lims[0], lims[1])
    axTemperature.set_ylim( lims[2], lims[3] )
    axTemperature.xaxis.set_major_locator(MaxNLocator(7))
    axTemperature.yaxis.set_major_locator(MaxNLocator(7))    

    
    
    mesh.set_rasterized(True) 

    plt.savefig(filename+'.jpeg')
    plt.savefig(filename+'.pdf')
    plt.close()






#Mass ratio Bias
def getq(ma,mb):
    m_ratio = [10,20,25,30,50]
    a_ratio = [0.47,0.53,0.56,0.59,0.72]
    ur = 0
    lr = len(m_ratio)-1
    
    for j in range(len(m_ratio)):
        if(ma > m_ratio[j]):
            ur = j
        if(ma <= m_ratio[j]):
            lr = j
            break

    if(ur == len(m_ratio)-1):
        lr = ur
        ur = ur-1
        
    if(lr == 0):
        ur = lr
        lr = lr+1
            
    if(ur == lr):
        print ("Warning!! division by zero",ur,lr)
            
    qslope = (a_ratio[ur]-a_ratio[lr])/(m_ratio[ur]-m_ratio[lr])*(ma - m_ratio[ur]) + a_ratio[ur]

    #print ma, mb, qslope, ur, lr
        
    return qslope
    


def count(mbh,mmax,mmin):
    cntbh = 0
    for i in range(len(mbh)):
        if(mbh[i] < mmax and mbh[i] > mmin):
          cntbh+=1
            
                
    return cntbh
            

def weighted_map(newbhmL,newqhmL,newbhmH,newqhmH,mgw1,qgw):
    # Define the locations for the axes
    left0,right = 0.15,0.55
    left, width = 0.15, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02
    
    # Set up the geometry of the three plot
    rect_temperature = [left0, bottom, right, height]   # dimensions of temp plot
    rect_histx       = [left, bottom_h, width, 0.25]    # dimensions of x-histogram
    rect_histy       = [left_h, bottom, 0.25, height]   # dimensions of y-histogram


    
    fig = plt.figure(1, figsize=(15,15))

    
    ax      = plt.axes(rect_temperature)
    axHistx = plt.axes(rect_histx) # x histogram
    axHisty = plt.axes(rect_histy) # y histogram

    
    xlabel = "$M_\mathrm{BH}$ [M$_\odot$]"
    ylabel = "$q$ [mass ratio]"


    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    #Plot the axes labels
    ax.set_xlabel(xlabel,fontsize=lsize,fontname='Times New Roman')    
    ax.set_ylabel(ylabel,fontsize=lsize,fontname='Times New Roman')
    
    #Make the tickmarks pretty
    ticklabels = ax.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    
    ticklabels = ax.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
        
    #Make the tickmarks pretty

    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
 
    #Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    

    ymin = 0.02
    ymax = 1
    xmin = 0.5
    xmax = 60.
    
    xlims = [xmin,xmax]
    ylims = [ymin,ymax]
    
    #Set up the plot limits
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    
    axHistx.set_xlim(xlims)
    axHisty.set_ylim(ylims)
    axHisty.set_yscale("log")
    
    #Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))
    
    axHisty.xaxis.set_major_formatter(FormatStrFormatter('%1.1f'))
    axHistx.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))

    axHisty.set_yticklabels("")
    
    
    ax.scatter(newbhmL,newqhmL,marker="o",s=400,c=gre2,alpha=0.2, edgecolor="black", label="Dynamical (metal-poor)")
    ax.scatter(newbhmH,newqhmH,marker="o",s=400,c=blu2,alpha=0.2, edgecolor="black", label="Dynamical (metal-rich)")

    # Make the plot

    ''''
    x = newbhmL
    y = newqhmL
    xcmp = newbhmH
    ycmp = newqhmH    
    nbins = 200
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    n = 100
    zi = zi/zi.sum()
    t = np.linspace(0, zi.max(), n)
    integral = ((zi >= t[:, None, None]) * zi).sum(axis=(1,2))
    t = t[::-1]
    integral = integral[::-1]
    f = interpolate.interp1d(integral, t)
    fname = 'Times New Roman'
    #t_contours = f(np.array([0.99,0.9,0.8,0.7,0.5,0.2]))
    t_contours = f(np.array([0.99,0.5,0.2]))
    fmt = {}
    #strs = ["99%","90%","80%","70%", "50%", "20%"]
    strs = ["99%", "50%", "20%"]
    for l, s in zip(t_contours, strs):
        fmt[l] = s    
    origin = 'lower'
    CS = ax.contour(xi, yi, zi.reshape(xi.shape), t_contours, colors=(gre2), linewidths=2, origin=origin)        
    ax.clabel(CS, inline=1, fmt=fmt, fontsize=fsize, fontname=fname)

    
    nbins = 200
    k = kde.gaussian_kde([xcmp,ycmp])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    n = 100
    zi = zi/zi.sum()
    t = np.linspace(0, zi.max(), n)
    integral = ((zi >= t[:, None, None]) * zi).sum(axis=(1,2))
    t = t[::-1]
    integral = integral[::-1]
    f = interpolate.interp1d(integral, t)
    fname = 'Times New Roman'
    #t_contours = f(np.array([0.99,0.9,0.8,0.7,0.5,0.2]))
    t_contours = f(np.array([0.99,0.5,0.2]))
    fmt = {}
    #strs = ["99%","90%","80%","70%", "50%", "20%"]
    strs = ["99%", "50%", "20%"]
    for l, s in zip(t_contours, strs):
        fmt[l] = s    
    origin = 'lower'
    CS = ax.contour(xi, yi, zi.reshape(xi.shape), t_contours, colors=(blu2), linewidths=2, origin=origin)        
    ax.clabel(CS, inline=1, fmt=fmt, fontsize=fsize, fontname=fname)
    '''

    
    value = []
    for i in range(len(newbhmH)):
        value.append(1./len(newbhmH))
       
    axHistx.hist(newbhmH,15,weights=value,color=blu2,alpha = 0.7, range=(xmin,xmax))



    value = []
    for i in range(len(newqhmH)):
        value.append(1./len(newqhmH))

    axHisty.hist(newqhmH,bins=np.logspace(np.log10(ymin),np.log10(ymax), 20),weights=value,orientation='horizontal',color=blu2,alpha = 0.7, range=(ymin,ymax))


    value = []
    for i in range(len(newbhmL)):
        value.append(1./len(newbhmL))
    

    axHistx.hist(newbhmL,15,weights=value,color=gre2,alpha = 0.7, range=(xmin,xmax))
        

    value = []
    for i in range(len(newqhmL)):
        value.append(1./len(newqhmL))

    axHisty.hist(newqhmL,bins=np.logspace(np.log10(ymin),np.log10(ymax), 20),weights=value,orientation='horizontal',color=gre2,alpha = 0.7, range=(ymin,ymax), rasterized=True)

    axHistx.set_xlim(xlims)

    rect = patches.Rectangle((mgw1*0.7,qgw*0.7),0.6*mgw1,0.6*qgw,linewidth=3,edgecolor="black",facecolor='none')
    ax.add_patch(rect)

    ax.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='on',direction="in", pad=10)
    axHistx.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='on',direction="in", pad=10)
    axHisty.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='on',direction="in", pad=10)
    

    axHisty.hlines(0.12,0.0,0.6,lw=3,color="black")
    axHisty.hlines(0.112,0.0,0.6,lw=3,linestyle="--",color="black")
    axHisty.hlines(0.103,0.0,0.6,lw=3,color="black")
    axHisty.set_xlim(0,0.55)
    
    axHistx.axvline(24.3,lw=3,color="black")
    axHistx.axvline(23.2,lw=3,linestyle="--",color="black")
    axHistx.axvline(22.2,lw=3,color="black")
    
    #axHistx.set_ylim(0,0.95)
    
    ax.set_yscale("log")

    #ax.set_ylim(0.0,0.6)
    #axHisty.set_ylim(0.0,0.6)
    
    
    mgw2 = qgw*mgw1
    
    ax.scatter(mgw1,qgw,marker="*",s=500,c="white",edgecolor="black",alpha=1,rasterized=True)
    ax.text(32.,0.1, 'GW190814')
    epsdx_q = np.sqrt(pow(1./mgw1 * 0.09, 2.) + pow(mgw2/(mgw1*mgw1) * 1.0,2.))
    epssx_q = np.sqrt(pow(1./mgw1 * 0.08, 2.) + pow(mgw2/(mgw1*mgw1) * 1.1,2.))


    mn1 = [5.7,8.9,5.7]
    edx1= [2.3,1.2,2.1]
    esx1= [4.0,1.5,1.8]
    mn2 = [1.5,1.9,1.5]
    edx2= [0.5,0.2,0.3]       
    esx2= [0.8,0.3,0.7]
    
    sname=['GW190426','GW200105','GW200115']
    
    for kkk in range(len(mn1)):
        qn1 = mn2[kkk] / mn1[kkk]
        qdx = np.sqrt(pow(1./mn1[kkk] * edx2[kkk],2.) + pow(mn2[kkk]/(mn1[kkk]*mn1[kkk]) * edx1[kkk],2.))
        qsx = np.sqrt(pow(1./mn1[kkk] * esx2[kkk],2.) + pow(mn2[kkk]/(mn1[kkk]*mn1[kkk]) * esx1[kkk],2.))
    
        axHisty.hlines(qn1-qdx,0.0,0.6,lw=3,color="black")
        axHisty.hlines(qn1,0.0,0.6,lw=3,linestyle="--",color="black")
        axHisty.hlines(qn1+qsx,0.0,0.6,lw=3,color="black")
        
        axHistx.axvline(mn1[kkk]-edx1[kkk],lw=3,color="black")
        axHistx.axvline(mn1[kkk],lw=3,linestyle="--",color="black")
        axHistx.axvline(mn1[kkk]+esx1[kkk],lw=3,color="black")
    
        ax.scatter(mn1[kkk],qn1,marker="*",s=2000,c="white",edgecolor="black",alpha=1)
        #rect = patches.Rectangle((mn1-edx1,qn1-qdx),mn1+abs(edx1-edx2),qn1+abs(qdx-qsx),linewidth=3,edgecolor="black",facecolor='none',linestyle="--")
        #ax.add_patch(rect)
        rect = patches.Rectangle((0.7*mn1[kkk],0.7*qn1),0.6*mn1[kkk],0.6*qn1,linewidth=3,edgecolor="black",facecolor='none')
        ax.add_patch(rect)
        ax.text(15,0.25+0.1*(len(mn1)-kkk), sname[kkk])
        
    
    


    
    #rect = patches.Rectangle( (mgw1-1,qgw-epssx_q), 1.+1.1, epssx_q+epsdx_q , linewidth=3,edgecolor="black",facecolor='none')
    #ax.add_patch(rect)

    
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    t = "volume weighted"
    leg = ax.legend(loc="upper right",title = t,  framealpha=0.0, prop={"size":"30"},markerfirst=False)
    plt.setp(leg.get_title(),fontsize='32')
    
    plt.savefig("mratioBHNS_weighted_new.jpeg")
    plt.savefig("mratioBHNS_weighted_new.pdf")
    plt.close()



def srt(list1,list2):
    zipped_lists = zip(list1, list2)
    sorted_pairs = sorted(zipped_lists)

    tuples = zip(*sorted_pairs)
    list1, list2 = [ list(tuple) for tuple in  tuples]

    xold = list1
    yold = list2

    return xold,yold



def colors():
    dpur = '#440154' # dark purple
    purp = '#472c7a' # purple
    blu2 = '#3b518b' # blue
    blu3 = '#2c718e' # blue
    blgr = '#21908d' # blue-green
    gre1 = '#27ad81' # green
    gre2 = '#5cc863' # green
    lgre = '#aadc32' # lime green
    yell = '#fde725' # yellow

    return dpur,purp,blu2,blu3,blgr,gre1,gre2,lgre,yell
    

dpur,purp,blu2,blu3,blgr,gre1,gre2,lgre,yell = colors()
colors = [dpur, blu2, blu3, blgr, gre1 , gre2, yell]
markers= ["o" , "^" , "s" , "*" , "p" , "H" , "D" ]
ls = ["-","--",":","-.","-","--",":","-.","-"]

def simple_map_bella(x,y,lims,xlab,ylab, colmap, filename):

    mx=lims[0]
    Mx=lims[1]
    my=lims[2]
    My=lims[3]
            
    fsize = 43
    lsize = 43
    family = 'Times New Roman'
    matplotlib.rcParams.update({'font.size': fsize})
    matplotlib.rcParams.update({'font.family': family})
    fig   = plt.figure(1, figsize=(18.0,15.0))    

 
    xlabel = xlab
    ylabel = ylab
    cmap   = cm.get_cmap(colmap)
    boxcol1 = 'white'
    boxcol2 = 'white'
    edgecol = '#440154FF'

    # Make the plot
    ymin =  my
    ymax =  My
    xmin =  mx
    xmax =  Mx
    
    xlims = [xmin,xmax]
    ylims = [ymin,ymax]
    
    
    nbins = 200
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    grid_coords = np.append(xi.reshape(-1,1),yi.reshape(-1,1),axis=1)
    
    
    # Define the locations for the axes
    left0,right = 0.12,0.55
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02
    
    # Set up the geometry of the three plot
    rect_temperature = [left0, bottom, right, height]   # dimensions of temp plot
    rect_histx       = [left, bottom_h, width, 0.25]    # dimensions of x-histogram
    rect_histy       = [left_h, bottom, 0.25, height]   # dimensions of y-histogram
    

    fig = plt.figure(1, figsize=(15,15))
    

    axTemperature = plt.axes(rect_temperature)
    axHistx       = plt.axes(rect_histx) # x histogram
    axHisty       = plt.axes(rect_histy) # y histogram
    
    #axTemperature.set_xscale('log')
    #axHistx.set_xscale('log')
    #axTemperature.set_yscale('log')
    #axHisty.set_yscale('log')
    
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    # Make the 'main' temperature plot
    
    # Define the number of bins
    nxbins = 20
    nybins = 20
    nbins = 50
    
    xbins   = linspace(start = xmin, stop = xmax, num = nxbins)
    ybins   = linspace(start = ymin, stop = ymax, num = nybins)
    xcenter = (xbins[0:-1]+xbins[1:])/2.0
    ycenter = (ybins[0:-1]+ybins[1:])/2.0
    
    # Find the min/max of the data
    xmin = min(xlims) #xmin #min(xlims)
    xmax = max(xlims) #xmax #max(xlims)
    ymin = min(ylims) #ymin #min(ylims)
    ymax = max(ylims) #ymax #max(ylims)
    
    H, yedges,xedges = np.histogram2d(y,x,bins=(ybins,xbins))
    
    X = xcenter
    Y = ycenter
    Z = H 
    
    #axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax],interpolation='nearest', origin='lower',aspect='auto')

    axTemperature.pcolormesh(xi, yi, zi.reshape(xi.shape),cmap=cmap)
    
    #transform zi to a 0-1 range
    Z = Z = (Z - Z.min())/(Z.max() - Z.min())
    
    n = 1000
    zi = zi/zi.sum()
    #zi = (zi - zi.min())/(zi.max() - zi.min())
    t = np.linspace(0, zi.max(), n)
    integral = ((zi >= t[:, None, None]) * zi).sum(axis=(1,2))


    t = t[::-1]
    integral = integral[::-1]

    f = interpolate.interp1d(integral, t)
    fname = 'Times New Roman'
    #t_contours = f(np.array([0.99,0.9,0.8,0.7,0.5,0.2]))
    t_contours = f(np.array([0.99,0.5,0.2]))
    fmt = {}
    #strs = ["99%","90%","80%","70%", "50%", "20%"]
    strs = ["99%", "50%", "20%"]
    for l, s in zip(t_contours, strs):
        fmt[l] = s
        
    origin = 'lower'
    CS = axTemperature.contour(xi, yi, zi.reshape(xi.shape), t_contours, colors=('white'), linewidths=2, origin=origin)        
    axTemperature.clabel(CS, inline=1, fmt=fmt, fontsize=fsize, fontname=fname)

    
    #Plot the axes labels
    axTemperature.set_xlabel(xlabel,fontsize=lsize,fontname='Times New Roman')    
    axTemperature.set_ylabel(ylabel,fontsize=lsize,fontname='Times New Roman')
    
    #Make the tickmarks pretty
    ticklabels = axTemperature.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    
    ticklabels = axTemperature.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')

    #Make the tickmarks pretty
    
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
 
    #Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    
    
    #Set up the plot limits
    axTemperature.set_xlim(xlims)
    axTemperature.set_ylim(ylims)
            
    
    #Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
    ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)
    
    #Plot the histograms
    axHistx.hist(x, bins=xbins, color = boxcol1, edgecolor=edgecol, density=1, histtype='step', linewidth=2)
    #axHistx.hist(x, bins=np.logspace(np.log10(xmin),np.log10(xmax), 50), color = boxcol1, density=0)

    
    axHisty.hist(y, bins=ybins, orientation='horizontal', color = boxcol2, edgecolor=edgecol, density=1, histtype='step', linewidth=2)
    #axHisty.hist(y, bins=np.logspace(np.log10(ymin),np.log10(ymax), 50), orientation='horizontal', color = boxcol2, density=0)
    
    #Set up the histogram limits
    axHistx.set_xlim( xmin, xmax )
    axHisty.set_ylim( ymin, ymax )
    
 
    #Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))
    
    axHisty.xaxis.set_major_formatter(FormatStrFormatter('%1.1e'))
    axHistx.yaxis.set_major_formatter(FormatStrFormatter('%1.1e'))
    ltit= "--"
    textstr1 = ltit
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', alpha=0.8)
    # place a text box in upper left in axes coords
    axHistx.text(9.5, 1.0, textstr1, fontsize=fsize, color='black', fontname='Times New Roman', bbox=props)
    #axHistx.set_ylim(0,1.2)
    #axHisty.set_xlim(0,2.5)

    plt.savefig(filename+'.jpeg')
    plt.savefig(filename+'.pdf')
    plt.close()

    return xi,yi,zi,grid_coords.T
    
def simple_map_bella2(x,y,lims,xcmp,ycmp,xlab,ylab, colmap, filename):

    mx=lims[0]
    Mx=lims[1]
    my=lims[2]
    My=lims[3]
  
    fsize = 43
    lsize = 43
    family = 'Times New Roman'
    matplotlib.rcParams.update({'font.size': fsize})
    matplotlib.rcParams.update({'font.family': family})
    fig   = plt.figure(1, figsize=(18.0,15.0))    

 
    xlabel = xlab
    ylabel = ylab
    cmap   = cm.get_cmap(colmap)
    boxcol1 = 'white'
    boxcol2 = 'white'
    edgecol = '#440154FF'

    # Make the plot
    xmin =  mx
    xmax =  Mx
    ymin =  my
    ymax =  My

    xlims = [xmin,xmax]
    ylims = [ymin,ymax]
    
    nbins = 200
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    grid_coords = np.append(xi.reshape(-1,1),yi.reshape(-1,1),axis=1)
    
    
    # Define the locations for the axes
    left0,right = 0.12,0.55
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02
    
    # Set up the geometry of the three plot
    rect_temperature = [left0, bottom, right, height]   # dimensions of temp plot
    rect_histx       = [left, bottom_h, width, 0.25]    # dimensions of x-histogram
    rect_histy       = [left_h, bottom, 0.25, height]   # dimensions of y-histogram
    

    fig = plt.figure(1, figsize=(15,15))
    

    axTemperature = plt.axes(rect_temperature)
    axHistx       = plt.axes(rect_histx) # x histogram
    axHisty       = plt.axes(rect_histy) # y histogram
    
    #axTemperature.set_xscale('log')
    #axHistx.set_xscale('log')
    #axTemperature.set_yscale('log')
    #axHisty.set_yscale('log')
    
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    # Make the 'main' temperature plot
    
    # Define the number of bins
    nxbins = 20
    nybins = 20
    nbins = 50
    
    xbins   = linspace(start = xmin, stop = xmax, num = nxbins)
    ybins   = linspace(start = ymin, stop = ymax, num = nybins)
    xcenter = (xbins[0:-1]+xbins[1:])/2.0
    ycenter = (ybins[0:-1]+ybins[1:])/2.0
    
    # Find the min/max of the data
    xmin = min(xlims) #xmin #min(xlims)
    xmax = max(xlims) #xmax #max(xlims)
    ymin = min(ylims) #ymin #min(ylims)
    ymax = max(ylims) #ymax #max(ylims)
    
    H, yedges,xedges = np.histogram2d(y,x,bins=(ybins,xbins))
    
    X = xcenter
    Y = ycenter
    Z = H 
    
    #axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax],interpolation='nearest', origin='lower',aspect='auto')

    axTemperature.pcolormesh(xi, yi, zi.reshape(xi.shape),cmap=cmap)
    
    #transform zi to a 0-1 range
    Z = Z = (Z - Z.min())/(Z.max() - Z.min())
    
    n = 1000
    zi = zi/zi.sum()
    #zi = (zi - zi.min())/(zi.max() - zi.min())
    t = np.linspace(0, zi.max(), n)
    integral = ((zi >= t[:, None, None]) * zi).sum(axis=(1,2))


    t = t[::-1]
    integral = integral[::-1]

    f = interpolate.interp1d(integral, t)
    fname = 'Times New Roman'
    #t_contours = f(np.array([0.99,0.9,0.8,0.7,0.5,0.2]))
    t_contours = f(np.array([0.99,0.5,0.2]))
    fmt = {}
    #strs = ["99%","90%","80%","70%", "50%", "20%"]
    strs = ["99%", "50%", "20%"]
    for l, s in zip(t_contours, strs):
        fmt[l] = s
        
    origin = 'lower'
    CS = axTemperature.contour(xi, yi, zi.reshape(xi.shape), t_contours, colors=('white'), linewidths=2, origin=origin)        
    axTemperature.clabel(CS, inline=1, fmt=fmt, fontsize=fsize, fontname=fname)

    
    #Plot the axes labels
    axTemperature.set_xlabel(xlabel,fontsize=lsize,fontname='Times New Roman')    
    axTemperature.set_ylabel(ylabel,fontsize=lsize,fontname='Times New Roman')
    
    #Make the tickmarks pretty
    ticklabels = axTemperature.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    
    ticklabels = axTemperature.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')

    #Make the tickmarks pretty
    
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
 
    #Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    
    
    #Set up the plot limits
    axTemperature.set_xlim(xlims)
    axTemperature.set_ylim(ylims)
            
    
    #Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
    ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)
    
    #Plot the histograms
    axHistx.hist(x, bins=xbins, color = boxcol1, edgecolor=edgecol, density=1, histtype='step', linewidth=2)
    axHistx.hist(xcmp, bins=xbins, color = "grey", edgecolor="black", density=1, linewidth=2,histtype="stepfilled",alpha = 0.7)
    
    #axHistx.hist(x, bins=np.logspace(np.log10(xmin),np.log10(xmax), 50), color = boxcol1, density=0)

    
    axHisty.hist(y, bins=ybins, orientation='horizontal', color = boxcol2, edgecolor=edgecol, density=1, histtype='step', linewidth=2)
    axHisty.hist(ycmp, bins=ybins, orientation='horizontal',  color = "grey", edgecolor="black", density=1,  linewidth=2,histtype="stepfilled",alpha = 0.7)
    #axHisty.hist(y, bins=np.logspace(np.log10(ymin),np.log10(ymax), 50), orientation='horizontal', color = boxcol2, density=0)
    
    #Set up the histogram limits
    axHistx.set_xlim( xmin, xmax )
    axHisty.set_ylim( ymin, ymax )
    
 
    #Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))
    
    axHisty.xaxis.set_major_formatter(FormatStrFormatter('%1.1e'))
    axHistx.yaxis.set_major_formatter(FormatStrFormatter('%1.1e'))
    ltit= "--"
    textstr1 = ltit
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', alpha=0.8)
    # place a text box in upper left in axes coords
    #axHistx.text(9.5, 1.0, textstr1, fontsize=fsize, color='black', fontname='Times New Roman', bbox=props)
    #axHistx.set_ylim(0,1.2)
    #axHisty.set_xlim(0,2.5)
    
    plt.savefig(filename+'.jpeg')
    plt.savefig(filename+'.pdf')
    plt.close()

    return xi,yi,zi,grid_coords.T
    
def simple_map_bella3(x,y,lims,xcmp,ycmp,xlab,ylab, colmap, filename):


    ax,fig=initplot()
    mx=lims[0]#0.9*lims[0]
    Mx=lims[1]#1.2*lims[1]
    my=lims[2]#0.9*lims[2]
    My=lims[3]#1.1*lims[3]
  

    xlabel = xlab
    ylabel = ylab
    cmap   = cm.get_cmap(colmap)
    boxcol1 = 'white'
    boxcol2 = 'white'
    edgecol = '#440154FF'

    # Make the plot
    xmin =  mx
    xmax =  Mx
    ymin =  my
    ymax =  My

    xlims = [xmin,xmax]
    ylims = [ymin,ymax]


    '''
    nbins = 50
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    n = 1000
    zi = zi/zi.sum()
    t = np.linspace(0, zi.max(), n)
    integral = ((zi >= t[:, None, None]) * zi).sum(axis=(1,2))
    t = t[::-1]
    integral = integral[::-1]
    f = interpolate.interp1d(integral, t)
    fname = 'Times New Roman'
    #t_contours = f(np.array([0.99,0.9,0.8,0.7,0.5,0.2]))
    t_contours = f(np.array([0.95,0.5,0.2]))
    fmt = {}
    #strs = ["99%","90%","80%","70%", "50%", "20%"]
    strs = ["95%", "50%", "20%"]
    for l, s in zip(t_contours, strs):
        fmt[l] = s    
    origin = 'lower'
    CS = ax.contour(xi, yi, zi.reshape(xi.shape), t_contours, colors=('black'), linewidths=2, origin=origin)        
    ax.clabel(CS, inline=1, fmt=fmt, fontsize=fsize, fontname=fname)


    nbins = 50
    k = kde.gaussian_kde([xcmp,ycmp])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    n = 1000
    zi = zi/zi.sum()
    t = np.linspace(0, zi.max(), n)
    integral = ((zi >= t[:, None, None]) * zi).sum(axis=(1,2))
    t = t[::-1]
    integral = integral[::-1]
    f = interpolate.interp1d(integral, t)
    fname = 'Times New Roman'
    #t_contours = f(np.array([0.95,0.9,0.8,0.7,0.5,0.2]))
    t_contours = f(np.array([0.95,0.5,0.2]))
    fmt = {}
    #strs = ["99%","90%","80%","70%", "50%", "20%"]
    strs = ["95%", "50%", "20%"]
    for l, s in zip(t_contours, strs):
        fmt[l] = s    
    origin = 'lower'
    CS = ax.contour(xi, yi, zi.reshape(xi.shape), t_contours, colors=('red'), linewidths=2, origin=origin)        
    ax.clabel(CS, inline=1, fmt=fmt, fontsize=fsize, fontname=fname)
    '''
    
    ax.scatter(x,y,alpha=0.1)
    ax.scatter(xcmp,ycmp,alpha=0.1)
    ax.set_xlim(0.1*xmin,1.1*xmax)
    #ax.set_yscale("log")
    plt.savefig(filename+".jpeg")
    plt.close()



    '''
    mx=lims[0]
    Mx=lims[1]
    my=lims[2]
    My=lims[3]
  
    fsize = 43
    lsize = 43
    family = 'Times New Roman'
    matplotlib.rcParams.update({'font.size': fsize})
    matplotlib.rcParams.update({'font.family': family})
    fig   = plt.figure(1, figsize=(18.0,15.0))    

 
    xlabel = xlab
    ylabel = ylab
    cmap   = cm.get_cmap(colmap)
    boxcol1 = 'white'
    boxcol2 = 'white'
    edgecol = '#440154FF'

    # Make the plot
    xmin =  mx
    xmax =  Mx
    ymin =  my
    ymax =  My

    xlims = [xmin,xmax]
    ylims = [ymin,ymax]
    
    nbins = 200
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    grid_coords = np.append(xi.reshape(-1,1),yi.reshape(-1,1),axis=1)
    
    
    # Define the locations for the axes
    left0,right = 0.12,0.55
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left+width+0.02
    
    # Set up the geometry of the three plot
    rect_temperature = [left0, bottom, right, height]   # dimensions of temp plot
    rect_histx       = [left, bottom_h, width, 0.25]    # dimensions of x-histogram
    rect_histy       = [left_h, bottom, 0.25, height]   # dimensions of y-histogram
    

    fig = plt.figure(1, figsize=(15,15))
    

    axTemperature = plt.axes(rect_temperature)
    axHistx       = plt.axes(rect_histx) # x histogram
    axHisty       = plt.axes(rect_histy) # y histogram
    
    #axTemperature.set_xscale('log')
    #axHistx.set_xscale('log')
    axTemperature.set_yscale('log')
    axHisty.set_yscale('log')
    
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    # Make the 'main' temperature plot
    
    # Define the number of bins
    nxbins = 20
    nybins = 20
    nbins = 50
    
    xbins   = linspace(start = xmin, stop = xmax, num = nxbins)
    ybins   = linspace(start = ymin, stop = ymax, num = nybins)
    xcenter = (xbins[0:-1]+xbins[1:])/2.0
    ycenter = (ybins[0:-1]+ybins[1:])/2.0
    
    # Find the min/max of the data
    xmin = min(xlims) #xmin #min(xlims)
    xmax = max(xlims) #xmax #max(xlims)
    ymin = min(ylims) #ymin #min(ylims)
    ymax = max(ylims) #ymax #max(ylims)
    
    H, yedges,xedges = np.histogram2d(y,x,bins=(ybins,xbins))
    
    X = xcenter
    Y = ycenter
    Z = H 
    
    #axTemperature.imshow(H, extent=[xmin,xmax,ymin,ymax],interpolation='nearest', origin='lower',aspect='auto')

    #axTemperature.pcolormesh(xi, yi, zi.reshape(xi.shape),cmap=cmap)
    
    #transform zi to a 0-1 range
    Z = Z = (Z - Z.min())/(Z.max() - Z.min())
    
    n = 1000
    zi = zi/zi.sum()
    #zi = (zi - zi.min())/(zi.max() - zi.min())
    t = np.linspace(0, zi.max(), n)
    integral = ((zi >= t[:, None, None]) * zi).sum(axis=(1,2))


    t = t[::-1]
    integral = integral[::-1]

    f = interpolate.interp1d(integral, t)
    fname = 'Times New Roman'
    #t_contours = f(np.array([0.99,0.9,0.8,0.7,0.5,0.2]))
    t_contours = f(np.array([0.99,0.5,0.2]))
    fmt = {}
    #strs = ["99%","90%","80%","70%", "50%", "20%"]
    strs = ["99%", "50%", "20%"]
    for l, s in zip(t_contours, strs):
        fmt[l] = s
        
    origin = 'lower'
    CS = axTemperature.contour(xi, yi, zi.reshape(xi.shape), t_contours, colors=('white'), linewidths=2, origin=origin)        
    axTemperature.clabel(CS, inline=1, fmt=fmt, fontsize=fsize, fontname=fname)

    
    #Plot the axes labels
    axTemperature.set_xlabel(xlabel,fontsize=lsize,fontname='Times New Roman')    
    axTemperature.set_ylabel(ylabel,fontsize=lsize,fontname='Times New Roman')
    
    #Make the tickmarks pretty
    ticklabels = axTemperature.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    
    ticklabels = axTemperature.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')

    #Make the tickmarks pretty
    
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
 
    #Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(fsize)
        label.set_family('Times New Roman')
    
    
    #Set up the plot limits
    axTemperature.set_xlim(xlims)
    axTemperature.set_ylim(ylims)
            
    
    #Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax-xmin)/nbins)
    ybins = np.arange(ymin, ymax, (ymax-ymin)/nbins)
    
    #Plot the histograms
    axHistx.hist(x, bins=xbins, color = boxcol1, edgecolor=edgecol, density=1, histtype='step', linewidth=2)
    axHistx.hist(xcmp, bins=xbins, color = "grey", edgecolor="black", density=1, linewidth=2,histtype="stepfilled",alpha = 0.7)
    
    #axHistx.hist(x, bins=np.logspace(np.log10(xmin),np.log10(xmax), 50), color = boxcol1, density=0)

    
    axHisty.hist(y, bins=ybins, orientation='horizontal', color = boxcol2, edgecolor=edgecol, density=1, histtype='step', linewidth=2)
    axHisty.hist(ycmp, bins=ybins, orientation='horizontal',  color = "grey", edgecolor="black", density=1,  linewidth=2,histtype="stepfilled",alpha = 0.7)
    #axHisty.hist(y, bins=np.logspace(np.log10(ymin),np.log10(ymax), 50), orientation='horizontal', color = boxcol2, density=0)
    
    #Set up the histogram limits
    axHistx.set_xlim( xmin, xmax )
    axHisty.set_ylim( ymin, ymax )
    
 
    #Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))
    
    axHisty.xaxis.set_major_formatter(FormatStrFormatter('%1.1e'))
    axHistx.yaxis.set_major_formatter(FormatStrFormatter('%1.1e'))
    ltit= "--"
    textstr1 = ltit
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='white', alpha=0.8)
    # place a text box in upper left in axes coords
    #axHistx.text(9.5, 1.0, textstr1, fontsize=fsize, color='black', fontname='Times New Roman', bbox=props)
    #axHistx.set_ylim(0,1.2)
    #axHisty.set_xlim(0,2.5)
    
    plt.savefig(filename+'.jpeg')
    plt.savefig(filename+'.pdf')
    plt.close()

    return xi,yi,zi,grid_coords.T
    '''

    
def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))



def reading(fname, col, skip):
    with open(fname) as f:
        if(skip == "yes"):
            next(f)
            
        c = [float(line.split()[col]) for line in f]
    return c

def readings(fname, col, skip):
    with open(fname) as f:
        if(skip == "yes"):
            next(f)
            
        c = [str(line.split()[col]) for line in f]
    return c


def tMer(a,le,M1,M2):
    e    = le
    c    = 300000.E3;
    G    = 6.67E-11;
    Msun = 1.99E30;
    Rpc  = 1.5e11;
    Gyr  = 365.*24.*3600.;
    tgw  = (5./256.*pow(c,5.)*pow(a*Rpc,4.)*pow(1.-e*e,7./2.)/(G*G*G*M1*M2*(M1+M2)*Msun*Msun*Msun))/Gyr
    return tgw


def initplot():
    # Define the locations for the axes
    left,   right  = 0.14, 0.82
    bottom, height = 0.12, 0.80
    
    # Set up the geometry of the three plot
    rect = [left, bottom, right, height]   # dimensions of temp plot
    fig = plt.figure(1, figsize=(15,15))
    ax      = plt.axes(rect)
    ax.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='on',direction="in", pad=10)
    return ax,fig

def initplot_banner():
    # Define the locations for the axes
    left,   right  = 0.12, 0.82
    bottom, height = 0.12, 0.80
    
    # Set up the geometry of the three plot
    rect = [left, bottom, right, height]   # dimensions of temp plot
    fig = plt.figure(1, figsize=(20,12))
    ax      = plt.axes(rect)
    ax.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='on',direction="in", pad=10)
    return ax,fig




def extract(j,NBrx,NBry):

    NBx = []
    NBy = []
    #NBx.append(0.0)
    #NBy.append(0.0)
    for i in range(len(NBrx)):
        NBx.append(NBrx[i])
        NBy.append(NBry[i])


        
    id1=0
    id2=len(NBy) - 1
    for n in range(len(NBx)):
        if(j > NBy[n]):
            id1 = n
        else:
            id2 = n
            break
        
    if(id2 == 0):
        id1 = 0
        id2 = 1
    if(id1 == len(NBy) - 1):
        id1 = len(NBy) - 2
        id2 = len(NBy) - 1
        
    if(id1 == 0):
        id2 = id1 + 1
    if(id2 == len(NBy) - 1):
        id1 = id2 - 1


    if(NBy[id1] - NBy[id2] == 0.0):
        m1 = max(NBx[id1],NBx[id2])
        m2 = min(NBx[id1],NBx[id2])
        mbhext = m1+(m2-m1)*random.random()
    else:
        mbhext = (NBx[id1] - NBx[id2])/(NBy[id1] - NBy[id2]) * (j - NBy[id2]) + NBx[id2]
   

    #print id1,id2,NBy[id1],NBy[id2],j
        
    return mbhext

def resampling(NBx,NBy,NQx,NQy,Nsource):

    slope = 0.0
    
    maxns = 3.0
    minns = 1.3
    
    maxbh = 60.
    minbh = 5.

    
    newbhm = []
    newqhm = []

    for i in range(Nsource):
        j = random.random()
        if(j<0 or j>1):
            while True:
                j = random.random()
                if(j>0 and j<=1):
                    break
            
        mmm = extract(j,NBx,NBy)
        if(mmm < maxns):
            while True:
                j = random.random()
                if(j<0 or j>1):
                    while True:
                        j = random.random()
                        if(j>0 and j<=1):
                            break
                mmm = extract(j,NBx,NBy)
                if(mmm > minbh and mmm < maxbh):
                    break

        newbhm.append(mmm)

        black = mmm
        
        neutron = -10
        while True:
            j = random.random()
            if(j<0 or j>1):
                while True:
                    j = random.random()
                    if(j>0 and j<=1):
                        break

            mmm = extract(j,NQx,NQy)
            if(mmm > 1.):
                while True:
                    mmm = extract(j,NQx,NQy)
                    if(mmm < 1.):
                        break

            neutron = black * mmm

            if(neutron > minns and neutron < maxns):
                break
        
        newqhm.append(mmm)

                

    
    return newbhm,newqhm
        

def newextract(mbhttZh,coeff, mmin, mmax, mmint):

    num   = 10000
    kdem  = stats.gaussian_kde(mbhttZh,coeff)

    x_eval = np.linspace(mmin,mmax,num=num)
    y_eval = kdem(x_eval)    
    c_eval = []
        
    cnt = 0
    for i in range(len(x_eval)):
        cnt += y_eval[i]
        c_eval.append(cnt)

    for i in range(len(x_eval)):
        c_eval[i] = c_eval[i] / cnt
        y_eval[i] = y_eval[i] 

    #print cnt,num,cnt*(mmax-mmin)/num,mmax-mmin
        
    xbhbin = x_eval
    ybhbin = y_eval
    cmbin  = c_eval
    
    return xbhbin,ybhbin,cmbin



def newsample(mbhttZh,mnsttZh,qraZh):

    coeff = 0.02
    mmin  = min(mbhttZh)
    mmax  = max(mbhttZh)

    x_eval,y_eval,c_eval = newextract(mbhttZh,coeff,mmin,mmax,0.0)

    xbhbin = x_eval
    ybhbin = y_eval
    cmbin  = c_eval
    

    coeff = 0.05
    mmin = min(mnsttZh)
    mmax = max(mnsttZh)
    x_eval,y_eval,c_eval = newextract(mnsttZh,coeff,mmin,mmax,0.0)
    xnsbin = x_eval
    ynsbin = y_eval
    cnsbin = c_eval

        
    coeff = 0.08
    mmin = min(qraZh)
    mmax = max(qraZh)
    x_eval,y_eval,c_eval = newextract(qraZh,coeff,mmin,mmax,0.0)
    xqbin = x_eval
    yqbin = y_eval
    cqbin = c_eval
   
        
    #    pltin(mbhttZh,mnsttZh,qraZh, xbhbin, xnsbin, xqbin,ybhbin, ynsbin, yqbin,cmbin,  cnsbin, cqbin)
    
    return  ybhbin,xbhbin,cmbin,ynsbin,xnsbin,cnsbin,yqbin,xqbin,cqbin



def pltin(mbh,mns,q, xbhbin, xnsbin, xqbin,ybhbin, ynsbin, yqbin,cmbin,  cnsbin, cqbin,ii,name):
    ### Plotting ###
    ax,fig = initplot()
    fvec  = [mbh, mns, q]
    fxvec = [xbhbin, xnsbin, xqbin]
    fyvec = [ybhbin, ynsbin, yqbin]
    fzvec = [cmbin,  cnsbin, cqbin]
    mmmax = [50., 4. , 1.]

    f  = fvec[ii]
    fx = fxvec[ii]
    fy = fyvec[ii]
    mmm = mmmax[ii] 

    wgh = [ 1./len(f) for i in range(len(f))] 

    ax.hist(f, 50, range=[0,mmm], alpha=0.2, edgecolor="black", facecolor="red", cumulative=False, density=True)
    yedge,xedge = np.histogram(f,50,range=[0,mmm],  density=True)
    cnt = 0
    for i in range(len(yedge)):
        cnt  +=  yedge[i]

    x = [0.5*(xedge[i]+xedge[i+1]) for i in range(len(xedge)-1)]
    y = [yedge[i] *1./ cnt for i in range(len(yedge))]

    ax.plot(x,y,alpha=1,lw=3,linestyle = "--",color="black")
    ax.plot(fx,fy,alpha=1.0, lw = 2,color="blue")
    plt.savefig("over_check"+name+".jpeg")
    plt.close()
    return
    
def interp(x,Nx,Ny):

    id1 = 0
    id2 = len(Nx)-1
    if(Nx[0] < Nx[1]):
        for i in range(len(Nx)-1):
            if(x>Nx[i]):
                id1 = i
            elif(x<=Nx[i]):
                id2 = i
                break
    else:
        for i in range(len(Nx)-1):
            if(x<Nx[i]):
                id1 = i
            elif(x>=Nx[i]):
                id2 = i
                break


    if(id1 >= len(Nx)-1):
        #print "Warning, we arrived at the end ", x, Nx[id1], Nx[id2], id1, id2
        id2 = id1
        id1 = id2-1
        
    if(id2 == 0):
        #print "Warning, we never leave the homerun", x, Nx[id1], Nx[id2], id1, id2
        id1 = id2
        id2 = id1+1
     
    y = (Ny[id2]-Ny[id1])/(Nx[id2] - Nx[id1])*(x-Nx[id2]) + Ny[id2]
    return y


def check(mbh, mbh_end, Mn, Mx, ftyp, cumul):
    dpur = '#440154' # dark purple
    purp = '#472c7a' # purple
    blu2 = '#3b518b' # blue
    blu3 = '#2c718e' # blue
    blgr = '#21908d' # blue-green
    gre1 = '#27ad81' # green
    gre2 = '#5cc863' # green
    lgre = '#aadc32' # lime green
    yell = '#fde725' # yellow
    ax,fig=initplot()
    wg = [1./len(mbh) for i in range(len(mbh))]
    ax.hist(mbh,    50,range=[Mn,Mx],histtype="step",lw=3,weights=wg,label="original",cumulative=cumul)
    wg = [1./len(mbh_end) for i in range(len(mbh_end))]
    ax.hist(mbh_end,50,range=[Mn,Mx],histtype="step",lw=3,weights=wg,label="resample",cumulative=cumul)

    if(cumul == True):
        legloc = "upper left"
        
    if(cumul == False):
        legloc = "upper right"
        ax.set_ylim(0,0.4)

    plt.legend(loc=legloc)
    plt.savefig("over_"+ftyp+".jpeg")
    plt.close()


def check2(mbh, mbh_res, mbh_end, Mn, Mx, ftyp, cumul):
    dpur = '#440154' # dark purple
    purp = '#472c7a' # purple
    blu2 = '#3b518b' # blue
    blu3 = '#2c718e' # blue
    blgr = '#21908d' # blue-green
    gre1 = '#27ad81' # green
    gre2 = '#5cc863' # green
    lgre = '#aadc32' # lime green
    yell = '#fde725' # yellow

    ax,fig=initplot()
    wg = [1./len(mbh) for i in range(len(mbh))]
    ax.hist(mbh,    20,range=[Mn,Mx], alpha = 0.8, facecolor = dpur, edgecolor="black", histtype="stepfilled",lw=3,weights=wg,label="original",cumulative=cumul)
    wg = [1./len(mbh_res) for i in range(len(mbh_end))]
    ax.hist(mbh_res,20,range=[Mn,Mx], alpha = 0.6, facecolor = blu2, edgecolor="black", histtype="stepfilled",lw=3,weights=wg,label="resample",cumulative=cumul)
    wg = [1./len(mbh_end) for i in range(len(mbh_end))]
    ax.hist(mbh_end,20,range=[Mn,Mx], alpha = 0.4, facecolor = lgre, edgecolor="black", histtype="stepfilled",lw=3,weights=wg,label="accereje",cumulative=cumul)

    legloc = "upper left"
    if(cumul == False):
        legloc = "upper right"
    
    ax.set_ylim(0,0.4)
    plt.legend(loc=legloc)
    plt.savefig("over_"+ftyp+".jpeg")
    plt.close()



def massspec(NSH,BHH,NSL,BHL):
    
    # Define the locations for the axes
    left,   right  = 0.14, 0.86
    bottom, height = 0.12, 0.80
    
    # Set up the geometry of the three plot
    rect = [left, bottom, right, height]   # dimensions of temp plot
    fig = plt.figure(1, figsize=(18,15))
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    # Turn off axis lines and ticks of the big subplot
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(top='off', bottom='off', left='off', right='off')
    
    # plot the same data on both axes
    
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    labels = []
    
    ax.set_ylim( 0.,0.5)
    ax1.set_ylim(0.,0.5)
    ax2.set_ylim(0.,0.5)
    ax1.set_xlim(1.0,3.0)
    ax2.set_xlim(4.0,60)
    
    ax2.axvline(23.2,ymin=0.0,ymax=1.,color="red",lw=3,linestyle="--",alpha=1)
    ax2.axvline(24.3,ymin=0.0,ymax=1.,color="red",lw=3,linestyle="-",alpha=1)
    ax2.axvline(22.2,ymin=0.0,ymax=1.,color="red",lw=3,linestyle="-",alpha=1)
    
    ax1.axvline(2.59,ymin=0.0,ymax=1.,color="blue",lw=3,linestyle="--",alpha=1)
    ax1.axvline(2.67,ymin=0.0,ymax=1.,color="blue",lw=3,linestyle="-",alpha=1)
    ax1.axvline(2.50,ymin=0.0,ymax=1.,color="blue",lw=3,linestyle="-",alpha=1)
    

    WH = [1./len(NSH) for i in range(len(NSH))]
    WL = [1./len(NSL) for i in range(len(NSL))]

    
    ax1.hist(NSH, 10, range=[0.8,3], histtype="stepfilled", weights=WH, facecolor="blue" , alpha = 0.4, rasterized=True, label="$Z = 0.02$")
    ax2.hist(BHH, 10, range=[4,30],  histtype="stepfilled", weights=WH, facecolor="red"  , alpha = 0.4, rasterized=True, label="$Z = 0.02$")
    ax1.hist(NSL, 10, range=[0.8,3], histtype="step"      , weights=WL, edgecolor="blue" , alpha = 1.0, lw = 5, linestyle="--",rasterized=True, label="$Z = 0.0002$")
    ax2.hist(BHL, 20, range=[4,60],  histtype="step"      , weights=WL, edgecolor="red"  , alpha = 1.0, lw = 5, linestyle="--",rasterized=True, label="$Z = 0.0002$")
    


    ax.set_ylabel("$\mathrm{d}N/\mathrm{d}M$",labelpad = 80)
    ax1.set_xlabel("$M_\mathrm{NS}$ [M$_\odot$]",labelpad=10)
    ax2.set_xlabel("$M_\mathrm{BH}$ [M$_\odot$]",labelpad=10)
    
    
    # hide the spines between ax and ax2
    ax1.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax1.yaxis.tick_left()
    ax1.tick_params(labelright='off')
    ax2.tick_params(labelleft='off')

    ax1.tick_params(axis="both",which="both",top='on', bottom='on', left='on', right='off',direction="in", pad=10)
    ax2.tick_params(axis="both",which="both",top='on', bottom='on', left='off', right='on',direction="in", pad=10)

    d = .015 # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((1-d,1+d), (-d,+d), **kwargs)
    ax1.plot((1-d,1+d),(1-d,1+d), **kwargs)
    
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d,+d), (1-d,1+d), **kwargs)
    ax2.plot((-d,+d), (-d,+d), **kwargs)
    
    
    
    #ax1.hist(allns ,10, range=[0.8,3], histtype="step",lw=2,color="black",linestyle="-.",weights=weigh_allns,label="Z = 0.0002\n NSMF")
    #ax1.hist(allnsS,10, range=[0.8,3], histtype="step",lw=2,color="black",linestyle=":",weights=weigh_allnsS,label="Z = 0.02\n NSMF")
    
    fig.subplots_adjust(bottom=0.11, top=0.89, left=0.15, right=0.90,
                        wspace=0.1, hspace=0.01)
    
    leg1 = ax1.legend(loc="upper left",title="NS",framealpha=1)
    leg2 = ax2.legend(loc="upper left",title="BH",framealpha=1)
    #ax.add_artist(leg1)
    #ax.add_artist(leg2)
    
    
    fig = plt.gcf()
    fig.suptitle("Volume Weighted",y=0.98)
    
    plt.savefig("MassSpectr.jpeg")
    plt.savefig("MassSpectr.pdf")
    
    plt.close()



def inter(x, X, Y, N):    
    id1=0
    id2=N-1
    if(X[0] > X[1]):
        for i in range(N):
            if(x < X[i]):
                id1 = i
            if(x >= X[i]):
                id2 = i
                break
    elif(X[0] <= X[1]):
        for i in range(N):
            if(x >= X[i]):
                id1 = i
            if(x < X[i]):
                id2 = i
                break

    if(id2 == 0):
        id1 = 0
        id2 = 1
        
    if(id1 == N-1):
        id1 = N-2
        id2 = N-1

    
    
    y = (Y[id1]-Y[id2])/(X[id1]-X[id2])*(x-X[id1]) + Y[id1]
  
    
    return y



def package(testo, simID):

    

    redfile = "/home/manuel/Scrivania/ACTIVE_PROJECTS/redshift_time.txt"
    reds = reading(redfile, 0, "no")
    tage = reading(redfile, 1, "no")

    
    m1 = reading(testo, 5, "no")
    m2 = reading(testo, 6, "no")
    tf = reading(testo, 13,"no")
    tg = reading(testo, 14,"no")
    st = readings(testo, 18, "no")
    zgw = reading(testo, 22, "no")



    mbin= np.array(m1)+np.array(m2)
    qbin= [min(m1[i],m2[i])/max(m1[i],m2[i]) for i in range(len(m1))]
    
    zre= [inter(tg[i]/1.E9, tage, reds, len(reds)) for i in range(len(tg))]    
    
    ZRE = np.logspace(-2.0,1.0,50)
    MRE = [[] for i in range(len(ZRE)-1)]
    NME = [0.0 for i in range(len(ZRE)-1)]
    for i in range(len(zre)):
        if(m1[i]+m2[i] > 500.):
            continue        
        for k in range(len(ZRE)-1):
            if(zre[i] >= ZRE[k] and zre[i] < ZRE[k+1]):
                MRE[k].append(m1[i]+m2[i])
                NME[k] += 1.0
                break

    NTOT = sum(NME)
    print ("Nmergers = ",NTOT)
    if(NTOT == 0.0):
        NTOT = 1
    
    for i in range(len(NME)):
        NME[i] /= NTOT
        
    redi = [0.5*(ZRE[i]+ZRE[i+1]) for i in range(len(ZRE)-1)]   
    medi = [np.percentile(MRE[k],50) if len(MRE[k]) > 0 else 0.0 for k in range(len(ZRE)-1)]
    seve = [np.percentile(MRE[k],70) if len(MRE[k]) > 0 else 0.0 for k in range(len(ZRE)-1)]
    thir = [np.percentile(MRE[k],30) if len(MRE[k]) > 0 else 0.0 for k in range(len(ZRE)-1)]
    nine = [np.percentile(MRE[k],90) if len(MRE[k]) > 0 else 0.0 for k in range(len(ZRE)-1)]
    teni = [np.percentile(MRE[k],10) if len(MRE[k]) > 0 else 0.0 for k in range(len(ZRE)-1)]
    
    ax,fig = initplot()
    ax.scatter(zre, mbin, s=200, alpha=0.3, edgecolor="black", color="red")
    ax.plot(redi, medi, lw = 5, linestyle="-", c="black")
    ax.plot(redi, seve, lw = 5, linestyle="--", c="black")
    ax.plot(redi, thir, lw = 5, linestyle="--", c="black")
    ax.plot(redi, nine, lw = 5, linestyle=":", c="black")
    ax.plot(redi, teni, lw = 5, linestyle=":", c="black")

    ax.set_ylim(2,500)
    ax.set_xlim(1.E-3,10.)
    ax.set_xlabel("z")
    ax.set_ylabel("$M_\mathrm{BBH}$ [M$_\odot$]")
    ax.set_yscale("log")
    ax.set_xscale("log")
    tit = "ID"+simID #"$M_\mathrm{cl} = (0.5-5)\\times 10^6$ M$_\odot$ \n $R_h = (1-2)$ pc \n $Z = 0.01$ Z$_\odot$"
    plt.legend(loc="upper left", title=tit)
    plt.savefig("redshift_evol_mass_test"+simID+".jpeg")
    plt.close()



    ax,fig = initplot()
    ax.scatter(redi, NME, s=200, alpha=0.8, edgecolor="black", color="red")
    ax.set_ylim(1.E-3,0.1)
    ax.set_xlabel("z")
    ax.set_ylabel("$n_\mathrm{GW}$ ")
    ax.set_xscale("log")
    ax.set_yscale("log")
    tit = "ID"+simID #"$M_\mathrm{cl} = (0.5-5)\\times 10^6$ M$_\odot$ \n $R_h = (1-2)$ pc \n $Z = 0.01$ Z$_\odot$"
    plt.legend(loc="lower left", title=tit)
    plt.savefig("redshift_evol_test"+simID+".jpeg")
    plt.close()


    ZRE = np.linspace(min(zre),max(zre),50)
    MRE = [[] for i in range(len(ZRE)-1)]
    for i in range(len(zre)):
        if(m1[i]+m2[i] > 500.):
            continue
        for k in range(len(ZRE)-1):
            if(zre[i] >= ZRE[k] and zre[i] < ZRE[k+1]):
                MRE[k].append(qbin[i])
                break

    redi = [0.5*(ZRE[i]+ZRE[i+1]) for i in range(len(ZRE)-1)]   
    medi = [np.percentile(MRE[k],50) if len(MRE[k]) > 0 else 0.0 for k in range(len(ZRE)-1)]

    ax,fig = initplot()
    ax.scatter(zre, qbin, s=200, alpha=0.3, edgecolor="black", color="red")
    ax.plot(redi, medi, lw = 5, linestyle="-", c="black")
    ax.set_ylim(0,1)
    ax.set_xlabel("z")
    ax.set_ylabel("$q_\mathrm{BBH}$")
    ax.set_xscale("log")
    tit = "ID"+simID #"$M_\mathrm{cl} = (0.5-5)\\times 10^6$ M$_\odot$ \n $R_h = (1-2)$ pc \n $Z = 0.01$ Z$_\odot$"
    plt.legend(loc="upper left", title=tit)
    plt.savefig("redshift_evol_ratio_test"+simID+".jpeg")
    plt.close()



    
