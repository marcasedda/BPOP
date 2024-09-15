import matplotlib
import matplotlib.colors as colors
import matplotlib.pyplot as plt

import numpy as np
import seaborn as sns
import pandas as pd
from sklearn.neighbors import KernelDensity
import scipy
from scipy.interpolate import splev, splrep, BSpline, interp1d
from scipy import signal
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
from random import random
import os
from os import path

def init_dir(fdyn, fgc, fyc, fnc, spinI, spinD, clurad, bhseed, umg, highg, mix, vms):
    f = "/SIM_Fdyn"+str(fdyn)+"_Ngc"+str(fgc)+"_Nyc"+str(fyc)+"_Nnc"+str(fnc)+"isolS_"+spinI+"dynaS_"+spinD+"_MetalDivi_0_8_logflat_logflat_Correction_no_SFR___sfronly_no"+vms+"CluRh_"+clurad+"_0.3BHseed"+bhseed+"_UMG"+umg+"_HighG"+highg+"_mix_"+mix+"/"
    return f


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


def psi(typ, z, CFE, B, Zn, Sn):

    if(typ == "mf17"):
        sfr = CFE * 1.0E-2 * pow(1.+z,2.6) / (1. + pow((1.+z)/3.2,6.2))        
    elif(typ == "kr13"):
        if(2 <= z <= 8):
            sfr = CFE * 0.003            
        else:
            sfr = 0.0
    else:
        sfr =  psiNC(z, B, Zn, Sn)

    return sfr

def psiNC(z, B, Zn, Sn):
    p = B * np.exp(-pow(z-Zn,2.) / (2.*Sn*Sn))
    return p

def GSS_cdf_f(x, xmean, xdisp):
  return 0.5*(1.+math.erf((x - xmean)/(np.sqrt(2.)*xdisp)))

def tgw(m1,m2,a,e):

    cl = 3.E8
    G  = 6.67E-11
    Ms = 1.99E30
    Rs = 1.5E11
    
    tgw = 256./5. * pow(cl,5.) * pow(a * Rs,4.) * pow(G*Ms,-3.) /(m1*m2*(m1+m2)) * pow(1.-e*e,3.5)
    tgw = tgw /(365.*24.*3600.)
   
    return tgw


def DtDz(zz, H0now, Omatt, Olamb):
    dtdz_now = (1.+zz) * H0now * ( 1.E3 / 3.08E22 * 365.*24.*3600. ) * np.sqrt(pow(1.+zz,3.)*Omatt + Olamb)
    return 1./dtdz_now



def test(gpc_to_mpc3, H0, Om, Ol):
    #TESTS#
    z = np.linspace(0,10,100)
    z = np.array(z)
    I = [None] * len(z)
    dtdz = DtDz(z, H0, Om, Ol)
    I = dtdz
    Thubb = scipy.integrate.simpson(I, z)
    print("Age of the Universe=",Thubb*1.E-9," Gyr")
    
    I = [None] * len(z)
    #### Local Volume mass density ####
    typ = ["mf17", "kr13", "eb19"]
    mty = [3.E5, 3.E5, 3.E5]
    cty = [3.e5 * 200./6.e10, 1.0, 1.0]
    for i in range(len(typ)):
        
        dtdz = DtDz(z, H0, Om, Ol)
        I = [psi(typ[i], z[j], cty[i], 2.E-4, 3.2, 1.5)*dtdz[j] for j in range(len(z))]
        
        R = scipy.integrate.simpson(I, z)

    print("GCSF-{0:s}".format(typ[i]), "Number density :", "{0:1.2f}".format(R/mty[i]), " Mpc^{-3}")

    return 



def comdis(z, H0, OmegaM, OmegaL, c):
    A = c/H0 * 1. / np.sqrt(pow(1.+z,3.)*OmegaM + OmegaL) 
    return A


def volume(Da, z, H0, OmegaM, OmegaL, c):  
    dVdz = (4.*np.pi) * (c/H0)  * (Da*Da) / (np.sqrt(pow(1.+z,3.)*OmegaM + OmegaL))
    return dVdz

def GlobRate(simname, maxZ, clt, lab, Hub0, OmegaM, OmegaL):
    H0 = Hub0 * 1000.
    c  = 3.E8
    
    ##### Calculate the cosmological volume ######
    nbin = 1000
    zmin = 0.0
    zmax = 15.0
    dz = (zmax - zmin)/nbin
    Rat=[0.0] * len(clt)
    
    for kt in range(len(clt)):
        
        X = np.loadtxt("Merger_rates_"+simname+"_"+clt[kt]+".txt", usecols=[0,1], unpack=True)
        Xz = X[0]
        Xr = X[1]
    

        rex = [0.0]*nbin
        Ivo = [0.0]*nbin

        for i in range(nbin):
            
            z = zmin + i*dz
            
            Rate = interp(z, Xz, Xr)
            
            if(z == 0):
                rex[i] = 0.0
                Ivo[i] = 0.0
                continue
    

        
            Ivol = [0.0]*nbin
            z_int= [0.0]*nbin
        
            #Calculating the comoving distance within a given z//
            for ii in range(nbin):
                z_int[ii] = z*ii*1./nbin
                Ivol[ii] = comdis(z_int[ii], H0, OmegaM, OmegaL, c)
            
            IV = scipy.integrate.simpson(Ivol, z_int)
            DL = (1.+z)*IV*1.E6;

            #//Cosmological volume//
            #//Checked, this works -- this is in Mpc^3!!
            cosmovol = volume(IV, z, H0, OmegaM, OmegaL, c)
            
            Ivo[i] = Rate * (cosmovol/1.e9) / (1.+z)
            rex[i] = z
            
    
        Vol = scipy.integrate.simpson(Ivo, rex)
        Rat[kt] = Vol
        #print(clt[kt], Rat[kt])

    Total = np.sum(Rat)
    Dynamical = Rat[0]+Rat[1]+Rat[2]
    print("f_dyn = ", "{0:1.4f}".format(Dynamical*1.0/Total))

   
    
    for kt in range(len(clt[kt])-1,-1,-1):
        if(clt[kt] != "none"):
            lst = "f_"+lab[kt]
            print(lst, "= ", "{0:1.4f}".format(Rat[kt]/Dynamical))
    
    
    #cFR. NED WRIGHTS - (z,V) = 0.5 - 30.848 Gpc^3)!!!! Works!
    
    return
