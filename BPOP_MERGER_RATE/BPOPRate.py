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
import time

from lib_BPOPRate import *

start = time.time()

gpc_to_mpc3 = 1.e9
H0 = 67.74
Om = 0.3089
Ol = 0.6911

plt.rc('text', usetex=True)
cmap = "viridis"
fsize = 41
lsize = 40
family = 'Times New Roman'
matplotlib.rcParams.update({'font.size': fsize})
matplotlib.rcParams.update({'font.family': family})


print("\n===================  BPOPRate package  =====================")


mstar   = 1.0
fbh     = 0.0026
fbin_is = 1.0
fbin_yc = 1.0
fbin_gc = 1.0
fbin_nc = 1.0

    
## User-defined quantities: type of cluster, CSFH, binary fraction, typical mass, max redshift, and repository path and name ##

IDmod = 1
a = "b"
maxz = 15.0
massive = "no"
max_mas = 1.E5

fdyn  = [0.85,0.85,0.85]
fgc   = [0.3 ,0.3 ,0.3]
fyc   = [0.4 ,0.4 ,0.4]
fnc   = [0.3 ,0.3 ,0.3]
spinI = ["maxwellian02","maxwellian02","maxwellian02"]
spinD = ["maxwellian02","maxwellian02","maxwellian02"]
clurad= ["AS20", "AS20", "AS20"]
bhseed= ["no"  , "no"  , "density"]
umg   = ["no"  , "yes"  , "yes"]
highg = ["no"  , "no"  , "no"]
mix   = ["0.5G", "0.5G", "0.5G"]
vms   = [""    , ""    , "VMSrho0.2"]

#-------------------------------------------------------------
clt  = ["nuclear", "globular", "young", "none"]
lab  = ["nuc", "glo", "you", "iso"]
if(a == "a"):
    sfrt = ["mf17", "eb19", "mf17", "mf17"]
elif(a == "b"):
    sfrt = ["eb19", "eb19", "mf17", "mf17"]
elif(a == "z"):
    sfrt = ["kr13", "kr13", "mf17", "mf17"]

if(clurad[IDmod] == "Mapelli20"):
    mtyp = [pow(10.,6.18), pow(10.,5.6), pow(10., 4.3), 1.0]
else:
    mtyp = [pow(10.,6.36), pow(10.,5.01), pow(10., 3.53), 1.0]
    
fbin = [fbin_nc, fbin_gc, fbin_yc, fbin_is]

#-------------------------------------------------------------


simname = "SET"+str(IDmod) 
home = "../PRODUCTION_RUNS_SEP2024/"+simname+"/"
codename = [simname+"_Continuous"]


fl = init_dir(fdyn[IDmod], fgc[IDmod], fyc[IDmod], fnc[IDmod], spinI[IDmod], spinD[IDmod], clurad[IDmod], bhseed[IDmod], umg[IDmod], highg[IDmod], mix[IDmod],vms[IDmod])

print("\n File analysed: ",fl)
print(simname+a, "\n")

print("CLType: ", clt)
print("CSFH:   ", sfrt)
print("fbin:   ", fbin)
print("m_typ:  ", np.log10(mtyp))
print("z_max:  ", maxz)
print("------------------------------------------------------------\n")

#### Calculations start here ####

Z = [0.0002,
     0.0004,
     0.0008,
     0.0012,
     0.0016,
     0.002,
     0.004,
     0.006,
     0.008,
     0.012,
     0.016,
     0.02]

eta_eff  = [[1.0]*len(Z), [1.0]*len(Z), [1.0]*len(Z), [1.0]*len(Z)]


#Merger efficiency obtained directly from the catalogues! 
X = np.loadtxt("eta_isolated_binaries",usecols=[3,2,1],unpack=True)
zmer_mobse = X[0]
pmer_mobse = X[1]
nmer_mobse = X[2]
if(zmer_mobse[0] != Z[0]):
    print("Error eta_isolated_binaries")
    exit()

# Directory and files to be analysed #

directory = [fl]
             
directory = [home + codename[i] + directory[i] for i in range(len(directory))]

file_all = "Larger_than_tH.txt"
file_mer = "Catalogue.txt"
file_esc = "retention_Z"


fmix = 0.5
        

z = np.linspace(0,15.0,200)
z = np.array(z)

Bf   = [[1.0]*len(z), [1.0]*len(z), [1.0]*len(z), [1.0]*len(z)]

def msca(CLU, SFR):
    fsca = [1.0,1.0,1.0,1.0]
    
    if((CLU == "nuclear" or CLU == "globular") and (SFR != "kr13")):
        if(SFR == "mf17"):
            if(CLU == "nuclear"):
                fsca[0] = 3.E7/6.E10
            elif(CLU == "globular"):
                fsca[0] = 6.E-5 * 1.E12 / 6.E10                
        elif(SFR == "eb19"):
            if(CLU == "nuclear"):
                fsca[1] = 1.E-5
                fsca[2] = 3.2
                fsca[3] = 1.5
            elif(CLU == "globular"):
                fsca[1] = 2.E-4
                fsca[2] = 3.2
                fsca[3] = 1.5

    if(SFR == "kr13"):
        if(CLU=="globular"):
            fsca[0] = 0.1
        if(CLU=="nuclear"):
            fsca[0] = 0.1 * 3.E7 / (200. * 3.E5) 
            
    return fsca
    

        

## Preliminarly, we need to evaluate the scaling between the SFR among different channels ##


sfr   = [None]*len(z)
sfr_N = [None]*len(z)
sfr_G = [None]*len(z)
sfr_Y = [None]*len(z)
sfr_I = [None]*len(z)
fysc_mw = 0.02 #(e.g. Larsen09, Piskunov+08)

for i in range(len(z)):

    sfr_tot =  psi("mf17", z[i], 1.0, 1.0, 1.0, 1.0)
    sfr[i] = sfr_tot
    
    ## Evaluate the scale factor w.r.t. the SFR ##
    fsca   = msca(clt[0], sfrt[0])
    sfr_nc = psi(sfrt[0], z[i], fsca[0], fsca[1], fsca[2], fsca[3])
    sfr_N[i] = sfr_nc    
    
    fsca   = msca(clt[1], sfrt[1])
    sfr_gc = psi(sfrt[1], z[i], fsca[0], fsca[1], fsca[2], fsca[3])
    sfr_G[i] = sfr_gc
    
    B_yc = max(0, min(fysc_mw, 1.-(sfr_nc+sfr_gc)/sfr_tot))
    sfr_yc = B_yc * psi(sfrt[2], z[i], 1.0, 1.0, 1.0, 1.0)
    sfr_Y[i] = sfr_yc    
    
    B_is = max (0, 1.0 - (sfr_yc+sfr_gc+sfr_nc)/sfr_tot)
    sfr_is = B_is * psi(sfrt[3], z[i], 1.0, 1.0, 1.0, 1.0) 
    sfr_I[i] = sfr_is
    
    Bf[2][i] = B_yc
    Bf[3][i] = B_is

    
ax, fig = initplot()
ax.plot(z, sfr, lw = 5, label="tot")
ax.plot(z, sfr_I, lw = 5,label="IS")
ax.plot(z, sfr_Y, lw = 5,label="YC")
ax.plot(z, sfr_G, lw = 5,label="GC")
ax.plot(z, sfr_N, lw = 5,label="NC")
ax.set_xlabel("$z$")
ax.set_ylabel("CSFR")
ax.set_yscale("log")
ax.set_ylim(1.E-7,0.5)
ax.set_xlim(0.0,15.0)
plt.legend(loc="upper right")
plt.savefig("SFR_uni.jpeg")
plt.close()
#exit()

print("--------------- Local denstity in Mpc^{-3} -----------------")

for kt in range(len(clt)):
    
    cluster = clt[kt]
    typ = sfrt[kt]
    fsca   = msca(cluster, typ)

    I = [None]*len(z)
    dtdz = DtDz(z, H0, Om, Ol)
    
    for j in range(len(z)):
        sfr = Bf[kt][j] * psi(typ, z[j], fsca[0], fsca[1], fsca[2], fsca[3])
        I[j] = sfr * dtdz[j]
        
    Cv = scipy.integrate.simpson(I,z)

    prtstr = "rho_"+lab[kt] + " (SFR " + sfrt[kt] + ") = "
    prtstr += "{0:1.5f}".format(Cv / mtyp[kt])


    print (prtstr)

print("------------------------------------------------------------")

## Now we are ready to calculate the merger rate, taking into account that the above loop calculate the star formation rate of different environments ##
rt = np.loadtxt("./redshift_time.txt", usecols=[0,2], unpack=True)
redshift = rt[0]
lookback = rt[1]

## Reading the file ## 

Fn = directory[0] + file_mer
X = np.loadtxt(Fn, dtype={'names':('zmer', 'zfor', 'met', 'clus', 'mclu', 'vesc','mas', 'tmer', 'nrep', 'msec','rclu'), 'formats':(float, float, float, 'U8',float,float,float,float,float,float,float)}, usecols=[24, 25, 1, 3, 15, 17, 5, 14, 2, 6,16], unpack=True)

Ln = directory[0] + file_all
Y = np.loadtxt(Ln, dtype={'names':('met','clus','nrecy'), 'formats':(float,'U8',float)}, usecols=[13,14,22], unpack=True)                

zmer_a = X[0]
zfor_a = X[1]
Zmet_a = X[2]
ctyp_a = X[3]
mclu_a = X[4]
vesc_a = X[5]
mpri_a = X[6]
msec_a = X[9]
tmer_a = X[7]
nrep_a = X[8]
rclu_a = X[10]

mthi_a = [max(mpri_a[i], msec_a[i]) for i in range(len(mpri_a))]

mpri_a = np.array(mthi_a)


if(massive == "no"):
    idx = np.where((zmer_a > 0.) & (mpri_a < max_mas))
else:
    idx = np.where((zmer_a > 0.) & (mpri_a > max_mas))
    
tmer = tmer_a[idx]
zmer = zmer_a[idx]
zfor = zfor_a[idx]
Zmet = Zmet_a[idx]
ctype= ctyp_a[idx]
nrep = nrep_a[idx]

Mmer = mclu_a[idx]
Rmer = rclu_a[idx]
Vmer = [0.0]*len(Rmer)
for i in range(len(Mmer)):
    if(Mmer[i] == 0.0):
        continue
    
    Vmer[i] = np.sqrt(Mmer[i] / Rmer[i] * 6.67E-11*1.99E30/3.08E16) * 0.001
    
Vmer = np.array(Vmer)

Zloss= Y[0]
Closs= Y[1]
nloss= Y[2]
        
f = open("efficiency_Z_SET"+str(IDmod)+a+".txt","w")
fstr = ""
for kt in range(len(clt)):         
            
    cluster = clt[kt]
    typ = sfrt[kt]

    ## Loss fraction -- to calculate the eta(Z) parameter ##

    dZ = 0.0001

    for i in range(len(Z)):
        dcen = Z[i]

        if(cluster == "none"):
            eta = 0.5 * 0.285 * pmer_mobse[i]  

        else:
            dfile_esc = directory[0]+file_esc+str(Z[i])+".dat"
            if(os.path.isfile(dfile_esc) == False):
                eta_eff[kt][i] = 0.0
                continue
            
            Xesc = np.loadtxt(dfile_esc, usecols=[0,1,2],unpack=True)
        
        
            if(i < len(Z)-1):
                dup = Z[i+1]
            else:
                dup = Z[i] * 1.1
        
            if(i > 0):
                dlow = Z[i-1]
            else:
                dlow = 0.9 * Z[i]
            
            dZsx = (dcen-dlow)/2.0
            dZdx = (dup-dcen)/2.0
        
            ## How many mergers per metallicity bin ? ##
            ## This is the denominator in the calculation of eta(Z) ##

            idx_Z = np.where((Zmet >= dcen - dZsx) & (Zmet < dcen + dZdx) & (ctype == cluster))
            idx_N = np.where((Zloss >= dcen - dZsx) & (Zloss < dcen + dZdx) & (Closs == cluster) & (nloss==0))
            Zmet_Z = Zmet_a[idx_Z]
            Zmet_no= Zloss[idx_N]

            ## Multiple mergers should also be counted here? ##
            nrep_Z = nrep[idx_Z]
            nrep_sum = np.sum(nrep_Z)

            Mmer_mer = Mmer[idx_Z]
            Vmer_mer = Vmer[idx_Z]
            fret_nat = np.ones(len(Vmer_mer))
            fret_mix = np.ones(len(Vmer_mer))
            for j in range(len(Vmer_mer)):
                if(str(Vmer_mer[j]) == "nan" or str(Vmer_mer[j]) == "-nan"):
                    print("Error")
                    exit()
                fret_nat[j] = interp(Vmer_mer[j],Xesc[0],Xesc[1]) 
                fret_mix[j] = interp(Vmer_mer[j],Xesc[0],Xesc[2])            

            fret = (1.-fmix) * fret_nat + fmix * fret_mix
            Mweight = np.sum(fret * Mmer_mer)
            Msummed = np.sum(Mmer_mer)
            
            NZgw = len(Zmet_Z) 
            NZno = len(Zmet_no)
            eta = 0.0
            if(mstar != 0.0 and (NZgw+NZno != 0) and Msummed!=0):
                eta = fbh / mstar * (NZgw + nrep_sum) / (NZgw + NZno) * fbin[kt] #* Mweight/Msummed
                
            #if(clt[kt] == "young"):                
            #    eta = 0.5 * eta                
            #    if(Z[i] > 0.01):                    
            #        eta = 0.2 * eta
            

        eta_eff[kt][i] = eta        

        
        
        fstr += "{0:1.3e} {1:1.4e} {2:s}\n".format(dcen, eta_eff[kt][i], cluster)
    fstr += "\n\n"

f.write(fstr)
f.close()



lenz = len(z)-1
zfin_f = [[None]*lenz,[None]*lenz,[None]*lenz,[None]*lenz]
Rate_f = [[None]*lenz,[None]*lenz,[None]*lenz,[None]*lenz]
Numb_f = [0,0,0,0]
Int = [0,0,0,0]

for kt in range(len(clt)):

    cluster = clt[kt]
    typ = sfrt[kt]

    ## mergers in given environment ##
    idx_all = np.where(ctype == cluster)

    zmer_all= zmer[idx_all]
    zfor_all= zfor[idx_all]
    Zmet_all= Zmet[idx_all]
    tmer_all= tmer[idx_all]
       
    Numb_f[kt] = len(zmer_all)

    zmob_Ztot=[None]*len(Z)
    tmob_Ztot=[None]*len(Z)
    zmgw_Ztot=[None]*len(Z)
    Nztot = [0] * len(Z)
    
    for k in range(len(Z)):
        dcen = Z[k]
        if(k < len(Z)-1):
            dup = Z[k+1]
        else:
            dup = Z[k] * 1.1
            
        if(k > 0):
            dlow = Z[k-1]
        else:
            dlow = 0.9 * Z[k]
                
        dZsx = (dcen-dlow)/2.0
        dZdx = (dup-dcen)/2.0
        
        ### All mergers with a given metallicity    
        idx_Z = np.where((Zmet_all >= Z[k] - dZsx) & (Zmet_all < Z[k] + dZdx))    
        zmob_Ztot[k] = zfor_all[idx_Z]
        tmob_Ztot[k] = tmer_all[idx_Z]
        zmgw_Ztot[k] = zmer_all[idx_Z]
        Nztot[k] = len(tmob_Ztot[k])
        ###
 

    
    fsca = msca(cluster, typ)
    
    N = 0
    Rz = [None] * (len(z) - 1)
    zcen = [None] * (len(z)-1)
    tcen = [interp(z[i], redshift, lookback)*1.E9 for i in range(len(z))]

       

    for i in range(len(z)-1):
        zcen[i] = 0.5*(z[i]+z[i+1])
        tlook = tcen[i+1] - tcen[i]
        Int = 0.0
        fsum = [0.0]*len(Z)
        Iz = [0.0]*len(Z)        
        
        zp = np.linspace(z[i], maxz, 100)
        tp = [np.interp(zp[j], redshift[::-1],lookback[::-1]) for j in range(len(zp))]
        tp = np.array(tp)
        
        zp_cen = [0.5*(zp[j]+zp[j+1]) for j in range(len(zp)-1)]
        zp_cen = np.array(zp_cen)
        dtdz = DtDz(zp_cen, H0, Om, Ol)
        
        Rp_int = [None]*len(zp_cen)
        
        for k in range(len(Z)):
            dcen = Z[k]
            if(k < len(Z)-1):
                dup = Z[k+1]
            else:
                dup = Z[k] * 1.1
                
            if(k > 0):
                dlow = Z[k-1]
            else:
                dlow = 0.9 * Z[k]
                
            dZsx = (dcen-dlow)/2.0
            dZdx = (dup-dcen)/2.0

            ### All mergers with a given metallicity                
            zmob_Z = zmob_Ztot[k]
            tmob_Z = tmob_Ztot[k]
            zmgw_Z = zmgw_Ztot[k]
            #Nz = Nztot[k]
            ###
            
            I = [0.0]*(len(zp)-1)
            for j in range(len(zp)-1):

                idx_p = np.where((zfor_all >= zp[j]) & (zfor_all < zp[j+1]))
                tmer_p = tmer_all[idx_p]
                Nzp = len(tmer_p)
                #Nzp = interp(zp_cen[j], zcen, NzpT)
                
                ### All mergers formed in the zp bin with metallicity Z 
                idx = np.where((zmob_Z >= zp[j]) & (zmob_Z < zp[j+1]))        
                tmer_z = tmob_Z[idx]
                zmer_z = zmgw_Z[idx]
                ###

                ### All mergers merging in the z bin, formed at zp bin, with metallicity Z
                id_H = np.where((zmer_z > z[i]) & (zmer_z < z[i+1]))
                tmer_H = tmer_z[id_H]
                ###
                
                dN = len(tmer_H)*1.0
                Nz = len(zmer_z)*1.0

                #### Normalization a la Mapelli et al, but without the metallicity check ####
                #Nz = Nztot[k]
                #### Normalization a la Arca Sedda et al ####
                #Nz = Nzp
                
                frac = 0.0
                Nz = Nzp
                if(Nz > 0):
                    frac = dN / Nz 
                                   
                fsum[k] = frac * eta_eff[kt][k]

                Bf_fd = interp(zp_cen[j], z, Bf[kt])
                sfr = Bf_fd * psi(typ, zp_cen[j], fsca[0], fsca[1], fsca[2], fsca[3])                             
                
                I[j] = sfr * dtdz[j] * fsum[k] * gpc_to_mpc3 / tlook
        
            Iz[k] = scipy.integrate.simpson(I, zp_cen) 
            Int += Iz[k]
            
        Rz[i] = Int        

        zfin_f[kt][i] = zcen[i]
        Rate_f[kt][i] = Rz[i]

        #if(kt == 3):
        #    print( zcen[i], Rz[i] )


    #stri += "\n\n"

    #print(cluster, scipy.integrate.simpson(Rz, zcen))

totN = sum(list(Numb_f))
print("\ntotal number of sources = ",totN)


print("\n-------  Merger rate at z = 0.2 in yr^{-1} Gpc^{-3} --------\n")
print("LIGO estimate: (17.9 - 44.0)\n")
LocRate = 0.0
for kt in range(len(clt)):
    LocR = interp(0.2, zfin_f[kt], Rate_f[kt])
    LocRate += LocR
    print("BPOP "+clt[kt]+":", LocR)
    
Ratestr = "{0:1.3f}".format(LocRate)

print("\nBPOP estimate:", Ratestr, "\n")
print("------------------------------------------------------------")


LVK_min = 17.9
LVK_max = 44.0
LVK = 0.5*(17.9 + 44)
LVK_c = [0]*len(clt)
RBpop = [0]*len(clt)
lvkcorr = "no"

for kt in range(len(clt)):    
    f = open("Merger_rates_"+simname+a+"_"+clt[kt]+".txt", "w")
    stri = ""
    LVK_c[kt] = LVK * Numb_f[kt]/totN ## This tells us who contributes to the rate
    if(len(Rate_f[kt]) < 2):
        RBpop[kt] = Rate_f[kt]
    else:
        RBpop[kt] = interp(0.2, zfin_f[kt], Rate_f[kt]) ## This is the rate we get from sim
    
    for i in range(len(Rate_f[kt])):
        if(lvkcorr == "yes"):
            Rate_f[kt][i] = Rate_f[kt][i] * LVK_c[kt]/RBpop[kt]
            
        stri += "{0:1.3f} {1:1.4e}\n".format(zfin_f[kt][i], Rate_f[kt][i])


    f.write(stri)
    f.close()


print("Fraction of mergers from different channels        \n ")
fname = simname+a
GlobRate(fname, maxz, clt, lab, H0, Om, Ol)

end = time.time()

print("Elapsed time [s] = ", "{0:1.6f}".format(end-start))

print("\n=========================  DONE  ===========================\n")

