from astropy.cosmology import Planck18 as Planck
import astropy.units as u
from astropy.cosmology import z_at_value
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
import random
import os, shutil, glob
from os import path
import time
from math import erf
import sys
sys.path.insert(0, '/mnt/sdb/cugolini/hier/mix05_full_pop/noseed/BPOP_LIBPY/')
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from lib_BPOPRate import *

start = time.time()
Zsun = 0.019 #0.017

gpc_to_mpc3 = 1.e9

H0 = Planck.H0/u.km*u.Mpc*u.s
Om = Planck.Om0
Ol = 1.-Om
tHubble = Planck.age(0)/u.Gyr*1e9

test(gpc_to_mpc3, H0, Om, Ol)

plt.rc('text', usetex=False)
cmap = "viridis"
fsize = 41
lsize = 40
# family = 'Times New Roman'
matplotlib.rcParams.update({'font.size': fsize})
# matplotlib.rcParams.update({'font.family': family})
matplotlib.rcParams.update({'font.family': 'DejaVu Serif'})

print("\n===================  BPOPRate package  =====================")
## User-defined quantities: type of cluster, CSFH, binary fraction, typical mass, max redshift, and repository path and name ##
## SEE Notes below. Additional notes: add parameters for sfr coefficients (e.g. el-badry gaussians) and for additional sfr (e.g. elbadry for nuclear)
##

IDmod = 2
a = "a"
name_file = "../"
IBfrac   = 0.4
fysc_mw = 0.01
f_NCoc = 1.0
f_GCsc = 1.0
Tobs = 10.0

clt  = ["nuclear", "globular", "young", "none"]
sigma_Z = 0.2


maxz = 15
massive = "no"
repetita = "no"
evapora  = "no"
adjust_red = "no"
adjust_typ = "tmed" #tmed, zmed, rnd, lst, tave, zave, rate
pre = "SEVN"
aCE = ""
IBonly=False

#**************************************************************#

sfrt, metdis, finp = models(a)
print(sfrt, metdis, finp)

if(pre=="MOBSE"):
    IBcorr = 0.285
elif(pre=="SEVN" or pre=="SEVN_all"):
    IBcorr = 0.185
elif(pre=="IORIO"):
    IBcorr = 0.255
else:
    print("please select SEVN or MOBSE and retry")
    exit()
    

choices = ["tmed", "zmed", "rnd", "lst", "tave", "zave", "rate"]
if adjust_typ not in choices:
    print("adjust_typ",adjust_typ," not recognised")
    exit()

if(massive == "no"):
    max_mas = 1.E8
else:
    max_mas = 1.E2

fmix = -1
if(IDmod < 5):
    fmix = (1.0 * IDmod)/4.
elif(IDmod == 5):
    fmix = (1.0 * IDmod)/10.
elif(IDmod < 0):
    fmix = 0
else:
    print("Uknown type of model, check fmix")
    exit()
 
#-------------------------------------------------------------

mstar   = 1.0
fbh     = 0.0026
fbin_is = 1.0
fbin_yc = 0.6
fbin_gc = 0.2
fbin_nc = 0.2




'''
if(clurad[IDmod] == "Mapelli20"):
    mtyp = [pow(10.,6.18), pow(10.,5.6), pow(10., 4.3), 1.0]
else:
    mtyp = [pow(10.,6.36), pow(10.,5.01), pow(10., 3.53), 1.0]
'''

mtyp = [pow(10.,6.36), pow(10.,5.01), pow(10., 3.53), 1.0]
fbin = [fbin_nc, fbin_gc, fbin_yc, fbin_is]

#-------------------------------------------------------------

home=""
codename=[""]
simname=""
fl = name_file


print("\n File analysed: ",codename[0], fl)
print(simname+a, "\n")
print("Model ID", IDmod)
print("CLType: ", clt)
print("CSFH:   ", sfrt)
print("fbin:   ", fbin)
print("m_typ:  ", np.log10(mtyp))
print("z_max:  ", maxz)
print("f_ysc:  ", fysc_mw)
print("f_ib :  ", IBfrac)
print("f_mix:  ", fmix)
print("massive: ", massive, " - mmax = ", max_mas)
print("------------------------------------------------------------\n")

#### Calculations start here ####
if(pre=="SEVN"):
    Z = [
        0.0002,
        0.0003,
        0.0004,
        0.0007,
        0.0010,
        0.0014,
        0.002,
        0.004,
        0.007,
        0.010,
        0.014,
        0.02]
elif(pre=="IORIO"):
     Z = [
        0.0001,
        0.0002,
        0.0004,
        0.0006,
        0.0008,
        0.0010,
        0.002,
        0.004,
        0.006,
        0.008,
        0.010,
        0.014,
        0.017,
        0.02,
        0.03]    
elif(pre=="SEVN_all"):
     Z = [
        0.0002,
        0.0003,
        0.0004,
        0.0005,
        0.0007,
        0.0010,
        0.0014,
        0.002,
        0.003,
        0.004,
        0.005,
        0.007,
        0.010,
        0.014,
        0.02]    
else:        
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
app=""
if(pre=="SEVN" or pre=="SEVN_all"):
    app = "_SEVN"
elif(pre=="IORIO"):
    app = "_IORIO"
    
X = np.loadtxt("./eta_isolated_binaries"+app+aCE,usecols=[3,2],unpack=True)
zmer_mobse = X[0]
pmer_mobse = X[1]

if(zmer_mobse[0] != Z[0]):
    print("Error eta_isolated_binaries")
    exit()

# Directory and files to be analysed #

directory = [fl]
             
directory = [home + codename[i] + directory[i] for i in range(len(directory))]

file_all = "Larger_than_tH.txt"
file_mer = "Catalogue.txt"
file_esc = "retention_Z"
file_hie = "Catalogue_multiple_dyn.txt"

lab  = ["nuc", "glo", "you", "iso"]

z = np.linspace(0,maxz,500) #ogspace(-2, np.log10(maxz), 100) #0.01,maxz,100)
z = np.array(z)

# The following vector contains scaling for NC, GC, YC, IB, respectively
Bf   = [[f_NCoc]*len(z), [f_GCsc]*len(z), [1.0]*len(z), [1.0]*len(z)]
if(IBonly):
    Bf =  [[0.0]*len(z), [0.0]*len(z), [0.0]*len(z), [1.0]*len(z)]

## Preliminarly, we need to evaluate the scaling between the SFR among different channels ##


sfr   = [None]*len(z)
sfr_N = [None]*len(z)
sfr_G = [None]*len(z)
sfr_Y = [None]*len(z)
sfr_I = [None]*len(z)

for i in range(len(z)):

    sfr_tot =  psi("mf17", z[i], 1.0, 1.0, 1.0, 1.0)
    sfr[i] = sfr_tot
    
    ## Evaluate the scale factor w.r.t. the SFR ##
    fsca   = msca(clt[0], sfrt[0], finp[0])
    sfr_nc = psi(sfrt[0], z[i], fsca[0], fsca[1], fsca[2], fsca[3])
    sfr_N[i] = sfr_nc    
    
    fsca   = msca(clt[1], sfrt[1], finp[1])
    sfr_gc = psi(sfrt[1], z[i], fsca[0], fsca[1], fsca[2], fsca[3])
    sfr_G[i] = sfr_gc
    
    B_yc = max(0, min(fysc_mw, 1.-(sfr_nc+sfr_gc)/sfr_tot))
    sfr_yc = B_yc * psi(sfrt[2], z[i], 1.0, 1.0, 1.0, 1.0)
    sfr_Y[i] = sfr_yc    
        
    B_is = max(0, min(IBfrac, 1.0 - (sfr_yc+sfr_gc+sfr_nc)/sfr_tot))            
    sfr_is = B_is * psi(sfrt[3], z[i], 1.0, 1.0, 1.0, 1.0) 
    sfr_I[i] = sfr_is
    
    Bf[2][i] = B_yc
    Bf[3][i] = B_is


#for i in range(len(z)):
#    print("star formation rate = ",sfr[i], sfr_I[i]+sfr_Y[i]+sfr_G[i]+sfr_N[i])


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

print("--------------- Local denstity in Mpc^{-3} -----------------")

for kt in range(len(clt)):
    
    cluster = clt[kt]
    typ = sfrt[kt]
    fsca   = msca(cluster, typ, finp[kt])
    
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
rt = np.loadtxt("redshift_time.txt", usecols=[0,1,2], unpack=True)
redshift = rt[0]
tage     = rt[1]
lookback = rt[2]

for i in range(len(redshift)):
    tage[i] = Planck.age(redshift[i])/u.Gyr
    lookback[i]=Planck.lookback_time(redshift[i])/u.Gyr



## Reading the file ## 

Fn = directory[0] + file_mer
# # X = np.loadtxt(Fn, dtype={'names':('zmer', 'zfor', 'met', 'clus', 'mclu', 'vesc','mas', 'tmer', 'nrep', 'msec','rclu', 'stat','id', 'tfor'), 'formats':(float, float, float, 'U8',float,float,float,float,float,float,float,'U8',float,float)}, usecols=[24, 25, 1, 3, 15, 17, 5, 14, 2, 6, 16, 18, 0, 13], unpack=True)
# X = np.loadtxt(Fn, dtype={'names': ('zmer','zfor','met','clus','mclu','vesc','mas','tmer','nrep','msec','rclu','stat','id','tfor'), 'formats': (float, float, float, 'U8', float, float, float, float, float, float, float, 'U8', float, float)}, usecols=[33, 34, 1, 4, 16, 18, 6, 15, 2, 7, 17, 19, 0, 14], unpack=True)

# zmer_a = X[0]
# zfor_a = X[1]
# Zmet_a = X[2]
# ctyp_a = X[3]
# mclu_a = X[4]
# vesc_a = X[5]
# mpri_a = X[6]
# msec_a = X[9]
# tmer_a = X[7]
# nrep_a = X[8]
# rclu_a = X[10]
# stat_a = X[11]
# id_a   = X[12]
# tfor_a = X[13]

# cols_final = "#ID Metal Nrec Nsec EnvType lab m1[Msun] m2[Msun] a1 a2 Mfin[Msun] afin xeff vGW[km/s] tfor[yr] tlast_mer[yr] Mclu_t0[Ms] Rclu_t0[pc] Vesc[km/s] BinaryStatus aeje[AU] aGW[AU] n_BHs_t0 nBHs_tot nBHs_1g nBHs_2g nBHs_3g nBHs_4g nBHs_5g nBHs_>5g nmer_tot Mcore_th[Ms] rcore_th[pc] redshift_merger redshift_formation tSMBH[yr] redshiftSMBH mprog[Ms] eccentricity semimajoraxis[AU] acrit[AU] tmerger[yr] cos(angle_s1s2) cos(angle_s1L) cos(angle_s2L)".split()

# with open(Fn, 'r') as f:
#     lines = [l for l in f if not l.startswith('#') and l.strip()]

# max_cols = len(cols_final)
# rows = []
# for line in lines:
#     fields = line.split()
#     fields += [np.nan] * (max_cols - len(fields))
#     rows.append(fields[:max_cols])

# df = pd.DataFrame(rows, columns=cols_final)
# for col in df.columns:
#     try:
#         df[col] = pd.to_numeric(df[col])
#     except (ValueError, TypeError):
#         pass

df = read_catalogue(Fn)

zmer_a = df['zgw'].values
zfor_a = df['zfo'].values
Zmet_a = df['Z'].values
ctyp_a = df['ctype'].values
mclu_a = df['mcl'].values
vesc_a = df['vesc'].values
mpri_a = df['m1'].values
msec_a = df['m2'].values
tmer_a = df['tmer'].values
nrep_a = df['nrec'].values
nsec_a = df['Nsec'].values
rclu_a = df['rcl'].values
stat_a = df['status'].values
id_a   = df['id'].values
tfor_a = df['tfor'].values

Ln = directory[0] + file_all
Y = np.loadtxt(Ln, dtype={'names':('met','clus','nrecy','tmer','tfor'), 'formats':(float,'U8',float,float,float)}, usecols=[13,14,22,23,24], unpack=True)                



if(adjust_red == "yes"):

    stri = ""
    
    Gn = directory[0] + file_hie
    # K = np.loadtxt(Gn, usecols = [0,1,2,3,15,15,30], unpack=True)
    K = np.loadtxt(Gn, usecols=[0, 1, 2, 3, 15, 15, 31], unpack=True)
    
    m1rep = K[0]
    m2rep = K[1]
    a1rep = K[2]
    a2rep = K[3]
    trep = K[4]
    zrep = np.array([interp(trep[i]/1.E9, tage, redshift) for i in range(len(trep))])
    # irep = K[6]+1
    irep = K[6]
    
    idx_rep = np.where(nrep_a > 1)[0]
    
    tmedian = 0
    zmedian = 0

    
    for i in range(len(idx_rep)):
        idx_cat = np.where(irep == id_a[idx_rep[i]])[0]

       
        if(adjust_typ == "tmed"):
            tmedian = np.median(trep[idx_cat])        
            iclose = (np.abs(trep[idx_cat] - tmedian)).argmin()
        elif(adjust_typ == "zmed"):
            zmedian = np.median(zrep[idx_cat])
            iclose = (np.abs(zrep[idx_cat] - zmedian)).argmin()
        elif(adjust_typ == "tave"):
            tmedian = np.mean(trep[idx_cat])        
            iclose = (np.abs(trep[idx_cat] - tmedian)).argmin()
        elif(adjust_typ == "zave"):
            zmedian = np.mean(zrep[idx_cat])
            iclose = (np.abs(zrep[idx_cat] - zmedian)).argmin()
        elif(adjust_typ == "rnd"):
            iclose = np.random.randint(0, len(idx_cat)-1)
        elif(adjust_typ == "lst"):
            iclose = len(idx_cat)-1
        elif(adjust_typ == "rate"):
            iclose = np.argmin(np.ediff1d(trep[idx_cat]))
        else:
            print("ERROR")
            exit()

        
    
        mpri_a[idx_rep[i]] = m1rep[idx_cat[iclose]]
        msec_a[idx_rep[i]] = m2rep[idx_cat[iclose]]
        zmer_a[idx_rep[i]]  = zrep[idx_cat[iclose]]
        tmer_a[idx_rep[i]]  = trep[idx_cat[iclose]]

        
        stri += "{0:1.4e} {1:1.5f} {2:1.5e} {3:1.5e} {4:1.5f} {5:1.5e}".format(tmedian, zmedian, mpri_a[idx_rep[i]], msec_a[idx_rep[i]], zmer_a[idx_rep[i]], tmer_a[idx_rep[i]])

    fadj_nm = adjust_red + adjust_typ + "_SET"+str(IDmod)+".txt"    
    fadj=open(fadj_nm, "w")
    fadj.write(stri)
    fadj.close()
    

mthi_a = [max(mpri_a[i], msec_a[i]) for i in range(len(mpri_a))]
mpri_a = np.array(mthi_a)

if(evapora == "no"):
    if(massive == "no"):
        idx = np.where((zmer_a > 0.) & (mpri_a < max_mas)) 
    elif(massive == "gwtc3"):
        idx = np.where((zmer_a > 0.) & (mpri_a < 100.) & (mpri_a > 50.)) 
    else:
        idx = np.where((zmer_a > 0.) & (mpri_a > max_mas)) 
else:
    if(massive == "no"):
        idx = np.where((zmer_a > 0.) & (mpri_a < max_mas) & (stat_a != "inevap"))
    elif(massive == "gwtc3"):
        idx = np.where((zmer_a > 0.) & (mpri_a < 100.) & (mpri_a > 50.) & (stat_a != "inevap")) 
    else:
        idx = np.where((zmer_a > 0.) & (mpri_a > max_mas) & (stat_a != "inevap"))

tmer = tmer_a[idx]
zmer = zmer_a[idx]
zfor = zfor_a[idx]
Zmet = Zmet_a[idx]
ctype= ctyp_a[idx]
nrep = nrep_a[idx]
nsec = nsec_a[idx]
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
tloss= Y[3]
tgwloss= Y[4]
zloss= [interp(tloss[i]/1.E9, tage, redshift) for i in range(len(tloss))]
zloss= np.array(zloss)



f = open("efficiency_Z_SET"+str(IDmod)+a+"_massive_"+massive+".txt","w")
fstr = ""
for kt in range(len(clt)):         
            
    cluster = clt[kt]
    typ = sfrt[kt]

    ## Loss fraction -- to calculate the eta(Z) parameter ##

    dZ = 0.0001

    for i in range(len(Z)):
        dcen = Z[i]

        if(cluster == "none"):
            eta = IBcorr * interp(dcen, zmer_mobse, pmer_mobse)

        else:
            dfile_esc = directory[0]+file_esc+str(Z[i])+".dat"
            if(os.path.isfile(dfile_esc) == False):
                eta_eff[kt][i] = 0.0
                print(dfile_esc)
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

            
            
            if(repetita == "no"):
                nrep_sum = 0.0
            
            if(mstar != 0.0 and (NZgw+NZno != 0) and Msummed!=0):
                eta = fbh / mstar * (NZgw + nrep_sum) / (NZgw + NZno) * fbin[kt] * Mweight/Msummed
            

        eta_eff[kt][i] = eta        

        
        
        fstr += "{0:1.3e} {1:1.4e} {2:s}\n".format(dcen, eta_eff[kt][i], cluster)
    fstr += "\n\n"

f.write(fstr)
f.close()


def prob(Zcx, xmean, xsigma):
    return 0.5*(1.+erf((Zcx - xmean)/(np.sqrt(2.)*xsigma)))

fname = np.array(["Merger_rates_"+simname+a+"_"+clt[kt]+"_massive_"+massive+".txt" for kt in range(len(clt))])

lenz   = len(z)-1
zfin_f = [[None]*lenz,[None]*lenz,[None]*lenz,[None]*lenz]
Rate_f = [[None]*lenz,[None]*lenz,[None]*lenz,[None]*lenz]
Numb_f = [0,0,0,0]
Int    = [0,0,0,0]

for kt in range(len(clt)):

    if(os.path.isfile(fname[kt])):
        print("Warning: merger rate file for ", clt[kt], " already exist, skipping this")
        continue
    

    
    cluster = clt[kt]
    typ = sfrt[kt]

    ## mergers in given environment ##
    idx_all = np.where(ctype == cluster)

    print("Analysis of ", cluster)
    
    zmer_all= zmer[idx_all]    
    zfor_all= zfor[idx_all]
    Zmet_all= Zmet[idx_all]
    tmer_all= tmer[idx_all]
    Numb_f[kt] = len(zmer_all)


    ## binaries that don't merge for a given environment ##
    ## we need to know how many don't merge at a formation redshift ##
    idx_non = np.where(Closs == cluster)
    zlos_non= zloss[idx_non]
    Zlos_non= Zloss[idx_non]
    
    zmob_Ztot=[None]*len(Z)
    tmob_Ztot=[None]*len(Z)
    zmgw_Ztot=[None]*len(Z)
    Nztot = [0] * len(Z)

    znon_Ztot=[None]*len(Z)
    
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

        ### All non-mergers with a given metallicity 
        idn_Z = np.where((Zlos_non >= Z[k] - dZsx) & (Zlos_non < Z[k] + dZdx))
        znon_Ztot[k] = zlos_non[idn_Z]

        
    fsca = msca(cluster, typ, finp[kt])
    
    N = 0
    Rz = [None] * (len(z) - 1)
    zcen = [None] * (len(z)-1)
    tcen = [interp(z[i], redshift, tage)*1.E9 for i in range(len(z))]       

    if(len(zmer_all) == 0):
        for i in range(len(z)-1):
            zfin_f[kt][i] = 0.0
            Rate_f[kt][i] = 0.0
        continue

    fout = open("Merger_rates_"+simname+a+"_"+clt[kt]+"_massive_"+massive+".txt", "w")

    stri="#redshift Rate "
    for k in range(len(Z)):
        stri+="{0:1.5f} ".format(Z[k])
    stri+="\n"
    fout.write(stri)
    stri=""
    
    for i in range(len(z)-1):
        zcen[i] = 0.5*(z[i]+z[i+1])
        
        tlook = tcen[i] - tcen[i+1]
        Int = 0.0
        fsum = [0.0]*len(Z)
        Iz = [0.0]*len(Z)        
        
        zp = np.linspace(zcen[i], maxz, 500)
        
        
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

            
            ###
            ### All non-mergers with a given metallicity
            znon_Z = znon_Ztot[k]
            ###

            I = [0.0]*(len(zp)-1)
            for j in range(len(zp)-1):

                ### All mergers formed in the zp bin 
                #idx_p = np.where((zfor_all >= zp[j]) & (zfor_all < zp[j+1]))
                #tmer_p = tmer_all[idx_p]
                #Nz = float(len(tmer_p))
                ### All non-mergers formed in the zp bin
                #idx_r = np.where((znon_Z >= zp[j]) & (znon_Z < zp[j+1]))
                #Nz += float(len(znon_Z[idx_r]))                
                ###
          
                ### All mergers formed in the zp bin with metallicity Z 
                idx = np.where((zmob_Z >= zp[j]) & (zmob_Z < zp[j+1]))        
                zmer_z = zmgw_Z[idx]
                tmob_z = tmob_Z[idx]
                Nz = len(zmer_z)

                idsx = np.where((znon_Z >= zp[j]) & (znon_Z < zp[j+1]))
                Nz+=len(znon_Z[idsx])
                ###
                
                ### All mergers merging in the z bin, formed at zp bin, with metallicity Z
                id_H = np.where((tmob_z >= tcen[i+1]) & (tmob_z < tcen[i])) #np.where((zmer_z > z[i]) & (zmer_z < z[i+1]))
                dN = float(len(zmer_z[id_H]))                
                ###       

                ### This is valid only if metcor is a la Mapelli+20 --- Add El-Badry2019 option for GCs (and for NCs?) ###
                ### Also, need to consider redshift with metallicity outside boundaries
                ### Something like:
                ### let's use as 'string' below the variable used to determine the sfr
                ### important NOTE: there must be coeherence between this and BPOP (or not?)
                ### if(clt[k] == "globular" and string == "EB19"): ### Mstar_mass = 1.E10
                ### lZ =  0.40 * (np.log10(Mstar_mass) - 10) + 0.67 * np.exp(-0.50*zp_cen[j]) - 1.04 *(1. + pow(zp_cen[j]/15,3)) [ElBadry+19, star phase]

                sZ = sigma_Z
                lZ = metcor(metdis[kt], zp_cen[j], sZ)
                        
                pzZsx = 0.5 * (1. + erf((np.log10((Z[k]-dZsx)/Zsun)-lZ)/(sZ*np.sqrt(2.)))) 
                pzZdx = 0.5 * (1. + erf((np.log10((Z[k]+dZdx)/Zsun)-lZ)/(sZ*np.sqrt(2.)))) 

                pdf = 1./(np.sqrt(2.*np.pi)*sZ) * np.exp(-pow(np.log10(Z[k]/Zsun)-lZ,2.)/(2.*sZ*sZ))

                pzZ = pzZdx - pzZsx

                frac = 0.0
                
                if(Nz > 0):
                    #mergers in z, born in zp, with Z, divided by
                    #total number of sources born in zp, multiplied by
                    #Delta_CDF from log-normal distri, divided by
                    #number of sources formed in zp with Z over number of sources formed in zp
                    #This corresponds to dN / Nz * pzZ

                    fcl = pzZ
                    frac =  (dN / Nz) * fcl
                    
                    

                fsum[k] = frac 

                
                Bf_fd = interp(zp_cen[j], z, Bf[kt])
                sfr = Bf_fd * psi(typ, zp_cen[j], fsca[0], fsca[1], fsca[2], fsca[3])

                

                            
                I[j] = sfr * eta_eff[kt][k] * fsum[k] * gpc_to_mpc3 / tlook  * dtdz[j]               
                
            Iz[k] = scipy.integrate.simpson(I, zp_cen) 
            

            
        Rz[i] = np.sum(Iz)

        zfin_f[kt][i] = zcen[i]
        Rate_f[kt][i] = Rz[i]

        stri = "{0:1.3f} {1:1.4e} ".format(zfin_f[kt][i], Rate_f[kt][i])
        for k in range(len(Z)):
            stri += "{0:1.5e} ".format(Iz[k])
        stri+="\n"
        
    
        fout.write(stri)
        fout.flush()

        Iz = [0.0]*len(Z)

    fout.close()      


totN = sum(list(Numb_f))
print("\ntotal number of sources = ",totN)


print("\n-------  Merger rate at z = 0.2 in yr^{-1} Gpc^{-3} --------\n")
print("LIGO estimate: (17.9 - 44.0)\n")
LocRate = 0.0
LocR = 0.0
for kt in range(len(clt)):
    if(np.count_nonzero(Rate_f[kt]) > 0):
        LocR = interp(0.2, zfin_f[kt], Rate_f[kt])
        
    LocRate += LocR
    print("BPOP "+clt[kt]+":", LocR)
    
Ratestr = "{0:1.3f}".format(LocRate)

print("\nBPOP estimate:", Ratestr, "\n")
print("------------------------------------------------------------")

global_rate=False
if(global_rate):
    print("Fraction of mergers from different channels        \n ")    
    GlobRate(fname, maxz, clt, lab, H0, Om, Ol, massive)


end = time.time()

print("MRD calculation - Elapsed time [s] = ", "{0:1.6f}".format(end-start))


#Need to pass path to the catalogue
#Need to pass path to the merger rate
#bsmplr(gpc_to_mpc3, H0, Om, Ol, IDmod, a, clt, massive, pre, IBonly, fdyn, fgc, fyc, fnc, spinI, spinD, clurad, bhseed, umg, highg, mix, vms)

start = time.time()

cat    = Fn
MRDfile= fname
bsmplr(cat, MRDfile, gpc_to_mpc3, H0, Om, Ol, clt, Z, massive, Tobs, max_mas)

strig = "SET{0:1.0f}{1:s}_IB{2:1.2f}_YC{3:1.0e}_GC{4:1.1f}_NC{5:1.1f}_stev{6:s}_Tobs{7:1.1f}".format(IDmod, a, IBfrac, fysc_mw, f_GCsc, f_NCoc, pre+aCE, Tobs)
if not(os.path.isdir(strig)):
        create_dir(strig)


print("\n====================== copying output directories ================\n")
files = glob.iglob(os.path.join("./", "SMP*"))
for file in files:
    if os.path.isdir(file):
        shutil.move(file, strig)
                   

end = time.time()
print("Catalogue creation - Elapsed time [s] = ", "{0:1.6f}".format(end-start))
print("\n====================== DONE ================\n")

'''
NOTES
a) Add an options to skip the recreation of MRD files if they already exists, specify they already exist, would you overwrite them?
b) Add an options to specify the path to the directories both in the main and in bsmplr
c) Remove the creation of directories
d) Add info about the selected stellar evolution
e) Add additional input about alphaCE
'''
