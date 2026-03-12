## <span style="color:red"> B-POP: Binary black hole POPulation synthesis

B-POP is a population synthesis code able to create synthetic Universes filled up with black hole binaries pairing and merging in different environments, namely the galactic field (as isolated binaries), 
and dynamical environments (young, globular, and nuclear star clusters). 

**1) Before compilation**

The code in its current version takes all parameters from a header file (soon to be changed into an ASCII input file). The parameter file is stored in the /include/ directory. The structure is shown in
section A) below

The main parameters must be set before compilation. As a best practice, for the time being, it is advisable to copy the input parameter.h file in the directory where BPOP 
is ultimately run, so as to keep track of the simulated setup. 
<u> Note that in the current version, the code needs to be recompiled whenever one parameter in the file must be changed. <u> 

**2) Compiling and running B-POP**
```
make clean && make
```

The executable is saved in /build/BPOP.x

To run the code, copy BPOP.x into the chosen directory and type 
```
./BPOP.x > Output.txt &
```

The Output.txt file will contain all printouts of the code. 
All files with *txt and *dat files produced by BPOP are ultimately stored in a new directory whose naming summarises some of the main parameters of the simulation.

The main files produced by a single run are:
- ["Catalogue.txt"](https://github.com/marcasedda/BPOP/blob/main/struct_Catalogue.md) This file contains all BBH mergers. It consists of 36 columns which represent:

- ["Catalogue_multiple_dyn.txt"](https://github.com/marcasedda/BPOP/blob/main/struct_CatMul.md) This file contains the merger trees of all multiple mergers. 

- ["Larger_than_tH.txt"](https://github.com/marcasedda/BPOP/blob/main/struct_LargtH.md) contains info about binaries that do not merge within a Hubble time


**3) Utilities and tools**

3.1) Merger rate calculation and mock merger catalogue
NOTE: The user can set a desired redshift distribution or metallicity distribution for the sample, but the code alone does not take into account the possible impact of 
merger efficiency related to different channels or metallicity variations. To facilitate the inclusion of these effects, we have created the BPOPRate package.

The directory /utils/ contains two sub-directories:
- BPOP_LIBPY Contains libraries

- BPOP_MERGER_RATE Contains the Python script (BPOPRate_v2025.py) that calculate the MRD and retrieve the corresponding source sample.

**Use of BPOPRate package**
The directory must be copied into the directory where BPOP is run. The package creates a synthetic Universe and fills it with BBH mergers drawn from BPOP catalogues.
The following parameters must be adjusted:

```
IDmod = 2 
a = "a"
name_file = "../"

IBfrac   = 0.4
fysc_mw = 0.005
f_NCoc = 1.0
f_GCsc = 1.0
Tobs = 1.0

maxz = 15
massive = "no"
repetita = "no"
evapora  = "no"
adjust_red = "no"
adjust_typ = "tmed" #tmed, zmed, rnd, lst, tave, zave, rate
pre = "SEVN"
aCE = ""
IBonly=False


```

3.2) Directory tools/ contains two bash scripts:
bpop_chunks.sh 



**A) The structure of the input file**

```
//input_params.h

//Guard
#ifndef PARAMETERS_H
#define PARAMETERS_H

#define CHUNKS "no"
#define Hubble 13.99E9 //13.803E9

// DATAFILES (Metal. distri, Single BHs, Binary BHs)
#define PREDIR "/home/manuel/Scrivania/ACTIVE_PROJECTS/BPOP/"
#define SINGPTH  "../SEVN_Catalog/rapid_M20" //"../MOBSE_Catalog/A5" //"../SEVN_Catalog/rapid_M20_a1/" 
#define PATH     "../SEVN_Catalog/rapid_M20" //"../MOBSE_Catalog/A5" //"../SEVN_Catalog/rapid_M20_a1/" 
#define SEVN "yes"
#define MERGER_EFF_CORR "no"

// GLOBAL
#define N        1000000
#define mmax     150.
#define mmin     18.5
#define mslope  -2.35
#define Zsun     0.019

//DYNAMICAL FRACTIONS
#define DynOvTot 0.95
#define pYC      0.8
#define pGC      0.1
#define pNC      0.1

#define uppergap "yes"
#define bhseed   "bifrost"
#define bifZ     0.001
#define bhpisn   270.
#define fupgp    0.15
#define mass_gap  60.0
#define upgtp   "dicarlo"
#define SFRTYPE_ISO "continuous" //"MF17" //"continuous" //"bigbang" //
#define SFRTYPE_CLU "continuous" //"EB18_MF17" //"KR13" //"EB18_MF17" //"continuous" //"bigbang" //


#define mixing  0.5
#define fbin    0.5

//YC mass-size relation
//Marks12, Rantala24, Mapelli20, AS20
//if None, sizes and masses are rescaled from GCs
#define TagR "AS20"

//Sigma of Gaussian sma distribution
#define SSMA 0.3

//CLUSTER EVOLUTION
#define CLfill       "GG23"
#define MonoZ        "no"
#define mono_Z       0.02
#define cluster_test "no"
#define cluster_test_env "NC"
#define msmbhmax 5.E7

//BH SEED
#define f_seed    0.2
#define maxseed   1000.0
#define minseed   100.0
#define seedslope 3.0

//HIGH MASS BHs
#define a_gp 2.0

//HIGH GENERATION SECONDARY BH
#define highgen "no"

//METALLLICITY SPREAD
#define sigma_metal 0.2
#define sigma_distri "Mapelli"

//SPINS
#define spinlb  "maxwellian02"
#define obslope  0.0

//STAR FORMATION
#define sfr_only "no"

//DYNAMICAL SCENARIO
#define zdyn       "logflat"
#define zslope     -1.5

//DEPRECATED
#define MRATIO "nouniform" //"uniform","pwl","gss","nouniform"
#define MRATIO_SLOPE -1.0


//ISOLATED SCENARIO
#define zdist       "logflat" //"gallazzi"
#define zdivi       0.0
#define correction  "no" 
#define slope       -1.5

#define INDEX_ALIGN 8.0

//FILE STRUCTURE
#define kpar 20
#define spar 5
#define numZ 12

//SIZE OF GENERAL VECTORS
#define bin_st 50
#define nsize 80
#define tsize 80
#define vsize 80

#endif

```
