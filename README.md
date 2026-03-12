## <span style="color:red"> B-POP: Binary black hole POPulation synthesizer

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

- ["retention_Z0_#.dat"](https://github.com/marcasedda/BPOP/blob/main/struct_retZ.md) contains the cumulative distribution of natal kicks of BHs from single and binary stars. The distribution is calculated for different values of the metallicity

**3) Tools**

 Directory tools/ contains two bash scripts:
- bpop_chunks.sh\
the script runs B-POP on multiple threads, and store results from each threads in separated directories (chnk#)\
to perform multiple runs on multiple threads, initialise the ```N``` variable in the script, then simply type
```
sh bpop_chunks.sh
```

- bpop_catcatalogues.sh\
after a succesfull B-POP run on multiple threads, your directory will contain N directories named chnk#\
this script collects and collates useful files from chnk# directories in the main directory\
```
sh bpop_catcatalogues.sh
```

**4) Utilities**

**Merger rate calculation and mock merger catalogue**
NOTE: The user can set a desired redshift distribution or metallicity distribution for the sample, but the code alone does not take into account the possible impact of 
merger efficiency related to different channels or metallicity variations. To facilitate the inclusion of these effects, we have created the BPOPRate package.

The directory /utils/ contains two sub-directories:
- BPOP_LIBPY Contains libraries

- BPOP_MERGER_RATE Contains the Python script (BPOPRate_v2025.py) that calculate the MRD and retrieve the corresponding source sample.

*Use of BPOPRate package*

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
Once the initial parameters are set, type
```
python BPOPRate.py > Output.txt &
```

**NOTE**: Make sure that the directory above BPOP_MERGER_RATE contains the following files from B-POP: Catalogue.txt, Catalogue_multiple_dyn.txt, Larger_than_tH.txt, and all
retention_Z0_#.dat files

BPOPRate produces the following files (# identifies the adopted model):
- "efficiency_Z_SET#_massive_no.txt" contains the metallicity (Col. 1) and merger efficiency (Col. 2) (in units of Msun^-1) of isolated, young, globular, and nuclear clusters, separated by white blocks. 
- "Merger_rates_#_##_massive_no.txt" contains the redshift (Col. 1), merger rate density (Col. 2), and merger rate density for different values of the metallicity (Col. 3-14) for different environments, either isolated binaries (##="none"), YCs (##="young"), GCs (##="globular"), or NCs (##="nuclear")
- "Catalogue_##.txt" contains BBH mergers for the whole cosmic population (##="global") and different environments as above. The number of BBHs depends on the selected time window. Files' header:\
```Metallicity zGW zfor m1 m2 a1 a2 mrem arem xeff vGW nrecy cluster mclu rclu ecce semiax status cos1 cos2 cos12 id```

**NOTE**: to downsample properly a catalogue, candidates should be extracted separately from different environments' files. \
Example: to obtain 1yr of events from a catalogue containing 10yr of events, you must extract 1/10 of all BBHs in each catalogue (Catalogue_none.txt, Catalogue_globular.txt, Catalogue_nuclear.txt, and Catalogue_young.txt)

**A) The structure of the input file**

Check [input](https://github.com/marcasedda/BPOP/blob/main/input.md) for an example

```CHUNKS``` Number of parallel simulations\
```Hubble``` Hubble time in yr\
```PREDIR``` path to the BPOP executable\
```SINGPTH``` path to the stellar evolution tracks for single BHs\
```PATH``` path to the stellar evolution tracks for isolated BBH mergers\
```SEVN``` if "yes" it takes SEVN catalogues and metallicity range, if "no" it uses MOBSE catalogues\
```MERGER_EFF_CORR``` Leave "no"\

```N``` Number of BBH to be simulated (NOTE: not all binaries merge)\
```mmax``` Max stellar mass  (UNUSED)\
```mmin``` Min stellar mass for BH progenitors (UNUSED)\
```mslope``` slope of the mass function (UNUSED)\
```Zsun``` Solar metallicity\

```DynOvTot``` Fraction of dynamical binaries over the total\
```pYC``` Fraction of dynamical binaries from YCs\
```pGC``` Fraction of dynamical binaries from GCs\
```pNC``` Fraction of dynamical binaries from NCs\

```uppergap``` If "yes", B-POP allows for the formation of BHs with a mass > ```mass_gap``` from dynamical processes, depending on the ```upgtp``` parameter. It extracts BHs according to a power-law with a probability ```fupgp``` and slope ```a_gp```\
```mass_gap``` Lower bound of the PISN gap\
```fupgp```    Probability that a BH mass is extracted from a power-law in the upper mass gap\
```a_gp```     Power-law index for mass sampling of upper mass-gap BHs from dynamical processes\
```upgtp``` if "dicarlo" it bypasses ```uppergap``` and samples BHs with a mass in the gap, combining Di Carlo et al (2021) and Arca Sedda et al (2024) results

```bhseed``` Enables formation of IMBH seeds from stellar collisions according to different models\
```bifZ``` maximum metallicity below which IMBH seeds from stellar collisions are included\
```bhpisn``` minimum stellar mass above which an IMBH from direct collapse can form (outer edge of the PISN gap)\
```f_seed``` probability for IMBH seed formation (depending on ```bhseed```)\
```maxseed```  maximum IMBH seed mass allowed (depending on ```bhseed```)\
```minseed```  minimum IMBH seed mass allowed (depending on ```bhseed```)\
```seedslope``` the IMBH seed mass is sampled from a power-law with slope ```seedslope``` (depending on ```bhseed```)\

```SFRTYPE_ISO``` distribution of formation redshift for isolated BBHs\
```SFRTYPE_CLU``` distribution of formation redshift for dynamical BBHs

```mixing``` fraction of BHs from single stars or binaries\
```fbin``` binary fraction (it sets the probability that a BBH is formed through three-body or binary--single interactions)

```TagR``` model for cluster mass-radius relation

```SSMA``` dispersion of Maxwellian distribution used to sample BBH semimajor axes


```CLfill``` Cluster evolution model\
```MonoZ``` if "yes", a single metallicity is considered\
```mono_Z``` metallicity adopted for ```MonoZ "yes"```\
```cluster_test``` if "yes" only clusters of ```cluster_test_env``` are simulated
```cluster_test_env``` cluster type to be simulated when ```cluster_test``` is enabled
```msmbhmax``` maximum BH mass allowed, integration stops if a BH mass exceeds this threshold

```highgen``` if "yes", higher generation secondaries are allowed (NOT IMPLEMENTED YET)

```sigma_metal``` spread of the cosmological spin distribution
```sigma_distri``` model for metallicity distribution 

```spinlb``` model for natal spins

```INDEX_ALIGN``` index of the power-law function adopted to sample the angle between spins and angular momentum vectors

All other parameters are deprecated and should not be changed
