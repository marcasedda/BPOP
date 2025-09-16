//input_params.h

//Guard
#ifndef PARAMETERS_H
#define PARAMETERS_H

#define CHUNKS "yes"
//#define Hubble 13.99E9 //13.803E9

// DATAFILES (Metal. distri, Single BHs, Binary BHs)
#define PREDIR "/home/manuel/Scrivania/ACTIVE_PROJECTS/BPOP/"
#define SINGPTH  "../SEVN_Catalog/rapid_M20/"
#define PATH     "../SEVN_Catalog/rapid_M20/" 
#define SEVN "yes"
#define MERGER_EFF_CORR "no"

// GLOBAL
#define N        1000000
#define mmax     150.
#define mmin     18.5
#define mslope  -2.35
#define Zsun     0.019

//DYNAMICAL FRACTIONS
#define DynOvTot 0.9
#define pYC      0.75
#define pGC      0.15
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
#define fbin    0.0

//YC mass-size relation
//Marks12, Rantala24, Mapelli20, AS20
//if None, sizes and masses are rescaled from GCs
#define TagR "AS20"

//Sigma of Gaussian sma distribution
#define SSMA 0.3

//CLUSTER EVOLUTION
#define CLfill       "GG23"
#define MonoZ        "no"
#define mono_Z       0.002
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
//NOTE: The choice Elbadry19 requires a sigma_metal_clu = 0.5 to embrace galaxy masses between ~ 10^8 - 10^12
#define sigma_metal_iso  0.2
#define sigma_distri_iso "Mapelli"
#define sigma_metal_clu  0.2
#define sigma_distri_clu "Elbadry19"

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
//#define numZ 12

//SIZE OF GENERAL VECTORS
#define bin_st 50
#define nsize 80
#define tsize 80
#define vsize 80

#endif
