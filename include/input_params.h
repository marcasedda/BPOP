#pragma once
#include <string>
#include "config.h"
using namespace std;

// === SWITCH / FLAGS ===
extern string CHUNKS;

// === DATAFILES ===
extern string PREDIR, SINGPTH, PATH;
extern string SEVN, MERGER_EFF_CORR;

// === GLOBAL ===
extern int    N;
extern double mmax, mmin, mslope, Zsun;

// === DYNAMICAL FRACTIONS ===
extern double DynOvTot, pYC, pGC, pNC;

// === UPPER GAP / SEEDS / SFR TYPES ===
extern string uppergap, bhseed, upgtp;
extern double min_hier;
extern double bifZ, bhpisn, fupgp, mass_gap, a_gp;
extern string SFRTYPE_ISO, SFRTYPE_CLU;

// === MIXING / BIN FRAC ===
extern double mixing, fbin;

// === YC MASS-SIZE RELATION ===
extern string TagR;
extern double SSMA;

// === CLUSTER EVOLUTION ===
extern string CLfill, MonoZ, cluster_test, cluster_test_env;
extern double mono_Z, msmbhmax;

// === BH SEED BLOCK ===
extern double f_seed, maxseed, minseed, seedslope;

// === HIGH GENERATION SECONDARY BH ===
extern string highgen;

// === METALLICITY SPREAD ===
extern double sigma_metal_iso, sigma_metal_clu;
extern string sigma_distri_iso, sigma_distri_clu;

// === SPINS ===
extern string spinlb;
extern double obslope;

// === STAR FORMATION ===
extern string sfr_only;

// === DYNAMICAL SCENARIO ===
extern string zdyn;
extern double zslope;

// === DEPRECATED (compatibilità) ===
extern string MRATIO;
extern double MRATIO_SLOPE;

// === ISOLATED SCENARIO ===
extern string zdist, correction;
extern double zdivi, slope;

// === MISC ===
extern double INDEX_ALIGN;

// === SEEDS - POP III ===
extern double MCL_popIII, RCL_popIII, MBH_popIII, aBH_popIII;

// === FILE STRUCTURE ===
extern int kpar, spar;

// === SIZE OF GENERAL VECTORS ===
extern int bin_st, nsize, tsize, vsize;
