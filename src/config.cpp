#include "config.h"
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <cctype>
using namespace std;

// ==========================
// Default (copiati dal tuo input_params.h)
// ==========================
string CHUNKS = "yes";

string PREDIR   = "/home/cugolini/BPOP_new/dev/BPOP/";
string SINGPTH  = "../../../SEVN_Catalogs/rapid_M20_a1/";
string PATH     = "../../../SEVN_Catalogs/rapid_M20_a1/";
string SEVN     = "yes";
string MERGER_EFF_CORR = "no";

int    N     = 1000000;
double mmax  = 150.;
double mmin  = 18.5;
double mslope= -2.35;
double Zsun  = 0.019;

double DynOvTot = 0.9;
double pYC = 0.75, pGC = 0.15, pNC = 0.1;

string uppergap = "yes";
string bhseed   = "bifrost";
double bifZ     = 0.001;
double bhpisn   = 270.;
double fupgp    = 0.15;
double mass_gap = 60.0;
string upgtp    = "dicarlo";
string SFRTYPE_ISO = "continuous";
string SFRTYPE_CLU = "continuous";

double mixing = 0.5;
double fbin   = 0.0;

string TagR = "AS20";
double SSMA = 0.3;

string CLfill = "GG23";
string MonoZ  = "no";
double mono_Z = 0.0002;
string cluster_test = "no";
string cluster_test_env = "NC";
double msmbhmax = 5.0e7;

double f_seed = 0.2;
double maxseed = 1000.0;
double minseed = 100.0;
double seedslope = 3.0;

double a_gp = 2.0;

string highgen = "no";

double sigma_metal_iso = 0.2;
string sigma_distri_iso = "Mapelli";
double sigma_metal_clu = 0.2;
string sigma_distri_clu = "Elbadry19";

string spinlb = "maxwellian02";
double obslope = 0.0;

string sfr_only = "no";

string zdyn = "logflat";
double zslope = -1.5;

string MRATIO = "nouniform";
double MRATIO_SLOPE = -1.0;

string zdist = "logflat";
double zdivi = 0.0;
string correction = "no";
double slope = -1.5;

double INDEX_ALIGN = 8.0;

double MCL_popIII = 1.0e6;
double RCL_popIII = 1.0;
double MBH_popIII = 300.0;
double aBH_popIII = 0.0001;

int kpar = 20;
int spar = 5;

int bin_st = 50;
int nsize  = 80;
int tsize  = 80;
int vsize  = 80;

// ==========================
// Utilità parsing key=value
// ==========================
static inline string trim(string s){
    auto ws = [](int c){ return std::isspace(c); };
    s.erase(s.begin(), std::find_if_not(s.begin(), s.end(), ws));
    s.erase(std::find_if_not(s.rbegin(), s.rend(), ws).base(), s.end());
    if (!s.empty() && s.front()=='"' && s.back()=='"' && s.size()>1)
        s = s.substr(1, s.size()-2);
    return s;
}

static unordered_map<string,string> read_kv(const string& path){
    ifstream in(path);
    unordered_map<string,string> kv;
    string line;
    if (!in) {
        cerr << "Config \"" << path << "\" non trovata: uso i default.\n";
        return kv;
    }
    while (getline(in, line)){
        auto comment = line.find('#');
        if (comment != string::npos) line.resize(comment);
        auto eq = line.find('=');
        if (eq == string::npos) continue;
        auto key = trim(line.substr(0, eq));
        auto val = trim(line.substr(eq + 1));
        if (!key.empty()) kv[key] = val;
    }
    return kv;
}

// ==========================
// Caricamento
// ==========================
void load_config(const string& filename){
    auto kv = read_kv(filename);

    auto S = [&](const char* k, string& dst){ auto it=kv.find(k); if(it!=kv.end()) dst = it->second; };
    auto D = [&](const char* k, double& dst){ auto it=kv.find(k); if(it!=kv.end()) dst = stod(it->second); };
    auto I = [&](const char* k, int&    dst){ auto it=kv.find(k); if(it!=kv.end()) dst = stoi(it->second); };

    // SWITCH / MISC
    S("CHUNKS", CHUNKS);

    // DATAFILES
    S("PREDIR", PREDIR); S("SINGPTH", SINGPTH); S("PATH", PATH);
    S("SEVN", SEVN); S("MERGER_EFF_CORR", MERGER_EFF_CORR);

    // GLOBAL
    I("N", N); D("mmax", mmax); D("mmin", mmin); D("mslope", mslope); D("Zsun", Zsun);

    // DYNAMICAL FRACTIONS
    D("DynOvTot", DynOvTot); D("pYC", pYC); D("pGC", pGC); D("pNC", pNC);

    // UPPER GAP / SEEDS / SFR
    S("uppergap", uppergap); S("bhseed", bhseed); D("bifZ", bifZ);
    D("bhpisn", bhpisn); D("fupgp", fupgp); D("mass_gap", mass_gap);
    S("upgtp", upgtp); S("SFRTYPE_ISO", SFRTYPE_ISO); S("SFRTYPE_CLU", SFRTYPE_CLU);
    D("a_gp", a_gp);

    // MIXING / BIN FRAC
    D("mixing", mixing); D("fbin", fbin);

    // YC mass-size
    S("TagR", TagR); D("SSMA", SSMA);

    // CLUSTER EVOLUTION
    S("CLfill", CLfill); S("MonoZ", MonoZ); D("mono_Z", mono_Z);
    S("cluster_test", cluster_test); S("cluster_test_env", cluster_test_env);
    D("msmbhmax", msmbhmax);

    // BH SEED block
    D("f_seed", f_seed); D("maxseed", maxseed); D("minseed", minseed); D("seedslope", seedslope);

    // HIGH GENERATION SECONDARY BH
    S("highgen", highgen);

    // METALLICITY SPREAD
    D("sigma_metal_iso", sigma_metal_iso); S("sigma_distri_iso", sigma_distri_iso);
    D("sigma_metal_clu", sigma_metal_clu); S("sigma_distri_clu", sigma_distri_clu);

    // SPINS
    S("spinlb", spinlb); D("obslope", obslope);

    // STAR FORMATION
    S("sfr_only", sfr_only);

    // DYNAMICAL SCENARIO
    S("zdyn", zdyn); D("zslope", zslope);

    // DEPRECATED (compatibilità)
    S("MRATIO", MRATIO); D("MRATIO_SLOPE", MRATIO_SLOPE);

    // ISOLATED SCENARIO
    S("zdist", zdist); D("zdivi", zdivi); S("correction", correction); D("slope", slope);

    // MISC
    D("INDEX_ALIGN", INDEX_ALIGN);

    // Seeds - pop III
    D("MCL_popIII", MCL_popIII); D("RCL_popIII", RCL_popIII);
    D("MBH_popIII", MBH_popIII); D("aBH_popIII", aBH_popIII);

    // FILE STRUCTURE
    I("kpar", kpar); I("spar", spar);

    // SIZES
    I("bin_st", bin_st); I("nsize", nsize); I("tsize", tsize); I("vsize", vsize);

    cout << "Config caricata da \"" << filename << "\"\n";
}

// ==========================
// Autoload all'avvio
// ==========================
struct _ConfigAutoLoader {
    _ConfigAutoLoader(){
        ifstream here("run_bpop.txt");
        if (here.good()) load_config("run_bpop.txt");
        else {
            ifstream up("../run_bpop.txt");
            if (up.good()) load_config("../run_bpop.txt");
            else cerr << "Nessun run_bpop.txt trovato; uso i default.\n";
        }
    }
} _bpop_autoload;
