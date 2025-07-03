#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <string>
#include <optional>

using namespace std;

class Functions
{
 public:
  double geo16Rnc(double mlog, double c1, double c2, double a, double b, double e);
  double LOGSMP(double mean, double sigma);
  double metcor(string metal_dis, double sigmaZ, double red_del);
  double GSS_pdf_f(double x,double xm,double sx);
  double GSS_cdf_f(double x, double xmean, double xdisp);
  double GSS_smpl(double minZ, double maxZ, double xm, double sx);
  double mevol(double t, double rh, double mh, double trel, string type, string tclus);
  double revol(double t, double rh, double mh, double trel, string type, string tclus);

  std::optional<size_t> search_closest(const std::vector<double> & sorted_array, double x);
  double finterp(vector<double> &Zeta, vector<double> &Eeta, double lz);


  void test();
  double inter(double x, double *X, double *Y, int N);

  double sfr_red(string sfr);

  double maxwell(double a);
  double rnd();
  double rnd_old();
  double rndG();
  double rndG7();
  double rndgen(double pp, double spp);
  double LogGaussian(double pp, double spp);
  double Gaussian_normal(double min, double max, double pp, double spp);
  double Gss_weight(vector<double>& Zeta, vector<double>& Eeta, double Zsun, double zmean, double zsigma, double zmin, double zmax);
  double Gaussian(double pp, double spp);

  double tred(double z);
  double zred(double t);

  double phiP(double x,double M,double a);
  double phiD(double x,double M,double a,double g);
  double Rh(double a,double g);
  double ah(double Rh,double g);
  double vescape(double r, double m, string type);
  
  void histo(double* X, int N, int nbin, string binning, string name);

  double findbhmin(double metal, string singpth);
  
  double singBHt(double metal, string singpath, double obslope, double mslope);
  

  void singBHt_mix_old(vector<double> zams_mix, vector<double> remn_mix, vector<double> tdel_mix, vector<double> kick_mix, double obslope, double mslope, double *sing_out, double msmax, double msmin, double mbmax, double mbmin, double vescape);
  //void singBHt_new(vector<double> zams_sin, vector<double> remn_sin, vector<double> tdel_sin, vector<double> kick_sin, double obslope, double mslope, double *sing_out, double msmax, double msmin, double mbmax, double mbmin, double vescape);
  void singBHt_new(vector<double>& zams_sin, vector<double>& remn_sin, vector<double>& tdel_sin, vector<double>& kick_sin, double *sing_out, double vescape);


  double mratio(double mpri, double MRATIO_SLOPE, string stri_mrat);

  double sevntab(double M0, double metal, string path);

  int bbh_UIBfits_setup(double mbh1, double mbh2, double chibh1, double chibh2, double *out);

  double bbh_final_mass_non_precessing_UIB2016(double m1, double m2, double chi1, double chi2, string version);

  double bbh_final_spin_non_precessing_UIB2016(double m1, double m2, double chi1, double chi2, string version);

  double angmom(double a, double atot, double v, int Nmet);

  double angmom2(double a1, double a2, double v, double q, double cosa, double cosb, double cosg);

  void SREM2(double ndx, double a1, double a2, double m1, double m2, string align, double *spins);

  double kicks(double m1, double m2, double a1, double a2, double q, double cosa, double cosb, double cosg);

  double spin(double mass, string spinning);

  string print(double mass, double mmax, double mmin,double wgh);

  Functions();

};
#endif
