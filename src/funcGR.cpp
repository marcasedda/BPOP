#include "funcGR.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <cstring>
#include <stdlib.h> 
#include <sys/types.h>  // For stat().
#include <sys/stat.h>   // For stat().
#include <string>
#include <random>
#include <algorithm>
#include <optional>

Functions::Functions(){
  srand(time(0));
}


optional<size_t> Functions::search_closest(const std::vector<double> & sorted_array, double x) {

    if (sorted_array.empty())
        return std::nullopt;

    const auto iter_geq = std::lower_bound(sorted_array.begin(), sorted_array.end(), x);

    if (iter_geq == sorted_array.begin())
        return 0;

    if (iter_geq == sorted_array.end())
        return sorted_array.size() - 1;

    const auto & a = *(iter_geq - 1);
    const auto & b = *(iter_geq);

    if (std::fabs(x - a) < std::fabs(x - b))
        return iter_geq - sorted_array.begin() - 1;

    return iter_geq - sorted_array.begin();
}

double Functions::finterp(vector<double> &Zeta, vector<double> &Eeta, double lz){
  int eta_id, eta_id2;
  eta_id = search_closest(Zeta, pow(10., lz)).value();
  double eta_m = Eeta[eta_id];
  eta_id2= eta_id;
  int lng = Zeta.size();
  if(Zeta[eta_id] > pow(10., lz))
    eta_id2 = std::max(0, eta_id2-1);      
  else
    eta_id2 = std::min(lng, eta_id2+1);
  
  double interp = (Eeta[eta_id] - Eeta[eta_id2])/(Zeta[eta_id] - Zeta[eta_id2]) * (pow(10., lz) - Zeta[eta_id2]) + Eeta[eta_id2];
  
  return interp;
}

double Functions::geo16Rnc(double mlog, double c1, double c2, double a, double b, double e){

  double lR = log10(c1) + a * (mlog - log10(c2)) + b ;

  lR = lR * (1.+ e * (1.-2.*rnd()));
  
  return lR;

}

double Functions::LOGSMP(double mean, double sigma){
  double logpdf, logrnd;
  logpdf = -1;
  logrnd = 0;
  double xmint = -1.0;
  double pdfmax = 1./(mean * sigma * sqrt(2.*M_PI));
  double xxmin = mean - 3.*sigma;
  double xxmax = mean + 3.*sigma;
  do{
    xmint = xxmin + (xxmax - xxmin)*rnd();
    logpdf = 1./(pow(10.,xmint)*sigma*sqrt(2.*M_PI))*exp(-pow(xmint - mean,2.)/(2.*sigma*sigma));
    logrnd = pdfmax * rnd();
  }while(logpdf < logrnd);
  
  return xmint;
}



double Functions::metcor(string metal_dis, double sigmaZ, double red_del){
      double logz;
      double logz_me;
      if(metal_dis == "Bouffanais")  
	logz_me = 0.153 - 0.074 * pow(red_del,1.34); // Bouffanais 2021
      else if(metal_dis == "Mapelli")
	logz_me = 0.153 - 0.074 * pow(red_del,1.34) - (log(10.0) * sigmaZ*sigmaZ/2.0); // Bavera et al 2020, Mapelli et al 2022
      else if(metal_dis == "Rafelski") //Giacobbo et al 2018, Rafelski et al 2012
	if(red_del <= 1.5)
	  logz_me = -0.19 * red_del ;
	else
	  logz_me = -0.22 * red_del ;
      else if(metal_dis == "Giacobbo20a")
	logz_me = -0.24*red_del;
      else if(metal_dis == "Giacobbo20b")
	logz_me = -0.24*red_del - 0.18;
      else if(metal_dis == "Elbadry19"){
	//Eq. 13 in EB19
	double Mstar_mass = 1.e10;
	//logz_me =  0.35 * (log10(Mstar_mass) - 10) + 0.93 * exp(-0.43*red_del) - 1.05 *(1. + pow(red_del/15,3)); //gas phase (see Ma et al 2016)
	logz_me =  0.40 * (log10(Mstar_mass) - 10) + 0.67 * exp(-0.50*red_del) - 1.04 *(1. + pow(red_del/15,3)); //star phase (see Ma et al 2016)
      }
      else{
	cout<<"Metallicity distribution not found, please retry"<<endl;
	exit(0);
      }

      return logz_me;
      
}

double Functions::GSS_pdf_f(double x,double xm,double sx){
  return 1./sqrt(2.*M_PI*sx*sx) * exp(-pow(x-xm,2.)/(2.*sx*sx));
}

double Functions::GSS_cdf_f(double x, double xmean, double xdisp){
  return 0.5*(1.+erf((x - xmean)/(sqrt(2.)*xdisp)));
}

double Functions::GSS_smpl(double minZ, double maxZ, double xm, double sx){

  double Zsun = 0.017;
  double logz_me = xm;
  double sigmaZ = sx;

  double pmax = GSS_pdf_f(logz_me, logz_me, sigmaZ);

  double pmax_pro = pmax;
  double ppro = 0.1*pmax_pro;
  double logz = 0.0;
  double logzmin=log10(1.E-6/Zsun);
  double logzmax=log10(0.3/Zsun);
  double logz_norm = logz;
  do{
    pmax_pro = pmax * rnd();
    logz_norm = logzmin + (logzmax - logzmin)*rnd();
    ppro = GSS_pdf_f(logz_norm,logz_me,sigmaZ);
  }while(ppro < pmax_pro);

  logz = logz_norm + log10(Zsun);

  if(logz < log10(minZ))
    logz = log10(minZ);
  if(logz > log10(maxZ))
    logz = log10(maxZ);

  
  return logz;
  
}

double Functions::Gaussian_normal(double min, double max, double pp, double spp){
  std::random_device rd;
  std::mt19937 mt(rd());
  std::normal_distribution<double> distG(pp,spp);
  
  double rndG = distG(mt);
  rndG = distG(mt);


  if(rndG < min)
    rndG = min;
  if(rndG > max)
    rndG = max;  
  
  return rndG;
}


double Functions::Gss_weight(vector<double>& Zeta, vector<double>& Eeta, double Zsun, double zmean, double zsigma, double zmin,double zmax){

  vector<double> cdf,pdf;
  for(int j=0;j<Zeta.size();j++)
    pdf.push_back(GSS_pdf_f(log10(Zeta[j]/Zsun), zmean, zsigma)*Eeta[j]);      
  
  cdf.push_back(pdf[0]);
  for(int j=1;j<Zeta.size();j++)
    cdf.push_back(cdf[j-1]+pdf[j]);
  
  double maxc = cdf[cdf.size()-1];
  for(int j=1;j<Zeta.size();j++){
    cdf[j] /= maxc;
  }
  
  int ln=cdf.size()-1;
  double p = rnd();
  int ida, ida2,idb,idb2;
  ida2 = search_closest(cdf, p).value();

  if(cdf[ida2] > p)
    idb2 = max(0,ida2-1);
  else
    idb2 = min(ln,ida2+1);

  ida = min(ida2, idb2);
  idb = max(ida2, idb2);
  
  double itp;
  if(ida==idb)
    itp = Zeta[ida];
  else
    itp = (Zeta[ida]-Zeta[idb])/(cdf[ida]-cdf[idb])*(p-cdf[ida]) + Zeta[ida];

  if(itp > zmax)
    itp = zmax;
  if(itp < zmin)
    itp = zmin;
    
  cdf.erase(cdf.begin(), cdf.end());
  
  return log10(itp);
}


double Functions::mevol(double t, double rh, double mh, double trel, string type, string tclus){

  double tcc = 0.138*mh/(150. * log(0.11*mh/150.))*sqrt(pow(rh*3.08E16,3.)/(6.67E-11*1.99E30*mh))/(365.*24.*3600.);        
  if(tcc < 0)
    tcc = 0.2 * 4.2E9 * pow(rh/4.0,1.5) * sqrt(mh/1.E7) ;

  
  double mfact_h ,  mfact;
  if(type == "GG23"){
    mfact_h = mhalf_evo( t,  rh,  mh,  trel,  type,  tclus ); 
    mfact = mfact_h * 0.03;
    
    double rcoll = revol(5.*tcc, rh, mh, trel, type, tclus);
    
    if(t < 5.*tcc){
      double rcoll = revol(t, rh, mh, trel, type, tclus);
      if(t < 2.*tcc)
	mfact *= pow(rcoll/0.15,0.4);
      else
	mfact *= pow(rcoll/0.15,1./1.3);
    }
    if(t >= 5.*tcc){
      mfact *= pow(rcoll/0.15,1./1.3);
    }
    
    
    if(mfact <= 0.0 || mfact_h <= 0.0)
      mfact = 0.0;
    
  }
  else if(type == "noevo"){
    mfact = 1.0;
  }
  else{
    cout<<"Select noevo, under, over, critical, GG23, or mix"<<endl;
    exit(0);    
  }

  return mfact;
}

double Functions::revol(double t, double rh, double mh, double trel, string type, string tclus){
  
  double tcc = 0.138*mh/(150. * log(0.11*mh/150.))*sqrt(pow(rh*3.08E16,3.)/(6.67E-11*1.99E30*mh))/(365.*24.*3600.);
  if(tcc < 0.0)
    tcc = 0.2 * 4.2E9 * pow(rh/4.0,1.5) * sqrt(mh/1.E7) ;

  double rfact_h;
  double rfact;
  double rmin = pow(5000. / mh , 1./0.43);
  
  if(type == "GG23"){
    rfact_h = rhalf_evo( t,  rh,  mh,  trel,  type,  tclus);
    double gamma = .25;//7;
    double keep = 5.;
    double fct = (1. + 0.25*(1.-2.*rnd()));
    double alfct = 0.2;
    
    if(t > keep*tcc){
      double beta = 0.05;  //pow(10., -3. + 2.*rnd());
      double r10 =  rmin + alfct*pow(keep/2. - 1.,gamma);
      rfact =  r10 + pow(t/(2.*tcc)-1.,beta) - pow(keep/2. - 1.,beta); //*rnd());
    }
    else if(t>2.*tcc && t<=keep*tcc){
      rfact =  (rmin + alfct*pow(t/(2.*tcc)-1.,gamma));//*rnd());
    }
    else{
      rfact = rfact_h * 0.2 * pow(1.-t/(2.*tcc),0.53); //* (0.05 + 0.25*rnd()) ;      
      //if(rfact < rmin)
      //rfact = rmin;
      rfact = max(rmin, rfact);		
    }

    rfact *= fct;
    
    stringstream chk;
    chk<<rfact;
    if(chk.str() == "nan" || rfact == 0.0){
      cout<< "Rhalf == 0: " <<rfact<<" "<<rfact_h<<" "<<t<<" "<<tcc<<" "<<mh<<" "<<rh<<endl;
      exit(0);
    }
    
  }
  else if(type == "noevo"){
    rfact = 1.0;
  }
  else{
    cout<<"Select noevo, under, over, critical, GG23, or mix"<<endl;
    exit(0);
  }

  
  //if(mfact == 0.0)
  //  rfact = 0.0;

  return rfact;
}



void Functions::test(){
  double prova[1000], provaold[1000];
  for(int i=0;i<1000;i++){
    prova[i] = 2.+4.*rnd();
    provaold[i]=2.+4.*rnd_old();
  }
  
  int ns = 20;
  double step = 4./(double)ns;
  double R[ns],N[ns],No[ns];
  for(int i=0; i<ns; i++){
    R[i] = 2.0 + i*step;
    N[i] = 0.0;
    No[i]= 0.0;
  }
  
  for(int j=0;j<1000;j++){
    for(int i=0; i<ns;i++)    {
      if(prova[j] >= R[i] && prova[j] < R[i+1]){
	N[i]++;
	break;
      }
    }
    for(int i=0; i<ns;i++)    {
      if(provaold[j] >= R[i] && provaold[j] < R[i+1]){
	No[i]++;
	break;
      }
    }      
  }
  
  for(int i=0; i<ns-1 ;i++)    
    cout<<0.5*(R[i]+R[i+1])<<" "<<N[i]<<" "<<No[i]<<endl;

  return;
  
}

double cor(double zx, double omM){
  return 1./((1.+zx)*sqrt(pow(1.+zx,3.)*omM + (1.-omM)));
}
double cov(double zx, double omM){
  double a = 4632.04;//          +/- 41.74        (0.9012%)
  double b = 0.379966;//        +/- 0.00475      (1.25%)

  double Dcom = a*zx/(1.+b*zx);
  double dVdz = Dcom * Dcom * cor(zx, omM) * (1.+zx);

  return dVdz;
}
    

  
double Functions::sfr_red(string sfrtype){

  string sfr = sfrtype;
  double zred = -1;
  double psisfr = 1.E30;
  double zredmax   = 1.08; //MF17 2.036;
  double Om = 0.3609;
  
  
  double psisfrmax = 0.01 * pow(1+zredmax,2.6) / (1. + pow((1+zredmax)/3.2,6.2)) ;/// cor(zredmax, Om);
  double psirnd;

  double reds_max = 20.0;

  
  if(sfr=="katz13" || sfr == "KR13"){
    zred = 2. + 8.*rnd(); //KATZ E RICOTTI 2013
  }
  else if(sfr=="madau17" || sfr == "MF17"){    

    do{
      zred = reds_max*rnd();
      psisfr = 0.01 * pow(1+zred,2.6) / (1. + pow((1+zred)/3.2,6.2)) ;/// cor(zred, Om);
      psirnd = psisfrmax * rnd();
      if(psisfr > psirnd)
	break;
    }while(psisfr < psirnd);

  }
  else if(sfr=="continuous" || sfr == "constant"){
    zred = 15.0*rnd(); //REFERENCE?
  }
  else if(sfr=="grid"){
    zred = 15.0*rnd();
  }
  else if(sfr=="burst"){
    zred = 20.0;
  }
  else if(sfr=="bigbang"){
    zred = 20.;
  }
  else if(sfr=="single"){
    zred = 20.;
  }
  else if(sfr=="elba18" || sfr=="EB18"){
    double Zn,Sn;
    
    Zn = 3.2;
    Sn = 1.5;
    
    double pup = 1.0;
    double pgas= 2.0;
    double zup;
    do{
      pup = rnd();
      zup = reds_max * rnd();
      pgas = exp(-pow(zup-Zn,2.) / (2.*Sn*Sn));
      if(pgas > pup)
	break;
    }while(pgas < pup);
    
    zred = zup;
    
  }
  else{
    cout<<"Please select Katz and Ricotti 2013 or Madau and Fragos 2017"<<endl;
    exit(0);
  }  
  
  
  if(zred > 1100){  
    cout<<"ATTENZIONE: ERRORE REDSHIFT - "<<zred<<" "<<endl;
    exit(0);
  }
  
  
  return zred;
}


double Functions::maxwell(double a){

  double mum = rnd();

  double pi = M_PI;
  double x = 0.;
  double cum = erf(x/sqrt(2.*a*a)) - sqrt(2./pi) * (x*exp(-x*x/(2.*a*a))) / a;

  double x1 = -1;
  double x2 = -1;
  double c1 = -1;
  double c2 = -1;
  do{

    if(cum < mum){
      x1 = x;
      c1 = cum;
    }
    x += 0.05;
    cum = erf(x/sqrt(2.*a*a)) - sqrt(2./pi) * (x*exp(-x*x/(2.*a*a))) / a;
    if(cum > mum){      
      x2 = x;
      c2 = cum;
      break;
    }
    
  }while(cum < mum);


  double x_f = (x1 - x2) / (c1 - c2) * (mum-c2) + x2;

  return x_f;
  
  
}

double Functions::rndG(){
  std::random_device rd;
  std::mt19937 mt(rd());
  std::normal_distribution<double> distG(0.5,0.2);

  double rndG = distG(mt);
  do{
    rndG = distG(mt);
  }while(rndG < 0. || rndG > 1.);
  return rndG;
}

double Functions::rndG7(){
  std::random_device rd;
  std::mt19937 mt(rd());
  std::normal_distribution<double> distG(0.7,0.2);

  double rndG = distG(mt);
  do{
    rndG = distG(mt);
  }while(rndG < 0. || rndG > 1.);
  return rndG;
}

double Functions::rndgen(double pp, double spp){
  std::random_device rd;
  std::mt19937 mt(rd());
  std::normal_distribution<double> distG(pp,spp);
  
  double rndG = distG(mt);

  do{
    rndG = distG(mt);
  }while(rndG < 0. || rndG > 1.0);
  
  return rndG;
}



double Functions::Gaussian(double pp, double spp){
  std::random_device rd;
  std::mt19937 mt(rd());
  std::normal_distribution<double> distG(pp,spp);
  
  double rndG = distG(mt);

  do{
    rndG = distG(mt);
  }while(rndG < 0.);
  
  return rndG;
}

double Functions::LogGaussian(double pp, double spp){
  std::random_device rd;
  std::mt19937 mt(rd());
  std::lognormal_distribution<double> distG(pp,spp);
  
  double rndG = distG(mt);

  do{
    rndG = distG(mt);
  }while(rndG < 0.);
  
  return rndG;
}


double Functions::inter(double x, double *X, double *Y, int N){
  double y;
  int id1=0;
  int id2=N;
  if(X[0] > X[1]){	     
    for (int i=0;i<N;i++){
      if(x < X[i])
	id1 = i;
      if(x >= X[i]){
	id2 = i;
	break;
      }    
    }
  }
  else{
    for (int i=0;i<N;i++){
      if(x >= X[i])
	id1 = i;
      if(x < X[i]){
	id2 = i;
	break;
      }    
    }

  }
  
  if(id2 == 0){
    id1 = 0;
    id2 = 1;
  }
  if(id1 == N){
    id1 = N-2;
    id2 = N-1;
  }
     
  y = (Y[id1]-Y[id2])/(X[id1]-X[id2])*(x-X[id1]) + Y[id1];
  
  
  return y;
}

double Functions::rnd(){
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(0,1);

  double rnd = dist(mt);
  //double rnd = rand()/(double)RAND_MAX;
  
  return rnd;
}
double Functions::rnd_old(){
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(0,1);

  //double rnd = dist(mt);
  double rnd = rand()/(double)RAND_MAX;
  
  return rnd;
}

double Functions::tred(double z){
  //// redshift - time relation:

  double tH = 27.0164323345212; // Gyr
  double alpha = 1.82262661377844;
  double scatter =  0.3; //Gyr;

  double toff = tH / (1.+pow(1.+z, alpha)) ; //+ (-scatter + 2.*scatter*rnd());
    
  return toff;
}

double Functions::zred(double toff){
  //// redshift - time relation:
  /// only if t > 1.35e10
   
  double tH = 27.7263309054165; // Gyr
  double alpha = 2.886178235975;
  double scatter =  0.1; //Gyr;
  double zr = (pow(tH/toff - 1,1./alpha)-1.); // + (-scatter + 2.*scatter*rnd());

  if(zr < 0.0 && toff < 13.803)
    cout<<"Attention, time = "<<toff<<" "<<zr<<endl;
  
  return zr;
}

double Functions::vphiP(double x,double M,double a){
  return sqrt(2.*M/sqrt(x*x+a*a)*6.67E-11*1.99E30/3.08E16)*1.E-3;
}
double Functions::vphiD(double x,double M,double a,double g){
  return sqrt(2.*M/((2-g)*a)*(1.-pow(x/(x+a),2.-g))*6.67E-11*1.99E30/3.08E16)*1.E-3;
}

double Functions::Rh(double a,double g){
  return a/(pow(2.,1./(3.-g))-1);
}

double Functions::ah(double Rh,double g){
  return Rh*(pow(2.,1./(3.-g))-1);
}
double Functions::vescape(double g, double r, double m, string type){
  double vesc;


  double rn = rnd();

  vesc = vphiP(0.0,m,r/1.3);
     
  if(g > 0.){
    double rcore = 0.2;
    double rndmu = rnd();
    double bD    = pow(rndmu,1./(3.-g));
    double pos = bD/(1.-bD);        
    double rpos  = min(rcore, pos);
    double a     = ah(r,g);    

    vesc = vphiD(rpos, m, a, g);
  }  

  double rho = 3.*m/(8.*M_PI*r*r*r);
  vesc = 40. * pow(m/1.E5,1./3.) * pow(rho/1.E5, 1./6.);

  
  return vesc;

}


double Functions::mhalf_evo(double t, double rh, double mh, double trel, string type, string tclus){

  double tcc = 0.138*mh/(150. * log(0.11*mh/150.))*sqrt(pow(rh*3.08E16,3.)/(6.67E-11*1.99E30*mh))/(365.*24.*3600.);        
  if(tcc < 0)
    tcc = 0.2 * 4.2E9 * pow(rh/4.0,1.5) * sqrt(mh/1.E7) ;

  
  double trlx = trel ; //0.78E9 / log(0.11*mh) * pow(mh/1.E5,0.5) * pow(rh,1.5);
    
  double fplu = 1.0; 
  
  double mfact = 1.0;
  double mfact_h = 1.0;
    

  double fcohn = 300.0;

  double tstev = 1.E7;


  if(type == "GG23"){

  
    /*if(type == "under" || type == "over" || type == "critical"){
    if(type == "under")
      fcohn = 300.0;
    else if(type == "over" && tclus != "young")
      fcohn = 60.0;
    else if(type == "critical" || (type == "over" && tclus == "young"))
      fcohn = 5.0;
  }
  else if(type == "noevo"){
    mfact = 1.0;
  }
  else if(type == "mix"){
    if(tclus == "young")
      fcohn = 5.0;
    else if (tclus == "globular")
      fcohn = 300.0;
    else if (tclus == "nuclear")
      fcohn = 1000.0;    
  }
  else if(type == "GG23"){*/

    
    //Note: in Gieles&Gnedin23 they derive the lifetime of GCs and compare with observations. I used the lifetime to adjust the parameters below, keeping for GCs and NCs calculations related to a distance >= 8 kpc, whilst for YCs we adopt a distance ~ 1 kpc.
    
    if(tclus == "young")
      fcohn = 100. * pow(mh/1.E6,0.6);
    else
      fcohn = 1000. * pow(mh/1.E6,0.6);

  }
  else{
    cout<<"Select GG23, or implement your choice"<<endl;
    exit(0);
  }


  mfact_h = pow(1. + t/tstev,-0.1) * (1.-(t/(fcohn*trlx)));

  return mfact_h;
  
  
}

double Functions::rhalf_evo(double t, double rh, double mh, double trel, string type, string tclus){

  double tcc = 0.138*mh/(150. * log(0.11*mh/150.))*sqrt(pow(rh*3.08E16,3.)/(6.67E-11*1.99E30*mh))/(365.*24.*3600.);
  if(tcc < 0.0)
    tcc = 0.2 * 4.2E9 * pow(rh/4.0,1.5) * sqrt(mh/1.E7) ;

  double rfact, alpha, t0, rfact_h;
  double rmin = pow(5000. / mh , 1./0.43);

  if(type == "GG23"){
    t0 = tcc * 5. * (rh + 1.);
    alpha = 0.15 + 0.1 * rnd();
    rfact_h = pow(1. + t/t0,alpha);    
  }
  else{
     cout<<"Select GG23, or implement your choice"<<endl;
     exit(0);
  }
  return rfact_h;
  
}

double Functions::vevol(double t, double rh, double mh, double trel, string type, string tclus){
  double r = rhalf_evo(t, rh, mh, trel, type, tclus);
  double m = mhalf_evo(t, rh, mh, trel, type, tclus);
  double ve = sqrt(m/r);
  if(r <= 0.0 || m <= 0.0)
    ve = 0.0;

  return ve;
}


void Functions::histo(double* XX, int N, int nbin, string binning, string name){

  
  double *X;
  X = new double [N];
  
  for(int i=0;i<N;i++){
    X[i]=XX[i];
    if(X[i]<1.E-10)
      X[i] = 0.0;
  }
  
  double xmin =  1.E30;
  double xmax = -1.E30;
  if(binning=="log"){
    for(int i=0;i<N;i++){
      if(X[i]<xmin && X[i]>0.0) xmin = X[i];
      if(X[i]>xmax) xmax = X[i];
    }
  }
  else{
    for(int i=0;i<N;i++){      
      if(X[i]<xmin && X[i]>0.0) xmin = X[i];
      if(X[i]>xmax) xmax = X[i];
    }
  }
  
  double dx ;
  if(binning == "log"){
    xmax = log10(xmax);
    xmin = log10(xmin);
  }

  //cout<<xmin<<endl;
  
  dx = (xmax - xmin)/(double) nbin;
  
  double x0;
  double x1;
  
  int Ncnt = 0;
  int Ncum = 0;
  ofstream out;
  out.open(name.c_str());
  for(int k=0;k<nbin;k++){
    x0 = xmin + k*dx;
    x1 = xmin + (k+1)*dx;
    for(int i=0;i<N;i++){
      if(binning=="log"){
	if(X[i]<pow(10.,x1) && X[i]>=pow(10.,x0))
	  Ncnt += 1;	
      } 
      else{
	if(X[i]<x1 && X[i]>=x0)
	  Ncnt += 1;
      }
    }

    Ncum += Ncnt;
    
    if(binning=="log")
      out<<0.5*(pow(10.,x1)+pow(10.,x0))<<" "<<Ncnt<<" "<<Ncum<<" "<<N<<endl;
    else
      out<<0.5*(x1+x0)<<" "<<Ncnt<<" "<<Ncum<<" "<<N<<endl;


    Ncnt = 0;
  }
  out.close();
  
  delete [] X;
  
}


double Functions::findbhmin(double metal, string singpath){
  
  double mmaxBH = -1.;
  double mminBH = 1.E30;       
  double mmax   = -1.;
  double mmin   = 1.E30;
  double P;
  
  double mpri;
  string metastr;
  if     (metal == 0.0002)
    metastr="0002";
  else if(metal == 0.0004)
    metastr="0004";    
  else if(metal == 0.0008)
    metastr="0008";
  else if(metal == 0.0012)
    metastr="0012";
  else if(metal == 0.0016)
    metastr="0016";
  else if(metal == 0.002)
    metastr="002";
  else if(metal == 0.004)
    metastr="004";
  else if(metal == 0.006)
    metastr="006";
  else if(metal == 0.008)
    metastr="008";
  else if(metal == 0.012)
    metastr="012";
  else if(metal == 0.016)
    metastr="016";
  else if(metal == 0.02)
    metastr="02";
  else{
    cout<<"Metallicity for dynamical binaries out of range ..."<<endl;
    exit(0);      
  }
  
  string fname = singpath + "spectrum" + metastr + ".txt";
  ifstream isin;

  //cout<<"Our BH has metallicity "<<fname<<endl;
  
  isin.open(fname.c_str());
  if(!isin.is_open()){
    cout<<"File "<<fname<<" not found"<<endl;
    exit(0);
  }
  string line;  
  getline(isin,line);
  int lcnt = 0;
  double mzam, mfin, ktype;
  do{
    isin>>mzam>>mfin>>ktype;
    if(mfin>mmaxBH && ktype == 14)mmaxBH = mfin;
    if(mfin<mminBH && ktype == 14)mminBH = mfin;
    if(mzam>mmax   && ktype == 14)mmax   = mzam;
    if(mzam<mmin   && ktype == 14)mmin   = mzam;
    lcnt+=1;
  }while(!isin.eof());
  isin.close();
  lcnt -= 1;


  return mminBH;


}


double Functions::singBHt(double metal, string singpath, double obslope, double mslope){

  
  double mmaxBH = -1.;
  double mminBH = 1.E30;       
  double mmax   = -1.;
  double mmin   = 1.E30;
  double P;
  
  double mpri;
  string metastr;
  
  if     (metal == 0.0002)
    metastr="0002";
  else if(metal == 0.0004)
    metastr="0004";    
  else if(metal == 0.0008)
    metastr="0008";
  else if(metal == 0.0012)
    metastr="0012";
  else if(metal == 0.0016)
    metastr="0016";
  else if(metal == 0.002)
    metastr="002";
  else if(metal == 0.004)
    metastr="004";
  else if(metal == 0.006)
    metastr="006";
  else if(metal == 0.008)
    metastr="008";
  else if(metal == 0.012)
    metastr="012";
  else if(metal == 0.016)
    metastr="016";
  else if(metal == 0.02)
    metastr="02";
  else{
    cout<<"Metallicity for dynamical binaries out of range ..."<<endl;
    exit(0);      
  }
  
  string fname = singpath + "spectrum" + metastr + ".txt";
  ifstream isin;

  //cout<<"Our BH has metallicity "<<fname<<endl;
  
  isin.open(fname.c_str());
  if(!isin.is_open()){
    cout<<"File "<<fname<<" not found"<<endl;
    exit(0);
  }
  string line;
  getline(isin,line);
  int lcnt = 0;
  double mzam, mfin, ktype;
  do{
    isin>>mzam>>mfin>>ktype;
    //cout<<mzam<<" "<<mfin<<" "<<ktype<<endl;
    if(mfin>mmaxBH && ktype == 14)mmaxBH = mfin;
    if(mfin<mminBH && ktype == 14)mminBH = mfin;
    if(mzam>mmax   && ktype == 14)mmax   = mzam;
    if(mzam<mmin   && ktype == 14)mmin   = mzam;
    lcnt+=1;
  }while(!isin.eof());
  isin.close();
  lcnt -= 1;

  
  double mstar,mobs;
  double mpre = -1; double mbhpre = -1;
  double msuc = -1; double mbhsuc = -1;
  double mzam0,mfin0;
 
  do{
    mpre   = -1;
    mbhpre = -1;
    msuc   = -1;
    mbhsuc = -1;
    P = rnd();   
    mstar = pow( pow(mmax,1.+mslope)*P + (1.-P)*pow(mmin,1.+mslope)   ,1./(1.+mslope));
    isin.open(fname.c_str());
    getline(isin,line);
    for(int k=0;k<lcnt;k++){
      isin>>mzam>>mfin>>ktype;
      if(k==0){
	mzam0 = mzam;
	mfin0 = mfin;
      }
      if     (mstar > mzam && ktype==14) {
	mbhpre = mfin;
	mpre   = mzam;
      }
      else if(mstar < mzam && ktype==14){
	mbhsuc = mfin;
	msuc   = mzam;
	break;
      }
    }
    isin.close();

    if(mbhsuc == -1){
      mbhsuc = mfin0;
      msuc   = mzam0;
    }
    if(mbhsuc == -1){
      mbhsuc = mfin;
      msuc   = mzam;
    }
    
    mobs = (mstar - mpre)*(mbhsuc-mbhpre)/(msuc-mpre) + mbhpre;


    double acc_rej = ( pow(mobs,1.+obslope)-pow(mminBH,1.+obslope) )/( pow(mmaxBH,1.+obslope)-pow(mminBH,1.+obslope) )  ;
    P = rnd();

    if(P<acc_rej) break;

  }while(1>0);    

  mpri = mobs;
  
  return mpri;

}



void Functions::singBHt_mix_old(vector<double> zams_mix, vector<double> remn_mix, vector<double> tdel_mix, vector<double> kick_mix, double obslope, double mslope, double *sing_out, double msmax, double msmin, double mbmax, double mbmin, double vescape){


  double mpri, tpri, kick;
  double mmax = msmax;
  double mmin = msmin;

  double P = rnd();   
  double mstar = pow( pow(mmax,1.+mslope)*P + (1.-P)*pow(mmin,1.+mslope)   ,1./(1.+mslope));

  vector<int> idix;
  
  for(int i=0;i<zams_mix.size();i++){
    if(abs(mstar - zams_mix[i])/zams_mix[i] < 0.1 && remn_mix[i] > 3.0)
      idix.push_back(i);
  }
    
  int ind = (int) (idix.size()*rnd());
  
  mpri = remn_mix[ind]; 
  tpri = tdel_mix[ind];
  kick = kick_mix[ind];
  
  sing_out[0] = mpri;
  sing_out[1] = tpri;
  sing_out[2] = kick;

  idix.erase(idix.begin(),idix.end());
  
  return;
  
}




void Functions::singBHt_new(vector<double>& zams_sin,
		 vector<double>& remn_sin,
		 vector<double>& tdel_sin,
		 vector<double>& kick_sin,
		 double *sing_out,
		 double vescape){

  
  double mblack, vblack, tblack;
  sing_out[0] = 0.0;
  sing_out[1] = 0.0;
  sing_out[2] = 0.0;

  int cat_size = zams_sin.size() ;
  int id;
  int nctn = 0;
  do{
    id = static_cast<int>(cat_size * rnd());
    nctn++;    
  }while(kick_sin[id] > vescape);
  

  
  sing_out[0] = remn_sin[id];
  sing_out[1] = tdel_sin[id];
  sing_out[2] = kick_sin[id];
  
  
  return ;
}


/*void Functions::singBHt_new_no(vector<double>& zams_sin, vector<double>& remn_sin, vector<double>& tdel_sin, vector<double>& kick_sin, double obslope, double mslope, double *sing_out, double msmax, double msmin, double mbmax, double mbmin, double vescape){

  double mpri;

  double mstar,mobs;
  double mpre = -1; double mbhpre = -1; double tpre = -1; double kpre = -1;
  double msuc = -1; double mbhsuc = -1; double tsuc = -1; double ksuc = -1;
  double mzam0,mfin0,t0;

  double time_for, kick_for;
  double P;

  double mmax = msmax;
  double mmin = msmin;
  double mmaxBH = mbmax;
  double mminBH = mbmin;
  
  int nsafe = 0;
  kick_for = 1.E30;
  
  vector<double> goodB, goodM, goodV, goodT;

  if(kick_sin.size() > 1)
    for(int i=0 ; i < kick_sin.size(); i++){
      if(kick_sin[i] < vescape){
	goodM.push_back(zams_sin[i]);
	goodB.push_back(remn_sin[i]);
	goodV.push_back(kick_sin[i]);
	goodT.push_back(tdel_sin[i]);
      }
    }
    
  int flag = 0;
  sing_out[3] = 0.0;
  if(goodB.size() == 0){
    //    cout<<"No BHs form here "<<vescape<<endl;
    sing_out[0] = 0.0;
    sing_out[1] = 0.0;
    sing_out[2] = 1.E3;
    sing_out[3] = 1.0;
    flag = 1;    
  }
  
  if(flag == 0){
    P = rnd();   
    mstar = pow( pow(mmax,1.+mslope)*P + (1.-P)*pow(mmin,1.+mslope)   ,1./(1.+mslope));
    int id1 = -1;
    int id2 = -1;
    for(int i=0;i<goodM.size();i++){
      if(mstar > goodM[i])
	id1 = i;
      if(mstar < goodM[i]){
	id2 = i;
	break;
      }            
    }
    
    msuc = goodM[id2];
    mpre = goodM[id1];
    mbhsuc=goodB[id2];
    mbhpre=goodB[id1];
    
    tsuc = goodT[id2];
    tpre = goodT[id1];
    
    ksuc = goodV[id2];
    kpre = goodV[id1];
    
    mobs = (mstar - msuc)*(mbhsuc-mbhpre)/(msuc-mpre) + mbhsuc;
    time_for = (tpre - tsuc)/(mpre-msuc) * (mstar - mpre) + tpre;
    kick_for = (kpre - ksuc)/(mpre-msuc) * (mstar - mpre) + kpre;
    
    
    sing_out[0] = mobs;
    sing_out[1] = time_for;
    sing_out[2] = kick_for;
    sing_out[3] = flag;

    
    //cout<<mstar<<" "<<mobs<<" "<<time_for<<" "<<kick_for<<" "<<goodM[id1]<<" "<<goodM[id2]<<" "<<goodB[id1]<<" "<<goodB[id2]<<endl;
    
    
    goodM.erase(goodM.begin(),goodM.end());
    goodB.erase(goodB.begin(),goodB.end());
    goodT.erase(goodT.begin(),goodT.end());
    goodV.erase(goodV.begin(),goodV.end());
  }
  
  return;

  }*/


double Functions::mratio(double mpri, double MRATIO_SLOPE, string MRATIO_STR){
  double q;

  if(MRATIO_STR=="uniform")
    q = rnd();
  else if(MRATIO_STR=="powerlaw" || MRATIO_STR=="pwl"){
    double qrnd = rnd();
    double ratslp = MRATIO_SLOPE;
    q = pow(qrnd,1./(1.-ratslp));
  }
  else if(MRATIO_STR=="gaussian" || MRATIO_STR=="gss"){
    q = rndG7();
    do{
      q = rndG7();
    }while(q<0. || q>1.);
  }
  else{     
    cout<<"mratio distribution not recognized -- retry "<<endl;
    exit(0);
  }   	



  return q;
}






double Functions::sevntab(double M0, double metal, string path){

    //reading tables
    string tabname;
    if(metal==0.02){
      tabname = path+"SperaMapelli17_z002.txt";
    }
    else if(metal==0.002){
      tabname = path+"SperaMapelli17_z0002.txt";
    }	
    else if(metal==0.0002){
      tabname = path+"SperaMapelli17_z00002.txt";
    }       
    else{
      cout<<"Metallicity available only Z=2x10^{-4,-3,-2} "<<metal<<endl;
      exit(0);
    } 

    //cout<<"Metallicity investigated "<<metal<<" tables "<<tabname<<endl;

    ifstream sprin;
    int clin = 0;
    string sline;
    sprin.open(tabname.c_str());
    if(!sprin.is_open()){
      cout<<"File "<<tabname<<" not found"<<endl;
      exit(0);
    }
    getline(sprin,sline);
    double aa;
    do{
      getline(sprin,sline);
      clin+=1;
    }while(!sprin.eof());
    sprin.close();
    //RE-READING AND SAVE THINGS OR STOP WHEN YOU FIND WHAT YOU'RE LOOKING FOR!!
    sprin.open(tabname.c_str());
    getline(sprin,sline);
    double preMI, postMI, preMF, postMF;
    preMI = -1;
    postMI = -1;
    preMF = -1;
    postMF = -1;
    double msevn0[clin-1],msevnF[clin-1];
    for(int ll=0;ll<clin-1;ll++){
      sprin>>msevn0[ll]>>aa>>aa>>aa>>aa>>msevnF[ll];	    	
    }      	
    sprin.close();
    
    for(int ll=0;ll<clin-1;ll++){
      if(msevn0[ll]<M0){
	preMI = msevn0[ll];
	preMF = msevnF[ll];
      }
      else if(msevn0[ll]>M0){
	postMI = msevn0[ll];
	postMF = msevnF[ll];
	break;
      }
      
    }

    //cout<<preMI<<" "<<preMF<<" "<<M0<<endl;
    double realmass = (preMF-postMF)/(preMI-postMI)*(M0-preMI) + preMF;
    return realmass;

}


int Functions::bbh_UIBfits_setup(double mbh1, double mbh2, double chibh1, double chibh2, double *out){

  double m1 = mbh1;
  double m2 = mbh2;
  double chi1=chibh1;
  double chi2=chibh2;

  if(m2>m1){m2 = mbh1; chi2 = chibh1; m1 = mbh2; chi1 = chibh2;}

  double m    = m1+m2;

  double msq  = m*m;
  double m1sq = m1*m1;
  double m2sq = m2*m2;

  double eta  = m1*m2/msq;
  if(eta>0.25)eta = min(eta,0.25);
  if(eta<0.0) eta = 0.0;

  double eta2 = eta*eta;
  double eta3 = eta2*eta;
  double eta4 = eta2*eta2;

  // spin variables (in m = 1 units)
  double S1    = chi1*m1sq/msq; //# spin angular momentum 1
  double S2    = chi2*m2sq/msq; //# spin angular momentum 2
  double Stot  = S1+S2 ;        //# total spin
  double Shat  = (chi1*m1sq+chi2*m2sq)/(m1sq+m2sq); //# effective spin, = msq*Stot/(m1sq+m2sq)
  double Shat2 = Shat*Shat;
  double Shat3 = Shat2*Shat;
  double Shat4 = Shat2*Shat2;

  double chidiff  = chi1 - chi2;
  double chidiff2 = chidiff*chidiff;

  //    # typical squareroots and functions of eta
  double sqrt2 = pow(2.,0.5);
  double sqrt3 = pow(3.,0.5);
  double sqrt1m4eta = pow(1. - 4.*eta,0.5);

  /*  out.push_back(m);
  out.push_back(eta);
  out.push_back(eta2);
  out.push_back(eta3);
  out.push_back(eta4);
  out.push_back(Stot);
  out.push_back(Shat);
  out.push_back(Shat2);
  out.push_back(Shat3);
  out.push_back(Shat4);
  out.push_back(chidiff);
  out.push_back(chidiff2);
  out.push_back(sqrt2);
  out.push_back(sqrt3);
  out.push_back(sqrt1m4eta);*/

  out[0] = m;
  out[1] = eta;
  out[2] = eta2;
  out[3] = eta3;
  out[4] = eta4;
  out[5] = Stot;
  out[6] = Shat;
  out[7] = Shat2;
  out[8] = Shat3;
  out[9] = Shat4;
  out[10] = chidiff;
  out[11] = chidiff2;
  out[12] = sqrt2;
  out[13] = sqrt3;
  out[14] = sqrt1m4eta;

  return 0;


}

double Functions::bbh_final_mass_non_precessing_UIB2016(double m1, double m2, double chi1, double chi2, string version){
  
  double mfin;

  double out[15];
  bbh_UIBfits_setup(m1, m2, chi1, chi2, out);

  double m = out[0];
  double eta= out[1];
  double eta2= out[2];
  double eta3= out[3];
  double eta4= out[4];
  double Stot = out[5];
  double Shat= out[6];
  double Shat2 = out[7];
  double Shat3= out[8];
  double Shat4= out[9];
  double chidiff= out[10]; 
  double chidiff2 = out[11]; 
  double sqrt2= out[12];
  double sqrt3 = out[13];
  double sqrt1m4eta = out[14];


  if(version!="v2"){
    cout<<"choose version v2"<<endl;exit(0);
  }

  //# rational-function Pade coefficients (exact) from Eq. (22) of 1611.00332v2
  double b10 = 0.346;
  double b20 = 0.211;
  double b30 = 0.128;
  double b50 = -0.212;
//# fit coefficients from Tables VII-X of 1611.00332v2
//# values at increased numerical precision copied from
//# https://git.ligo.org/uib-papers/finalstate2016/blob/master/LALInference/EradUIB2016v2_pyform_coeffs.txt
//# git commit f490774d3593adff5bb09ae26b7efc6deab76a42
   double     a2 = 0.5609904135313374;
   double     a3 = -0.84667563764404;
   double    a4 = 3.145145224278187;
   double     b1 = -0.2091189048177395;
   double     b2 = -0.19709136361080587;
   double     b3 = -0.1588185739358418;
   double     b5 = 2.9852925538232014;
   double     f20 = 4.271313308472851;
   double    f30 = 31.08987570280556;
   double    f50 = 1.5673498395263061;
   double    f10 = 1.8083565298668276;
   double     f21 = 0.;
   double     d10 = -0.09803730445895877;
   double     d11 = -3.2283713377939134;
   double     d20 = 0.01118530335431078;
   double     d30 = -0.01978238971523653;
   double     d31 = -4.91667749015812;
   double    f11 = 15.738082204419655;
   double     f31 = -243.6299258830685;
   double     f51 = -0.5808669012986468;


   double Erad = (((1. + -2.0/3.0*sqrt2)*eta + a2*eta2 + a3*eta3 + a4*eta4)*(1. + b10*b1*Shat*(f10 + f11*eta + (16. - 16.*f10 - 4.*f11)*eta2) + b20*b2*Shat2*(f20 + f21*eta + (16. - 16.*f20 - 4.*f21)*eta2) + b30*b3*Shat3*(f30 + f31*eta + (16. - 16.*f30 - 4.*f31)*eta2)))/(1. + b50*b5*Shat*(f50 + f51*eta + (16. - 16.*f50 - 4.*f51)*eta2)) + d10*sqrt1m4eta*eta2*(1. + d11*eta)*chidiff + d30*Shat*sqrt1m4eta*eta*(1. + d31*eta)*chidiff + d20*eta3*chidiff2;

     //    # Convert to actual final mass
     double Mf = m*(1.-Erad);

     mfin = Mf;



  return mfin;

}

double Functions::bbh_final_spin_non_precessing_UIB2016(double m1, double m2, double chi1, double chi2, string version){
  
  double out[15];
  bbh_UIBfits_setup(m1, m2, chi1, chi2, out);

  double m = out[0];
  double eta= out[1];
  double eta2= out[2];
  double eta3= out[3];
  double eta4= out[4];
  double Stot = out[5];
  double Shat= out[6];
  double Shat2 = out[7];
  double Shat3= out[8];
  double Shat4= out[9];
  double chidiff= out[10]; 
  double chidiff2 = out[11]; 
  double sqrt2= out[12];
  double sqrt3 = out[13];
  double sqrt1m4eta = out[14];


  double sfin;
  if(version!="v2"){
    cout<<"choose version v2"<<endl;exit(0);
  }

  //# rational-function Pade coefficients (exact) from Eqs. (7) and (8) of 1611.00332v2
  double      a20 = 5.24;
  double      a30 = 1.3;
  double      a50 = 2.88;
  double      b10 = -0.194;
  double      b20 = 0.0851;
  double      b30 = 0.00954;
  double      b50 = -0.579;
  //# fit coefficients from Tables I-IV of 1611.00332v2
  //# values at increased numerical precision copied from
  //# https://git.ligo.org/uib-papers/finalstate2016/blob/master/LALInference/FinalSpinUIB2016v2_pyform_coeffs.txt
  //# git commit f490774d3593adff5bb09ae26b7efc6deab76a42
  double      a2 = 3.8326341618708577;
  double      a3 = -9.487364155598392;
  double      a5 = 2.5134875145648374;
  double      b1 = 1.0009563702914628;
  double      b2 = 0.7877509372255369;
  double      b3 = 0.6540138407185817;
  double      b5 = 0.8396665722805308;
   double     f21 = 8.77367320110712;
  double      f31 = 22.830033250479833;
  double      f50 = 1.8804718791591157;
  double      f11 = 4.409160174224525;
  double      f52 = 0.;
  double      d10 = 0.3223660562764661;
  double      d11 = 9.332575956437443;
  double      d20 = -0.059808322561702126;
  double      d30 = 2.3170397514509933;
  double      d31 = -3.2624649875884852;
  double      f12 = 0.5118334706832706;
  double      f22 = -32.060648277652994;
  double      f32 = -153.83722669033995;
  double      f51 = -4.770246856212403;

  double  Lorb = (2.*sqrt3*eta + a20*a2*eta2 + a30*a3*eta3)/(1. + a50*a5*eta) + (b10*b1*Shat*(f11*eta + f12*eta2 + (64. - 16.*f11 - 4.*f12)*eta3) + b20*b2*Shat2*(f21*eta + f22*eta2 + (64. - 16.*f21 - 4.*f22)*eta3) + b30*b3*Shat3*(f31*eta + f32*eta2 + (64. - 16.*f31 - 4.*f32)*eta3))/(1. + b50*b5*Shat*(f50 + f51*eta + f52*eta2 + (64. - 64.*f50 - 16.*f51 - 4.*f52)*eta3)) + d10*sqrt1m4eta*eta2*(1. + d11*eta)*chidiff + d30*Shat*sqrt1m4eta*eta3*(1. + d31*eta)*chidiff + d20*eta3*chidiff2;

    //    # Convert to actual final spin
  double chif = Lorb + Stot;

  sfin = chif;

  return sfin;

}
double Functions::angmom(double a, double atot, double v, int Nmet){
  

  double Z1 = 1.+pow(1.-a*a, 1./3.)*(pow(1.+a, 1./3.)+pow(1.-a, 1./3.))  ;
  double Z2 = sqrt(3.*a*a+Z1*Z1);
  
  double sgn = 1.;
  if(a<0.)sgn = -1.;
  double  risco = 3.+Z2 - sgn*sqrt((3.-Z1)*(3.+Z1+2.*Z2))  ;
  
  double  Lisco = 2./(3.*sqrt(3.))*(1.+2.*sqrt(3.*risco-2.));
  
  double  Eisco = sqrt(1.-2./(3.*risco));
  
  double k[Nmet][Nmet];
  if(Nmet == 4){
  
    k[0][0] =  -5.9;
    k[0][1] =  2.87025;
    k[0][2] = -1.53315;
    k[0][3] = -3.78893;
    
    k[1][0] =  32.9127 ;
    k[1][1] = -62.9901;
    k[1][2] =  10.0068;
    k[1][3] =  56.1926;
    
    k[2][0] = -136.832;
    k[2][1] =  329.32 ;
    k[2][2] = -13.2034 ;
    k[2][3] = -252.27 ;
    
    k[3][0] =  210.075; 
    k[3][1] = -545.35 ;
    k[3][2] = -3.97509 ;
    k[3][3] =  368.405 ;
  }
  else{
    k[0][0] = -5.9;
    k[0][1] = 3.39221;
    k[0][2] = 4.48865;
    k[0][3] = -5.77101;
    k[0][4] = -13.0459;
    
    k[1][0] = 35.1278;
    k[1][1] = -72.9336;
    k[1][2] = -86.0036;
    k[1][3] = 93.7371;
    k[1][4] = 200.975;
    
    k[2][0] = -146.822;
    k[2][1] = 387.184;
    k[2][2] = 447.009;
    k[2][3] = -467.383;
    k[2][4] = -884.339;
    
    k[3][0] = 223.911;
    k[3][1] = -648.502;
    k[3][2] = -697.177;
    k[3][3] =  753.738; 
    k[3][4] = 1166.89;  
  }

  double sum = 0.0;

  for(int i=0;i<Nmet;i++){
    for(int j=0;j<Nmet;j++){
      sum += k[i][j]*pow(v,i+1)*pow(a,j);
    }
  }
  
 

  double l = abs(Lisco - 2.*atot*(Eisco-1.) + sum);
  
  //cout<<"Angular mom "<<Lisco<<" "<<atot<<" "<<Eisco<<" "<<sum<<endl;
  
  
  return l;
  
}

double Functions::angmom2(double a1, double a2, double v, double q, double cosa, double cosb, double cosg){

  double l;

  double s4 = -0.1229;
  double s5 =  0.4537;
  double t0 = -2.8904;
  double t3 =  2.5763;
  double t2 = -3.5171;

  l = 2.*sqrt(3.)+    t2*v+    t3*v*v+    s4/((1.+q*q)*(1.+q*q))*(a1*a1+a2*a2*q*q*q*q+2.*a1*a2*q*q*cosa)+ (s5*v+t0+2.)/(1.+q*q)*(a1*cosb+a2*q*q*cosg);

  return l;
}



void Functions::SREM2(double ndx, double a1, double a2, double m1, double m2, string align, double *spins){

  double originalm1 = m1;
  double originalm2 = m2;
    
  double afin = 2.0;
  double afin1,afin2,afin3;
  double CHI;

  if(a1>1.09 || a2>1.09 || a1 <0.0 || a2 < 0.0){
    cout<<"Critical spin error "<<a1<<" "<<a2<<endl;
    if(a1>1.0)a1=1.0;
    if(a2>1.0)a2=1.0;
  }
  if(m1<0.0 || m2<0.0){
    cout<<"Critical mass error "<<endl;exit(0);
  }
 
  
  double q   = m1/m2;
  if(m1>m2) q = m2/m1;
  if(m1==0.0) q = 0.0;
  if(m2==0.0) q = 0.0;


  if(m1 > m2){
    double useful = m1;
    m1 = m2;
    m2 = useful;

    useful = a1;
    a1 = a2;
    a2 = useful;

  }
  

  
  double v = q/((1.+q)*(1.+q));

  //cosalpha = angle between the spins;
  //cosbeta  = angle with the angular momentum;
  //cosgamma = angle with the angular momentum;
  

  double l;
  double cosa,cosb,cosg;
  double chkneg;

  do{
    double ccc  = rnd();
    double cosbeta = 2.*pow(ccc,1./(ndx+1.))-1.;    

    double ccc2 = rnd();
    double cosgamma = 2.*pow(ccc2,1./(ndx+1.))-1.;

    double cosalpha = -1.+2.*rnd();


    if(align=="align"){
      cosalpha = 1.0;
      cosbeta  = cosgamma;
    }
    else if(align=="antialign"){
      cosalpha = -1.0;
      cosbeta  = cosgamma;      
    }
    else if(align=="dynamical"){
      cosalpha = -1.+2.*rnd();
      cosbeta  = -1.+2.*rnd();
      cosgamma = -1.+2.*rnd();      
    }
    else if(align!="whatever"){
      cout<<"Please choose align, antialign, dynamical, or whatever and retry"<<endl;
      exit(0);
    }
    
    cosa = cosalpha;
    cosb = cosbeta;
    cosg = cosgamma;
    
    /*cosa = 1.0;
      cosb = 1.0;
      cosg = 1.0;*/
    
    
    if(a1>1. || a2>1.){cout<<"BH spins too large "<<a1<<" "<<a2<<endl;exit(0);}
    
    
    int Nmet = 4;
    double csi;
    if(Nmet == 5) csi= 0.474046; 
    else if(Nmet == 4) csi = 0.463926;
    else{
      cout<<"choose a method please "<<endl;exit(0);
    }
    double atot = (a1*cosb + a2*cosg*q*q)/((1+q)*(1+q));
    double aeff = atot + csi*v*(a1*cosb+a2*cosg);
    
    //Hofmann et al. 2016
    l = angmom(aeff, atot, v, Nmet); 
    afin1 = 1./((1.+q)*(1.+q))*sqrt(a1*a1 + a2*a2*q*q*q*q + 2.*a2*a1*q*q*cosa + 2.*(a1*cosb + a2*q*q*cosg)*l*q +l*l*q*q);
    
    //Barausse & Rezzolla 2009 if(afin>1.0){ 
    l = angmom2(a1, a2, v, q, cosa, cosb, cosg);    
    afin2 = 1./((1.+q)*(1.+q))*sqrt(a1*a1 + a2*a2*q*q*q*q + 2.*a2*a1*q*q*cosa + 2.*(a1*cosb + a2*q*q*cosg)*l*q +l*l*q*q);
    //}
     chkneg = a1*a1 + a2*a2*q*q*q*q + 2.*a2*a1*q*q*cosa + 2.*(a1*cosb + a2*q*q*cosg)*l*q +l*l*q*q ;
  }while(chkneg<0.);
  
  double a1par = a1*cosa*cosb;
  double a2par = a2*cosg;
  double Mfin = bbh_final_mass_non_precessing_UIB2016( m1,  m2,  a1,  a2, "v2");
  double Sfin = bbh_final_spin_non_precessing_UIB2016( m1,  m2,  a1par,  a2par, "v2");
  
  double a1per = a1*sin(acos(cosb))*sin(acos(cosa)); 
  double a2per = a2*sin(acos(cosg));
  double Sper = (m1*m1*a1per + m2*m2*a2per)/((m1+m2)*(m1+m2));
  
  afin3 = sqrt(Sfin*Sfin + Sper*Sper);

  if(a1<0. || a2 <0.){
    cout<<"PD2 = PDDD"<<endl;
    exit(0);
  }

  //  -- MOD APRIL 2021--> USE ONLY JIMENEZ
  /*if(afin1<=1.0 && afin1>0.0){
    afin = afin1;
  }
  else{*/
  
  //cout<<"Using Jimenez-Forteza 2018 "<<q<<" "<<a1<<" "<<a2<<" "<<afin1<<" "<<afin2<<" "<<afin3<<endl;
  afin = afin3;
  if(afin > 1)
    afin = afin1;
  
    //}

  //if(afin>1.0){cout<<"PD"<<endl;exit(0);}


  stringstream aaa;
  aaa<<afin;
  if(aaa.str()=="-nan" || aaa.str()=="nan"){
    cout<<"Got a nan"<<m1<<" "<<m2<<" "<<a1<<" "<<a2<<" "<<afin<<endl;
    exit(0);
  }



  //// RECOILING KICK VELOCITIES.... Ref: Campanelli+06, Lousto & Zlochower (2008), Antonini+2019, Gonzalez et al. 2007a (numerical values)
  double vkick = kicks(m1, m2, a1, a2, q, cosa, cosb, cosg); 

  double cos_primary   = cosg;
  double cos_secondary = cosb;

  if(originalm1 < originalm2){
    cos_primary = cosb;
    cos_secondary= cosg;
  }
  
  
  spins[0] = afin;
  spins[1] = (m1*a1*cosb + m2*a2*cosg)/(m1+m2);
  spins[2] = Mfin;
  spins[3] = vkick;

  spins[4] = cosa;
  spins[5] = cos_primary;
  spins[6] = cos_secondary;

  return;

}


double Functions::kicks(double m1, double m2, double a1, double a2, double q, double cosa, double cosb, double cosg){ 
  
  double sina = sin(acos(cosa));
  double sinb = sin(acos(cosb));
  double sing = sin(acos(cosg));

    
  double a2par = a2 * cosg;
  double a2per1= a2 * sing; 
  double a2per2= 0.0;
  
  double a1par = a1 * cosb;
  double a1per1= a1 * sinb*cosa;
  double a1per2= a1 * sinb*sina;
  
  double Deltapar  = (m1+m2)*(m1+m2) / (1+q) * (a2par  - q*a1par);
  double Deltaper1 = (m1+m2)*(m1+m2) / (1+q) * (a2per1 - q*a1per1);
  double Deltaper2 = (m1+m2)*(m1+m2) / (1+q) * (a2per2 - q*a1per2);
  
  double Deltaper  = sqrt(Deltaper1*Deltaper1 + Deltaper2*Deltaper2);
  
  double dir1 = rnd();
  double pm   = 1.-2.*rnd();
  double dir2 = abs(pm)/pm * sqrt(1.-dir1*dir1);
  double dir  = sqrt(dir1*dir1 + dir2*dir2);


  double SEpar = 2.*(a2par + q*q*a1par)/((1.+q)*(1.+q));

  
  double PHI = (Deltaper1*dir1 + Deltaper2*dir2)/(Deltaper * dir);
  double PHI1= 2.*M_PI*rnd();
  double A = 1.2E4; // km/s  
  double B = -0.93; // adim
  double H = 6.9E3; //km/s
  double XI= 145.*M_PI/180.; //degrees 
  double V11 = 3677.76; // km/s
  double VA  = 2481.21; // km/s
  double VB  = 1792.45; // km/s
  double VC  = 1506.52; // km/s

  double ETA= q/((1.+q)*(1.+q)); //asimmetric mass ratio
  
  
  double VM   = A*ETA*ETA*sqrt(1.-4.*ETA)*(1.+B*ETA);


  double VPER = H*ETA*ETA/(1.+q) * (a2par - q*a1par);

    
  double VTER = V11 + VA*SEpar + VB*SEpar*SEpar + VC*SEpar*SEpar*SEpar;
  double SPPE = sqrt((a2per1 - q*a1per1)*(a2per1 - q*a1per1) + (a2per2 - q*a1per2)*(a2per2 - q*a1per2))*cos(PHI-PHI1);
  double VPAR = 16.*ETA*ETA/(1.+q) * VTER * SPPE;
  
  double vkiper1 = VM + VPER*cos(XI);
  double vkiper2 = VPER*sin(XI);
  double vkipar  = VPAR;

  double vkick = sqrt(vkiper1*vkiper1 + vkiper2*vkiper2 + vkipar*vkipar);

  return vkick;
}
    
				      
				      
double Functions::spin(double mass, string spinning){  
  double A = 0.893;double B=8.0;double Mo=60.0;double C = 0.21;
  double spi;
  if(spinning=="chris")
    spi   = A/pow(1.+pow(mass/Mo,B),C);   
  else if(spinning=="random")
    spi = rnd();
  else if(spinning=="gaussian")
    spi = rndG();
  else if(spinning=="gaussian02")
    spi = rndgen(0.2,0.1);  
  else if(spinning=="maxwellian02"){
    spi = maxwell(0.2);
    do{
     spi = maxwell(0.2);
    }while(spi > 1.);
  }
  else if(spinning=="maxwellian01"){
    spi = maxwell(0.1);
    do{
     spi = maxwell(0.1);
    }while(spi > 1.);
  }
  else if(spinning=="fuller")
    spi = 0.001;
  else if(spinning=="no")
    spi = 1.E-6;
  else if(spinning=="zero")
    if(mass<30.)
      spi = rnd();
    else
      spi = 30./mass*(0.075 + 0.45*rnd());
  else{
    cout<<"Spinning either 'chris', 'random' or 'zero'"<<endl;
    exit(0);
  }
  return spi;


}




string Functions::print(double mass, double mmax, double mmin,double wgh){
  double P = pow(mass,wgh);  //(pow(mmax,1.+wgh) - pow(mass,1.+wgh))/(pow(mmax,1.+wgh) - pow(mmin,1.+wgh));
  double O = mmin + mmax*rnd();
  O = pow(O,wgh);
  string succ;
  if(O<=P) succ = "success";
  else    succ = "failed";
  
  succ = "success";
  return succ;
}


double Functions::GWeff(string pcluster, double met){

  //This function computes the efficiency of formation of merging BBHs in different environments
  //as a function of the metallicity of the progenitor stars and the type of cluster
  //The formula is a fit from the catalogs

  Functions func;

  double eps;
  double a,c;
  if(pcluster == "young"){
    a = 0.000134696 + 1.543E-5*(1.-2.*func.rnd());
    c = -4.25476 + 1.221 *(1.-2.*func.rnd());
  }
  else if(pcluster == "globular"){
    a = 0.00053 + 2.46E-5*(1.-2.*func.rnd());
    c = -3.1727 + 0.3818*(1.-2.*func.rnd());
  }
  else if(pcluster == "nuclear"){
    a = 0.0011 + 2.2E-5 * (1.-2.*func.rnd());
    c = -1.84 + 0.11 * (1.-2.*func.rnd());  
  }

  eps = a * pow(1.+met/0.02,c);
   
  return eps;
}



void Functions::singBHt_mix(vector<double>& zams_mix,
//This function selects a single BH from the mixed population
      vector<double>& remn_mix,
      vector<double>& tdel_mix,
      vector<double>& kick_mix,
      double *sing_out,
      double vescape){

        
    Functions func;
    double mblack, vblack, tblack;
    sing_out[0] = 0.0;
    sing_out[1] = 0.0;
    sing_out[2] = 0.0;

    int cat_size = zams_mix.size() ;
    int id;
    int nctn = 0;
    do{
      id = static_cast<int>(cat_size * func.rnd());
      nctn++;    
    }while(kick_mix[id] > vescape);
    
    //cout<<id<<" "<<zams_mix[id]<<" "<<remn_mix[id]<<" "<<tdel_mix[id]<<" "<<kick_mix[id]<<" "<<vescape<<endl;
    
    sing_out[0] = remn_mix[id];
    sing_out[1] = tdel_mix[id];
    sing_out[2] = kick_mix[id];
    
    //cout<<sing_out[0]<<" "<<sing_out[1]<<" "<<sing_out[2]<<" "<<vescape<<" "<<nctn<<endl;
    return ;
  }

  double Functions::inter_rate(double m1, double m2, double vesc, double m_hg, int gen,
            double nbhs, double mhalf, double mcore, double rcore, double n_bin,
            double trelax, double t12capt, double tbbhform, double tcc,
            string pcluster){

  //  This function computes the interaction rate considering:
  //  - the density of stellar BHs in the core
  //  - the density of hierarchical BHs in the core --> computed considering the influence radius and the wandering radius of the hierarchical BHs
  //  - the average mass of BHs
  //  - the mass of the binary
  //  - the mass of the hierarchical BH
  //  - the escape velocity of the cluster
  //  - the relaxation time of the cluster
  //  - the time of formation of the BBH
  //  - the time of core collapse of the cluster

  //We compute the number densities
  // rho_star = M_core/R_core^3 ==> n_star = rho_star/mstar_avg
  // R_hier = max(R_inf, R_wand) 
  // with the wandering radius R_wand being the region of the core where we can find the BH
  // and the influence radius being the region of the core where the BHs dominates the interactions
  //                   | alpha, with alpha = m_hier/m_core
  // R_hier = R_core x |
  //                   | mu exp[-0.5], with mu = m_hier/mstar_avg

  // Thus, we define the number density of hierarchical BHs as:
  // rho_hier = n_hier / (R_hier)^3 from Di Cinto et al. 2023
  // n_hier = rho_hier / (R_hier)^3 from Di Cinto et al. 2023

  // with m_hier = max( 2 * mstar_avg, s[2] )

  // in this way we account for the average mass of the 1-g BHs in the cluster for low gens and for the higher gen mergers    Functions func;
    
    double IR;
    
    // Let's define the quantities intervining in the interaction rate
    double mstar_avg = 17.4 - 4.0 * log10(tbbhform/trelax);  // Dragon-II paper II, eq. 8
    
    // Masses involved in the interaction rate
    double m_bin = m1 + m2;
    double m_hier = gen * mstar_avg;

    // Let's define the influence sphere of the hierarchical BH, see Di Cinto et al. 2023
    double alpha =  m_hier/mcore; 
    double mu = m_hier/mstar_avg;

    double r_inf = rcore * alpha;
    double r_wand = rcore * pow(mu, -0.5);
    double r_hier = max(r_inf, r_wand);

    // Let's define the densities of stars and hierarchical BHs
    double n_star = mcore / (mstar_avg * pow(rcore, 3)); // Number density of stellar BHs in the core
    double n_hier = nbhs / pow(r_hier, 3); // Number density of hierarchical BHs in the core
    
    // Now we can compute the interaction rate
    IR = (n_hier / n_star) * pow(m_hg / mstar_avg, 3. / 2.) * (m_bin + m_hg) / (m_bin + mstar_avg);

    if (IR>1.0) IR = 1.0; // Cap the interaction rate to 1.0
    
    return IR;

    }


     void Functions::evolve_bhs(vector<double>& nbhs,           // total of bhs
                  double n_bin,                   // fraction of bhs in binaries
                  const vector<double>& gwK,      // kick distribution
                  const vector<double>& gwK_cdf,  // kick cdf
                  double vesc,                     // escape velocity of the cluster
                  int gen2){                       // generation of the hierarchical merger (to set the max gen to update)

      // This function evolves the population of BHs in the cluster considering:
      // - the fraction of BHs in binaries (n_bin)
      // - the kick distribution (gwK and gwK_cdf)
      // - the escape velocity of the cluster (vesc)
      // - the generation of the hierarchical merger (gen2)

      // We'll append one zero slot to allow the formation of a next generation
      nbhs.push_back(0.0);                // reserve the next generation (zero)

      // old snapshot
      const vector<double> old = nbhs;
      const size_t oldS = old.size();     // e.g. 2 when we have [total, gen1]

      // ---- retention fraction ret_fract = P(v_kick <= vesc) ----
      double ret_fract = 0.0;
      
      auto it = upper_bound(gwK.begin(), gwK.end(), vesc);
      size_t idx = (it == gwK.begin()) ? 0
                : (it == gwK.end())   ? gwK_cdf.size() - 1
                                      : size_t(it - gwK.begin() - 1);
      ret_fract = gwK_cdf[idx];

      // Update gen-1 explicitly (loses a fraction n_bin because of mergers)
      nbhs[1] = old[1] * (1.0 - n_bin);

      // Update gens 2..oldS
      // (flow within and from previous gen, aka BHs not merging of the current gen and merging ones of the previous)
      for (size_t g = 2; g < oldS; ++g) {
          const double stay      = old[g]     * (1.0 - n_bin);
          const double from_prev = old[g - 1] * 0.5 * n_bin * ret_fract;
          //cout << "gen " << g << ": stay " << stay << ", from_prev " << from_prev << endl;
          nbhs[g] = stay + from_prev;
      }

      // Recompute total as sum over all actual generations (exclude the newly appended zero at the end)
      double total = 0.0;
      for (size_t g = 1; g < oldS; ++g) total += nbhs[g];
      nbhs[0] = total;
      
      // Debug output to check for negative burnt BHs
      double burnt_bhs = old[0] - nbhs[0];
      if (burnt_bhs < 0.0) {
          cout << "Warning: negative burnt BHs: " << burnt_bhs << endl;
          cout << "Cluster properties: vesc " << vesc << ", n_bin " << n_bin << ", ret_fract " << ret_fract << ", gen2 " << gen2 << endl;
          cout << " old total: " << old[0] << ", new total: " << nbhs[0] << endl;
          // Debug output for nbhs and old arrays before exiting
          cout << "nbhs array contents:" << endl;
          for (size_t i = 0; i < nbhs.size(); ++i) {
            cout << "nbhs[" << i << "] = " << nbhs[i] << endl;
          }
          cout << "old array contents:" << endl;
          for (size_t i = 0; i < old.size(); ++i) {
            cout << "old[" << i << "] = " << old[i] << endl;
          }
      }

      return;
  }

  void Functions::DiCarlo_BHs(double* mpri, double* msec, double Z, bool processed, string uppergap, double fupgp, double a_gp, double mass_gap, string upgtp){
    //This section serves for the binary component masses --- need to be added also in the hierarchical merger chain
    //We follow Di Carlo+2020 for the upper mass gap treatment
    // The function modifies the masses of the binary components considering:
    // - the metallicity of the cluster (Z)
    // - if the BHs are processed or not (processed), i.e. if they have been already modified in the hierarchical merger chain
    // - the user choice for the upper mass gap treatment (uppergap)
    // - the fraction of mergers in the upper mass gap (fupgp)
    // - the slope of the power-law distribution for the masses in the gap (a_gp)
    // - the mass of the upper mass gap (mass_gap)

    //We assume a power-law distribution for the masses in the gap, with slope a_gp
    Functions func;

    //From DRAGON: distribution of mergers with no compo. in the gap, both compo. in the gap, one compo. in the gap
    double pbelow = 45./78. ;
    double pupper = 12./78. ;
    double pbelup = 21./78. ;

    double prob_ugp = func.rnd();
    double prob_fgp = func.rnd();
    
    string UP = upgtp;
    double fUP= fupgp;
          
    if(UP == "dicarlo")
      fUP = 0.01 * (1. + 5.797 * exp(Z / 0.0002 * log(5./5.797)));
    
    if(uppergap == "yes" && prob_ugp >= pbelow && prob_fgp < fUP){
        
      double p_gp;
      p_gp = func.rnd();
      double m1_gp = pow(p_gp * pow(100.,1.-a_gp) + (1.-p_gp)*pow(50., 1.-a_gp), 1./(1.-a_gp));
      p_gp = func.rnd();
      double m2_gp;

      //The following make results DRAGON-II like!
      if(m1_gp < 100.)
        m2_gp = m1_gp * (0.4 + 0.6*func.rnd()); 
      else
        m2_gp = m1_gp * pow(m1_gp / 65., -1.78) * (1. + 0.2*(-1. + 2.*func.rnd()));

      if(m2_gp > m1_gp){
        double temp_gp = m1_gp;
        m1_gp = m2_gp;
        m2_gp = temp_gp;
      }
        
      // We need to implement something for the spins too ...
      if(prob_ugp >= pbelow && prob_ugp < pbelow + pbelup){

        if(*mpri < mass_gap && !processed) *mpri = m1_gp;

        else if (*msec < mass_gap) *msec = m2_gp;
        
        else if (!processed){
          //cout<<"Both objects above the gap!"<<endl; //This may happen in multiple mergers?
        }
      }	      
      
      else if(prob_ugp > pbelow + pbelup){

        if(*mpri < mass_gap && !processed) *mpri = m1_gp;

        if(*msec < mass_gap) *msec = m2_gp;

      }
    }
    return;
  }

  