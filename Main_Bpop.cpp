#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <vector>
#include "funcGR.h"
#include <ctime>
#include <stdlib.h>
#include <omp.h>
#include <chrono>

//#define THREAD 4
#define Hubble 14.E9 //13.803E9

// DATAFILES (Metal. distri, Single BHs, Binary BHs)
#define PREDIR "../"
#define zPATH   "gallazzi05ZDATA.ttt"
#define SINGPTH "rapid_M20/" // "DATI_SingleBH/"
#define PATH    "rapid_M20/" // "DATI_GiaMap18/"
#define PATHSIN "DATI_SingleBH/"

// GLOBAL
#define N        1000000
#define mmax     150.
#define mmin     18.5
#define mslope  -2.35
#define Zsun     0.02

//DYNAMICAL FRACTIONS
#define DynOvTot  0.9
#define pYC 0.75
#define pGC 0.15
#define pNC 0.10

#define uppergap "yes"
#define bhseed   "bifrost"
#define bifZ     0.001
#define bhpisn   270.
#define fupgp    0.15
#define mass_gap  60.0
#define upgtp   "dicarlo"
#define SFRTYPE_ISO "bigbang" //"continuous" //"MF17" //
#define SFRTYPE_CLU "bigbang" //"continuous" //"EB18_MF17" //


#define mixing  0.5
#define fbin    1.0

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

using namespace std;

void singBHt_mix(double mssx[], double msdx[], double mbsx[], double mbdx[], double tbsx[], double tbdx[], double vbsx[], double vbdx[], double mbhmix[][nsize], double tbhmix[][tsize], double vbhmix[][vsize], double mslp, double *sing_out, double saximus_mix, double sinimus_mix, double maximus_mix, double minimus_mix, double vescape){

  Functions func;

  int nsafe;
  
  double P = func.rnd();
  double mzams = pow((P*pow(saximus_mix,1.+mslp) + (1.-P)*pow(sinimus_mix,1.+mslp)),1./(1.+mslp));
  int idx = -1;
  for(int ii=0;ii<bin_st;ii++){
    if(mzams >= mssx[ii] && mzams < msdx[ii]){
      idx = ii;
      break;
    }
  }


  
  double maxprob;
  double maxtest;
  double maxpoint;
  int idy;

  maxtest = 1.E30;
  maxpoint= -1;
  idy = -1;
  nsafe = 0;
  maxprob = -1;     
  for(int jj=0;jj<nsize;jj++){    
    if(maxprob < mbhmix[idx][jj])
      maxprob = mbhmix[idx][jj];
  }
  do{
    double ddy = 1.*nsize*func.rnd();
    idy = ddy;
    maxtest = maxprob * func.rnd();
    maxpoint= mbhmix[idx][idy];
    /*if(nsafe > 100){
	cout<<"Warning mass. : "<<vsize_max<<" "<<vbsx[vsize_max]<<" "<<vbdx[vsize_max]<<" "<<vescape<<endl;
	}*/
    nsafe ++;
  }while(maxtest > maxpoint);	
  
  double mblack = mbsx[idy] + (mbdx[idy] - mbsx[idy]) * func.rnd();
  
  
  int vsize_max = -1;
  for(int jj=0;jj<vsize;jj++){
    if(vbsx[jj] <=vescape)
      vsize_max = jj;	    
    if(vbsx[jj] > vescape)
      break;    
  }
  if(vsize_max <= 0){
    //cout<<"WARNING: no BHs can be retained in this cluster due to natal kick "<<vsize_max<<" "<<vbsx[vsize_max]<<" "<<vbdx[vsize_max]<<" "<<vescape<<endl;
    vsize_max = 0;    
  }

  //vsize_max = vsize;
  
  //cout<<vsize_max<<" "<<vbsx[vsize_max]<<" "<<vbdx[vsize_max]<<" "<<vescape<<endl;
  //exit(0);
  nsafe = 0;
  maxprob = -1;
  for(int jj=0;jj<vsize_max;jj++){
    if(maxprob < vbhmix[idy][jj])
      maxprob = vbhmix[idy][jj];
  }
  double vblack;
  if(maxprob > 0){
    maxtest = 1.E30;
    maxpoint= -1;
    int idz = -1;
    do{
      double ddz = 1.*vsize_max*func.rnd();
      idz = ddz;
      maxtest = maxprob * func.rnd();
      maxpoint= vbhmix[idy][idz];
      nsafe ++;
      /*if(nsafe > 100){
	cout<<"Warning vel. : "<<vsize_max<<" "<<vbsx[vsize_max]<<" "<<vbdx[vsize_max]<<" "<<vescape<<endl;
	}*/
    }while(maxtest > maxpoint);
    
    vblack = vbsx[idz] + (vbdx[idz] - vbsx[idz]) * func.rnd();
  }
  else
    vblack = 0.0;
  
  

  nsafe = 0;
  maxprob = -1;
  for(int jj=0;jj<tsize;jj++){
    if(maxprob < tbhmix[idy][jj])
      maxprob = tbhmix[idy][jj];
  }
  double tblack;
  if(maxprob > 0){
    maxtest = 1.E30;
    maxpoint= -1;
    int idz = -1;
    nsafe = 0;
    do{
      double ddz = 1.*tsize*func.rnd();
      idz = ddz;
      maxtest = maxprob * func.rnd();
      maxpoint= tbhmix[idy][idz];
      /*if(nsafe > 100){
	cout<<"Warning time. : "<<vsize_max<<" "<<vbsx[vsize_max]<<" "<<vbdx[vsize_max]<<" "<<vescape<<endl;
	}*/
      nsafe ++;
    }while(maxtest > maxpoint);
    
    tblack = tbsx[idz] + (tbdx[idz] - tbsx[idz]) * func.rnd();
  }
  else
    tblack = 0.0;
  
       
  sing_out[0] = mblack;
  sing_out[1] = tblack;
  sing_out[2] = vblack;

  //cout<<mzams<<" "<<mblack<<" "<<tblack<<" "<<vblack<<endl;
  
  return ;
}
    
int main(){
  srand(time(0));
  Functions func;

  //int NTH = THREAD;
  using clock = std::chrono::system_clock;
  using sec = std::chrono::duration<double>;
  const auto before = clock::now();

  double sigmaZ = sigma_metal;
  string metal_dis = sigma_distri;
  string predir = PREDIR;

  double sigma_sma = SSMA;
  
  int Niso_real = 0;
  int Ndyn_real = 0;
  int Nyou_real = 0;
  int Nglo_real = 0;
  int Nnuc_real = 0;
  

  if(CLfill != "under" && CLfill != "over" && CLfill != "critical" && CLfill != "mix" && CLfill != "noevo" && CLfill != "GG23"){
    cout<<"Select CLfill = under/over/critical/mix/GG23/noevo "<<endl;
    cout<<"You've chosen: "<<CLfill<<endl;
    exit(0);
  }

  //From DRAGON: distribution of mergers with no compo. in the gap, both compo. in the gap, one compo. in the gap
  double pbelow = 45./78. ;
  double pupper = 12./78. ;
  double pbelup = 21./78. ;


  double probGC = pGC;
  double probYC = pYC;
  double probNC = pNC;

  if( cluster_test == "yes" ){
    if( cluster_test_env == "YC" ) {
      probYC = 1.0;
      probGC = 0.0;
      probNC = 0.0;
    }
    else if (cluster_test_env == "GC" ) {
      probYC = 0.0;
      probGC = 1.0;
      probNC = 0.0;
    }
    else if (cluster_test_env == "NC" ) {
      probYC = 0.0;
      probGC = 0.0;
      probNC = 1.0;      
    }				  
  }
  
  //  func.test();
  //  exit(0);
  
  if(MonoZ=="yes" || cluster_test=="yes"){
    cout<<"WARNING: cluster test mono Z mode enabled ";
    if(DynOvTot != 1.0 && probGC != 1.0){
      cout<<"please set DynOvTot and probGC to 1";
    }
    cout<<endl;
      
  }

  ifstream in;
  
  int redline = 100;
  double reds[redline],age[redline],lkbk[redline];
  in.open(predir+"redshift_time.ttt");
  if(!in.is_open()){
    cout<<"File redshift vs time not found"<<endl;
    exit(0);
  }
  double d1, d2, d3;
  for(int i=0;i<redline;i++){
    in>>d1>>d2>>d3;
    reds[i] = d1;
    age[i]  = d2;
    lkbk[i] = d3;
    //cout<<reds[i]<<" "<<age[i]<<" "<<lkbk[i]<<endl;
  }
  in.close();
  
  string isospin = spinlb;  
  string dynaS   = spinlb;
  
  if(spinlb == "mixed"){
    isospin = "maxwellian02";
    dynaS   = "fuller";
  }
  else if(spinlb== "mixedhigh"){
    isospin = "gaussian";
    dynaS   = "fuller";
  }
  
  string isolS = isospin;
  

  string cmdstr = "rm Catalogue*txt total_mremn_distri*txt";
  char *cmd;
  cmd = new char [cmdstr.length()+1];
  strcpy(cmd,cmdstr.c_str());
  int ck;
  delete [] cmd;


  string sfr_iso = SFRTYPE_ISO;
  string sfr_clu = SFRTYPE_CLU;
  string sfr_clu_gc = sfr_clu;
  string sfr_clu_nc = sfr_clu;
  if(sfr_clu == "EB18_MF17"){
    sfr_clu_gc = "EB18";
    sfr_clu_nc = "MF17";
  }
  

  if(cluster_test == "yes"){
    if(cluster_test_env == "NC")
      sfr_clu = "bigbang";
    else
      if(cluster_test_env == "YC")
	sfr_clu = "MF17";
      else
	sfr_clu = "KR13";
  }
  int Nsrc = N;
    
  double *Spinning;
  Spinning = new double [7];
 
  Spinning[3] = -1.E30;
  
  string path   = predir+PATH;
  string pathse = "";
  string pathsp = "";
  
 /* 
  //METALLICITY AVAILABLES 0.0002 - 0.0004 - 0.0008 - 0.0012 - 0.0016 - 0.002 - 0.004 - 0.006 - 0.008 - 0.012 - 0.016 - 0.02 
  double met[13];
  met[0]  = 0.0002;
  met[1]  = 0.0004;
  met[2]  = 0.0008;
  met[3]  = 0.0012;
  met[4]  = 0.0016;
  met[5]  = 0.002;
  met[6]  = 0.004;
  met[7]  = 0.006;
  met[8]  = 0.008;
  met[9]  = 0.012;
  met[10] = 0.016;
  met[11] = 0.02;
  met[12] = 0.03;

  double mis[13];
  for(int i = 0;i<13;i++)mis[i] = 0.0;
  
  ofstream out;

  double metdyn[13];
  metdyn[0]  = 0.0002;
  metdyn[1]  = 0.0004;
  metdyn[2]  = 0.0008;
  metdyn[3]  = 0.0012;
  metdyn[4]  = 0.0016;
  metdyn[5]  = 0.002;
  metdyn[6]  = 0.004;
  metdyn[7]  = 0.006;
  metdyn[8]  = 0.008;
  metdyn[9]  = 0.012;
  metdyn[10] = 0.016;
  metdyn[11] = 0.02;
  metdyn[12] = 0.03;
  */
    //METALLICITY AVAILABLES 0.0002 - 0.0004 - 0.0008 - 0.0012 - 0.0016 - 0.002 - 0.004 - 0.006 - 0.008 - 0.012 - 0.016 - 0.02 
	double met[13];
	met[0]  = 0.0002;
	met[1]  = 0.0003;
	met[2]  = 0.0004;
	met[3]  = 0.0007;
	met[4]  = 0.0010;
	met[5]  = 0.0014;
	met[6]  = 0.002;
	met[7]  = 0.004;
	met[8]  = 0.007;
	met[9]  = 0.010;
	met[10] = 0.014;
	met[11] = 0.02;
	met[12] = 0.03;
  
	double mis[13];
	for(int i = 0;i<13;i++)mis[i] = 0.0;
	
	ofstream out;
  
	double metdyn[13];
	for (int i = 0; i < 13; i++) {
		metdyn[i] = met[i];
	}

  string singpthA = predir+SINGPTH;  
  double ndx;
  double apri;
  double asec;
  double mpri;
  double msec;
  double tpri;
  double tsec;
  double kpri;
  double ksec;
  double smaiso;
  
  double vpri;
  double vsec;
  
  string align;

  double rnd;
  string label = "0";

  double qmin, recy;
  string cluster = "none";
  vector<double> mpost;
 
  align = "whatever";

  // CREATING ESCAPE VEL ARRAYS //

  vector<double> YCrx,YCry,YCmx,YCmy;
  vector<double> GCrx,GCry,GCmx,GCmy;
  vector<double> NCrx,NCry,NCmx,NCmy;

  for(int ciname=0;ciname<3;ciname++){
    string cltype;
    if(ciname == 0)
      cltype = "YC";
    if(ciname == 1)
      cltype = "GC";
    if(ciname == 2)
      cltype = "NC";
    
    for(int cjname=0;cjname<2;cjname++){
      
      string clvari;
      if(cjname == 0)
	clvari = "radi";
      else if(cjname == 1)
	clvari = "mass";

      
      string clufile = predir+"ClustersRhM/scottrule_"+clvari+cltype+".ttt";
      cout<<"Opening "<<clufile<<endl;
      ifstream cluIO;
      cluIO.open(clufile);
      if(!cluIO.is_open()){
	cout<<"File "<<clufile<<" not found"<<endl;
	exit(0);    
      }
      do{
	double x,y;
	cluIO>>x>>y;
	if(cltype == "YC"){
	  if(clvari=="radi"){
	    YCrx.push_back(x);
	    YCry.push_back(y);
	  }
	  else if(clvari=="mass"){
	    YCmx.push_back(x-2.5);
	    YCmy.push_back(y);
	  }
	  
	}
	else if(cltype=="GC"){
	  if(clvari=="radi"){
	    GCrx.push_back(x);
	    GCry.push_back(y);
	  }
	  else if(clvari=="mass"){
	    GCmx.push_back(x);
	    GCmy.push_back(y);
	  }
	  
	}
	else if(cltype=="NC"){
	  if(clvari=="radi"){
	    NCrx.push_back(x);
	    NCry.push_back(y);
	  }
	  else if(clvari=="mass"){
	    NCmx.push_back(x);
	    NCmy.push_back(y);	    
	  }
	}
	else{
	  cout<<"Error Cl esc "<<endl;
	  exit(0);
	}
      }while(!cluIO.eof());    
    }
  }



  cout<<"Real stuff starting "<<endl;

  double par[20];

  
  double *Mrem;
  Mrem = new double[Nsrc];
  double *Srem;
  Srem = new double[Nsrc];
  double *Xrem;
  Xrem = new double[Nsrc];
  double *Krem;
  Krem = new double[Nsrc];
  double *Cosa;
  Cosa = new double[Nsrc];
  double *Cosb;
  Cosb = new double[Nsrc];
  double *Cosg;
  Cosg = new double[Nsrc];

   
  double *Z;
  Z = new double [Nsrc];
  double *Zi;
  Zi = new double [Nsrc];

  int *Npar;
  Npar = new int [numZ];
  for(int i=0;i<numZ;i++) Npar[i] = 0.0;



  
  vector<double> miso1,miso2;
  vector<double> aiso1,aiso2;
  vector<double> sma_iso;
  double P;

  //CREATE THE CUMULATIVE AND SELECT FROM THAT... (Gallazzi+05)
  string zpth = predir+zPATH;
  ifstream gin;
  gin.open(zpth);
  if(!gin.is_open()){
    cout<<"DATA FILE FROM GALLAZZI DOESN'T FOUND..."<<endl;
    exit(0);
  }
  vector<double> zgal;
  vector<double> ngal;
  double cnt = 0;
  string line;
  do{
    getline(gin,line);
    cnt += 1;
  }while(!gin.eof());
  cnt -= 2;
  gin.close();
  gin.open(zpth);
  double a1,a2;
  double a3 = 0.0;
  for(int ll=0;ll<cnt;ll++){
    gin>>a1>>a2;
    if(zdivi != 0)
      a1 = a1 - log10(zdivi);
    if(correction == "no" || correction == "No" || correction == "NO"){
      a3+=a2;
    }
    else{
      a3 += (a2*pow(Zsun*pow(10.,a1),slope));
    }
    zgal.push_back(a1);
    ngal.push_back(a2);
  }
  gin.close();
  
  vector<double> Nz;
  double a4 = 0; double a5=0;
  ofstream gout;gout.open("gallazzi_cumul.txt");
  for(int ll=0;ll<zgal.size();ll++){
    double ncorr;
    if (correction == "no" || correction == "No" || correction == "NO"){
      ncorr = ngal[ll]/a3;
    }
    else{
      ncorr = ngal[ll]/a3*pow(Zsun*pow(10.,zgal[ll]),slope);
    }
    
    a4 += ncorr;
    a5  = ncorr;

    Nz.push_back(a4);
    gout<<pow(10.,zgal[ll])*Zsun<<" "<<Nz[ll]<<" "<<a5<<endl;

  }
  gout.close();

  
  int Ndyn = 0;
  int Nyc  = 0;
  int Ngc  = 0;
  int Nnc  = 0;
  int Niso = 0;
  int Ncnt = 0;
  
  /*do{
    double csu = func.rnd();
    if(csu>DynOvTot) Niso += 1;
    else{
      double csv = func.rnd();
      if(csv < probGC)                            Ngc += 1;
      else if(csv>probGC && csv<probGC+probYC)    Nyc += 1;
      else                                        Nnc += 1;
    }
    Ncnt += 1;

  }while(Ncnt<Nsrc);  
  
  Nsrc = Niso + Nyc + Ngc + Nnc;*/

  Niso = Nsrc * (1.-DynOvTot);
  Nyc  = Nsrc * DynOvTot * probYC;
  Ngc  = Nsrc * DynOvTot * probGC;
  Nnc  = Nsrc - (Niso + Nyc + Ngc);

  cout<<"Creating N = "<<Nsrc<<" "<<Niso<<" "<<Nyc<<" "<<Ngc<<" "<<Nnc<<endl;

  double Z0 = 1.E30;
  double Z1 = -1.E30;


  //FORMATION TIMES

  vector<double> tfor;
  for(int i=0;i<Niso+Nyc;i++){
    double zred0 = func.sfr_red(sfr_iso);      
    double toff = 1.E9*func.inter(zred0, reds, age, redline);
    tfor.push_back(toff);
  }

  for(int i=0;i<Ngc;i++){
    double zred0 = func.sfr_red(sfr_clu_gc);      
    double toff  = 1.E9*func.inter(zred0, reds, age, redline);
    tfor.push_back(toff);
  }
  for(int i=0;i<Nnc;i++){
    double zred0 = func.sfr_red(sfr_clu_nc);      
    double toff  = 1.E9*func.inter(zred0, reds, age, redline);
    tfor.push_back(toff);
  }
  
  //ASSIGNING METALLICITIES to ISO and CLU BBHs
  double z;
  cout<<"Assigning metallicities"<<endl;

  
  for(int i=0;i<Niso+Nyc;i++){
    
    double tof = tfor[i] / 1.E9;      
    double zinit = z;
    double red_del = func.inter(tof, age, reds, redline);
    
    double logz;
    double logz_me;            
    logz_me = func.metcor(metal_dis, sigmaZ, red_del);
    double logz1 = func.Gaussian_normal(log10(met[0]/Zsun),log10(met[11]/Zsun), logz_me, sigmaZ) + log10(Zsun);
    
    logz = logz1;

    z = pow(10.,logz);  
    if(z < 0.9*met[0])
      z = met[0];
    if(z > 1.2*met[11])
      z = met[11];

    Z[i] = z;
    if(Z[i]>Z1) Z1 = Z[i];
    if(Z[i]<Z0) Z0 = Z[i];

  }

  int nout=0;
  for(int i=Niso+Nyc;i<Niso+Nyc+Ngc+Nnc;i++){

    double tof = tfor[i] / 1.E9;      
    double zinit = z;
    
    double red_del = func.inter(tof, age, reds, redline);
    double logz;
    double logz_me;            
    logz_me = func.metcor(metal_dis, sigmaZ, red_del);
    double logz1 = func.Gaussian_normal(log10(met[0]/Zsun),log10(met[11]/Zsun), logz_me, sigmaZ) + log10(Zsun);

    logz = logz1;
    
    z = pow(10.,logz);  
    if(z < 0.9*met[0])
      z = met[0];
    if(z > 1.2*met[11])
      z = met[11];

    Z[i] = z;
    if(Z[i]>Z1) Z1 = Z[i];
    if(Z[i]<Z0) Z0 = Z[i];
      
  }
  //----------------------------------------------------------
  
  
  if(MonoZ == "yes")
    for(int i=0;i<Nsrc;i++){
      Z[i] = mono_Z;
      Z0 = Z[i];
      Z1 = Z[i];
    }

  
  //   CLOSEST DISCRETE Z FOR FILES READING?
  
  int NzIS = 0;
  int NzDY = 0;
  int NzGC = 0;
  int NzYC = 0;
  int NzNC = 0;
  out.open("met_dist.txt");

  Z1 = log10(1.2*Z1);
  Z0 = log10(0.9*Z0);
  int numbin_Z = numZ;
  double dZ = (Z1-Z0)/(double)numbin_Z;
  double Zmin = Z0;
  for(int k=0;k<numbin_Z;k++){
    Z0 = Zmin + k*dZ;
    Z1 = Zmin + (k+1)*dZ;
    for(int i=0;i<Nsrc;i++){  
      if(Z[i]<pow(10.,Z1) && Z[i]>=pow(10.,Z0)){
	if(i<Niso) NzIS += 1;
	else{
	  NzDY += 1;
	  if(i < Niso + Nyc)
	    NzYC += 1;
	  else if(i >= Niso + Nyc && i < Niso + Nyc + Ngc)
	    NzGC += 1;
	  else
	    NzNC += 1;
	}
	
	Npar[k] += 1;
      }
      /*else if(Z[i]==met[numZ-1]){
	if(i<Niso+Nyc) NzIS += 1;
	else           NzGC += 1;	
	Npar[numZ-1] += 1;	
	} */     
    }
    out<<0.5*(pow(10.,Z0)+pow(10.,Z1))<<" "<<Npar[k]<<" "<<NzIS<<" "<<NzDY<<" "<<pow(10.,Z0)<<" "<<pow(10.,Z1)<<" "<<NzYC<<" "<<NzGC<<" "<<NzNC<<endl; //Printout metallicity bins
    NzIS = 0;
    NzDY = 0;
    NzGC = 0;
    NzYC = 0;
    NzNC = 0;

    
  }
  out.close();

  int itot = 0;


  zgal.erase(zgal.begin(),zgal.end());
  ngal.erase(ngal.begin(),ngal.end());
  Nz.erase(Nz.begin(),Nz.end());
  
  //cout<<"Metallicity and Number of sources"<<endl;
  //for(int k=0;k<numZ;k++)cout<<met[k]<<" "<<Npar[k]<<endl;
  
  
  // ----------------------------------------------------------------------------  

  cout<<"    PREPARING THE CATALOGUE   "<<endl;
  //  exit(0);

  double zmer, zfor, zsmbh;
  
  int Ngap = 0;
  string REC = "RecNO";
  cout<<"Header: src No., metallicity, scenario, cluster type"<<endl;
  ofstream out2;

  string outname = "Catalogue.txt";
  out.open(outname.c_str());
  out<<"#ID Metal Nrec EnvType lab m1[Msun] m2[Msun] a1 a2 Mfin[Msun] afin xeff vGW[km/s] tfor[yr] tlast_mer[yr] Mclu_t0[Ms] Rclu_t0[pc] Vesc[km/s] BinaryStatus aeje[AU] aGW[AU] nBHs Mcore_th[Ms] rcore_th[pc] redshift_merger redshift_formation tSMBH[yr] redshiftSMBH mprog[Ms] eccentricity semimajoraxis[AU] acrit[AU] tmerger[yr] cos(angle_s1s2) cos(angle_s1L) cos(angle_s2L)"<<endl;

  out2.open("Catalogue_clean.txt");

  cout<<"    ISOLATED BINARIES   "<<endl;

  ofstream hout;
  hout.open("Larger_than_tH.txt");
 

  
  for(int k=0;k<numZ;k++)Npar[k] = 0;
  for(int i=0;i<Niso;i++)Zi[i]=Z[i];
  

  vector<double> arr[numZ];
  
  for(int i=0;i<Niso;i++){
    if(Z[i] <= met[0]){
      Z[i] = met[0];
      Npar[0] += 1;
      arr[0].push_back(i);      
      continue;
    }
    if(Z[i] >= met[numZ-1]){
      Z[i] = met[numZ-1];
      Npar[numZ-1] += 1;
      arr[numZ-1].push_back(i);
      continue;
    }
    
    for(int k=0;k<numZ;k++){
      double dcen, dup, dlow, dZsx, dZdx;
      dcen = met[k];
      if(k < numZ-1)
	dup = met[k+1];
      else
	dup = met[k] * 1.1;
            
      if(k > 0)
	dlow = met[k-1];
      else
	dlow = 0.9 * met[k];
            
      dZsx = (dcen-dlow)/2.0;
      dZdx = (dup-dcen)/2.0;

      
      if(Z[i] >= met[k] - dZsx && Z[i] < met[k] + dZdx){
	Npar[k] += 1;
	Z[i] = met[k];	
	arr[k].push_back(i);
	break;
    
      }
      /*else if(Z[i]==met[numZ-1]){
	Npar[numZ-1] += 1;	
	} */
    }    
  }

  double tdel;
  //cout<<"Metallicity and Number of sources for isolated channel "<<endl;
  //for(int k=0;k<numZ;k++)cout<<met[k]<<" "<<Npar[k]<<endl;
  double mobs;

  
  for(int k=0;k<numZ;k++){
    cout<<"      -- Metallicity "<<met[k]<<" "<<Npar[k]<<endl;
    if(Npar[k]==0) continue;


    //READING MICHELA'S DATA FILES 
    double maxbhiso = -1.0;
    double minbhiso = 10000.0;
    
    stringstream s;
    s<<met[k];

    vector<double> tdel_iso;


    string filepath = path + pathse + pathsp + "/data_" + s.str() + pathsp + ".dat";

    in.open(filepath.c_str());
    if(!in.is_open()){
      cout<<"Data file "<<filepath<<" not found"<<endl;
      exit(0);
    }
    int chk = 0;
    do{	
      for(int kk=0;kk<kpar;kk++)in>>par[kk];

      
      if(par[1]>par[2]){
	miso1.push_back(par[1]);
	miso2.push_back(par[2]);
	aiso1.push_back(par[8]);
	aiso2.push_back(par[9]);	
      }
      else{
	miso1.push_back(par[2]);
	miso2.push_back(par[1]);
	aiso1.push_back(par[9]);
	aiso2.push_back(par[8]);
      }
      sma_iso.push_back(par[7]*6.9E8/1.5E11);
      tdel_iso.push_back(par[6]*1.E6);

      
      if(par[2] > maxbhiso) maxbhiso = par[2];
      if(par[1] > maxbhiso) maxbhiso = par[1];
      if(par[2] < minbhiso) minbhiso = par[2];
      if(par[1] < minbhiso) minbhiso = par[1];
	 
      chk += 1;
    }while(!in.eof());
    in.close();

    int safe = 0;
    double deltabh = 1.E30;
    double bhcandi = 1.E30;
    
    REC = "RecNO";
    mpri = -10000.0;

    //int manuel = NTH;
    //#pragma omp parallel for num_threads(4) 
    for(int i=0; i<Npar[k];i++){
      int jump = 0;
      REC = "RecIS";
      ndx = INDEX_ALIGN;
      align = "whatever";
      double extract = chk * func.rnd();
      int ext = extract;
      P = func.rnd();
      mobs = pow( pow(maxbhiso,1.+obslope)*(1.-P) + P*pow(minbhiso,1.+obslope)   ,1./(1.+obslope));

      
      safe = 0;
      deltabh = 1.E30;
      bhcandi = 1.E30;
      //BH SPINS SELECTION (TO BE DECIDED)    
      if(isospin != "giacobbo" || isospin != "Giacobbo"){
	apri   = func.spin(mpri,isospin);
	asec   = func.spin(msec,isospin);	
      }
      else{
	apri = aiso1[ext];
	asec = aiso2[ext];
      }


      
      do{
	extract = chk * func.rnd();
	ext = extract;
	if(abs(mpri-mobs) < deltabh){
	  deltabh = abs(mpri-mobs);
	  bhcandi = mpri;
	}
	
	safe += 1;
	if(safe > 1000){
	  mpri = bhcandi;
	  mis[k] += 1./Npar[k];
	  break;
	}

	if(miso1[ext] < mobs*1.2 && miso1[ext] > mobs*0.8 && tdel_iso[ext]+tfor[i] > Hubble){
	  hout<<itot<<" 0.0 0.0 0.0 1.E30 1.E30 0.0 0.0 0.0 0.0 0.0 "<<tdel_iso[ext]<<" "<<safe<<" "<<Z[k]<<" none "<<miso1[ext]<<" "<<miso2[ext]<<" "<<apri<<" "<<asec<<" 0.0 "<<sma_iso[ext]<<" 0.0 0"<<" "<<tdel_iso[ext]<<endl;
	  jump = 1;
	  break;
	  
	}

	  
      }while((miso1[ext] > mobs*1.2 || miso1[ext] < mobs*0.8) || tdel_iso[ext]+tfor[i] > Hubble);

      
      if(jump == 1){
	itot++;       
	continue;
      }
      
      
      mpri = miso1[ext];
      msec = miso2[ext];
      tdel = tdel_iso[ext];


      Spinning[0] = 0.0;
      Spinning[1] = 0.0;
      Spinning[2] = 0.0;
      Spinning[3] = 0.0;
      Spinning[4] = 0.0;
      Spinning[5] = 0.0;
      Spinning[6] = 0.0;
      if(mpri > 0.0 && msec > 0.0)
	func.SREM2(ndx, apri, asec, mpri, msec, align, Spinning);
	
      Srem[i] = Spinning[0];
      Xrem[i] = Spinning[1];
      Mrem[i] = Spinning[2];
      Krem[i] = Spinning[3];
      Cosa[i] = Spinning[4];
      Cosb[i] = Spinning[5];
      Cosg[i] = Spinning[6];
      
      if(mpri < msec){
	double mpri_sec = mpri;
	mpri = msec;
	msec = mpri;
      }

      //cout<<"Check3"<<endl;

      double tform = tfor[arr[k][i]];
      
      tdel += tform;
      smaiso = sma_iso[ext];
	
      if(tdel < 1.35E10){
	zmer = func.inter(tdel / 1.E9, age, reds, redline);
	zfor = func.inter(tform / 1.E9, age, reds, redline);
      }
      else{
	zmer = func.zred(tdel/1.E9);
	zfor = func.zred(tform/1.E9);
      }
      
      


      
      //Dec 23: a space was missing between 0.0 and zmer, resulting in error for isolated binaries ...
      out<<itot<<" "<<met[k]<<" "<<label<<" "<<cluster<<" "<<REC<<" "<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<Mrem[i]<<" "<<Srem[i]<<" "<<Xrem[i]<<" "<<Krem[i]<<" "<<tform<<" "<<tdel;
      out<<" 0.0 0.0 0.0 none "<<smaiso<<" "<<smaiso<<" 0.0 0.0 0.0 "<<zmer<<" "<<zfor<<" -1 -1 " <<mpri<<" 0.0 "<<smaiso<<" "<<smaiso<<" "<<tdel<<" "<<Cosa[i]<<" "<<Cosb[i]<<" "<<Cosg[i]<<endl;
      

      Niso_real++;
      
      if(mpri!=0.0)	
	out2<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<Mrem[i]<<" "<<Srem[i]<<" "<<Xrem[i]<<" "<<Krem[i]<<" "<<tform<<" "<<tdel<<"-1 -1"<<endl;
	  
      mpri = -10000.0;

      itot++;

    }
    for(int i = 0;i<numZ;i++)
      arr[i].erase(arr[i].begin(), arr[i].end());
    miso1.erase(miso1.begin(),miso1.end());
    aiso1.erase(aiso1.begin(),aiso1.end());
    miso2.erase(miso2.begin(),miso2.end());
    aiso2.erase(aiso2.begin(),aiso2.end());
    tdel_iso.erase(tdel_iso.begin(),tdel_iso.end());
    sma_iso.erase(sma_iso.begin(),sma_iso.end());
    
  }
  out.close();
  out2.close();


  hout.close();
  
  outname = "Catalogue.txt";
  out.open(outname.c_str(),ios::app);
  out2.open("Catalogue_clean.txt");      
  ofstream out3;
  out3.open("Catalogue_multiple_dyn.txt");

  
  hout.open("Larger_than_tH.txt",ios::app);
  
  double semi_ej,semi_gw;


  for(int k=0;k<numZ;k++)  cout<<"      -- Misint. BHs (isolated only) "<<met[k]<<" "<<mis[k]<<endl;

  ofstream clout;
  clout.open("dynsrc_vescape.txt");

  cout<<"     DYNAMICAL BINARIES    "<<endl;
  REC = "RecNO";
  int miso = 0;
  int nrecy= 0;
  double cpri ;

  double *single_bh;
  single_bh = new double [4];
  single_bh[0] = 0.0;
  single_bh[1] = -1.0;
  single_bh[2] = 0.0;
  single_bh[3] = 0.0;
  
  for(int k=0;k<numZ;k++)Npar[k] = 0;
  for(int i=Niso;i<Nsrc;i++){
    for(int k=0;k<numZ;k++){
      Z0 = Zmin + k*dZ;
      Z1 = Zmin + (k+1)*dZ;
      if(Z[i]<pow(10.,Z1) && Z[i]>=pow(10.,Z0))
	Npar[k] += 1;
    }
   
  }

  if(MonoZ=="yes"){
    int idN = -1;
    for(int k=0;k<numZ;k++)
      if(Npar[k] != 0){
	idN = k;
	break;
      }
    for(int k=0;k<numZ;k++)      
      if(Z[Niso+1] == met[k]){
	Npar[k] = Npar[idN];
	Npar[idN] = 0;
	break;
      }         
  }
	
  double total_sources = 0;
  /*
  for(int k=0;k<numbin_Z;k++){
    total_sources += Npar[k];
    cout<<"      -- Metallicity "<<met[k]<<" "<<Npar[k]<<" "<<endl;
  }
  */

  
  double fbin_rnd;
    
  total_sources = 0.0;
  int numme = numZ;
  if(Nyc + Ngc + Nnc == 0)
    numme = 1;    


  int nreal = 0;
  
  if(MonoZ == "yes")
      Zmin = log10(0.0002);
    
  for(int k=0;k<numme;k++){
	auto start = std::chrono::high_resolution_clock::now();
    cout<<"N. of dyn binaries with Z = "<<met[k]<<" "<<Npar[k]<<endl;
    if(Npar[k] == 0){
      continue;
    }
    // qui dobbiamo aprire single vecchi e nuovi...
    vector<double> zams_mix,remn_mix;
    vector<double> zams_sin,remn_sin;
    vector<double> tdel_mix,tdel_sin;
    vector<double> kick_mix,kick_sin;
    
    stringstream s;
    s<<met[k];
    string lol = s.str();
    char *buf;
    buf = new char [lol.length()+1];
    strcpy(buf, lol.c_str());

    stringstream sss;
    for(int j=2;j<lol.length();j++){
      sss<<buf[j];
    }

    double maximus = -1.0;
    double minimus = 1.E30;
    double saximus = -1.0;
    double sinimus = 1.E30;
    double maximus_mix = -1.0;
    double minimus_mix = 1.E30;
    double saximus_mix = -1.0;
    double sinimus_mix = 1.E30;
    
    string filepath;


    double vnat_dx[50];
    double nnat_dx[50];
    double nmix_dx[50];
    double dvnat = pow(5000. / 0.1, 1./50.);
    vnat_dx[0] = 0.1;
    nnat_dx[0] = 0.0;
    nmix_dx[0] = 0.0;
    for(int im=1;im<50;im++){
      vnat_dx[im] = 0.1 * pow(dvnat,im);
      nnat_dx[im] = 0.0;
      nmix_dx[im] = 0.0;
      //cout<<vnat_dx[im-1]<<" "<<vnat_dx[im]<<endl;
    }
    //exit(0);
    
    //DIVORCED BHs
    filepath = path + pathse + pathsp + "/spectrum" + sss.str() + ".txt";
    //cout<<"      -- -- -- opening "<<filepath<<endl;
    in.open(filepath.c_str());
    if(!in.is_open()){
      cout<<"Data file "<<filepath<<" not found"<<endl;
      exit(0);
    }
    do{
      double par[spar];
      for(int jj=0;jj<spar;jj++)
	in>>par[jj];
		       
      zams_mix.push_back(par[0]);
      remn_mix.push_back(par[1]);
      tdel_mix.push_back(1.E6 * par[3]);
      kick_mix.push_back(par[4]);
      
      if(par[1] > maximus_mix && par[1] < 200.0)
	maximus_mix = par[1];
      if(par[1] < minimus_mix)
	minimus_mix = par[1];

      if(par[0] > saximus_mix)
	saximus_mix = par[0];
      if(par[0] < sinimus_mix)
	sinimus_mix = par[0];

      
    }while(!in.eof());
    in.close();


    
    //SINGLE BHs
    filepath = path + pathse + pathsp + "/spectrum_cleaned" + sss.str() + ".txt";
    //cout<<"      -- -- -- opening "<<filepath<<endl;
    in.open(filepath.c_str());
    if(!in.is_open()){
      cout<<"Data file "<<filepath<<" not found"<<endl;
      exit(0);
    }
    do{
      double par[spar];
      for(int jj=0;jj<spar;jj++)
	in>>par[jj];

      /*for(int jj=0;jj<spar;jj++)
	cout<<par[jj]<<" ";
	cout<<endl;*/     

      zams_sin.push_back(par[0]);
      remn_sin.push_back(par[1]);
      tdel_sin.push_back(1.E6 * par[3]);
      kick_sin.push_back(par[4]);
      if(par[1] > maximus)
	maximus = par[1];
      if(par[1] < minimus)
	minimus = par[1];
      if(par[0] > saximus)
	saximus = par[0];
      if(par[0] < sinimus)
	sinimus = par[0];
      
    }while(!in.eof());
    in.close();



    
    double msdx[bin_st],mssx[bin_st];
    double mbdx[nsize],mbsx[nsize];
    double tbdx[nsize],tbsx[nsize];
    double vbdx[nsize],vbsx[nsize];

    double mbhmix[bin_st][nsize];
    double tbhmix[nsize][tsize];
    double vbhmix[nsize][vsize];
       
    //  cout<<"Filling the matrix"<<endl;
    double dm = (saximus_mix - sinimus_mix)/(double)bin_st;
    double db = (maximus_mix - minimus_mix)/(double)nsize;
    double dtb= 500.E6  / (double) tsize;
    double dvb= 1.E3 / (double) vsize;
    
    for(int jj=0;jj<bin_st;jj++){
      mssx[jj] = sinimus_mix + jj*dm;
      msdx[jj] = sinimus_mix +(jj+1)*dm;
    }
    for(int jj=0;jj<nsize;jj++){
      mbsx[jj] = minimus_mix + jj*db;
      mbdx[jj] = minimus_mix +(jj+1)*db;
    }
    for(int jj=0;jj<tsize;jj++){
      tbsx[jj] = jj*dtb;
      tbdx[jj] = (jj+1.)*dtb;
    }
    for(int jj=0;jj<vsize;jj++){
      vbsx[jj] = jj*dvb;
      vbdx[jj] = (jj+1.)*dvb;
    }
    
    //cout<<"Initialize matrix"<<endl;
    for(int jj = 0; jj < bin_st ; jj ++){
      for(int ii=0;ii<nsize;ii++)
	mbhmix[jj][ii] = 0.0;
      for(int ii=0;ii<tsize;ii++)
	tbhmix[jj][ii] = 0.0;
      for(int ii=0;ii<vsize;ii++)
	vbhmix[jj][ii] = 0.0;      
    }


    for(int ii=0;ii<zams_mix.size();ii++){      
      for(int jj=0;jj<bin_st;jj++){
	if(zams_mix[ii] > mssx[jj] && zams_mix[ii] <= msdx[jj]){	  
	  for(int kk=0;kk<nsize;kk++){
	    if(remn_mix[ii] > mbsx[kk] && remn_mix[ii] <= mbdx[kk]){
	      mbhmix[jj][kk] += 1;	     
	    }
	  }
	  break;
	}
      }
      for(int jj=0;jj<nsize;jj++){
	if(remn_mix[ii] > mbsx[jj] && remn_mix[ii] <= mbdx[jj]){	  
	  for(int kk=0;kk<tsize;kk++){
	    if(tdel_mix[ii] > tbsx[kk] && tdel_mix[ii] <= tbdx[kk]){
	      tbhmix[jj][kk] += 1;	      
	    }	      
	  }
	  
	  for(int kk=0;kk<vsize;kk++){
	    if(kick_mix[ii] > vbsx[kk] && kick_mix[ii] <= vbdx[kk]){
	      vbhmix[jj][kk] += 1;	 
	    }	      
	  }
	  
	  break;
	}
      }
          
    }

    /*for(int i=0;i<1000;i++)
      func.singBHt_new(zams_sin, remn_sin, tdel_sin, kick_sin, obslope, mslope, single_bh, saximus,sinimus,maximus,minimus, 1.E30);
      exit(0);*/
    
    double MSLP = mslope;

    for(int i=0;i<1000;i++){
      double kikki;
      func.singBHt_new(zams_sin, remn_sin, tdel_sin, kick_sin, obslope, 0.0, single_bh, saximus,sinimus,maximus,minimus, 1.E30);
      kikki = single_bh[2];
      for(int im=0;im<50;im++){
	if(kikki < vnat_dx[im]){
	  nnat_dx[im] += 1.0;	  
	}	
      }
      singBHt_mix(mssx, msdx, mbsx, mbdx, tbsx, tbdx, vbsx, vbdx, mbhmix, tbhmix, vbhmix, 0.0, single_bh, saximus_mix, sinimus_mix, maximus_mix, minimus_mix, 1.E30);
      kikki = single_bh[2];
      for(int im=0;im<50;im++){
	if(kikki < vnat_dx[im]){
	  nmix_dx[im] += 1.0;	  
	}	
      }
      
      
    }
    stringstream zstring;
    zstring << metdyn[k];
    string fzname = "retention_Z"+zstring.str()+".dat";
    ofstream escout;
    escout.open(fzname.c_str());
    for(int im=0;im<50;im++){
      escout<<vnat_dx[im]<<" "<<double(nnat_dx[im])*1./1000.0 <<" "<<double(nmix_dx[im])*1./1000.0 <<endl;
    }
    escout.close();


    
    
    ndx = 1000.;
    align="dynamical";
    
    cout<<"      -- Metallicity "<<met[k]<<" "<<Npar[k]<<" "<<minimus<<" "<<maximus<<endl;        
    int npar_runtime = 0;
    int nH = 0;

    
    for(int i=Niso; i<Nsrc; i++){
      
      string cluster_stat="none";
	
      
      Zi[i] = Z[i];
      
      Z0 = Zmin + k*dZ;
      Z1 = Zmin + (k+1)*dZ;

      //cout<<pow(10.,Z0)<<" "<<pow(10.,Z1)<<" "<<Z[i]<<" "<<Z[i] - pow(10.,Z0) << endl;

      if( (Z[i] < pow(10.,Z1) && Z[i] >= pow(10.,Z0)) || (MonoZ == "yes" && Z[i] == met[k])) {

	//cout<<"ENTRATO"<<endl;
	
	npar_runtime ++;
	if(i>=Niso     && i<Niso+Nyc)     cluster = "young";
	if(i>=Niso+Nyc && i<Niso+Nyc+Ngc) cluster = "globular";
	if(i>=Niso+Nyc+Ngc)               cluster = "nuclear";

	
	//NB: the primary shouldn't have zero mass!!!	
	// CLUSTER'S PROPERTIES //

	double vthre = 0.0;
	double nbhs  = 0.0;
	
	vector<double> mx,my;
	vector<double> rx,ry;
	
	int ipre,ipos;
	string pcluster = cluster;
	string label = "retained";
	
	if(cluster == "young"){
	  //vthre = 3.0;
	  mx = YCmx;
	  my = YCmy;
	  rx = YCrx;
	  ry = YCry;
	}
	if(cluster == "globular"){
	  //vthre = 15.0;
	  mx = GCmx;
	  my = GCmy;
	  rx = GCrx;
	  ry = GCry;
	}
	if(cluster == "nuclear"){
	  //vthre = vmaxNC;
	  mx = NCmx;
	  my = NCmy;
	  rx = NCrx;
	  ry = NCry;
	}

	int kmeta;
	double esc,mint,rint,lsig,sig_clu,rho_clu,mstar,rho_cubicpc;
	double vthre_in, sig_clu_in,rho_clu_in,rho_cubicpc_in;
	double metal, zita, logL, rhalf, mhalf, fb;
	double tSNe, mper, trelax, trelax0, tcc, tdf, semihard, kappa, semi, ecc, t12, acrit, tbbh, tmer, mu_red, sma;
	double t3bb, t12capt, tbbhform;
	double time,nsafe_glob,nsafe,mass_ratio, mixer, nhigen;
	string stri_mrat;

	double mclcorr = 1.0;
	double rclcorr = 1.0;

	zita = 1.0;
	logL = 15.;

	//Here we need to set a specific binary fraction per cluster type (and perhaps time evolution?)
	//fb   = 0.01; 	
	if(cluster == "young")
	  fb = 0.4;
	else if(cluster == "globular")
	  fb = 0.1;
	else if(cluster == "nuclear")
	  fb = 0.01;
	//
	

	mper = 15.0;
	
	mixer = func.rnd();

	double vlimiting = 2.;


	if(mixer > mixing){
	  if(met[k] == 0.02)
	    vlimiting = 10.;
	}
	
	

	int nsafe_cal = 0;
	//do{ removed this to take into account also ICs that do not produce merging BHs
	
	  nsafe_cal ++;
	  //SELECTCLUSTER PROPERTIES
	  mint = 0;
	  do{	    
	    
	    esc = func.rnd();
	    ipre=mx.size();
	    ipos=0;
	    for(int escid=0;escid<mx.size();escid++){
	      if(esc>my[escid])
		ipre = escid;
	      if(esc<my[escid]){
		ipos = escid;
		break;
	      }	  	    
	    }
	    
	    if(ipos==0){
	      ipre = 0;
	      ipos = 1;
	    }
	    if(ipre==mx.size()-1){
	      ipre = mx.size()-2;
	      ipos = mx.size()-1;
	    }
	    
	    
	    mint = (mx[ipre] - mx[ipos])/(my[ipre] - my[ipos])*(esc - my[ipre])+mx[ipre];
	    if(cluster_test == "yes"){
	      
	      mint = 4.0 + 3.0*func.rnd();

	      if(cluster_test_env == "NC")
		mint = 6.0+4.0*func.rnd();

	      if(cluster_test_env == "YC")
		mint = 3.0 + 2.0*func.rnd();
	      
	      
	    }
	   

	      
	    do{
	      esc = func.rnd();
	      ipre=rx.size();
	      ipos=rx[0];
	      for(int escid=0;escid<rx.size();escid++){
		if(esc>ry[escid])
		  ipre = escid;
		if(esc<ry[escid]){
		  ipos = escid;
		  break;
		}	  	    
	      }
	      
	      
	      rint = (rx[ipre] - rx[ipos])/(ry[ipre] - ry[ipos])*(esc - ry[ipre])+rx[ipre];

	      double mclu_int = pow(10., mint);
	      if(TagR == "Marks12" && cluster == "young"){

		/*
		double mu_mrk12 = func.rnd();
		double smrk12 = 2.0;
		mint = log10(pow(mu_mrk12 * pow(1.E6,1.-smrk12) + (1.-mu_mrk12)*pow(3.E2,1.-smrk12),1./(1.-smrk12)));

		mclu_int = pow(10.,mint);
		*/
		
		double alphaMK12 = 0.1 * (1. + (-0.04+0.11*func.rnd()));		
		double betaMK12 = 0.13 * (1. + 0.04*(1.-2.*func.rnd()));
		rint = log10(alphaMK12 * pow(mclu_int,betaMK12));

	      }	     
	      else if(TagR == "Rantala24" && cluster == "young"){

		double alphaBG = 2.365 * (1. + 0.106*(1.-2.*func.rnd()));
		double betaBG = 0.180 * (1. + 0.028 * (1.-2.*func.rnd()));
		rint = log10(0.125/1.3 * alphaBG * pow(mclu_int/1.E4,betaBG));

	      }
	      else if(TagR == "Mapelli20"){

		double mean, sigma;
		sigma = log(pow(10.,0.4));
		if(cluster=="young")
		  mean = log(pow(10.,4.3));
		else if(cluster == "globular")
		  mean = log(pow(10.,5.6));
		else if(cluster == "nuclear")
		  mean = log(pow(10.,6.18));

		//mint = func.LOGSMP(mean, sigma);
		mint = func.LogGaussian(mean, sigma);		
		
		mint = log10(mint);
		double rhoint = -1;
		
		if(cluster=="young")
		  mean = log(pow(10.,3.3));
		else if(cluster == "globular")
		  mean = log(pow(10.,3.7));
		else if(cluster == "nuclear")
		  mean = log(pow(10.,5.0));


		//rhoint = func.LOGSMP(mean, sigma);
		rhoint = func.LogGaussian(mean, sigma);
		rhoint = log10(rhoint);
		
		rint = (mint - rhoint)/3.0;

	      }
	      else if(TagR == "AS20" && cluster == "nuclear"){
		//Select if it's late type (50%) or early type (50%)
		//calculate rint from mint using a 3sigma dispersion
		//assign rint

		double Plate = func.rnd();
		double c1_nc, c2_nc, a_nc, b_nc, e_nc;
		if(Plate > 0.5){
		  c1_nc = 3.31;
		  c2_nc = 3.6E6;
		  a_nc  = 0.321;
		  b_nc  = -0.011;
		  e_nc =  5.*0.133;
		}		  
		else{
		  c1_nc = 6.27;
		  c2_nc = 1.95E6;
		  a_nc  = 0.347;
		  b_nc  = -0.024;
		  e_nc  = 5.*0.131;
		}

		rint = func.geo16Rnc(mint, c1_nc, c2_nc, a_nc, b_nc, e_nc);

	      }
	      else{
		if(TagR != "AS20" && cluster=="young"){
		  cout<<"Please choose a recipe for initial cluster radii"<<endl;
		  exit(0);
		}
	      }
	      
	      
	      
	      if(TagR != "AS20" && (cluster=="nuclear" || cluster_test == "yes")){
		double lnrhalf;
		
		double A_r = 0.222887 ; // +/- 0.03415      (15.32%)
		double B_r = -0.909231; // +/- 0.2196       (24.15%)
		
		if(cluster_test == "yes" && cluster_test_env != "NC")
		  lnrhalf = -1.5 + 2.5*func.rnd();
		else if(cluster == "nuclear" || cluster_test_env == "NC"){
		  lnrhalf = A_r * (1.+ 3.*0.03415*(1.-2.*func.rnd())) * mint + B_r * (1. + 5.*0.2196 * (1.-2.*func.rnd()));

		}
		
		rint = lnrhalf ;
		
	      }
	    }while( pow(10.,rint) < 0.0 ); //10.*0.8*pow(pow(10.,mint)/1000.,0.2));
	    
	    
	    
	    vthre = func.vescape(pow(10.,rint),pow(10.,mint),pcluster);
	  }while(vthre < vlimiting);// || mint < 3.0 || pow(10.,mint)/pow(pow(10.,rint),3.) < 100.);
	  	  
	  //the limiting values above lead to 20-50% loss of sources and ensure that the time remain t < 13.5 Gyr//
	  //particularly the M/R^3 < 100 Msun/pc^3 region has no sources with t < 13.5 Gyr with single BHs  !!//
	  
	  lsig = 0.5 * (-1.14 + log10(6.67E-11*1.99E30/3.08E16 * pow(10.,mint)/pow(10.,rint) ));			 
	  sig_clu = 0.001*pow(10.,lsig);
	  
	  rho_clu = 3. *pow(10.,mint) / ( 4. * M_PI * pow(pow(10.,rint),3.));//(3.-2.*func.rnd())*pow(10.,mint) / ( 4. * M_PI * pow(pow(10.,rint),3.));
	  
	  mstar = 1.0;
	  rho_cubicpc = rho_clu / mstar;

	  //nbhs = 0.01*pow(10.,mint)/30.;
	  double alpha_AS = 0.903;
	  double beta_AS  = -0.79; 
	  double alpha1_AS= 1.14;
	  double beta1_AS = -0.62;
	  double Mcluster = pow(10.,mint);

	  
	  double fblack   = 0.083; 
	  
	  
	  rhalf = pow(10.,rint);
	  mhalf = pow(10.,mint);	  	  

	  
	  trelax =  0.78E9 /log(0.11 * mhalf) * pow(mhalf/1.E5, 0.5) * pow(rhalf, 1.5); //4.2E9 * (15./logL) * pow(rhalf/4.0,1.5) * sqrt(mhalf/1.E7) ;
	  trelax0 = trelax;

	  tcc = 0.138*mhalf/(150. * log(0.11*mhalf/150.))*sqrt(pow(rhalf*3.08E16,3.)/(6.67E-11*1.99E30*mhalf))/(365.*24.*3600.);

	  if(tcc < 0. || tcc > trelax)
	    tcc = 0.2 * trelax; //Portegies-Zwart & McMillan 2002/2004
	   
	  
	  vthre_in = vthre;
	  sig_clu_in = sig_clu;
	  rho_clu_in = rho_clu;
	  rho_cubicpc_in = rho_cubicpc;
	  
	  // HERE WE HAVE THE CLUSTER INITIAL PROPERTIES
	  /*

	    here we can plan the cluster evolution, considering the following examples:
	    mcl_time = pow(10.,mint) * exp(ttt / (0.45*trelax));
	    rcl_time = pow(10.,rint) * (1. + log10(1. + ttt/(0.45*trelax))); //REVISE THIS

	    In this regards, we should go as follows: 
	    - calculate bh birth time and df time, use maximus to calculate mcl_time and rcl_time
	    - calculate t3bb and tbs with the updated mcl_time and rcl_time, use the min(t3bb,tbs) to update mcl_time, rcl_time
	    - calculate t12 with updates and update mcl_time, rcl_time

	   */

	  // Cluster scale radius (from Dehnen) //
	  double g_cl;
	  if(cluster == "young")
	    g_cl = 1.0;
	  else if(cluster == "globular")
	    g_cl = 1.5;
	  else if(cluster == "nuclear")
	    g_cl = 1.9;
	  
	  double a_cl = rhalf*(pow(2.,1./(3.-g_cl))-1);

	  // DOUBLE CHECK THE FOLLOWING //
	  
	  // Max radius from which BHs can spiral-in over a Hubble time via DF //	  
	  double radius = rhalf * pow( (Hubble - tfor[i]*1.E9) / (0.42E9 * (10.*mstar/(mpri+msec)) * (trelax / 4.2E9)) , 1./1.74);

	  // Fraction of mass enclosed within the infall radius above //
	  double fencl = (radius / (radius + a_cl),3.-g_cl) * (1. + 0.2*(1.-2.*func.rnd()));

	  // retention fraction freten //
	  double freten = 0.5 * (1. + 0.3*(1.-2.*func.rnd()));

	  // Fraction number of BHs in a power-law IMF between 0.08 and 150 Msun//
	  double fraBH = 0.0008 * (1. + 0.1*(1.-2.*func.rnd()));
	  
	  // This depends on the number fraction of BHs in the cluster, we're also assuming mint == N_* 
	  nbhs = fraBH * pow(10.,mint) * freten * fencl;


	  
	  if(nbhs < 2)
	    nbhs = 2;
	  
	  
	  if(nbhs < 1){
	    cout<<"Warning -- no BHs in the centre!? "<<nbhs<<" "<<pow(10.,mint)<<" "<<pow(10.,rint)<<" "<<vthre<<" "<<pcluster<<endl;
	    exit(0);
	  }
	  
	  
	  stri_mrat = MRATIO;
	  nsafe_glob = 0;
	  nsafe = 0;	 
	  nhigen = 0;
	  
	  if(mixer > mixing){		
	    mpri = -1;
	    kpri = 1.E30;
	    do{
	      func.singBHt_new(zams_sin, remn_sin, tdel_sin, kick_sin, obslope, mslope, single_bh, saximus,sinimus,maximus,minimus, vthre);
	      mpri = single_bh[0];	 
	      tpri = single_bh[1];
	      kpri = single_bh[2];
	      if(mpri > 0.0 && kpri > vthre)
		break;

	      nsafe ++;
	      
	    }while(mpri <= 0.0 || kpri > vthre);
	  }
	  else{	  
	    MSLP = mslope;
	    mpri = -1;
	    kpri = 1.E30;
	    int nsafe = 0;
	    do{
	      singBHt_mix(mssx, msdx, mbsx, mbdx, tbsx, tbdx, vbsx, vbdx, mbhmix, tbhmix, vbhmix, MSLP, single_bh, saximus_mix, sinimus_mix, maximus_mix, minimus_mix, vthre);
	      mpri = single_bh[0];	  	  
	      tpri = single_bh[1];
	      kpri = single_bh[2];	      
	      if(nsafe > 1000)
		break;
	      
	      nsafe ++;
	    }while(mpri <= 0.0 || kpri > vthre);
	  }	  
	  apri   = func.spin(mpri,dynaS);	

	  if(nsafe == 1000)
	    cout<<"Wrong BH"<<endl;
	  
	  nsafe_glob += nsafe;
	  
	  mass_ratio = -1;
	  nsafe = 0;       
	  if(stri_mrat != "nouniform"){
	    do{
	      mass_ratio = func.mratio(mpri, MRATIO_SLOPE, stri_mrat);
	      msec = mpri * mass_ratio;
	      for(int iii=0; iii<remn_sin.size()-1;iii++){
		if(msec > remn_sin[iii] && msec < remn_sin[iii+1]){
		  tsec = 0.5*(tdel_sin[iii] + tdel_sin[iii+1]) ;
		}
	      }

	      if(nsafe > 1000){
		cout<<"mratio fails"<<endl;
		exit(0);
	      }
	      nsafe += 1;
	    }while(msec < 1. || msec > 500.);
	    
	  }
	  else{
	    msec = -1;
	    ksec = 1.E30;
	    nsafe = 0;
	    if(mixer > mixing){
	      do{
		
		func.singBHt_new(zams_sin, remn_sin, tdel_sin, kick_sin, obslope, mslope, single_bh,saximus,sinimus,maximus,minimus,vthre);	  
		msec = single_bh[0];	 
		tsec = single_bh[1];
		ksec = single_bh[2];	    

		if(msec > minimus && ksec < vthre)
		  break;
		
		nsafe ++;
		if(single_bh[3] == 1)
		  break;
		
		else if(nsafe > 500){
		  cout<<npar_runtime<<" "<<i<<" "<<Z[i]<<" "<<msec<<" "<<tsec<<" "<<ksec<<endl;
		}
		
		if(nsafe > 1000){
		  cout<<"Something wrong " <<msec<<" "<<tsec<<" "<<ksec<<" "<<minimus<<" "<<vthre<<endl;
		  exit(0);
		}
		
	      }while(msec <= minimus || ksec > vthre);
	      
	      if(msec < minimus){
		cout<<"Second BH mass below mmin = "<<minimus<<" "<<msec<<endl;
	      }
	      
	    }
	    else{
	      double MSLP = mslope;
	      msec = -1;
	      ksec = 1.E30;
	      do{
		singBHt_mix(mssx, msdx, mbsx, mbdx, tbsx, tbdx, vbsx, vbdx, mbhmix, tbhmix, vbhmix, MSLP, single_bh, saximus_mix, sinimus_mix, maximus_mix, minimus_mix, vthre);
		msec = single_bh[0];	  	  
		tsec = single_bh[1];
		ksec = single_bh[2];	    
		if(nsafe > 1000)
		  break;
		nsafe ++;

	      }while(msec <= 0.0 || ksec > vthre);
	      
	    }

	    if(msec == 0){
	      time = 1.e30;
	      break;
	    }
	    
	    if(highgen == "yes"){
	      //This section needs to be improved, especially for what concern the probability of repeated mergers
	      double p_BPOP;
	      if(cluster == "globular")
		p_BPOP = 295. / 1.E4; //this is valid for GCs with Z = 0.0002 (would it change with different metallicity?)
	      else if(cluster == "nuclear")
		p_BPOP = 2740./ 1.E4;
	      else if(cluster == "young")
		p_BPOP = 0.0;
	    
	      double p_hgen = func.rnd();
	      double mav = 30.;
	      if(p_hgen < p_BPOP){
		double coef = 0.2 + 0.8*func.rnd();
		msec = mav * (1. + msec/mav * (1.+mav/6.0));	       
	      }
	      nhigen += int(msec / mav);

	      asec = 0.7; // Need adjustments
	      
	    }


	    if(msec > mpri){
	      double dum1;
	      dum1 = msec;
	      msec = mpri;
	      mpri = dum1;
	      dum1 = tsec;
	      tsec = tpri;
	      tpri = dum1;
	    }
	  }

	  //This section serves for the binary component masses --- need to be added also in the hierarchical merger chain
	  double prob_ugp = func.rnd();
	  double prob_fgp = func.rnd();
	  string UP = upgtp;
	  double fUP= fupgp;
			     
	  if(UP == "dicarlo")
	    fUP = 0.01 * (1. + 5.797 * exp(Z[i] / 0.0002 * log(5./5.797)));
	  
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
	      double mex_gp = m1_gp;
	      m1_gp = m2_gp;
	      m2_gp = m1_gp;
	    }
	      
	    // We need to implement something for the spins too ...
	    if(prob_ugp >= pbelow && prob_ugp < pbelow + pbelup){
	      if(mpri < mass_gap)
		mpri = m1_gp;		  	      
	      else if (msec < mass_gap)
		msec = m2_gp;
	      else{
		cout<<"Both objects above the gap!"<<endl; //This may happen in multiple mergers?
	      }
	      
	    }	      
	    else if(prob_ugp > pbelow + pbelup){
	      if(mpri < mass_gap)
		mpri = m1_gp;
	      if(msec < mass_gap)
		msec = m2_gp;
	    }		
	    
	  }

	  
	  nsafe_glob += nsafe;
	  double dmy = func.rnd();
	  if(bhseed == "bifrost" && dmy < f_seed && mint > log10(5.3E3) && pow(10., mint - 3.*rint) > 1.E5 && Z[i] < bifZ){ //min mass to form a seed of at least 150 Msun

	    double mprimin = min(150., pow(pow(10.,mint), 0.53)*pow(10., -0.23));
	    double mprimax = min(2.E4, 0.02 * pow(10.,mint));
	      
	    double mprit = mprimin + (mprimax - mprimin)*func.rnd();

	    if(mprit > bhpisn)
	      mpri = 0.9*mprit;

	    
	    	    
	  }
	  else if(bhseed == "vms"){	    
	    if(dmy < f_seed){
	      double A_vms, B_vms;

	      double provms = func.rnd();
	      if(provms < 0.95){
		A_vms = 309.;
		B_vms = 0.015;
	      }
	      else if(provms > 0.95 && provms < 0.99){
		A_vms = 648.;
		B_vms = 0.05;
	      }
	      else{
		A_vms = 1284.;
		B_vms = 0.1;
	      }
	      
	      double rho_vms = pow(10.,mint - 3.*rint);
	      mpri = A_vms * pow(rho_vms / 1.E9, B_vms) * pow(10., -0.05 + 0.1 * func.rnd());
	      
	    }
	  }
	  else if(bhseed == "yes"){
	    if(dmy < f_seed){
	      double mimmi = func.rnd();
	      mpri = pow(mimmi * pow(maxseed, 1.-seedslope) + (1.-mimmi) * pow(minseed,1.-seedslope) , 1./(1.-seedslope));	      
	    }	    
	  }
	  else if(bhseed == "density"){
	    double rho_vms = pow(10.,mint - 3.*rint);
	    if(dmy < f_seed && rho_vms > 3.E5){
		double mimmi = func.rnd();
		mpri = pow(mimmi * pow(maxseed, 1.-seedslope) + (1.-mimmi) * pow(minseed,1.-seedslope) , 1./(1.-seedslope));	      
	    }
	    	  
	  }



	  
	  //Also, we should separate between the density stuff and the upper mass-gap stuff, that is more related to the binary fraction indeed

	  

	  
	  asec   = func.spin(msec,dynaS);
	  
	  
	  //FIRST MERGER
	  Spinning[0] = 0.0;
	  Spinning[1] = 0.0;
	  Spinning[2] = 0.0;
	  Spinning[3] = 0.0;
	  Spinning[4] = 0.0;
	  Spinning[5] = 0.0;
	  Spinning[6] = 0.0;

	  if(mpri>0.0 && msec>0.0)
	    func.SREM2(ndx, apri, asec, mpri, msec, align, Spinning);	  
	  Srem[i] = Spinning[0];
	  Xrem[i] = Spinning[1];
	  Mrem[i] = Spinning[2];
	  Krem[i] = Spinning[3];
	  Cosa[i] = Spinning[4];
	  Cosb[i] = Spinning[5];
	  Cosg[i] = Spinning[6];

	  
	  
	  time = tfor[i];
	  tSNe = max(tpri,tsec);

	  if(bhseed != "no")
	    tSNe = max(5.0e6, tcc);
	    
	  
	  rhalf = pow(10.,rint);
	  mhalf = pow(10.,mint);
	   
	  trelax = 0.78E9 /log(0.11 * mhalf) * pow(mhalf/1.E5, 0.5) * pow(rhalf, 1.5); //4.2E9 * (15./logL) * pow(rhalf/4.0,1.5) * sqrt(mhalf/1.E7) ;		  
	  //tDF to sink	
	  tdf = 0.42E9 * (10.*mstar/(mpri+msec)) * (trelax / 4.2E9);	  	  
	  
	  if(tdf > tSNe)
	    time += tdf;	  
	  else
	    time += tSNe;

	  
	  //if(CLevo == "yes"){
	    
	  mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  	    
	  sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	  rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	  rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	  
	  vthre = vthre_in * sqrt(mclcorr/rclcorr);
	  
	  semihard = 1./pow(sig_clu/30., 2.); //2.*6.67E-11*1.99E30/(1.E6*1.5E11) * (mpri+msec) / (sig_clu*sig_clu);
	  
	  t12capt = 3.E9 * (0.01 / fb) * (1.E6 / rho_cubicpc) * (sig_clu / 30.) * 10./(mpri + msec + mper) * 1./semihard;
	  t3bb = 4.E9 * pow(1.E6 / rho_clu,2.) * pow(sig_clu/(zita * 30.),9.) * pow(mstar / 30. * 10., 9./2.) * pow(10./30.,-5.);
	  
	  fbin_rnd = func.rnd();
	  if(fbin_rnd < fbin)
	    t3bb = 0.0;
	  
	  if(t3bb < 0){
	    cout<<t3bb<<" "<<rho_clu<<" "<<sig_clu<<" "<<mstar<<" "<<zita<<endl;
	    exit(0);
	  }
	  tbbhform = min(t12capt,t3bb);
	  tbbhform *= func.rndgen(1.0, 0.1);
	  
	  if(mclcorr * mhalf < 1.E1)
	    tbbhform = 1.E12;
	  
	    
	  if(time + tbbhform > tfor[i]+tcc && time < tfor[i]+tcc){
	    time = tfor[i]+tcc;
	    
	    
	    mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	    rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	    	    
	    sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	    rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	    rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	    
	    vthre = vthre_in * sqrt(mclcorr/rclcorr);

	    if(mclcorr == 0.0 && rclcorr == 0.0)
	      time += 1.E12;

	  }	    
	    
	    
	  //}

	  
	 	  
	  //Binary formation and BH captures

	  
	  semihard = 1./pow(sig_clu/30., 2.); //2.*6.67E-11*1.99E30/(1.E6*1.5E11) * (mpri+msec) / (sig_clu*sig_clu);
	  
	  t12capt = 3.E9 * (0.01 / fb) * (1.E6 / rho_cubicpc) * (sig_clu / 30.) * 10./(mpri + msec + mper) * 1./semihard;
	  t3bb = 4.E9 * pow(1.E6 / rho_clu,2.) * pow(sig_clu/(zita * 30.),9.) * pow(mstar / 30. * 10., 9./2.) * pow(10./30.,-5.);

	  fbin_rnd = func.rnd();
	  if(fbin_rnd < fbin)
	    t3bb = 0.0;
	  
	  if(t3bb < 0){
	    cout<<t3bb<<" "<<rho_clu<<" "<<sig_clu<<" "<<mstar<<" "<<zita<<endl;
	    exit(0);
	  }
	  tbbhform = min(t12capt,t3bb);
	  tbbhform *= func.rndgen(1.0, 0.1);
	  
	  if(mclcorr * mhalf < 1.E1)
	    tbbhform = 1.E12;
	  
	  
	  time += tbbhform;	  

	  
	  
	  //if(CLevo == "yes"){
	    
	  mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  
	  sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	  rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	  rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	  
	  vthre = vthre_in * sqrt(mclcorr/rclcorr);
	  
	  t12 = 3.E8 / zita * (0.01/fb) * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(mstar / (mpri+msec)) * (30./(mpri+msec+mper)) * (1./semihard);
	  t12 *= func.rndgen(1.0, 0.1);
	    
	  if(mclcorr * mhalf < 1.E1)
	    t12 = 1.E12;
	    
	  
	  if(time < tfor[i]+tcc && time + t12 > tfor[i]+tcc){
	    time = tfor[i]+tcc;
	    
	    mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	    rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	    
	    sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	    rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	    rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	    
	    vthre = vthre_in * sqrt(mclcorr/rclcorr);
	  }
	  //}
	  

	  //Time for BHs to pair --> assume it's already hard //see Antonini&Rasio2016
	  
	  t12 = 3.E8 / zita * (0.01/fb) * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(mstar / (mpri+msec)) * (30./(mpri+msec+mper)) * (1./semihard);
	  t12 *= func.rndgen(1.0, 0.1);

	  if(mclcorr * mhalf < 1.E1)
	    t12 = 1.E12;

	  time += t12;

	 
	  
	  
	  //if(CLevo == "yes"){

	  mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  
	  sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	  rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	  rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	  
	  
	  vthre = vthre_in * sqrt(mclcorr/rclcorr);
	    //}

	  //Calculate hard binary separation and assign binary sma
	  semihard = 1./pow(sig_clu/30., 2.); 
	  
	  sma = semihard * func.Gaussian(1.,sigma_sma);

	  if(sigma_sma < 0)
	    sma = 0.1 + 0.1*func.rnd();

	  
	  
	  //calculate aej vs agw to merge inside or outside --> inside (tgwpeters), outside (5/q3 * t2-1)	  
	  semi_ej = 0.07 * (mu_red * mper / (mpri + msec + mper) * mper/(mpri+msec)) * pow(vthre/50.,-2.);
	  semi_gw = 0.05 * pow((mpri+msec)/20,3./5.)*pow( (msec/mpri) / pow(1+(msec/mpri),2.) , 1./5.) * pow(sig_clu/30., 1./5.) * pow(1.E6 / rho_clu,1./5.);
	  
	  //acrit = max(semi_ej,semi_gw);
	  //acrit = min(acrit,0.01); removed on 13-03-2024
	  //new recipe for calculations

	  double afin;
	  if(semi_ej > semi_gw)
	    afin = min(sma, semi_ej);
	  else
	    afin = min(sma, semi_gw);
	  
	  acrit = afin;

	  tbbh = 2.E7/zita * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(10. * mstar / (mpri+msec)) * (0.05 / acrit) * (20./(mpri+msec));
	  tbbh *= func.rndgen(1.0, 0.1);	  
	  
	  if(mclcorr * mhalf < 1.E1){
	    tbbh = 1.E12;
	    semi_ej = -2.0;
	    semi_gw = -1.0;
	  }
	  
	  tmer = 5. * (mpri + msec)/mper * tbbh;						 	  

	  mclcorr = func.mevol(time+tmer-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  rclcorr = func.revol(time+tmer-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);	  
	  sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	  rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	  rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);	  	 
	  vthre = vthre_in * sqrt(mclcorr/rclcorr);

	  if(mclcorr*mhalf <= 1.E1)
	    tmer = 1.E12;
	  
	  
	  if((time > 1.E9 || time+tmer > 1.E9) && cluster == "young" && CLfill == "critical"){
	    time = 1.E30;
	    label = "YCdisrupt";
	  }


	  
	  time += tmer;


	  
	  /*if(nsafe_cal > 10000)
	    break;*/ //this was in the do{}while statement, which has been removed, see below

	  
	// }while(time > Hubble); I removed this condition to allow for all clusters that do not produce a merger 

	
	   
	  /* Insert a statistical consideration to check whether binaries with semi_ej > semi_gw can merge before ejection owing to a sufficiently large value of ecc -- May 2024 */
	  // At 0th order, we can get

	  double clight = 3.E8;
	  double Grav   = 6.67E-11;
	  double Msun   = 1.99E30;
	  double Rau    = 1.5E11;	    
	  
	  /* Prescription suggested by F. Antonini, here however overestimate the amount of in-cluster mergers 
	    double elim = sqrt(1 - pow(256./5.*((tbbh+tmer) * (365.*24*3600.) * pow(Grav * Msun, 3.) * (mpri*msec*(mpri+msec)))/(pow(clight,5.) * pow(Rau,4.)  * pow(acrit,4.)), 2./7.));
	  double Pelim = 1. - (elim * elim);
	  double Prnd = func.rnd();
	  if(Prnd < Pelim){
	    label = "inside";
	    double Pecc = Pelim + (1.-Pelim)*func.rnd();
	    ecc = sqrt(1.-Pecc);	    
	  }
	  else{
	    ecc = sqrt(func.rnd());
	    if(semi_ej < semi_gw)
	      label="inside";
	    else{
	      label="ejected";
	      acrit = 2. * afin / (2. + afin/semihard); 
	    }
	  }
	  */
	//Old prescription -- we must ensure that tight binaries merge inside if acrit < semi_ej//
	  if(semi_ej < semi_gw)
	    label="inside";
	  else
	    if(acrit < semi_gw)
	      label="inside";
	    else
	      label="ejected";
	    
	  ecc = sqrt(func.rnd());

	  double teje = tmer;
	
	  tmer = 5./256. * pow(clight,5.) / pow(Grav * Msun, 3.) * pow(Rau,4.)  * pow(acrit,4.)/(mpri*msec*(mpri+msec)) * pow(1.-pow(ecc,2.),3.5) / (1. + 73./24. * ecc*ecc + 37./96.*ecc*ecc*ecc*ecc) * 1./ (365.*24.*3600.*1.E0);

	  if(tmer < teje && label == "ejected")
	    label = "inside";

	  
	  
	/*REMOVING THE FOLLOWING TO TAKE INTO ACCOUNT NON-MERGING BINARIES 
	  if(time + tmer > Hubble && time < Hubble && mclcorr*mhalf > 1.E1){
	  do{
	    ecc = sqrt(func.rnd());
	    tmer = 5./256. * pow(clight,5.) / pow(Grav * Msun, 3.) * pow(Rau,4.)  * pow(acrit,4.)/(mpri*msec*(mpri+msec)) * pow(1.-pow(ecc,2.),3.5) / (1. + 73./24. * ecc*ecc + 37./96.*ecc*ecc*ecc*ecc) * 1./ (365.*24.*3600.*1.E0);

	    //count how many throws needed to have one more merger
	    nsafe_glob ++;
	    
	  }while(time + tmer > Hubble);
	}
	else */
	if(mclcorr*mhalf < 1.E1) 
	  tmer = 1.E12;
	
	time += tmer;

	
	//ADDED PRINTOUT OF BINARY PROPERTIES      
	if(time > Hubble){ //The other condition is not needed if nsafe_cal is just 1 --> || nsafe_cal > 1){
	  hout<<itot<<" "<<pow(10.,mint)<<" "<<pow(10.,rint)<<" "<<mclcorr<<" "<<rclcorr<<" "<<trelax0<<" "<<vthre<<" "<<sig_clu<<" "<<tdf<<" "<<tbbhform<<" "<<t12<<" "<<tmer<<" "<<nsafe_cal+nsafe_glob<<" "<<Z[i]<<" "<<cluster<<" "<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<ecc<<" "<<sma<<" "<<acrit<<" "<<nrecy<<" "<<time<<endl;
	}

	stringstream nan;
	nan<<time;
	
	if(nan.str() == "nan" || nan.str() == "-nan"){
	  cout<<"Step error "<<mpri<<" "<<msec<<" "<<tfor[i]<<" "<<tSNe<<" "<<tdf<<" "<<t12<<" "<<tbbh<<" "<<tmer<<" "<<sig_clu<<" "<<vthre<<" "<<mclcorr<<" "<<rclcorr<<endl;
	  exit(0);
	}

	
	if(time < Hubble)
	  nH ++;
	
	clout<<pow(10.,mint)<<" "<<pow(10.,rint)<<" "<<vthre<<" "<<sig_clu<<" "<<rho_clu<<" "<<pcluster<<" "<<nbhs<<" ";
	clout<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<kpri<<" "<<ksec<<" "<<Mrem[i]<<" "<<Srem[i]<<" "<<Xrem[i]<<" "<<Krem[i]<<" "<<time<<" "<<tdf<<" "<<t12<<" "<<tbbh<<" "<<tmer<<" "<<(double) nH / (double) npar_runtime;
	clout<<endl;	


	//MULTIPLE MERGER CHAIN//
	Spinning[0] = 0.0;
	Spinning[1] = 0.0;	
	Spinning[2] = 0.0;
	Spinning[3] = 0.0;
	Spinning[4] = 0.0;
	Spinning[5] = 0.0;
	Spinning[6] = 0.0;

	if(mpri>0.0 && msec>0.0)
	  func.SREM2(ndx, apri, asec, mpri, msec, align, Spinning);	  

	Srem[i] = Spinning[0];
	Xrem[i] = Spinning[1];
	Mrem[i] = Spinning[2];
	Krem[i] = Spinning[3];
	Cosa[i] = Spinning[4];
	Cosb[i] = Spinning[5];
	Cosg[i] = Spinning[6];

	vpri = Spinning[3];

	double mzero = mpri;
	
	//if(CLevo == "yes"){
	  

	mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	
	sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	
	
	vthre = vthre_in * sqrt(mclcorr/rclcorr);
	
	if(mclcorr * mhalf < 1.E1)
	  cluster_stat="evaporated";
	//}

	double rinfinite = pow(10.,rint) * sqrt(pow(vthre,4.) / pow(vthre*vthre - vpri*vpri,2.) - 1.);       
	//if(CLevo == "yes")
	rinfinite *= rclcorr;

	
	int nrecy = nhigen;
	double trecy = 0.0;
	double tsmbh = 0.0;

	double sig_clu0 = sig_clu;

	
	
	do{

	  if(time > Hubble)
	    break;
	  
	  stringstream nan;
	  nan<<time;
	  
	  if(nan.str() == "nan" || nan.str() == "-nan"){
	    cout<<"Step error mean"<<tfor[i]<<" "<<tSNe<<" "<<tdf<<" "<<t12<<" "<<tbbh<<" "<<tmer<<endl;
	    exit(0);
	  }
	  
	  if(time+trecy > Hubble){
	    label = "highlander";
	    
	    Srem[i] = Spinning[0];
	    Xrem[i] = Spinning[1];
	    Mrem[i] = Spinning[2];
	    Krem[i] = Spinning[3];
	    Cosa[i] = Spinning[4];
	    Cosb[i] = Spinning[5];
	    Cosg[i] = Spinning[6];

	    hout<<itot<<" "<<pow(10.,mint)<<" "<<pow(10.,rint)<<" "<<mclcorr<<" "<<rclcorr<<" "<<trelax0<<" "<<vthre<<" "<<sig_clu<<" "<<tdf<<" "<<tbbhform<<" "<<t12<<" "<<tmer<<" "<<nsafe_cal+nsafe_glob<<" "<<Z[i]<<" "<<cluster<<" "<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<ecc<<" "<<sma<<" "<<acrit<<" "<<nrecy<<" "<<time<<endl;
	  

	    break;
	  }


	  if(Spinning[3] > vthre && label != "ejected"){
	    label = "GWrecoiled";
	    Srem[i] = Spinning[0];
	    Xrem[i] = Spinning[1];
	    Mrem[i] = Spinning[2];
	    Krem[i] = Spinning[3];
	    Cosa[i] = Spinning[4];
	    Cosb[i] = Spinning[5];
	    Cosg[i] = Spinning[6];
	    break;
	  }
	  
	  if(cluster_stat == "evaporated"){
	    if(label != "ejected")
	      label = "inevap";
	    break;
	  }
	  
	  time += trecy;
	  nan.str("");
	  nan<<time;	  
	  if(nan.str() == "nan" || nan.str() == "-nan"){
	    cout<<"Step error Final "<<tfor[i]<<" "<<tSNe<<" "<<tdf<<" "<<t12<<" "<<tbbh<<" "<<tmer<<endl;
	    exit(0);
	  }

	  trecy = 0.0;
	  
	  Srem[i] = Spinning[0];
	  Xrem[i] = Spinning[1];
	  Mrem[i] = Spinning[2];
	  Krem[i] = Spinning[3];
	  Cosa[i] = Spinning[4];
	  Cosb[i] = Spinning[5];
	  Cosg[i] = Spinning[6];


	  double cj = pow(vthre,4.) / pow(vthre*vthre - Krem[i]*Krem[i],2.) - 1.;
	  if(cj < 0.0 && abs(cj) < 1.E-10)
	    cj = 0.0;
	  
	  rinfinite = pow(10.,rint) * sqrt(cj);
	  //if(CLevo == "yes")
	  rinfinite *= rclcorr;


	  
	  if(vthre < Krem[i] ||  (cj < 0.0 && abs(cj) > 1.E-10))
	    rinfinite = 1.E10;

	  out3<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<semi<<" "<<semi_ej<<" "<<semi_gw<<" "<<tfor[i]<<" "<<tSNe<<" "<<t12capt<<" "<<t3bb<<" "<<tdf<<" "<<t12<<" "<<tbbh<<" "<<tmer<<" "<<time<<" "<<nrecy<<" "<<pow(10., mint)*mclcorr<<" "<<rhalf*rclcorr<<" "<<pow(10.,mint)<<" "<<pow(10.,rint)<<" "<<tcc<<" "<<i<<" "<<label<<" "<<cluster<<" "<<" "<<Mrem[i]<<" "<<Srem[i]<<" "<<Xrem[i]<<" "<<Krem[i]<<" "<<vthre<<" "<<itot<<endl;	
	  
	  if(mpri > msmbhmax && tsmbh == 0.0){
	    tsmbh = time;
	    break;
	  }
	  if(nrecy > nbhs)
	    break;
	  
	  if(label == "ejected" || label == "GWrecoiled")
	    break;

	  
	  mpri = Mrem[i];
	  apri = Srem[i];

	  mass_ratio = -1;
	  nsafe = 0;       

	  double msec_prec = msec;
	  
	  if(stri_mrat != "nouniform"){
	    do{
	      mass_ratio = func.mratio(mpri, MRATIO_SLOPE, stri_mrat);
	      msec = mpri * mass_ratio;
	      for(int iii=0; iii<remn_sin.size()-1;iii++){
		if(msec > remn_sin[iii] && msec < remn_sin[iii+1]){
		  tsec = 0.5*(tdel_sin[iii] + tdel_sin[iii+1]) ;
		}
	      }
	      nsafe += 1;
	      if(nsafe > 1000){
		cout<<"mratio fails"<<endl;
		exit(0);
	      }
	    }while(msec < 1. || msec > 500.);
	    
	  }
	  else{
	    msec = -1;
	    ksec = 1.E30;
	    nsafe = 0;
	    if(mixer > mixing){
	      do{
		
		func.singBHt_new(zams_sin, remn_sin, tdel_sin, kick_sin, obslope, mslope, single_bh,saximus,sinimus,maximus,minimus,vthre);	  
		msec = single_bh[0];	 
		tsec = single_bh[1];
		ksec = single_bh[2];	    
		
		nsafe ++;
		if(single_bh[3] == 1)
		  break;
		else if(nsafe > 500){
		  cout<<npar_runtime<<" "<<i<<" "<<Z[i]<<" "<<msec<<" "<<tsec<<" "<<ksec<<endl;
		}
		if(nsafe > 1000){
		  cout<<"Something wrong " <<msec<<" "<<tsec<<" "<<ksec<<" "<<minimus<<" "<<vthre<<endl;
		  exit(0);
		}
		
	      }while(msec <= 0.0 || ksec > vthre);
	      
	      if(msec < minimus){
		cout<<"Second BH mass below mmin = "<<minimus<<" "<<msec<<endl;
	      }

	    }
	    else{
	      double MSLP = mslope;
	      msec = -1;
	      ksec = 1.E30;	      
	      do{
		singBHt_mix(mssx, msdx, mbsx, mbdx, tbsx, tbdx, vbsx, vbdx, mbhmix, tbhmix, vbhmix, MSLP, single_bh, saximus_mix, sinimus_mix, maximus_mix, minimus_mix, vthre);
		msec = single_bh[0];	  	  
		tsec = single_bh[1];
		ksec = single_bh[2];
		nsafe ++;
		if(nsafe > 1000)
		  break;
	      }while(msec <= 0.0 || ksec > vthre);
	      
	    }
	  }

	  if(msec == 0){
	    msec = msec_prec;
	    break;
	  }
	  
	  if(highgen == "yes" && msec < 65.){
	    double p_BPOP;
	    if(cluster == "globular")
	      p_BPOP = 295. / 1.E4; //this is valid for GCs with Z = 0.0002 (would it change with different metallicity?)
	    else if(cluster == "nuclear")
	      p_BPOP = 2740./ 1.E4;
	    else if(cluster == "young")
	      p_BPOP = 0.0;
	    
	    double p_hgen = func.rnd();
	    double mav = 30.;
	    if(p_hgen < p_BPOP){
	      double coef = 0.2 + 0.8*func.rnd();
	      msec = mav * (1. + msec/mav * (1.+mav/6.0));	       
	    }
	    nhigen += int(msec / mav);
	  }

	  
	  asec   = func.spin(msec,dynaS);

	  if(apri > 1. || asec > 1.){
	    cout<<"Critical error "<<mpri<<" "<<apri<<" "<<msec<<" "<<asec<<" "<<nrecy<<endl;
	    exit(0);
	  }


	  mu_red = mpri * msec / (mpri+msec);

	  	 
	  //tDF to sink, [(t3bb or t1-2)+t1-2bh] for BBH formation	
	  tdf = 0.42E9 * (10.*mstar/(mpri+msec)) * (trelax / 4.2E9) * pow(rinfinite/pow(10.,rint),1.74);
	  //if(CLevo=="yes"){
	  if(mclcorr*mhalf > 1.E1)
	    tdf *= (pow(rclcorr,3./2.-1.74)*sqrt(mclcorr));
	  else
	    tdf  = 1.E12;
	    
	  //}
	  
	  nan.str("");
	  nan<<tdf;
	  if(nan.str()=="nan"||nan.str()=="-nan"){
	    cout<<"DF ERROR: "<<tdf<<endl;
	    exit(0);
	  }
	  trecy += tdf;

	  //Time for BHs to pair --> assume it's already hard
	  //if(CLevo == "yes"){
	  
	  mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  
	  sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	  rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	  rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	  
	  
	  vthre = vthre_in * sqrt(mclcorr/rclcorr);
	  
	  semihard = 1./pow(sig_clu/30., 2.);
	  
	  
	  t12 = 3.E8 / zita * (0.01/fb) * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(mstar / (mpri+msec)) * (30./(mpri+msec+mper)) * (1./semihard);
	  t12 *= func.rndgen(1.0, 0.1);

	  
	  if(time + trecy < tfor[i] +tcc && time + trecy + t12 > tfor[i] +tcc){
	    time = tfor[i] + tcc;
	    
	    mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	    rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	    
	    sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	    rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	    rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	    
	    
	    vthre = vthre_in * sqrt(mclcorr/rclcorr);
	  }
	    
	    
	    //}
	  
	  
	  semihard = 1./pow(sig_clu/30., 2.); 

	  
	  t12 = 3.E8 / zita * (0.01/fb) * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(mstar / (mpri+msec)) * (30./(mpri+msec+mper)) * (1./semihard);
	  t12 *= func.rndgen(1.0, 0.1);

	  if(mclcorr*mhalf < 1.E1) // && CLevo == "yes")
	    t12 = 1.E12;
	    
	  
	  trecy += t12;

	  //if(CLevo == "yes"){
	    
	  mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  
	  sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	  rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	  rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	  
	  
	  vthre = vthre_in * sqrt(mclcorr/rclcorr);
	    
	  //}
	  

	  semihard = 1./pow(sig_clu/30., 2.); 
	  sma = semihard * func.Gaussian(1.,sigma_sma);
	  if(sigma_sma < 0)
	    sma = 0.1 + 0.1*func.rnd();
	  
	  //calculate aej vs agw to merge inside or outside --> inside (tgwpeters), outside (5/q3 * t2-1)
	  semi_ej = 0.07 * (mu_red * mper / (mpri + msec + mper) * mper/(mpri+msec)) * pow(vthre/50.,-2.);
	  semi_gw = 0.05 * pow((mpri+msec)/20,3./5.)*pow( (msec/mpri) / pow(1+(msec/mpri),2.) , 1./5.) * pow(sig_clu/30., 1./5.) * pow(1.E6 / rho_clu,1./5.);

	  //acrit = max(semi_ej,semi_gw); changed with the part below
	  double afin;
	  if(semi_ej > semi_gw){
	    afin = min(sma, semi_ej);	    
	  }
	  else{
	    afin = min(sma, semi_gw);
	  }
	  acrit = afin;
	  
	  
	  tbbh = 2.E7/zita * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(10. * mstar / (mpri+msec)) * (0.05 / acrit) * (20./(mpri+msec));
	  tbbh *= func.rndgen(1.0, 0.1);	  
	  
	  if(semi_ej > semi_gw)
	    acrit = 2. * afin / (2. + afin/semihard); 

	  
	  
	  tmer = 5. * (mpri + msec)/mper * tbbh;						 	  

	  if(mclcorr*mhalf < 1.E1){
	    semi_ej = -2.0;
	    semi_gw = -1.0;
	    acrit   = 1.E12;
	    tmer    = 1.E12;
	  }
	  
	  trecy+= tmer;
	  
	  //if(CLevo == "yes"){
	    

	  mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  
	  
	  sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	  rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	  rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	  
	  
	  vthre = vthre_in * sqrt(mclcorr/rclcorr);
	  
	  //}
	 	  
	  if(semi_ej < 0.0 && semi_gw < 0.0)
	    label="dissolved";
	  else if(semi_ej < semi_gw)	    
	    label="inside";
	  else	   
	    if(acrit < semi_gw)
	      label="inside";
	    else
	      label="ejected";

	  double teje = tmer;
       

	  ecc = sqrt(func.rnd());
	 	  
	  double clight = 3.E8;
	  double Grav   = 6.67E-11;
	  double Msun   = 1.99E30;
	  double Rau    = 1.5E11;	    
	  tmer = 5./256. * pow(clight,5.) / pow(Grav * Msun, 3.) * pow(Rau,4.)  * pow(acrit,4.)/(mpri*msec*(mpri+msec)) * pow(1.-pow(ecc,2.),3.5) / (1. + 73./24. * ecc*ecc + 37./96.*ecc*ecc*ecc*ecc) * 1./ (365.*24.*3600.*1.E0);
	  
	  if(tmer < teje && label == "ejected")
	    label = "inside";
	  
	  if(mclcorr*mhalf < 1.E1)
	    tmer = 1.E12;
	  
	  trecy += tmer;

	  nan.str("");
	  nan<<time;
	  if(nan.str()=="-nan" || nan.str()=="nan"){
	    cout<<"Critical error "<<tmer<<" "<<tbbh<<" "<<tdf<<" "<<t12<<" "<<time<<" "<<trecy<<endl;
	    exit(0);
	  }
	  
	  Spinning[0] = 0.0;
	  Spinning[1] = 0.0;
	  Spinning[2] = 0.0;
	  Spinning[3] = 0.0;
	  Spinning[4] = 0.0;
	  Spinning[5] = 0.0;
	  Spinning[6] = 0.0;

	  if(mpri>0.0 && msec>0.0)
	    func.SREM2(ndx, apri, asec, mpri, msec, align, Spinning);	  	 


	  if(mpri > 200.)
	    sig_clu = sqrt(sig_clu0*sig_clu0 + mpri / (0.1*pow(10.,rint)));
	  
	  nrecy += 1;
	  nrecy += nhigen;

	  //This will include all repeated mergers into the main catalogue ... 
	  if(time < 1.35e10){
	    zmer = func.inter(time / 1.E9, age, reds, redline);
	    zfor = func.inter(tfor[i] / 1.E9, age, reds, redline);
	    zsmbh= func.inter(tsmbh/1.E9, age, reds, redline);
	  }
	  else{
	    zmer = func.zred(time/1.E9);
	    zfor = func.zred(tfor[i]/1.E9);	  
	    zsmbh = func.zred(tsmbh/1.E9);
	  }
	 

	}while(1>0);

	
	//if(CLevo == "yes"){
	  
	mclcorr = func.mevol(Hubble-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	rclcorr = func.revol(Hubble-tfor[i], rhalf, mhalf, trelax0, CLfill, cluster);
	  
	//}
	

	
	

	if(time < 1.35e10){
	  zmer = func.inter(time / 1.E9, age, reds, redline);
	  zfor = func.inter(tfor[i] / 1.E9, age, reds, redline);
	  zsmbh= func.inter(tsmbh/1.E9, age, reds, redline);
	}
	else{
	  zmer = func.zred(time/1.E9);
	  zfor = func.zred(tfor[i]/1.E9);
	  zsmbh = func.zred(tsmbh/1.E9);
	}

	itot++;
	
	if(time < Hubble){
	  out<<itot<<" "<<Z[i]<<" "<<nrecy<<" "<<cluster<<" "<<REC<<" "<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<Mrem[i]<<" "<<Srem[i]<<" "<<Xrem[i]<<" "<<Krem[i]<<" "<<tfor[i]<<" "<<time<<" ";
	  out<<pow(10.,mint)<<" "<<pow(10.,rint)<<" "<<vthre<<" "<<label<<" "<<semi_ej<<" "<<semi_gw<<" "<<nbhs<<" "<<mhalf*mclcorr<<" "<<rhalf*rclcorr<<" "<<zmer<<" "<<zfor<<" "<<tsmbh<<" "<<zsmbh<<" "<<mzero<<" "<<ecc<<" "<<sma<<" "<<acrit<<" "<<tmer<<" "<<Cosa[i]<<" "<<Cosb[i]<<" "<<Cosg[i]<<endl; //Eccentricity added to output


	  
	  Ndyn_real++;
	  if(cluster == "young")
	    Nyou_real++;
	  else if(cluster == "globular")
	    Nglo_real++;
	  else if(cluster == "nuclear")
	    Nnuc_real++;
	  else{
	    cout<<"Error"<<endl;
	    exit(0);
	  }

	  
	  if(mpri!=0.0)      
	    out2<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<Mrem[i]<<" "<<Srem[i]<<" "<<Xrem[i]<<" "<<Krem[i]<<" "<<tfor[i]<<" "<<time<<endl;
	}



	mx.erase(mx.begin(),mx.end());
	my.erase(my.begin(),my.end());
	rx.erase(rx.begin(),rx.end());
	ry.erase(ry.begin(),ry.end());
      }
      
    }

    //    cout<<" ----> Actual BBH found "<<npar_runtime<<endl;


    
    zams_mix.erase(zams_mix.begin(),zams_mix.end());
    remn_mix.erase(remn_mix.begin(),remn_mix.end());
    zams_sin.erase(zams_sin.begin(),zams_sin.end());
    remn_sin.erase(remn_sin.begin(),remn_sin.end());
    tdel_mix.erase(tdel_mix.begin(),tdel_mix.end());
    tdel_sin.erase(tdel_sin.begin(),tdel_sin.end());
    kick_mix.erase(kick_mix.begin(),kick_mix.end());
    kick_sin.erase(kick_sin.begin(),kick_sin.end());

    
    // Cristiano 16/04/2025
	// Let's check how long it take per metallicity
	cout << "time for the metallicity function: " 
	<< chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count() 
	<< " seconds." << endl;
    
  }
  clout.close();

  out.close();
  out2.close();
  out3.close();
  hout.close();
  
  
  

  for(int i=0;i<100;i++){

    stringstream s;
    s<<i;

    string inname = "Catalogue"+s.str()+".txt";
    in.open(inname.c_str());
    if(!in.is_open())break;
    int Nsss = 0;
    do{
      string line;
      getline(in,line);
      Nsss += 1;
    }while(!in.eof());
    in.close();
    Nsss-=1;
    in.open(inname.c_str());
    
    string par2, par3;
    double apa[9];
    for(int k=0;k<Nsss;k++){
      in>>apa[0]>>apa[1]>>par2>>par3;
      for(int ll = 2;ll<9;ll++)in>>apa[ll];

      mpost.push_back(apa[6]);
      //cout<<mpost[k]<<endl;
    }
    in.close();
  }

  if(mpost.size() > 0){
    double *X;
    X = new double [mpost.size()];
    for(int i=0;i<mpost.size();i++)X[i]=mpost[i];
    func.histo(X,mpost.size(),30,"linear","Many_catalogues.txt");
    delete [] X;
  }

  stringstream Fcl;
  
  stringstream ggc;  stringstream gyc;  stringstream gnc;
  stringstream qgc;  stringstream qyc;  stringstream qnc;
  stringstream fgc;  stringstream fyc;  stringstream fnc;
  stringstream alg;
  

  Fcl<<DynOvTot;
  ggc<<probGC;  gyc<<probYC;  gnc<<probNC;
  alg<<INDEX_ALIGN;

  string ZDIS = zdist;
  string ZDYN = zdyn;
  
  stringstream metaldivi;
  metaldivi<<zdivi;

  string corr = correction;
  string type_mrat = MRATIO;

  //FINAL FILE MOVING
  string cmdstr_zero =  "./SIM_Fdyn"+Fcl.str()+"_Ngc"+ggc.str()+"_Nyc"+gyc.str()+"_Nnc"+gnc.str()+"isolS_"+isolS+"dynaS_"+dynaS+"_"+"MetalDivi_"+metaldivi.str()+"_"+alg.str()+"_"+ZDIS+"_"+ZDYN+"_Correction_"+corr;
  
  string SFR;
  if(sfr_iso == "katz13" || sfr_iso == "KR13")
    SFR = "KR13";
  else if(sfr_iso == "madau17" || sfr_iso == "MF17")
    SFR = "MF17";
  SFR += "_";
  if(sfr_clu == "katz13" || sfr_clu == "KR13")
    SFR += "KR13";
  else if(sfr_clu == "madau17" || sfr_clu == "MF17")
    SFR += "MF17";
  else if(sfr_clu == "burst")
    SFR += "burst";
  else if(sfr_clu == "elba18" || sfr_clu == "EB18")
    SFR += "EB18";
  else if(sfr_clu == "EB18_MF17")
    SFR += "EB18MF17";
  
  cmdstr_zero += "_SFR_"+SFR;
  cmdstr_zero += "_sfronly_";
  cmdstr_zero += sfr_only;

  if(bhseed == "yes"){
    stringstream strifseed;
    strifseed<<f_seed;
    cmdstr_zero += "HSBH"+strifseed.str();
  }
  if(bhseed == "vms"){
    stringstream strifseed;
    strifseed<<f_seed;
    cmdstr_zero += "VMSBH"+strifseed.str();
  }
  if(bhseed == "density"){
    stringstream strifseed;
    strifseed<<f_seed;
    cmdstr_zero += "VMSrho"+strifseed.str();
  }

  string TagRC = TagR;
  cmdstr_zero += "CluRh_";
  cmdstr_zero += TagRC;
  stringstream smastr;
  smastr << SSMA;
  cmdstr_zero += "_";
  cmdstr_zero += smastr.str();
  cmdstr_zero += "BHseed";
  cmdstr_zero += bhseed;
  cmdstr_zero += "_UMG";
  if(upgtp == "dicarlo")
    cmdstr_zero += upgtp;
  else
    cmdstr_zero += uppergap;
  
  cmdstr_zero += "_HighG";
  cmdstr_zero += highgen;

  
  stringstream smix;
  smix<<mixing;

  cmdstr_zero += "_mix_"+smix.str();

  if(highgen == "yes")
    cmdstr_zero += "_HGN";
  if(MonoZ == "yes")
    cmdstr_zero += "_MonoZ";
  if(cluster_test == "yes")
    if(cluster_test_env=="NC")
      cmdstr_zero += "_clteNC";
    else if(cluster_test_env=="YC")
      cmdstr_zero += "_clteYC";
    else      
      cmdstr_zero += "_cltest";
  
  if(CLfill == "under")
    cmdstr_zero += "U";
  else if(CLfill == "over")
    cmdstr_zero += "O";
  else if(CLfill == "critical")
    cmdstr_zero += "C";
  else if(CLfill == "noevo")
    cmdstr_zero += "N";
  else if(CLfill == "mix")
    cmdstr_zero += "M";
  else if(CLfill == "GG23")
    cmdstr_zero += "G";
  else
    cout<<"wrong selection of cluster evolution parameter"<<endl;

  
  
  cmdstr = "mkdir "+cmdstr_zero;  
  cmd = new char [cmdstr.length()+1];
  strcpy(cmd,cmdstr.c_str());
  ck = system(cmd);
  delete [] cmd;
  
  cmdstr = "mv *.txt *.dat "+cmdstr_zero;
  cmd = new char [cmdstr.length()+1];
  strcpy(cmd,cmdstr.c_str());
  ck = system(cmd);
  delete [] cmd;
  cmdstr = "mv *jpeg "+cmdstr_zero;
 
  cmd = new char [cmdstr.length()+1];
  strcpy(cmd,cmdstr.c_str());
  ck = system(cmd);
  cout<<"Copied everything into "<<cmd<<endl;
  delete [] cmd;
  

  
  
  YCrx.erase(YCrx.begin(),YCrx.end());
  YCry.erase(YCry.begin(),YCry.end());
  YCmx.erase(YCmx.begin(),YCmx.end());
  YCmy.erase(YCmy.begin(),YCmy.end());

  GCrx.erase(GCrx.begin(),GCrx.end());
  GCry.erase(GCry.begin(),GCry.end());
  GCmx.erase(GCmx.begin(),GCmx.end());
  GCmy.erase(GCmy.begin(),GCmy.end());

  NCrx.erase(NCrx.begin(),NCrx.end());
  NCry.erase(NCry.begin(),NCry.end());
  NCmx.erase(NCmx.begin(),NCmx.end());
  NCmy.erase(NCmy.begin(),NCmy.end());

  
  mpost.erase(mpost.begin(),mpost.end());
  
  delete [] Spinning;
  delete [] Npar;

  delete [] single_bh;
  
  delete [] Mrem;
  delete [] Srem;
  delete [] Xrem;
  delete [] Krem;
  delete [] Cosa;
  delete [] Cosb;
  delete [] Cosg;
  delete [] Z;
  delete [] Zi;

  cout<<"============="<<endl;
  cout<<"Final remarks"<<endl;
  cout<<"Number of sources requested = "<<N<<endl;
  cout<<"Number of simulated mergers = "<<Niso_real + Ndyn_real<<endl;
  cout<<"Actual number of sources "<<Niso_real<<" "<<Ndyn_real<<" "<<Nyou_real<<" "<<Nglo_real<<" "<<Nnuc_real<<endl;
  cout<<"f_dyn = "<<Ndyn_real * 1./(Ndyn_real + Niso_real)<<endl;
  if(Ndyn_real > 0)
    cout<<"F_YC, GC, NC / Dyn = "<<Nyou_real * 1./Ndyn_real <<" "<<Nglo_real * 1./Ndyn_real<<" "<<Nnuc_real * 1./Ndyn_real<<endl;
  cout<<"============="<<endl;
  
  const sec duration = clock::now() - before;

  
  cout<<"Elapsed time = "<<duration.count()<<"s"<<endl;

  
  

  return 0;
  
  
}

