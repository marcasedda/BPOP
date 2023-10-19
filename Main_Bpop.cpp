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

#define Hubble 1.35E10

// DATAFILES (Metal. distri, Single BHs, Binary BHs)
#define zPATH   "./gallazzi05ZDATA.ttt"
#define SINGPTH "./A5/" // "DATI_SingleBH/"
#define PATH    "./A5/" // "DATI_GiaMap18/"
#define PATHSIN "DATI_SingleBH/"

// GLOBAL
#define N        500
#define mmax     150.
#define mmin     18.5
#define mslope  -2.35
#define Zsun     0.02

// GENERAL PARAMS
#define delaytime "yes"
#define kick      "yes"

#define DynOvTot  1.0
#define mixing    0.0
#define fbin      1.0

//CLUSTER EVOLUTION
#define CLevo        "no"
#define CLfill       "under"
#define MonoZ        "yes"
#define cluster_test "no"
#define cluster_test_env "NC"
#define msmbhmax 1.E10

//BH SEED
#define bhseed   "yes"
#define vms      "75"
#define f_seed    0.1
#define maxseed   1500.0
#define minseed   100.0
#define seedslope 2.0

//SPINS
#define spinlb  "gaussian02"
#define obslope  0.0

//STAR FORMATION
#define SFRTYPE_ISO "MF17"
#define SFRTYPE_CLU "KR13"
#define sfr_only "no"

//DYNAMICAL SCENARIO
#define zdyn       "logflat"
#define zslope     -1.5

//DEPRECATED
#define MRATIO "nouniform" //"uniform","pwl","gss","nouniform"
#define MRATIO_SLOPE -1.0

//DYNAMICAL FRACTIONS
#define probGC 0.0
#define probYC 0.0
#define probNC 1.0

//DEPRECATED
#define MinQGC 0.0
#define recyGC 1.0
#define MinQYC 0.0
#define recyYC 1.0
#define MinQNC 0.0
#define recyNC 1.0

//HIGHLY DEPRECATED
#define vmaxNC 100.E0
#define crazy  "no"

//ISOLATED SCENARIO
#define zdist       "gallazzi"
#define zdivi       0.0
#define correction  "yes" 
#define slope       -1.5

#define INDEX_ALIGN 8.0

//DEPRECATED
#define PATHse "" //"CC15alpha1"
#define PATHsp "" //"_spin01"

//FILE STRUCTURE
#define kpar 20
#define spar 5
#define numZ 12

//DEPRECATED
// REDUCING CALCULATIONS
#define onlyrepetita "no"
#define bin_st 50
#define nsize 200
#define tsize 200
#define vsize 200

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

  if(CLevo == "yes" && (CLfill!= "filling" && CLfill != "under" && CLfill != "postcol" && CLfill != "nuclear" && CLfill != "nuclearF")){
    cout<<"Select CLfill = filling/under/postcol/nuclear "<<endl;
    cout<<"You've chosen: "<<CLfill<<endl;
    exit(0);
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
  in.open("redshift_time.ttt");
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
  //if(yes=="Y") ck = system(cmd);
  delete [] cmd;
  //delete [] yes;

  string sfr_iso = SFRTYPE_ISO;
  string sfr_clu = SFRTYPE_CLU;
  if(cluster_test == "yes"){
    if(cluster_test_env == "NC")
      sfr_clu = "bigbang";
    else
      sfr_clu = "burst";
  }
  int Nsrc = N;
    
  double *Spinning;
  Spinning = new double [4];

  /*
  double *sss;
  sss = new double [1000];
  for(int i=0;i<1000;i++)
    sss[i] = 700.*func.rndgen(1.0, 0.1);
  func.histo(sss,1000,30,"linear","test_general_GSS.txt");
  delete [] sss;
  exit(0);
  */
  
  Spinning[3] = -1.E30;
  
  string path   = PATH;
  string pathse = PATHse;
  string pathsp = PATHsp;

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

  string singpthA = SINGPTH;  
  double ndx;
  double apri;
  double asec;
  double mpri;
  double msec;
  double tpri;
  double tsec;
  double kpri;
  double ksec;
  
  double vpri;
  double vsec;
  
  string align;

  double rnd;
  string label = "0";

  double qmin, recy;
  string cluster = "none";

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

      
      string clufile = "ClustersRhM/scottrule_"+clvari+cltype+".ttt";
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
  
  double *Z;
  Z = new double [Nsrc];
  double *Zi;
  Zi = new double [Nsrc];

  int *Npar;
  Npar = new int [numZ];
  for(int i=0;i<numZ;i++) Npar[i] = 0.0;



  
  vector<double> miso1,miso2;
  vector<double> aiso1,aiso2;
  double P;

  //CREATE THE CUMULATIVE AND SELECT FROM THAT... (Gallazzi+05)
  string zpth = zPATH;
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
  
  do{
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
  
  Nsrc = Niso + Nyc + Ngc + Nnc;


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
  for(int i=0;i<Ngc+Nnc;i++){
    double zred0 = func.sfr_red(sfr_clu);      
    double toff  = 1.E9*func.inter(zred0, reds, age, redline);
    tfor.push_back(toff);
  }

  //ASSIGNING METALLICITIES to ISO and CLU BBHs
  double z;

  for(int i=0;i<Niso+Nyc;i++){

    if(zdist == "logflat"){
      z = log10(met[0])+(log10(met[12])-log10(met[0]))*func.rnd();
      z = pow(10.,z);
    }
    else if(zdist == "powerlaw"){
      P = func.rnd();
      z = pow( (1.-P)*pow(met[12],1.+zslope) + P*pow(met[0],1.+zslope) , 1./(1.+zslope));
      if(zslope == -1.0){
	z = (1.-P)*log(met[12]) + P*log(met[0]);
	z = exp(z);
      }
      if(zslope == 0.0){
	z = met[12]*func.rnd();
	do{z = met[12]*func.rnd();}while(z>met[11]);
      }
    }
    else if(zdist == "gallazzi" || zdist == "Gallazzi" || zdist == "Gallazzi05" || zdist == "gallazzi05"){
	
      double Pzrnd = func.rnd();	
      double zfrs, zsec;
      zfrs = -3;
      zsec = -3;
      for(int ll=0;ll<zgal.size();ll++){
	if(Pzrnd>Nz[ll])
	  zfrs = zgal[ll];	
	else{
	  zsec = zgal[ll];
	  break;
	}	 	 
      }
      if(zfrs == -3){
	//cout<<"out of range (below) "<<Pzrnd<<" "<<Nz[0]<<endl;
	zfrs = zgal[0] + (Pzrnd-Nz[0])/(Nz[1]-Nz[0])*(zgal[1]-zgal[0]);
	if(zfrs<0.0) zfrs = 0.000001;
      }	
      if(zsec == -3){
	//cout<<"out of range (above) "<<Pzrnd<<" "<<Nz[zgal.size()-1]<<endl;
	zsec = zgal[zgal.size()-1];
      }
      
      double lnzt = (zsec-zfrs)*func.rnd();      
      double lnz = zfrs+lnzt;
      
      z = Zsun*pow(10.,lnz);
      //cout<<zfrs<<" "<<zsec<<" "<<lnzt<<" "<<z<<endl;
    }   
    else{
      cout<<"Z distribution not recognized "<<endl;
      exit(0);
    }


    if(delaytime=="yes" && (zdist == "gallazzi" || zdist == "Gallazzi" || zdist == "Gallazzi05" || zdist == "gallazzi05") ){

      ///INTRODUCING A CORRECTION FOR THE DELAY TIME...///

      double tof = tfor[i] / 1.E9;      
      double zinit = z;
      double red_del = func.inter(tof, age, reds, redline);
      double logz = log10(zinit) - 0.074*pow(red_del,1.34);
      z = pow(10.,logz);

      
      /*PUT THIS BACK*/
      if(z<0.9*met[0]){
	do{
	  double zred0 = func.sfr_red(sfr_iso);      
	  double toff = 1.E9*func.inter(zred0, reds, age, redline);
	  tfor[i] = toff;
	  tof = tfor[i] / 1.E9;      
	  red_del = zred0;

	  logz = log10(zinit) - 0.074*pow(red_del,1.34);

	  z = pow(10.,logz);
	  
	  if(z > 0.9*met[0])
	    break;
	  
	  
	}while(1>0);
	}
      
      
      if(z<0.9*met[0]){
	cout<<"--> WARNING"<<endl;
	exit(0);
	}

    }

      
    Z[i] = z;
    if(Z[i]>Z1) Z1 = Z[i];
    if(Z[i]<Z0) Z0 = Z[i];

  }

  cout<<"isolated binaries Z assigned"<<endl;


  int nout=0;
  for(int i=Niso+Nyc;i<Niso+Nyc+Ngc+Nnc;i++){

    double Zgcmin = 0.074*pow(2.,1.34) + log10(0.0005);
    double Zgcmax = 0.074*pow(2.,1.34) + log10(0.001);

    if(zdyn == "logflat"){
      z = Zgcmin+(Zgcmax-Zgcmin)*func.rnd();
      z = pow(10.,z);
      if(z > met[12]){
	cout<<"argh: "<<z<<" "<<Zgcmin<<" "<<Zgcmax<<endl;
      }
	
    }
    else if(zdyn == "powerlaw"){
      P = func.rnd();
      z = pow( (1.-P)*pow(met[12],1.+zslope) + P*pow(met[0],1.+zslope) , 1./(1.+zslope));
      if(zslope == -1.0){
	z = (1.-P)*log(met[12]) + P*log(met[0]);
	z = exp(z);
      }
      if(zslope == 0.0){
	z = met[12]*func.rnd();
	do{z = met[12]*func.rnd();}while(z>met[11]);
      }
    }
    else if(zdyn== "gallazzi" || zdyn == "Gallazzi" || zdyn == "Gallazzi05" || zdyn == "gallazzi05"){
	
      double Pzrnd = func.rnd();	
      double zfrs, zsec;
      zfrs = -3;
      zsec = -3;
      for(int ll=0;ll<zgal.size();ll++){
	if(Pzrnd>Nz[ll])
	  zfrs = zgal[ll];	
	else{
	  zsec = zgal[ll];
	  break;
	}	 	 
      }
      if(zfrs == -3){
	//cout<<"out of range (below) "<<Pzrnd<<" "<<Nz[0]<<endl;
	zfrs = zgal[0] + (Pzrnd-Nz[0])/(Nz[1]-Nz[0])*(zgal[1]-zgal[0]);
	if(zfrs<0.0) zfrs = 0.000001;
      }	
      if(zsec == -3){
	//cout<<"out of range (above) "<<Pzrnd<<" "<<Nz[zgal.size()-1]<<endl;
	zsec = zgal[zgal.size()-1];
      }
      
      double lnzt = (zsec-zfrs)*func.rnd();      
      double lnz = zfrs+lnzt;
      
      z = Zsun*pow(10.,lnz);

      /*double Pzrnd = func.rnd();	
      double zfrs, zsec;
      zfrs = -3;
      zsec = -3;
      for(int ll=0;ll<zgal.size();ll++){
	if(Pzrnd>Nz[ll])
	  zfrs = zgal[ll];	
	else{
	  zsec = zgal[ll];
	  break;
	}	 	 
	}*/
    }
    else{
      cout<<"Metallicity distribution for dynamical unknown"<<endl;
      exit(0);
    }

    if((delaytime=="yes" && (zdist == "gallazzi" || zdist == "Gallazzi" || zdist == "Gallazzi05" || zdist == "gallazzi05" || zdist=="logflat")) && MonoZ!="yes"){
      ///INTRODUCING A CORRECTION FOR THE DELAY TIME...///

      double tof = tfor[i] / 1.E9;      
      double zinit = z;

      double red_del = func.inter(tof, age, reds, redline);

      double logz = log10(zinit) - 0.074*pow(red_del,1.2);
      z = pow(10.,logz);
      
      if(z<0.9*metdyn[0] || z>metdyn[1]){
	do{
	  double zred0 = func.sfr_red(sfr_clu);      
	  double toff = 1.E9*func.inter(zred0, reds, age, redline);
	  tfor[i] = toff;
	  tof = tfor[i] / 1.E9;      
	  red_del = zred0;

	  logz = log10(zinit) - 0.074*pow(red_del,1.2);

	  z = pow(10.,logz);

	  /*cout<<" -- > WARNING: "<<zred0<<" "<<toff<<" "<<tof<<" "<<logz<<" "<<z<<" "<<zinit<<endl;
	    exit(0);*/
	  if(z > 0.9*metdyn[0] && z < metdyn[12])
	    break;
	  
	}while(1>0);
     
      }
      
      //cout<<"Checking metallicity "<<tof<<" "<<red_del<<" "<<logz<<" "<<log10(zinit)<<" "<<zinit<<endl;
      
      if(z<0.9*metdyn[0]){
	cout<<"--> WARNING"<<endl;
	exit(0);
      }
    }

    Z[i] = z;
    if(Z[i]>Z1) Z1 = Z[i];
    if(Z[i]<Z0) Z0 = Z[i];
      
  }
  //cout<<"Source outside metallicity range ... "<<nout<<endl;
  cout<<"dynamical binaries Z assigned"<<endl;

  if(MonoZ == "yes")
    for(int i=0;i<Nsrc;i++){
      Z[i] = 0.0002;
      Z0 = Z[i];
      Z1 = Z[i];
    }
  
  
  //   CLOSEST DISCRETE Z FOR FILES READING?
  
  int NzIS = 0;
  int NzGC = 0;
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
	else       NzGC += 1;
	
	Npar[k] += 1;
      }
      /*else if(Z[i]==met[numZ-1]){
	if(i<Niso+Nyc) NzIS += 1;
	else           NzGC += 1;	
	Npar[numZ-1] += 1;	
	} */     
    }
    out<<0.5*(pow(10.,Z0)+pow(10.,Z1))<<" "<<Npar[k]<<" "<<NzIS<<" "<<NzGC<<endl;
    NzIS = 0;
    NzGC = 0;
    
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
  
  int Ngap = 0;
  string REC = "RecNO";
  cout<<"Header: src No., metallicity, scenario, cluster type"<<endl;
  ofstream out2;

  string outname = "Catalogue.txt";
  out.open(outname.c_str());
  out2.open("Catalogue_clean.txt");

  cout<<"    ISOLATED BINARIES   "<<endl;

  for(int k=0;k<numZ;k++)Npar[k] = 0;
  for(int i=0;i<Niso;i++)Zi[i]=Z[i];
    
  for(int i=0;i<Niso;i++){  

    for(int k=0;k<numZ;k++){
      if(Z[i]<met[k+1] && Z[i]>=met[k]){
	Npar[k] += 1;
	double delta1 = abs(Z[i] - met[k+1]);
	double delta2 = abs(Z[i] - met[k]);
	if(delta1<delta2)       Z[i] = met[k+1];
	else if(delta1>delta2)  Z[i] = met[k];
	else if(delta1==delta2){
	  P=func.rnd();
	  if(P<0.5) Z[i] = met[k+1];
	  else      Z[i] = met[k];
	}
	break;
      }
      else if(Z[i]==met[numZ-1]){
	Npar[numZ-1] += 1;	
      } 
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
	aiso1.push_back(par[7]);
	aiso2.push_back(par[13]);
      }
      else{
	miso1.push_back(par[2]);
	miso2.push_back(par[1]);
	aiso1.push_back(par[13]);
	aiso2.push_back(par[7]);
      }
      
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

    for(int i=0; i<Npar[k];i++){
    
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
      }while(miso1[ext] > mobs*1.2 || miso1[ext] < mobs*0.8);

      mpri = miso1[ext];
      msec = miso2[ext];
      tdel = tdel_iso[ext];
      
      //BH SPINS SELECTION (TO BE DECIDED)    
      if(isospin != "giacobbo" || isospin != "Giacobbo"){
	apri   = func.spin(mpri,isospin);
	asec   = func.spin(msec,isospin);	
      }
      else{
	apri = aiso1[ext];
	asec = aiso2[ext];
      }

      Spinning[0] = 0.0;
      Spinning[1] = 0.0;
      Spinning[2] = 0.0;
      if(mpri > 0.0 && msec > 0.0)
	func.SREM2(ndx, apri, asec, mpri, msec, align, Spinning);
	
      Srem[i] = Spinning[0];
      Xrem[i] = Spinning[1];
      Mrem[i] = Spinning[2];
      Krem[i] = Spinning[3];
      

      if(mpri < msec){
	double mpri_sec = mpri;
	mpri = msec;
	msec = mpri;
      }


      tdel += tfor[i];
      
      double zmer = func.inter(tdel / 1.E9, age, reds, redline);
      double zfor = func.inter(tfor[i] / 1.E9, age, reds, redline);

      out<<itot<<" "<<met[k]<<" "<<label<<" "<<cluster<<" "<<REC<<" "<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<Mrem[i]<<" "<<Srem[i]<<" "<<Xrem[i]<<" "<<Krem[i]<<" "<<tfor[i]<<" "<<tdel;
      out<<" 0.0 0.0 0.0 none 0.0 0.0 0.0 0.0 0.0"<<zmer<<" "<<zfor<<" -1 -1" <<endl;
      if(mpri!=0.0)	
	out2<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<Mrem[i]<<" "<<Srem[i]<<" "<<Xrem[i]<<" "<<Krem[i]<<" "<<tfor[i]<<" "<<tdel<<"-1 -1"<<endl;
	  
      mpri = -10000.0;


    }
    
    miso1.erase(miso1.begin(),miso1.end());
    aiso1.erase(aiso1.begin(),aiso1.end());
    miso2.erase(miso2.begin(),miso2.end());
    aiso2.erase(aiso2.begin(),aiso2.end());
    tdel_iso.erase(tdel_iso.begin(),tdel_iso.end());
    
    
  }
  out.close();
  out2.close();


  outname = "Catalogue.txt";
  out.open(outname.c_str(),ios::app);
  out2.open("Catalogue_clean.txt");      
  ofstream out3;
  out3.open("Catalogue_multiple_dyn.txt");

  ofstream hout;
  hout.open("Larger_than_tH.txt");
  
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
  single_bh = new double [3];
  single_bh[0] = 0.0;
  single_bh[1] = -1.0;
  single_bh[2] = 0.0;
  
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
  
  total_sources = 0.0;
  int numme = numZ;
  if(Nyc + Ngc + Nnc == 0)
    numme = 1;    


  int nreal = 0;
  
  if(MonoZ == "yes")
      Zmin = log10(Z[Niso+1]);
    
  for(int k=0;k<numme;k++){
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

      Zi[i] = Z[i];
      
      Z0 = Zmin + k*dZ;
      Z1 = Zmin + (k+1)*dZ;

      //cout<<pow(10.,Z0)<<" "<<pow(10.,Z1)<<" "<<Z[i]<<" "<<Z[i] - pow(10.,Z0) << endl;

      if( (Z[i] < pow(10.,Z1) && Z[i] >= pow(10.,Z0)) || (MonoZ == "yes" && abs(Z[i] - pow(10.,Z0))<1.E-10)) {

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
	double tSNe, mper, trelax, tcc, tdf, semihard, kappa, semi, ecc, t12, acrit, tbbh, tmer, mu_red;
	double t3bb, t12capt, tbbhform;
	double time,nsafe_glob,nsafe,mass_ratio, mixer;
	string stri_mrat;

	double mclcorr = 1.0;
	double rclcorr = 1.0;

	zita = 1.0;
	logL = 15.;
	fb   = 0.01;

	mper = 15.0;
	
	mixer = func.rnd();

	double vlimiting = 2.;

	if(mixer > mixing){
	  if(met[k] == 0.02)
	    vlimiting = 10.;
	}

	string manuel = "daje";

	int nsafe_cal = 0;
	do{
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
	      
	      
	      if(cluster=="nuclear" || cluster_test == "yes"){
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
	    }while( pow(10.,rint) > 10.*0.8*pow(pow(10.,mint)/1000.,0.2));
	    
	    
	    
	    vthre = func.vescape(pow(10.,rint),pow(10.,mint),pcluster);
	  }while(vthre < vlimiting || mint < 3.0 || pow(10.,mint)/pow(pow(10.,rint),3.) < 100.);
	  

	  
	  
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

	  
	  trelax = 4.2E9 * (15./logL) * pow(rhalf/4.0,1.5) * sqrt(mhalf/1.E7) ;

	  tcc = 0.138*mhalf/(150. * log(0.11*mhalf/150.))*sqrt(pow(rhalf*3.08E16,3.)/(6.67E-11*1.99E30*mhalf))/(365.*24.*3600.);
	  
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
	  
	  
	  if(mixer > mixing){		
	    mpri = -1;
	    kpri = 1.E30;
	    do{
	      func.singBHt_new(zams_sin, remn_sin, tdel_sin, kick_sin, obslope, mslope, single_bh, saximus,sinimus,maximus,minimus, vthre);
	      mpri = single_bh[0];	 
	      tpri = single_bh[1];
	      kpri = single_bh[2];
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
	      nsafe ++;
	      if(nsafe > 1000)
		break;
	    }while(mpri <= 0.0 || kpri > vthre);
	  }
	  apri   = func.spin(mpri,dynaS);	
	  
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
		if(nsafe > 500){
		  cout<<npar_runtime<<" "<<i<<" "<<Z[i]<<" "<<msec<<" "<<tsec<<" "<<ksec<<endl;
		}
		if(nsafe > 1000){
		  cout<<"Something wrong " <<endl;
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
	      int nsafe = 0;
	      do{
		singBHt_mix(mssx, msdx, mbsx, mbdx, tbsx, tbdx, vbsx, vbdx, mbhmix, tbhmix, vbhmix, MSLP, single_bh, saximus_mix, sinimus_mix, maximus_mix, minimus_mix, vthre);
		msec = single_bh[0];	  	  
		tsec = single_bh[1];
		ksec = single_bh[2];	    
		if(nsafe > 1000)
		  break;
	      }while(msec <= 0.0 || ksec > vthre);
	      
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

	  if(bhseed == "vms"){
	    
	    double A_vms, B_vms;
	    if(vms == "max"){
	      A_vms = 2.1041788;
	      B_vms = 0.08170471;
	    }
	    else if(vms == "75"){
	      A_vms = 2.22932958;
	      B_vms = 0.04713929;
	    }
	    else if(vms == "low"){
	      A_vms = 2.31804413;
	      B_vms = 0.01607658;
	    }
	    else{
	      cout<<"Please choose vms == low, 75, max"<<endl;
	      cout<<"you have chosen "<<vms<<endl;
	      exit(0);
	    }
	    double rho_vms = mint - 3.*rint;
	    mpri = pow(10., A_vms + B_vms * rho_vms );
	  }
	  if(bhseed == "yes"){
	    double dmy = func.rnd();
	    if(dmy < f_seed){
	      double mimmi = func.rnd();
	      mpri = pow(mimmi * pow(maxseed, 1.-seedslope) + (1.-mimmi) * pow(minseed,1.-seedslope) , 1./(1.-seedslope));	      
	    }
	    
	  }
	  
	  asec   = func.spin(msec,dynaS);
	  
	  
	  //FIRST MERGER
	  Spinning[0] = 0;
	  Spinning[1] = 0;
	  Spinning[2] = 0;
	  Spinning[3] = 0;
	  if(mpri>0.0 && msec>0.0)
	    func.SREM2(ndx, apri, asec, mpri, msec, align, Spinning);	  
	  Srem[i] = Spinning[0];
	  Xrem[i] = Spinning[1];
	  Mrem[i] = Spinning[2];
	  Krem[i] = Spinning[3];
	  
	  
	  
	  time = tfor[i];
	  tSNe = max(tpri,tsec);

	  if(bhseed == "vms")
	    tSNe = max(5.0e6, tcc);
	    
	  
	  rhalf = pow(10.,rint);
	  mhalf = pow(10.,mint);
	   
	  trelax = 4.2E9 * (15./logL) * pow(rhalf/4.0,1.5) * sqrt(mhalf/1.E7) ;		  
	  //tDF to sink	
	  tdf = 0.42E9 * (10.*mstar/(mpri+msec)) * (trelax / 4.2E9);	  	  



	  
	  if(tdf > tSNe)
	    time += tdf;	  
	  else
	    time += tSNe;

	  if(CLevo == "yes"){
	    
	    mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	    rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax, CLfill);

	    //cout<<"M(t), R(t) = "<<time-tfor[i]<<" "<<mclcorr*mhalf<<" "<<rclcorr*rhalf<<endl;
	    
	    sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	    rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	    rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);

	    vthre = vthre_in * sqrt(mclcorr/rclcorr);
	    
	    
	    semihard = 1./pow(sig_clu/30., 2.); //2.*6.67E-11*1.99E30/(1.E6*1.5E11) * (mpri+msec) / (sig_clu*sig_clu);
	  
	    t12capt = 3.E9 * (0.01 / fb) * (1.E6 / rho_cubicpc) * (sig_clu / 30.) * 10./(mpri + msec + mper) * 1./semihard;
	    t3bb = 4.E9 * pow(1.E6 / rho_clu,2.) * pow(sig_clu/(zita * 30.),9.) * pow(mstar / 30. * 10., 9./2.) * pow(10./30.,-5.);
	    
	    double fbin_rnd = func.rnd();
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

	    
	    if(time + tbbhform > tcc && time < tcc){
	      time = tcc;
	    
	      
	      mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	      rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	      

	      
	      sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	      rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	      rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	      
	      vthre = vthre_in * sqrt(mclcorr/rclcorr);
	    }	    
	    
	    
	  }

	 	  
	  //Binary formation and BH captures

	  
	  semihard = 1./pow(sig_clu/30., 2.); //2.*6.67E-11*1.99E30/(1.E6*1.5E11) * (mpri+msec) / (sig_clu*sig_clu);
	  
	  t12capt = 3.E9 * (0.01 / fb) * (1.E6 / rho_cubicpc) * (sig_clu / 30.) * 10./(mpri + msec + mper) * 1./semihard;
	  t3bb = 4.E9 * pow(1.E6 / rho_clu,2.) * pow(sig_clu/(zita * 30.),9.) * pow(mstar / 30. * 10., 9./2.) * pow(10./30.,-5.);

	  double fbin_rnd = func.rnd();
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

	  

	  
	  if(CLevo == "yes"){
	    
	    mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	    rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax, CLfill);

	    sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	    rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	    rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);

	    vthre = vthre_in * sqrt(mclcorr/rclcorr);
	    
	    t12 = 3.E8 / zita * (0.01/fb) * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(mstar / (mpri+msec)) * (30./(mpri+msec+mper)) * (1./semihard);
	    t12 *= func.rndgen(1.0, 0.1);
	    
	    if(mclcorr * mhalf < 1.E1)
	      t12 = 1.E12;
	    
	    
	    if(time < tcc && time + t12 > tcc){
	      time = tcc;

	      mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	      rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	      
	      sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	      rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	      rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	      
	      vthre = vthre_in * sqrt(mclcorr/rclcorr);
	    }
	  }
	    

	  //Time for BHs to pair --> assume it's already hard //see Antonini&Rasio2016
	  
	  t12 = 3.E8 / zita * (0.01/fb) * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(mstar / (mpri+msec)) * (30./(mpri+msec+mper)) * (1./semihard);
	  t12 *= func.rndgen(1.0, 0.1);

	  if(mclcorr * mhalf < 1.E1)
	    t12 = 1.E12;

	  time += t12;

	 
	  
	  
	  if(CLevo == "yes"){

	    mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	    rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax, CLfill);

	    sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	    rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	    rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	  

	    vthre = vthre_in * sqrt(mclcorr/rclcorr);
	  }

	  
	  //calculate aej vs agw to merge inside or outside --> inside (tgwpeters), outside (5/q3 * t2-1)

	  semi_ej = 0.07 * (mu_red * mper / (mpri + msec + mper) * mper/(mpri+msec)) * pow(vthre/50.,-2.);
	  semi_gw = 0.05 * pow((mpri+msec)/20,3./5.)*pow( (msec/mpri) / pow(1+(msec/mpri),2.) , 1./5.) * pow(sig_clu/30., 1./5.) * pow(1.E6 / rho_clu,1./5.);
	  
	  acrit = max(semi_ej,semi_gw);
	  acrit = min(acrit,0.01);
	  
	  tbbh = 2.E7/zita * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(10. * mstar / (mpri+msec)) * (0.05 / acrit) * (20./(mpri+msec));
	  tbbh *= func.rndgen(1.0, 0.1);	  
	  
	  if(mclcorr * mhalf < 1.E1){
	    tbbh = 1.E12;
	    semi_ej = -2.0;
	    semi_gw = -1.0;
	  }
	  
	  tmer = 5. * (mpri + msec)/mper * tbbh;						 	  
	    
	  time += tmer;
	  
	  if(nsafe_cal > 10000)
	    break;
	  
	}while(time > Hubble);

	
	if(time > Hubble){  // || nsafe_cal > 1){
	  hout<<"Hubble time -- "<<pow(10.,mint)<<" "<<pow(10.,rint)<<" "<<mclcorr<<" "<<rclcorr<<" "<<trelax<<" "<<vthre<<" "<<sig_clu<<" "<<tdf<<" "<<tbbhform<<" "<<t12<<" "<<tmer<<" "<<nsafe_cal<<endl;
          //exit(0);
	}
	  
	
	if(semi_ej < semi_gw)	    
	  label="inside";
	else
	  label="ejected";
	
	
	ecc = sqrt(func.rnd());
	
	double clight = 3.E8;
	double Grav   = 6.67E-11;
	double Msun   = 1.99E30;
	double Rau    = 1.5E11;	    
	tmer = 5./256. * pow(clight,5.) / pow(Grav * Msun, 3.) * pow(Rau,4.)  * pow(acrit,4.)/(mpri*msec*(mpri+msec)) * pow(1.-pow(ecc,2.),3.5) / (1. + 73./24. * ecc*ecc + 37./96.*ecc*ecc*ecc*ecc) * 1./ (365.*24.*3600.*1.E0);
	
	if(time + tmer > Hubble && time < Hubble && mclcorr*mhalf > 1.E1){
	  do{
	    ecc = sqrt(func.rnd());
	    tmer = 5./256. * pow(clight,5.) / pow(Grav * Msun, 3.) * pow(Rau,4.)  * pow(acrit,4.)/(mpri*msec*(mpri+msec)) * pow(1.-pow(ecc,2.),3.5) / (1. + 73./24. * ecc*ecc + 37./96.*ecc*ecc*ecc*ecc) * 1./ (365.*24.*3600.*1.E0);	 
	  }while(time + tmer < Hubble);
	}
	else if(mclcorr*mhalf < 1.E1)
	  tmer = 1.E12;
	
	time += tmer;
	
	/*if(time > Hubble){
	  cout<<"Hubble time -- tGW "<<tmer<<" "<<time<<" "<<cluster<<" "<<label<<endl;
	  									  //exit(0);
										  }*/


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
	Spinning[0] = 0;
	Spinning[1] = 0;
	Spinning[2] = 0;
	Spinning[3] = 0;
	if(mpri>0.0 && msec>0.0)
	  func.SREM2(ndx, apri, asec, mpri, msec, align, Spinning);	  

	Srem[i] = Spinning[0];
	Xrem[i] = Spinning[1];
	Mrem[i] = Spinning[2];
	Krem[i] = Spinning[3];
	
	vpri = Spinning[3];


	
	if(CLevo == "yes"){
	  

	  mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	  rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax, CLfill);

	  sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	  rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	  rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	  
	  
	  vthre = vthre_in * sqrt(mclcorr/rclcorr);

	  if(mclcorr * mhalf < 1.E1)
	    manuel = "none";
	}

	double rinfinite = pow(10.,rint) * sqrt(pow(vthre,4.) / pow(vthre*vthre - vpri*vpri,2.) - 1.);       
	if(CLevo == "yes")
	  rinfinite *= rclcorr;

	
	int nrecy = 0;
	double trecy = 0.0;
	double tsmbh = 0.0;

	double sig_clu0 = sig_clu;

	do{
	  
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
	    break;
	  }


	  if(Spinning[3] > vthre && label != "ejected"){
	    label = "recoiled";
	    Srem[i] = Spinning[0];
	    Xrem[i] = Spinning[1];
	    Mrem[i] = Spinning[2];
	    Krem[i] = Spinning[3];
	    break;
	  }
	  
	  if(manuel == "none")
	    break;

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

	  rinfinite = pow(10.,rint) * sqrt(pow(vthre,4.) / pow(vthre*vthre - Krem[i]*Krem[i],2.) - 1.);
	  if(CLevo == "yes")
	    rinfinite *= rclcorr;
	  
	  if(vthre < Krem[i] || pow(vthre,4.) / pow(vthre*vthre - Krem[i]*Krem[i],2.) - 1. < 0.0)
	    rinfinite = 1.E10;
	  
	  out3<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<semi<<" "<<semi_ej<<" "<<semi_gw<<" "<<tfor[i]<<" "<<tSNe<<" "<<t12capt<<" "<<t3bb<<" "<<tdf<<" "<<t12<<" "<<tbbh<<" "<<tmer<<" "<<time<<" "<<nrecy<<" "<<pow(10., mint)*mclcorr<<" "<<rhalf*rclcorr<<" "<<pow(10.,mint)<<" "<<pow(10.,rint)<<" "<<tcc<<" "<<i<<" "<<label<<" "<<cluster<<endl;	

	  if(mpri > msmbhmax && tsmbh == 0.0){
	    tsmbh = time;
	    break;
	  }
	  if(nrecy > nbhs)
	    break;
	  
	  if(label == "ejected" || label == "recoiled")
	    break;
	  
	  mpri = Mrem[i];
	  apri = Srem[i];

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
		if(nsafe > 500){
		  cout<<npar_runtime<<" "<<i<<" "<<Z[i]<<" "<<msec<<" "<<tsec<<" "<<ksec<<endl;
		}
		if(nsafe > 1000){
		  cout<<"Something wrong " <<endl;
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

	  asec   = func.spin(msec,dynaS);

	  if(apri > 1. || asec > 1.){
	    cout<<"Critical error "<<mpri<<" "<<apri<<" "<<msec<<" "<<asec<<" "<<nrecy<<endl;
	    exit(0);
	  }

	  
	  mu_red = mpri * msec / (mpri+msec);

	  	 
	  //tDF to sink, [(t3bb or t1-2)+t1-2bh] for BBH formation	
	  tdf = 0.42E9 * (10.*mstar/(mpri+msec)) * (trelax / 4.2E9) * pow(rinfinite/pow(10.,rint),1.74);
	  if(CLevo=="yes"){
	    if(mclcorr*mhalf > 1.E1)
	      tdf *= (pow(rclcorr,3./2.-1.74)*sqrt(mclcorr));
	    else
	      tdf  = 1.E12;
	    
	  }
	  
	  nan.str("");
	  nan<<tdf;
	  if(nan.str()=="nan"||nan.str()=="-nan"){
	    cout<<"DF ERROR: "<<tdf<<endl;
	    exit(0);
	  }
	  trecy += tdf;

	  //Time for BHs to pair --> assume it's already hard
	  if(CLevo == "yes"){
	  
	    /*mclcorr = exp(-(time+trecy-tfor[i])/ (0.45*trelax));
	      rclcorr = (1. + log10(1. + (time+trecy-tfor[i])/(0.45*trelax))); //REVISE THIS*/

	    mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	    rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax, CLfill);

	    sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	    rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	    rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	    
	    
	    vthre = vthre_in * sqrt(mclcorr/rclcorr);

	    semihard = 1./pow(sig_clu/30., 2.); //2.*6.67E-11*1.99E30/(1.E6*1.5E11) * (mpri+msec) / (sig_clu*sig_clu);

	    
	    t12 = 3.E8 / zita * (0.01/fb) * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(mstar / (mpri+msec)) * (30./(mpri+msec+mper)) * (1./semihard);
	    t12 *= func.rndgen(1.0, 0.1);

	    
	    if(time + trecy < tcc && time + trecy + t12 > tcc){
	      time = tcc;
	      
	      mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	      rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax, CLfill);

	      sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	      rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	      rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	      
	      
	      vthre = vthre_in * sqrt(mclcorr/rclcorr);
	    }
	    
	    
	  }
	  
	  
	  semihard = 1./pow(sig_clu/30., 2.); //2.*6.67E-11*1.99E30/(1.E6*1.5E11) * (mpri+msec) / (sig_clu*sig_clu);

	  
	  t12 = 3.E8 / zita * (0.01/fb) * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(mstar / (mpri+msec)) * (30./(mpri+msec+mper)) * (1./semihard);
	  t12 *= func.rndgen(1.0, 0.1);

	  if(mclcorr*mhalf < 1.E1 && CLevo == "yes")
	    t12 = 1.E12;
	    
	  
	  trecy += t12;

	  if(CLevo == "yes"){
	    
	    /*mclcorr = exp(-(time+trecy-tfor[i])/ (0.45*trelax));
	      rclcorr = (1. + log10(1. + (time+trecy-tfor[i])/(0.45*trelax))); //REVISE THIS*/

	    mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	    rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax, CLfill);

	    sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	    rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	    rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	    
	    
	    vthre = vthre_in * sqrt(mclcorr/rclcorr);
	    
	  }
	  
	  
	  //calculate aej vs agw to merge inside or outside --> inside (tgwpeters), outside (5/q3 * t2-1)
	  semi_ej = 0.07 * (mu_red * mper / (mpri + msec + mper) * mper/(mpri+msec)) * pow(vthre/50.,-2.);
	  semi_gw = 0.05 * pow((mpri+msec)/20,3./5.)*pow( (msec/mpri) / pow(1+(msec/mpri),2.) , 1./5.) * pow(sig_clu/30., 1./5.) * pow(1.E6 / rho_clu,1./5.);

	  acrit = max(semi_ej,semi_gw);
	  
	  tbbh = 2.E7/zita * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(10. * mstar / (mpri+msec)) * (0.05 / acrit) * (20./(mpri+msec));
	  tbbh *= func.rndgen(1.0, 0.1);	  
	  
	  tmer = 5. * (mpri + msec)/mper * tbbh;						 	  

	  if(mclcorr*mhalf < 1.E1){
	    semi_ej = -2.0;
	    semi_gw = -1.0;
	    acrit   = 1.E12;
	    tmer    = 1.E12;
	  }
	  
	  trecy+= tmer;
	  
	  if(CLevo == "yes"){
	    
	    /*mclcorr = exp(-(time+trecy-tfor[i])/ (0.45*trelax));
	      rclcorr = (1. + log10(1. + (time+trecy-tfor[i])/(0.45*trelax))); //REVISE THIS*/

	    mclcorr = func.mevol(time-tfor[i], rhalf, mhalf, trelax, CLfill);
	    rclcorr = func.revol(time-tfor[i], rhalf, mhalf, trelax, CLfill);

	    
	    sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
	    rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
	    rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
	    
	    
	    vthre = vthre_in * sqrt(mclcorr/rclcorr);
	    
	  }
	 	  
	  if(semi_ej < 0.0 && semi_gw < 0.0)
	    label="dissolved";
	  else if(semi_ej < semi_gw)	    
	    label="inside";
	  else
	    label="ejected";


	  ecc = sqrt(func.rnd());
	 	  
	  double clight = 3.E8;
	  double Grav   = 6.67E-11;
	  double Msun   = 1.99E30;
	  double Rau    = 1.5E11;	    
	  tmer = 5./256. * pow(clight,5.) / pow(Grav * Msun, 3.) * pow(Rau,4.)  * pow(acrit,4.)/(mpri*msec*(mpri+msec)) * pow(1.-pow(ecc,2.),3.5) / (1. + 73./24. * ecc*ecc + 37./96.*ecc*ecc*ecc*ecc) * 1./ (365.*24.*3600.*1.E0);

	  if(mclcorr*mhalf < 1.E1)
	    tmer = 1.E12;
	  
	  trecy += tmer;

	  nan.str("");
	  nan<<time;
	  if(nan.str()=="-nan" || nan.str()=="nan"){
	    cout<<"Critical error Jesus "<<tmer<<" "<<tbbh<<" "<<tdf<<" "<<t12<<" "<<time<<" "<<trecy<<endl;
	    exit(0);
	  }
	  
	  Spinning[0] = 0;
	  Spinning[1] = 0;
	  Spinning[2] = 0;
	  Spinning[3] = 0;
	  if(mpri>0.0 && msec>0.0)
	    func.SREM2(ndx, apri, asec, mpri, msec, align, Spinning);	  	 


	  if(mpri > 200.)
	    sig_clu = sqrt(sig_clu0*sig_clu0 + mpri / (0.1*pow(10.,rint)));
	  
	  nrecy += 1;
	  	  
	}while(1>0);

	
	if(CLevo == "yes"){
	  
	  mclcorr = func.mevol(Hubble-tfor[i], rhalf, mhalf, trelax, CLfill);
	  rclcorr = func.revol(Hubble-tfor[i], rhalf, mhalf, trelax, CLfill);
	  
	}
	

	
	
	int itot = i;

	double zmer = func.inter(time / 1.E9, age, reds, redline);
	double zfor = func.inter(tfor[i] / 1.E9, age, reds, redline);	
	double zsmbh= func.inter(tsmbh/1.E9, age, reds, redline);
	
	if(time < Hubble){
	  out<<itot<<" "<<Z[i]<<" "<<nrecy<<" "<<cluster<<" "<<REC<<" "<<mpri<<" "<<msec<<" "<<apri<<" "<<asec<<" "<<Mrem[i]<<" "<<Srem[i]<<" "<<Xrem[i]<<" "<<Krem[i]<<" "<<tfor[i]<<" "<<time<<" ";
	  out<<pow(10.,mint)<<" "<<pow(10.,rint)<<" "<<vthre<<" "<<label<<" "<<semi_ej<<" "<<semi_gw<<" "<<nbhs<<" "<<mhalf*mclcorr<<" "<<rhalf*rclcorr<<" "<<zmer<<" "<<zfor<<" "<<tsmbh<<" "<<zsmbh<<endl;
	  
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

    
    
    
  }
  clout.close();

  out.close();
  out2.close();
  out3.close();
  hout.close();
  
  
  vector<double> mpost;

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

  double *X;
  X = new double [mpost.size()];
  for(int i=0;i<mpost.size();i++)X[i]=mpost[i];
  func.histo(X,mpost.size(),30,"linear","Many_catalogues.txt");
  delete [] X;


  stringstream Fcl;
  
  stringstream ggc;  stringstream gyc;  stringstream gnc;
  stringstream qgc;  stringstream qyc;  stringstream qnc;
  stringstream fgc;  stringstream fyc;  stringstream fnc;
  stringstream alg;
  

  Fcl<<DynOvTot;
  ggc<<probGC;  gyc<<probYC;  gnc<<probNC;
  qgc<<MinQGC;  qyc<<MinQYC;  qnc<<MinQNC;
  fgc<<recyGC;  fyc<<recyYC;  fnc<<recyNC;
  alg<<INDEX_ALIGN;

  string ZDIS = zdist;
  string ZDYN = zdyn;
  
  stringstream metaldivi;
  metaldivi<<zdivi;

  string corr = correction;
  string type_mrat = MRATIO;

  //FINAL FILE MOVING
  string cmdstr_zero =  "./SIM_Fdyn"+Fcl.str()+"_Ngc"+ggc.str()+"_GCq_"+qgc.str()+"_GCf"+fgc.str()+"_Nyc"+gyc.str()+"_YCq"+qyc.str()+"_YCf"+fyc.str()+"_Nnc"+gnc.str()+"_NCq"+qnc.str()+"_NCf"+fnc.str()+"isolS_"+isolS+"dynaS_"+dynaS+"_"+"MetalDivi_"+metaldivi.str()+"_"+alg.str()+"_"+ZDIS+"_"+ZDYN+"_Correction_"+corr;  
  if(kick=="yes")
    cmdstr_zero += "_kick_Yes";
  else
    cmdstr_zero += "_kick_No";
  
  if(crazy=="yes")
    cmdstr_zero += "_crazy_Yes";    
  else
    cmdstr_zero += "_crazy_No";    

  cmdstr_zero += "_mratio";
  cmdstr_zero += type_mrat;
  if(type_mrat=="pwl"){
    stringstream mrslp_str;
    mrslp_str<<MRATIO_SLOPE;
    cmdstr_zero +=  mrslp_str.str();
  }
  
  if(delaytime=="yes")
    cmdstr_zero += "_delaytimes_Yes";
  else
    cmdstr_zero += "_delaytimes_No";
  
  cmdstr_zero += "_primslope_";
  stringstream msl;
  msl<<obslope;
  cmdstr_zero += msl.str();


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
  
  stringstream smix;
  smix<<mixing;

  cmdstr_zero += "_mixing_"+smix.str();

  if(MonoZ == "yes")
    cmdstr_zero += "_MonoZ";
  if(cluster_test == "yes")
    if(cluster_test_env=="NC")
      cmdstr_zero += "_clteNC";
    else      
      cmdstr_zero += "_cltest";
  if(CLevo == "yes"){
    if(CLfill == "under")
      cmdstr_zero += "U";
    else if(CLfill == "filling") 
      cmdstr_zero += "F";
    else if(CLfill == "postcol")
      cmdstr_zero += "P";
    else if(CLfill == "nuclear")
      cmdstr_zero += "N";
    else if(CLfill == "nuclearF")
      cmdstr_zero += "NF";
  }
  
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
  delete [] Z;
  delete [] Zi;


  return 0;
  
  
  }

