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
#include <chrono>
#include <algorithm>
#include "funcGR.h"

#define CLfill       "GG23"

using namespace std;

double extM(){
  Functions func;
  double mslp = -2.0;
  double P = func.rnd();
  double m1 = pow((P*pow(80.,1.+mslp) + (1.-P)*pow(5.0,1.+mslp)),1./(1.+mslp));

  return m1;
}

void singBHt_mix(vector<double>& zams_mix,
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

double GWeff(string pcluster, double met){
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

void hgen(double eps, double m1, double a1, double m2, double a2, double vesc, string stype, // star variables
          vector<double>& zams_mix, vector<double>& remn_mix, vector<double>& tdel_mix, vector<double>& kick_mix,// catalog variables
          double *c, double *s, double nbhs, double nrecy, double nmerg, int id, // storage vectors
          double mhalf, double mcore, double rcore, double n_bin, vector<double>& gwK, vector<double>& gwK_cdf, // core properties
          double trelax, double t12capt, double tbbhform, double tcc, string pcluster){ //timescales
  //What about the probability for a higher gen merger to occur? interaction rate?
  //What about the delay time for high-gen merger to develop?
  

  Functions func;
  double tau, interaction_rate;
  double mstar_avg;
  int cnt=0;
  
  string success="no";

  // Let's define the timescale for a secondary merger
  // mean value of Tau = tGW/tcc = alpha * log10(Mcl/Ms) + betha + 1.5 dex gaussian scatter
  tau = -0.53 * log10(mhalf) + 5.6 + 1.5* func.rndgen(0.0, 1.0);
  
  // Let's define the quantities intervining in the interaction rate
  mstar_avg = 17.4 - 4.0 * log10(tbbhform/trelax);
  double m_bin = m1 + m2;
  
  //Let's infer the retention fraction of hierarchical from this cluster
  double ret_fract;
  int i=0;

  do{
    ret_fract = gwK_cdf[i];
    i++;
    //cout << "ret_fract: " << ret_fract << endl;

    }while(gwK[i]<vesc);
  
  if(nmerg < 1){
    double P = func.rnd();
    if(P < nmerg) //We have more than one merger, we likely end mergers here
      success = "yes";       
  }
  // If I have >1 merger in the cluster AND the timescale is enough to have the formation of a merging BBH
  // I can have a hierarchical secondary BH
  else if(nmerg >= 1){
    success = "no";
    // Cristiano 21/05/2025:
    // Changed to a while loop:
    // If the time condition is not satisfied (i.e., not enough time for a BBH merger)
    // OR I exhausted the BH reservoir (i.e., I have no more BHs to merge)
    // I cannot have a hierarchical merger

    interaction_rate = 0.0; //in this way I can check if the while loop is verified at least once

    while (cnt < nmerg && max(tbbhform/tcc, t12capt/tcc) >= tau) {
      success = "yes";
      // Let's initialize the star merging with the secondary BH
      double *single_bh;
      single_bh = new double [4];
      single_bh[0] = 0.0;
      single_bh[1] = -1.0;
      single_bh[2] = 0.0;
      single_bh[3] = 0.0;

      // Let's compute the 0-g BH properties
      singBHt_mix(zams_mix, remn_mix, tdel_mix, kick_mix, single_bh, vesc);
      double m2b = single_bh[0];	 
      double a2b = func.spin(m2b, stype);

      //Let's compute POTENTIAL the merger remnant and its natal kick
      func.SREM2(8.0, a2, a2b, m2, m2b, "dynamical", s);

      // Let's compute the interaction rate
      // rho_star = M_core/R_core^3 ==> n_star = rho_star/mstar_avg
      // R_hier = max(R_inf, R_wand) 
      // with the wandering radius R_wand being the region of the core where we can find the BH
      // and the influence radius being the region of the vore where the BHs dominates the interactions
      //                   | alpha, with alpha = m_hier/m_core
      // R_hier = R_core x |
      //                   | mu exp[-0.5], with mu = m_hier/mstar_avg

      // Thus, we define:
      // rho_hier = n_hier * m_hier/(R_hier)^3 from Di Cinto et al. 2023
      // with m_hier = max( 2 * mstar_avg, s[2] )
      // in this way we account for the average mass of the 1-g BHs in the cluster for low gens and for the higher gen mergers
      
      //double m_hier = max(2*mstar_avg, s[2]);
      double m_hier = 2*mstar_avg;

      double alpha =  m_hier/mcore; 
      double mu = m_hier/mstar_avg;

      double r_inf = rcore * alpha;
      double r_wand = rcore * pow(mu, -0.5);

      double r_hier = max(r_inf, r_wand);

      double rho_star = mcore / (mstar_avg * pow(rcore, 3));

      // for the 1-g BH, we assume to generate n_hier = n_bin * nBHs/2
      // Thus, rho_hier = n_hier * m_hier / V_hier 
      // We need to fix this number with the fraction of retained 1-g BHs, as GW recoils will eject part of them
      // but in the interaction rate formula we have rho/m
      // Thus, we can write:
      double n_ret = ret_fract * pow((n_bin /2 ) * nbhs, cnt+1);
      double n_hier = max(1.0, n_ret);
      // We are also considering the retention fraction by multiplying for the retention fraction
      
      double rho_hier = n_hier / pow(r_hier, 3);
      
      // Now we can compute the interaction rate
      interaction_rate = (rho_hier / rho_star) * pow(s[2] / mstar_avg, 3. / 2.) * (m_bin + s[2]) / (m_bin + mstar_avg);
      
      // if (interaction_rate < 0){
      //  // Print out the components for debugging the interaction rate calculation:
      //   cout << "Interaction rate <0 details:" << endl;
      //   cout << "mcore: " << mcore << endl;
      //   cout << "r_core: " << rcore << endl;
      //   cout << "r_inf: " << r_inf << endl;
      //   cout << "r_wand: " << r_wand << endl;
      //   cout << "rho_hier: " << rho_hier << endl;
      //   cout << "rho_star: " << rho_star << endl;
      //   cout << "m2: " << m2 << endl;
      //   cout << "s[2]: " << s[2] << endl;
      //   cout << "mstar_avg: " << mstar_avg << endl;
      //   cout << "m_hier: " << m_hier << endl;
      //   cout << "m_bin: " << m_bin << endl;
      //   cout << "n_bhs: " << nbhs << endl;
      //   cout << "retention fraction: " << ret_fract << endl;
      //   cout << "n_hiers: " << n_hier << endl;
      //   cout << "interaction_rate: " << interaction_rate << endl;
      // }

      // Now I trow a dice, if the interaction rate is high enough, I can have a merger
      double dice = func.rnd();
      //If the interaction rate is high enough, we find a companion BH coming from a hierarchical merger

      if(dice > interaction_rate){
        // If the interaction rate is low, I keep the stellar secondary
        //cout << "interaction rate too low => 0-th gen secondary" << endl;
        break;
      }
      
      //If the 1g BH is not expelled from the cluster it pairs with the primary BH
      double vrec = s[3];

      if(vrec > vesc){
        //cout << " ejected => 0-th gen secondary - vrec: " << vrec << " vesc: " << vesc << endl;
        break;
      }
      
      // cout << "#####################################" << endl;
      // cout << "ID: " << id << " interacted" << endl;
      // cout <<  "Interaction rate: " << interaction_rate << " dice: " << dice << endl;
      // cout << "#####################################" << endl;

      // I can have a merger and the new secondary is the result of the merger of 
      // the 0-g stellar BH and the secondary stellar BH giving birth to a 1-g secondary BH
      m2 = s[2];
      a2 = s[0];

      // We have one BBH merger
      cnt++;

      //We need to adjust the timescale for the next merger
      tau += -0.53 * log10(mhalf) + 5.6 + 1.5* func.rndgen(0.0, 1.0); //See notion notes/plot

      break; // I can have only one hierarchical merger per generation

    }
    
  }
  // If I have no hierarchical merger, I can still have a 0-g merger
  
  c[0] = m2;
  c[1] = a2;
  c[2] = cnt;
  c[3] = interaction_rate;

  //cout << "hgen output for the " << cnt << "-th generation" << endl;
  //cout << "ID: " << id << " m2: " << c[0] << " a2: " << c[1] << " cnt: " << c[2] << " interaction_rate: " << c[3] << endl;
  
  return ;
  
}



int main(){
  srand(time(0));
  
  using clock = std::chrono::system_clock;
  using sec = std::chrono::duration<double>;
  const auto before = clock::now();

  int ssize = 7;
  double *Spinning;
  Spinning = new double [ssize];
  
  Functions func;

  double rhalf, mhalf;
  double t12capt, t3bb, tbbhform, tbbh, tmer, tau;
  double interaction_rate;
  double m_core, r_core;
  double vthre, vthre_in;
  double sig_clu, sig_clu_in;
  double rho_clu, rho_clu_in;
  double rho_cubicpc, rho_cubicpc_in;
  double mclcorr, rclcorr;
  double mstar, mper;
  double m1, a1, m2, a2;

  // Let's define the parameters for the cluster
  double Z = 0.0003;
  string pcluster = "nuclear";

  //Here we need to set a specific binary fraction per cluster type (and perhaps time evolution?)
  double fb; 	
  if(pcluster == "young")
    fb = 0.4;
  else if(pcluster == "globular")
    fb = 0.1;
  else if(pcluster == "nuclear")
    fb = 0.01;

  // Let's draw the catalog of initial masses of the binary components
  
  // Let's input the catalogs
  vector<double> zams_mix, remn_mix;
  vector<double> tdel_mix, kick_mix;
  vector<double> gw_recoil, gw_recoil_cdf;
  vector<double> gw_recoil_hg, gw_recoil_hg_cdf;

  string cat_name = "spectrum0003.txt";
  ifstream in;
  in.open(cat_name.c_str());
  double par[20];
  int spar = 5; // number of parameters in the catalog
  
  do{
    double par[spar];

    for(int jj=0;jj<spar;jj++)
      in>>par[jj];
          
      zams_mix.push_back(par[0]);
      remn_mix.push_back(par[1]);
      tdel_mix.push_back(1.E6 * par[3]);
      kick_mix.push_back(par[4]);

    }while(!in.eof());

    in.close();

  // Let's read the GW kick cdf
  string gw_cat_name = "kick_velocity_cdf.csv";
  in.open(gw_cat_name.c_str());
  //cout << "reading kicks" << endl;

  int gw_kick_col=2;

  do{
    double bpar[spar];

    for(int jj=0; jj<gw_kick_col; jj++) in>>par[jj];
      gw_recoil.push_back(par[0]);
      gw_recoil_cdf.push_back(par[1]);
    
    }while (!in.eof());

  in.close();

  // Let's read the GW kick cdf for generation above the first
  string gw_hg_cat_name = "kick_velocity_hg_cdf.csv";
  in.open(gw_cat_name.c_str());
  //cout << "reading kicks" << endl;

  do{
    double bpar[spar];

    for(int jj=0; jj<gw_kick_col; jj++) in>>par[jj];
      gw_recoil_hg.push_back(par[0]);
      gw_recoil_hg_cdf.push_back(par[1]);
    
    }while (!in.eof());

  in.close();
  
  //cout << "read kicks" << endl;
  
  ofstream outfile("output_hg.txt");
  //outfile << "i mhalf rhalf mcore m1 m2 a1 a2 interaction_rate number" << endl;
  outfile << "ID time mhalf rhalf mcore m1 m2 a1 a2 interaction_rate nhier N" << endl;
  
  for(int ID=0; ID<100000; ID++){
    
    // Cluster scale radius (from Dehnen) & initial properties //
	  double g_cl, mint, mean;
	  if(pcluster == "young"){
	    g_cl = 1.0;
      mint= 3.0 + 2.0*func.rnd();
      mean = log(pow(10.,3.3));
    }
	  else if(pcluster == "globular"){
	    g_cl = 1.5;
      mint = 4.0 + 3.0*func.rnd();
      mean = log(pow(10.,3.7));
  }
	  else if(pcluster == "nuclear"){
	    g_cl = 1.9;
      mint = 6.0+4.0*func.rnd(); 
      mean = log(pow(10.,5.0));
  }
  
    double sigma = log(pow(10.,0.4));
    double rhoint = func.LogGaussian(mean, sigma);
    rhoint = log10(rhoint);
    
    double rint = (mint - rhoint)/3.0;
    
    double Nrecy = 0.0;
    double Nbhs;
    double zita = 1.0;
    
    //cout << "Cicle: "<< ID << endl;
    //cout << "Cicle: "<< i << " mhalf: " << mhalf <<  " rhalf: " << rhalf << endl;
    
    double eps = GWeff(pcluster,Z);

    double nmerg = min(Nbhs-Nrecy, eps*mhalf);
    // cout << "nrecy: " << Nrecy << "; "
    //     << "nbhs: " << Nbhs << "; "
    //     << "eps: " << eps << "; "
    //     << "Mc: " << mhalf << "; "
    //     << "eps*Mc: " << eps*mhalf << "; "
    //     << "nmerg: " << nmerg << "\n"
    //     << endl;

    // Let's define the velocity dispersion and the semimajor axis for an hard binary
    // scaling relation in paper 2021 BHNS mergers (Arca Sedda 2020)
    // relazione tra massa, raggio e sigma di nuclear clusters  --> empirica

    double lsig = 0.5 * (-1.14 + log10(6.67E-11*1.99E30/3.08E16 * pow(10.,mint)/pow(10.,rint) ));			 
    sig_clu = 0.001*pow(10.,lsig);
    
    rho_clu = 3. *pow(10.,mint) / ( 4. * M_PI * pow(pow(10.,rint),3.));
    rho_clu_in = rho_clu;
        
    // mass of the stars aka the perturbers
    mstar = 1.0;
    mper = 15.0;
    rho_cubicpc = rho_clu / mstar;

    // Let's compute the relaxation time and the core-collapse time
    double trelax =  0.78E9 /log(0.11 * mhalf) * pow(mhalf/1.E5, 0.5) * pow(rhalf, 1.5);
    double trelax0 = trelax;
    double tcc = 0.138*mhalf/(150. * log(0.11*mhalf/150.))*sqrt(pow(rhalf*3.08E16,3.)/(6.67E-11*1.99E30*mhalf))/(365.*24.*3600.);

    vthre = func.vescape(pow(10.,rint),pow(10.,mint),pcluster);

    // The values defined up to now are the init cond of the cluster
    vthre_in = vthre;
    sig_clu_in = sig_clu;
    rho_clu_in = rho_clu;
    rho_cubicpc_in = rho_cubicpc;

    // Let's initialize the primary 0-gen (i.e. not hierarchically formed) 
    // string spintype = "maxwellian02";
    // m1 = extM();
    // a1 = func.spin(m1, spintype);

    // Catalog implementation
    double *primary_bh;
    primary_bh = new double [4];
    primary_bh[0] = 0.0;
    primary_bh[1] = -1.0;
    primary_bh[2] = 0.0;
    primary_bh[3] = 0.0;


    double *secondary_bh;
    secondary_bh = new double [4];
    secondary_bh[0] = 0.0;
    secondary_bh[1] = -1.0;
    secondary_bh[2] = 0.0;
    secondary_bh[3] = 0.0;

    singBHt_mix(zams_mix, remn_mix, tdel_mix, kick_mix, primary_bh, vthre);
    singBHt_mix(zams_mix, remn_mix, tdel_mix, kick_mix, secondary_bh, vthre);

    double tpri, tsec, kpri, ksec;
    if(primary_bh[0]>secondary_bh[0]){
     
      m1 = primary_bh[0];	 
      m2 = secondary_bh[0];	 
      tpri = primary_bh[1];
      kpri = primary_bh[2];
      tsec = secondary_bh[1];
      ksec = secondary_bh[2];
    }
    else{
      m1 = secondary_bh[0];	 
      m2 = primary_bh[0];	 
      tpri = secondary_bh[1];
      kpri = secondary_bh[2];
      tsec = primary_bh[1];
      ksec = primary_bh[2];
    }
    
    // Let's define the spin type for the primary and secondary BHs
    // lines from BPOP
    string spintype = "maxwellian02"; 
    a1   = func.spin(m2,spintype);
	  a2   = func.spin(m2,spintype);	


    for(int i=0; i<ssize; i++) Spinning[i] = 0.0;
    Spinning[3] = -1.E30;

    int elem = 4;
    double *Comp;
    Comp = new double [elem];
    for(int i=0;i<elem;i++) Comp[i] = 0.0;
    
    // Let's compute the number of BHs in the cluster
    double a_cl = rhalf*(pow(2.,1./(3.-g_cl))-1);

	  // DOUBLE CHECK THE FOLLOWING //

    // Max radius from which BHs can spiral-in over a Hubble time via DF //	  
    double Hubble = 13.99E9;
	  double radius = rhalf * pow( Hubble / (0.42E9 * (10.*mstar/(m1+m2)) * (trelax0 / 4.2E9)) , 1./1.74);
	  
    // Fraction of mass enclosed within the infall radius above //
	  double fencl = (radius / (radius + a_cl),3.-g_cl) * (1. + 0.2*(1.-2.*func.rnd()));

	  // retention fraction freten //
	  double freten = 0.5 * (1. + 0.3*(1.-2.*func.rnd()));

	  // Fraction number of BHs in a power-law IMF between 0.08 and 150 Msun//
	  double fraBH = 0.0008 * (1. + 0.1*(1.-2.*func.rnd()));
	  
	  // This depends on the number fraction of BHs in the cluster, we're also assuming mint == N_* 
	  //Nbhs = fraBH * pow(10.,mint) * freten * fencl;


    //Nbhs = 10*4 * func.rndgen(0.0, 1.0);


    rhalf = pow(10.,rint);
    mhalf = pow(10.,mint);
    
    Nbhs = 0.01*pow(10.,mint)/30.;
    //cout << "ID: " << ID << " nBHs: " << Nbhs << " mhalf: " << mhalf << endl;

    tbbhform = 0.0; // Let's initialize the time for the BBH formation
    double t_star = max(tpri, tsec); //minimum time for BHs
    double time = max(t_star, tcc); // Let's initialize the time of the cluster evolution
    
    while(time < 13.5E9 && Nbhs>0) { // Let's evolve the cluster until 13.5 Gyr
      // Let's simulate the evolution of the cluster
  
      if (nmerg == 0) break;

      //let's recompute them after the evolution of the cluster
      mclcorr = func.mevol(time, rhalf, mhalf, trelax0, CLfill, pcluster);
      rclcorr = func.revol(time, rhalf, mhalf, trelax0, CLfill, pcluster);
            
      sig_clu    = sig_clu_in*sqrt(mclcorr/rclcorr);
      rho_clu    = rho_clu_in*mclcorr/pow(rclcorr,3.);
      rho_cubicpc = rho_cubicpc_in*mclcorr/pow(rclcorr,3.);
      
      vthre = vthre_in * sqrt(mclcorr/rclcorr);

      m_core = mclcorr*mhalf;
      r_core = rclcorr*rhalf;

      //cout<< "ID:" << ID << " time: " << time << " m_core: " << m_core << endl;
      
      // Let's compute the binary formation timescales

      double semihard = 1./pow(sig_clu/30., 2.); // Antonini&Rasio 2016

      t12capt = 3.E9 * (0.01 / fb) * (1.E6 / rho_cubicpc) * (sig_clu / 30.) * 10./(m1 + m2 + mper) * 1./semihard;
      t3bb = 4.E9 * pow(1.E6 / rho_clu,2.) * pow(sig_clu/(zita * 30.),9.) * pow(mstar / 30. * 10., 9./2.) * pow(10./30.,-5.);
      
      tbbhform = min(t12capt,t3bb);
      tbbhform *= func.rndgen(1.0, 0.1); //resampling gaussiano intorno al valore medio
      
      time += tbbhform;
      hgen(eps, m1, a1, m2, a2, vthre, spintype,zams_mix, remn_mix, tdel_mix, kick_mix, Comp, Spinning, Nbhs, Nrecy, nmerg, ID, mhalf, m_core, r_core, fb, gw_recoil, gw_recoil_cdf, trelax0, t12capt, time, tcc, pcluster);
        
      m2 = Comp[0];
      a2 = Comp[1];
      int ng = int(Comp[2]);
      double interaction_rate = Comp[3];
      
      nmerg -= ng; // Let's adjust the number of mergers
      //cout << "BH reservoir: " << nmerg << endl;

      // cout << "routine output for the " << i << "-th iteration" << endl;
      // cout << "m2: " << Comp[0] << " a2: " << Comp[1] << " cnt: " << Comp[2] << " interaction_rate: " << Comp[3] << endl;


      //outfile << ID << " " << time << " " << mhalf << " " << rhalf << " " << m_core << " "<< m1 << " " << m2 << " " << a1 << " " << a2 << " " << interaction_rate << " " << ng << " "<< Nrecy << endl;
      

      //cout<<"m1"<<" "<<"m2"<<" "<<"a1"<<" "<<"a2"<<" "<<"number"<<" "<<endl;
      //cout<<m1<<" "<<m2<<" "<<a1<<" "<<a2<<" "<<ng<<" "<<endl;

      
      for(int i=0; i<ssize; i++) Spinning[i] = 0.0;
      Spinning[3] = -1.E30;


      tbbh = 2.E7/zita * (1.E6/rho_cubicpc) * (sig_clu/30.) * sqrt(10. * mstar / (m1+m2)) * (0.05 / semihard) * (20./(m1+m2));
	    tbbh *= func.rndgen(1.0, 0.1);	  
      tmer = 5. * (m1+m2)/mper * tbbh;						 	  

      time += tmer;
      //cout <<"ID: "<< ID << " time: " << time << " tbbh: " << tbbh << " tmer: " << tmer << endl;
      
      func.SREM2(8.0, a1, a2, m1, m2, "dynamical", Spinning);
      Nrecy ++;
      outfile << ID << " " << time << " " << mhalf << " " << rhalf << " " << m_core << " "<< m1 << " " << m2 << " " << a1 << " " << a2 << " " << interaction_rate << " " << ng << " "<< Nrecy << endl;
      m1 = Spinning[2];
      a1 = Spinning[0];

      if(Spinning[3] > vthre) {
        //cout << "Primary ejected => end of growth - kick: " << Spinning[3] << " vesc: " << vthre << endl;
        break; // If the primary is ejected, we stop the evolution
      }
      else {
        //cout << "not ejected" << endl;
        Nbhs--;
        nmerg --;
        //cout << "Nbhs: " << Nbhs << " Nrecy: " << Nrecy << endl;
        
      }
       
    }
    
  }
  outfile.close();

  const sec duration = clock::now() - before;

  
  cout<<"Elapsed time = "<<duration.count()<<"s"<<endl;

  return 0;
}

  
