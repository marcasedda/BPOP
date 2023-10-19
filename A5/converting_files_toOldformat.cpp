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
#include <ctime>
#include <stdlib.h>

#define MODEL "0"

using namespace std;

int main(){

  string model = MODEL;
  
  char *nmodel;
  nmodel = new char [model.length()+1];

  strcpy(nmodel, model.c_str());
  
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


  for(int i=0;i<13;i++){

    stringstream ss;
    ss<<met[i];
    
    string fname = "DATI_ORIGINALI/data_BBHs_"+ss.str()+".txt";
    ifstream in;
    in.open(fname);
    if(!in.is_open()){
      cout<<"File "<<fname<<" not found"<<endl;
      break;
    }
    string line;
    getline(in,line);
    getline(in,line);
    getline(in,line);
    int cnl = 0;
    do{
      getline(in,line);
      cnl += 1;
    }while(!in.eof());
    in.close();
    
    string oname = "data_"+ss.str()+".dat";
    ofstream out;
    out.open(oname);
    in.open(fname);
    getline(in,line);
    getline(in,line);
    getline(in,line);
    int flag = 0;
    for(int j=0;j<cnl-1;j++){
      int kpar = 20;
      string dmy1;
      double par[kpar];
      double m1, m2, mbh1, mbh2, sma, ecc, delay;
      in>>dmy1;
      for(int k=0;k<kpar;k++)
	in>>par[k];


      
      out<<par[2] + par[3]<<" "<<par[2]<<" "<<par[3]<<" "<<par[0]<<" "<<par[1]<<" "<<par[4]<<" "<<par[5]<<" "<<par[6];
      for(int k=0;k<kpar - 8;k++)
	out<<" 0.0 ";
      out<<endl;
      flag = 1;
      
    }
    in.close();
    out.close();

    vector<double> zams_new,remn_new;
    vector<double> zams_old,remn_old;
    vector<double> time_new, kick_new;
    vector<int>    type_old;
    
    fname = "DATI_ORIGINALI/single_BH_"+ss.str()+".txt";
    in.open(fname);
    if(!in.is_open()){
      cout<<"File "<<fname<<" not found"<<endl;
      break;
    }
    getline(in,line);
    getline(in,line);
    getline(in,line);
    cnl = 0;
    do{
      getline(in,line);
      cnl += 1;
    }while(!in.eof());
    in.close();

    string lol = ss.str();
    char *buf;
    buf = new char [lol.length()+1];
    strcpy(buf, lol.c_str());

    stringstream sss;
    for(int j=2;j<lol.length();j++){
      sss<<buf[j];
    }

    flag = 0;
    oname = "spectrum"+sss.str()+".txt";
    out.open(oname);
    in.open(fname);
    getline(in,line);
    getline(in,line);
    getline(in,line);
    for(int j=0;j<cnl-1;j++){
      int kpar = 5;
      string dmy1;
      double par[kpar];
      in>>dmy1;
      for(int k=0;k<kpar;k++)
	in>>par[k];

      
      flag = 1;
      out<<par[0]<<" "<<par[1]<<" "<<14<<" "<<par[2]<<" "<<par[3];
      out<<endl;
      
      zams_new.push_back(par[0]);
      remn_new.push_back(par[1]);
      time_new.push_back(par[2]);
      kick_new.push_back(par[3]);
      
    } 
    in.close();
    out.close();


    delete [] buf;

    
    oname = "../DATI_SingleBH/spectrum"+sss.str()+".txt";
    in.open(oname);
    getline(in,line);
    do{
      int kpar = 3;
      double par[kpar];
      for(int k=0;k<kpar;k++)
	in>>par[k];

      zams_old.push_back(par[0]);
      remn_old.push_back(par[1]);
      type_old.push_back(par[2]);
      
    }while(!in.eof());
    in.close();

    vector<double> ratios;
    for(int j=0;j<zams_old.size()-2;j++){
      ratios.push_back(remn_old[j] / zams_old[j]);
    }
    
    oname = "spectrum_full"+sss.str()+".txt";
    ofstream outu;
    outu.open(oname);

    vector<double> tform,nform, vform;
    for(int j=0;j<zams_old.size();j++){
      tform.push_back(0.0);
      vform.push_back(0.0);
      nform.push_back(0.0);
    }
    
    for(int j=0;j<zams_new.size();j++){
      int id1 = -1;
      int id2 = -1;
      for(int k = 0; k<zams_old.size()-2; k++){
	if(zams_new[j] > zams_old[k])
	  id1 = k;
	if(zams_new[j] < zams_old[k]){
	  id2 = k;
	  break;
	}

	if(id1 == zams_old.size()-1){
	  id1 = zams_old.size()-2;
	  id2 = zams_old.size()-1;
	}
	if(id2 == -1){
	  id2 = id1 + 1;
	}
	if(id1 == -1){
	  id1+= 1;
	  id2+= 2;	    
	}
	
      }
      
           
      if( abs(remn_new[j] / zams_new[j] - ratios[id1]) / ratios[id1] < 0.02 ||
	  abs(remn_new[j] / zams_new[j] - ratios[id2]) / ratios[id2] < 0.02 ){

	tform[id1] += time_new[j];
	tform[id2] += time_new[j];
	vform[id1] += kick_new[j];
	vform[id2] += kick_new[j];
	nform[id1] += 1;
	nform[id2] += 1;

	outu<<zams_new[j]<<" "<<remn_new[j]<<" "<<14<<" "<<time_new[j]<<" "<<kick_new[j]<<endl;
	
      }
      
      
    }

    outu.close();
    
    oname = "spectrum_cleaned"+sss.str()+".txt";
    outu.open(oname);
    for(int i=0;i<remn_old.size();i++){
      if(type_old[i] < 14)
	continue;
           
      
      double taver = tform[i]/nform[i];
      double kaver = vform[i]/nform[i];
      if(nform[i] <= 0.0){
	taver = 0.0;
	kaver = 0.0;	
      }
      
      stringstream pd;
      pd<<taver;
      if(pd.str() == "nan" || pd.str() == "-nan"){
	cout<<pd.str()<<" "<<tform[i]<<" "<<nform[i]<<endl;
	exit(0);
      }
      if(nform[i] > 0.0)
	outu<<zams_old[i]<<" "<<remn_old[i]<<" "<<type_old[i]<<" "<<taver<<" "<<kaver<<endl;
    }
    outu.close();
    zams_old.erase(zams_old.begin(),zams_old.end());
    remn_old.erase(remn_old.begin(),remn_old.end());
    type_old.erase(type_old.begin(),type_old.end());

    nform.erase(nform.begin(),nform.end());
    tform.erase(tform.begin(),tform.end());
    vform.erase(vform.begin(),vform.end());

    zams_new.erase(zams_new.begin(),zams_new.end());
    remn_new.erase(remn_new.begin(),remn_new.end());
    time_new.erase(time_new.begin(),time_new.end());
    kick_new.erase(kick_new.begin(),kick_new.end());
  }


  delete [] nmodel;


  return 0;
}
