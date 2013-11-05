#include <TROOT.h>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TAxis.h"

struct PullInfo
{
  double mean;
  double std;
  double strength;
  double strength_m1;
  double strength_p1;
};

std::string itoa(int i) 
{
  char res[4];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

PullInfo infoFromToy(std::string mass, std::string prodXsec, std::string toy)
{
  PullInfo info;
  info.mean=1;
  info.std=1;
  info.strength=1;
  info.strength_m1=1;
  info.strength_p1=1;
  
  // std::string filename="BgShapeGaussExp_SigShapeExpGaussExp_19fb_KinFit_"+prodXsec+"pb/HbbHbb_19invfb_mX"+mass+"_toy"+toy+"_MaxLikelihood.log";
  std::string filename="HbbHbb_19invfb_mX"+mass+"_toy"+toy+"_MaxLikelihood.log";
  std::ifstream file(filename.c_str());
  // std::cout<<"Opened file "<<filename<<std::endl;
  std::string line;
  
  while (line!="    Floating Parameter  InitialValue    FinalValue (+HiError,-LoError)    GblCorr." && !file.eof()) getline(file, line);
  if (file.eof()) {file.close(); return info;}
  while (line.find("signal_p0")==std::string::npos && !file.eof()) getline(file, line);
  if (file.eof()) {file.close(); return info;}
  std::string mean_string, std_string;
  double mean_old, std_old;
  std::stringstream lineStream(line);
  lineStream>>mean_string>>mean_old>>info.mean;
  getline(file, line);
  lineStream.str(line);
  lineStream>>std_string>>std_old>>info.std;
  
  while (line!=" --- MaxLikelihoodFit ---" && !file.eof()) getline(file, line);
  if (file.eof()) {file.close(); return info;}
  getline(file, line);
  int pos1=line.find("r:");
  int pos2=line.find("-", pos1);
  int pos3=line.find("+", pos2);
  int pos4=line.find(" ", pos3);
  info.strength=atof(line.substr(pos1+3, (pos2-pos1-5)).c_str());
  info.strength_m1=atof(line.substr(pos2+1, (pos3-pos2-2)).c_str());
  info.strength_p1=atof(line.substr(pos3+1, (pos4-pos3-1)).c_str());
  
  file.close();
  return info;
}

void PullMassMean()
{
  std::string mass[]={"400", "450", "500", "550", "600", "650", "700", "800", "900", "1000", "1100"};
  double cols[]={kRed-2, kRed, kRed+2, kGreen, kBlue, kOrange, kBlack, kViolet, kCyan, kMagenta, kBlue+2};
  int nPoints=sizeof(mass)/sizeof(mass[0]);
  
  std::string injections[]={"0.01", "0.05", "0.1", "0.5", "1", "5", "10"};
  int nInjections=sizeof(injections)/sizeof(injections[0]);
  
  int toys=100;
  
  std::cout<<"nPoints = "<<nPoints<<std::endl;
  std::cout<<"nInjections = "<<nInjections<<std::endl;
  
  TH1F *h_strength[11];
  TH1F *h_mean[11];
  TH1F *h_std[11];
  
  // for (unsigned int i_inj=0; i_inj<nInjections; ++i_inj)
  // {
    for (int i_mass=0; i_mass<nPoints; ++i_mass)
    {
      h_strength[i_mass]=new TH1F(("h_strength_"+mass[i_mass]).c_str(), "Spread of the Signal Strength", 40, 0.6, 1.4);
      // double mean_mean=4.+(20.-4.)*(atof((mass[i_mass].c_str()))-400.)/(1100.-400.);
      // std::cout<<"mass = "<<mass[i_mass]<<" mean_mean = "<<mean_mean<<std::endl;
      // h_mean[i_mass]=new TH1F(("h_mean_"+mass[i_mass]).c_str(), "Signal Mean", 200, mean_mean-5., mean_mean+5.);
      h_mean[i_mass]=new TH1F(("h_mean_"+mass[i_mass]).c_str(), "Spread of the Signal Mean", 100, 5., 20.);
      h_std[i_mass]=new TH1F(("h_std_"+mass[i_mass]).c_str(), "Spread of the Signal Std Dev", 400, 0., 40.);
      for (unsigned int i_toy=0; i_toy<toys; ++i_toy)
      {
        PullInfo info=infoFromToy(mass[i_mass], "1", std::string(itoa(i_toy)));
        h_strength[i_mass]->Fill(info.strength);
        h_mean[i_mass]->Fill(info.mean-atof(mass[i_mass].c_str()));
        h_std[i_mass]->Fill(info.std);
      }
    }
  // }
  
  gROOT->SetStyle("Plain");
  
  TCanvas *c_strength;
  TCanvas *c_mean;
  TCanvas *c_std;
  for (unsigned int i_mass=0; i_mass<nPoints; ++i_mass)
  {
    c_strength=new TCanvas("c_strength", "Spread of the Strength", 700, 700);
    h_strength[i_mass]->SetLineColor(cols[i_mass]);
    h_strength[i_mass]->Draw("e");
    h_strength[i_mass]->Fit("gaus");
    c_strength->SaveAs(("c_strength_"+mass[i_mass]+".png").c_str());
    
    c_mean=new TCanvas("c_mean", "Spread of the Mean", 700, 700);
    h_mean[i_mass]->SetLineColor(cols[i_mass]);
    h_mean[i_mass]->Draw("e");
    h_mean[i_mass]->Fit("gaus");
    c_mean->SaveAs(("c_mean_"+mass[i_mass]+".png").c_str());
    
    c_std=new TCanvas("c_std", "Spread of the Std Dev", 700, 700);
    h_std[i_mass]->SetLineColor(cols[i_mass]);
    h_std[i_mass]->Draw("C");
    // h_std[i_mass]->Fit("gaus");
    c_std->SaveAs(("h_std_"+mass[i_mass]+".png").c_str());
  }
  
  // delete h_strength;
  // delete h_mean;
  // delete h_std;

}
