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
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"

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
  info.mean=-1;
  info.std=-1;
  info.strength=-1;
  info.strength_m1=-1;
  info.strength_p1=-1;
  
  std::string filename="BgShapeGaussExp_SigShapeExpGaussExp_19fb_KinFit_"+prodXsec+"pb/HbbHbb_19invfb_mX"+mass+"_toy"+toy+"_MaxLikelihood.log";
  std::ifstream file(filename.c_str());
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
  
  int toys=200;
  
  std::cout<<"nPoints = "<<nPoints<<std::endl;
  std::cout<<"nInjections = "<<nInjections<<std::endl;
  
  TH1F *h_strength[11], *h_strength_p1[11], *h_strength_m1[11], *h_pull[11];
  TH1F *h_mean[11];
  TH1F *h_std[11];
  
  gROOT->SetStyle("Plain");
  TCanvas *c_StrengthBias=new TCanvas("c_StrengthBias", "c_StrengthBias", 700, 700);
  for (unsigned int i_inj=0; i_inj<nInjections; ++i_inj)
  {
    double injection=atof((injections[i_inj]).c_str());
    std::cout<<"Injection Strength = "<<injection<<std::endl;
    
    std::vector<double> masses;
    std::vector<double> strength_mean;
    std::vector<double> strength_m1, strength_p1;
    std::vector<double> strength_0;
    std::vector<double> pull_mean;
    std::vector<double> pull_rms;
    std::vector<double> mean_center;
    std::vector<double> mean_err;
    for (int i_mass=0; i_mass<nPoints; ++i_mass)
    {
      h_strength[i_mass]=new TH1F(("h_strength_"+mass[i_mass]+"_"+injections[i_inj]).c_str(), "Reconstructed Signal Strength", 200, 0, 3*injection);
      h_strength_p1[i_mass]=new TH1F(("h_strength_p1"+mass[i_mass]+"_"+injections[i_inj]).c_str(), "Reconstructed Signal Strength", 200, 0, 2*injection);
      h_strength_m1[i_mass]=new TH1F(("h_strength_m1"+mass[i_mass]+"_"+injections[i_inj]).c_str(), "Reconstructed Signal Strength", 200, 0, 2*injection);
      h_pull[i_mass]=new TH1F(("h_pull_"+mass[i_mass]+"_"+injections[i_inj]).c_str(), "Pull of the Reconstructed Signal Strength; (Strength-Injection)/Uncertainty", 40, -1/injection, 1/injection);
      h_mean[i_mass]=new TH1F(("h_mean_"+mass[i_mass]+"_"+injections[i_inj]).c_str(), "Spread of the Signal Mean", 1200, 0., 1200.);
      h_std[i_mass]=new TH1F(("h_std_"+mass[i_mass]+"_"+injections[i_inj]).c_str(), "Spread of the Signal Std Dev", 50, 0., 50.);
      for (unsigned int i_toy=0; i_toy<toys; ++i_toy)
      {
        PullInfo info=infoFromToy(mass[i_mass], injections[i_inj], std::string(itoa(i_toy)));
        h_strength[i_mass]->Fill(info.strength);
        h_strength_p1[i_mass]->Fill(info.strength_p1);
        h_strength_m1[i_mass]->Fill(info.strength_m1);
        double pull;
        if (info.strength-injection>0) pull=(info.strength-injection)/info.strength_m1;
        else pull=(info.strength-injection)/info.strength_p1;
        h_pull[i_mass]->Fill(pull);
        h_mean[i_mass]->Fill(info.mean);
        h_std[i_mass]->Fill(info.std);
      }
      masses.push_back(atof(mass[i_mass].c_str()));
      strength_mean.push_back(h_strength[i_mass]->GetMean());
      strength_p1.push_back(h_strength_p1[i_mass]->GetMean());
      strength_m1.push_back(h_strength_m1[i_mass]->GetMean());
      strength_0.push_back(0);
      pull_mean.push_back(h_pull[i_mass]->GetMean());
      pull_rms.push_back(h_pull[i_mass]->GetRMS());
      mean_center.push_back(h_mean[i_mass]->GetMean());
      mean_err.push_back(h_mean[i_mass]->GetMeanError());
    }
    
    c_StrengthBias->cd();
    TGraphAsymmErrors *g_strength=new TGraphAsymmErrors(nPoints, &masses[0], &strength_mean[0], &strength_0[0], &strength_0[0], &strength_m1[0], &strength_p1[0]);
    g_strength->SetTitle("Reconstructed Signal Strength for Different Injected Signal Strengths and Various m_{X}; m_{X} (GeV); Signal Strength");
    g_strength->SetMinimum(0.002);
    g_strength->SetMaximum(20);
    if (i_inj==0) g_strength->Draw("AL*");
    else g_strength->Draw("L*");
    // g_strength->SetFillColor(kGreen);
    // g_strength->Draw("3");
    // g_strength->Draw("L*");
    TLine *line=new TLine(400, injection, 1100, injection); line->SetLineColor(kRed); line->Draw();
    /*
    TCanvas *c_PullBias=new TCanvas("c_PullBias", "c_PullBias", 700, 700);
    TGraphErrors *g_pull=new TGraphErrors(nPoints, &masses[0], &pull_mean[0], &strength_0[0], &pull_rms[0]);
    // g_pull->SetMinimum(-1);
    // g_pull->SetMaximum(+1);
    g_pull->SetTitle("Signal Bias for Various m_{X}; m_{X} (GeV); Signal Bias/#sigma");
    g_pull->Draw("AL*");
    TLine *line1=new TLine(400, 0, 1100, 0); line1->SetLineColor(kRed); line1->Draw();
    c_PullBias->SaveAs(("c_PullBias_"+injections[i_inj]+"pb.png").c_str());
    
    TCanvas *c_MeanBias=new TCanvas("c_MeanBias", "c_MeanBias", 700, 700);
    TGraphErrors *g_mean=new TGraphErrors(nPoints, &masses[0], &mean_center[0], &strength_0[0], &mean_err[0]);
    g_mean->SetTitle(("Bias in Signal Mean with "+injections[i_inj]+" pb injection; m_{X} (GeV); Reconstructed Signal Mean (GeV)").c_str());
    g_mean->Draw("A*");
    g_mean->Fit("pol1");
    c_MeanBias->SaveAs(("c_MeanBias_"+injections[i_inj]+"pb.png").c_str());
    */
  }
  c_StrengthBias->SetLogy();
  c_StrengthBias->SaveAs("c_StrengthBias.png");
  
  

}
