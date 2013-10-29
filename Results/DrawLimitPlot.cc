#include <TROOT.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"

std::string itoa(int i) 
{
  char res[4];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

void DrawLimitPlot()
{
  const unsigned int nPoints=10;
  double xsec[nPoints], xsecNeg1[nPoints], xsecPos1[nPoints], xsecNeg2[nPoints], xsecPos2[nPoints];
  double obs[nPoints];
  double expNeg2[nPoints], expNeg1[nPoints], expPos1[nPoints], expPos2[nPoints];
  // double mass[nPoints]={300, 400, 500, 600, 700, 800};
  double mass[nPoints]={450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100};
  
  for (unsigned int i=0; i<nPoints; ++i)
  {
    std::string mass_string=itoa(mass[i]);
    std::string filename="HbbHbb_19invfb_mX"+mass_string+"_Asymptotic.log";
    std::ifstream file(filename.c_str(), ios::in);
    std::cout<<"Opened file "<<filename<<std::endl;
    std::string line;
    getline(file, line);
    getline(file, line);
    getline(file, line);
    
    getline(file, line); obs[i]=atof(line.substr(line.find("<")+1).c_str());
    getline(file, line); xsecNeg2[i]=atof(line.substr(line.find("<")+1).c_str());
    getline(file, line); xsecNeg1[i]=atof(line.substr(line.find("<")+1).c_str());
    getline(file, line); xsec[i]=atof(line.substr(line.find("<")+1).c_str());
    getline(file, line); xsecPos1[i]=atof(line.substr(line.find("<")+1).c_str());
    getline(file, line); xsecPos2[i]=atof(line.substr(line.find("<")+1).c_str());
    
    expNeg2[i]=xsec[i]-xsecNeg2[i];
    expNeg1[i]=xsec[i]-xsecNeg1[i];
    expPos1[i]=xsecPos1[i]-xsec[i];
    expPos2[i]=xsecPos2[i]-xsec[i];
    
    // std::cout<<"obs="<<obs[i]<<", exp="<<xsec[i]<<", -1sigma = "<<expNeg1[i]<<", +1sigma = "<<expPos1[i]<<std::endl;
    std::cout<<"exp="<<xsec[i]<<", -1sigma = "<<expNeg1[i]<<", +1sigma = "<<expPos1[i]<<std::endl;
    
  }
  
  gROOT->SetStyle("Plain");
  // gStyle->SetOptStat(000000000);
  
  TGraph *g_xsec=new TGraph(nPoints, mass, xsec);
  g_xsec->SetTitle("95% CL upper limits on a di-Higgs resonance decaying to 4 b-quarks; m_{X} (GeV); #sigma(X) #times Br(X#rightarrowH(b#bar{b}) H(b#bar{b})) (pb)");
  g_xsec->SetLineWidth(2);
  g_xsec->SetLineStyle(2);
  TGraphAsymmErrors *g_xsec_1sigma=new TGraphAsymmErrors(nPoints, mass, xsec, 0, 0, expNeg1, expPos1);
  g_xsec_1sigma->SetLineColor(kGreen);
  g_xsec_1sigma->SetFillColor(kGreen);
  TGraphAsymmErrors *g_xsec_2sigma=new TGraphAsymmErrors(nPoints, mass, xsec, 0, 0, expNeg2, expPos2);
  g_xsec_2sigma->SetLineColor(kYellow);
  g_xsec_2sigma->SetFillColor(kYellow);
  TGraph *g_obs=new TGraph(nPoints, mass, obs);
  g_obs->SetLineWidth(2);
  g_obs->SetLineStyle(1);
  TCanvas *c_xsec=new TCanvas("c_xsec", "c_xsec", 1000, 700);
  c_xsec->SetLogy();
  g_xsec->SetMaximum(5.0); g_xsec->SetMinimum(0.005);
  g_xsec->Draw("AL*");
  g_xsec_2sigma->Draw("3");
  g_xsec_1sigma->Draw("3");
  g_xsec->Draw("L*");
  // g_obs->Draw("L* SAME");
  TLegend *leg=new TLegend(0.45, 0.6, 0.9, 0.85);
  // TLegend *leg=new TLegend(0.45, 0.5, 0.9, 0.7);
  leg->SetFillStyle(0);
  leg->AddEntry((TObject*)0, "CMS Experiment. #sqrt{s} = 8 TeV, L = 18.6 fb^{-1}", "");
  leg->AddEntry((TObject*)0, "X#rightarrowH(b#bar{b}) H(b#bar{b})", "");
  leg->AddEntry((TObject*)0, "Jet Tagging: MMMM", "");
  leg->AddEntry(g_xsec, "Expected Limit", "LP");
  leg->AddEntry(g_xsec_1sigma, "Expected #pm 1 #sigma");
  leg->AddEntry(g_xsec_2sigma, "Expected #pm 2 #sigma");
  // leg->AddEntry(g_obs, "Observed Limit with fake data (background + 1 pb signal)", "LP");
  leg->Draw();
  c_xsec->SaveAs("UpperLimit.png");
  c_xsec->SaveAs("UpperLimit.eps");
  
  TFile *file=new TFile("UpperLimits_xsec.root", "RECREATE");
  g_obs->Write("g_obs");
  g_xsec->Write("g_xsec");
  g_xsec_1sigma->Write("g_xsec_1sigma");
  g_xsec_2sigma->Write("g_xsec_2sigma");
  file->Close();
  
}
    
