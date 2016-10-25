#include <TROOT.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TLatex.h"

#include "/Users/souvik/HbbHbb/Analysis/TDRStyle.h"

bool compareATLAS=true;

std::string itoa(int i) 
{
  char res[4];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

void DrawLimitPlot()
{
  // double xsec[nPoints], xsecNeg1[nPoints], xsecPos1[nPoints], xsecNeg2[nPoints], xsecPos2[nPoints];
  // double obs[nPoints];
  // double expNeg2[nPoints], expNeg1[nPoints], expPos1[nPoints], expPos2[nPoints];
  // double mass[nPoints]={450, 500, 550, 600, 650, 700, 800};
  // double mass[nPoints]={400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100};
  
  std::vector<double> mass;
  std::vector<double> xsec;
  std::vector<double> obs;
  std::vector<double> expNeg2, expNeg1, expPos1, expPos2;
  
  for (double imass=700; imass<=1100; imass+=10)
  {
    mass.push_back(imass);
    
    std::string mass_string=itoa(imass);
    std::string filename="f_datacard_"+mass_string+".log";
    std::ifstream file(filename.c_str(), ios::in);
    // std::cout<<"Opened file "<<filename<<std::endl;
    std::string line;
    getline(file, line);
    getline(file, line);
    getline(file, line);
    
    getline(file, line); double d_obs=atof(line.substr(line.find("<")+1).c_str())*100; obs.push_back(d_obs);
    getline(file, line); double d_xsecNeg2=atof(line.substr(line.find("<")+1).c_str())*100;
    getline(file, line); double d_xsecNeg1=atof(line.substr(line.find("<")+1).c_str())*100;
    getline(file, line); double d_xsec=atof(line.substr(line.find("<")+1).c_str())*100; xsec.push_back(d_xsec);
    getline(file, line); double d_xsecPos1=atof(line.substr(line.find("<")+1).c_str())*100;
    getline(file, line); double d_xsecPos2=atof(line.substr(line.find("<")+1).c_str())*100;
    
    expNeg2.push_back(d_xsec-d_xsecNeg2);
    expNeg1.push_back(d_xsec-d_xsecNeg1);
    expPos1.push_back(d_xsecPos1-d_xsec);
    expPos2.push_back(d_xsecPos2-d_xsec);
    
    if (imass==700 ||
        imass==800 ||
        imass==900 ||
        imass==1000 ||
        imass==1100)
    std::cout<<"mX = "<<mass_string<<", obs = "<<d_obs<<", exp limit ="<<d_xsec<<", -1 sigma = "<<d_xsec-d_xsecNeg1<<", +1 sigma = "<<d_xsecPos1-d_xsec<<std::endl;
    
  }
  
  gROOT->SetStyle("Plain");
  // gStyle->SetOptStat(000000000);
  
  // ATLAS curve
  double masses_ATLAS[8]={500, 600, 700, 800, 900, 1000, 1100};
  double limit_ATLAS[8]={80, 25, 15, 12, 10, 9, 8};
  TGraph *g_ATLAS=new TGraph(7, masses_ATLAS, limit_ATLAS); g_ATLAS->SetLineWidth(2); g_ATLAS->SetLineColor(kRed);
  
  TStyle *tdrStyle=setTDRStyle();
  tdrStyle->cd();
  
  TGraph *g_xsec=new TGraph(mass.size(), &mass[0], &xsec[0]);
  g_xsec->SetTitle("; m_{X} (GeV); #sigma(pp#rightarrowX) #times Br(X#rightarrowH(b#bar{b}) H(b#bar{b})) (fb)");
  g_xsec->SetLineWidth(2);
  g_xsec->SetLineStyle(2);
  TGraphAsymmErrors *g_xsec_1sigma=new TGraphAsymmErrors(mass.size(), &mass[0], &xsec[0], 0, 0, &expNeg1[0], &expPos1[0]);
  g_xsec_1sigma->SetLineColor(kGreen);
  g_xsec_1sigma->SetFillColor(kGreen);
  TGraphAsymmErrors *g_xsec_2sigma=new TGraphAsymmErrors(mass.size(), &mass[0], &xsec[0], 0, 0, &expNeg2[0], &expPos2[0]);
  g_xsec_2sigma->SetLineColor(kYellow);
  g_xsec_2sigma->SetFillColor(kYellow);
  TGraph *g_obs=new TGraph(mass.size(), &mass[0], &obs[0]);
  g_obs->SetLineWidth(2);
  g_obs->SetLineStyle(1);
  g_obs->SetMarkerStyle(20);
  TCanvas *c_xsec=new TCanvas("c_xsec", "c_xsec", 1000, 700);
  c_xsec->SetLogy();
  c_xsec->SetGridx(); c_xsec->SetGridy();
  // g_xsec->SetMaximum(500); g_xsec->SetMinimum(1);
  g_xsec->SetMaximum(1000); g_xsec->SetMinimum(1);
  g_xsec->Draw("AL");
  g_xsec_2sigma->Draw("3");
  g_xsec_1sigma->Draw("3");
  g_xsec->Draw("L");
  g_obs->Draw("LP SAME");
  if (compareATLAS) g_ATLAS->Draw("L same");
  TLegend *leg=new TLegend(0.45, 0.6, 0.9, 0.85);
  leg->SetFillStyle(1); leg->SetFillColor(kWhite);
  leg->AddEntry(g_xsec, "Expected Upper Limit", "L");
  leg->AddEntry(g_xsec_1sigma, "Expected #pm 1 #sigma", "F");
  leg->AddEntry(g_xsec_2sigma, "Expected #pm 2 #sigma", "F");
  leg->AddEntry(g_obs, "Observed Upper Limit", "LP");
  if (compareATLAS) leg->AddEntry(g_ATLAS, "ATLAS expected limit", "L");
  leg->Draw();
  TLatex * tPrel = new TLatex();
  tPrel->SetTextSize(0.05);
  tPrel->DrawLatexNDC(0.1, 0.94, "CMS Preliminary; #sqrt{s} = 8 TeV, L = 17.928 fb^{-1}");
  c_xsec->Update();
  c_xsec->SaveAs("UpperLimit.png");
  c_xsec->SaveAs("UpperLimit.pdf");
  
  TFile *file=new TFile("UpperLimits_Interpol.root", "RECREATE");
  g_obs->Write("g_obs");
  g_xsec->Write("g_xsec");
  g_xsec_1sigma->Write("g_xsec_1sigma");
  g_xsec_2sigma->Write("g_xsec_2sigma");
  file->Close();
  
}
    
