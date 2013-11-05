#include <TROOT.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TF1.h"
#include "TAxis.h"

std::string itoa(int i) 
{
  char res[4];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

std::string dirName(int j)
{
  std::string dir="Step2V8_BgShapeSyst_SigShapeGaussianSyst_19fb_KinFit";
  if (j==0) dir=dir+"_1pb";
  else if (j==1) dir=dir+"_p5pb";
  else if (j==2) dir=dir+"_p1pb";
  else if (j==3) dir=dir+"_p05pb";
  else if (j==4) dir=dir+"_p01pb";
  return dir;
}

void ClosureTest_SignalStrength()
{
  
  std::string mass[]={"400", "450", "500", "550", "600", "650", "700", "800", "900", "1000", "1100"};
  int nMasses=sizeof(mass)/sizeof(mass[0]);
  
  std::string injections[]={"0.01", "0.05", "0.1", "0.5", "1", "5", "10"};
  int nInjections=sizeof(injections)/sizeof(injections[0]);
  
  gROOT->SetStyle("Plain");
  
  TGraphAsymmErrors *g_strength[nMasses];
  TCanvas *c_strength[nMasses];
  
  for (unsigned int i=0; i<nMasses; ++i)
  {
    std::string mass_string=mass[i];
    
    double xsec[nInjections];
  
    double strengths_Mass[nInjections];
    double strengths_m1_Mass[nInjections];
    double strengths_p1_Mass[nInjections];
      
    for (unsigned int j=0; j<nInjections; ++j)
    {
      xsec[j]=atof((injections[j]).c_str());
      
      std::string filename="BgShapeGaussExp_SigShapeExpGaussExp_19fb_KinFit_"+injections[j]+"pb/HbbHbb_19invfb_mX"+mass_string+"_MaxLikelihood.log";
      std::ifstream file(filename.c_str(), ios::in);
      std::cout<<"Opened file "<<filename<<std::endl;
      std::string line;
      while (line!=" --- MaxLikelihoodFit ---")
      {
        getline(file, line);
      }
      getline(file, line);
      // getline(file, line);
      int pos1=line.find("r:");
      int pos2=line.find("-", pos1);
      int pos3=line.find("+", pos2);
      int pos4=line.find(" ", pos3);
      
      strengths_Mass[j]=atof(line.substr(pos1+3, (pos2-pos1-5)).c_str());
      strengths_m1_Mass[j]=atof(line.substr(pos2+1, (pos3-pos2-2)).c_str());
      strengths_p1_Mass[j]=atof(line.substr(pos3+1, (pos4-pos3-1)).c_str());
      
      //  std::cout<<"Mass = "<<mass[i]<<", xsec = "<<xsec[j]<<std::endl;
      //  std::cout<<strengths_Mass[j]<<std::endl;
      //  std::cout<<strengths_m1_Mass[j]<<std::endl;
      //  std::cout<<strengths_p1_Mass[j]<<std::endl;
    }
    
    g_strength[i]=new TGraphAsymmErrors(nInjections, xsec, strengths_Mass, 0, 0, strengths_m1_Mass, strengths_p1_Mass);
    // g_strength[i]->SetLineColor(cols[i]);
    g_strength[i]->SetTitle(("Signal Strength of Maximum Likelihood for mass "+mass_string+"; true x-sec (pb); returned x-sec (pb)").c_str());
    c_strength[i]=new TCanvas("c_strength", "c_strength", 700, 700);
    g_strength[i]->Draw("A*");
    c_strength[i]->SetLogy(); c_strength[i]->SetLogx();
    g_strength[i]->GetXaxis()->SetLimits(0.005, 50.);
    g_strength[i]->GetYaxis()->SetRangeUser(0.005, 50.);
    TF1 *f_ratio=new TF1("f_ratio", "pol1");
    g_strength[i]->Fit(f_ratio);
    c_strength[i]->SaveAs(("c_strength_"+mass_string+".png").c_str());
    
  }
  
  
}
