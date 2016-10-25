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

bool compareATLAS=false;
bool graviton=false;
bool radion=true;

std::string itoa(int i) 
{
  char res[4];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

void DrawLimitPlot_LMR_MMR_HMR()
{
  std::vector<double> mass;
  std::vector<double> xsec;
  std::vector<double> obs;
  std::vector<double> expNeg2, expNeg1, expPos1, expPos2;
  
  for (double imass=270; imass<450;)
  {
    mass.push_back(imass);
    
    std::string mass_string=itoa(imass);
    std::string filename="InterpolLMR/InterpolatedLogs_4b/InterpolLMR/HbbHbb_19invfb_mX"+mass_string+"_Asymptotic.log";
    std::ifstream file(filename.c_str(), ios::in);
    // std::cout<<"Opened file "<<filename<<std::endl;
    std::string line;
    getline(file, line);
    getline(file, line);
    getline(file, line);
    
    getline(file, line); double d_obs=atof(line.substr(line.find("<")+1).c_str())*1000; obs.push_back(d_obs);
    getline(file, line); double d_xsecNeg2=atof(line.substr(line.find("<")+1).c_str())*1000;
    getline(file, line); double d_xsecNeg1=atof(line.substr(line.find("<")+1).c_str())*1000;
    getline(file, line); double d_xsec=atof(line.substr(line.find("<")+1).c_str())*1000; xsec.push_back(d_xsec);
    getline(file, line); double d_xsecPos1=atof(line.substr(line.find("<")+1).c_str())*1000;
    getline(file, line); double d_xsecPos2=atof(line.substr(line.find("<")+1).c_str())*1000;
    
    expNeg2.push_back(d_xsec-d_xsecNeg2);
    expNeg1.push_back(d_xsec-d_xsecNeg1);
    expPos1.push_back(d_xsecPos1-d_xsec);
    expPos2.push_back(d_xsecPos2-d_xsec);
    
    // if (imass==350)
    // std::cout<<"mX = "<<mass_string<<", obs = "<<d_obs<<", exp limit ="<<d_xsec<<", -1 sigma = "<<d_xsec-d_xsecNeg1<<", +1 sigma = "<<d_xsecPos1-d_xsec<<std::endl;
    std::cout<<mass_string<<" & "<<d_obs<<" & "<<d_xsec<<" & "<<d_xsecNeg1<<" & "<<d_xsecPos1<<" & "<<d_xsecNeg2<<" & "<<d_xsecPos2<<" \\\\ "<<std::endl;
    
    
    if (imass<300) imass+=3;
    if (imass>=300 && imass<350) imass+=5;
    if (imass>=350) imass+=10;
    
  }
  
  for (double imass=450; imass<=730; imass+=10)
  {
    mass.push_back(imass);
    
    std::string mass_string=itoa(imass);
    std::string filename="InterpolMMR/InterpolDatacards_10GeV/f_datacard_"+mass_string+".log";
    std::ifstream file(filename.c_str(), ios::in);
    // std::cout<<"Opened file "<<filename<<std::endl;
    std::string line;
    getline(file, line);
    getline(file, line);
    getline(file, line);
    
    getline(file, line); double d_obs=atof(line.substr(line.find("<")+1).c_str())*1000; obs.push_back(d_obs);
    getline(file, line); double d_xsecNeg2=atof(line.substr(line.find("<")+1).c_str())*1000;
    getline(file, line); double d_xsecNeg1=atof(line.substr(line.find("<")+1).c_str())*1000;
    getline(file, line); double d_xsec=atof(line.substr(line.find("<")+1).c_str())*1000; xsec.push_back(d_xsec);
    getline(file, line); double d_xsecPos1=atof(line.substr(line.find("<")+1).c_str())*1000;
    getline(file, line); double d_xsecPos2=atof(line.substr(line.find("<")+1).c_str())*1000;
    
    expNeg2.push_back(d_xsec-d_xsecNeg2);
    expNeg1.push_back(d_xsec-d_xsecNeg1);
    expPos1.push_back(d_xsecPos1-d_xsec);
    expPos2.push_back(d_xsecPos2-d_xsec);
    
    // if (imass==700)
    // std::cout<<"mX = "<<mass_string<<", obs = "<<d_obs<<", exp limit ="<<d_xsec<<", -1 sigma = "<<d_xsec-d_xsecNeg1<<", +1 sigma = "<<d_xsecPos1-d_xsec<<std::endl;
    std::cout<<mass_string<<" & "<<d_obs<<" & "<<d_xsec<<" & "<<d_xsecNeg1<<" & "<<d_xsecPos1<<" & "<<d_xsecNeg2<<" & "<<d_xsecPos2<<" \\\\ "<<std::endl;
    
    
  }
  
  for (double imass=740; imass<=1100; imass+=10)
  {
    mass.push_back(imass);
    
    std::string mass_string=itoa(imass);
    std::string filename="InterpolHMR/InterpolDatacards_10GeV/f_datacard_"+mass_string+".log";
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
    
    // if (imass==900)
    // std::cout<<"mX = "<<mass_string<<", obs = "<<d_obs<<", exp limit ="<<d_xsec<<", -1 sigma = "<<d_xsec-d_xsecNeg1<<", +1 sigma = "<<d_xsecPos1-d_xsec<<std::endl;
    std::cout<<mass_string<<" & "<<d_obs<<" & "<<d_xsec<<" & "<<d_xsecNeg1<<" & "<<d_xsecPos1<<" & "<<d_xsecNeg2<<" & "<<d_xsecPos2<<" \\\\ "<<std::endl;
   
    
  }
  
  gROOT->SetStyle("Plain");
  // gStyle->SetOptStat(000000000);
  
  // ATLAS curve
  double masses_ATLAS[8]={500, 600, 700, 800, 900, 1000, 1100};
  double limit_ATLAS[8]={80, 25, 15, 12, 10, 9, 8};
  TGraph *g_ATLAS=new TGraph(7, masses_ATLAS, limit_ATLAS); g_ATLAS->SetLineWidth(2); g_ATLAS->SetLineColor(kRed);
  
  // Graviton curve
  double masses_graviton[19]={260,     300,     400,     450,     500,     550,     600,     650,    700,     750,     800,     850,     900,   950,     1000,    1050,    1100,    1150};
  double x_graviton[19]=     {10.7403, 185.658, 265.788, 199.824, 141.978, 97.9176, 68.5976, 47.067, 32.7501, 22.8192, 16.2399, 11.5632, 8.363, 6.01585, 4.44136, 3.24999, 2.42234, 1.8144};
  TGraph *g_graviton=new TGraph(18, masses_graviton, x_graviton); g_graviton->SetLineWidth(2); g_graviton->SetLineColor(kMagenta); g_graviton->SetFillColor(kWhite);
  
  // Radion curves
  // L_R = 3 TeV
  /*
  double masses_radion[11]={300, 350, 400, 450, 500, 550, 600, 700, 800, 900, 1000};
  double x_radion_3TeV[11]= {2272.0725413375076, 1446.7421962329545,  830.6271984524935,  581.3817599746898,  426.34584069784813,320.7551137323525, 246.3303275668544,152.53424031123356, 100.55254575140823, 69.89048380004162, 50.7888474680315};
  double x_radion_1TeV[11], x_radion_2TeV[11];
  for(int i=0; i<12; i++ )
  {
    x_radion_3TeV[i]=x_radion_3TeV[i]*0.25; //BR X->HH 
    x_radion_3TeV[i]=x_radion_3TeV[i]*0.57*0.57;
    
    x_radion_1TeV[i]=x_radion_3TeV[i]*9.;
    x_radion_2TeV[i]=x_radion_1TeV[i]*0.25;
  }
  TGraph *g_radion_3TeV=new TGraph(11, masses_radion, x_radion_3TeV); g_radion_3TeV->SetLineWidth(2); g_radion_3TeV->SetLineColor(kOrange+2); g_radion_3TeV->SetFillColor(kWhite);
  TGraph *g_radion_1TeV=new TGraph(11, masses_radion, x_radion_1TeV); g_radion_1TeV->SetLineWidth(2); g_radion_1TeV->SetLineColor(kRed); g_radion_1TeV->SetFillColor(kWhite);
  TGraph *g_radion_2TeV=new TGraph(11, masses_radion, x_radion_2TeV); g_radion_2TeV->SetLineWidth(2); g_radion_2TeV->SetLineColor(kBlue); g_radion_2TeV->SetFillColor(kWhite);
  */
  double masses_radion[]={250,
                          252,
                          254,
                          256,
                          258,
                          260,
                          280,
                          300,
                          320,
                          340,
                          360,
                          380,
                          400,
                          450,
                          500,
                          550,
                          600,
                          650,
                          700,
                          750,
                          800,
                          850,
                          900,
                          950,
                          1000,
                          1100};
  double x_radion_1TeV[]={0,
                          1447.244932,
                          1858.8143,
                          2102.391868,
                          2261.981552,
                          2370.797056,
                          2474.378497,
                          2154.631125,
                          1816.67727,
                          1524.826584,
                          1123.600715,
                          857.1077551,
                          690.037282,
                          446.399517,
                          312.8601625,
                          229.3429377,
                          173.5143952,
                          134.3418404,
                          106.0147116,
                          85.25396084,
                          69.62412993,
                          57.64151745,
                          48.37995955,
                          41.04551271,
                          35.19099514,
                          24.23219521};
  std::vector<double> v_masses_radion(masses_radion, masses_radion+sizeof(masses_radion)/sizeof(double));
  std::vector<double> v_x_radion_1TeV(x_radion_1TeV, x_radion_1TeV+sizeof(x_radion_1TeV)/sizeof(double));
  
  for (unsigned int i=0; i<v_masses_radion.size(); ++i)
  {
    std::cout<<"masses_radion = "<<masses_radion[i]<<", x_radion_1TeV = "<<x_radion_1TeV[i]<<std::endl;
  }
  
  // std::vector<double> v_x_radion_1TeV(x_radion_1TeV, x_radion_1TeV+sizeof(x_radion_1TeV)/sizeof(double));
  TGraph *g_radion_1TeV=new TGraph(v_masses_radion.size(), &v_masses_radion[0], &v_x_radion_1TeV[0]); g_radion_1TeV->SetLineWidth(2); g_radion_1TeV->SetLineColor(kRed); g_radion_1TeV->SetFillColor(kWhite);
  
  
  TStyle *tdrStyle=setTDRStyle();
  tdrStyle->cd();
  
  TGraph *g_xsec=new TGraph(mass.size(), &mass[0], &xsec[0]);
  g_xsec->SetTitle("; m_{X} [GeV]; 95\% CL limit #sigma(pp #rightarrow X #rightarrowHH #rightarrowb #bar{b}b#bar{b}) [fb]");
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
  TCanvas *c_xsec=new TCanvas("c_xsec", "c_xsec", 1000, 1000);
  c_xsec->SetLogy();
  // c_xsec->SetGridx(); c_xsec->SetGridy();
  g_xsec->SetMaximum(20000); g_xsec->SetMinimum(2);
  g_xsec->Draw("AL");
  g_xsec_2sigma->Draw("3");
  g_xsec_1sigma->Draw("3");
  g_xsec->Draw("L");
  g_obs->Draw("LP SAME");
  if (compareATLAS) g_ATLAS->Draw("L same");
  if (radion)
  {
    g_radion_1TeV->SetTitle("; m_{X} [GeV]; #sigma(pp#rightarrowX) #times Br(X#rightarrowHH#rightarrowb#bar{b}b#bar{b}) [fb]");
    g_radion_1TeV->SetMaximum(20000); g_radion_1TeV->SetMinimum(2);
    g_radion_1TeV->Draw("C same");
  }
  if (graviton) g_graviton->Draw("C same");
  TLegend *leg=new TLegend(0.40, 0.55, 0.9, 0.85);
  leg->SetFillStyle(1); leg->SetFillColor(kWhite); leg->SetLineColor(kWhite); leg->SetTextSize(0.04);
  /*leg->AddEntry(g_xsec, "Expected Upper Limit", "L");
  leg->AddEntry(g_xsec_1sigma, "Expected #pm 1 #sigma", "F");
  leg->AddEntry(g_xsec_2sigma, "Expected #pm 2 #sigma", "F");
  leg->AddEntry(g_obs, "Observed Upper Limit", "LP");*/
  if (compareATLAS) leg->AddEntry(g_ATLAS, "ATLAS expected limit", "L");
  if (graviton) leg->AddEntry(g_graviton, "RS1 KK-Graviton, k/m_{P}=0.1");
  if (radion) leg->AddEntry(g_radion_1TeV, "RS1 Radion, #Lambda_{R}=1 TeV");
  leg->Draw();
  TLatex * tPrel = new TLatex();
  tPrel->SetTextSize(0.03);
  // tPrel->DrawLatexNDC(0.1, 0.94, "CMS Preliminary; #sqrt{s} = 8 TeV, L = 17.93 fb^{-1}");
  // tPrel->DrawLatexNDC(0.1, 0.94, "Theoretical Cross Sections of the Radion");
  // tPrel->SetTextSize(0.03); tPrel->DrawLatexNDC(0.75, 0.65, "kL=35, no R/H mixing");
  c_xsec->Update();
  c_xsec->SaveAs("UpperLimit.png");
  c_xsec->SaveAs("UpperLimit.pdf");
  c_xsec->SaveAs("UpperLimit.root");
  
  TFile *file=new TFile("UpperLimits_Interpol.root", "RECREATE");
  g_obs->Write("g_obs");
  g_xsec->Write("g_xsec");
  g_xsec_1sigma->Write("g_xsec_1sigma");
  g_xsec_2sigma->Write("g_xsec_2sigma");
  if (graviton) g_graviton->Write("graviton");
  if (radion)
  {
    g_radion_1TeV->Write("g_radion_1TeV");
  }
  file->Close();
  
}
    
