#include <TH1F.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>

#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooRealVar.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooDataHist.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooWorkspace.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooPlot.h"

bool fillRooFit=false;
bool drawPlots=true;

double SR_lo=320.;
double SR_hi=1200.;

double lnN(double b, double a, double c)
{
  double err=0;
  if (b>0) err=1.+fabs(a-c)/(2.*b);
  return err;
}

int CreateRooFit_ttbar()
{

  // Set constants

  double totalLuminosity=17928; // /pb
  
  double xsec_ttbar_fulllept=24.56;
  double xsec_ttbar_semilept=103.12;
  double xsec_ttbar_hadronic=106.32;
  
  int rebin=4;
  
  // Open files
  
  TFile *ttbar_fulllept=new TFile("MMMM_nominal/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept=new TFile("MMMM_nominal/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic=new TFile("MMMM_nominal/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_JECp1=new TFile("KinSel_JECp1/MMMM_nominal/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_JECp1=new TFile("KinSel_JECp1/MMMM_nominal/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_JECp1=new TFile("KinSel_JECp1/MMMM_nominal/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_JECm1=new TFile("KinSel_JECm1/MMMM_nominal/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_JECm1=new TFile("KinSel_JECm1/MMMM_nominal/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_JECm1=new TFile("KinSel_JECm1/MMMM_nominal/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_JERp1=new TFile("KinSel_JERp1/MMMM_nominal/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_JERp1=new TFile("KinSel_JERp1/MMMM_nominal/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_JERp1=new TFile("KinSel_JERp1/MMMM_nominal/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_JERm1=new TFile("KinSel_JERm1/MMMM_nominal/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_JERm1=new TFile("KinSel_JERm1/MMMM_nominal/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_JERm1=new TFile("KinSel_JERm1/MMMM_nominal/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_upBC=new TFile("MMMM_upBC/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_upBC=new TFile("MMMM_upBC/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_upBC=new TFile("MMMM_upBC/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_downBC=new TFile("MMMM_downBC/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_downBC=new TFile("MMMM_downBC/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_downBC=new TFile("MMMM_downBC/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_upL=new TFile("MMMM_upL/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_upL=new TFile("MMMM_upL/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_upL=new TFile("MMMM_upL/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_downL=new TFile("MMMM_downL/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_downL=new TFile("MMMM_downL/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_downL=new TFile("MMMM_downL/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_upTrigCSV=new TFile("MMMM_upTrigCSV/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_upTrigCSV=new TFile("MMMM_upTrigCSV/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_upTrigCSV=new TFile("MMMM_upTrigCSV/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_downTrigCSV=new TFile("MMMM_downTrigCSV/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_downTrigCSV=new TFile("MMMM_downTrigCSV/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_downTrigCSV=new TFile("MMMM_downTrigCSV/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_upTrigpT=new TFile("MMMM_upTrigpT/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_upTrigpT=new TFile("MMMM_upTrigpT/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_upTrigpT=new TFile("MMMM_upTrigpT/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  TFile *ttbar_fulllept_downTrigpT=new TFile("MMMM_downTrigpT/a_KinFit/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept_downTrigpT=new TFile("MMMM_downTrigpT/a_KinFit/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic_downTrigpT=new TFile("MMMM_downTrigpT/a_KinFit/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  // Extract histograms
  
  TH1F *h_mX_SR_ttbar_fulllept=(TH1F*)ttbar_fulllept->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept=(TH1F*)ttbar_semilept->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic=(TH1F*)ttbar_hadronic->Get("h_mX_SR");
  
  TH1F *h_mX_SR_ttbar_fulllept_JECp1=(TH1F*)ttbar_fulllept_JECp1->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_JECp1=(TH1F*)ttbar_semilept_JECp1->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_JECp1=(TH1F*)ttbar_hadronic_JECp1->Get("h_mX_SR");
  
  TH1F *h_mX_SR_ttbar_fulllept_JECm1=(TH1F*)ttbar_fulllept_JECm1->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_JECm1=(TH1F*)ttbar_semilept_JECm1->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_JECm1=(TH1F*)ttbar_hadronic_JECm1->Get("h_mX_SR");
  
  TH1F *h_mX_SR_ttbar_fulllept_JERp1=(TH1F*)ttbar_fulllept_JERp1->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_JERp1=(TH1F*)ttbar_semilept_JERp1->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_JERp1=(TH1F*)ttbar_hadronic_JERp1->Get("h_mX_SR");

  TH1F *h_mX_SR_ttbar_fulllept_JERm1=(TH1F*)ttbar_fulllept_JERm1->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_JERm1=(TH1F*)ttbar_semilept_JERm1->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_JERm1=(TH1F*)ttbar_hadronic_JERm1->Get("h_mX_SR");
  
  TH1F *h_mX_SR_ttbar_fulllept_upBC=(TH1F*)ttbar_fulllept_upBC->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_upBC=(TH1F*)ttbar_semilept_upBC->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_upBC=(TH1F*)ttbar_hadronic_upBC->Get("h_mX_SR");

  TH1F *h_mX_SR_ttbar_fulllept_downBC=(TH1F*)ttbar_fulllept_downBC->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_downBC=(TH1F*)ttbar_semilept_downBC->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_downBC=(TH1F*)ttbar_hadronic_downBC->Get("h_mX_SR");
  
  TH1F *h_mX_SR_ttbar_fulllept_upL=(TH1F*)ttbar_fulllept_upL->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_upL=(TH1F*)ttbar_semilept_upL->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_upL=(TH1F*)ttbar_hadronic_upL->Get("h_mX_SR");

  TH1F *h_mX_SR_ttbar_fulllept_downL=(TH1F*)ttbar_fulllept_downL->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_downL=(TH1F*)ttbar_semilept_downL->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_downL=(TH1F*)ttbar_hadronic_downL->Get("h_mX_SR");
  
  TH1F *h_mX_SR_ttbar_fulllept_upTrigCSV=(TH1F*)ttbar_fulllept_upTrigCSV->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_upTrigCSV=(TH1F*)ttbar_semilept_upTrigCSV->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_upTrigCSV=(TH1F*)ttbar_hadronic_upTrigCSV->Get("h_mX_SR");

  TH1F *h_mX_SR_ttbar_fulllept_downTrigCSV=(TH1F*)ttbar_fulllept_downTrigCSV->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_downTrigCSV=(TH1F*)ttbar_semilept_downTrigCSV->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_downTrigCSV=(TH1F*)ttbar_hadronic_downTrigCSV->Get("h_mX_SR");
  
  TH1F *h_mX_SR_ttbar_fulllept_upTrigpT=(TH1F*)ttbar_fulllept_upTrigpT->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_upTrigpT=(TH1F*)ttbar_semilept_upTrigpT->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_upTrigpT=(TH1F*)ttbar_hadronic_upTrigpT->Get("h_mX_SR");

  TH1F *h_mX_SR_ttbar_fulllept_downTrigpT=(TH1F*)ttbar_fulllept_downTrigpT->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept_downTrigpT=(TH1F*)ttbar_semilept_downTrigpT->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic_downTrigpT=(TH1F*)ttbar_hadronic_downTrigpT->Get("h_mX_SR");
  
  // Set scales
  
  double init_ttbar_fulllept=((TH1F*)ttbar_fulllept->Get("CountWithPU"))->GetBinContent(1);
  double init_ttbar_semilept=((TH1F*)ttbar_semilept->Get("CountWithPU"))->GetBinContent(1);
  double init_ttbar_hadronic=((TH1F*)ttbar_hadronic->Get("CountWithPU"))->GetBinContent(1);
  
  std::cout<<"init_ttbar_fulllept = "<<init_ttbar_fulllept<<std::endl;
  std::cout<<"init_ttbar_semilept = "<<init_ttbar_semilept<<std::endl;
  std::cout<<"init_ttbar_hadronic = "<<init_ttbar_hadronic<<std::endl;
  
  double scale_ttbar_fulllept=totalLuminosity*xsec_ttbar_fulllept/init_ttbar_fulllept;
  double scale_ttbar_semilept=totalLuminosity*xsec_ttbar_semilept/init_ttbar_semilept;
  double scale_ttbar_hadronic=totalLuminosity*xsec_ttbar_hadronic/init_ttbar_hadronic;
  
  h_mX_SR_ttbar_fulllept->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_JECp1->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_JECp1->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_JECp1->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_JECm1->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_JECm1->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_JECm1->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_JERp1->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_JERp1->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_JERp1->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_JERm1->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_JERm1->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_JERm1->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_upBC->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_upBC->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_upBC->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_downBC->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_downBC->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_downBC->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_upL->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_upL->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_upL->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_downL->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_downL->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_downL->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_upTrigCSV->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_upTrigCSV->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_upTrigCSV->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_downTrigCSV->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_downTrigCSV->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_downTrigCSV->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_upTrigpT->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_upTrigpT->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_upTrigpT->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept_downTrigpT->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept_downTrigpT->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic_downTrigpT->Scale(scale_ttbar_hadronic);
  
  // Add up the three ttbar components
  
  TH1F *h_mX_SR_ttbar=(TH1F*)h_mX_SR_ttbar_fulllept->Clone("h_mX_SR_ttbar");
  h_mX_SR_ttbar->Add(h_mX_SR_ttbar_semilept);
  h_mX_SR_ttbar->Add(h_mX_SR_ttbar_hadronic);
  
  TH1F *h_mX_SR_ttbar_JECp1=(TH1F*)h_mX_SR_ttbar_fulllept_JECp1->Clone("h_mX_SR_ttbar_JECp1");
  h_mX_SR_ttbar_JECp1->Add(h_mX_SR_ttbar_semilept_JECp1);
  h_mX_SR_ttbar_JECp1->Add(h_mX_SR_ttbar_hadronic_JECp1);
  
  TH1F *h_mX_SR_ttbar_JECm1=(TH1F*)h_mX_SR_ttbar_fulllept_JECm1->Clone("h_mX_SR_ttbar_JECm1");
  h_mX_SR_ttbar_JECm1->Add(h_mX_SR_ttbar_semilept_JECm1);
  h_mX_SR_ttbar_JECm1->Add(h_mX_SR_ttbar_hadronic_JECm1);
  
  TH1F *h_mX_SR_ttbar_JERp1=(TH1F*)h_mX_SR_ttbar_fulllept_JERp1->Clone("h_mX_SR_ttbar_JERp1");
  h_mX_SR_ttbar_JERp1->Add(h_mX_SR_ttbar_semilept_JERp1);
  h_mX_SR_ttbar_JERp1->Add(h_mX_SR_ttbar_hadronic_JERp1);
  
  TH1F *h_mX_SR_ttbar_JERm1=(TH1F*)h_mX_SR_ttbar_fulllept_JERm1->Clone("h_mX_SR_ttbar_JERm1");
  h_mX_SR_ttbar_JERm1->Add(h_mX_SR_ttbar_semilept_JERm1);
  h_mX_SR_ttbar_JERm1->Add(h_mX_SR_ttbar_hadronic_JERm1);
  
  TH1F *h_mX_SR_ttbar_upBC=(TH1F*)h_mX_SR_ttbar_fulllept_upBC->Clone("h_mX_SR_ttbar_upBC");
  h_mX_SR_ttbar_upBC->Add(h_mX_SR_ttbar_semilept_upBC);
  h_mX_SR_ttbar_upBC->Add(h_mX_SR_ttbar_hadronic_upBC);
  
  TH1F *h_mX_SR_ttbar_downBC=(TH1F*)h_mX_SR_ttbar_fulllept_downBC->Clone("h_mX_SR_ttbar_downBC");
  h_mX_SR_ttbar_downBC->Add(h_mX_SR_ttbar_semilept_downBC);
  h_mX_SR_ttbar_downBC->Add(h_mX_SR_ttbar_hadronic_downBC);
  
  TH1F *h_mX_SR_ttbar_upL=(TH1F*)h_mX_SR_ttbar_fulllept_upL->Clone("h_mX_SR_ttbar_upL");
  h_mX_SR_ttbar_upL->Add(h_mX_SR_ttbar_semilept_upL);
  h_mX_SR_ttbar_upL->Add(h_mX_SR_ttbar_hadronic_upL);
  
  TH1F *h_mX_SR_ttbar_downL=(TH1F*)h_mX_SR_ttbar_fulllept_downL->Clone("h_mX_SR_ttbar_downL");
  h_mX_SR_ttbar_downL->Add(h_mX_SR_ttbar_semilept_downL);
  h_mX_SR_ttbar_downL->Add(h_mX_SR_ttbar_hadronic_downL);
  
  TH1F *h_mX_SR_ttbar_upTrigCSV=(TH1F*)h_mX_SR_ttbar_fulllept_upTrigCSV->Clone("h_mX_SR_ttbar_upTrigCSV");
  h_mX_SR_ttbar_upTrigCSV->Add(h_mX_SR_ttbar_semilept_upTrigCSV);
  h_mX_SR_ttbar_upTrigCSV->Add(h_mX_SR_ttbar_hadronic_upTrigCSV);
  
  TH1F *h_mX_SR_ttbar_downTrigCSV=(TH1F*)h_mX_SR_ttbar_fulllept_downTrigCSV->Clone("h_mX_SR_ttbar_downTrigCSV");
  h_mX_SR_ttbar_downTrigCSV->Add(h_mX_SR_ttbar_semilept_downTrigCSV);
  h_mX_SR_ttbar_downTrigCSV->Add(h_mX_SR_ttbar_hadronic_downTrigCSV);
  
  TH1F *h_mX_SR_ttbar_upTrigpT=(TH1F*)h_mX_SR_ttbar_fulllept_upTrigpT->Clone("h_mX_SR_ttbar_upTrigpT");
  h_mX_SR_ttbar_upTrigpT->Add(h_mX_SR_ttbar_semilept_upTrigpT);
  h_mX_SR_ttbar_upTrigpT->Add(h_mX_SR_ttbar_hadronic_upTrigpT);
  
  TH1F *h_mX_SR_ttbar_downTrigpT=(TH1F*)h_mX_SR_ttbar_fulllept_downTrigpT->Clone("h_mX_SR_ttbar_downTrigpT");
  h_mX_SR_ttbar_downTrigpT->Add(h_mX_SR_ttbar_semilept_downTrigpT);
  h_mX_SR_ttbar_downTrigpT->Add(h_mX_SR_ttbar_hadronic_downTrigpT);
  
  int bin1=h_mX_SR_ttbar->FindBin(SR_lo);
  int bin2=h_mX_SR_ttbar->FindBin(SR_hi)-1;
  
  // Print out the lnN systematics of normalization
  double ttbar=h_mX_SR_ttbar->Integral(bin1, bin2);
  double JECp1=h_mX_SR_ttbar_JECp1->Integral(bin1, bin2);
  double JECm1=h_mX_SR_ttbar_JECm1->Integral(bin1, bin2);
  double JERp1=h_mX_SR_ttbar_JERp1->Integral(bin1, bin2);
  double JERm1=h_mX_SR_ttbar_JERm1->Integral(bin1, bin2);
  double upBC=h_mX_SR_ttbar_upBC->Integral(bin1, bin2);
  double downBC=h_mX_SR_ttbar_downBC->Integral(bin1, bin2);
  double upL=h_mX_SR_ttbar_upL->Integral(bin1, bin2);
  double downL=h_mX_SR_ttbar_downL->Integral(bin1, bin2);
  double upTrigCSV=h_mX_SR_ttbar_upTrigCSV->Integral(bin1, bin2);
  double downTrigCSV=h_mX_SR_ttbar_downTrigCSV->Integral(bin1, bin2);
  double upTrigpT=h_mX_SR_ttbar_upTrigpT->Integral(bin1, bin2);
  double downTrigpT=h_mX_SR_ttbar_downTrigpT->Integral(bin1, bin2);
  
  std::cout<<"h_mX_SR_ttbar->Integral(bin1, bin2) = "<<h_mX_SR_ttbar->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_JECp1->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_JECp1->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_JECm1->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_JECm1->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_JERp1->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_JERp1->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_JERm1->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_JERm1->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_upBC->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_upBC->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_downBC->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_downBC->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_upL->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_upL->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_downL->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_downL->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_upTrigCSV->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_upTrigCSV->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_downTrigCSV->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_downTrigCSV->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_upTrigpT->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_upTrigpT->Integral(bin1, bin2)<<std::endl;
  std::cout<<"h_mX_SR_ttbar_downTrigpT->Integral(bin1, bin2) = "<<h_mX_SR_ttbar_downTrigpT->Integral(bin1, bin2)<<std::endl;
  
  std::cout<<"JEC lnN  - - "<<lnN(ttbar, JECp1, JECm1)<<std::endl;
  std::cout<<"JER lnN  - - "<<lnN(ttbar, JERp1, JERm1)<<std::endl;
  std::cout<<"bTagSFbc - - "<<lnN(ttbar, upBC, downBC)<<std::endl;
  std::cout<<"bTagSFl  - - "<<lnN(ttbar, upL, downL)<<std::endl;
  std::cout<<"trigSFCSV - - "<<lnN(ttbar, upTrigCSV, downTrigCSV)<<std::endl;
  std::cout<<"trigSFpT - - "<<lnN(ttbar, upTrigpT, downTrigpT)<<std::endl;
  
  // Put the histograms into a RooFit workspace and save it
  
  if (fillRooFit)
  {/*
    RooRealVar x("x", "m_{X} (GeV)", SR_lo, SR_hi);
    x.setBins(88);
  
    RooDataHist ttbar("ttbar", "ttbar", RooArgList(x), h_mX_SR_ttbar);
  
    RooDataHist ttbar_JECUp("ttbar_JECUp", "ttbar_JECUp", RooArgList(x), h_mX_SR_ttbar_JECp1);
    RooDataHist ttbar_JECDown("ttbar_JECDown", "ttbar_JECDown", RooArgList(x), h_mX_SR_ttbar_JECm1);
  
    RooDataHist ttbar_JERUp("ttbar_JERUp", "ttbar_JERUp", RooArgList(x), h_mX_SR_ttbar_JERp1);
    RooDataHist ttbar_JERDown("ttbar_JERDown", "ttbar_JERDown", RooArgList(x), h_mX_SR_ttbar_JERm1);
  
    RooDataHist ttbar_BCUp("ttbar_BCUp", "ttbar_BCUp", RooArgList(x), h_mX_SR_ttbar_upBC);
    RooDataHist ttbar_BCDown("ttbar_BCDown", "ttbar_BCDown", RooArgList(x), h_mX_SR_ttbar_downBC);
  
    RooDataHist ttbar_LUp("ttbar_LUp", "ttbar_LUp", RooArgList(x), h_mX_SR_ttbar_upL);
    RooDataHist ttbar_LDown("ttbar_LDown", "ttbar_LDown", RooArgList(x), h_mX_SR_ttbar_downL);
  
    RooDataHist ttbar_TrigCSVUp("ttbar_TrigCSVUp", "ttbar_TrigCSVUp", RooArgList(x), h_mX_SR_ttbar_upTrigCSV);
    RooDataHist ttbar_TrigCSVDown("ttbar_TrigCSVDown", "ttbar_TrigCSVDown", RooArgList(x), h_mX_SR_ttbar_downTrigCSV);
  
    RooDataHist ttbar_TrigpTUp("ttbar_TrigpTUp", "ttbar_TrigpTUp", RooArgList(x), h_mX_SR_ttbar_upTrigpT);
    RooDataHist ttbar_TrigpTDown("ttbar_TrigpTDown", "ttbar_TrigpTDown", RooArgList(x), h_mX_SR_ttbar_downTrigpT);
  
    RooWorkspace *w=new RooWorkspace("HbbHbb");
    w->import(ttbar);
    w->import(ttbar_JECUp);
    w->import(ttbar_JECDown);
    w->import(ttbar_JERUp);
    w->import(ttbar_JERDown);
    w->import(ttbar_BCUp);
    w->import(ttbar_BCDown);
    w->import(ttbar_LUp);
    w->import(ttbar_LDown);
    w->import(ttbar_TrigCSVUp);
    w->import(ttbar_TrigCSVDown);
    w->import(ttbar_TrigpTUp);
    w->import(ttbar_TrigpTDown);
    w->SaveAs("w_ttbar.root");
  */}
  
  // Draw the histograms on a TCanvas
  
  if (drawPlots)
  {
  
    h_mX_SR_ttbar->Rebin(rebin);
    h_mX_SR_ttbar_JECp1->Rebin(rebin);
    h_mX_SR_ttbar_JECm1->Rebin(rebin);
    h_mX_SR_ttbar_JERp1->Rebin(rebin);
    h_mX_SR_ttbar_JERm1->Rebin(rebin);
    h_mX_SR_ttbar_upBC->Rebin(rebin);
    h_mX_SR_ttbar_downBC->Rebin(rebin);
    h_mX_SR_ttbar_upL->Rebin(rebin);
    h_mX_SR_ttbar_downL->Rebin(rebin);
    h_mX_SR_ttbar_upTrigCSV->Rebin(rebin);
    h_mX_SR_ttbar_downTrigCSV->Rebin(rebin);
    h_mX_SR_ttbar_upTrigpT->Rebin(rebin);
    h_mX_SR_ttbar_downTrigpT->Rebin(rebin);
    
    h_mX_SR_ttbar->SetLineColor(kBlack); h_mX_SR_ttbar->SetLineWidth(2);
    h_mX_SR_ttbar_JECp1->SetLineColor(kRed);
    h_mX_SR_ttbar_JECm1->SetLineColor(kRed+2);
    h_mX_SR_ttbar_JERp1->SetLineColor(kMagenta);
    h_mX_SR_ttbar_JERm1->SetLineColor(kMagenta+2);
    h_mX_SR_ttbar_upBC->SetLineColor(kGreen);
    h_mX_SR_ttbar_downBC->SetLineColor(kGreen+3);
    h_mX_SR_ttbar_upL->SetLineColor(kCyan);
    h_mX_SR_ttbar_downL->SetLineColor(kCyan+2);
    h_mX_SR_ttbar_upTrigCSV->SetLineColor(kBlue);
    h_mX_SR_ttbar_downTrigCSV->SetLineColor(kBlue+2);
    h_mX_SR_ttbar_upTrigpT->SetLineColor(kBlue-1);
    h_mX_SR_ttbar_downTrigpT->SetLineColor(kBlue-2);
    
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(000000000);
    
    TCanvas *c_ttbar=new TCanvas("c_ttbar", "c_ttbar", 700, 700);
    h_mX_SR_ttbar->Draw("hist"); h_mX_SR_ttbar->SetMaximum(10); h_mX_SR_ttbar->GetXaxis()->SetRangeUser(200, 1200);
    h_mX_SR_ttbar->Draw("ep3 same");
    h_mX_SR_ttbar_JECp1->Draw("hist same");
    h_mX_SR_ttbar_JECm1->Draw("hist same");
    h_mX_SR_ttbar_JERp1->Draw("hist same");
    h_mX_SR_ttbar_JERm1->Draw("hist same");
    h_mX_SR_ttbar_upBC->Draw("hist same");
    h_mX_SR_ttbar_downBC->Draw("hist same");
    h_mX_SR_ttbar_upL->Draw("hist same");
    h_mX_SR_ttbar_downL->Draw("hist same");
    h_mX_SR_ttbar_upTrigCSV->Draw("hist same");
    h_mX_SR_ttbar_downTrigCSV->Draw("hist same");
    h_mX_SR_ttbar_upTrigpT->Draw("hist same");
    h_mX_SR_ttbar_downTrigpT->Draw("hist same");
    TLegend *leg=new TLegend(0.5, 0.5, 0.9, 0.9);
    leg->AddEntry(h_mX_SR_ttbar, "Baseline");
    leg->AddEntry(h_mX_SR_ttbar_JECp1, "JEC +1 #sigma");
    leg->AddEntry(h_mX_SR_ttbar_JECm1, "JEC -1 #sigma");
    leg->AddEntry(h_mX_SR_ttbar_JERp1, "JER +1 #sigma");
    leg->AddEntry(h_mX_SR_ttbar_JERm1, "JER -1 #sigma");
    leg->AddEntry(h_mX_SR_ttbar_upBC, "SF_{bc} +1 #sigma");
    leg->AddEntry(h_mX_SR_ttbar_downBC, "SF_{bc} -1 #sigma");
    leg->AddEntry(h_mX_SR_ttbar_upL, "SF_{l} +1 #sigma");
    leg->AddEntry(h_mX_SR_ttbar_downL, "SF_{l} -1 #sigma");
    leg->AddEntry(h_mX_SR_ttbar_upTrigCSV, "Trig SF (CSV) +1 #sigma");
    leg->AddEntry(h_mX_SR_ttbar_downTrigCSV, "Trig SF (CSV) -1 #sigma");
    leg->AddEntry(h_mX_SR_ttbar_upTrigpT, "Trig SF (p_{T}) +1 #sigma");
    leg->AddEntry(h_mX_SR_ttbar_downTrigpT, "Trig SF (p_{T}) -1 #sigma");
    leg->Draw();
    c_ttbar->SaveAs("c_ttbar.png");
  }
  
  return 1;
}

