#include <TH1F.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <algorithm>

#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooRealVar.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooDataHist.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooWorkspace.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooPlot.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooFitResult.h"

// #include "/Users/souvik/HbbHbb/Analysis/PDFs/GaussExp.h"

bool doHistograms=false;
bool doFit=true;

double SR_lo=320.;
double SR_hi=1200.;

std::string itoa(int i) 
{
  char res[4];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

double quad(double a, double b, double c=0, double d=0, double e=0, double f=0, double g=0, double h=0, double i=0, double j=0, double k=0)
{
  return pow(a*a+b*b+c*c+d*d+e*e+f*f+g*g+h*h+i*i+j*j+k*k, 0.5);
}

struct Params
{
  double sg_p0;
  double sg_p1;
  double sg_p2;
  double sg_p3;
  double sg_p0_err;
  double sg_p1_err;
  double sg_p2_err;
  double sg_p3_err;
};

double lnN(double b, double a, double c)
{
  double err=0;
  if (b>0) err=1.+fabs(a-c)/(2.*b);
  return err;
}

RooPlot* fitttbar(RooDataHist r, int color, Params &params, RooWorkspace *w)
{
  RooRealVar *x=new RooRealVar("x", "m_{X} (GeV)", SR_lo, SR_hi);
  // RooRealVar *x=(RooRealVar*)r.FindObject("x");
  std::cout<<"x="<<x<<std::endl;
  
  RooRealVar *tt_p0=new RooRealVar("tt_p0", "tt_p0", 400, 550);
  RooRealVar *tt_p1=new RooRealVar("tt_p1", "tt_p1", 30, 90);
  RooRealVar *tt_p2=new RooRealVar("tt_p2", "tt_p2", 0.0, 1.0);
  
  GaussExp ttbar("ttbar", "ttbar", *x, *tt_p0, *tt_p1, *tt_p2);
  RooFitResult *r_ttbar=ttbar.fitTo(r, RooFit::Range(SR_lo, SR_hi), RooFit::Save());
  
  params.sg_p0=tt_p0->getVal(); params.sg_p0_err=tt_p0->getError();
  params.sg_p1=tt_p1->getVal(); params.sg_p1_err=tt_p1->getError();
  params.sg_p2=tt_p2->getVal(); params.sg_p2_err=tt_p2->getError();
  RooPlot *plot=x->frame();
  r.plotOn(plot);
  if (color==kBlack) 
  {
    // ttbar.plotOn(plot, RooFit::VisualizeError(*r_ttbar, 2), RooFit::FillColor(kGreen));
    ttbar.plotOn(plot, RooFit::VisualizeError(*r_ttbar, 1), RooFit::FillColor(kCyan));
  }
  ttbar.plotOn(plot, RooFit::LineColor(color));
  r.plotOn(plot, RooFit::MarkerColor(color));
  if (color==kBlack)
  {
    RooRealVar ttbar_p0("ttbar_p0", "ttbar_p0", tt_p0->getVal());
    RooRealVar ttbar_p1("ttbar_p1", "ttbar_p1", tt_p1->getVal());
    RooRealVar ttbar_p2("ttbar_p2", "ttbar_p2", tt_p2->getVal());
    GaussExp ttbar_fixed("ttbar", "ttbar", *x, ttbar_p0, ttbar_p1, ttbar_p2);
    w->import(ttbar_fixed);
  }
  
  return plot;
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
  
  h_mX_SR_ttbar->SetTitle("t#bar{t} Distribution in Signal Region; m_{X} (GeV); Events / 40 GeV");
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  
  // Put the histograms into a RooFit workspace and save it
  RooRealVar x("x", "m_{X} (GeV)", SR_lo, SR_hi);
  x.setBins(88);
  RooWorkspace *w=new RooWorkspace("HbbHbb");
  
  RooDataHist r_ttbar("ttbar", "ttbar", RooArgList(x), h_mX_SR_ttbar);
  
  RooDataHist r_ttbar_JECUp("ttbar_ttbar_JECUp", "ttbar_ttbar_JECUp", RooArgList(x), h_mX_SR_ttbar_JECp1);
  RooDataHist r_ttbar_JECDown("ttbar_ttbar_JECDown", "ttbar_ttbar_JECDown", RooArgList(x), h_mX_SR_ttbar_JECm1);
  
  RooDataHist r_ttbar_JERUp("ttbar_ttbar_JERUp", "ttbar_ttbar_JERUp", RooArgList(x), h_mX_SR_ttbar_JERp1);
  RooDataHist r_ttbar_JERDown("ttbar_ttbar_JERDown", "ttbar_ttbar_JERDown", RooArgList(x), h_mX_SR_ttbar_JERm1);
  
  RooDataHist r_ttbar_BCUp("ttbar_ttbar_BCUp", "ttbar_ttbar_BCUp", RooArgList(x), h_mX_SR_ttbar_upBC);
  RooDataHist r_ttbar_BCDown("ttbar_ttbar_BCDown", "ttbar_ttbar_BCDown", RooArgList(x), h_mX_SR_ttbar_downBC);
  
  RooDataHist r_ttbar_LUp("ttbar_ttbar_LUp", "ttbar_ttbar_LUp", RooArgList(x), h_mX_SR_ttbar_upL);
  RooDataHist r_ttbar_LDown("ttbar_ttbar_LDown", "ttbar_ttbar_LDown", RooArgList(x), h_mX_SR_ttbar_downL);
  
  RooDataHist r_ttbar_TrigCSVUp("ttbar_ttbar_TrigCSVUp", "ttbar_ttbar_TrigCSVUp", RooArgList(x), h_mX_SR_ttbar_upTrigCSV);
  RooDataHist r_ttbar_TrigCSVDown("ttbar_ttbar_TrigCSVDown", "ttbar_ttbar_TrigCSVDown", RooArgList(x), h_mX_SR_ttbar_downTrigCSV);
  
  RooDataHist r_ttbar_TrigpTUp("ttbar_ttbar_TrigpTUp", "ttbar_ttbar_TrigpTUp", RooArgList(x), h_mX_SR_ttbar_upTrigpT);
  RooDataHist r_ttbar_TrigpTDown("ttbar_ttbar_TrigpTDown", "ttbar_ttbar_TrigpTDown", RooArgList(x), h_mX_SR_ttbar_downTrigpT);
    
  if (doHistograms)
  {
    w->import(r_ttbar);
    w->import(r_ttbar_JECUp);
    w->import(r_ttbar_JECDown);
    w->import(r_ttbar_JERUp);
    w->import(r_ttbar_JERDown);
    w->import(r_ttbar_BCUp);
    w->import(r_ttbar_BCDown);
    w->import(r_ttbar_LUp);
    w->import(r_ttbar_LDown);
    w->import(r_ttbar_TrigCSVUp);
    w->import(r_ttbar_TrigCSVDown);
    w->import(r_ttbar_TrigpTUp);
    w->import(r_ttbar_TrigpTDown);
    
    for (int i=0; i<h_mX_SR_ttbar->GetNbinsX(); ++i)
    {
      double binContent=h_mX_SR_ttbar->GetBinContent(i);
      double binError=h_mX_SR_ttbar->GetBinError(i);
      // std::cout<<"bin = "<<i<<", binContent = "<<binContent<<", binError = "<<binError<<std::endl;
      
      if (binContent>0)
      {
        TH1F *h_mX_SR_ttbar_statUp=(TH1F*)h_mX_SR_ttbar->Clone("h_mX_SR_ttbar_statUp");
        TH1F *h_mX_SR_ttbar_statDown=(TH1F*)h_mX_SR_ttbar->Clone("h_mX_SR_ttbar_statDown");
        h_mX_SR_ttbar_statUp->SetBinContent(i, binContent+binError);
        h_mX_SR_ttbar_statDown->SetBinContent(i, binContent-binError);
        
        std::string histname="ttbar_StatBin"+itoa(i);
        std::cout<<"StatBin"+itoa(i)<<"   shape    -    -    1"<<std::endl;
        RooDataHist r_ttbar_StatBinUp((histname+"Up").c_str(), (histname+"Up").c_str(), RooArgList(x), h_mX_SR_ttbar_statUp);
        RooDataHist r_ttbar_StatBinDown((histname+"Down").c_str(), (histname+"Down").c_str(), RooArgList(x), h_mX_SR_ttbar_statDown);
        w->import(r_ttbar_StatBinUp);
        w->import(r_ttbar_StatBinDown);
      }
    }
    
    for (int i=0; i<h_mX_SR_ttbar->GetNbinsX(); ++i)
    {
      double binContent=h_mX_SR_ttbar->GetBinContent(i);
      if (binContent>0) std::cout<<"StatBin"+itoa(i)<<"   shape    -    -    1"<<std::endl;
    }
    
    TCanvas *c_ttbar=new TCanvas("c_ttbar", "c_ttbar", 700, 700);
    h_mX_SR_ttbar->Draw("hist"); h_mX_SR_ttbar->SetMaximum(h_mX_SR_ttbar->GetMaximum()*1.5); h_mX_SR_ttbar->GetXaxis()->SetRangeUser(200, 1200);
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
  
  if (doFit)
  {
    gSystem->Load("/Users/souvik/HbbHbb/Analysis/PDFs/GaussExp_cxx.so");
    
    Params par, par_JECp1, par_JECm1, par_JERp1, par_JERm1, par_upBC, par_downBC, par_upL, par_downL, par_upTrigCSV, par_downTrigCSV, par_upTrigpT, par_downTrigpT;
    RooPlot *plot=fitttbar(r_ttbar, kBlack, par, w);
    RooPlot *plot_JECp1=fitttbar(r_ttbar_JECUp, kRed, par_JECp1, w);
    RooPlot *plot_JECm1=fitttbar(r_ttbar_JECDown, kRed+2, par_JECm1, w);
    RooPlot *plot_JERp1=fitttbar(r_ttbar_JERUp, kMagenta, par_JERp1, w);
    RooPlot *plot_JERm1=fitttbar(r_ttbar_JERDown, kMagenta+2, par_JERm1, w);
    RooPlot *plot_upBC=fitttbar(r_ttbar_BCUp, kGreen, par_upBC, w);
    RooPlot *plot_downBC=fitttbar(r_ttbar_BCDown, kGreen+3, par_downBC, w);
    RooPlot *plot_upL=fitttbar(r_ttbar_LUp, kOrange, par_upL, w);
    RooPlot *plot_downL=fitttbar(r_ttbar_LDown, kOrange+2, par_downL, w);
    RooPlot *plot_upTrigCSV=fitttbar(r_ttbar_TrigCSVUp, kBlue, par_upTrigCSV, w);
    RooPlot *plot_downTrigCSV=fitttbar(r_ttbar_TrigCSVDown, kBlue+2, par_downTrigCSV, w);
    RooPlot *plot_upTrigpT=fitttbar(r_ttbar_TrigpTUp, kBlue+3, par_upTrigpT, w);
    RooPlot *plot_downTrigpT=fitttbar(r_ttbar_TrigpTDown, kBlue+4, par_downTrigpT, w);
    TCanvas *c_ttbar_fit=new TCanvas("c_ttbar_fit", "c_ttbar_fit", 700, 700);
    plot->Draw();
    plot_JECp1->Draw("same");
    plot_JECm1->Draw("same");
    plot_JERp1->Draw("same");
    plot_JERm1->Draw("same");
    plot_upBC->Draw("same");
    plot_downBC->Draw("same");
    plot_upL->Draw("same");
    plot_downL->Draw("same");
    plot_upTrigCSV->Draw("same");
    plot_downTrigCSV->Draw("same");
    plot_upTrigpT->Draw("same");
    plot_downTrigpT->Draw("same");
    c_ttbar_fit->SaveAs("ttbar/c_ttbar_fit.png");
    
    ofstream outfile;
    outfile.open("ttbar/ttbarSystematics.html");
    outfile<<"<html>"<<std::endl;
    outfile<<"<head>"<<std::endl;
    outfile<<"</head>"<<std::endl;
    outfile<<"<body>"<<std::endl;
    outfile<<"<h1> ttbar with Statistical and Systematic Uncertainties </h1>"<<std::endl;
    outfile<<"<br/><hr/>"<<std::endl;
    outfile<<"<img src='c_ttbar_fit.png'/> <br/>"<<std::endl;
    outfile<<"=== Baseline plot === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par.sg_p0<<" +- "<<par.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par.sg_p1<<" +- "<<par.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par.sg_p2<<" +- "<<par.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === JEC +1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_JECp1->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_JECp1.sg_p0<<" +- "<<par_JECp1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_JECp1.sg_p1<<" +- "<<par_JECp1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_JECp1.sg_p2<<" +- "<<par_JECp1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === JEC -1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_JECm1->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_JECm1.sg_p0<<" +- "<<par_JECm1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_JECm1.sg_p1<<" +- "<<par_JECm1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_JECm1.sg_p2<<" +- "<<par_JECm1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === JER +1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_JERp1->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_JERp1.sg_p0<<" +- "<<par_JERp1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_JERp1.sg_p1<<" +- "<<par_JERp1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_JERp1.sg_p2<<" +- "<<par_JERp1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === JER -1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_JERm1->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_JERm1.sg_p0<<" +- "<<par_JERm1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_JERm1.sg_p1<<" +- "<<par_JERm1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_JERm1.sg_p2<<" +- "<<par_JERm1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(bc) +1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_upBC->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_upBC.sg_p0<<" +- "<<par_upBC.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_upBC.sg_p1<<" +- "<<par_upBC.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_upBC.sg_p2<<" +- "<<par_upBC.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(bc) -1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_downBC->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_downBC.sg_p0<<" +- "<<par_downBC.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_downBC.sg_p1<<" +- "<<par_downBC.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_downBC.sg_p2<<" +- "<<par_downBC.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(l) +1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_upL->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_upL.sg_p0<<" +- "<<par_upL.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_upL.sg_p1<<" +- "<<par_upL.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_upL.sg_p2<<" +- "<<par_upL.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(l) -1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_downL->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_downL.sg_p0<<" +- "<<par_downL.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_downL.sg_p1<<" +- "<<par_downL.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_downL.sg_p2<<" +- "<<par_downL.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (CSV) +1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_upTrigCSV->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_upTrigCSV.sg_p0<<" +- "<<par_upTrigCSV.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_upTrigCSV.sg_p1<<" +- "<<par_upTrigCSV.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_upTrigCSV.sg_p2<<" +- "<<par_upTrigCSV.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (CSV) -1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_downTrigCSV->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_downTrigCSV.sg_p0<<" +- "<<par_downTrigCSV.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_downTrigCSV.sg_p1<<" +- "<<par_downTrigCSV.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_downTrigCSV.sg_p2<<" +- "<<par_downTrigCSV.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (pT) +1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_upTrigpT->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_upTrigpT.sg_p0<<" +- "<<par_upTrigpT.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_upTrigpT.sg_p1<<" +- "<<par_upTrigpT.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_upTrigpT.sg_p2<<" +- "<<par_upTrigpT.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (pT) -1 sigma === <br/>"<<std::endl;
    outfile<<" norm = "<<h_mX_SR_ttbar_downTrigpT->Integral(bin1, bin2)<<" <br/>"<<std::endl;
    outfile<<" ttbar_p0 = "<<par_downTrigpT.sg_p0<<" +- "<<par_downTrigpT.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p1 = "<<par_downTrigpT.sg_p1<<" +- "<<par_downTrigpT.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<" ttbar_p2 = "<<par_downTrigpT.sg_p2<<" +- "<<par_downTrigpT.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"=== === <br/>"<<std::endl;
    
    double ttbar_p0_errStat=par.sg_p0_err;
    double ttbar_p0_errSyst[]={par.sg_p0,
                               par_JECp1.sg_p0, par_JECm1.sg_p0,
                               par_JERp1.sg_p0, par_JERm1.sg_p0,
                               par_upBC.sg_p0, par_downBC.sg_p0,
                               par_upL.sg_p0, par_downL.sg_p0,
                               par_upTrigCSV.sg_p0, par_downTrigCSV.sg_p0,
                               par_upTrigpT.sg_p0, par_downTrigpT.sg_p0};
    double ttbar_p0_errSyst_min=par.sg_p0-(*std::min_element(ttbar_p0_errSyst, ttbar_p0_errSyst+13));
    double ttbar_p0_errSyst_max=(*std::max_element(ttbar_p0_errSyst, ttbar_p0_errSyst+13))-par.sg_p0;
    outfile<<" Uncertainty on ttbar_p0 = "<<par.sg_p0<<" +- "<<ttbar_p0_errStat/2.<<" (stat) - "<<ttbar_p0_errSyst_min<<" + "<<ttbar_p0_errSyst_max<<" (syst); -"<<quad(ttbar_p0_errStat/2., ttbar_p0_errSyst_min)<<"/+"<<quad(ttbar_p0_errStat/2., ttbar_p0_errSyst_max)<<" (total) <br/>"<<std::endl;
    
    double ttbar_p1_errStat=par.sg_p1_err;
    double ttbar_p1_errSyst[]={par.sg_p1,
                               par_JECp1.sg_p1, par_JECm1.sg_p1,
                               par_JERp1.sg_p1, par_JERm1.sg_p1,
                               par_upBC.sg_p1, par_downBC.sg_p1,
                               par_upL.sg_p1, par_downL.sg_p1,
                               par_upTrigCSV.sg_p1, par_downTrigCSV.sg_p1,
                               par_upTrigpT.sg_p1, par_downTrigpT.sg_p1};
    double ttbar_p1_errSyst_min=par.sg_p1-(*std::min_element(ttbar_p1_errSyst, ttbar_p1_errSyst+13));
    double ttbar_p1_errSyst_max=(*std::max_element(ttbar_p1_errSyst, ttbar_p1_errSyst+13))-par.sg_p1;
    outfile<<" Uncertainty on ttbar_p1 = "<<par.sg_p1<<" +- "<<ttbar_p1_errStat/2.<<" (stat) - "<<ttbar_p1_errSyst_min<<" + "<<ttbar_p1_errSyst_max<<" (syst); -"<<quad(ttbar_p1_errStat/2., ttbar_p1_errSyst_min)<<"/+"<<quad(ttbar_p1_errStat/2., ttbar_p1_errSyst_max)<<" (total) <br/>"<<std::endl;
    
    double ttbar_p2_errStat=par.sg_p2_err;
    double ttbar_p2_errSyst[]={par.sg_p2,
                               par_JECp1.sg_p2, par_JECm1.sg_p2,
                               par_JERp1.sg_p2, par_JERm1.sg_p2,
                               par_upBC.sg_p2, par_downBC.sg_p2,
                               par_upL.sg_p2, par_downL.sg_p2,
                               par_upTrigCSV.sg_p2, par_downTrigCSV.sg_p2,
                               par_upTrigpT.sg_p2, par_downTrigpT.sg_p2};
    double ttbar_p2_errSyst_min=par.sg_p2-(*std::min_element(ttbar_p2_errSyst, ttbar_p2_errSyst+13));
    double ttbar_p2_errSyst_max=(*std::max_element(ttbar_p2_errSyst, ttbar_p2_errSyst+13))-par.sg_p2;
    outfile<<" Uncertainty on ttbar_p2 = "<<par.sg_p2<<" +- "<<ttbar_p2_errStat/2.<<" (stat) - "<<ttbar_p2_errSyst_min<<" + "<<ttbar_p2_errSyst_max<<" (syst); -"<<quad(ttbar_p2_errStat/2., ttbar_p2_errSyst_min)<<"/+"<<quad(ttbar_p2_errStat/2., ttbar_p2_errSyst_max)<<" (total) <br/>"<<std::endl;
    
    outfile<<"=== === <br/>"<<std::endl;
    outfile<<"   Copy this into the datacard: <br/>"<<std::endl;
    outfile<<" <pre>"<<std::endl;
    outfile<<"JEC       lnN     "<<lnN(h_mX_SR_ttbar->Integral(bin1, bin2), h_mX_SR_ttbar_JECp1->Integral(bin1, bin2), h_mX_SR_ttbar_JECm1->Integral(bin1, bin2))<<    "     -"<<std::endl;
    outfile<<"JER       lnN     "<<lnN(h_mX_SR_ttbar->Integral(bin1, bin2), h_mX_SR_ttbar_JERp1->Integral(bin1, bin2), h_mX_SR_ttbar_JERm1->Integral(bin1, bin2))<<    "     -"<<std::endl;
    outfile<<"bTagSFbc  lnN     "<<lnN(h_mX_SR_ttbar->Integral(bin1, bin2), h_mX_SR_ttbar_upBC->Integral(bin1, bin2), h_mX_SR_ttbar_downBC->Integral(bin1, bin2))<<    "     -"<<std::endl;
    outfile<<"bTagSFl   lnN     "<<lnN(h_mX_SR_ttbar->Integral(bin1, bin2), h_mX_SR_ttbar_upL->Integral(bin1, bin2), h_mX_SR_ttbar_downL->Integral(bin1, bin2))<<      "     -"<<std::endl;
    outfile<<"trigSFCSV lnN     "<<lnN(h_mX_SR_ttbar->Integral(bin1, bin2), h_mX_SR_ttbar_upTrigCSV->Integral(bin1, bin2), h_mX_SR_ttbar_downTrigCSV->Integral(bin1, bin2))<<"     -"<<std::endl;
    outfile<<"trigSFpT  lnN     "<<lnN(h_mX_SR_ttbar->Integral(bin1, bin2), h_mX_SR_ttbar_upTrigpT->Integral(bin1, bin2), h_mX_SR_ttbar_downTrigpT->Integral(bin1, bin2))<<"     -"<<std::endl;
    outfile<<"ttbar_p0 param   "<<par.sg_p0<<" -"<<quad(ttbar_p0_errStat/2., ttbar_p0_errSyst_min)<<"/+"<<quad(ttbar_p0_errStat/2., ttbar_p0_errSyst_max)<<std::endl;
    outfile<<"ttbar_p1 param   "<<par.sg_p1<<" -"<<quad(ttbar_p1_errStat/2., ttbar_p1_errSyst_min)<<"/+"<<quad(ttbar_p1_errStat/2., ttbar_p1_errSyst_max)<<std::endl;
    outfile<<"ttbar_p2 param   "<<par.sg_p2<<" -"<<quad(ttbar_p2_errStat/2., ttbar_p2_errSyst_min)<<"/+"<<quad(ttbar_p2_errStat/2., ttbar_p2_errSyst_max)<<std::endl;
    outfile<<" </pre>"<<std::endl;
    outfile<<"=== === <br/>"<<std::endl;
    
    outfile<<"</body>"<<std::endl;
    outfile<<"</html>"<<std::endl;
  }
  
  w->SaveAs("w_ttbar.root");
  
  return 0;
}

