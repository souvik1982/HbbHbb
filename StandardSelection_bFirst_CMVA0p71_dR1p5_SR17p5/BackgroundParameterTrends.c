#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <iostream>
#include <TStyle.h>
#include <TPaveText.h>
#include <THStack.h>
#include <TF1.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooRealVar.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooArgList.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooChebychev.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooDataHist.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooExtendPdf.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooWorkspace.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooPlot.h"

double rebin=4;

std::string tags="MMMM_nominal"; // MMMM

double quad(double a, double b, double c=0, double d=0, double e=0, double f=0)
{
  return pow(a*a+b*b+c*c+d*d+e*e+f*f, 0.5);
}

struct Params
{
  double bg_p0;
  double bg_p1;
  double bg_p2;
  double bg_p3;
  double bg_p0_err;
  double bg_p1_err;
  double bg_p2_err;
  double bg_p3_err;
};

TH1F* returnSR(TFile *f)
{
  return (TH1F*)f->Get("h_mX_SR");
}

TH1F* returnSB(TFile *f)
{
  TH1F *h_mX_CR=(TH1F*)f->Get("h_mX_CR2");
  h_mX_CR->Add((TH1F*)f->Get("h_mX_CR4"));
  TH1F *h_mX_SR=(TH1F*)f->Get("h_mX_SR");
  h_mX_CR->Scale(h_mX_SR->GetSumOfWeights()/h_mX_CR->GetSumOfWeights());
  return h_mX_CR;
}

void fitBackground(TH1F* h_SR, TH1F*h_SB, std::string id, Params &params_SR, Params &params_SB)
{
  RooRealVar *x;
  RooRealVar *p0_SR, *p1_SR, *p2_SR;
  RooRealVar *p0_SB, *p1_SB, *p2_SB;
  
  std::string mH;
  double rangeLo=1, rangeHi=-1;
  if (id=="b")
  {
    mH="90";
    rangeLo=250; rangeHi=1200;
    
    p0_SR=new RooRealVar("p0_SR", "p0_SR", 300., 700.);
    p1_SR=new RooRealVar("p1_SR", "p1_SR", 40., 100.1);
    p2_SR=new RooRealVar("p2_SR", "p2_SR", 0.1, 10.1);
    
    p0_SB=new RooRealVar("p0_SB", "p0_SB", 300., 700.);
    p1_SB=new RooRealVar("p1_SB", "p1_SB", 40., 100.1);
    p2_SB=new RooRealVar("p2_SB", "p2_SB", 0.1, 10.1);
  }
  else if (id=="c")
  {
    mH="160";
    rangeLo=350; rangeHi=1200;
    
    p0_SR=new RooRealVar("p0_SR", "p0_SR", 500., 800.);
    p1_SR=new RooRealVar("p1_SR", "p1_SR", 40., 150.1);
    p2_SR=new RooRealVar("p2_SR", "p2_SR", 0.1, 10.1);
    
    p0_SB=new RooRealVar("p0_SB", "p0_SB", 500., 800.);
    p1_SB=new RooRealVar("p1_SB", "p1_SB", 40., 150.1);
    p2_SB=new RooRealVar("p2_SB", "p2_SB", 0.1, 10.1);
  }
  else if (id=="d")
  {
    mH="107.5";
    rangeLo=250; rangeHi=1200;
    
    p0_SR=new RooRealVar("p0_SR", "p0_SR", 300., 700.);
    p1_SR=new RooRealVar("p1_SR", "p1_SR", 40., 100.1);
    p2_SR=new RooRealVar("p2_SR", "p2_SR", 0.1, 10.1);
    
    p0_SB=new RooRealVar("p0_SB", "p0_SB", 300., 700.);
    p1_SB=new RooRealVar("p1_SB", "p1_SB", 40., 100.1);
    p2_SB=new RooRealVar("p2_SB", "p2_SB", 0.1, 10.1);
  }
  else if (id=="e")
  {
    mH="142.5";
    rangeLo=350; rangeHi=1200;
    
    p0_SR=new RooRealVar("p0_SR", "p0_SR", 450., 700.);
    p1_SR=new RooRealVar("p1_SR", "p1_SR", 40., 100.1);
    p2_SR=new RooRealVar("p2_SR", "p2_SR", 0.1, 10.1);
    
    p0_SB=new RooRealVar("p0_SB", "p0_SB", 450., 700.);
    p1_SB=new RooRealVar("p1_SB", "p1_SB", 40., 100.1);
    p2_SB=new RooRealVar("p2_SB", "p2_SB", 0.1, 10.1);
  }
  x=new RooRealVar("x", ("m_{X} at m_{H} = "+mH+" GeV").c_str(), rangeLo-50, rangeHi+50);
  GaussExp fit_SR("fit_SR", ("SR fit in Sideband "+id).c_str(), *x, *p0_SR, *p1_SR, *p2_SR);
  GaussExp fit_SB("fit_SB", ("SB fit in Sideband "+id).c_str(), *x, *p0_SB, *p1_SB, *p2_SB);
  RooDataHist hist_SR("hist_SR", ("SR data in Sideband "+id).c_str(), RooArgList(*x), h_SR);
  RooDataHist hist_SB("hist_SB", ("SB data in Sideband "+id).c_str(), RooArgList(*x), h_SB);
  RooFitResult *r_fit_SR=fit_SR.fitTo(hist_SR, RooFit::Range(rangeLo, rangeHi), RooFit::Save());
  RooFitResult *r_fit_SB=fit_SB.fitTo(hist_SB, RooFit::Range(rangeLo, rangeHi), RooFit::Save());
  params_SR.bg_p0=p0_SR->getVal(); params_SR.bg_p0_err=p0_SR->getError();
  params_SR.bg_p1=p1_SR->getVal(); params_SR.bg_p1_err=p1_SR->getError();
  params_SR.bg_p2=p2_SR->getVal(); params_SR.bg_p2_err=p2_SR->getError();
  params_SB.bg_p0=p0_SB->getVal(); params_SB.bg_p0_err=p0_SB->getError();
  params_SB.bg_p1=p1_SB->getVal(); params_SB.bg_p1_err=p1_SB->getError();
  params_SB.bg_p2=p2_SB->getVal(); params_SB.bg_p2_err=p2_SB->getError();
  RooPlot *plot=x->frame();
  hist_SR.plotOn(plot);
  hist_SB.plotOn(plot);
  fit_SR.plotOn(plot, RooFit::VisualizeError(*r_fit_SR, 1), RooFit::FillColor(kCyan));
  fit_SB.plotOn(plot, RooFit::VisualizeError(*r_fit_SB, 1), RooFit::FillColor(kOrange));
  hist_SR.plotOn(plot, RooFit::LineColor(kBlue), RooFit::MarkerColor(kBlue));
  hist_SB.plotOn(plot, RooFit::LineColor(kRed), RooFit::MarkerColor(kRed));
  fit_SR.plotOn(plot, RooFit::LineColor(kBlue));
  fit_SB.plotOn(plot, RooFit::LineColor(kRed));
  TCanvas *c_mX_Sideband=new TCanvas("c_mX_Sideband", "c_mX_Sideband", 700, 700);
  plot->Draw();
  c_mX_Sideband->SaveAs(("c_mX_Sideband_"+id+".png").c_str());
  delete c_mX_Sideband;
  delete plot;
  delete r_fit_SR;
  delete r_fit_SB;
  delete x;
  delete p0_SR; delete p1_SR; delete p2_SR;
  delete p0_SB; delete p1_SB; delete p2_SB;
} 

void BackgroundParameterTrends()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  gSystem->Load("../PDFs/GaussExp_cxx.so");
  
  TFile *f_data_a=new TFile((tags+"/a_KinFit/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TFile *f_data_b=new TFile((tags+"/b_KinFit/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TFile *f_data_c=new TFile((tags+"/c_KinFit/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TFile *f_data_d=new TFile((tags+"/d_KinFit/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TFile *f_data_e=new TFile((tags+"/e_KinFit/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  
  TFile *f_ttbar_a=new TFile((tags+"/a_KinFit/Histograms_TTJets_mX.root").c_str());
  TFile *f_ttbar_b=new TFile((tags+"/b_KinFit/Histograms_TTJets_mX.root").c_str());
  TFile *f_ttbar_c=new TFile((tags+"/c_KinFit/Histograms_TTJets_mX.root").c_str());
  TFile *f_ttbar_d=new TFile((tags+"/d_KinFit/Histograms_TTJets_mX.root").c_str());
  TFile *f_ttbar_e=new TFile((tags+"/e_KinFit/Histograms_TTJets_mX.root").c_str());
  
  TH1F *h_data_b_SR=returnSR(f_data_b);
  TH1F *h_data_b_SB=returnSB(f_data_b);
  TH1F *h_data_c_SR=returnSR(f_data_c);
  TH1F *h_data_c_SB=returnSB(f_data_c);
  TH1F *h_data_d_SR=returnSR(f_data_d);
  TH1F *h_data_d_SB=returnSB(f_data_d);
  TH1F *h_data_e_SR=returnSR(f_data_e);
  TH1F *h_data_e_SB=returnSB(f_data_e);
  
  TH1F *h_ttbar_b_SR=(TH1F*)f_ttbar_b->Get("h_mX_SR_ttbar");
  TH1F *h_ttbar_b_SB=(TH1F*)f_ttbar_b->Get("h_mX_SB_ttbar");
  TH1F *h_ttbar_c_SR=(TH1F*)f_ttbar_c->Get("h_mX_SR_ttbar");
  TH1F *h_ttbar_c_SB=(TH1F*)f_ttbar_c->Get("h_mX_SB_ttbar");
  TH1F *h_ttbar_d_SR=(TH1F*)f_ttbar_d->Get("h_mX_SR_ttbar");
  TH1F *h_ttbar_d_SB=(TH1F*)f_ttbar_d->Get("h_mX_SB_ttbar");
  TH1F *h_ttbar_e_SR=(TH1F*)f_ttbar_e->Get("h_mX_SR_ttbar");
  TH1F *h_ttbar_e_SB=(TH1F*)f_ttbar_e->Get("h_mX_SB_ttbar");
  /*
  h_data_b_SR->Add(h_ttbar_b_SR, -1);
  h_data_b_SB->Add(h_ttbar_b_SB, -1);
  h_data_c_SR->Add(h_ttbar_c_SR, -1);
  h_data_c_SB->Add(h_ttbar_c_SB, -1);
  h_data_d_SR->Add(h_ttbar_d_SR, -1);
  h_data_d_SB->Add(h_ttbar_d_SB, -1);
  h_data_e_SR->Add(h_ttbar_e_SR, -1);
  h_data_e_SB->Add(h_ttbar_e_SB, -1);
  */
  h_data_b_SR->Rebin(rebin);
  h_data_b_SB->Rebin(rebin);
  h_data_c_SR->Rebin(6);
  h_data_c_SB->Rebin(6);
  h_data_d_SR->Rebin(rebin);
  h_data_d_SB->Rebin(rebin);
  h_data_e_SR->Rebin(rebin);
  h_data_e_SB->Rebin(rebin);
  
  Params par_b_SR, par_b_SB;
  Params par_c_SR, par_c_SB;
  Params par_d_SR, par_d_SB;
  Params par_e_SR, par_e_SB;
  
  fitBackground(h_data_b_SR, h_data_b_SB, "b", par_b_SR, par_b_SB);
  fitBackground(h_data_c_SR, h_data_c_SB, "c", par_c_SR, par_c_SB);
  fitBackground(h_data_d_SR, h_data_d_SB, "d", par_d_SR, par_d_SB);
  fitBackground(h_data_e_SR, h_data_e_SB, "e", par_e_SR, par_e_SB);
  
  const unsigned int nPoints=4;
  double mass[nPoints]={90., 107.5, 142.5, 160.};
  double p0_SR[nPoints],  p1_SR[nPoints], p2_SR[nPoints];
  double p0_SB[nPoints],  p1_SB[nPoints], p2_SB[nPoints];
  double p0_SR_err[nPoints],  p1_SR_err[nPoints], p2_SR_err[nPoints];
  double p0_SB_err[nPoints],  p1_SB_err[nPoints], p2_SB_err[nPoints];
  double p0_ratio[nPoints], p1_ratio[nPoints], p2_ratio[nPoints];
  double p0_errorY[nPoints], p1_errorY[nPoints], p2_errorY[nPoints];
  double p0_errorX[nPoints], p1_errorX[nPoints], p2_errorX[nPoints];
  
  p0_SR[0]=par_b_SR.bg_p0; p0_SR_err[0]=par_b_SR.bg_p0_err; 
  p0_SR[1]=par_d_SR.bg_p0; p0_SR_err[1]=par_d_SR.bg_p0_err;
  p0_SR[2]=par_e_SR.bg_p0; p0_SR_err[2]=par_e_SR.bg_p0_err;
  p0_SR[3]=par_c_SR.bg_p0; p0_SR_err[3]=par_c_SR.bg_p0_err;
  
  p0_SB[0]=par_b_SB.bg_p0; p0_SB_err[0]=par_b_SB.bg_p0_err;
  p0_SB[1]=par_d_SB.bg_p0; p0_SB_err[1]=par_d_SB.bg_p0_err;
  p0_SB[2]=par_e_SB.bg_p0; p0_SB_err[2]=par_e_SB.bg_p0_err;
  p0_SB[3]=par_c_SB.bg_p0; p0_SB_err[3]=par_c_SB.bg_p0_err;
  
  p1_SR[0]=par_b_SR.bg_p1; p1_SR_err[0]=par_b_SR.bg_p1_err; 
  p1_SR[1]=par_d_SR.bg_p1; p1_SR_err[1]=par_d_SR.bg_p1_err;
  p1_SR[2]=par_e_SR.bg_p1; p1_SR_err[2]=par_e_SR.bg_p1_err;
  p1_SR[3]=par_c_SR.bg_p1; p1_SR_err[3]=par_c_SR.bg_p1_err;
  
  p1_SB[0]=par_b_SB.bg_p1; p1_SB_err[0]=par_b_SB.bg_p1_err;
  p1_SB[1]=par_d_SB.bg_p1; p1_SB_err[1]=par_d_SB.bg_p1_err;
  p1_SB[2]=par_e_SB.bg_p1; p1_SB_err[2]=par_e_SB.bg_p1_err;
  p1_SB[3]=par_c_SB.bg_p1; p1_SB_err[3]=par_c_SB.bg_p1_err;
  
  p2_SR[0]=par_b_SR.bg_p2; p2_SR_err[0]=par_b_SR.bg_p2_err; 
  p2_SR[1]=par_d_SR.bg_p2; p2_SR_err[1]=par_d_SR.bg_p2_err;
  p2_SR[2]=par_e_SR.bg_p2; p2_SR_err[2]=par_e_SR.bg_p2_err;
  p2_SR[3]=par_c_SR.bg_p2; p2_SR_err[3]=par_c_SR.bg_p2_err;
  
  p2_SB[0]=par_b_SB.bg_p2; p2_SB_err[0]=par_b_SB.bg_p2_err;
  p2_SB[1]=par_d_SB.bg_p2; p2_SB_err[1]=par_d_SB.bg_p2_err;
  p2_SB[2]=par_e_SB.bg_p2; p2_SB_err[2]=par_e_SB.bg_p2_err;
  p2_SB[3]=par_c_SB.bg_p2; p2_SB_err[3]=par_c_SB.bg_p2_err;
  
  
  for (unsigned int i=0; i<nPoints; ++i)
  {
    p0_ratio[i]=p0_SR[i]/p0_SB[i];
    p0_errorY[i]=p0_ratio[i]*quad(p0_SR_err[i]/p0_SR[i], p0_SB_err[i]/p0_SB[i]);
    p0_errorX[i]=0;
    
    p1_ratio[i]=p1_SR[i]/p1_SB[i];
    p1_errorY[i]=p1_ratio[i]*quad(p1_SR_err[i]/p1_SR[i], p1_SB_err[i]/p1_SB[i]);
    p1_errorX[i]=0;
    
    p2_ratio[i]=p2_SR[i]/p2_SB[i];
    p2_errorY[i]=p2_ratio[i]*quad(p2_SR_err[i]/p2_SR[i], p2_SB_err[i]/p2_SB[i]);
    p2_errorX[i]=0;
  }
  
  TGraphErrors *g_p0_ratio=new TGraphErrors(nPoints, mass, p0_ratio, p0_errorX, p0_errorY);
  g_p0_ratio->SetTitle("SR/SB ratio for the mean of the Gaussian core; m_{H} (GeV); SR/SB Ratio of Mean");
  TCanvas *c_p0_ratio=new TCanvas("c_p0_ratio", "c_p0_ratio", 700, 700);
  g_p0_ratio->SetMinimum(0.); g_p0_ratio->SetMaximum(2.);
  g_p0_ratio->Draw("A*");
  TF1 *f_p0_ratio=new TF1("f_p0_ratio", "pol1");
  g_p0_ratio->Fit(f_p0_ratio);
  c_p0_ratio->SaveAs("c_p0_ratio.png");
  double p0_ratioAt125=f_p0_ratio->Eval(125.);
  double p0_errorAt125=(p0_errorY[0]+p0_errorY[1]+p0_errorY[2]+p0_errorY[3])/4.;
  std::cout<<"p0_ratioAt125 = "<<p0_ratioAt125<<" +- "<<p0_errorAt125<<std::endl;
  
  TGraphErrors *g_p1_ratio=new TGraphErrors(nPoints, mass, p1_ratio, p1_errorX, p1_errorY);
  g_p1_ratio->SetTitle("SR/SB ratio for the std dev of the Gaussian core; m_{H} (GeV); SR/SB Ratio of Std Dev");
  TCanvas *c_p1_ratio=new TCanvas("c_p1_ratio", "c_p1_ratio", 700, 700);
  g_p1_ratio->SetMinimum(0.); g_p1_ratio->SetMaximum(2.);
  g_p1_ratio->Draw("A*");
  TF1 *f_p1_ratio=new TF1("f_p1_ratio", "pol1");
  g_p1_ratio->Fit(f_p1_ratio);
  c_p1_ratio->SaveAs("c_p1_ratio.png");
  double p1_ratioAt125=f_p1_ratio->Eval(125.);
  double p1_errorAt125=(p1_errorY[0]+p1_errorY[1]+p1_errorY[2]+p1_errorY[3])/4.;
  std::cout<<"p1_ratioAt125 = "<<p1_ratioAt125<<" +- "<<p1_errorAt125<<std::endl;
  
  TGraphErrors *g_p2_ratio=new TGraphErrors(nPoints, mass, p2_ratio, p2_errorX, p2_errorY);
  g_p2_ratio->SetTitle("SR/SB ratio for the exponent of the tail; m_{H} (GeV); SR/SB Ratio of Exponent");
  TCanvas *c_p2_ratio=new TCanvas("c_p2_ratio", "c_p2_ratio", 700, 700);
  g_p2_ratio->SetMinimum(0.); g_p2_ratio->SetMaximum(2.);
  g_p2_ratio->Draw("A*");
  TF1 *f_p2_ratio=new TF1("f_p2_ratio", "pol1");
  g_p2_ratio->Fit(f_p2_ratio);
  c_p2_ratio->SaveAs("c_p2_ratio.png");
  double p2_ratioAt125=f_p2_ratio->Eval(125.);
  double p2_errorAt125=(p2_errorY[0]+p2_errorY[1]+p2_errorY[2]+p2_errorY[3])/4.;
  std::cout<<"p2_ratioAt125 = "<<p2_ratioAt125<<" +- "<<p2_errorAt125<<std::endl;
  
}
  

