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
// #include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooFitResult.h"

// #include "/Users/souvik/HbbHbb/Analysis/PDFs/GaussExp.h"

double H_mass=125.0;
double mH_diff_cut=40.;
double mH_mean_cut=20.;

double rebin=1;
bool useRatioFit=false;
bool bReg=false;

std::string tags="MMMM_nominal"; // MMMM

double VR_lo=250.;
double VR_hi=1200.;
double SR_lo=350.; // 350 for MMMM_nominal and 400 for MMMMbar
double SR_hi=1200.;

double quad(double a, double b, double c=0, double d=0, double e=0, double f=0)
{
  return pow(a*a+b*b+c*c+d*d+e*e+f*f, 0.5);
}

std::string itoa(int i) 
{
  char res[10];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

TCanvas* comparePlots2(RooPlot *plot_bC, RooPlot *plot_bS, TH1F *data, TH1F *qcd, std::string title)
{
  TCanvas *c=new TCanvas(("c_RooFit_"+title).c_str(), "c", 700, 700);
  TPad *p_2=new TPad("p_2", "p_2", 0, 0, 1, 0.35);
  TPad *p_1=new TPad("p_1", "p_1", 0, 0.35, 1, 1);
  p_1->SetBottomMargin(0.05);
  p_1->SetFillStyle(4000);
  p_1->SetFrameFillColor(0);
  p_2->SetFillStyle(4000);
  p_2->SetFrameFillColor(0);
  p_1->Draw();
  p_2->Draw();
  p_1->cd();
  double maxdata=data->GetMaximum();
  double maxqcd=qcd->GetMaximum();
  double maxy=(maxdata>maxqcd) ? maxdata : maxqcd;
  title="; m_{X} (GeV); Events / ("+itoa(data->GetBinWidth(1))+" GeV)";
  p_1->DrawFrame(VR_lo-100, 0, VR_hi+100, maxy*1., title.c_str());
  plot_bC->Draw("same");
  plot_bS->Draw("same");
  p_2->cd();
  p_2->SetGridy();
  TH1F *h_ratio=(TH1F*)data->Clone("h_ratio");
  h_ratio->SetTitle("; m_{X} (GeV); VR/VB Ratio");
  h_ratio->Divide(qcd);
  h_ratio->SetLineColor(1);
  h_ratio->SetMarkerStyle(20);
  h_ratio->GetXaxis()->SetRangeUser(VR_lo-100, VR_hi+100);                         
  h_ratio->SetMinimum(h_ratio->GetMinimum()*1.2); h_ratio->SetMaximum(h_ratio->GetMaximum()*1.1);                  
  h_ratio->Draw();
  TLine *m_one_line=new TLine(VR_lo, 1, VR_hi, 1);
  m_one_line->Draw("same");
  p_1->cd(); 
  return c;                          
} 

void BackgroundPrediction_Kinematic_GaussExp_Chebychev()
{

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  
  if (bReg) tags=tags+"_bReg";
  
  const unsigned int nPoints=4;
  double mass[nPoints]={90., 107.5, 142.5, 160.};
  double n_SB[nPoints], n_SR[nPoints];
  double ratio[nPoints];
  double errorsY[nPoints], errorsX[nPoints];
  double ratioAt125=-1, errorAt125=-1;
  
  // TFile *f_MMMM_a=new TFile((tags+"/a/Histograms_8TeVData2012BCD_Skim.root").c_str());
  // TFile *f_MMMM_a=new TFile((tags+"/a/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TFile *f_MMMM_a=new TFile((tags+"/a_KinFit/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  // TFile *f_MMMM_a=new TFile((tags+"/a_KinFit/Histograms_BJetPlusX_Run2012BCD_Skim_selected_bTagged_.root").c_str());
  TH1F *h_mX_CR2_a=(TH1F*)f_MMMM_a->Get("h_mX_CR2");
  TH1F *h_mX_CR4_a=(TH1F*)f_MMMM_a->Get("h_mX_CR4");
  TH1F *h_mX_SR_a=(TH1F*)f_MMMM_a->Get("h_mX_SR");
  ratioAt125=h_mX_SR_a->GetSumOfWeights()/(h_mX_CR2_a->GetSumOfWeights()+h_mX_CR4_a->GetSumOfWeights());
  
  // === MMMM/b ===
  // TFile *f_MMMM_b=new TFile((tags+"/b/Histograms_8TeVData2012BCD_Skim.root").c_str());
  // TFile *f_MMMM_b=new TFile((tags+"/b/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  TFile *f_MMMM_b=new TFile((tags+"/b_KinFit/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  // TFile *f_MMMM_b=new TFile((tags+"/b_KinFit/Histograms_BJetPlusX_Run2012BCD_Skim_selected_bTagged_.root").c_str());
  TH1F *h_mX_CR2_b=(TH1F*)f_MMMM_b->Get("h_mX_CR2");
  TH1F *h_mX_CR4_b=(TH1F*)f_MMMM_b->Get("h_mX_CR4");
  TH1F *h_mX_SR_b=(TH1F*)f_MMMM_b->Get("h_mX_SR");
  n_SB[0]=(h_mX_CR2_b->GetSumOfWeights()+h_mX_CR4_b->GetSumOfWeights());
  n_SR[0]=h_mX_SR_b->GetSumOfWeights();
  
  if (useRatioFit)
  {
    // === MMMM/d ===
    // TFile *f_MMMM_d=new TFile((tags+"/d/Histograms_8TeVData2012BCD_Skim.root").c_str());
    TFile *f_MMMM_d=new TFile((tags+"/d/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
    // TFile *f_MMMM_d=new TFile((tags+"/d/Histograms_BJetPlusX_Run2012BCD_Skim_selected_bTagged_.root").c_str());
    TH1F *h_mX_CR2_d=(TH1F*)f_MMMM_d->Get("h_mX_CR2");
    TH1F *h_mX_CR4_d=(TH1F*)f_MMMM_d->Get("h_mX_CR4");
    TH1F *h_mX_SR_d=(TH1F*)f_MMMM_d->Get("h_mX_SR");
    n_SB[1]=(h_mX_CR2_d->GetSumOfWeights()+h_mX_CR4_d->GetSumOfWeights());
    n_SR[1]=h_mX_SR_d->GetSumOfWeights();
  
    // === MMMM/e ===
    // TFile *f_MMMM_e=new TFile((tags+"/e/Histograms_8TeVData2012BCD_Skim.root").c_str());
    TFile *f_MMMM_e=new TFile((tags+"/e/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
    // TFile *f_MMMM_e=new TFile((tags+"/e/Histograms_BJetPlusX_Run2012BCD_Skim_selected_bTagged_.root").c_str());
    TH1F *h_mX_CR2_e=(TH1F*)f_MMMM_e->Get("h_mX_CR2");
    TH1F *h_mX_CR4_e=(TH1F*)f_MMMM_e->Get("h_mX_CR4");
    TH1F *h_mX_SR_e=(TH1F*)f_MMMM_e->Get("h_mX_SR");
    n_SB[2]=(h_mX_CR2_e->GetSumOfWeights()+h_mX_CR4_e->GetSumOfWeights());
    n_SR[2]=h_mX_SR_e->GetSumOfWeights();
  
    // === MMMM/c ===
    // TFile *f_MMMM_c=new TFile((tags+"/c/Histograms_8TeVData2012BCD_Skim.root").c_str());
    TFile *f_MMMM_c=new TFile((tags+"/c/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
    // TFile *f_MMMM_c=new TFile((tags+"/c/Histograms_BJetPlusX_Run2012BCD_Skim_selected_bTagged_.root").c_str());
    TH1F *h_mX_CR2_c=(TH1F*)f_MMMM_c->Get("h_mX_CR2");
    TH1F *h_mX_CR4_c=(TH1F*)f_MMMM_c->Get("h_mX_CR4");
    TH1F *h_mX_SR_c=(TH1F*)f_MMMM_c->Get("h_mX_SR");
    n_SB[3]=(h_mX_CR2_c->GetSumOfWeights()+h_mX_CR4_c->GetSumOfWeights());
    n_SR[3]=h_mX_SR_c->GetSumOfWeights();
  
    for (unsigned int i=0; i<nPoints; ++i)
    {
      ratio[i]=n_SR[i]/n_SB[i];
      errorsY[i]=ratio[i]*pow(1./n_SR[i]+1./n_SB[i], 0.5);
      errorsX[i]=0.;
    }
  
    TGraphErrors *g_ratio=new TGraphErrors(nPoints, mass, ratio, errorsX, errorsY);
    g_ratio->SetTitle("SR/SB ratio");
    TCanvas *c_ratio=new TCanvas("c_ratio", "c_ratio", 700, 700);
    g_ratio->SetMinimum(0.); g_ratio->SetMaximum(2.);
    g_ratio->Draw("A*");
    TF1 *f_ratio=new TF1("f_ratio", "pol1");
    g_ratio->Fit(f_ratio);
    c_ratio->SaveAs(("c_ratio_"+tags+".png").c_str());
  
    ratioAt125=f_ratio->Eval(125.);
    errorAt125=(errorsY[0]+errorsY[1]+errorsY[2]+errorsY[3])/4.;
  }
  
  std::cout<<"ratioAt125 = "<<ratioAt125<<" +- "<<errorAt125<<std::endl;
  std::cout<<"bgFloat   lnN     -   "<<1.+errorAt125/ratioAt125<<std::endl;
  
  std::cout<<" = MMMM b ======================================== "<<std::endl;
  TH1F *h_mMMMMb_3Tag_CR2=(TH1F*)f_MMMM_b->Get("h_mX_CR2");
  TH1F *h_mMMMMb_3Tag_CR4=(TH1F*)f_MMMM_b->Get("h_mX_CR4");
  TH1F *h_mMMMMb_3Tag_SR=(TH1F*)f_MMMM_b->Get("h_mX_SR");
  double bS=h_mMMMMb_3Tag_SR->GetSumOfWeights();
  std::cout<<"Number of events in MMMM b signal region = "<<bS<<std::endl;
  TH1F *h_mMMMMb_3Tag_CR24=(TH1F*)h_mMMMMb_3Tag_CR2->Clone("h_mX_CR24");
  h_mMMMMb_3Tag_CR24->Add(h_mMMMMb_3Tag_CR4);
  h_mMMMMb_3Tag_CR24->Rebin(rebin);
  h_mMMMMb_3Tag_SR->Rebin(rebin);
  h_mMMMMb_3Tag_CR24->SetLineColor(kRed);
  h_mMMMMb_3Tag_SR->SetLineColor(kBlue);
  double bC=h_mMMMMb_3Tag_CR24->GetSumOfWeights();
  std::cout<<"bC = "<<bC<<", bS = "<<bS<<std::endl;
  h_mMMMMb_3Tag_CR24->SetMaximum(h_mMMMMb_3Tag_CR24->GetMaximum()*1.3);
  h_mMMMMb_3Tag_CR24->SetTitle(("Kinematic Extrapolation in "+tags+" Validation Region; m_{X} GeV").c_str());
  h_mMMMMb_3Tag_SR->Scale(bC/bS);
  // Do the fits using RooFit
  // gSystem->Load("../PDFs/GaussExp_cxx.so");
  gSystem->Load("../PDFs/WrappedChebychev3_cxx.so");
  RooRealVar x("x", "m_{X} (GeV)", VR_lo-100., VR_hi+100.);
  // bC
  RooRealVar bC_p0("bC_p0", "bC_p0", 250.);
  RooRealVar bC_p1("bC_p1", "bC_p1", 800.);
  RooRealVar bC_p2("bC_p2", "bC_p2", -1., 1.);
  RooRealVar bC_p3("bC_p3", "bC_p3", -1., 1.);
  RooRealVar bC_p4("bC_p4", "bC_p4", -1., 1.);
  RooRealVar bC_p5("bC_p5", "bC_p5", -1., 1.);
  WrappedChebychev3 bC_fit("bC_fit", "bC Fit", x, bC_p0, bC_p1, bC_p2, bC_p3, bC_p4, bC_p5);
  h_mMMMMb_3Tag_CR24->GetXaxis()->SetRangeUser(VR_lo-100., VR_hi+100.);
  RooDataHist bC_data("bC_data", "bC Data", RooArgList(x), h_mMMMMb_3Tag_CR24);
  RooFitResult *r_bC_fit=bC_fit.fitTo(bC_data, RooFit::Range(VR_lo, VR_hi), RooFit::Save());
  RooPlot *bC_plot=x.frame();
  bC_data.plotOn(bC_plot);
  bC_fit.plotOn(bC_plot, RooFit::VisualizeError(*r_bC_fit, 1), RooFit::FillColor(kOrange));
  bC_fit.plotOn(bC_plot, RooFit::LineColor(kRed));
  bC_data.plotOn(bC_plot, RooFit::LineColor(kRed), RooFit::MarkerColor(kRed));
  TCanvas *c_bC=new TCanvas("c_bC", "c_bC", 700, 700);
  bC_plot->Draw();
  c_bC->SaveAs(("c_compareData_"+tags+"_VR_RooFit_Chebychev3.png").c_str());
  
  return;
  
  // bS
  RooRealVar bS_p0("bS_p0", "bS_p0", 300., 500.);
  RooRealVar bS_p1("bS_p1", "bS_p1", 40., 100.1);
  RooRealVar bS_p2("bS_p2", "bS_p2", 0.1, 10.1);
  GaussExp bS_fit("bS_fit", "bS Fit", x, bS_p0, bS_p1, bS_p2);
  h_mMMMMb_3Tag_SR->GetXaxis()->SetRangeUser(VR_lo-100., VR_hi+100.);
  RooDataHist bS_data("bS_data", "bS Data", RooArgList(x), h_mMMMMb_3Tag_SR);
  RooFitResult *r_bS_fit=bS_fit.fitTo(bS_data, RooFit::Range(VR_lo, VR_hi), RooFit::Save()); // RooFit::SumW2Error(kTRUE), 
  RooPlot *bS_plot=x.frame();
  bS_data.plotOn(bS_plot);
  bS_fit.plotOn(bS_plot,RooFit::VisualizeError(*r_bS_fit), RooFit::FillColor(kCyan), RooFit::FillStyle(3001));
  bS_fit.plotOn(bS_plot, RooFit::LineColor(kBlue));
  bS_data.plotOn(bS_plot, RooFit::LineColor(kBlue), RooFit::MarkerColor(kBlue));
  std::cout<<" === === "<<std::endl;
  // std::cout<<"chi^2/ndof of bC = "<<bC_plot->chiSquare()<<std::endl;
  // std::cout<<"chi^2/ndof of bS = "<<bS_plot->chiSquare()<<std::endl;
  std::cout<<" === === "<<std::endl;
  TCanvas *c_bC=comparePlots2(bC_plot, bS_plot, h_mMMMMb_3Tag_SR, h_mMMMMb_3Tag_CR24, "Kinematic Extrapolation in "+tags+" Validation Region of Data; m_{X} GeV");
  double x_mean_bC=bC_p0.getVal();
  double x_k_bC=bC_p0.getVal()+bC_p2.getVal()*bC_p1.getVal();
  TLine *l_mean_bC=new TLine(x_mean_bC, 0, x_mean_bC, h_mMMMMb_3Tag_CR24->GetMaximum()*0.8); l_mean_bC->SetLineColor(kRed); l_mean_bC->Draw();
  TLine *l_k_bC=new TLine(x_k_bC, 0, x_k_bC, h_mMMMMb_3Tag_CR24->GetMaximum()*0.8); l_k_bC->SetLineColor(kRed); l_k_bC->SetLineStyle(9); l_k_bC->Draw();
  double x_mean_bS=bS_p0.getVal();
  double x_k_bS=bS_p0.getVal()+bS_p2.getVal()*bS_p1.getVal();
  TLine *l_mean_bS=new TLine(x_mean_bS, 0, x_mean_bS, h_mMMMMb_3Tag_SR->GetMaximum()); l_mean_bS->SetLineColor(kBlue); l_mean_bS->Draw();
  TLine *l_k_bS=new TLine(x_k_bS, 0, x_k_bS, h_mMMMMb_3Tag_SR->GetMaximum()); l_k_bS->SetLineColor(kBlue); l_k_bS->SetLineStyle(9); l_k_bS->Draw();
  c_bC->SaveAs(("c_compareData_"+tags+"_VR_RooFit_Chebyshev3.png").c_str());
  
  // Calculate Pi and DPi and dPi -- for shape systematics
  double PbC_0=bC_p0.getVal();
  double PbC_1=bC_p1.getVal();
  double PbC_2=bC_p2.getVal();
  double dPbC_0=bC_p0.getError();
  double dPbC_1=bC_p1.getError();
  double dPbC_2=bC_p2.getError();
  double PbS_0=bS_p0.getVal();
  double PbS_1=bS_p1.getVal();
  double PbS_2=bS_p2.getVal();
  double dPbS_0=bS_p0.getError();
  double dPbS_1=bS_p1.getError();
  double dPbS_2=bS_p2.getError();
  
  std::cout<<" = MMMM Background Prediction ==== "<<std::endl;
  TH1F *h_mMMMMa_3Tag_CR2=(TH1F*)f_MMMM_a->Get("h_mX_CR2");
  TH1F *h_mMMMMa_3Tag_CR4=(TH1F*)f_MMMM_a->Get("h_mX_CR4");
  TH1F *h_mMMMMa_3Tag_SR;
  if (tags!="MMMM") h_mMMMMa_3Tag_SR=(TH1F*)f_MMMM_a->Get("h_mX_SR");
  TH1F *h_mMMMMa_3Tag_CR24=(TH1F*)h_mMMMMa_3Tag_CR2->Clone("h_mX_CR24");
  h_mMMMMa_3Tag_CR24->Add(h_mMMMMa_3Tag_CR4);
  h_mMMMMa_3Tag_CR24->Rebin(rebin);
  h_mMMMMa_3Tag_CR24->SetLineColor(kBlack);
  if (tags!="MMMM") h_mMMMMa_3Tag_SR->Rebin(rebin);
  if (tags!="MMMM") h_mMMMMa_3Tag_SR->SetLineColor(kBlue);
  TH1F *h_mMMMMa_3Tag_SR_Prediction=(TH1F*)h_mMMMMa_3Tag_CR24->Clone("h_mMMMMa_3Tag_SR_Prediction");
  // double aC=h_mMMMMa_3Tag_CR24->GetSumOfWeights();
  // Get the scale of the prediction right
  std::cout<<"bS/bC = "<<bS/bC<<std::endl;
  std::cout<<"(bC) bgFloat   lnN     -    "<<1.+sqrt(1./bC)<<std::endl;
  std::cout<<"(bC + bS) bgFloat   lnN     -    "<<1.+sqrt(1./bC+1./bS)<<std::endl;
  // std::cout<<"ratioAt125 = "<<ratioAt125<<", +- "<<errorAt125<<" (fract unc.) = "<<1.+errorAt125/ratioAt125<<std::endl;
  // h_mMMMMa_3Tag_SR_Prediction->Scale(ratioAt125);
  std::cout<<"Number of predicted events in 17.928 /fb = "<<h_mMMMMa_3Tag_SR_Prediction->GetSumOfWeights()*ratioAt125<<std::endl;
  std::cout<<"Number of predicted events in 17.928 /fb around mX=509(+-26) GeV = "<<(h_mMMMMa_3Tag_SR_Prediction->Integral(h_mMMMMa_3Tag_SR_Prediction->FindBin(509.-26.), h_mMMMMa_3Tag_SR_Prediction->FindBin(509.+26.)))*249./h_mMMMMa_3Tag_SR_Prediction->GetSumOfWeights()<<std::endl;
  std::cout<<"Number of predicted events in 17.928 /fb around mX=714(+-40) GeV = "<<(h_mMMMMa_3Tag_SR_Prediction->Integral(h_mMMMMa_3Tag_SR_Prediction->FindBin(714.-40.), h_mMMMMa_3Tag_SR_Prediction->FindBin(714.+40.)))*249./h_mMMMMa_3Tag_SR_Prediction->GetSumOfWeights()<<std::endl;
  std::cout<<"Number of predicted events in 17.928 /fb around mX=450-550 GeV = "<<(h_mMMMMa_3Tag_SR_Prediction->Integral(h_mMMMMa_3Tag_SR_Prediction->FindBin(450.), h_mMMMMa_3Tag_SR_Prediction->FindBin(550.)))*249./h_mMMMMa_3Tag_SR_Prediction->GetSumOfWeights()<<std::endl;
  std::cout<<"Number of predicted events in 17.928 /fb around mX=600-800 GeV = "<<(h_mMMMMa_3Tag_SR_Prediction->Integral(h_mMMMMa_3Tag_SR_Prediction->FindBin(600.), h_mMMMMa_3Tag_SR_Prediction->FindBin(800.)))*249./h_mMMMMa_3Tag_SR_Prediction->GetSumOfWeights()<<std::endl;
  // RooFit fit to background prediction
  // RooRealVar bg_p0("bg_p0", "bg_p0", 400., 600.);
  // RooRealVar bg_p1("bg_p1", "bg_p1", 50., 100.1);
  // RooRealVar bg_p2("bg_p2", "bg_p2", 0.1, 10.1);
  // For mX300
  RooRealVar bg_p0("bg_p0", "bg_p0", 350., 600.);
  RooRealVar bg_p1("bg_p1", "bg_p1", 30., 100.1);
  RooRealVar bg_p2("bg_p2", "bg_p2", 0.01, 10.1);
  GaussExp bg("bg", "Background Prediction PDF", x, bg_p0, bg_p1, bg_p2);
  RooDataHist pred("pred", "Prediction from SB", RooArgList(x), h_mMMMMa_3Tag_SR_Prediction);
  RooFitResult *r_bg=bg.fitTo(pred, RooFit::Range(SR_lo, SR_hi), RooFit::Save());
  RooPlot *aC_plot=x.frame();
  pred.plotOn(aC_plot, RooFit::MarkerColor(kRed));
  bg.plotOn(aC_plot, RooFit::VisualizeError(*r_bg, 1), RooFit::FillColor(kOrange));
  bg.plotOn(aC_plot, RooFit::LineColor(kRed));
  pred.plotOn(aC_plot, RooFit::LineColor(kRed), RooFit::MarkerColor(kRed));
  
  TCanvas *c_rooFit=new TCanvas("c_rooFit", "c_rooFit", 700, 700);
  TPad *p_2=new TPad("p_2", "p_2", 0, 0, 1, 0.35);
  TPad *p_1=new TPad("p_1", "p_1", 0, 0.35, 1, 1);
  p_1->SetBottomMargin(0.05);
  p_1->SetFillStyle(4000);
  p_1->SetFrameFillColor(0);
  p_2->SetFillStyle(4000);
  p_2->SetFrameFillColor(0);
  p_1->Draw();
  p_2->Draw();
  p_1->cd();
  if (tags!="MMMM") h_mMMMMa_3Tag_SR->Draw("Ep9 SAME");
  aC_plot->SetTitle("; m_{X} (GeV); Events / (10 GeV)");
  aC_plot->Draw();
  double x_mean_aC=bg_p0.getVal();
  double x_k_aC=bg_p0.getVal()+bg_p2.getVal()*bg_p1.getVal();
  TLine *l_mean_aC=new TLine(x_mean_aC, 0, x_mean_aC, h_mMMMMa_3Tag_SR_Prediction->GetMaximum()); l_mean_aC->SetLineColor(kRed); l_mean_aC->Draw();
  TLine *l_k_aC=new TLine(x_k_aC, 0, x_k_aC, h_mMMMMa_3Tag_SR_Prediction->GetMaximum()); l_k_aC->SetLineColor(kRed); l_k_aC->SetLineStyle(9); l_k_aC->Draw();
  
  // Prediction Curve with Shape Systematics
  double PaC_0=bg_p0.getVal();
  double PaC_1=bg_p1.getVal();
  double PaC_2=bg_p2.getVal();
  double dPaC_0=bg_p0.getError();
  double dPaC_1=bg_p1.getError();
  double dPaC_2=bg_p2.getError();
  double PaS_0=PaC_0*PbS_0/PbC_0;
  double PaS_1=PaC_1*PbS_1/PbC_1;
  double PaS_2=PaC_2*PbS_2/PbC_2;
  double dPaS_0=PaS_0*quad((dPaC_0/PaC_0), (dPbS_0/PbS_0), (dPbC_0/PbC_0));
  double dPaS_1=PaS_1*quad((dPaC_1/PaC_1), (dPbS_1/PbS_1), (dPbC_1/PbC_1));
  double dPaS_2=PaS_2*quad((dPaC_2/PaC_2), (dPbS_2/PbS_2), (dPbC_2/PbC_2));
  std::cout<<"(dPaC_0/PaC_0) = ("<<dPaC_0<<"/"<<PaC_0<<") = "<<(dPaC_0/PaC_0)<<"; (dPbS_0/PbS_0) = ("<<dPbS_0<<"/"<<PbS_0<<") = "<<(dPbS_0/PbS_0)<<"; (dPbC_0/PbC_0) = ("<<dPbC_0<<"/"<<PbC_0<<") = "<<(dPbC_0/PbC_0)<<std::endl;
  std::cout<<"(dPaC_1/PaC_1) = ("<<dPaC_1<<"/"<<PaC_1<<") = "<<(dPaC_1/PaC_1)<<"; (dPbS_1/PbS_1) = ("<<dPbS_1<<"/"<<PbS_1<<") = "<<(dPbS_1/PbS_1)<<"; (dPbC_1/PbC_1) = ("<<dPbC_1<<"/"<<PbC_1<<") = "<<(dPbC_1/PbC_1)<<std::endl; 
  std::cout<<"(dPaC_2/PaC_2) = ("<<dPaC_2<<"/"<<PaC_2<<") = "<<(dPaC_2/PaC_2)<<"; (dPbS_2/PbS_2) = ("<<dPbS_2<<"/"<<PbS_2<<") = "<<(dPbS_2/PbS_2)<<"; (dPbC_2/PbC_2) = ("<<dPbC_2<<"/"<<PbC_2<<") = "<<(dPbC_2/PbC_2)<<std::endl; 
  std::cout<<" Predicted PaS_0 = "<<PaS_0<<" +- "<<dPaS_0<<std::endl;
  std::cout<<" Predicted PaS_1 = "<<PaS_1<<" +- "<<dPaS_1<<std::endl;
  std::cout<<" Predicted PaS_2 = "<<PaS_2<<" +- "<<dPaS_2<<std::endl;
  RooRealVar *bg_pred0;
  RooRealVar *bg_pred1;
  RooRealVar *bg_pred2;
  // Let the mean and std dev float a bit
  // dPaS_0=dPaS_0*5.;
  // dPaS_1=dPaS_1*5.;
  // dPaS_2=dPaS_2*5.;
  std::cout<<"Relaxed ranges:"<<std::endl;
  std::cout<<" Predicted PaS_0 = "<<PaS_0<<" +- "<<dPaS_0<<std::endl;
  std::cout<<" Predicted PaS_1 = "<<PaS_1<<" +- "<<dPaS_1<<std::endl;
  std::cout<<" Predicted PaS_2 = "<<PaS_2<<" +- "<<dPaS_2<<std::endl;
  if (tags!="MMMM_nominal")
  { 
    bg_pred0=(new RooRealVar("bg_pred0", "bg_pred0", PaS_0-dPaS_0/2., PaS_0+dPaS_0/2.));
    bg_pred1=(new RooRealVar("bg_pred1", "bg_pred1", PaS_1-dPaS_1/2., PaS_1+dPaS_1/2.));
    bg_pred2=(new RooRealVar("bg_pred2", "bg_pred2", PaS_2-dPaS_2/2., PaS_2+dPaS_2/2.));
  }
  else
  {
    bg_pred0=(new RooRealVar("bg_pred0", "bg_pred0", PaS_0));  bg_pred0->setError(dPaS_0);
    bg_pred1=(new RooRealVar("bg_pred1", "bg_pred1", PaS_1));  bg_pred1->setError(dPaS_1);
    bg_pred2=(new RooRealVar("bg_pred2", "bg_pred2", PaS_2));  bg_pred2->setError(dPaS_2);
  }
  GaussExp bg_pred_init("background_init", "Background Predicted for Signal Region", x, *bg_pred0, *bg_pred1, *bg_pred2);
  GaussExp bg_pred("background", "Background Predicted for Signal Region", x, *bg_pred0, *bg_pred1, *bg_pred2);
  RooPlot *aS_plot=x.frame();
  if (tags!="MMMM_nominal")
  {
    RooDataHist unblind("unblind", "Signal Region", RooArgList(x), h_mMMMMa_3Tag_SR);
    unblind.plotOn(aS_plot, RooFit::LineColor(kBlue), RooFit::MarkerColor(kBlue));
    // bg_pred_init.plotOn(aS_plot, RooFit::LineColor(kGreen), RooFit::Range(SR_lo, SR_hi));
    RooFitResult *r_bg_pred=bg_pred.fitTo(unblind, RooFit::Range(SR_lo, SR_hi), RooFit::Save());
    bg_pred.plotOn(aS_plot, RooFit::VisualizeError(*r_bg_pred, 1), RooFit::FillColor(kCyan), RooFit::FillStyle(3001));
    bg_pred.plotOn(aS_plot, RooFit::LineColor(kBlue));
    bg_pred.plotOn(aS_plot, RooFit::Name("r_bg_prediction"));
    unblind.plotOn(aS_plot, RooFit::LineColor(kBlue), RooFit::MarkerColor(kBlue));
    aS_plot->Draw("same");
  }
  else
  {
    bg_pred.plotOn(aC_plot, RooFit::LineColor(kGreen), RooFit::Range(SR_lo, SR_hi));
    aC_plot->Draw("same");
  }
  double x_mean_aS=bg_pred0->getVal();
  double x_k_aS=bg_pred0->getVal()+bg_pred2->getVal()*bg_pred1->getVal();
  TLine *l_mean_aS=new TLine(x_mean_aS, 0, x_mean_aS, h_mMMMMa_3Tag_SR_Prediction->GetMaximum()); l_mean_aS->SetLineColor(kBlue); l_mean_aS->Draw();
  TLine *l_k_aS=new TLine(x_k_aS, 0, x_k_aS, h_mMMMMa_3Tag_SR_Prediction->GetMaximum()); l_k_aS->SetLineColor(kBlue); l_k_aS->SetLineStyle(9); l_k_aS->Draw();
  
  std::cout<<" === === "<<std::endl;
  // std::cout<<"chi^2/ndof of bC = "<<bC_plot->chiSquare()<<std::endl;
  // std::cout<<"chi^2/ndof of bS = "<<bS_plot->chiSquare()<<std::endl;
  // std::cout<<"chi^2/ndof of aC = "<<aC_plot->chiSquare()<<std::endl;
  // std::cout<<"chi^2/ndof of aS = "<<aS_plot->chiSquare()<<std::endl;
  std::cout<<" === === "<<std::endl;
  
  TLatex * tPrel = new TLatex();
  tPrel->SetNDC();
  tPrel->SetTextColor(kBlack);
  tPrel->SetTextSize(0.04);
  tPrel->DrawLatex(0.1,0.95,"CMS Preliminary; #sqrt{s} =  8 TeV, L=17.928 fb^{-1}");
  
  p_2->cd();
  p_2->SetGridy();
  RooHist *hpull;
  hpull=aC_plot->pullHist();
  RooPlot *frameP=x->frame();
  frameP->SetTitle("; m_{X} (GeV); (data-fit)/fit");
  frameP->addPlotable(hpull, "P");
  frameP->Draw();
  TLine *m_one_line=new TLine(SR_lo, 0, SR_hi, 0); m_one_line->Draw();
   
  c_rooFit->SaveAs(("c_compareData_"+tags+"_SR_RooFit_GaussExp.png").c_str());
  
  RooWorkspace *w=new RooWorkspace("HbbHbb");
  w->import(bg_pred);
  w->SaveAs("w_background_GaussExp.root");
  
}
  
  
  
