// Creates the images and HTML 
// for displaying changes in Signal MC
// due to JEC+1-1, and JER+1-1

#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <iostream>
#include <TFractionFitter.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <THStack.h>
#include <TArrow.h>
#include <TColor.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooRealVar.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooArgList.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooChebychev.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooDataHist.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooAbsPdf.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooWorkspace.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooPlot.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooFitResult.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooCBShape.h"
#include "/Applications/FWLite/CMSSW_6_1_1_FWLITE/osx108_amd64_gcc472/lcg/roofit/5.34.04-cms2/include/RooGaussian.h"

// #include "/Users/souvik/HbbHbb/Analysis/PDFs/ExpGaussExp.h"

bool focus=false;
bool writeDataCard=false;

int rebin=1;
ofstream outfile;

std::string tostr(float t)
{
  std::ostringstream os; 
  os<<t; 
  return os.str(); 
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

TLegend* twoStatBoxes(TH1F* h1, TH1F* h2)
{
  TLegend *leg=new TLegend(0.6, 0.6, 0.9, 0.9);
  double wherepeak=(h2->GetMean())/(h2->GetXaxis()->GetXmax());
  if (wherepeak>0.63) leg=new TLegend(0.1, 0.9, 0.4, 0.6);
  leg->AddEntry(h1, "Baseline");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h1->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h1->GetRMS())).c_str(), "");
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry(h2, "JEC +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h2->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h2->GetRMS())).c_str(), "");
  leg->SetFillColor(10);
  return leg;
}

TLegend* threeStatBoxes(TH1F* h1, TH1F* h2, TH1F *h3)
{
  TLegend *leg=new TLegend(0.6, 0.6, 0.9, 0.9);
  double wherepeak=(h2->GetMean())/(h2->GetXaxis()->GetXmax());
  if (wherepeak>0.63) leg=new TLegend(0.1, 0.9, 0.4, 0.6);
  leg->AddEntry(h1, "Baseline");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h1->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h1->GetRMS())).c_str(), "");
  leg->AddEntry(h2, "JEC +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h2->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h2->GetRMS())).c_str(), "");
  leg->AddEntry(h3, "JEC -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h3->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h3->GetRMS())).c_str(), "");
  leg->SetFillColor(10);
  return leg;
}

TLegend* fiveStatBoxes(TH1F* h1, TH1F* h2, TH1F *h3, TH1F *h4, TH1F *h5)
{
  TLegend *leg=new TLegend(0.6, 0.4, 0.9, 0.9);
  double wherepeak=(h2->GetMean())/(h2->GetXaxis()->GetXmax());
  if (wherepeak>0.63) leg=new TLegend(0.1, 0.9, 0.4, 0.6);
  leg->AddEntry(h1, "Baseline");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h1->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h1->GetRMS())).c_str(), "");
  leg->AddEntry(h2, "JEC +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h2->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h2->GetRMS())).c_str(), "");
  leg->AddEntry(h3, "JEC -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h3->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h3->GetRMS())).c_str(), "");
  leg->AddEntry(h4, "JER +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h4->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h4->GetRMS())).c_str(), "");
  leg->AddEntry(h5, "JER -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h5->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h5->GetRMS())).c_str(), "");
  leg->SetFillColor(10);
  return leg;
}

TLegend* nineStatBoxes(TH1F* h1, TH1F* h2, TH1F *h3, TH1F *h4, TH1F *h5, TH1F *h6, TH1F *h7, TH1F *h8, TH1F *h9)
{
  TLegend *leg=new TLegend(0.7, 0.1, 0.9, 0.9);
  double wherepeak=(h2->GetMean())/(h2->GetXaxis()->GetXmax());
  if (wherepeak>0.63) leg=new TLegend(0.1, 0.9, 0.3, 0.1);
  leg->AddEntry(h1, "Baseline");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h1->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h1->GetRMS())).c_str(), "");
  leg->AddEntry(h2, "JEC +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h2->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h2->GetRMS())).c_str(), "");
  leg->AddEntry(h3, "JEC -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h3->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h3->GetRMS())).c_str(), "");
  leg->AddEntry(h4, "JER +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h4->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h4->GetRMS())).c_str(), "");
  leg->AddEntry(h5, "JER -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h5->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h5->GetRMS())).c_str(), "");
  leg->AddEntry(h6, "Trig SF (CSV) +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h6->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h6->GetRMS())).c_str(), "");
  leg->AddEntry(h7, "Trig SF (CSV) -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h7->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h7->GetRMS())).c_str(), "");
  leg->AddEntry(h8, "Trig SF (p_{T}) +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h8->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h8->GetRMS())).c_str(), "");
  leg->AddEntry(h9, "Trig SF (p_{T}) -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h9->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h9->GetRMS())).c_str(), "");
  leg->SetFillColor(10);
  return leg;
}

TLegend* elevenStatBoxes(TH1F* h1, TH1F* h2, TH1F *h3, TH1F *h4, TH1F *h5, TH1F *h6, TH1F *h7, TH1F *h8, TH1F *h9, TH1F *h10, TH1F *h11)
{
  TLegend *leg=new TLegend(0.7, 0.1, 0.9, 0.9);
  double wherepeak=(h2->GetMean())/(h2->GetXaxis()->GetXmax());
  if (wherepeak>0.63) leg=new TLegend(0.1, 0.9, 0.3, 0.1);
  leg->AddEntry(h1, "Baseline");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h1->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h1->GetRMS())).c_str(), "");
  leg->AddEntry(h2, "JEC +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h2->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h2->GetRMS())).c_str(), "");
  leg->AddEntry(h3, "JEC -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h3->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h3->GetRMS())).c_str(), "");
  leg->AddEntry(h4, "JER +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h4->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h4->GetRMS())).c_str(), "");
  leg->AddEntry(h5, "JER -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h5->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h5->GetRMS())).c_str(), "");
  leg->AddEntry(h6, "SF_{bc} +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h6->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h6->GetRMS())).c_str(), "");
  leg->AddEntry(h7, "SF_{bc} -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h7->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h7->GetRMS())).c_str(), "");
  leg->AddEntry(h8, "SF_{l} +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h8->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h8->GetRMS())).c_str(), "");
  leg->AddEntry(h9, "SF_{l} -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h9->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h9->GetRMS())).c_str(), "");
  leg->AddEntry(h10, "Trig SF +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h10->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h10->GetRMS())).c_str(), "");
  leg->AddEntry(h11, "Trig SF -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h11->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h11->GetRMS())).c_str(), "");
  leg->SetFillColor(10);
  return leg;
}

TLegend* thirteenStatBoxes(TH1F* h1, TH1F* h2, TH1F *h3, TH1F *h4, TH1F *h5, TH1F *h6, TH1F *h7, TH1F *h8, TH1F *h9, TH1F *h10, TH1F *h11, TH1F *h12, TH1F *h13)
{
  TLegend *leg=new TLegend(0.7, 0.1, 0.9, 0.9);
  double wherepeak=(h2->GetMean())/(h2->GetXaxis()->GetXmax());
  if (wherepeak>0.63) leg=new TLegend(0.1, 0.9, 0.3, 0.1);
  leg->AddEntry(h1, "Baseline");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h1->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h1->GetRMS())).c_str(), "");
  leg->AddEntry(h2, "JEC +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h2->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h2->GetRMS())).c_str(), "");
  leg->AddEntry(h3, "JEC -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h3->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h3->GetRMS())).c_str(), "");
  leg->AddEntry(h4, "JER +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h4->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h4->GetRMS())).c_str(), "");
  leg->AddEntry(h5, "JER -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h5->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h5->GetRMS())).c_str(), "");
  leg->AddEntry(h6, "SF_{bc} +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h6->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h6->GetRMS())).c_str(), "");
  leg->AddEntry(h7, "SF_{bc} -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h7->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h7->GetRMS())).c_str(), "");
  leg->AddEntry(h8, "SF_{l} +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h8->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h8->GetRMS())).c_str(), "");
  leg->AddEntry(h9, "SF_{l} -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h9->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h9->GetRMS())).c_str(), "");
  leg->AddEntry(h10, "Trig SF (CSV) +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h10->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h10->GetRMS())).c_str(), "");
  leg->AddEntry(h11, "Trig SF (CSV) -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h11->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h11->GetRMS())).c_str(), "");
  leg->AddEntry(h12, "Trig SF (p_{T}) +1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h12->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h12->GetRMS())).c_str(), "");
  leg->AddEntry(h13, "Trig SF (p_{T}) -1 #sigma");
  leg->AddEntry((TObject*)0, ("mean="+tostr(h13->GetMean())).c_str(), "");
  leg->AddEntry((TObject*)0, ("rms="+tostr(h13->GetRMS())).c_str(), "");
  leg->SetFillColor(10);
  return leg;
}

Double_t crystalBall(Double_t *x, Double_t *par)
{
  Double_t std=(x[0]-par[2])/par[3];
  Double_t A=pow(par[1]/par[0], par[1])*exp(-0.5*pow(par[0], 2));
  Double_t B=par[1]/par[0]-par[0];
  Double_t result=0.;
  
  if (std<par[0]) // Gaussian region
  {
    result=exp(-0.5*pow(std, 2));
  }
  else // Power Law region
  {
    result=A/pow(B+std, par[1]);
  }
  
  result=result*par[4];
  
  return result;
}

THStack* drawCombinatorics(std::string mass, TH1F *h_right, TH1F *h_wrong, TH1F *h_no)
{
  std::string s_name=h_right->GetName();
  s_name.replace(0, 1, "s");
  
  double rangeLo=-1, rangeHi=-1;
  if (mass=="300") {rangeLo=200., rangeHi=600.;}
  else if (mass=="400") {rangeLo=200., rangeHi=600.;}
  else if (mass=="500") {rangeLo=470., rangeHi=650.;}
  else if (mass=="600") {rangeLo=580., rangeHi=670.;}
  else if (mass=="700") {rangeLo=650., rangeHi=870.;}
  else if (mass=="800") {rangeLo=750., rangeHi=990.;}
  
  h_right->GetXaxis()->SetRangeUser(rangeLo, rangeHi);
  h_wrong->GetXaxis()->SetRangeUser(rangeLo, rangeHi);
  h_no->GetXaxis()->SetRangeUser(rangeLo, rangeHi);
  
  THStack *stack=new THStack(s_name.c_str(), s_name.c_str());
  stack->Add(h_no);
  stack->Add(h_wrong);
  stack->Add(h_right);
  
  return stack;
}

RooPlot* fitSignal(TH1F *h, std::string mass, int color, TLegend *leg, Params &params, bool kinFit=false)
{
  
  RooRealVar *x, *sg_p0, *sg_p1, *sg_p2, *sg_p3;
  x=new RooRealVar("x", "m_{X} (GeV)", 300., 1100.);
  // x=new RooRealVar("x", "m_{X} (GeV)", 300., 800.);
  double rangeLo=-1, rangeHi=-1;
  if (!kinFit)
  {
    if (mass=="270")
    {
      rangeLo=200., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 280., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 1.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 1.);
    }
    else if (mass=="300")
    {
      rangeLo=200., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 280., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 1.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 1.);
    }
    else if (mass=="350")
    {
      rangeLo=200., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 300., 370.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 1.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 1.);
    }
    else if (mass=="400")
    {
      rangeLo=300., rangeHi=550.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 380., 420.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 15.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="450")
    {
      rangeLo=370., rangeHi=540.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 430., 470.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 25.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="500")
    {
      rangeLo=420., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 460., 520.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="550")
    {
      rangeLo=450., rangeHi=650.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 510., 570.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="600")
    {
      rangeLo=490., rangeHi=700.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 560., 620.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="650")
    {
      rangeLo=540., rangeHi=750.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 610., 670.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 20., 35.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="700")
    {
      rangeLo=580., rangeHi=800.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 660., 720.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 20., 35.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="800")
    {
      rangeLo=650., rangeHi=900.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 760., 820.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 25., 40.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="900")
    {
      rangeLo=720., rangeHi=1000.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 870., 930.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 25., 40.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="1000")
    {
      rangeLo=800., rangeHi=1150.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 950., 1050.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 25., 40.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="1100")
    {
      rangeLo=850., rangeHi=1250.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 1080., 1150.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 30., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 5.);
    }
  }
  else
  {
    if (mass=="270")
    {
      rangeLo=200., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 280., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 1.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 1.);
    }
    else if (mass=="300")
    {
      rangeLo=200., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 280., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 1.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 1.);
    }
    else if (mass=="350")
    {
      rangeLo=200., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 300., 370.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 1.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 1.);
    }
    else if (mass=="400")
    {
      rangeLo=360., rangeHi=490.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 380., 420.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 15.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="450")
    {
      rangeLo=400., rangeHi=540.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 430., 480.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 20.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="500")
    {
      rangeLo=450., rangeHi=600.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 480., 530.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 25.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="550")
    {
      rangeLo=490., rangeHi=650.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 530., 580.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 25.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="600")
    {
      rangeLo=540., rangeHi=700.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 580., 630.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 25.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="650")
    {
      rangeLo=580., rangeHi=750.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 630., 680.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 25.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="700")
    {
      rangeLo=620., rangeHi=810.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 690., 740.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 15., 25.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="800")
    {
      rangeLo=720., rangeHi=920.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 780., 850.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 15., 30.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="900")
    {
      rangeLo=790., rangeHi=1030.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 870., 950.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 20., 35.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="1000")
    {
      rangeLo=850., rangeHi=1150.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 980., 1050.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 25., 40.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="1100")
    {
      rangeLo=950., rangeHi=1250.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 1080., 1170.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 25., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
  }
  x=new RooRealVar("x", "m_{X} (GeV)", rangeLo-50., rangeHi+50.);
  // RevCrystalBall signal("signal", "Signal Prediction", *x, *sg_p0, *sg_p1, *sg_p2, *sg_p3);
  ExpGaussExp signal("signal", "Signal Prediction", *x, *sg_p0, *sg_p1, *sg_p2, *sg_p3);
  RooDataHist signalHistogram("signalHistogram", "Signal Histogram", RooArgList(*x), h);
  signal.fitTo(signalHistogram, RooFit::Range(rangeLo, rangeHi), RooFit::Save());
  params.sg_p0=sg_p0->getVal(); params.sg_p0_err=sg_p0->getError();
  params.sg_p1=sg_p1->getVal(); params.sg_p1_err=sg_p1->getError();
  params.sg_p2=sg_p2->getVal(); params.sg_p2_err=sg_p2->getError();
  params.sg_p3=sg_p3->getVal(); params.sg_p3_err=sg_p3->getError();
  RooPlot *plot=x->frame();
  if (color==kBlack)
  {
    signalHistogram.plotOn(plot, RooFit::MarkerColor(color), RooFit::MarkerSize(1.2));
    signal.plotOn(plot, RooFit::LineColor(color), RooFit::LineWidth(3));
  }
  else 
  {
    signalHistogram.plotOn(plot, RooFit::MarkerColor(color));
    signal.plotOn(plot, RooFit::LineColor(color), RooFit::LineWidth(0));
  }
  // leg->AddEntry((TObject*)0, ("norm="+tostr(h->GetSumOfWeights())).c_str(), "");
  // leg->AddEntry((TObject*)0, ("mean="+tostr(sg_p2->getVal())+"#pm"+tostr(sg_p2->getError())).c_str(), "");
  // leg->AddEntry((TObject*)0, ("rms="+tostr(sg_p3->getVal())+"#pm"+tostr(sg_p3->getError())).c_str(), "");
  // std::cout<<"chi2/dof = "<<plot->chiSquare()<<std::endl;
  
  // Save modified signal shape to workspace
  if (color==kBlack)
  {
    RooRealVar signal_p0("signal_p0", "signal_p0", sg_p0->getVal());
    RooRealVar signal_p1("signal_p1", "signal_p1", sg_p1->getVal());
    RooRealVar signal_p2("signal_p2", "signal_p2", sg_p2->getVal());
    RooRealVar signal_p3("signal_p3", "signal_p3", sg_p3->getVal());
    // RevCrystalBall signal_fixed("signal", "Signal Prediction Fixed", *x, signal_p0, signal_p1, signal_p2, signal_p3);
    ExpGaussExp signal_fixed("signal", "Signal Prediction Fixed", *x, signal_p0, signal_p1, signal_p2, signal_p3);
    RooWorkspace *w=new RooWorkspace("HbbHbb");
    w->import(signal_fixed);
    if (!kinFit) w->SaveAs(("SignalShapes/w_signal_"+mass+".root").c_str());
    if (kinFit) w->SaveAs(("SignalShapes_KinFit/w_signal_"+mass+".root").c_str());
  }
  return plot;
}

RooPlot* fitSignal_Gaussian(TH1F *h, std::string mass, int color, TLegend *leg, Params &params, bool kinFit=false)
{ 
  RooRealVar *x, *sg_p0, *sg_p1;
  double rangeLo=-1, rangeHi=-1;
  if (!kinFit)
  {
    if (mass=="270")
    {
      rangeLo=170., rangeHi=470.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 250., 290.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 60.);
    }
    else if (mass=="300")
    {
      rangeLo=200., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 280., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 60.);
    }
    else if (mass=="350")
    {
      rangeLo=250., rangeHi=650.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 330., 370.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 60.);
    }
    else if (mass=="400")
    {
      rangeLo=300., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 380., 420.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
    }
    else if (mass=="450")
    {
      rangeLo=350., rangeHi=650.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 430., 470.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
    }
    else if (mass=="500")
    {
      rangeLo=300., rangeHi=700.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 470., 530.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
    }
    else if (mass=="550")
    {
      rangeLo=350., rangeHi=750.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 530., 580.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
    }
    else if (mass=="600")
    {
      rangeLo=400., rangeHi=800.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 580., 620.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 20., 30.);
    }
    else if (mass=="650")
    {
      rangeLo=450., rangeHi=850.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 630., 680.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 20., 30.);
    }
    else if (mass=="700")
    {
      rangeLo=500., rangeHi=900.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 680., 720.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 25., 35.);
    }
    else if (mass=="800")
    {
      rangeLo=600., rangeHi=1000.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 780., 850.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 25., 40.);
    }
    else if (mass=="900")
    {
      rangeLo=700., rangeHi=1100.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 880., 950.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 25., 40.);
    }
    else if (mass=="1000")
    {
      rangeLo=800., rangeHi=1200.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 980., 1050.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 25., 40.);
    }
    else if (mass=="1100")
    {
      rangeLo=900., rangeHi=1300.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 1080., 1150.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 25., 40.);
    }
  }
  else
  {
    if (mass=="270")
    {
      rangeLo=200., rangeHi=400.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 250., 290.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 60.);
    }
    else if (mass=="300")
    {
      rangeLo=250., rangeHi=350.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 280., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 7., 60.);
    }
    else if (mass=="350")
    {
      rangeLo=300., rangeHi=400.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 330., 370.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 60.);
    }
    else if (mass=="400")
    {
      // rangeLo=300., rangeHi=600.;
      rangeLo=380., rangeHi=440.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 390., 410.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 15.);
    }
    else if (mass=="450")
    {
      rangeLo=430., rangeHi=480.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 430., 470.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
    }
    else if (mass=="500")
    {
      // rangeLo=300., rangeHi=700.;
      rangeLo=475., rangeHi=540.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 470., 530.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
    }
    else if (mass=="550")
    {
      rangeLo=520., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 530., 580.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
    }
    else if (mass=="600")
    {
      // rangeLo=400., rangeHi=800.;
      rangeLo=575., rangeHi=650.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 600., 620.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
    }
    else if (mass=="650")
    {
      rangeLo=630., rangeHi=700.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 630., 680.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
    }
    else if (mass=="700")
    {
      // rangeLo=500., rangeHi=900.;
      rangeLo=675., rangeHi=750.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 680., 750.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 15., 30.);
    }
    else if (mass=="800")
    {
      // rangeLo=600., rangeHi=1000.; 
      rangeLo=770., rangeHi=860.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 780., 850.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 15., 30.);
    }
    else if (mass=="900")
    {
      rangeLo=860., rangeHi=970.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 880., 950.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 20., 40.);
    }
    else if (mass=="1000")
    {
      rangeLo=950., rangeHi=1080.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 980., 1050.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 25., 40.);
    }
    else if (mass=="1100")
    {
      rangeLo=1040., rangeHi=1180.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 1080., 1150.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 30., 50.);
    }
  }
  x=new RooRealVar("x", "m_{X} (GeV)", rangeLo-50, rangeHi+50);
  RooGaussian signal("signal", "Signal Prediction", *x, *sg_p0, *sg_p1);
  RooDataHist signalHistogram("signalHistogram", "Signal Histogram", RooArgList(*x), h);
  signal.fitTo(signalHistogram, RooFit::Range(rangeLo, rangeHi), RooFit::Save());
  params.sg_p0=sg_p0->getVal(); params.sg_p0_err=sg_p0->getError();
  params.sg_p1=sg_p1->getVal(); params.sg_p1_err=sg_p1->getError();
  params.sg_p2=-1;
  params.sg_p3=-1;
  RooPlot *plot=x->frame();
  if (color==kBlack)
  {
    signalHistogram.plotOn(plot, RooFit::MarkerColor(color), RooFit::MarkerSize(1.2));
    signal.plotOn(plot, RooFit::LineColor(color), RooFit::LineWidth(3));
  }
  else 
  {
    signalHistogram.plotOn(plot, RooFit::MarkerColor(color));
    signal.plotOn(plot, RooFit::LineColor(color), RooFit::LineWidth(0));
  }
  // leg->AddEntry((TObject*)0, ("norm="+tostr(h->GetSumOfWeights())).c_str(), "");
  // leg->AddEntry((TObject*)0, ("mean="+tostr(sg_p0->getVal())+"#pm"+tostr(sg_p0->getError())).c_str(), "");
  // leg->AddEntry((TObject*)0, ("rms="+tostr(sg_p1->getVal())+"#pm"+tostr(sg_p1->getError())).c_str(), "");
  // std::cout<<"chi2/dof = "<<plot->chiSquare()<<std::endl;
  
  // Save modified signal shape to workspace
  if (color==kBlack)
  {
    RooRealVar signal_p0("signal_p0", "signal_p0", sg_p0->getVal());
    RooRealVar signal_p1("signal_p1", "signal_p1", sg_p1->getVal());
    RooGaussian signal_fixed("signal", "Signal Prediction Fixed", *x, signal_p0, signal_p1);
    RooWorkspace *w=new RooWorkspace("HbbHbb");
    w->import(signal_fixed);
    if (!kinFit) w->SaveAs(("SignalShapes/w_signal_Gaussian_"+mass+".root").c_str());
    if (kinFit) w->SaveAs(("SignalShapes_KinFit/w_signal_Gaussian_"+mass+".root").c_str());
  }
  return plot;
}

double lnN(double b, double a, double c)
{
  // std::cout<<"a = "<<a<<", b = "<<b<<", c = "<<c<<std::endl;
  // std::cout<<"1.+(a-c)/(2.*b) = "<<1.+fabs(a-c)/(2.*b)<<std::endl;
  double err=0;
  if (b>0) err=1.+fabs(a-c)/(2.*b);
  return err;
}

int DisplayJECJERDistributions_FULLSIM()
{

  std::cout<<"hello"<<std::endl;
  
  std::vector<std::string> masses;
  masses.push_back("270");
  masses.push_back("300");
  masses.push_back("350");
  masses.push_back("400");
  masses.push_back("450");
  masses.push_back("500");
  masses.push_back("550");
  masses.push_back("600");
  masses.push_back("650");
  masses.push_back("700");
  masses.push_back("800");
  masses.push_back("900");
  masses.push_back("1000");
  masses.push_back("1100");
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  // gStyle->SetPalette(1);
  // gSystem->Load("/Users/souvik/HbbHbb/Analysis/PDFs/RevCrystalBall_cxx.so");
  gSystem->Load("/Users/souvik/HbbHbb/Analysis/PDFs/ExpGaussExp_cxx.so");
  
  // Calculate nSignal events given production cross section, branching fractions and efficiency
  // double totalLumi=18600.0; // /pb
  double totalLumi=17928.0; // /pb
  double prodXsec_1=1.; // pb
  
  // Interpolation Plots
  std::vector<double> v_sg_p0, v_sg_p0_err;
  std::vector<double> v_sg_p1, v_sg_p1_err;
  std::vector<double> v_sg_p2, v_sg_p2_err;
  std::vector<double> v_sg_p3, v_sg_p3_err;
  std::vector<double> v_zero;
  
  // Write to an HTML File
  outfile.open("SignalSystematics/SignalSystematics.html");
  outfile<<"<html>"<<std::endl;
  outfile<<"<head>"<<std::endl;
  // outfile<<"<base href=\"https://cmslpcweb.fnal.gov/uscms_data/souvik/SignalSystematics\" target=\"_blank\">"<<std::endl;
  outfile<<"<script type=\"text\/javascript\">"<<std::endl;
  outfile<<"function toggleMe(a){"<<std::endl;
  outfile<<"var e=document.getElementById(a);"<<std::endl;
  outfile<<"if(!e)return true;"<<std::endl;
  outfile<<"if(e.style.display==\"none\"){"<<std::endl;
  outfile<<"e.style.display=\"block\""<<std::endl;
  outfile<<"}"<<std::endl;
  outfile<<"else{"<<std::endl;
  outfile<<"e.style.display=\"none\""<<std::endl;
  outfile<<"}"<<std::endl;
  outfile<<"return true;"<<std::endl;
  outfile<<"}"<<std::endl;
  outfile<<"</script>"<<std::endl;
  outfile<<"</head>"<<std::endl;
  outfile<<"<body>"<<std::endl;
  outfile<<"<h1> Signal with JEC, JER Systematic Uncertainties</h1>"<<std::endl;
  
  for (unsigned int i=0; i<masses.size(); ++i)
  // for (unsigned int i=10; i<11; ++i)
  {
    v_zero.push_back(0);
    
    TFile *file=new TFile(("MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT=(TH1F*)file->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT=(TH1F*)file->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged=(TH1F*)file->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged->Add((TH1F*)file->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged=(TH1F*)file->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged->Add((TH1F*)file->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR=(TH1F*)file->Get("h_mX_SR");
    TH1F *h_mX_SR_rightComb=(TH1F*)file->Get("h_mX_SR_rightComb");
    TH1F *h_mX_SR_wrongComb=(TH1F*)file->Get("h_mX_SR_wrongComb");
    TH1F *h_mX_SR_noComb=(TH1F*)file->Get("h_mX_SR_noComb");
    TH1F *h_CountWithPU=(TH1F*)file->Get("CountWithPU");
    double nSignal_init=h_CountWithPU->GetBinContent(1);
    std::cout<<"nSignal_init = "<<nSignal_init<<std::endl;
    
    TFile *file_JECp1;
    if (focus) file_JECp1=new TFile(("MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_JECp1=new TFile(("KinSel_JECp1/MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_JECp1=(TH1F*)file_JECp1->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_JECp1=(TH1F*)file_JECp1->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_JECp1=(TH1F*)file_JECp1->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_JECp1->Add((TH1F*)file_JECp1->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_JECp1=(TH1F*)file_JECp1->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_JECp1->Add((TH1F*)file_JECp1->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_JECp1=(TH1F*)file_JECp1->Get("h_mX_SR");
    
    TFile *file_JECm1;
    if (focus) file_JECm1=new TFile(("MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_JECm1=new TFile(("KinSel_JECm1/MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_JECm1=(TH1F*)file_JECm1->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_JECm1=(TH1F*)file_JECm1->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_JECm1=(TH1F*)file_JECm1->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_JECm1->Add((TH1F*)file_JECm1->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_JECm1=(TH1F*)file_JECm1->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_JECm1->Add((TH1F*)file_JECm1->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_JECm1=(TH1F*)file_JECm1->Get("h_mX_SR");
    
    TFile *file_JERp1;
    if (focus) file_JERp1=new TFile(("MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_JERp1=new TFile(("KinSel_JERp1/MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_JERp1=(TH1F*)file_JERp1->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_JERp1=(TH1F*)file_JERp1->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_JERp1=(TH1F*)file_JERp1->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_JERp1->Add((TH1F*)file_JERp1->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_JERp1=(TH1F*)file_JERp1->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_JERp1->Add((TH1F*)file_JERp1->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_JERp1=(TH1F*)file_JERp1->Get("h_mX_SR");
    
    TFile *file_JERm1;
    if (focus) file_JERm1=new TFile(("MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_JERm1=new TFile(("KinSel_JERm1/MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_JERm1=(TH1F*)file_JERm1->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_JERm1=(TH1F*)file_JERm1->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_JERm1=(TH1F*)file_JERm1->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_JERm1->Add((TH1F*)file_JERm1->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_JERm1=(TH1F*)file_JERm1->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_JERm1->Add((TH1F*)file_JERm1->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_JERm1=(TH1F*)file_JERm1->Get("h_mX_SR");
    
    TFile *file_upTrigCSV;
    if (focus) file_upTrigCSV=new TFile(("MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_upTrigCSV=new TFile(("MMMM_upTrigCSV/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_upTrigCSV=(TH1F*)file_upTrigCSV->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_upTrigCSV=(TH1F*)file_upTrigCSV->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_upTrigCSV=(TH1F*)file_upTrigCSV->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_upTrigCSV->Add((TH1F*)file_upTrigCSV->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_upTrigCSV=(TH1F*)file_upTrigCSV->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_upTrigCSV->Add((TH1F*)file_upTrigCSV->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_upTrigCSV=(TH1F*)file_upTrigCSV->Get("h_mX_SR");
    
    TFile *file_downTrigCSV;
    if (focus) file_downTrigCSV=new TFile(("MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_downTrigCSV=new TFile(("MMMM_downTrigCSV/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_downTrigCSV=(TH1F*)file_downTrigCSV->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_downTrigCSV=(TH1F*)file_downTrigCSV->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_downTrigCSV=(TH1F*)file_downTrigCSV->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_downTrigCSV->Add((TH1F*)file_downTrigCSV->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_downTrigCSV=(TH1F*)file_downTrigCSV->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_downTrigCSV->Add((TH1F*)file_downTrigCSV->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_downTrigCSV=(TH1F*)file_downTrigCSV->Get("h_mX_SR");
    
    TFile *file_upTrigpT;
    if (focus) file_upTrigpT=new TFile(("MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_upTrigpT=new TFile(("MMMM_upTrigpT/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_upTrigpT=(TH1F*)file_upTrigpT->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_upTrigpT=(TH1F*)file_upTrigpT->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_upTrigpT=(TH1F*)file_upTrigpT->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_upTrigpT->Add((TH1F*)file_upTrigpT->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_upTrigpT=(TH1F*)file_upTrigpT->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_upTrigpT->Add((TH1F*)file_upTrigpT->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_upTrigpT=(TH1F*)file_upTrigpT->Get("h_mX_SR");
    
    TFile *file_downTrigpT;
    if (focus) file_downTrigpT=new TFile(("MMMM_nominal/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_downTrigpT=new TFile(("MMMM_downTrigpT/a/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_downTrigpT=(TH1F*)file_downTrigpT->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_downTrigpT=(TH1F*)file_downTrigpT->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_downTrigpT=(TH1F*)file_downTrigpT->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_downTrigpT->Add((TH1F*)file_downTrigpT->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_downTrigpT=(TH1F*)file_downTrigpT->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_downTrigpT->Add((TH1F*)file_downTrigpT->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_downTrigpT=(TH1F*)file_downTrigpT->Get("h_mX_SR");
    
    // Kinematically Constrained
    TFile *file_KinFit=new TFile(("MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_KinFit=(TH1F*)file_KinFit->Get("h_mX_SR");
    TH1F *h_mX_SR_rightComb_KinFit=(TH1F*)file_KinFit->Get("h_mX_SR_rightComb");
    TH1F *h_mX_SR_wrongComb_KinFit=(TH1F*)file_KinFit->Get("h_mX_SR_wrongComb");
    TH1F *h_mX_SR_noComb_KinFit=(TH1F*)file_KinFit->Get("h_mX_SR_noComb");
    
    TFile *file_JECp1_KinFit;
    if (focus) file_JECp1_KinFit=new TFile(("MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_JECp1_KinFit=new TFile(("KinSel_JECp1/MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_JECp1_KinFit=(TH1F*)file_JECp1_KinFit->Get("h_mX_SR");
   
    TFile *file_JECm1_KinFit;
    if (focus) file_JECm1_KinFit=new TFile(("MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_JECm1_KinFit=new TFile(("KinSel_JECm1/MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_JECm1_KinFit=(TH1F*)file_JECm1_KinFit->Get("h_mX_SR");
    
    TFile *file_JERp1_KinFit;
    if (focus) file_JERp1_KinFit=new TFile(("MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_JERp1_KinFit=new TFile(("KinSel_JERp1/MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_JERp1_KinFit=(TH1F*)file_JERp1_KinFit->Get("h_mX_SR");
   
    TFile *file_JERm1_KinFit;
    if (focus) file_JERm1_KinFit=new TFile(("MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_JERm1_KinFit=new TFile(("KinSel_JERm1/MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_JERm1_KinFit=(TH1F*)file_JERm1_KinFit->Get("h_mX_SR");
    
    TFile *file_upTrigCSV_KinFit;
    if (focus) file_upTrigCSV_KinFit=new TFile(("MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_upTrigCSV_KinFit=new TFile(("MMMM_upTrigCSV/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_upTrigCSV_KinFit=(TH1F*)file_upTrigCSV_KinFit->Get("h_mX_SR");
    
    TFile *file_downTrigCSV_KinFit;
    if (focus) file_downTrigCSV_KinFit=new TFile(("MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_downTrigCSV_KinFit=new TFile(("MMMM_downTrigCSV/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_downTrigCSV_KinFit=(TH1F*)file_downTrigCSV_KinFit->Get("h_mX_SR");
    
    TFile *file_upTrigpT_KinFit;
    if (focus) file_upTrigpT_KinFit=new TFile(("MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_upTrigpT_KinFit=new TFile(("MMMM_upTrigpT/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_upTrigpT_KinFit=(TH1F*)file_upTrigpT_KinFit->Get("h_mX_SR");
    
    TFile *file_downTrigpT_KinFit;
    if (focus) file_downTrigpT_KinFit=new TFile(("MMMM_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    else file_downTrigpT_KinFit=new TFile(("MMMM_downTrigpT/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_downTrigpT_KinFit=(TH1F*)file_downTrigpT_KinFit->Get("h_mX_SR");
    
    TCanvas *c_H1Jet1pT=new TCanvas("c_H1Jet1pT", "c_H1Jet1pT", 700, 700);
    h_H1Jet1pT->SetLineWidth(2);
    h_H1Jet1pT_JECp1->SetLineStyle(9); h_H1Jet1pT_JECp1->SetLineColor(kRed);
    h_H1Jet1pT_JECm1->SetLineStyle(9); h_H1Jet1pT_JECm1->SetLineColor(kRed+2);
    h_H1Jet1pT_JERp1->SetLineStyle(9); h_H1Jet1pT_JERp1->SetLineColor(kMagenta);
    h_H1Jet1pT_JERm1->SetLineStyle(9); h_H1Jet1pT_JERm1->SetLineColor(kMagenta+2);
    h_H1Jet1pT_upTrigCSV->SetLineStyle(9); h_H1Jet1pT_upTrigCSV->SetLineColor(kBlue);
    h_H1Jet1pT_downTrigCSV->SetLineStyle(9); h_H1Jet1pT_downTrigCSV->SetLineColor(kBlue+2);
    h_H1Jet1pT_upTrigpT->SetLineStyle(9); h_H1Jet1pT_upTrigpT->SetLineColor(kBlue-2);
    h_H1Jet1pT_downTrigpT->SetLineStyle(9); h_H1Jet1pT_downTrigpT->SetLineColor(kBlue+3);
    h_H1Jet1pT_upTrigCSV->Draw();
    h_H1Jet1pT->Draw("same");
    h_H1Jet1pT_JECp1->Draw("same");
    h_H1Jet1pT_JECm1->Draw("same");
    h_H1Jet1pT_JERp1->Draw("same");
    h_H1Jet1pT_JERm1->Draw("same");
    h_H1Jet1pT_upTrigCSV->Draw("same");
    h_H1Jet1pT_downTrigCSV->Draw("same");
    h_H1Jet1pT_upTrigpT->Draw("same");
    h_H1Jet1pT_downTrigpT->Draw("same");
    nineStatBoxes(h_H1Jet1pT, h_H1Jet1pT_JECp1, h_H1Jet1pT_JECm1, h_H1Jet1pT_JERp1, h_H1Jet1pT_JERm1, h_H1Jet1pT_upTrigCSV, h_H1Jet1pT_downTrigCSV, h_H1Jet1pT_upTrigpT, h_H1Jet1pT_downTrigpT)->Draw();
    c_H1Jet1pT->SaveAs(("SignalSystematics/c_H1Jet1pT_"+masses.at(i)+".png").c_str());
    
    TCanvas *c_H1Jet2pT=new TCanvas("c_H1Jet2pT", "c_H1Jet2pT", 700, 700);
    h_H1Jet2pT->SetLineWidth(2);
    h_H1Jet2pT_JECp1->SetLineStyle(9); h_H1Jet2pT_JECp1->SetLineColor(kRed);
    h_H1Jet2pT_JECm1->SetLineStyle(9); h_H1Jet2pT_JECm1->SetLineColor(kRed+2);
    h_H1Jet2pT_JERp1->SetLineStyle(9); h_H1Jet2pT_JERp1->SetLineColor(kMagenta);
    h_H1Jet2pT_JERm1->SetLineStyle(9); h_H1Jet2pT_JERm1->SetLineColor(kMagenta+2);
    h_H1Jet2pT_upTrigCSV->SetLineStyle(9); h_H1Jet2pT_upTrigCSV->SetLineColor(kBlue);
    h_H1Jet2pT_downTrigCSV->SetLineStyle(9); h_H1Jet2pT_downTrigCSV->SetLineColor(kBlue+2);
    h_H1Jet2pT_upTrigpT->SetLineStyle(9); h_H1Jet2pT_upTrigpT->SetLineColor(kBlue);
    h_H1Jet2pT_downTrigpT->SetLineStyle(9); h_H1Jet2pT_downTrigpT->SetLineColor(kBlue+2);
    h_H1Jet2pT_upTrigCSV->Draw();
    h_H1Jet2pT->Draw("same");
    h_H1Jet2pT_JECp1->Draw("same");
    h_H1Jet2pT_JECm1->Draw("same");
    h_H1Jet2pT_JERp1->Draw("same");
    h_H1Jet2pT_JERm1->Draw("same");
    h_H1Jet2pT_upTrigCSV->Draw("same");
    h_H1Jet2pT_downTrigCSV->Draw("same");
    h_H1Jet2pT_upTrigpT->Draw("same");
    h_H1Jet2pT_downTrigpT->Draw("same");
    nineStatBoxes(h_H1Jet2pT, h_H1Jet2pT_JECp1, h_H1Jet2pT_JECm1, h_H1Jet2pT_JERp1, h_H1Jet2pT_JERm1, h_H1Jet2pT_upTrigCSV, h_H1Jet2pT_downTrigCSV, h_H1Jet2pT_upTrigpT, h_H1Jet2pT_downTrigpT)->Draw();
    c_H1Jet2pT->SaveAs(("SignalSystematics/c_H1Jet2pT_"+masses.at(i)+".png").c_str());
    
    TCanvas *c_H1_mass_bTagged=new TCanvas("c_H1_mass_bTagged", "c_H1_mass_bTagged", 700, 700);
    h_H1_mass_bTagged->SetLineWidth(2);
    h_H1_mass_bTagged_JECp1->SetLineStyle(9); h_H1_mass_bTagged_JECp1->SetLineColor(kRed);
    h_H1_mass_bTagged_JECm1->SetLineStyle(9); h_H1_mass_bTagged_JECm1->SetLineColor(kRed+2);
    h_H1_mass_bTagged_JERp1->SetLineStyle(9); h_H1_mass_bTagged_JERp1->SetLineColor(kMagenta);
    h_H1_mass_bTagged_JERm1->SetLineStyle(9); h_H1_mass_bTagged_JERm1->SetLineColor(kMagenta+2);
    h_H1_mass_bTagged_upTrigCSV->SetLineStyle(9); h_H1_mass_bTagged_upTrigCSV->SetLineColor(kBlue);
    h_H1_mass_bTagged_downTrigCSV->SetLineStyle(9); h_H1_mass_bTagged_downTrigCSV->SetLineColor(kBlue+2);
    h_H1_mass_bTagged_upTrigpT->SetLineStyle(9); h_H1_mass_bTagged_upTrigpT->SetLineColor(kBlue);
    h_H1_mass_bTagged_downTrigpT->SetLineStyle(9); h_H1_mass_bTagged_downTrigpT->SetLineColor(kBlue+2);
    h_H1_mass_bTagged_upTrigCSV->Draw();
    h_H1_mass_bTagged->Draw("same");
    h_H1_mass_bTagged_JECp1->Draw("same");
    h_H1_mass_bTagged_JECm1->Draw("same");
    h_H1_mass_bTagged_JERp1->Draw("same");
    h_H1_mass_bTagged_JERm1->Draw("same");
    h_H1_mass_bTagged_upTrigCSV->Draw("same");
    h_H1_mass_bTagged_downTrigCSV->Draw("same");
    h_H1_mass_bTagged_upTrigpT->Draw("same");
    h_H1_mass_bTagged_downTrigpT->Draw("same");
    nineStatBoxes(h_H1_mass_bTagged, h_H1_mass_bTagged_JECp1, h_H1_mass_bTagged_JECm1, h_H1_mass_bTagged_JERp1, h_H1_mass_bTagged_JERm1, h_H1_mass_bTagged_upTrigCSV, h_H1_mass_bTagged_downTrigCSV, h_H1_mass_bTagged_upTrigpT, h_H1_mass_bTagged_downTrigpT)->Draw();
    c_H1_mass_bTagged->SaveAs(("SignalSystematics/c_H1_mass_bTagged_"+masses.at(i)+".png").c_str());
    
    TCanvas *c_H1_pT_bTagged=new TCanvas("c_H1_pT_bTagged", "c_H1_pT_bTagged", 700, 700);
    h_H1_pT_bTagged->SetLineWidth(2);
    h_H1_pT_bTagged_JECp1->SetLineStyle(9); h_H1_pT_bTagged_JECp1->SetLineColor(kRed);
    h_H1_pT_bTagged_JECm1->SetLineStyle(9); h_H1_pT_bTagged_JECm1->SetLineColor(kRed+2);
    h_H1_pT_bTagged_JERp1->SetLineStyle(9); h_H1_pT_bTagged_JERp1->SetLineColor(kMagenta);
    h_H1_pT_bTagged_JERm1->SetLineStyle(9); h_H1_pT_bTagged_JERm1->SetLineColor(kMagenta+2);
    h_H1_pT_bTagged_upTrigCSV->SetLineStyle(9); h_H1_pT_bTagged_upTrigCSV->SetLineColor(kBlue);
    h_H1_pT_bTagged_downTrigCSV->SetLineStyle(9); h_H1_pT_bTagged_downTrigCSV->SetLineColor(kBlue+2);
    h_H1_pT_bTagged_upTrigpT->SetLineStyle(9); h_H1_pT_bTagged_upTrigpT->SetLineColor(kBlue);
    h_H1_pT_bTagged_downTrigpT->SetLineStyle(9); h_H1_pT_bTagged_downTrigpT->SetLineColor(kBlue+2);
    h_H1_pT_bTagged_upTrigCSV->Draw();
    h_H1_pT_bTagged->Draw("same");
    h_H1_pT_bTagged_JECp1->Draw("same");
    h_H1_pT_bTagged_JECm1->Draw("same");
    h_H1_pT_bTagged_JERp1->Draw("same");
    h_H1_pT_bTagged_JERm1->Draw("same");
    h_H1_pT_bTagged_upTrigCSV->Draw("same");
    h_H1_pT_bTagged_downTrigCSV->Draw("same");
    h_H1_pT_bTagged_upTrigpT->Draw("same");
    h_H1_pT_bTagged_downTrigpT->Draw("same");
    nineStatBoxes(h_H1_pT_bTagged, h_H1_pT_bTagged_JECp1, h_H1_pT_bTagged_JECm1, h_H1_pT_bTagged_JERp1, h_H1_pT_bTagged_JERm1, h_H1_pT_bTagged_upTrigCSV, h_H1_pT_bTagged_downTrigCSV, h_H1_pT_bTagged_upTrigpT, h_H1_pT_bTagged_downTrigpT)->Draw();
    c_H1_pT_bTagged->SaveAs(("SignalSystematics/c_H1_pT_bTagged_"+masses.at(i)+".png").c_str());
    
    TCanvas *c_SignalmX=new TCanvas(("c_SignalmX"+masses.at(i)).c_str(), ("c_SignalmX"+masses.at(i)).c_str(), 700, 700);
    h_mX_SR->SetTitle(("m_{X} Peak in Signal MC (m_{X}="+masses.at(i)+" GeV); m_{X} (GeV)").c_str());
    h_mX_SR->Rebin(rebin);
    h_mX_SR_rightComb->Rebin(rebin);
    h_mX_SR_wrongComb->Rebin(rebin);
    h_mX_SR_noComb->Rebin(rebin);
    h_mX_SR_JECp1->Rebin(rebin);
    h_mX_SR_JECm1->Rebin(rebin);
    h_mX_SR_JERp1->Rebin(rebin);
    h_mX_SR_JERm1->Rebin(rebin);
    h_mX_SR_upTrigCSV->Rebin(rebin);
    h_mX_SR_downTrigCSV->Rebin(rebin);
    h_mX_SR_upTrigpT->Rebin(rebin);
    h_mX_SR_downTrigpT->Rebin(rebin);
    h_mX_SR_rightComb->SetFillColor(kRed);
    h_mX_SR_wrongComb->SetFillColor(kGreen);
    h_mX_SR_noComb->SetFillColor(kBlue);
    h_mX_SR_JECp1->SetLineColor(kRed);
    h_mX_SR_JECm1->SetLineColor(kRed+2);
    h_mX_SR_JERp1->SetLineColor(kMagenta);
    h_mX_SR_JERm1->SetLineColor(kMagenta+2);
    h_mX_SR_upTrigCSV->SetLineColor(kBlue);
    h_mX_SR_downTrigCSV->SetLineColor(kBlue+2);
    h_mX_SR_upTrigpT->SetLineColor(kBlue);
    h_mX_SR_downTrigpT->SetLineColor(kBlue+2);
    h_mX_SR->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_rightComb->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_wrongComb->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_noComb->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JECp1->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JECm1->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JERp1->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JERm1->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_upTrigCSV->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_downTrigCSV->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_upTrigpT->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_downTrigpT->GetXaxis()->SetRangeUser(0, 1200);
    TLegend *leg=new TLegend(0.7, 0.5, 0.9, 0.9);
    if (i>10) leg=new TLegend(0.1, 0.9, 0.3, 0.5);
    leg->AddEntry(h_mX_SR, "Baseline");
    Params par, par_JECp1, par_JECm1, par_JERp1, par_JERm1, par_upTrigCSV, par_downTrigCSV, par_upTrigpT, par_downTrigpT;
    RooPlot *plot=fitSignal(h_mX_SR, masses.at(i), kBlack, leg, par);
    RooPlot *plot_JECp1, *plot_JECm1, *plot_JERp1, *plot_JERm1, *plot_upTrigCSV, *plot_downTrigCSV, *plot_upTrigpT, *plot_downTrigpT;
    if (!focus)
    {
      leg->AddEntry(h_mX_SR_JECp1, "JEC +1 #sigma");
      leg->AddEntry(h_mX_SR_JECm1, "JEC -1 #sigma");
      leg->AddEntry(h_mX_SR_JERp1, "JER +1 #sigma");
      leg->AddEntry(h_mX_SR_JERm1, "JER -1 #sigma");
      leg->AddEntry(h_mX_SR_upTrigCSV, "Trig SF (CSV) +1 #sigma");
      leg->AddEntry(h_mX_SR_downTrigCSV, "Trig SF (CSV) -1 #sigma");
      leg->AddEntry(h_mX_SR_upTrigpT, "Trig SF (p_{T}) +1 #sigma");
      leg->AddEntry(h_mX_SR_downTrigpT, "Trig SF (p_{T}) -1 #sigma");
      
      plot_JECp1=fitSignal(h_mX_SR_JECp1, masses.at(i), kRed, leg, par_JECp1);
      plot_JECm1=fitSignal(h_mX_SR_JECm1, masses.at(i), kRed+2, leg, par_JECm1);
      plot_JERp1=fitSignal(h_mX_SR_JERp1, masses.at(i), kMagenta, leg, par_JERp1);
      plot_JERm1=fitSignal(h_mX_SR_JERm1, masses.at(i), kMagenta+2, leg, par_JERm1);
      plot_upTrigCSV=fitSignal(h_mX_SR_upTrigCSV, masses.at(i), kBlue, leg, par_upTrigCSV);
      plot_downTrigCSV=fitSignal(h_mX_SR_downTrigCSV, masses.at(i), kBlue+2, leg, par_downTrigCSV);
      plot_upTrigpT=fitSignal(h_mX_SR_upTrigpT, masses.at(i), kBlue, leg, par_upTrigpT);
      plot_downTrigpT=fitSignal(h_mX_SR_downTrigpT, masses.at(i), kBlue+2, leg, par_downTrigpT);
    }
    plot->SetMaximum(plot->GetMaximum()*1.1);
    plot->Draw();
    if (!focus)
    {
      plot_JECp1->Draw("same");
      plot_JECm1->Draw("same");
      plot_JERp1->Draw("same");
      plot_JERm1->Draw("same");
      plot_upTrigCSV->Draw("same");
      plot_downTrigCSV->Draw("same");
      plot_upTrigpT->Draw("same");
      plot_downTrigpT->Draw("same");
      leg->SetFillColor(0);
    }
    plot->Draw("same");
    leg->Draw();
    c_SignalmX->SaveAs(("SignalSystematics/c_SignalmX_"+masses.at(i)+".root").c_str());
    c_SignalmX->SaveAs(("SignalSystematics/c_SignalmX_"+masses.at(i)+".png").c_str());
    
    TCanvas *c_SignalmX_KinFit=new TCanvas(("c_SignalmX_KinFit"+masses.at(i)).c_str(), ("c_SignalmX_KinFit"+masses.at(i)).c_str(), 700, 700);
    h_mX_SR_KinFit->SetTitle(("m_{X} Peak in Signal MC (m_{X}="+masses.at(i)+" GeV); m_{X} (GeV)").c_str());
    h_mX_SR_KinFit->Rebin(rebin);
    h_mX_SR_rightComb_KinFit->Rebin(rebin);
    h_mX_SR_wrongComb_KinFit->Rebin(rebin);
    h_mX_SR_noComb_KinFit->Rebin(rebin);
    h_mX_SR_JECp1_KinFit->Rebin(rebin);
    h_mX_SR_JECm1_KinFit->Rebin(rebin);
    h_mX_SR_JERp1_KinFit->Rebin(rebin);
    h_mX_SR_JERm1_KinFit->Rebin(rebin);
    h_mX_SR_upTrigCSV_KinFit->Rebin(rebin);
    h_mX_SR_downTrigCSV_KinFit->Rebin(rebin);
    h_mX_SR_upTrigpT_KinFit->Rebin(rebin);
    h_mX_SR_downTrigpT_KinFit->Rebin(rebin);
    h_mX_SR_rightComb_KinFit->SetFillColor(kRed);
    h_mX_SR_wrongComb_KinFit->SetFillColor(kGreen);
    h_mX_SR_noComb_KinFit->SetFillColor(kBlue);
    h_mX_SR_JECp1_KinFit->SetLineColor(kRed);
    h_mX_SR_JECm1_KinFit->SetLineColor(kRed+2);
    h_mX_SR_JERp1_KinFit->SetLineColor(kMagenta);
    h_mX_SR_JERm1_KinFit->SetLineColor(kMagenta+2);
    h_mX_SR_upTrigCSV_KinFit->SetLineColor(kBlue);
    h_mX_SR_downTrigCSV_KinFit->SetLineColor(kBlue+2);
    h_mX_SR_upTrigpT_KinFit->SetLineColor(kBlue);
    h_mX_SR_downTrigpT_KinFit->SetLineColor(kBlue+2);
    h_mX_SR_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_rightComb_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_wrongComb_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_noComb_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JECp1_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JECm1_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JERp1_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JERm1_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_upTrigCSV_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_downTrigCSV_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_upTrigpT_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_downTrigpT_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    leg=new TLegend(0.7, 0.5, 0.9, 0.9);
    if (i>10) leg=new TLegend(0.1, 0.9, 0.3, 0.5);
    leg->AddEntry(h_mX_SR_KinFit, "Baseline");
    Params par_KinFit, par_JECp1_KinFit, par_JECm1_KinFit, par_JERp1_KinFit, par_JERm1_KinFit, par_upTrigCSV_KinFit, par_downTrigCSV_KinFit, par_upTrigpT_KinFit, par_downTrigpT_KinFit;
    RooPlot *plot_JECp1_KinFit, *plot_JECm1_KinFit, 
            *plot_JERp1_KinFit, *plot_JERm1_KinFit, 
            *plot_upTrigCSV_KinFit, *plot_downTrigCSV_KinFit,
            *plot_upTrigpT_KinFit, *plot_downTrigpT_KinFit;
    RooPlot *plot_KinFit=fitSignal(h_mX_SR_KinFit, masses.at(i), kBlack, leg, par_KinFit, true);
    v_sg_p0.push_back(par_KinFit.sg_p0); v_sg_p0_err.push_back(par_KinFit.sg_p0_err);
    v_sg_p1.push_back(par_KinFit.sg_p1); v_sg_p1_err.push_back(par_KinFit.sg_p1_err);
    v_sg_p2.push_back(par_KinFit.sg_p2); v_sg_p2_err.push_back(par_KinFit.sg_p2_err);
    v_sg_p3.push_back(par_KinFit.sg_p3); v_sg_p3_err.push_back(par_KinFit.sg_p3_err);
    if (!focus) {
      leg->AddEntry(h_mX_SR_JECp1_KinFit, "JEC +1 #sigma");
      leg->AddEntry(h_mX_SR_JECm1_KinFit, "JEC -1 #sigma");
      leg->AddEntry(h_mX_SR_JERp1_KinFit, "JER +1 #sigma");
      leg->AddEntry(h_mX_SR_JERm1_KinFit, "JER -1 #sigma");
      leg->AddEntry(h_mX_SR_upTrigCSV_KinFit, "Trig SF (CSV) +1 #sigma");
      leg->AddEntry(h_mX_SR_downTrigCSV_KinFit, "Trig SF (CSV) -1 #sigma");
      leg->AddEntry(h_mX_SR_upTrigpT_KinFit, "Trig SF (p_{T}) +1 #sigma");
      leg->AddEntry(h_mX_SR_downTrigpT_KinFit, "Trig SF (p_{T}) -1 #sigma");
    
      plot_JECp1_KinFit=fitSignal(h_mX_SR_JECp1_KinFit, masses.at(i), kRed, leg, par_JECp1_KinFit, true);
      plot_JECm1_KinFit=fitSignal(h_mX_SR_JECm1_KinFit, masses.at(i), kRed+2, leg, par_JECm1_KinFit, true);
      plot_JERp1_KinFit=fitSignal(h_mX_SR_JERp1_KinFit, masses.at(i), kMagenta, leg, par_JERp1_KinFit, true);
      plot_JERm1_KinFit=fitSignal(h_mX_SR_JERm1_KinFit, masses.at(i), kMagenta+2, leg, par_JERm1_KinFit, true);
      plot_upTrigCSV_KinFit=fitSignal(h_mX_SR_upTrigCSV_KinFit, masses.at(i), kBlue, leg, par_upTrigCSV_KinFit, true);
      plot_downTrigCSV_KinFit=fitSignal(h_mX_SR_downTrigCSV_KinFit, masses.at(i), kBlue+2, leg, par_downTrigCSV_KinFit, true);
      plot_upTrigpT_KinFit=fitSignal(h_mX_SR_upTrigpT_KinFit, masses.at(i), kBlue, leg, par_upTrigpT_KinFit, true);
      plot_downTrigpT_KinFit=fitSignal(h_mX_SR_downTrigpT_KinFit, masses.at(i), kBlue+2, leg, par_downTrigpT_KinFit, true);
    }
    plot_KinFit->SetMaximum(plot_KinFit->GetMaximum()*1.2);
    plot_KinFit->Draw();
    // drawCombinatorics(masses.at(i), h_mX_SR_rightComb_KinFit, h_mX_SR_wrongComb_KinFit, h_mX_SR_noComb_KinFit)->Draw("same");
    if (!focus) 
    {
      plot_JECp1_KinFit->Draw("same");
      plot_JECm1_KinFit->Draw("same");
      plot_JERp1_KinFit->Draw("same");
      plot_JERm1_KinFit->Draw("same");
      plot_upTrigCSV_KinFit->Draw("same");
      plot_downTrigCSV_KinFit->Draw("same");
      plot_upTrigpT_KinFit->Draw("same");
      plot_downTrigpT_KinFit->Draw("same");
    }
    plot_KinFit->Draw("same");
    leg->SetFillColor(0);
    leg->Draw();
    c_SignalmX_KinFit->SaveAs(("SignalSystematics/c_SignalmX_KinFit_"+masses.at(i)+".root").c_str());
    c_SignalmX_KinFit->SaveAs(("SignalSystematics/c_SignalmX_KinFit_"+masses.at(i)+".png").c_str());
    
    TCanvas *c_SignalmX_Gaussian_KinFit=new TCanvas(("c_SignalmX_Gaussian_KinFit"+masses.at(i)).c_str(), ("c_SignalmX_KinFit"+masses.at(i)).c_str(), 700, 700);
    leg=new TLegend(0.35, 0.4, 0.1, 0.9);
    if (i==0) leg=new TLegend(0.7, 0.9, 0.9, 0.4);
    leg->AddEntry(h_mX_SR_KinFit, "Baseline");
    Params par_Gaussian_KinFit, par_JECp1_Gaussian_KinFit, par_JECm1_Gaussian_KinFit, par_JERp1_Gaussian_KinFit, par_JERm1_Gaussian_KinFit, par_upTrigCSV_Gaussian_KinFit, par_downTrigCSV_Gaussian_KinFit, par_upTrigpT_Gaussian_KinFit, par_downTrigpT_Gaussian_KinFit;
    RooPlot *plot_JECp1_Gaussian_KinFit, *plot_JECm1_Gaussian_KinFit, 
            *plot_JERp1_Gaussian_KinFit, *plot_JERm1_Gaussian_KinFit,
            *plot_upTrigCSV_Gaussian_KinFit, *plot_downTrigCSV_Gaussian_KinFit,
            *plot_upTrigpT_Gaussian_KinFit, *plot_downTrigpT_Gaussian_KinFit;
    RooPlot *plot_Gaussian_KinFit=fitSignal_Gaussian(h_mX_SR_KinFit, masses.at(i), kBlack, leg, par_Gaussian_KinFit, true);
    if (!focus) {
      leg->AddEntry(h_mX_SR_JECp1_KinFit, "JEC +1 #sigma");
      leg->AddEntry(h_mX_SR_JECm1_KinFit, "JEC -1 #sigma");
      leg->AddEntry(h_mX_SR_JERp1_KinFit, "JER +1 #sigma");
      leg->AddEntry(h_mX_SR_JERm1_KinFit, "JER -1 #sigma");
      leg->AddEntry(h_mX_SR_upTrigCSV_KinFit, "Trig SF (CSV) +1 #sigma");
      leg->AddEntry(h_mX_SR_downTrigCSV_KinFit, "Trig SF (CSV) -1 #sigma");
      leg->AddEntry(h_mX_SR_upTrigpT_KinFit, "Trig SF (p_{T}) +1 #sigma");
      leg->AddEntry(h_mX_SR_downTrigpT_KinFit, "Trig SF (p_{T}) -1 #sigma");
      
      plot_JECp1_Gaussian_KinFit=fitSignal_Gaussian(h_mX_SR_JECp1_KinFit, masses.at(i), kRed, leg, par_JECp1_Gaussian_KinFit, true);
      plot_JECm1_Gaussian_KinFit=fitSignal_Gaussian(h_mX_SR_JECm1_KinFit, masses.at(i), kRed+2, leg, par_JECm1_Gaussian_KinFit, true);
      plot_JERp1_Gaussian_KinFit=fitSignal_Gaussian(h_mX_SR_JERp1_KinFit, masses.at(i), kMagenta, leg, par_JERp1_Gaussian_KinFit, true);
      plot_JERm1_Gaussian_KinFit=fitSignal_Gaussian(h_mX_SR_JERm1_KinFit, masses.at(i), kMagenta+2, leg, par_JERm1_Gaussian_KinFit, true);
      plot_upTrigCSV_Gaussian_KinFit=fitSignal_Gaussian(h_mX_SR_upTrigCSV_KinFit, masses.at(i), kBlue, leg, par_upTrigCSV_Gaussian_KinFit, true);
      plot_downTrigCSV_Gaussian_KinFit=fitSignal_Gaussian(h_mX_SR_downTrigCSV_KinFit, masses.at(i), kBlue+2, leg, par_downTrigCSV_Gaussian_KinFit, true);
      plot_upTrigpT_Gaussian_KinFit=fitSignal_Gaussian(h_mX_SR_upTrigpT_KinFit, masses.at(i), kBlue, leg, par_upTrigpT_Gaussian_KinFit, true);
      plot_downTrigpT_Gaussian_KinFit=fitSignal_Gaussian(h_mX_SR_downTrigpT_KinFit, masses.at(i), kBlue+2, leg, par_downTrigpT_Gaussian_KinFit, true);
    }
    plot_Gaussian_KinFit->SetMaximum(plot_Gaussian_KinFit->GetMaximum()*1.2);
    plot_Gaussian_KinFit->Draw();
    // drawCombinatorics(masses.at(i), h_mX_SR_rightComb_KinFit, h_mX_SR_wrongComb_KinFit, h_mX_SR_noComb_KinFit)->Draw("same");
    if (!focus) {
      plot_JECp1_Gaussian_KinFit->Draw("same");
      plot_JECm1_Gaussian_KinFit->Draw("same");
      plot_JERp1_Gaussian_KinFit->Draw("same");
      plot_JERm1_Gaussian_KinFit->Draw("same");
      plot_upTrigCSV_Gaussian_KinFit->Draw("same");
      plot_downTrigCSV_Gaussian_KinFit->Draw("same");
      plot_upTrigpT_Gaussian_KinFit->Draw("same");
      plot_downTrigpT_Gaussian_KinFit->Draw("same");
    }
    leg->SetFillColor(0);
    // leg->Draw();
    c_SignalmX_Gaussian_KinFit->SaveAs(("SignalSystematics/c_SignalmX_Gaussian_KinFit_"+masses.at(i)+".png").c_str());
    
    outfile<<"<br/><hr/>"<<std::endl;
    outfile<<"<h2> mX = "<<masses.at(i)<<" </h2>"<<std::endl;
    outfile<<"<table border='1'>"<<std::endl;
    outfile<<" <tr>"<<std::endl;
    outfile<<"  <td>"<<std::endl;
    outfile<<"   Higgs jet1 pT"<<std::endl;
    outfile<<"   <img src='"<<("c_H1Jet1pT_"+masses.at(i)+".png")<<"'/>"<<std::endl;
    outfile<<"  </td>"<<std::endl;
    outfile<<"  <td>"<<std::endl;
    outfile<<"   Higgs jet2 pT"<<std::endl;
    outfile<<"   <img src='"<<("c_H1Jet2pT_"+masses.at(i)+".png")<<"'/>"<<std::endl;
    outfile<<"  </td>"<<std::endl;
    outfile<<"  <td>"<<std::endl;
    outfile<<"   Higgs mass"<<std::endl;
    outfile<<"   <img src='"<<("c_H1_mass_bTagged_"+masses.at(i)+".png")<<"'/>"<<std::endl;
    outfile<<"  </td>"<<std::endl;
    outfile<<"  <td>"<<std::endl;
    outfile<<"   Higgs pT"<<std::endl;
    outfile<<"   <img src='"<<("c_H1_pT_bTagged_"+masses.at(i)+".png")<<"'/>"<<std::endl;
    outfile<<"  </td>"<<std::endl;
    outfile<<"  <td>"<<std::endl;
    outfile<<"   <img src='"<<("c_SignalmX_"+masses.at(i)+".png")<<"'/><br/>"<<std::endl;
    outfile<<"   <h2 align='center'>Without Kin-Fit. Fitted to an Exp-Gauss-Exp function.</h2><br/>"<<std::endl;
    outfile<<"   === Baseline plot === </br>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR->GetSumOfWeights()*totalLumi*prodXsec_1/nSignal_init<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par.sg_p0<<" +- "<<par.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par.sg_p1<<" +- "<<par.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par.sg_p2<<" +- "<<par.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par.sg_p3<<" +- "<<par.sg_p3_err<<" <br/>"<<std::endl;
    if (!focus) {
    outfile<<"   <center><input type='button' onclick=\"return toggleMe('"<<masses.at(i)<<"')\" value='Systematics'></center><br/>"<<std::endl;
    outfile<<"   <div id=\""<<masses.at(i)<<"\" style=\"display:none\">"<<std::endl;
    outfile<<"   === JEC +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECp1->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JECp1->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECp1.sg_p0<<" +- "<<par_JECp1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECp1.sg_p1<<" +- "<<par_JECp1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JECp1.sg_p2<<" +- "<<par_JECp1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JECp1.sg_p3<<" +- "<<par_JECp1.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JEC -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECm1->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JECm1->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECm1.sg_p0<<" +- "<<par_JECm1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECm1.sg_p1<<" +- "<<par_JECm1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JECm1.sg_p2<<" +- "<<par_JECm1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JECm1.sg_p3<<" +- "<<par_JECm1.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JER +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERp1->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JERp1->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERp1.sg_p0<<" +- "<<par_JERp1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERp1.sg_p1<<" +- "<<par_JERp1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JERp1.sg_p2<<" +- "<<par_JERp1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JERp1.sg_p3<<" +- "<<par_JERp1.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JER -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERm1->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JERm1->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERm1.sg_p0<<" +- "<<par_JERm1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERm1.sg_p1<<" +- "<<par_JERm1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JERm1.sg_p2<<" +- "<<par_JERm1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JERm1.sg_p3<<" +- "<<par_JERm1.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (CSV) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upTrigCSV->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_upTrig->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upTrigCSV.sg_p0<<" +- "<<par_upTrigCSV.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upTrigCSV.sg_p1<<" +- "<<par_upTrigCSV.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_upTrigCSV.sg_p2<<" +- "<<par_upTrigCSV.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_upTrigCSV.sg_p3<<" +- "<<par_upTrigCSV.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (CSV) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downTrigCSV->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_downTrig->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downTrigCSV.sg_p0<<" +- "<<par_downTrigCSV.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downTrigCSV.sg_p1<<" +- "<<par_downTrigCSV.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_downTrigCSV.sg_p2<<" +- "<<par_downTrigCSV.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_downTrigCSV.sg_p3<<" +- "<<par_downTrigCSV.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (pT) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upTrigpT->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_upTrig->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upTrigpT.sg_p0<<" +- "<<par_upTrigpT.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upTrigpT.sg_p1<<" +- "<<par_upTrigpT.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_upTrigpT.sg_p2<<" +- "<<par_upTrigpT.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_upTrigpT.sg_p3<<" +- "<<par_upTrigpT.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (pT) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downTrigpT->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_downTrig->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downTrigpT.sg_p0<<" +- "<<par_downTrigpT.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downTrigpT.sg_p1<<" +- "<<par_downTrigpT.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_downTrigpT.sg_p2<<" +- "<<par_downTrigpT.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_downTrigpT.sg_p3<<" +- "<<par_downTrigpT.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    outfile<<"   JEC lnN = "<<lnN(h_mX_SR->GetSumOfWeights(), h_mX_SR_JECp1->GetSumOfWeights(), h_mX_SR_JECm1->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   JER lnN = "<<lnN(h_mX_SR->GetSumOfWeights(), h_mX_SR_JERp1->GetSumOfWeights(), h_mX_SR_JERm1->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   Trigger CSV SF lnN = "<<lnN(h_mX_SR->GetSumOfWeights(), h_mX_SR_upTrigCSV->GetSumOfWeights(), h_mX_SR_downTrigCSV->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   Trigger pT SF lnN = "<<lnN(h_mX_SR->GetSumOfWeights(), h_mX_SR_upTrigpT->GetSumOfWeights(), h_mX_SR_downTrigpT->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    outfile<<"   </div>"<<std::endl;
    }
    outfile<<"  </td>"<<std::endl;
    outfile<<"  <td>"<<std::endl;
    outfile<<"   <img src='"<<("c_SignalmX_KinFit_"+masses.at(i)+".png")<<"'/><br/>"<<std::endl;
    outfile<<"   <h2 align='center'>With Kin-Fit. Fitted to an Exp-Gauss-Exp function.</h2><br/>"<<std::endl;
    outfile<<"   === Baseline plot === </br>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_KinFit->GetSumOfWeights()*totalLumi*prodXsec_1/nSignal_init<<" <br/>"<<std::endl;
    if (masses.at(i)=="1100") outfile<<"   Window 1000-1200 = "<<(h_mX_SR_KinFit->Integral(h_mX_SR_KinFit->FindBin(1000), h_mX_SR_KinFit->FindBin(1200)))*totalLumi*prodXsec_1/nSignal_init<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_KinFit.sg_p0<<" +- "<<par_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_KinFit.sg_p1<<" +- "<<par_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_KinFit.sg_p2<<" +- "<<par_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_KinFit.sg_p3<<" +- "<<par_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    if (!focus) {
    outfile<<"   <center><input type='button' onclick=\"return toggleMe('"<<masses.at(i)<<"_KinFit')\" value='Systematics'></center><br/>"<<std::endl;
    outfile<<"   <div id=\""<<masses.at(i)<<"_KinFit\" style=\"display:none\">"<<std::endl;
    outfile<<"   === JEC +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECp1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JECp1_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECp1_KinFit.sg_p0<<" +- "<<par_JECp1_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECp1_KinFit.sg_p1<<" +- "<<par_JECp1_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JECp1_KinFit.sg_p2<<" +- "<<par_JECp1_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JECp1_KinFit.sg_p3<<" +- "<<par_JECp1_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JEC -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECm1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JECm1_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECm1_KinFit.sg_p0<<" +- "<<par_JECm1_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECm1_KinFit.sg_p1<<" +- "<<par_JECm1_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JECm1_KinFit.sg_p2<<" +- "<<par_JECm1_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JECm1_KinFit.sg_p3<<" +- "<<par_JECm1_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JER +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERp1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JERp1_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERp1_KinFit.sg_p0<<" +- "<<par_JERp1_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERp1_KinFit.sg_p1<<" +- "<<par_JERp1_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JERp1_KinFit.sg_p2<<" +- "<<par_JERp1_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JERp1_KinFit.sg_p3<<" +- "<<par_JERp1_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JER -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERm1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JERm1_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERm1_KinFit.sg_p0<<" +- "<<par_JERm1_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERm1_KinFit.sg_p1<<" +- "<<par_JERm1_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JERm1_KinFit.sg_p2<<" +- "<<par_JERm1_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JERm1_KinFit.sg_p3<<" +- "<<par_JERm1_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (CSV) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upTrigCSV_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_upTrig_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upTrigCSV_KinFit.sg_p0<<" +- "<<par_upTrigCSV_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upTrigCSV_KinFit.sg_p1<<" +- "<<par_upTrigCSV_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_upTrigCSV_KinFit.sg_p2<<" +- "<<par_upTrigCSV_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_upTrigCSV_KinFit.sg_p3<<" +- "<<par_upTrigCSV_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (CSV) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downTrigCSV_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_downTrig_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downTrigCSV_KinFit.sg_p0<<" +- "<<par_downTrigCSV_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downTrigCSV_KinFit.sg_p1<<" +- "<<par_downTrigCSV_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_downTrigCSV_KinFit.sg_p2<<" +- "<<par_downTrigCSV_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_downTrigCSV_KinFit.sg_p3<<" +- "<<par_downTrigCSV_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (pT) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upTrigpT_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_upTrig_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upTrigpT_KinFit.sg_p0<<" +- "<<par_upTrigpT_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upTrigpT_KinFit.sg_p1<<" +- "<<par_upTrigpT_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_upTrigpT_KinFit.sg_p2<<" +- "<<par_upTrigpT_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_upTrigpT_KinFit.sg_p3<<" +- "<<par_upTrigpT_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (pT) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downTrigpT_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_downTrig_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downTrigpT_KinFit.sg_p0<<" +- "<<par_downTrigpT_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downTrigpT_KinFit.sg_p1<<" +- "<<par_downTrigpT_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_downTrigpT_KinFit.sg_p2<<" +- "<<par_downTrigpT_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_downTrigpT_KinFit.sg_p3<<" +- "<<par_downTrigpT_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    outfile<<"   JEC lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_JECp1_KinFit->GetSumOfWeights(), h_mX_SR_JECm1_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   JER lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_JERp1_KinFit->GetSumOfWeights(), h_mX_SR_JERm1_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   Trigger CSV SF lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upTrigCSV_KinFit->GetSumOfWeights(), h_mX_SR_downTrigCSV_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   Trigger pT SF lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upTrigpT_KinFit->GetSumOfWeights(), h_mX_SR_downTrigpT_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    double sg_p0_errStat=par_KinFit.sg_p0_err;
    double sg_p0_errSyst[]={par_KinFit.sg_p0,
                            par_JECp1_KinFit.sg_p0, par_JECm1_KinFit.sg_p0,
                            par_JERp1_KinFit.sg_p0, par_JERm1_KinFit.sg_p0,
                            par_upTrigCSV_KinFit.sg_p0, par_downTrigCSV_KinFit.sg_p0,
                            par_upTrigpT_KinFit.sg_p0, par_downTrigpT_KinFit.sg_p0};
    double sg_p0_errSyst_min=par_KinFit.sg_p0-(*std::min_element(sg_p0_errSyst, sg_p0_errSyst+9));
    double sg_p0_errSyst_max=(*std::max_element(sg_p0_errSyst, sg_p0_errSyst+9))-par_KinFit.sg_p0;
    outfile<<"   Uncertainty on sg_p0 = "<<par_KinFit.sg_p0<<" +- "<<sg_p0_errStat<<" (stat) - "<<sg_p0_errSyst_min<<" + "<<sg_p0_errSyst_max<<" (syst); -"<<quad(sg_p0_errStat/2., sg_p0_errSyst_min)<<"/+"<<quad(sg_p0_errStat/2., sg_p0_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p1_errStat=par_KinFit.sg_p1_err;
    double sg_p1_errSyst[]={par_KinFit.sg_p1,
                            par_JECp1_KinFit.sg_p1, par_JECm1_KinFit.sg_p1,
                            par_JERp1_KinFit.sg_p1, par_JERm1_KinFit.sg_p1,
                            par_upTrigCSV_KinFit.sg_p1, par_downTrigCSV_KinFit.sg_p1,
                            par_upTrigpT_KinFit.sg_p1, par_downTrigpT_KinFit.sg_p1};
    double sg_p1_errSyst_min=par_KinFit.sg_p1-(*std::min_element(sg_p1_errSyst, sg_p1_errSyst+9));
    double sg_p1_errSyst_max=(*std::max_element(sg_p1_errSyst, sg_p1_errSyst+9))-par_KinFit.sg_p1;
    outfile<<"   Uncertainty on sg_p1 = "<<par_KinFit.sg_p1<<" +- "<<sg_p1_errStat<<" (stat) - "<<sg_p1_errSyst_min<<" + "<<sg_p1_errSyst_max<<" (syst); -"<<quad(sg_p1_errStat/2., sg_p1_errSyst_min)<<"/+"<<quad(sg_p1_errStat/2., sg_p1_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p2_errStat=par_KinFit.sg_p2_err;
    double sg_p2_errSyst[]={par_KinFit.sg_p2,
                            par_JECp1_KinFit.sg_p2, par_JECm1_KinFit.sg_p2,
                            par_JERp1_KinFit.sg_p2, par_JERm1_KinFit.sg_p2,
                            par_upTrigCSV_KinFit.sg_p2, par_downTrigCSV_KinFit.sg_p2,
                            par_upTrigpT_KinFit.sg_p2, par_downTrigpT_KinFit.sg_p2};
    double sg_p2_errSyst_min=par_KinFit.sg_p2-(*std::min_element(sg_p2_errSyst, sg_p2_errSyst+9));
    double sg_p2_errSyst_max=(*std::max_element(sg_p2_errSyst, sg_p2_errSyst+9))-par_KinFit.sg_p2;
    outfile<<"   Uncertainty on sg_p2 = "<<par_KinFit.sg_p2<<" +- "<<sg_p2_errStat<<" (stat) - "<<sg_p2_errSyst_min<<" + "<<sg_p2_errSyst_max<<" (syst); -"<<quad(sg_p2_errStat/2., sg_p2_errSyst_min)<<"/+"<<quad(sg_p2_errStat/2., sg_p2_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p3_errStat=par_KinFit.sg_p3_err;
    double sg_p3_errSyst[]={par_KinFit.sg_p3,
                            par_JECp1_KinFit.sg_p3, par_JECm1_KinFit.sg_p3,
                            par_JERp1_KinFit.sg_p3, par_JERm1_KinFit.sg_p3,
                            par_upTrigCSV_KinFit.sg_p3, par_downTrigCSV_KinFit.sg_p3,
                            par_upTrigpT_KinFit.sg_p3, par_downTrigpT_KinFit.sg_p3};
    double sg_p3_errSyst_min=par_KinFit.sg_p3-(*std::min_element(sg_p3_errSyst, sg_p3_errSyst+9));
    double sg_p3_errSyst_max=(*std::max_element(sg_p3_errSyst, sg_p3_errSyst+9))-par_KinFit.sg_p3;
    outfile<<"   Uncertainty on sg_p3 = "<<par_KinFit.sg_p3<<" +- "<<sg_p3_errStat<<" (stat) - "<<sg_p3_errSyst_min<<" + "<<sg_p3_errSyst_max<<" (syst); -"<<quad(sg_p3_errStat/2., sg_p3_errSyst_min)<<"/+"<<quad(sg_p3_errStat/2., sg_p3_errSyst_max)<<" (total) <br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    outfile<<"   Copy this into the datacard: <br/>"<<std::endl;
    outfile<<"   <pre>"<<std::endl;
    outfile<<"JEC       lnN     "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_JECp1_KinFit->GetSumOfWeights(), h_mX_SR_JECm1_KinFit->GetSumOfWeights())<<    "     -"<<std::endl;
    outfile<<"JER       lnN     "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_JERp1_KinFit->GetSumOfWeights(), h_mX_SR_JERm1_KinFit->GetSumOfWeights())<<    "     -"<<std::endl;
    outfile<<"bTagSF    lnN     1.1402      -"<<std::endl;
    outfile<<"trigSFCSV lnN     "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upTrigCSV_KinFit->GetSumOfWeights(), h_mX_SR_downTrigCSV_KinFit->GetSumOfWeights())<<"     -"<<std::endl;
    outfile<<"trigSFpT  lnN     "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upTrigpT_KinFit->GetSumOfWeights(), h_mX_SR_downTrigpT_KinFit->GetSumOfWeights())<<"     -"<<std::endl;
    outfile<<"signal_p0 param   "<<par_KinFit.sg_p0<<" -"<<quad(sg_p0_errStat/2., sg_p0_errSyst_min)<<"/+"<<quad(sg_p0_errStat/2., sg_p0_errSyst_max)<<std::endl;
    outfile<<"signal_p1 param   "<<par_KinFit.sg_p1<<" -"<<quad(sg_p1_errStat/2., sg_p1_errSyst_min)<<"/+"<<quad(sg_p1_errStat/2., sg_p1_errSyst_max)<<std::endl;
    outfile<<"signal_p2 param   "<<par_KinFit.sg_p2<<" -"<<quad(sg_p2_errStat/2., sg_p2_errSyst_min)<<"/+"<<quad(sg_p2_errStat/2., sg_p2_errSyst_max)<<std::endl;
    outfile<<"signal_p3 param   "<<par_KinFit.sg_p3<<" -"<<quad(sg_p3_errStat/2., sg_p3_errSyst_min)<<"/+"<<quad(sg_p3_errStat/2., sg_p3_errSyst_max)<<std::endl;
    outfile<<"   </pre>"<<std::endl;
    }
    outfile<<"   </div>"<<std::endl;
    outfile<<"  </td>"<<std::endl;
    outfile<<"  <td>"<<std::endl;
    outfile<<"   <img src='"<<("c_SignalmX_Gaussian_KinFit_"+masses.at(i)+".png")<<"'/><br/>"<<std::endl;
    outfile<<"   <h2 align='center'>With Kin-Fit. Fitted to a Gaussian.</h2><br/>"<<std::endl;
    outfile<<"   === Baseline plot === </br>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_KinFit->GetSumOfWeights()*totalLumi*prodXsec_1/nSignal_init<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_Gaussian_KinFit.sg_p0<<" +- "<<par_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_Gaussian_KinFit.sg_p1<<" +- "<<par_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    if (!focus) {
    outfile<<"   <center><input type='button' onclick=\"return toggleMe('"<<masses.at(i)<<"_Gaussian_KinFit')\" value='Systematics'></center><br/>"<<std::endl;
    outfile<<"   <div id=\""<<masses.at(i)<<"_Gaussian_KinFit\" style=\"display:none\">"<<std::endl;
    outfile<<"   === JEC +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECp1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JECp1_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECp1_Gaussian_KinFit.sg_p0<<" +- "<<par_JECp1_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECp1_Gaussian_KinFit.sg_p1<<" +- "<<par_JECp1_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   === JEC -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECm1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JECm1_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECm1_Gaussian_KinFit.sg_p0<<" +- "<<par_JECm1_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECm1_Gaussian_KinFit.sg_p1<<" +- "<<par_JECm1_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   === JER +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERp1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JERp1_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERp1_Gaussian_KinFit.sg_p0<<" +- "<<par_JERp1_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERp1_Gaussian_KinFit.sg_p1<<" +- "<<par_JERp1_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   === JER -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERm1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_JERm1_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERm1_Gaussian_KinFit.sg_p0<<" +- "<<par_JERm1_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERm1_Gaussian_KinFit.sg_p1<<" +- "<<par_JERm1_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (CSV) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upTrigCSV_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_upTrig_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upTrigCSV_Gaussian_KinFit.sg_p0<<" +- "<<par_upTrigCSV_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upTrigCSV_Gaussian_KinFit.sg_p1<<" +- "<<par_upTrigCSV_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (CSV) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downTrigCSV_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_downTrig_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downTrigCSV_Gaussian_KinFit.sg_p0<<" +- "<<par_downTrigCSV_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downTrigCSV_Gaussian_KinFit.sg_p1<<" +- "<<par_downTrigCSV_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (pT) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upTrigpT_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_upTrig_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upTrigpT_Gaussian_KinFit.sg_p0<<" +- "<<par_upTrigpT_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upTrigpT_Gaussian_KinFit.sg_p1<<" +- "<<par_upTrigpT_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger (pT) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downTrigpT_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    // outfile<<"   chi2/ndof = "<<plot_downTrig_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downTrigpT_Gaussian_KinFit.sg_p0<<" +- "<<par_downTrigpT_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downTrigpT_Gaussian_KinFit.sg_p1<<" +- "<<par_downTrigpT_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    outfile<<"   JEC lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_JECp1_KinFit->GetSumOfWeights(), h_mX_SR_JECm1_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   JER lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_JERp1_KinFit->GetSumOfWeights(), h_mX_SR_JERm1_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   btag SF lnN = 1.1402 <br/>"<<std::endl;
    outfile<<"   Trigger CSV SF lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upTrigCSV_KinFit->GetSumOfWeights(), h_mX_SR_downTrigCSV_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   Trigger pT SF lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upTrigpT_KinFit->GetSumOfWeights(), h_mX_SR_downTrigpT_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    double sg_p0_errStat=par_Gaussian_KinFit.sg_p0_err;
    double sg_p0_errSyst[]={par_Gaussian_KinFit.sg_p0,
                            par_JECp1_Gaussian_KinFit.sg_p0, par_JECm1_Gaussian_KinFit.sg_p0,
                            par_JERp1_Gaussian_KinFit.sg_p0, par_JERm1_Gaussian_KinFit.sg_p0,
                            par_upTrigCSV_Gaussian_KinFit.sg_p0, par_downTrigCSV_Gaussian_KinFit.sg_p0,
                            par_upTrigpT_Gaussian_KinFit.sg_p0, par_downTrigpT_Gaussian_KinFit.sg_p0};
    double sg_p0_errSyst_min=par_Gaussian_KinFit.sg_p0-(*std::min_element(sg_p0_errSyst, sg_p0_errSyst+9));
    double sg_p0_errSyst_max=(*std::max_element(sg_p0_errSyst, sg_p0_errSyst+9))-par_Gaussian_KinFit.sg_p0;
    outfile<<"   Uncertainty on sg_p0 = "<<par_Gaussian_KinFit.sg_p0<<" +- "<<sg_p0_errStat<<" (stat) - "<<sg_p0_errSyst_min<<" + "<<sg_p0_errSyst_max<<" (syst); -"<<quad(sg_p0_errStat/2., sg_p0_errSyst_min)<<"/+"<<quad(sg_p0_errStat/2., sg_p0_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p1_errStat=par_Gaussian_KinFit.sg_p1_err;
    double sg_p1_errSyst[]={par_Gaussian_KinFit.sg_p1,
                            par_JECp1_Gaussian_KinFit.sg_p1, par_JECm1_Gaussian_KinFit.sg_p1,
                            par_JERp1_Gaussian_KinFit.sg_p1, par_JERm1_Gaussian_KinFit.sg_p1,
                            par_upTrigCSV_Gaussian_KinFit.sg_p1, par_downTrigCSV_Gaussian_KinFit.sg_p1,
                            par_upTrigpT_Gaussian_KinFit.sg_p1, par_downTrigpT_Gaussian_KinFit.sg_p1};
    double sg_p1_errSyst_min=par_Gaussian_KinFit.sg_p1-(*std::min_element(sg_p1_errSyst, sg_p1_errSyst+9));
    double sg_p1_errSyst_max=(*std::max_element(sg_p1_errSyst, sg_p1_errSyst+9))-par_Gaussian_KinFit.sg_p1;
    outfile<<"   Uncertainty on sg_p1 = "<<par_Gaussian_KinFit.sg_p1<<" +- "<<sg_p1_errStat<<" (stat) - "<<sg_p1_errSyst_min<<" + "<<sg_p1_errSyst_max<<" (syst); -"<<quad(sg_p1_errStat/2., sg_p1_errSyst_min)<<"/+"<<quad(sg_p1_errStat/2., sg_p1_errSyst_max)<<" (total) <br/>"<<std::endl;
    }
    outfile<<"   </div>"<<std::endl;
    outfile<<"  </td>"<<std::endl;
    
    outfile<<" </tr>"<<std::endl;
    outfile<<"</table>"<<std::endl;
    
    if (writeDataCard)
    {
    
    
    
    }
    
    // Close all files
    file->Close();
    file_JECp1->Close();
    file_JECm1->Close();
    file_JERp1->Close();
    file_JERm1->Close();
    file_upTrigCSV->Close();
    file_downTrigCSV->Close();
    file_upTrigpT->Close();
    file_downTrigpT->Close();
    file_KinFit->Close();
    file_JECp1_KinFit->Close();
    file_JECm1_KinFit->Close();
    file_JERp1_KinFit->Close();
    file_JERm1_KinFit->Close();
    file_upTrigCSV_KinFit->Close();
    file_downTrigCSV_KinFit->Close();
    file_upTrigpT_KinFit->Close();
    file_downTrigpT_KinFit->Close();
    
  }
  
  std::vector<double> masses_d;
  for (unsigned int i=0; i<masses.size(); ++i) masses_d.push_back(atof(masses.at(i).c_str()));
  TGraphErrors *g_sg_p0=new TGraphErrors(masses.size()-4, &masses_d[4], &v_sg_p0[4], &v_zero[0], &v_sg_p0_err[4]);
  TGraphErrors *g_sg_p1=new TGraphErrors(masses.size()-4, &masses_d[4], &v_sg_p1[4], &v_zero[0], &v_sg_p1_err[4]);
  TGraphErrors *g_sg_p2=new TGraphErrors(masses.size()-4, &masses_d[4], &v_sg_p2[4], &v_zero[0], &v_sg_p2_err[4]);
  TGraphErrors *g_sg_p3=new TGraphErrors(masses.size()-4, &masses_d[4], &v_sg_p3[4], &v_zero[0], &v_sg_p3_err[4]);
  
  // TSpline3 *sp_sg_p0=new TGraphErrors("sp_sg_p0", &masses_d[4], &v_sg_p0[4], masses.size()-4, "b2e2", 0, 0);
  
  TCanvas *c_sg_p0=new TCanvas("c_sg_p0", "c_sg_p0", 700, 700);
  g_sg_p0->SetTitle("Signal Mean Interpolation; m_{X} (GeV); Signal Mean");
  g_sg_p0->Draw("AC*");
  g_sg_p0->Fit("pol1");
  c_sg_p0->SaveAs("SignalSystematics/c_sg_p0.png");
  
  TCanvas *c_sg_p1=new TCanvas("c_sg_p1", "c_sg_p1", 700, 700);
  g_sg_p1->SetTitle("Signal RMS Interpolation; m_{X} (GeV); Signal RMS");
  g_sg_p1->Draw("AC*");
  g_sg_p1->Fit("pol1");
  c_sg_p1->SaveAs("SignalSystematics/c_sg_p1.png");
  
  TCanvas *c_sg_p2=new TCanvas("c_sg_p2", "c_sg_p2", 700, 700);
  g_sg_p2->SetTitle("Signal Right Exponential Interpolation; m_{X} (GeV); Signal k_{right}");
  g_sg_p2->Draw("AC*");
  g_sg_p2->Fit("pol1");
  c_sg_p2->SaveAs("SignalSystematics/c_sg_p2.png");
  
  TCanvas *c_sg_p3=new TCanvas("c_sg_p3", "c_sg_p3", 700, 700);
  g_sg_p3->SetTitle("Signal Left Exponential Interpolation; m_{X} (GeV); Signal k_{left}");
  g_sg_p3->Draw("AC*");
  g_sg_p3->Fit("pol1");
  c_sg_p3->SaveAs("SignalSystematics/c_sg_p3.png");
  
  outfile<<"<h1> Signal Mean Interpolation Plot </h1>"<<std::endl;
  outfile<<"<img src='c_sg_p0.png'/> <br/> <hr/>"<<std::endl;
  outfile<<"<h1> Signal RMS Interpolation Plot </h1>"<<std::endl;
  outfile<<"<img src='c_sg_p1.png'/> <br/> <hr/>"<<std::endl;
  outfile<<"<h1> Signal Right Exponential Interpolation Plot </h1>"<<std::endl;
  outfile<<"<img src='c_sg_p2.png'/> <br/> <hr/>"<<std::endl;
  outfile<<"<h1> Signal Left Exponential Interpolation Plot </h1>"<<std::endl;
  outfile<<"<img src='c_sg_p3.png'/> <br/> <hr/>"<<std::endl;
  outfile<<"</body>"<<std::endl;
  outfile<<"</html>"<<std::endl;
  
  return 0;
}
    
