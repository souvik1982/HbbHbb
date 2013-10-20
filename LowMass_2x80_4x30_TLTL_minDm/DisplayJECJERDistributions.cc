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

// #include "/Users/souvik/HbbHbb/Analysis/PDFs/RevCrystalBall.h"

bool focus=true;
bool genMatch=false;

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
  double sg_p4;
  double sg_p5;
  double sg_p6;
  double sg_p7;
  double sg_p8;
  double sg_p0_err;
  double sg_p1_err;
  double sg_p2_err;
  double sg_p3_err;
  double sg_p4_err;
  double sg_p5_err;
  double sg_p6_err;
  double sg_p7_err;
  double sg_p8_err;
};

TLegend* twoStatBoxes(TH1F* h1, TH1F* h2)
{
  TLegend *leg=new TLegend(0.6, 0.6, 0.9, 0.9);
  double wherepeak=(h2->GetMean())/(h2->GetXaxis()->GetXmax());
  // std::cout<<"wherepeak = "<<wherepeak<<std::endl;
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
  // std::cout<<"wherepeak = "<<wherepeak<<std::endl;
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
  // std::cout<<"wherepeak = "<<wherepeak<<std::endl;
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

TLegend* elevenStatBoxes(TH1F* h1, TH1F* h2, TH1F *h3, TH1F *h4, TH1F *h5, TH1F *h6, TH1F *h7, TH1F *h8, TH1F *h9, TH1F *h10, TH1F *h11)
{
  TLegend *leg=new TLegend(0.7, 0.1, 0.9, 0.9);
  double wherepeak=(h2->GetMean())/(h2->GetXaxis()->GetXmax());
  // std::cout<<"wherepeak = "<<wherepeak<<std::endl;
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
  double rangeLo=-1, rangeHi=-1;
  if (!kinFit)
  {
    if (mass=="270")
    {
      rangeLo=220., rangeHi=500.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 260., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="300")
    {
      rangeLo=220., rangeHi=550.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 260., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="350")
    {
      rangeLo=250., rangeHi=550.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 320., 380.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="400")
    {
      rangeLo=250., rangeHi=550.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 370., 430.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 15.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="450")
    {
      rangeLo=300., rangeHi=550.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 420., 480.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 25.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="500")
    {
      rangeLo=300., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 470., 530.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="550")
    {
      rangeLo=350., rangeHi=650.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 520., 580.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
  }
  else
  {
    if (mass=="270")
    {
      rangeLo=220., rangeHi=500.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 260., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="300")
    {
      rangeLo=220., rangeHi=550.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 260., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="350")
    {
      rangeLo=250., rangeHi=550.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 320., 380.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 50.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="400")
    {
      rangeLo=250., rangeHi=550.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 370., 430.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 15.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="450")
    {
      rangeLo=300., rangeHi=550.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 420., 480.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 25.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="500")
    {
      rangeLo=300., rangeHi=600.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 470., 530.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
    else if (mass=="550")
    {
      rangeLo=350., rangeHi=650.; 
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 520., 580.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 30.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", 0., 5.);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", 0., 7.);
    }
  }
  x=new RooRealVar("x", "m_{X} (GeV)", rangeLo-100., rangeHi+100.);
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

RooPlot* fitSignal_GaussianCombinatoric(TH1F *h, std::string mass, int color, TLegend *leg, Params &params, bool kinFit=false)
{

  RooRealVar *x, *sg_p0, *sg_p1;
  x=new RooRealVar("x", "m_{X} (GeV)", 300., 1100.);
  double rangeLo=-1, rangeHi=-1;
  if (!kinFit)
  {
    if (mass=="300")
    {
      rangeLo=200., rangeHi=600.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 280., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 60.);
    }
    else if (mass=="400")
    {
      rangeLo=200., rangeHi=600.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 380., 420.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 10., 50.);
    }
  }
  else
  {
    if (mass=="270")
    {
      // rangeLo=220., rangeHi=550.;
      rangeLo=220., rangeHi=500.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 260., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 15.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", rangeLo);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", rangeHi);
      sg_p4=new RooRealVar("sg_p4", "sg_p4", -5., 0.);
      sg_p5=new RooRealVar("sg_p5", "sg_p5", 5., 15.);
      sg_p6=new RooRealVar("sg_p6", "sg_p6", -10., -2.);
      sg_p7=new RooRealVar("sg_p7", "sg_p7", 1., 8.);
      sg_p8=new RooRealVar("sg_p8", "sg_p8", 0., 1.);
    }
    if (mass=="300")
    {
      // rangeLo=220., rangeHi=550.;
      rangeLo=220., rangeHi=550.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 260., 320.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 15.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", rangeLo);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", rangeHi);
      sg_p4=new RooRealVar("sg_p4", "sg_p4", -5., 0.);
      sg_p5=new RooRealVar("sg_p5", "sg_p5", 5., 15.);
      sg_p6=new RooRealVar("sg_p6", "sg_p6", -10., -2.);
      sg_p7=new RooRealVar("sg_p7", "sg_p7", 1., 8.);
      sg_p8=new RooRealVar("sg_p8", "sg_p8", 0., 1.);
    }
    if (mass=="350")
    {
      // rangeLo=220., rangeHi=550.;
      rangeLo=250., rangeHi=550.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 320., 380.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 5., 15.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", rangeLo);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", rangeHi);
      sg_p4=new RooRealVar("sg_p4", "sg_p4", -5., 0.);
      sg_p5=new RooRealVar("sg_p5", "sg_p5", 5., 15.);
      sg_p6=new RooRealVar("sg_p6", "sg_p6", -10., -2.);
      sg_p7=new RooRealVar("sg_p7", "sg_p7", 1., 8.);
      sg_p8=new RooRealVar("sg_p8", "sg_p8", 0., 1.);
    }
    else if (mass=="400")
    {
      // rangeLo=300., rangeHi=600.;
      rangeLo=250., rangeHi=550.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 370., 430.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 3., 15.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", rangeLo);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", rangeHi);
      sg_p4=new RooRealVar("sg_p4", "sg_p4", -5., 0.);
      sg_p5=new RooRealVar("sg_p5", "sg_p5", 5., 15.);
      sg_p6=new RooRealVar("sg_p6", "sg_p6", -10., -2.);
      sg_p7=new RooRealVar("sg_p7", "sg_p7", 1., 8.);
      sg_p8=new RooRealVar("sg_p8", "sg_p8", 0., 1.);
    }
    else if (mass=="450")
    {
      rangeLo=300., rangeHi=550.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 420., 480.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 3., 15.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", rangeLo);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", rangeHi);
      sg_p4=new RooRealVar("sg_p4", "sg_p4", -5., 0.);
      sg_p5=new RooRealVar("sg_p5", "sg_p5", 5., 15.);
      sg_p6=new RooRealVar("sg_p6", "sg_p6", -10., -2.);
      sg_p7=new RooRealVar("sg_p7", "sg_p7", 1., 8.);
      sg_p8=new RooRealVar("sg_p8", "sg_p8", 0., 1.);
    }
    else if (mass=="500")
    {
      // rangeLo=300., rangeHi=600.;
      rangeLo=300., rangeHi=600.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 470., 530.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 3., 15.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", rangeLo);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", rangeHi);
      sg_p4=new RooRealVar("sg_p4", "sg_p4", -5., 0.);
      sg_p5=new RooRealVar("sg_p5", "sg_p5", 5., 15.);
      sg_p6=new RooRealVar("sg_p6", "sg_p6", -10., -2.);
      sg_p7=new RooRealVar("sg_p7", "sg_p7", 1., 8.);
      sg_p8=new RooRealVar("sg_p8", "sg_p8", 0., 1.);
    }
    else if (mass=="550")
    {
      // rangeLo=300., rangeHi=600.;
      rangeLo=350., rangeHi=650.;
      sg_p0=new RooRealVar("sg_p0", "sg_p0", 520., 580.);
      sg_p1=new RooRealVar("sg_p1", "sg_p1", 3., 15.);
      sg_p2=new RooRealVar("sg_p2", "sg_p2", rangeLo);
      sg_p3=new RooRealVar("sg_p3", "sg_p3", rangeHi);
      sg_p4=new RooRealVar("sg_p4", "sg_p4", -5., 0.);
      sg_p5=new RooRealVar("sg_p5", "sg_p5", 5., 15.);
      sg_p6=new RooRealVar("sg_p6", "sg_p6", -10., -2.);
      sg_p7=new RooRealVar("sg_p7", "sg_p7", 1., 8.);
      sg_p8=new RooRealVar("sg_p8", "sg_p8", 0., 1.);
    }
  }
  // sg_p2->setConstant(true);
  // sg_p3->setConstant(true);
  x=new RooRealVar("x", "m_{X} (GeV)", rangeLo-100., rangeHi+100.);
  RooGaussian signalCore("signalCore", "Signal Prediction", *x, *sg_p0, *sg_p1);
  // RooChebychev polyn("polyn", "Combinatoric", *x, RooArgList(*sg_p2, *sg_p3, *sg_p4));
  WrappedChebychev3 polyn("polyn", "Combinatoric", *x, *sg_p2, *sg_p3, *sg_p4, *sg_p5, *sg_p6, *sg_p7);
  RooAddPdf signal("signal", "signal", RooArgList(signalCore, polyn), *sg_p8);
  RooDataHist signalHistogram("signalHistogram", "Signal Histogram", RooArgList(*x), h);
  signal.fitTo(signalHistogram, RooFit::Range(rangeLo, rangeHi), RooFit::Save());
  //polyn.fitTo(signalHistogram, RooFit::Range(rangeLo, rangeHi), RooFit::Save());
  params.sg_p0=sg_p0->getVal(); params.sg_p0_err=sg_p0->getError();
  params.sg_p1=sg_p1->getVal(); params.sg_p1_err=sg_p1->getError();
  params.sg_p2=sg_p2->getVal(); params.sg_p2_err=sg_p2->getError();
  params.sg_p3=sg_p3->getVal(); params.sg_p3_err=sg_p3->getError();
  params.sg_p4=sg_p4->getVal(); params.sg_p4_err=sg_p4->getError();
  params.sg_p5=sg_p5->getVal(); params.sg_p5_err=sg_p5->getError();
  params.sg_p6=sg_p6->getVal(); params.sg_p6_err=sg_p6->getError();
  params.sg_p7=sg_p7->getVal(); params.sg_p7_err=sg_p7->getError();
  params.sg_p8=sg_p8->getVal(); params.sg_p8_err=sg_p8->getError();
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
  // Save modified signal shape to workspace
  if (color==kBlack)
  {
    RooRealVar signal_p0("signal_p0", "signal_p0", sg_p0->getVal());
    RooRealVar signal_p1("signal_p1", "signal_p1", sg_p1->getVal());
    RooRealVar signal_p2("signal_p2", "signal_p2", sg_p2->getVal());
    RooRealVar signal_p3("signal_p3", "signal_p3", sg_p3->getVal());
    RooRealVar signal_p4("signal_p4", "signal_p4", sg_p4->getVal());
    RooRealVar signal_p5("signal_p5", "signal_p5", sg_p5->getVal());
    RooRealVar signal_p6("signal_p6", "signal_p6", sg_p6->getVal());
    RooRealVar signal_p7("signal_p7", "signal_p7", sg_p7->getVal());
    RooRealVar signal_p8("signal_p8", "signal_p8", sg_p8->getVal());
    // RooChebychev polyn_fixed("polyn", "Combinatoric Prediction Fixed", *x, RooArgList(signal_p2, signal_p3, signal_p4));
    RooGaussian signalCore_fixed("signalCore", "Signal Prediction Fixed", *x, signal_p0, signal_p1); 
    WrappedChebychev3 polyn_fixed("polyn", "Combinatoric Prediction Fixed", *x, signal_p2, signal_p3, signal_p4, signal_p5, signal_p6, signal_p7);
    RooAddPdf signal_fixed("signal", "Signal Prediction Fixed", signalCore_fixed, polyn_fixed, signal_p8);
    RooWorkspace *w=new RooWorkspace("HbbHbb");
    w->import(signal_fixed);
    if (!kinFit) w->SaveAs(("SignalShapes/w_signal_GaussianComb_"+mass+".root").c_str());
    if (kinFit) w->SaveAs(("SignalShapes_KinFit/w_signal_GaussianComb_"+mass+".root").c_str());
  }
  return plot;
}

double lnN(double b, double a, double c)
{
  // std::cout<<"a = "<<a<<", b = "<<b<<", c = "<<c<<std::endl;
  // std::cout<<"1.+(a-c)/(2.*b) = "<<1.+fabs(a-c)/(2.*b)<<std::endl;
  return 1.+fabs(a-c)/(2.*b);
}

int DisplayJECJERDistributions()
{
  std::vector<std::string> masses;
  masses.push_back("270");
  masses.push_back("300");
  masses.push_back("350");
  masses.push_back("400");
  masses.push_back("450");
  masses.push_back("500");
  masses.push_back("550");
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  // gStyle->SetPalette(1);
  // gSystem->Load("/Users/souvik/HbbHbb/Analysis/PDFs/RevCrystalBall_cxx.so");
  // gSystem->Load("/Users/souvik/HbbHbb/Analysis/PDFs/GaussExp_cxx.so");
  gSystem->Load("/Users/souvik/HbbHbb/Analysis/PDFs/ExpGaussExp_cxx.so"); 
  gSystem->Load("/Users/souvik/HbbHbb/Analysis/PDFs/WrappedChebychev3_cxx.so");
  // Calculate nSignal events given production cross section, branching fractions and efficiency
  double totalLumi=18600.0; // /pb
  double prodXsec_1=1.; // pb
  
  // Write to an HTML File
  outfile.open("SignalSystematics/SignalSystematics.html");
  outfile<<"<html>"<<std::endl;
  outfile<<"<head>"<<std::endl;
  // outfile<<"<base href=\"https://cmslpcweb.fnal.gov/uscms_data/souvik/SignalSystematics\" target=\"_blank\">"<<std::endl;
  outfile<<"</head>"<<std::endl;
  outfile<<"<body>"<<std::endl;
  outfile<<"<h1> Signal with JEC, JER Systematic Uncertainties</h1>"<<std::endl;
  
  for (unsigned int i=0; i<masses.size(); ++i)
  {
    TFile *file=new TFile(("TLTL_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT=(TH1F*)file->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT=(TH1F*)file->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged=(TH1F*)file->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged->Add((TH1F*)file->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged=(TH1F*)file->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged->Add((TH1F*)file->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR=(TH1F*)file->Get("h_mX_SR");
    if (genMatch)
    {
      TH1F *h_mX_SR_rightComb=(TH1F*)file->Get("h_mX_SR_rightComb");
      TH1F *h_mX_SR_wrongComb=(TH1F*)file->Get("h_mX_SR_wrongComb");
      TH1F *h_mX_SR_noComb=(TH1F*)file->Get("h_mX_SR_noComb");
    }
    TH1F *h_CountWithPU=(TH1F*)file->Get("CountWithPU");
    double nSignal_init=h_CountWithPU->GetBinContent(1);
    std::cout<<"nSignal_init = "<<nSignal_init<<std::endl;
    
    TFile *file_JECp1=new TFile(("KinSel_JECp1/TLTL_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_JECp1=(TH1F*)file_JECp1->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_JECp1=(TH1F*)file_JECp1->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_JECp1=(TH1F*)file_JECp1->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_JECp1->Add((TH1F*)file_JECp1->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_JECp1=(TH1F*)file_JECp1->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_JECp1->Add((TH1F*)file_JECp1->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_JECp1=(TH1F*)file_JECp1->Get("h_mX_SR");
    
    TFile *file_JECm1=new TFile(("KinSel_JECm1/TLTL_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_JECm1=(TH1F*)file_JECm1->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_JECm1=(TH1F*)file_JECm1->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_JECm1=(TH1F*)file_JECm1->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_JECm1->Add((TH1F*)file_JECm1->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_JECm1=(TH1F*)file_JECm1->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_JECm1->Add((TH1F*)file_JECm1->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_JECm1=(TH1F*)file_JECm1->Get("h_mX_SR");
    
    TFile *file_JERp1=new TFile(("KinSel_JERp1/TLTL_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_JERp1=(TH1F*)file_JERp1->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_JERp1=(TH1F*)file_JERp1->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_JERp1=(TH1F*)file_JERp1->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_JERp1->Add((TH1F*)file_JERp1->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_JERp1=(TH1F*)file_JERp1->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_JERp1->Add((TH1F*)file_JERp1->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_JERp1=(TH1F*)file_JERp1->Get("h_mX_SR");
    
    TFile *file_JERm1=new TFile(("KinSel_JERm1/TLTL_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_JERm1=(TH1F*)file_JERm1->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_JERm1=(TH1F*)file_JERm1->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_JERm1=(TH1F*)file_JERm1->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_JERm1->Add((TH1F*)file_JERm1->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_JERm1=(TH1F*)file_JERm1->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_JERm1->Add((TH1F*)file_JERm1->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_JERm1=(TH1F*)file_JERm1->Get("h_mX_SR");
    
    TFile *file_upBC=new TFile(("KinSel_upBC/TLTL_upBC/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_upBC=(TH1F*)file_upBC->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_upBC=(TH1F*)file_upBC->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_upBC=(TH1F*)file_upBC->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_upBC->Add((TH1F*)file_upBC->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_upBC=(TH1F*)file_upBC->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_upBC->Add((TH1F*)file_upBC->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_upBC=(TH1F*)file_upBC->Get("h_mX_SR");
    
    TFile *file_downBC=new TFile(("KinSel_downBC/TLTL_downBC/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_downBC=(TH1F*)file_downBC->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_downBC=(TH1F*)file_downBC->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_downBC=(TH1F*)file_downBC->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_downBC->Add((TH1F*)file_downBC->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_downBC=(TH1F*)file_downBC->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_downBC->Add((TH1F*)file_downBC->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_downBC=(TH1F*)file_downBC->Get("h_mX_SR");
    
    TFile *file_upL=new TFile(("KinSel_upL/TLTL_upL/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_upL=(TH1F*)file_upL->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_upL=(TH1F*)file_upL->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_upL=(TH1F*)file_upL->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_upL->Add((TH1F*)file_upL->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_upL=(TH1F*)file_upL->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_upL->Add((TH1F*)file_upL->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_upL=(TH1F*)file_upL->Get("h_mX_SR");
    
    TFile *file_downL=new TFile(("KinSel_downL/TLTL_downL/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_downL=(TH1F*)file_downL->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_downL=(TH1F*)file_downL->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_downL=(TH1F*)file_downL->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_downL->Add((TH1F*)file_downL->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_downL=(TH1F*)file_downL->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_downL->Add((TH1F*)file_downL->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_downL=(TH1F*)file_downL->Get("h_mX_SR");
    
    TFile *file_upTrig=new TFile(("TLTL_upTrig/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_upTrig=(TH1F*)file_upTrig->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_upTrig=(TH1F*)file_upTrig->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_upTrig=(TH1F*)file_upTrig->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_upTrig->Add((TH1F*)file_upTrig->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_upTrig=(TH1F*)file_upTrig->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_upTrig->Add((TH1F*)file_upTrig->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_upTrig=(TH1F*)file_upTrig->Get("h_mX_SR");
    
    TFile *file_downTrig=new TFile(("TLTL_downTrig/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_H1Jet1pT_downTrig=(TH1F*)file_downTrig->Get("h_H1Jet1pT");
    TH1F *h_H1Jet2pT_downTrig=(TH1F*)file_downTrig->Get("h_H1Jet2pT");
    TH1F *h_H1_mass_bTagged_downTrig=(TH1F*)file_downTrig->Get("h_H1_mass_bTagged");
    h_H1_mass_bTagged_downTrig->Add((TH1F*)file_downTrig->Get("h_H2_mass_bTagged"));
    TH1F *h_H1_pT_bTagged_downTrig=(TH1F*)file_downTrig->Get("h_H1_pT_bTagged");
    h_H1_pT_bTagged_downTrig->Add((TH1F*)file_downTrig->Get("h_H2_pT_bTagged"));
    TH1F *h_mX_SR_downTrig=(TH1F*)file_downTrig->Get("h_mX_SR");
    
    // Kinematically Constrained
    TFile *file_KinFit=new TFile(("TLTL_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_KinFit=(TH1F*)file_KinFit->Get("h_mX_SR");
    if (genMatch)
    {
      TH1F *h_mX_SR_rightComb_KinFit=(TH1F*)file->Get("h_mX_SR_rightComb");
      TH1F *h_mX_SR_wrongComb_KinFit=(TH1F*)file->Get("h_mX_SR_wrongComb");
      TH1F *h_mX_SR_noComb_KinFit=(TH1F*)file->Get("h_mX_SR_noComb");
    }
    
    TFile *file_JECp1_KinFit=new TFile(("KinSel_JECp1/TLTL_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_JECp1_KinFit=(TH1F*)file_JECp1_KinFit->Get("h_mX_SR");
   
    TFile *file_JECm1_KinFit=new TFile(("KinSel_JECm1/TLTL_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_JECm1_KinFit=(TH1F*)file_JECm1_KinFit->Get("h_mX_SR");
    
    TFile *file_JERp1_KinFit=new TFile(("KinSel_JERp1/TLTL_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_JERp1_KinFit=(TH1F*)file_JERp1_KinFit->Get("h_mX_SR");
   
    TFile *file_JERm1_KinFit=new TFile(("KinSel_JERm1/TLTL_nominal/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_JERm1_KinFit=(TH1F*)file_JERm1_KinFit->Get("h_mX_SR");
    
    TFile *file_upBC_KinFit=new TFile(("KinSel_upBC/TLTL_upBC/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_upBC_KinFit=(TH1F*)file_upBC_KinFit->Get("h_mX_SR");
    
    TFile *file_downBC_KinFit=new TFile(("KinSel_downBC/TLTL_downBC/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_downBC_KinFit=(TH1F*)file_downBC_KinFit->Get("h_mX_SR");
    
    TFile *file_upL_KinFit=new TFile(("KinSel_upL/TLTL_upL/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_upL_KinFit=(TH1F*)file_upL_KinFit->Get("h_mX_SR");
    
    TFile *file_downL_KinFit=new TFile(("KinSel_downL/TLTL_downL/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_downL_KinFit=(TH1F*)file_downL_KinFit->Get("h_mX_SR");
    
    TFile *file_upTrig_KinFit=new TFile(("TLTL_upTrig/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_upTrig_KinFit=(TH1F*)file_upTrig_KinFit->Get("h_mX_SR");
    
    TFile *file_downTrig_KinFit=new TFile(("TLTL_downTrig/a_KinFit/Histograms_RadionToHH_4b_M-"+masses.at(i)+"_TuneZ2star_8TeV_FULLSIM.root").c_str());
    TH1F *h_mX_SR_downTrig_KinFit=(TH1F*)file_downTrig_KinFit->Get("h_mX_SR");
    
    if (!focus) {
    TCanvas *c_H1Jet1pT=new TCanvas("c_H1Jet1pT", "c_H1Jet1pT", 700, 700);
    h_H1Jet1pT_JECp1->SetLineStyle(9); h_H1Jet1pT_JECp1->SetLineColor(kRed);
    h_H1Jet1pT_JECm1->SetLineStyle(9); h_H1Jet1pT_JECm1->SetLineColor(kRed+2);
    h_H1Jet1pT_JERp1->SetLineStyle(9); h_H1Jet1pT_JERp1->SetLineColor(kMagenta);
    h_H1Jet1pT_JERm1->SetLineStyle(9); h_H1Jet1pT_JERm1->SetLineColor(kMagenta+2);
    h_H1Jet1pT_upBC->SetLineStyle(9); h_H1Jet1pT_upBC->SetLineColor(kGreen);
    h_H1Jet1pT_downBC->SetLineStyle(9); h_H1Jet1pT_downBC->SetLineColor(kGreen+3);
    h_H1Jet1pT_upL->SetLineStyle(9); h_H1Jet1pT_upL->SetLineColor(kCyan);
    h_H1Jet1pT_downL->SetLineStyle(9); h_H1Jet1pT_downL->SetLineColor(kCyan+2);
    h_H1Jet1pT_upTrig->SetLineStyle(9); h_H1Jet1pT_upTrig->SetLineColor(kBlue);
    h_H1Jet1pT_downTrig->SetLineStyle(9); h_H1Jet1pT_downTrig->SetLineColor(kBlue+2);
    h_H1Jet1pT_upTrig->Draw();
    h_H1Jet1pT->Draw("same");
    h_H1Jet1pT_JECp1->Draw("same");
    h_H1Jet1pT_JECm1->Draw("same");
    h_H1Jet1pT_JERp1->Draw("same");
    h_H1Jet1pT_JERm1->Draw("same");
    h_H1Jet1pT_upBC->Draw("same");
    h_H1Jet1pT_downBC->Draw("same");
    h_H1Jet1pT_upL->Draw("same");
    h_H1Jet1pT_downL->Draw("same");
    h_H1Jet1pT_upTrig->Draw("same");
    h_H1Jet1pT_downTrig->Draw("same");
    elevenStatBoxes(h_H1Jet1pT, h_H1Jet1pT_JECp1, h_H1Jet1pT_JECm1, h_H1Jet1pT_JERp1, h_H1Jet1pT_JERm1, h_H1Jet1pT_upBC, h_H1Jet1pT_downBC, h_H1Jet1pT_upL, h_H1Jet1pT_downL, h_H1Jet1pT_upTrig, h_H1Jet1pT_downTrig)->Draw();
    c_H1Jet1pT->SaveAs(("SignalSystematics/c_H1Jet1pT_"+masses.at(i)+".png").c_str());
    
    TCanvas *c_H1Jet2pT=new TCanvas("c_H1Jet2pT", "c_H1Jet2pT", 700, 700);
    h_H1Jet2pT_JECp1->SetLineStyle(9); h_H1Jet2pT_JECp1->SetLineColor(kRed);
    h_H1Jet2pT_JECm1->SetLineStyle(9); h_H1Jet2pT_JECm1->SetLineColor(kRed+2);
    h_H1Jet2pT_JERp1->SetLineStyle(9); h_H1Jet2pT_JERp1->SetLineColor(kMagenta);
    h_H1Jet2pT_JERm1->SetLineStyle(9); h_H1Jet2pT_JERm1->SetLineColor(kMagenta+2);
    h_H1Jet2pT_upBC->SetLineStyle(9); h_H1Jet2pT_upBC->SetLineColor(kGreen);
    h_H1Jet2pT_downBC->SetLineStyle(9); h_H1Jet2pT_downBC->SetLineColor(kGreen+3);
    h_H1Jet2pT_upL->SetLineStyle(9); h_H1Jet2pT_upL->SetLineColor(kCyan);
    h_H1Jet2pT_downL->SetLineStyle(9); h_H1Jet2pT_downL->SetLineColor(kCyan+2);
    h_H1Jet2pT_upTrig->SetLineStyle(9); h_H1Jet2pT_upTrig->SetLineColor(kBlue);
    h_H1Jet2pT_downTrig->SetLineStyle(9); h_H1Jet2pT_downTrig->SetLineColor(kBlue+2);
    h_H1Jet2pT_upTrig->Draw();
    h_H1Jet2pT->Draw("same");
    h_H1Jet2pT_JECp1->Draw("same");
    h_H1Jet2pT_JECm1->Draw("same");
    h_H1Jet2pT_JERp1->Draw("same");
    h_H1Jet2pT_JERm1->Draw("same");
    h_H1Jet2pT_upBC->Draw("same");
    h_H1Jet2pT_downBC->Draw("same");
    h_H1Jet2pT_upL->Draw("same");
    h_H1Jet2pT_downL->Draw("same");
    h_H1Jet2pT_upTrig->Draw("same");
    h_H1Jet2pT_downTrig->Draw("same");
    elevenStatBoxes(h_H1Jet2pT, h_H1Jet2pT_JECp1, h_H1Jet2pT_JECm1, h_H1Jet2pT_JERp1, h_H1Jet2pT_JERm1, h_H1Jet2pT_upBC, h_H1Jet2pT_downBC, h_H1Jet2pT_upL, h_H1Jet2pT_downL, h_H1Jet2pT_upTrig, h_H1Jet2pT_downTrig)->Draw();
    c_H1Jet2pT->SaveAs(("SignalSystematics/c_H1Jet2pT_"+masses.at(i)+".png").c_str());
    
    TCanvas *c_H1_mass_bTagged=new TCanvas("c_H1_mass_bTagged", "c_H1_mass_bTagged", 700, 700);
    h_H1_mass_bTagged_JECp1->SetLineStyle(9); h_H1_mass_bTagged_JECp1->SetLineColor(kRed);
    h_H1_mass_bTagged_JECm1->SetLineStyle(9); h_H1_mass_bTagged_JECm1->SetLineColor(kRed+2);
    h_H1_mass_bTagged_JERp1->SetLineStyle(9); h_H1_mass_bTagged_JERp1->SetLineColor(kMagenta);
    h_H1_mass_bTagged_JERm1->SetLineStyle(9); h_H1_mass_bTagged_JERm1->SetLineColor(kMagenta+2);
    h_H1_mass_bTagged_upBC->SetLineStyle(9); h_H1_mass_bTagged_upBC->SetLineColor(kGreen);
    h_H1_mass_bTagged_downBC->SetLineStyle(9); h_H1_mass_bTagged_downBC->SetLineColor(kGreen+3);
    h_H1_mass_bTagged_upL->SetLineStyle(9); h_H1_mass_bTagged_upL->SetLineColor(kCyan);
    h_H1_mass_bTagged_downL->SetLineStyle(9); h_H1_mass_bTagged_downL->SetLineColor(kCyan+2);
    h_H1_mass_bTagged_upTrig->SetLineStyle(9); h_H1_mass_bTagged_upTrig->SetLineColor(kBlue);
    h_H1_mass_bTagged_downTrig->SetLineStyle(9); h_H1_mass_bTagged_downTrig->SetLineColor(kBlue+2);
    h_H1_mass_bTagged_upTrig->Draw();
    h_H1_mass_bTagged->Draw("same");
    h_H1_mass_bTagged_JECp1->Draw("same");
    h_H1_mass_bTagged_JECm1->Draw("same");
    h_H1_mass_bTagged_JERp1->Draw("same");
    h_H1_mass_bTagged_JERm1->Draw("same");
    h_H1_mass_bTagged_upBC->Draw("same");
    h_H1_mass_bTagged_downBC->Draw("same");
    h_H1_mass_bTagged_upL->Draw("same");
    h_H1_mass_bTagged_downL->Draw("same");
    h_H1_mass_bTagged_upTrig->Draw("same");
    h_H1_mass_bTagged_downTrig->Draw("same");
    elevenStatBoxes(h_H1_mass_bTagged, h_H1_mass_bTagged_JECp1, h_H1_mass_bTagged_JECm1, h_H1_mass_bTagged_JERp1, h_H1_mass_bTagged_JERm1, h_H1_mass_bTagged_upBC, h_H1_mass_bTagged_downBC, h_H1_mass_bTagged_upL, h_H1_mass_bTagged_downL, h_H1_mass_bTagged_upTrig, h_H1_mass_bTagged_downTrig)->Draw();
    c_H1_mass_bTagged->SaveAs(("SignalSystematics/c_H1_mass_bTagged_"+masses.at(i)+".png").c_str());
    
    TCanvas *c_H1_pT_bTagged=new TCanvas("c_H1_pT_bTagged", "c_H1_pT_bTagged", 700, 700);
    h_H1_pT_bTagged_JECp1->SetLineStyle(9); h_H1_pT_bTagged_JECp1->SetLineColor(kRed);
    h_H1_pT_bTagged_JECm1->SetLineStyle(9); h_H1_pT_bTagged_JECm1->SetLineColor(kRed+2);
    h_H1_pT_bTagged_JERp1->SetLineStyle(9); h_H1_pT_bTagged_JERp1->SetLineColor(kMagenta);
    h_H1_pT_bTagged_JERm1->SetLineStyle(9); h_H1_pT_bTagged_JERm1->SetLineColor(kMagenta+2);
    h_H1_pT_bTagged_upBC->SetLineStyle(9); h_H1_pT_bTagged_upBC->SetLineColor(kGreen);
    h_H1_pT_bTagged_downBC->SetLineStyle(9); h_H1_pT_bTagged_downBC->SetLineColor(kGreen+3);
    h_H1_pT_bTagged_upL->SetLineStyle(9); h_H1_pT_bTagged_upL->SetLineColor(kCyan);
    h_H1_pT_bTagged_downL->SetLineStyle(9); h_H1_pT_bTagged_downL->SetLineColor(kCyan+2);
    h_H1_pT_bTagged_upTrig->SetLineStyle(9); h_H1_pT_bTagged_upTrig->SetLineColor(kBlue);
    h_H1_pT_bTagged_downTrig->SetLineStyle(9); h_H1_pT_bTagged_downTrig->SetLineColor(kBlue+2);
    h_H1_pT_bTagged_upTrig->Draw();
    h_H1_pT_bTagged->Draw("same");
    h_H1_pT_bTagged_JECp1->Draw("same");
    h_H1_pT_bTagged_JECm1->Draw("same");
    h_H1_pT_bTagged_JERp1->Draw("same");
    h_H1_pT_bTagged_JERm1->Draw("same");
    h_H1_pT_bTagged_upBC->Draw("same");
    h_H1_pT_bTagged_downBC->Draw("same");
    h_H1_pT_bTagged_upL->Draw("same");
    h_H1_pT_bTagged_downL->Draw("same");
    h_H1_pT_bTagged_upTrig->Draw("same");
    h_H1_pT_bTagged_downTrig->Draw("same");
    elevenStatBoxes(h_H1_pT_bTagged, h_H1_pT_bTagged_JECp1, h_H1_pT_bTagged_JECm1, h_H1_pT_bTagged_JERp1, h_H1_pT_bTagged_JERm1, h_H1_pT_bTagged_upBC, h_H1_pT_bTagged_downBC, h_H1_pT_bTagged_upL, h_H1_pT_bTagged_downL, h_H1_pT_bTagged_upTrig, h_H1_pT_bTagged_downTrig)->Draw();
    c_H1_pT_bTagged->SaveAs(("SignalSystematics/c_H1_pT_bTagged_"+masses.at(i)+".png").c_str());
    
    TCanvas *c_SignalmX=new TCanvas(("c_SignalmX"+masses.at(i)).c_str(), ("c_SignalmX"+masses.at(i)).c_str(), 700, 700);
    h_mX_SR->SetTitle(("m_{X} Peak in Signal MC (m_{X}="+masses.at(i)+" GeV); m_{X} (GeV)").c_str());
    h_mX_SR->Rebin(rebin);
    if (genMatch)
    {
      h_mX_SR_rightComb->Rebin(rebin);
      h_mX_SR_wrongComb->Rebin(rebin);
      h_mX_SR_noComb->Rebin(rebin);
    }
    h_mX_SR_JECp1->Rebin(rebin);
    h_mX_SR_JECm1->Rebin(rebin);
    h_mX_SR_JERp1->Rebin(rebin);
    h_mX_SR_JERm1->Rebin(rebin);
    h_mX_SR_upBC->Rebin(rebin);
    h_mX_SR_downBC->Rebin(rebin);
    h_mX_SR_upL->Rebin(rebin);
    h_mX_SR_downL->Rebin(rebin);
    h_mX_SR_upTrig->Rebin(rebin);
    h_mX_SR_downTrig->Rebin(rebin);
    if (genMatch)
    {
      h_mX_SR_rightComb->SetFillColor(kRed);
      h_mX_SR_wrongComb->SetFillColor(kGreen);
      h_mX_SR_noComb->SetFillColor(kBlue);
    }
    h_mX_SR_JECp1->SetLineColor(kRed);
    h_mX_SR_JECm1->SetLineColor(kRed+2);
    h_mX_SR_JERp1->SetLineColor(kMagenta);
    h_mX_SR_JERm1->SetLineColor(kMagenta+2);
    h_mX_SR_upBC->SetLineColor(kGreen);
    h_mX_SR_downBC->SetLineColor(kGreen+3);
    h_mX_SR_upL->SetLineColor(kCyan);
    h_mX_SR_downL->SetLineColor(kCyan+2);
    h_mX_SR_upTrig->SetLineColor(kBlue);
    h_mX_SR_downTrig->SetLineColor(kBlue+2);
    if (genMatch)
    {
      h_mX_SR_rightComb->GetXaxis()->SetRangeUser(0, 1200);
      h_mX_SR_wrongComb->GetXaxis()->SetRangeUser(0, 1200);
      h_mX_SR_noComb->GetXaxis()->SetRangeUser(0, 1200);
    }
    h_mX_SR->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JECp1->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JECm1->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JERp1->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JERm1->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_upBC->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_downBC->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_upL->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_downL->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_upTrig->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_downTrig->GetXaxis()->SetRangeUser(0, 1200);
    // h_mX_SR->Scale(totalLumi*prodXsec_1/nSignal_init);
    if (genMatch)
    {
      // h_mX_SR_rightComb->Scale(totalLumi*prodXsec_1/nSignal_init);
      // h_mX_SR_wrongComb->Scale(totalLumi*prodXsec_1/nSignal_init);
      // h_mX_SR_noComb->Scale(totalLumi*prodXsec_1/nSignal_init);
    }
    // h_mX_SR_JECp1->Scale(totalLumi*prodXsec_1/nSignal_init);
    // h_mX_SR_JECm1->Scale(totalLumi*prodXsec_1/nSignal_init);
    // h_mX_SR_JERp1->Scale(totalLumi*prodXsec_1/nSignal_init);
    // h_mX_SR_JERm1->Scale(totalLumi*prodXsec_1/nSignal_init);
    // h_mX_SR_upBC->Scale(totalLumi*prodXsec_1/nSignal_init);
    // h_mX_SR_downBC->Scale(totalLumi*prodXsec_1/nSignal_init);
    // h_mX_SR_upL->Scale(totalLumi*prodXsec_1/nSignal_init);
    // h_mX_SR_downL->Scale(totalLumi*prodXsec_1/nSignal_init);
    // h_mX_SR_upTrig->Scale(totalLumi*prodXsec_1/nSignal_init);
    // h_mX_SR_downTrig->Scale(totalLumi*prodXsec_1/nSignal_init);
    TLegend *leg=new TLegend(0.7, 0.2, 0.9, 0.9);
    leg->AddEntry(h_mX_SR, "Baseline");
    Params par, par_JECp1, par_JECm1, par_JERp1, par_JERm1, par_upBC, par_downBC, par_upL, par_downL, par_upTrig, par_downTrig;
    RooPlot *plot=fitSignal(h_mX_SR, masses.at(i), kBlack, leg, par);
    leg->AddEntry(h_mX_SR_JECp1, "JEC +1 #sigma");
    RooPlot *plot_JECp1=fitSignal(h_mX_SR_JECp1, masses.at(i), kRed, leg, par_JECp1);
    leg->AddEntry(h_mX_SR_JECm1, "JEC -1 #sigma");
    RooPlot *plot_JECm1=fitSignal(h_mX_SR_JECm1, masses.at(i), kRed+2, leg, par_JECm1);
    leg->AddEntry(h_mX_SR_JERp1, "JER +1 #sigma");
    RooPlot *plot_JERp1=fitSignal(h_mX_SR_JERp1, masses.at(i), kMagenta, leg, par_JERp1);
    leg->AddEntry(h_mX_SR_JERm1, "JER -1 #sigma");
    RooPlot *plot_JERm1=fitSignal(h_mX_SR_JERm1, masses.at(i), kMagenta+2, leg, par_JERm1);
    leg->AddEntry(h_mX_SR_upBC, "SF_{bc} +1 #sigma");
    RooPlot *plot_upBC=fitSignal(h_mX_SR_upBC, masses.at(i), kGreen, leg, par_upBC);
    leg->AddEntry(h_mX_SR_downBC, "SF_{bc} -1 #sigma");
    RooPlot *plot_downBC=fitSignal(h_mX_SR_downBC, masses.at(i), kGreen+3, leg, par_downBC);
    leg->AddEntry(h_mX_SR_upL, "SF_{l} +1 #sigma");
    RooPlot *plot_upL=fitSignal(h_mX_SR_upL, masses.at(i), kCyan, leg, par_upL);
    leg->AddEntry(h_mX_SR_downL, "SF_{l} -1 #sigma");
    RooPlot *plot_downL=fitSignal(h_mX_SR_downL, masses.at(i), kCyan+2, leg, par_downL);
    leg->AddEntry(h_mX_SR_upTrig, "Trig SF +1 #sigma");
    RooPlot *plot_upTrig=fitSignal(h_mX_SR_upTrig, masses.at(i), kBlue, leg, par_upTrig);
    leg->AddEntry(h_mX_SR_downTrig, "Trig SF -1 #sigma");
    RooPlot *plot_downTrig=fitSignal(h_mX_SR_downTrig, masses.at(i), kBlue+2, leg, par_downTrig);
    plot->SetMaximum(plot->GetMaximum()*1.1);
    plot->Draw();
    plot_JECp1->Draw("same");
    plot_JECm1->Draw("same");
    plot_JERp1->Draw("same");
    plot_JERm1->Draw("same");
    plot_upBC->Draw("same");
    plot_downBC->Draw("same");
    plot_upL->Draw("same");
    plot_downL->Draw("same");
    plot_upTrig->Draw("same");
    plot_downTrig->Draw("same");
    leg->SetFillColor(0);
    leg->Draw();
    c_SignalmX->SaveAs(("SignalSystematics/c_SignalmX_"+masses.at(i)+".png").c_str());
    }
    
    std::cout<<"0"<<std::endl;
    
    TCanvas *c_SignalmX_KinFit=new TCanvas(("c_SignalmX_KinFit"+masses.at(i)).c_str(), ("c_SignalmX_KinFit"+masses.at(i)).c_str(), 700, 700);
    h_mX_SR_KinFit->SetTitle(("m_{X} Peak in Signal MC (m_{X}="+masses.at(i)+" GeV); m_{X} (GeV)").c_str());
    h_mX_SR_KinFit->Rebin(rebin);
    if (genMatch)
    {
      h_mX_SR_rightComb_KinFit->Rebin(rebin);
      h_mX_SR_wrongComb_KinFit->Rebin(rebin);
      h_mX_SR_noComb_KinFit->Rebin(rebin);
    }
    h_mX_SR_JECp1_KinFit->Rebin(rebin);
    h_mX_SR_JECm1_KinFit->Rebin(rebin);
    h_mX_SR_JERp1_KinFit->Rebin(rebin);
    h_mX_SR_JERm1_KinFit->Rebin(rebin);
    h_mX_SR_upBC_KinFit->Rebin(rebin);
    h_mX_SR_downBC_KinFit->Rebin(rebin);
    h_mX_SR_upL_KinFit->Rebin(rebin);
    h_mX_SR_downL_KinFit->Rebin(rebin);
    h_mX_SR_upTrig_KinFit->Rebin(rebin);
    h_mX_SR_downTrig_KinFit->Rebin(rebin);
    if (genMatch)
    {
      h_mX_SR_rightComb_KinFit->SetFillColor(kRed);
      h_mX_SR_wrongComb_KinFit->SetFillColor(kGreen);
      h_mX_SR_noComb_KinFit->SetFillColor(kBlue);
    }
    h_mX_SR_JECp1_KinFit->SetLineColor(kRed);
    h_mX_SR_JECm1_KinFit->SetLineColor(kRed+2);
    h_mX_SR_JERp1_KinFit->SetLineColor(kMagenta);
    h_mX_SR_JERm1_KinFit->SetLineColor(kMagenta+2);
    h_mX_SR_upBC_KinFit->SetLineColor(kGreen);
    h_mX_SR_downBC_KinFit->SetLineColor(kGreen+3);
    h_mX_SR_upL_KinFit->SetLineColor(kCyan);
    h_mX_SR_downL_KinFit->SetLineColor(kCyan+2);
    h_mX_SR_upTrig_KinFit->SetLineColor(kBlue);
    h_mX_SR_downTrig_KinFit->SetLineColor(kBlue+2);
    if (genMatch)
    {
      h_mX_SR_rightComb_KinFit->GetXaxis()->SetRangeUser(0, 1200);
      h_mX_SR_wrongComb_KinFit->GetXaxis()->SetRangeUser(0, 1200);
      h_mX_SR_noComb_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    }
    h_mX_SR_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JECp1_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JECm1_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JERp1_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_JERm1_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_upBC_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_downBC_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_upL_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_downL_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_upTrig_KinFit->GetXaxis()->SetRangeUser(0, 1200);
    h_mX_SR_downTrig_KinFit->GetXaxis()->SetRangeUser(0, 1200);
//     h_mX_SR_KinFit->Scale(totalLumi*prodXsec_1/nSignal_init);
//     h_mX_SR_JECp1_KinFit->Scale(totalLumi*prodXsec_1/nSignal_init);
//     h_mX_SR_JECm1_KinFit->Scale(totalLumi*prodXsec_1/nSignal_init);
//     h_mX_SR_JERp1_KinFit->Scale(totalLumi*prodXsec_1/nSignal_init);
//     h_mX_SR_JERm1_KinFit->Scale(totalLumi*prodXsec_1/nSignal_init);
//     h_mX_SR_upBC_KinFit->Scale(totalLumi*prodXsec_1/nSignal_init);
//     h_mX_SR_downBC_KinFit->Scale(totalLumi*prodXsec_1/nSignal_init);
//     h_mX_SR_upL_KinFit->Scale(totalLumi*prodXsec_1/nSignal_init);
//     h_mX_SR_downL_KinFit->Scale(totalLumi*prodXsec_1/nSignal_init);
//     h_mX_SR_upTrig_KinFit->Scale(totalLumi*prodXsec_1/nSignal_init);
//     h_mX_SR_downTrig_KinFit->Scale(totalLumi*prodXsec_1/nSignal_init);
    TLegend *leg=new TLegend(0.7, 0.2, 0.9, 0.9);
    leg->AddEntry(h_mX_SR_KinFit, "Kin-Constrained Baseline");
    Params par_KinFit, par_JECp1_KinFit, par_JECm1_KinFit, par_JERp1_KinFit, par_JERm1_KinFit, par_upBC_KinFit, par_downBC_KinFit, par_upL_KinFit, par_downL_KinFit, par_upTrig_KinFit, par_downTrig_KinFit;
    RooPlot *plot_KinFit=fitSignal(h_mX_SR_KinFit, masses.at(i), kBlack, leg, par_KinFit, true);
    if (!focus) {
    leg->AddEntry(h_mX_SR_JECp1_KinFit, "JEC +1 #sigma");
    RooPlot *plot_JECp1_KinFit=fitSignal(h_mX_SR_JECp1_KinFit, masses.at(i), kRed, leg, par_JECp1_KinFit, true);
    leg->AddEntry(h_mX_SR_JECm1_KinFit, "JEC -1 #sigma");
    RooPlot *plot_JECm1_KinFit=fitSignal(h_mX_SR_JECm1_KinFit, masses.at(i), kRed+2, leg, par_JECm1_KinFit, true);
    leg->AddEntry(h_mX_SR_JERp1_KinFit, "JER +1 #sigma");
    RooPlot *plot_JERp1_KinFit=fitSignal(h_mX_SR_JERp1_KinFit, masses.at(i), kMagenta, leg, par_JERp1_KinFit, true);
    leg->AddEntry(h_mX_SR_JERm1_KinFit, "JER -1 #sigma");
    RooPlot *plot_JERm1_KinFit=fitSignal(h_mX_SR_JERm1_KinFit, masses.at(i), kMagenta+2, leg, par_JERm1_KinFit, true);
    leg->AddEntry(h_mX_SR_upBC_KinFit, "SF_{bc} +1 #sigma");
    RooPlot *plot_upBC_KinFit=fitSignal(h_mX_SR_upBC_KinFit, masses.at(i), kGreen, leg, par_upBC_KinFit, true);
    leg->AddEntry(h_mX_SR_downBC_KinFit, "SF_{bc} -1 #sigma");
    RooPlot *plot_downBC_KinFit=fitSignal(h_mX_SR_downBC_KinFit, masses.at(i), kGreen+3, leg, par_downBC_KinFit, true);
    leg->AddEntry(h_mX_SR_upL_KinFit, "SF_{l} +1 #sigma");
    RooPlot *plot_upL_KinFit=fitSignal(h_mX_SR_upL_KinFit, masses.at(i), kCyan, leg, par_upL_KinFit, true);
    leg->AddEntry(h_mX_SR_downL_KinFit, "SF_{l} -1 #sigma");
    RooPlot *plot_downL_KinFit=fitSignal(h_mX_SR_downL_KinFit, masses.at(i), kCyan+2, leg, par_downL_KinFit, true);
    leg->AddEntry(h_mX_SR_upTrig_KinFit, "Trig SF +1 #sigma");
    RooPlot *plot_upTrig_KinFit=fitSignal(h_mX_SR_upTrig_KinFit, masses.at(i), kBlue, leg, par_upTrig_KinFit, true);
    leg->AddEntry(h_mX_SR_downTrig_KinFit, "Trig SF -1 #sigma");
    RooPlot *plot_downTrig_KinFit=fitSignal(h_mX_SR_downTrig_KinFit, masses.at(i), kBlue+2, leg, par_downTrig_KinFit, true);
    }
    plot_KinFit->SetMaximum(plot_KinFit->GetMaximum()*1.2);
    plot_KinFit->Draw();
    if (!focus) {
    plot_JECp1_KinFit->Draw("same");
    plot_JECm1_KinFit->Draw("same");
    plot_JERp1_KinFit->Draw("same");
    plot_JERm1_KinFit->Draw("same");
    plot_upBC_KinFit->Draw("same");
    plot_downBC_KinFit->Draw("same");
    plot_upL_KinFit->Draw("same");
    plot_downL_KinFit->Draw("same");
    plot_upTrig_KinFit->Draw("same");
    plot_downTrig_KinFit->Draw("same");
    }
    leg->SetFillColor(0);
    leg->Draw();
    c_SignalmX_KinFit->SaveAs(("SignalSystematics/c_SignalmX_KinFit_"+masses.at(i)+".png").c_str());
    
    std::cout<<"1"<<std::endl;
    
    TCanvas *c_SignalmX_Gaussian_KinFit=new TCanvas(("c_SignalmX_Gaussian_KinFit"+masses.at(i)).c_str(), ("c_SignalmX_Gaussian_KinFit"+masses.at(i)).c_str(), 700, 700);
    TLegend *leg=new TLegend(0.7, 0.2, 0.9, 0.9);
    leg->AddEntry(h_mX_SR_KinFit, "Kin-Constrained Baseline - Gaussian Fit");
    Params par_Gaussian_KinFit, par_JECp1_Gaussian_KinFit, par_JECm1_Gaussian_KinFit, par_JERp1_Gaussian_KinFit, par_JERm1_Gaussian_KinFit, par_upBC_Gaussian_KinFit, par_downBC_Gaussian_KinFit, par_upL_Gaussian_KinFit, par_downL_Gaussian_KinFit, par_upTrig_Gaussian_KinFit, par_downTrig_Gaussian_KinFit;
    RooPlot *plot_Gaussian_KinFit=fitSignal_GaussianCombinatoric(h_mX_SR_KinFit, masses.at(i), kBlack, leg, par_Gaussian_KinFit, true);
    std::cout<<"1.1"<<std::endl;
    if (!focus) 
    {
      leg->AddEntry(h_mX_SR_JECp1_KinFit, "JEC +1 #sigma");
      RooPlot *plot_JECp1_Gaussian_KinFit=fitSignal_GaussianCombinatoric(h_mX_SR_JECp1_KinFit, masses.at(i), kRed, leg, par_JECp1_Gaussian_KinFit, true);
      leg->AddEntry(h_mX_SR_JECm1_KinFit, "JEC -1 #sigma");
      RooPlot *plot_JECm1_Gaussian_KinFit=fitSignal_GaussianCombinatoric(h_mX_SR_JECm1_KinFit, masses.at(i), kRed+2, leg, par_JECm1_Gaussian_KinFit, true);
      leg->AddEntry(h_mX_SR_JERp1_KinFit, "JER +1 #sigma");
      RooPlot *plot_JERp1_Gaussian_KinFit=fitSignal_GaussianCombinatoric(h_mX_SR_JERp1_KinFit, masses.at(i), kMagenta, leg, par_JERp1_Gaussian_KinFit, true);
      leg->AddEntry(h_mX_SR_JERm1_KinFit, "JER -1 #sigma");
      RooPlot *plot_JERm1_Gaussian_KinFit=fitSignal_GaussianCombinatoric(h_mX_SR_JERm1_KinFit, masses.at(i), kMagenta+2, leg, par_JERm1_Gaussian_KinFit, true);
      leg->AddEntry(h_mX_SR_upBC_KinFit, "SF_{bc} +1 #sigma");
      RooPlot *plot_upBC_Gaussian_KinFit=fitSignal_GaussianCombinatoric(h_mX_SR_upBC_KinFit, masses.at(i), kGreen, leg, par_upBC_Gaussian_KinFit, true);
      leg->AddEntry(h_mX_SR_downBC_KinFit, "SF_{bc} -1 #sigma");
      RooPlot *plot_downBC_Gaussian_KinFit=fitSignal_GaussianCombinatoric(h_mX_SR_downBC_KinFit, masses.at(i), kGreen+3, leg, par_downBC_Gaussian_KinFit, true);
      leg->AddEntry(h_mX_SR_upL_KinFit, "SF_{l} +1 #sigma");
      RooPlot *plot_upL_Gaussian_KinFit=fitSignal_GaussianCombinatoric(h_mX_SR_upL_KinFit, masses.at(i), kCyan, leg, par_upL_Gaussian_KinFit, true);
      leg->AddEntry(h_mX_SR_downL_KinFit, "SF_{l} -1 #sigma");
      RooPlot *plot_downL_Gaussian_KinFit=fitSignal_GaussianCombinatoric(h_mX_SR_downL_KinFit, masses.at(i), kCyan+2, leg, par_downL_Gaussian_KinFit, true);
      leg->AddEntry(h_mX_SR_upTrig_KinFit, "Trig SF +1 #sigma");
      RooPlot *plot_upTrig_Gaussian_KinFit=fitSignal_GaussianCombinatoric(h_mX_SR_upTrig_KinFit, masses.at(i), kBlue, leg, par_upTrig_Gaussian_KinFit, true);
      leg->AddEntry(h_mX_SR_downTrig_KinFit, "Trig SF -1 #sigma");
      RooPlot *plot_downTrig_Gaussian_KinFit=fitSignal_GaussianCombinatoric(h_mX_SR_downTrig_KinFit, masses.at(i), kBlue+2, leg, par_downTrig_Gaussian_KinFit, true);
    }
    plot_Gaussian_KinFit->SetMaximum(plot_Gaussian_KinFit->GetMaximum()*1.2);
    plot_Gaussian_KinFit->Draw();
    if (genMatch) drawCombinatorics(masses.at(i), h_mX_SR_rightComb_KinFit, h_mX_SR_wrongComb_KinFit, h_mX_SR_noComb_KinFit)->Draw("same");
    if (!focus) {
      plot_JECp1_Gaussian_KinFit->Draw("same");
      plot_JECm1_Gaussian_KinFit->Draw("same");
      plot_JERp1_Gaussian_KinFit->Draw("same");
      plot_JERm1_Gaussian_KinFit->Draw("same");
      plot_upBC_Gaussian_KinFit->Draw("same");
      plot_downBC_Gaussian_KinFit->Draw("same");
      plot_upL_Gaussian_KinFit->Draw("same");
      plot_downL_Gaussian_KinFit->Draw("same");
      plot_upTrig_Gaussian_KinFit->Draw("same");
      plot_downTrig_Gaussian_KinFit->Draw("same");
    }
    leg->SetFillColor(0);
    leg->Draw();
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
    outfile<<"   No kinematic-constraint. Exp-Gauss-Exp fit"<<std::endl;
    outfile<<"   <img src='"<<("c_SignalmX_"+masses.at(i)+".png")<<"'/><br/>"<<std::endl;
    if (!focus) {
    outfile<<"   === Baseline plot === </br>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par.sg_p0<<" +- "<<par.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par.sg_p1<<" +- "<<par.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par.sg_p2<<" +- "<<par.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par.sg_p3<<" +- "<<par.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JEC +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECp1->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JECp1->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECp1.sg_p0<<" +- "<<par_JECp1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECp1.sg_p1<<" +- "<<par_JECp1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JECp1.sg_p2<<" +- "<<par_JECp1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JECp1.sg_p3<<" +- "<<par_JECp1.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JEC -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECm1->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JECm1->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECm1.sg_p0<<" +- "<<par_JECm1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECm1.sg_p1<<" +- "<<par_JECm1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JECm1.sg_p2<<" +- "<<par_JECm1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JECm1.sg_p3<<" +- "<<par_JECm1.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JER +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERp1->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JERp1->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERp1.sg_p0<<" +- "<<par_JERp1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERp1.sg_p1<<" +- "<<par_JERp1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JERp1.sg_p2<<" +- "<<par_JERp1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JERp1.sg_p3<<" +- "<<par_JERp1.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JER -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERm1->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JERm1->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERm1.sg_p0<<" +- "<<par_JERm1.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERm1.sg_p1<<" +- "<<par_JERm1.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JERm1.sg_p2<<" +- "<<par_JERm1.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JERm1.sg_p3<<" +- "<<par_JERm1.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(bc) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upBC->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_upBC->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upBC.sg_p0<<" +- "<<par_upBC.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upBC.sg_p1<<" +- "<<par_upBC.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_upBC.sg_p2<<" +- "<<par_upBC.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_upBC.sg_p3<<" +- "<<par_upBC.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(bc) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downBC->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_downBC->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downBC.sg_p0<<" +- "<<par_downBC.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downBC.sg_p1<<" +- "<<par_downBC.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_downBC.sg_p2<<" +- "<<par_downBC.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_downBC.sg_p3<<" +- "<<par_downBC.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(l) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upL->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_upL->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upL.sg_p0<<" +- "<<par_upL.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upL.sg_p1<<" +- "<<par_upL.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_upL.sg_p2<<" +- "<<par_upL.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_upL.sg_p3<<" +- "<<par_upL.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(l) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downL->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_downL->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downL.sg_p0<<" +- "<<par_downL.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downL.sg_p1<<" +- "<<par_downL.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_downL.sg_p2<<" +- "<<par_downL.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_downL.sg_p3<<" +- "<<par_downL.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upTrig->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_upTrig->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upTrig.sg_p0<<" +- "<<par_upTrig.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upTrig.sg_p1<<" +- "<<par_upTrig.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_upTrig.sg_p2<<" +- "<<par_upTrig.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_upTrig.sg_p3<<" +- "<<par_upTrig.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downTrig->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_downTrig->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downTrig.sg_p0<<" +- "<<par_downTrig.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downTrig.sg_p1<<" +- "<<par_downTrig.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_downTrig.sg_p2<<" +- "<<par_downTrig.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_downTrig.sg_p3<<" +- "<<par_downTrig.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    outfile<<"   JEC lnN = "<<lnN(h_mX_SR->GetSumOfWeights(), h_mX_SR_JECp1->GetSumOfWeights(), h_mX_SR_JECm1->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   JER lnN = "<<lnN(h_mX_SR->GetSumOfWeights(), h_mX_SR_JERp1->GetSumOfWeights(), h_mX_SR_JERm1->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   btag SF(bc) lnN = "<<lnN(h_mX_SR->GetSumOfWeights(), h_mX_SR_upBC->GetSumOfWeights(), h_mX_SR_downBC->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   btag SF(l) lnN = "<<lnN(h_mX_SR->GetSumOfWeights(), h_mX_SR_upL->GetSumOfWeights(), h_mX_SR_downL->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   Trigger SF lnN = "<<lnN(h_mX_SR->GetSumOfWeights(), h_mX_SR_upTrig->GetSumOfWeights(), h_mX_SR_downTrig->GetSumOfWeights())<<"<br/>"<<std::endl;
    }
    outfile<<"  </td>"<<std::endl;
    outfile<<"  <td>"<<std::endl;
    outfile<<"   Kinematic-constraint. Exp-Gauss-Exp fit"<<std::endl;
    outfile<<"   <img src='"<<("c_SignalmX_KinFit_"+masses.at(i)+".png")<<"'/><br/>"<<std::endl;
    outfile<<"   === Baseline plot === </br>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_KinFit->GetSumOfWeights()*totalLumi*prodXsec_1/nSignal_init<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_KinFit.sg_p0<<" +- "<<par_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_KinFit.sg_p1<<" +- "<<par_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_KinFit.sg_p2<<" +- "<<par_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_KinFit.sg_p3<<" +- "<<par_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    if (!focus) {
    outfile<<"   === JEC +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECp1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JECp1_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECp1_KinFit.sg_p0<<" +- "<<par_JECp1_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECp1_KinFit.sg_p1<<" +- "<<par_JECp1_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JECp1_KinFit.sg_p2<<" +- "<<par_JECp1_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JECp1_KinFit.sg_p3<<" +- "<<par_JECp1_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JEC -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECm1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JECm1_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECm1_KinFit.sg_p0<<" +- "<<par_JECm1_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECm1_KinFit.sg_p1<<" +- "<<par_JECm1_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JECm1_KinFit.sg_p2<<" +- "<<par_JECm1_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JECm1_KinFit.sg_p3<<" +- "<<par_JECm1_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JER +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERp1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JERp1_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERp1_KinFit.sg_p0<<" +- "<<par_JERp1_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERp1_KinFit.sg_p1<<" +- "<<par_JERp1_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JERp1_KinFit.sg_p2<<" +- "<<par_JERp1_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JERp1_KinFit.sg_p3<<" +- "<<par_JERp1_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === JER -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERm1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JERm1_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERm1_KinFit.sg_p0<<" +- "<<par_JERm1_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERm1_KinFit.sg_p1<<" +- "<<par_JERm1_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_JERm1_KinFit.sg_p2<<" +- "<<par_JERm1_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_JERm1_KinFit.sg_p3<<" +- "<<par_JERm1_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(bc) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upBC_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_upBC_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upBC_KinFit.sg_p0<<" +- "<<par_upBC_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upBC_KinFit.sg_p1<<" +- "<<par_upBC_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_upBC_KinFit.sg_p2<<" +- "<<par_upBC_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_upBC_KinFit.sg_p3<<" +- "<<par_upBC_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(bc) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downBC_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_downBC_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downBC_KinFit.sg_p0<<" +- "<<par_downBC_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downBC_KinFit.sg_p1<<" +- "<<par_downBC_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_downBC_KinFit.sg_p2<<" +- "<<par_downBC_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_downBC_KinFit.sg_p3<<" +- "<<par_downBC_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(l) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upL_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_upL_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upL_KinFit.sg_p0<<" +- "<<par_upL_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upL_KinFit.sg_p1<<" +- "<<par_upL_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_upL_KinFit.sg_p2<<" +- "<<par_upL_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_upL_KinFit.sg_p3<<" +- "<<par_upL_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(l) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downL_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_downL_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downL_KinFit.sg_p0<<" +- "<<par_downL_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downL_KinFit.sg_p1<<" +- "<<par_downL_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_downL_KinFit.sg_p2<<" +- "<<par_downL_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_downL_KinFit.sg_p3<<" +- "<<par_downL_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upTrig_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_upTrig_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upTrig_KinFit.sg_p0<<" +- "<<par_upTrig_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upTrig_KinFit.sg_p1<<" +- "<<par_upTrig_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_upTrig_KinFit.sg_p2<<" +- "<<par_upTrig_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_upTrig_KinFit.sg_p3<<" +- "<<par_upTrig_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downTrig_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_downTrig_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downTrig_KinFit.sg_p0<<" +- "<<par_downTrig_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downTrig_KinFit.sg_p1<<" +- "<<par_downTrig_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p2 = "<<par_downTrig_KinFit.sg_p2<<" +- "<<par_downTrig_KinFit.sg_p2_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p3 = "<<par_downTrig_KinFit.sg_p3<<" +- "<<par_downTrig_KinFit.sg_p3_err<<" <br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    outfile<<"   JEC lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_JECp1_KinFit->GetSumOfWeights(), h_mX_SR_JECm1_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   JER lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_JERp1_KinFit->GetSumOfWeights(), h_mX_SR_JERm1_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   btag SF(bc) lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upBC_KinFit->GetSumOfWeights(), h_mX_SR_downBC_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   btag SF(l) lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upL_KinFit->GetSumOfWeights(), h_mX_SR_downL_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   Trigger SF lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upTrig_KinFit->GetSumOfWeights(), h_mX_SR_downTrig_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    double sg_p0_errStat=par_KinFit.sg_p0_err;
    double sg_p0_errSyst[]={par_KinFit.sg_p0,
                            par_JECp1_KinFit.sg_p0, par_JECm1_KinFit.sg_p0,
                            par_JERp1_KinFit.sg_p0, par_JERm1_KinFit.sg_p0,
                            par_upBC_KinFit.sg_p0, par_downBC_KinFit.sg_p0,
                            par_upL_KinFit.sg_p0, par_downL_KinFit.sg_p0,
                            par_upTrig_KinFit.sg_p0, par_downTrig_KinFit.sg_p0};
    double sg_p0_errSyst_min=par_KinFit.sg_p0-(*std::min_element(sg_p0_errSyst, sg_p0_errSyst+11));
    double sg_p0_errSyst_max=(*std::max_element(sg_p0_errSyst, sg_p0_errSyst+11))-par_KinFit.sg_p0;
    outfile<<"   Uncertainty on sg_p0 = "<<par_KinFit.sg_p0<<" +- "<<sg_p0_errStat<<" (stat) - "<<sg_p0_errSyst_min<<" + "<<sg_p0_errSyst_max<<" (syst); - "<<quad(sg_p0_errStat/2., sg_p0_errSyst_min)<<" + "<<quad(sg_p0_errStat/2., sg_p0_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p1_errStat=par_KinFit.sg_p1_err;
    double sg_p1_errSyst[]={par_KinFit.sg_p1,
                            par_JECp1_KinFit.sg_p1, par_JECm1_KinFit.sg_p1,
                            par_JERp1_KinFit.sg_p1, par_JERm1_KinFit.sg_p1,
                            par_upBC_KinFit.sg_p1, par_downBC_KinFit.sg_p1,
                            par_upL_KinFit.sg_p1, par_downL_KinFit.sg_p1,
                            par_upTrig_KinFit.sg_p1, par_downTrig_KinFit.sg_p1};
    double sg_p1_errSyst_min=par_KinFit.sg_p1-(*std::min_element(sg_p1_errSyst, sg_p1_errSyst+11));
    double sg_p1_errSyst_max=(*std::max_element(sg_p1_errSyst, sg_p1_errSyst+11))-par_KinFit.sg_p1;
    outfile<<"   Uncertainty on sg_p1 = "<<par_KinFit.sg_p1<<" +- "<<sg_p1_errStat<<" (stat) - "<<sg_p1_errSyst_min<<" + "<<sg_p1_errSyst_max<<" (syst); - "<<quad(sg_p1_errStat/2., sg_p1_errSyst_min)<<" + "<<quad(sg_p1_errStat/2., sg_p1_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p2_errStat=par_KinFit.sg_p2_err;
    double sg_p2_errSyst[]={par_KinFit.sg_p2,
                            par_JECp1_KinFit.sg_p2, par_JECm1_KinFit.sg_p2,
                            par_JERp1_KinFit.sg_p2, par_JERm1_KinFit.sg_p2,
                            par_upBC_KinFit.sg_p2, par_downBC_KinFit.sg_p2,
                            par_upL_KinFit.sg_p2, par_downL_KinFit.sg_p2,
                            par_upTrig_KinFit.sg_p2, par_downTrig_KinFit.sg_p2};
    double sg_p2_errSyst_min=par_KinFit.sg_p2-(*std::min_element(sg_p2_errSyst, sg_p2_errSyst+11));
    double sg_p2_errSyst_max=(*std::max_element(sg_p2_errSyst, sg_p2_errSyst+11))-par_KinFit.sg_p2;
    outfile<<"   Uncertainty on sg_p2 = "<<par_KinFit.sg_p2<<" +- "<<sg_p2_errStat<<" (stat) - "<<sg_p2_errSyst_min<<" + "<<sg_p2_errSyst_max<<" (syst); - "<<quad(sg_p2_errStat/2., sg_p2_errSyst_min)<<" + "<<quad(sg_p2_errStat/2., sg_p2_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p3_errStat=par_KinFit.sg_p3_err;
    double sg_p3_errSyst[]={par_KinFit.sg_p3,
                            par_JECp1_KinFit.sg_p3, par_JECm1_KinFit.sg_p3,
                            par_JERp1_KinFit.sg_p3, par_JERm1_KinFit.sg_p3,
                            par_upBC_KinFit.sg_p3, par_downBC_KinFit.sg_p3,
                            par_upL_KinFit.sg_p3, par_downL_KinFit.sg_p3,
                            par_upTrig_KinFit.sg_p3, par_downTrig_KinFit.sg_p3};
    double sg_p3_errSyst_min=par_KinFit.sg_p3-(*std::min_element(sg_p3_errSyst, sg_p3_errSyst+11));
    double sg_p3_errSyst_max=(*std::max_element(sg_p3_errSyst, sg_p3_errSyst+11))-par_KinFit.sg_p3;
    outfile<<"   Uncertainty on sg_p3 = "<<par_KinFit.sg_p3<<" +- "<<sg_p3_errStat<<" (stat) - "<<sg_p3_errSyst_min<<" + "<<sg_p3_errSyst_max<<" (syst); - "<<quad(sg_p3_errStat/2., sg_p3_errSyst_min)<<" + "<<quad(sg_p3_errStat/2., sg_p3_errSyst_max)<<" (total) <br/>"<<std::endl;
    }
    outfile<<"  </td>"<<std::endl;
    outfile<<"  <td>"<<std::endl;
    outfile<<"   Kinematic-constraint. Gaussian + polynomial fit"<<std::endl;
    outfile<<"   <img src='"<<("c_SignalmX_Gaussian_KinFit_"+masses.at(i)+".png")<<"'/><br/>"<<std::endl;
    outfile<<"   === Baseline plot === </br>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_KinFit->GetSumOfWeights()*totalLumi*prodXsec_1/nSignal_init<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_Gaussian_KinFit.sg_p0<<" +- "<<par_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_Gaussian_KinFit.sg_p1<<" +- "<<par_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p4 = "<<par_Gaussian_KinFit.sg_p4<<" +- "<<par_Gaussian_KinFit.sg_p4_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p5 = "<<par_Gaussian_KinFit.sg_p5<<" +- "<<par_Gaussian_KinFit.sg_p5_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p6 = "<<par_Gaussian_KinFit.sg_p6<<" +- "<<par_Gaussian_KinFit.sg_p6_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p7 = "<<par_Gaussian_KinFit.sg_p7<<" +- "<<par_Gaussian_KinFit.sg_p7_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p8 = "<<par_Gaussian_KinFit.sg_p8<<" +- "<<par_Gaussian_KinFit.sg_p8_err<<" <br/>"<<std::endl;
    if (!focus) {
    outfile<<"   === JEC +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECp1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JECp1_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECp1_Gaussian_KinFit.sg_p0<<" +- "<<par_JECp1_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECp1_Gaussian_KinFit.sg_p1<<" +- "<<par_JECp1_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p4 = "<<par_JECp1_Gaussian_KinFit.sg_p4<<" +- "<<par_JECp1_Gaussian_KinFit.sg_p4_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p5 = "<<par_JECp1_Gaussian_KinFit.sg_p5<<" +- "<<par_JECp1_Gaussian_KinFit.sg_p5_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p6 = "<<par_JECp1_Gaussian_KinFit.sg_p6<<" +- "<<par_JECp1_Gaussian_KinFit.sg_p6_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p7 = "<<par_JECp1_Gaussian_KinFit.sg_p7<<" +- "<<par_JECp1_Gaussian_KinFit.sg_p7_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p8 = "<<par_JECp1_Gaussian_KinFit.sg_p8<<" +- "<<par_JECp1_Gaussian_KinFit.sg_p8_err<<" <br/>"<<std::endl;
    outfile<<"   === JEC -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JECm1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JECm1_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JECm1_Gaussian_KinFit.sg_p0<<" +- "<<par_JECm1_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JECm1_Gaussian_KinFit.sg_p1<<" +- "<<par_JECm1_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p4 = "<<par_JECm1_Gaussian_KinFit.sg_p4<<" +- "<<par_JECm1_Gaussian_KinFit.sg_p4_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p5 = "<<par_JECm1_Gaussian_KinFit.sg_p5<<" +- "<<par_JECm1_Gaussian_KinFit.sg_p5_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p6 = "<<par_JECm1_Gaussian_KinFit.sg_p6<<" +- "<<par_JECm1_Gaussian_KinFit.sg_p6_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p7 = "<<par_JECm1_Gaussian_KinFit.sg_p7<<" +- "<<par_JECm1_Gaussian_KinFit.sg_p7_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p8 = "<<par_JECm1_Gaussian_KinFit.sg_p8<<" +- "<<par_JECm1_Gaussian_KinFit.sg_p8_err<<" <br/>"<<std::endl;
    outfile<<"   === JER +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERp1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JERp1_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERp1_Gaussian_KinFit.sg_p0<<" +- "<<par_JERp1_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERp1_Gaussian_KinFit.sg_p1<<" +- "<<par_JERp1_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p4 = "<<par_JERp1_Gaussian_KinFit.sg_p4<<" +- "<<par_JERp1_Gaussian_KinFit.sg_p4_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p5 = "<<par_JERp1_Gaussian_KinFit.sg_p5<<" +- "<<par_JERp1_Gaussian_KinFit.sg_p5_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p6 = "<<par_JERp1_Gaussian_KinFit.sg_p6<<" +- "<<par_JERp1_Gaussian_KinFit.sg_p6_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p7 = "<<par_JERp1_Gaussian_KinFit.sg_p7<<" +- "<<par_JERp1_Gaussian_KinFit.sg_p7_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p8 = "<<par_JERp1_Gaussian_KinFit.sg_p8<<" +- "<<par_JERp1_Gaussian_KinFit.sg_p8_err<<" <br/>"<<std::endl;
    outfile<<"   === JER -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_JERm1_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_JERm1_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_JERm1_Gaussian_KinFit.sg_p0<<" +- "<<par_JERm1_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_JERm1_Gaussian_KinFit.sg_p1<<" +- "<<par_JERm1_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p4 = "<<par_JERm1_Gaussian_KinFit.sg_p4<<" +- "<<par_JERm1_Gaussian_KinFit.sg_p4_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p5 = "<<par_JERm1_Gaussian_KinFit.sg_p5<<" +- "<<par_JERm1_Gaussian_KinFit.sg_p5_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p6 = "<<par_JERm1_Gaussian_KinFit.sg_p6<<" +- "<<par_JERm1_Gaussian_KinFit.sg_p6_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p7 = "<<par_JERm1_Gaussian_KinFit.sg_p7<<" +- "<<par_JERm1_Gaussian_KinFit.sg_p7_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p8 = "<<par_JERm1_Gaussian_KinFit.sg_p8<<" +- "<<par_JERm1_Gaussian_KinFit.sg_p8_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(bc) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upBC_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_upBC_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upBC_Gaussian_KinFit.sg_p0<<" +- "<<par_upBC_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upBC_Gaussian_KinFit.sg_p1<<" +- "<<par_upBC_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p4 = "<<par_upBC_Gaussian_KinFit.sg_p4<<" +- "<<par_upBC_Gaussian_KinFit.sg_p4_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p5 = "<<par_upBC_Gaussian_KinFit.sg_p5<<" +- "<<par_upBC_Gaussian_KinFit.sg_p5_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p6 = "<<par_upBC_Gaussian_KinFit.sg_p6<<" +- "<<par_upBC_Gaussian_KinFit.sg_p6_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p7 = "<<par_upBC_Gaussian_KinFit.sg_p7<<" +- "<<par_upBC_Gaussian_KinFit.sg_p7_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p8 = "<<par_upBC_Gaussian_KinFit.sg_p8<<" +- "<<par_upBC_Gaussian_KinFit.sg_p8_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(bc) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downBC_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_downBC_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downBC_Gaussian_KinFit.sg_p0<<" +- "<<par_downBC_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downBC_Gaussian_KinFit.sg_p1<<" +- "<<par_downBC_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p4 = "<<par_downBC_Gaussian_KinFit.sg_p4<<" +- "<<par_downBC_Gaussian_KinFit.sg_p4_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p5 = "<<par_downBC_Gaussian_KinFit.sg_p5<<" +- "<<par_downBC_Gaussian_KinFit.sg_p5_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p6 = "<<par_downBC_Gaussian_KinFit.sg_p6<<" +- "<<par_downBC_Gaussian_KinFit.sg_p6_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p7 = "<<par_downBC_Gaussian_KinFit.sg_p7<<" +- "<<par_downBC_Gaussian_KinFit.sg_p7_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p8 = "<<par_downBC_Gaussian_KinFit.sg_p8<<" +- "<<par_downBC_Gaussian_KinFit.sg_p8_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(l) +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upL_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_upL_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upL_Gaussian_KinFit.sg_p0<<" +- "<<par_upL_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upL_Gaussian_KinFit.sg_p1<<" +- "<<par_upL_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p4 = "<<par_upL_Gaussian_KinFit.sg_p4<<" +- "<<par_upL_Gaussian_KinFit.sg_p4_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p5 = "<<par_upL_Gaussian_KinFit.sg_p5<<" +- "<<par_upL_Gaussian_KinFit.sg_p5_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p6 = "<<par_upL_Gaussian_KinFit.sg_p6<<" +- "<<par_upL_Gaussian_KinFit.sg_p6_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p7 = "<<par_upL_Gaussian_KinFit.sg_p7<<" +- "<<par_upL_Gaussian_KinFit.sg_p7_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p8 = "<<par_upL_Gaussian_KinFit.sg_p8<<" +- "<<par_upL_Gaussian_KinFit.sg_p8_err<<" <br/>"<<std::endl;
    outfile<<"   === btag SF(l) -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downL_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_downL_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downL_Gaussian_KinFit.sg_p0<<" +- "<<par_downL_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downL_Gaussian_KinFit.sg_p1<<" +- "<<par_downL_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p4 = "<<par_downL_Gaussian_KinFit.sg_p4<<" +- "<<par_downL_Gaussian_KinFit.sg_p4_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p5 = "<<par_downL_Gaussian_KinFit.sg_p5<<" +- "<<par_downL_Gaussian_KinFit.sg_p5_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p6 = "<<par_downL_Gaussian_KinFit.sg_p6<<" +- "<<par_downL_Gaussian_KinFit.sg_p6_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p7 = "<<par_downL_Gaussian_KinFit.sg_p7<<" +- "<<par_downL_Gaussian_KinFit.sg_p7_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p8 = "<<par_downL_Gaussian_KinFit.sg_p8<<" +- "<<par_downL_Gaussian_KinFit.sg_p8_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger +1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_upTrig_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_upTrig_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_upTrig_Gaussian_KinFit.sg_p0<<" +- "<<par_upTrig_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_upTrig_Gaussian_KinFit.sg_p1<<" +- "<<par_upTrig_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p4 = "<<par_upTrig_Gaussian_KinFit.sg_p4<<" +- "<<par_upTrig_Gaussian_KinFit.sg_p4_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p5 = "<<par_upTrig_Gaussian_KinFit.sg_p5<<" +- "<<par_upTrig_Gaussian_KinFit.sg_p5_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p6 = "<<par_upTrig_Gaussian_KinFit.sg_p6<<" +- "<<par_upTrig_Gaussian_KinFit.sg_p6_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p7 = "<<par_upTrig_Gaussian_KinFit.sg_p7<<" +- "<<par_upTrig_Gaussian_KinFit.sg_p7_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p8 = "<<par_upTrig_Gaussian_KinFit.sg_p8<<" +- "<<par_upTrig_Gaussian_KinFit.sg_p8_err<<" <br/>"<<std::endl;
    outfile<<"   === Trigger -1 sigma === <br/>"<<std::endl;
    outfile<<"   norm = "<<h_mX_SR_downTrig_KinFit->GetSumOfWeights()<<" <br/>"<<std::endl;
    outfile<<"   chi2/ndof = "<<plot_downTrig_Gaussian_KinFit->chiSquare()<<" <br/>"<<std::endl;
    outfile<<"   sg_p0 = "<<par_downTrig_Gaussian_KinFit.sg_p0<<" +- "<<par_downTrig_Gaussian_KinFit.sg_p0_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p1 = "<<par_downTrig_Gaussian_KinFit.sg_p1<<" +- "<<par_downTrig_Gaussian_KinFit.sg_p1_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p4 = "<<par_downTrig_Gaussian_KinFit.sg_p4<<" +- "<<par_downTrig_Gaussian_KinFit.sg_p4_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p5 = "<<par_downTrig_Gaussian_KinFit.sg_p5<<" +- "<<par_downTrig_Gaussian_KinFit.sg_p5_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p6 = "<<par_downTrig_Gaussian_KinFit.sg_p6<<" +- "<<par_downTrig_Gaussian_KinFit.sg_p6_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p7 = "<<par_downTrig_Gaussian_KinFit.sg_p7<<" +- "<<par_downTrig_Gaussian_KinFit.sg_p7_err<<" <br/>"<<std::endl;
    outfile<<"   sg_p8 = "<<par_downTrig_Gaussian_KinFit.sg_p8<<" +- "<<par_downTrig_Gaussian_KinFit.sg_p8_err<<" <br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    outfile<<"   JEC lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_JECp1_KinFit->GetSumOfWeights(), h_mX_SR_JECm1_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   JER lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_JERp1_KinFit->GetSumOfWeights(), h_mX_SR_JERm1_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   btag SF(bc) lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upBC_KinFit->GetSumOfWeights(), h_mX_SR_downBC_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   btag SF(l) lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upL_KinFit->GetSumOfWeights(), h_mX_SR_downL_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   Trigger SF lnN = "<<lnN(h_mX_SR_KinFit->GetSumOfWeights(), h_mX_SR_upTrig_KinFit->GetSumOfWeights(), h_mX_SR_downTrig_KinFit->GetSumOfWeights())<<"<br/>"<<std::endl;
    outfile<<"   === === <br/>"<<std::endl;
    double sg_p0_errStat=par_Gaussian_KinFit.sg_p0_err;
    double sg_p0_errSyst[]={par_Gaussian_KinFit.sg_p0,
                            par_JECp1_Gaussian_KinFit.sg_p0, par_JECm1_Gaussian_KinFit.sg_p0,
                            par_JERp1_Gaussian_KinFit.sg_p0, par_JERm1_Gaussian_KinFit.sg_p0,
                            par_upBC_Gaussian_KinFit.sg_p0, par_downBC_Gaussian_KinFit.sg_p0,
                            par_upL_Gaussian_KinFit.sg_p0, par_downL_Gaussian_KinFit.sg_p0,
                           par_upTrig_Gaussian_KinFit.sg_p0, par_downTrig_Gaussian_KinFit.sg_p0};
    double sg_p0_errSyst_min=par_Gaussian_KinFit.sg_p0-(*std::min_element(sg_p0_errSyst, sg_p0_errSyst+11));
    double sg_p0_errSyst_max=(*std::max_element(sg_p0_errSyst, sg_p0_errSyst+11))-par_Gaussian_KinFit.sg_p0;
    outfile<<"   Uncertainty on sg_p0 = "<<par_Gaussian_KinFit.sg_p0<<" +- "<<sg_p0_errStat<<" (stat) - "<<sg_p0_errSyst_min<<" + "<<sg_p0_errSyst_max<<" (syst); -"<<quad(sg_p0_errStat/2., sg_p0_errSyst_min)<<"/+"<<quad(sg_p0_errStat/2., sg_p0_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p1_errStat=par_Gaussian_KinFit.sg_p1_err;
    double sg_p1_errSyst[]={par_Gaussian_KinFit.sg_p1,
                            par_JECp1_Gaussian_KinFit.sg_p1, par_JECm1_Gaussian_KinFit.sg_p1,
                            par_JERp1_Gaussian_KinFit.sg_p1, par_JERm1_Gaussian_KinFit.sg_p1,
                            par_upBC_Gaussian_KinFit.sg_p1, par_downBC_Gaussian_KinFit.sg_p1,
                            par_upL_Gaussian_KinFit.sg_p1, par_downL_Gaussian_KinFit.sg_p1,
                            par_upTrig_Gaussian_KinFit.sg_p1, par_downTrig_Gaussian_KinFit.sg_p1};
    double sg_p1_errSyst_min=par_Gaussian_KinFit.sg_p1-(*std::min_element(sg_p1_errSyst, sg_p1_errSyst+11));
    double sg_p1_errSyst_max=(*std::max_element(sg_p1_errSyst, sg_p1_errSyst+11))-par_Gaussian_KinFit.sg_p1;
    outfile<<"   Uncertainty on sg_p1 = "<<par_Gaussian_KinFit.sg_p1<<" +- "<<sg_p1_errStat<<" (stat) - "<<sg_p1_errSyst_min<<" + "<<sg_p1_errSyst_max<<" (syst); -"<<quad(sg_p1_errStat/2., sg_p1_errSyst_min)<<"/+"<<quad(sg_p1_errStat/2., sg_p1_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p4_errStat=par_Gaussian_KinFit.sg_p4_err;
    double sg_p4_errSyst[]={par_Gaussian_KinFit.sg_p4,
                            par_JECp1_Gaussian_KinFit.sg_p4, par_JECm1_Gaussian_KinFit.sg_p4,
                            par_JERp1_Gaussian_KinFit.sg_p4, par_JERm1_Gaussian_KinFit.sg_p4,
                            par_upBC_Gaussian_KinFit.sg_p4, par_downBC_Gaussian_KinFit.sg_p4,
                            par_upL_Gaussian_KinFit.sg_p4, par_downL_Gaussian_KinFit.sg_p4,
                            par_upTrig_Gaussian_KinFit.sg_p4, par_downTrig_Gaussian_KinFit.sg_p4};
    double sg_p4_errSyst_min=par_Gaussian_KinFit.sg_p4-(*std::min_element(sg_p4_errSyst, sg_p4_errSyst+11));
    double sg_p4_errSyst_max=(*std::max_element(sg_p4_errSyst, sg_p4_errSyst+11))-par_Gaussian_KinFit.sg_p4;
    outfile<<"   Uncertainty on sg_p4 = "<<par_Gaussian_KinFit.sg_p4<<" +- "<<sg_p4_errStat<<" (stat) - "<<sg_p4_errSyst_min<<" + "<<sg_p4_errSyst_max<<" (syst); -"<<quad(sg_p4_errStat/2., sg_p4_errSyst_min)<<"/+"<<quad(sg_p4_errStat/2., sg_p4_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p5_errStat=par_Gaussian_KinFit.sg_p5_err;
    double sg_p5_errSyst[]={par_Gaussian_KinFit.sg_p5,
                            par_JECp1_Gaussian_KinFit.sg_p5, par_JECm1_Gaussian_KinFit.sg_p5,
                            par_JERp1_Gaussian_KinFit.sg_p5, par_JERm1_Gaussian_KinFit.sg_p5,
                            par_upBC_Gaussian_KinFit.sg_p5, par_downBC_Gaussian_KinFit.sg_p5,
                            par_upL_Gaussian_KinFit.sg_p5, par_downL_Gaussian_KinFit.sg_p5,
                            par_upTrig_Gaussian_KinFit.sg_p5, par_downTrig_Gaussian_KinFit.sg_p5};
    double sg_p5_errSyst_min=par_Gaussian_KinFit.sg_p5-(*std::min_element(sg_p5_errSyst, sg_p5_errSyst+11));
    double sg_p5_errSyst_max=(*std::max_element(sg_p5_errSyst, sg_p5_errSyst+11))-par_Gaussian_KinFit.sg_p5;
    outfile<<"   Uncertainty on sg_p5 = "<<par_Gaussian_KinFit.sg_p5<<" +- "<<sg_p5_errStat<<" (stat) - "<<sg_p5_errSyst_min<<" + "<<sg_p5_errSyst_max<<" (syst); -"<<quad(sg_p5_errStat/2., sg_p5_errSyst_min)<<"/+"<<quad(sg_p5_errStat/2., sg_p5_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p6_errStat=par_Gaussian_KinFit.sg_p6_err;
    double sg_p6_errSyst[]={par_Gaussian_KinFit.sg_p6,
                            par_JECp1_Gaussian_KinFit.sg_p6, par_JECm1_Gaussian_KinFit.sg_p6,
                            par_JERp1_Gaussian_KinFit.sg_p6, par_JERm1_Gaussian_KinFit.sg_p6,
                            par_upBC_Gaussian_KinFit.sg_p6, par_downBC_Gaussian_KinFit.sg_p6,
                            par_upL_Gaussian_KinFit.sg_p6, par_downL_Gaussian_KinFit.sg_p6,
                            par_upTrig_Gaussian_KinFit.sg_p6, par_downTrig_Gaussian_KinFit.sg_p6};
    double sg_p6_errSyst_min=par_Gaussian_KinFit.sg_p6-(*std::min_element(sg_p6_errSyst, sg_p6_errSyst+11));
    double sg_p6_errSyst_max=(*std::max_element(sg_p6_errSyst, sg_p6_errSyst+11))-par_Gaussian_KinFit.sg_p6;
    outfile<<"   Uncertainty on sg_p6 = "<<par_Gaussian_KinFit.sg_p6<<" +- "<<sg_p6_errStat<<" (stat) - "<<sg_p6_errSyst_min<<" + "<<sg_p6_errSyst_max<<" (syst); -"<<quad(sg_p6_errStat/2., sg_p6_errSyst_min)<<"/+"<<quad(sg_p6_errStat/2., sg_p6_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p7_errStat=par_Gaussian_KinFit.sg_p7_err;
    double sg_p7_errSyst[]={par_Gaussian_KinFit.sg_p7,
                            par_JECp1_Gaussian_KinFit.sg_p7, par_JECm1_Gaussian_KinFit.sg_p7,
                            par_JERp1_Gaussian_KinFit.sg_p7, par_JERm1_Gaussian_KinFit.sg_p7,
                            par_upBC_Gaussian_KinFit.sg_p7, par_downBC_Gaussian_KinFit.sg_p7,
                            par_upL_Gaussian_KinFit.sg_p7, par_downL_Gaussian_KinFit.sg_p7,
                            par_upTrig_Gaussian_KinFit.sg_p7, par_downTrig_Gaussian_KinFit.sg_p7};
    double sg_p7_errSyst_min=par_Gaussian_KinFit.sg_p7-(*std::min_element(sg_p7_errSyst, sg_p7_errSyst+11));
    double sg_p7_errSyst_max=(*std::max_element(sg_p7_errSyst, sg_p7_errSyst+11))-par_Gaussian_KinFit.sg_p7;
    outfile<<"   Uncertainty on sg_p7 = "<<par_Gaussian_KinFit.sg_p7<<" +- "<<sg_p7_errStat<<" (stat) - "<<sg_p7_errSyst_min<<" + "<<sg_p7_errSyst_max<<" (syst); -"<<quad(sg_p7_errStat/2., sg_p7_errSyst_min)<<"/+"<<quad(sg_p7_errStat/2., sg_p7_errSyst_max)<<" (total) <br/>"<<std::endl;
    double sg_p8_errStat=par_Gaussian_KinFit.sg_p8_err;
    double sg_p8_errSyst[]={par_Gaussian_KinFit.sg_p8,
                            par_JECp1_Gaussian_KinFit.sg_p8, par_JECm1_Gaussian_KinFit.sg_p8,
                            par_JERp1_Gaussian_KinFit.sg_p8, par_JERm1_Gaussian_KinFit.sg_p8,
                            par_upBC_Gaussian_KinFit.sg_p8, par_downBC_Gaussian_KinFit.sg_p8,
                            par_upL_Gaussian_KinFit.sg_p8, par_downL_Gaussian_KinFit.sg_p8,
                            par_upTrig_Gaussian_KinFit.sg_p8, par_downTrig_Gaussian_KinFit.sg_p8};
    double sg_p8_errSyst_min=par_Gaussian_KinFit.sg_p8-(*std::min_element(sg_p8_errSyst, sg_p8_errSyst+11));
    double sg_p8_errSyst_max=(*std::max_element(sg_p8_errSyst, sg_p8_errSyst+11))-par_Gaussian_KinFit.sg_p8;
    outfile<<"   Uncertainty on sg_p8 = "<<par_Gaussian_KinFit.sg_p8<<" +- "<<sg_p8_errStat<<" (stat) - "<<sg_p8_errSyst_min<<" + "<<sg_p8_errSyst_max<<" (syst); -"<<quad(sg_p8_errStat/2., sg_p8_errSyst_min)<<"/+"<<quad(sg_p8_errStat/2., sg_p8_errSyst_max)<<" (total) <br/>"<<std::endl;
    }
    outfile<<"  </td>"<<std::endl;
    
    outfile<<" </tr>"<<std::endl;
    outfile<<"</table>"<<std::endl;
    
  }
  outfile<<"</body>"<<std::endl;
  outfile<<"</html>"<<std::endl;
}
    

