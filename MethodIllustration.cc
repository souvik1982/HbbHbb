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
#include <TArc.h>

double mH_mean_cut=20.;
double marg=0.9*mH_mean_cut;

bool signal=false;
bool data=true;
bool ttbar=true;

double Normalize(TH1* h)
{
  double nEntries=h->GetSumOfWeights();
  h->Scale(1./nEntries);
  return nEntries;
}

void drawRegion(double H_mass, double r1, double r2, std::string sr, std::string cr, bool blind=false)
{
  TArc *circle1=new TArc(H_mass, H_mass, r1); circle1->SetLineWidth(2); circle1->SetLineColor(kBlue); circle1->SetFillStyle(0); 
  if (blind==true) 
  {
    circle1->SetFillStyle(1001);
    circle1->SetFillColor(kBlack);
  } 
  circle1->Draw();
  TArc *circle2=new TArc(H_mass, H_mass, r2, 90., 180.); circle2->SetLineWidth(2); circle2->SetNoEdges(); circle2->SetLineColor(kBlack); circle2->SetFillStyle(0); circle2->Draw();
  TArc *circle3=new TArc(H_mass, H_mass, r2, 270., 360.); circle3->SetLineWidth(3); circle3->SetNoEdges(); circle3->SetLineColor(kBlack); circle3->SetFillStyle(0); circle3->Draw();
  TLine *line1=new TLine(H_mass-r2, H_mass, H_mass-r1, H_mass); line1->SetLineWidth(2); line1->SetLineColor(kBlack); line1->Draw();
  TLine *line2=new TLine(H_mass+r2, H_mass, H_mass+r1, H_mass); line2->SetLineWidth(2); line2->SetLineColor(kBlack); line2->Draw();
  TLine *line3=new TLine(H_mass, H_mass-r2, H_mass, H_mass-r1); line3->SetLineWidth(2); line3->SetLineColor(kBlack); line3->Draw();
  TLine *line4=new TLine(H_mass, H_mass+r2, H_mass, H_mass+r1); line4->SetLineWidth(2); line4->SetLineColor(kBlack); line4->Draw();
  TArrow *arrow1=new TArrow(H_mass, H_mass+r2*5, H_mass, H_mass, 0.02); arrow1->SetLineWidth(2); arrow1->SetLineColor(kBlack); arrow1->Draw();
  TPaveText *mod1=new TPaveText(H_mass-marg, H_mass+r2*5-marg, H_mass+marg, H_mass+r2*5+marg);
  mod1->SetBorderSize(0); mod1->SetFillColor(0); mod1->AddText(sr.c_str()); mod1->SetLineColor(kBlue); mod1->Draw("ARC");
  TArrow *arrow2_1=new TArrow(H_mass+r2*4., H_mass, H_mass-r2/2., H_mass+r2/2., 0.02); arrow2_1->SetLineWidth(2); arrow2_1->SetLineColor(kBlack);     
  TArrow *arrow2_2=new TArrow(H_mass+r2*4., H_mass, H_mass+r2/2., H_mass-r2/2., 0.02); arrow2_2->SetLineWidth(2); arrow2_2->SetLineColor(kBlack);
  TLine *arrow2_3=new TLine(H_mass+r2*4., H_mass, H_mass+r2*5., H_mass); arrow2_3->SetLineWidth(2); arrow2_3->SetLineColor(kBlack);
  arrow2_1->Draw(); arrow2_2->Draw(); arrow2_3->Draw();
  TPaveText *mod2=new TPaveText(H_mass+r2*5.-marg, H_mass+marg, H_mass+r2*5.+marg, H_mass-marg);
  mod2->SetBorderSize(0); mod2->SetFillColor(0); mod2->AddText(cr.c_str()); mod2->SetLineColor(kBlack); mod2->Draw("ARC");
  
}

void MethodIllustration()
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  
  TLegend *leg=new TLegend(0.45, 0.75, 0.9, 0.9);
  leg->SetFillStyle(0);
  leg->AddEntry((TObject*)0, "SR: Signal Region", "");
  leg->AddEntry((TObject*)0, "SB: Sideband", "");
  leg->AddEntry((TObject*)0, "VR: Control Region ", "");
  leg->AddEntry((TObject*)0, "VR-SB: Control Region  Sideband", "");
  // leg->AddEntry((TObject*)0, "CR2: Control Region 2", "");
  // leg->AddEntry((TObject*)0, "CR2-SB: Control Region 2 Sideband", ""); 
  
  if (data)
  {
    TFile *data_8TeVData2012=new TFile("Histograms_BJetPlusX_Run2012BCD_Skim.root");
    TH2F *h_mH_mH_3Tag_data=(TH2F*)data_8TeVData2012->Get("h_mH_mH");
    TCanvas *c_mH_mH_3Tag=new TCanvas("c_mH_mH_3Tag", "c_mH_mH_3Tag", 700, 700);
    h_mH_mH_3Tag_data->SetTitle("m_{H1} vs m_{H2} for Data; m_{H1} [GeV]; m_{H2} [GeV]");
    h_mH_mH_3Tag_data->Draw("colz");
    drawRegion(125., 17.5, 35., "SR", "SB", true);
    drawRegion(90., 17.5, 35., "VR", "VB");
    leg->Draw();
    c_mH_mH_3Tag->SaveAs("Illustration_nTag_data.png");
  }
  
  if (signal)
  {
    TFile *signal_8TeVData2012=new TFile("Histograms_RadionToHH_4b_M-700_TuneZ2star_8TeV_FULLSIM.root");
    TH2F *h_mH_mH_3Tag_signal=(TH2F*)signal_8TeVData2012->Get("h_mH_mH");
    TCanvas *c_mH_mH_3Tag_signal=new TCanvas("c_mH_mH_3Tag_signal", "c_mH_mH_3Tag_signal", 700, 700);
    h_mH_mH_3Tag_signal->SetTitle("mH1 vs mH2 for Signal mX=600 GeV; mH1 (GeV); mH2 (GeV)");
    h_mH_mH_3Tag_signal->Draw("colz");
    drawRegion(125., 17.5, 35., "SR", "SB");
    drawRegion(90., 17.5, 35., "VR", "VB");
    leg->Draw();
    c_mH_mH_3Tag_signal->SaveAs("Illustration_nTag_signal.png");
  }
  
  if (ttbar)
  {
    TFile *ttbar=new TFile("Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
    TH2F *h_mH_mH_ttbar=(TH2F*)ttbar->Get("h_mH_mH");
    TCanvas *c_mH_mH_ttbar=new TCanvas("c_mH_mH_ttbar", "c_mH_mH_ttbar", 700, 700);
    h_mH_mH_ttbar->SetTitle("mH1 vs mH2 for hadronic ttbar; m_{H1} [GeV]; m_{H2} [GeV]");
    h_mH_mH_ttbar->Draw("colz");
    drawRegion(125., 17.5, 35., "SR", "SB");
    drawRegion(90., 17.5, 35., "VR", "VB");
    leg->Draw();
    c_mH_mH_ttbar->SaveAs("Illustration_nTag_ttbar.png");
  }
}
    
  
  
