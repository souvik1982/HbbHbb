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
#include <TGraphAsymmErrors.h>

struct Samples
{
  std::vector<TH1F *> v_Cuts;
  std::vector<TH1F *> v_CountWithPU;
} signals;

void pushBackHistograms(Samples &sample, TFile *file, bool MC=true)
{
  if (MC==true)
  {
    sample.v_CountWithPU.push_back((TH1F*)file->Get("CountWithPU"));
  }
  sample.v_Cuts.push_back((TH1F*)file->Get("h_Cuts"));
}

void CutFlow()
{
  int massPts[15]={270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100};
  
  TFile *glugluToX270=new TFile("Histograms_RadionToHH_4b_M-270_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX300=new TFile("Histograms_RadionToHH_4b_M-300_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX350=new TFile("Histograms_RadionToHH_4b_M-350_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX400=new TFile("Histograms_RadionToHH_4b_M-400_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX450=new TFile("Histograms_RadionToHH_4b_M-450_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX500=new TFile("Histograms_RadionToHH_4b_M-500_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX550=new TFile("Histograms_RadionToHH_4b_M-550_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX600=new TFile("Histograms_RadionToHH_4b_M-600_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX650=new TFile("Histograms_RadionToHH_4b_M-650_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX700=new TFile("Histograms_RadionToHH_4b_M-700_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX800=new TFile("Histograms_RadionToHH_4b_M-800_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX900=new TFile("Histograms_RadionToHH_4b_M-900_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX1000=new TFile("Histograms_RadionToHH_4b_M-1000_TuneZ2star_8TeV_FULLSIM.root");
  TFile *glugluToX1100=new TFile("Histograms_RadionToHH_4b_M-1100_TuneZ2star_8TeV_FULLSIM.root");
  
  pushBackHistograms(signals, glugluToX270);
  pushBackHistograms(signals, glugluToX300);
  pushBackHistograms(signals, glugluToX350);
  pushBackHistograms(signals, glugluToX400);
  pushBackHistograms(signals, glugluToX450);
  pushBackHistograms(signals, glugluToX500);
  pushBackHistograms(signals, glugluToX550);
  pushBackHistograms(signals, glugluToX600);
  pushBackHistograms(signals, glugluToX650);
  pushBackHistograms(signals, glugluToX700);
  pushBackHistograms(signals, glugluToX800);
  pushBackHistograms(signals, glugluToX900);
  pushBackHistograms(signals, glugluToX1000);
  pushBackHistograms(signals, glugluToX1100);
  
  TH1F *h_Init=new TH1F("h_Init", "h_Init", 88, 320, 1200);
  std::vector<TH1F *> h_Cut;
  
  TH1F *h_Events2=(TH1F*)h_Init->Clone("h_Events2");
  TH1F *h_Events4=(TH1F*)h_Init->Clone("h_Events4");
  TH1F *h_Events6=(TH1F*)h_Init->Clone("h_Events6");
  TH1F *h_Events8=(TH1F*)h_Init->Clone("h_Events8");
  TH1F *h_Events10=(TH1F*)h_Init->Clone("h_Events10");
  TH1F *h_Events14=(TH1F*)h_Init->Clone("h_Events14");
  TH1F *h_Events16=(TH1F*)h_Init->Clone("h_Events16");
  for (unsigned int i=0; i<signals.v_Cuts.size(); ++i)
  {
    h_Init->SetBinContent(h_Init->FindBin(massPts[i]), signals.v_CountWithPU.at(i)->GetBinContent(1));
    // std::cout<<"signals.v_CountWithPU.at(i)->GetBinContent(1) = "<<.signals.v_CountWithPU.at(i)->GetBinContent(1)<<std::endl;
    // h_Events2->SetBinContent(i, signals.v_Cuts.at(i)->GetBinContent(2));
    h_Events4->SetBinContent(h_Init->FindBin(massPts[i]), signals.v_Cuts.at(i)->GetBinContent(4));
    h_Events6->SetBinContent(h_Init->FindBin(massPts[i]), signals.v_Cuts.at(i)->GetBinContent(6));
    h_Events8->SetBinContent(h_Init->FindBin(massPts[i]), signals.v_Cuts.at(i)->GetBinContent(8));
    h_Events10->SetBinContent(h_Init->FindBin(massPts[i]), signals.v_Cuts.at(i)->GetBinContent(10));
    h_Events14->SetBinContent(h_Init->FindBin(massPts[i]), signals.v_Cuts.at(i)->GetBinContent(14));
    h_Events16->SetBinContent(h_Init->FindBin(massPts[i]), signals.v_Cuts.at(i)->GetBinContent(16));
    
    std::cout<<"m_X (GeV) = "<<massPts[i]<<", Efficiency = "<<signals.v_Cuts.at(i)->GetBinContent(16)/signals.v_CountWithPU.at(i)->GetBinContent(1)<<std::endl;
  }
  
  // TGraphAsymmErrors *g_Ae_2=new TGraphAsymmErrors(h_Events2, h_Init); g_Ae_2->SetTitle("Step 2 ntuplization efficiency");
  TGraphAsymmErrors *g_Ae_4=new TGraphAsymmErrors(h_Events4, h_Init); g_Ae_4->SetTitle("Signal Acceptance #times Efficiency");
  TGraphAsymmErrors *g_Ae_6=new TGraphAsymmErrors(h_Events6, h_Init); g_Ae_6->SetTitle("Vtype efficiency");
  TGraphAsymmErrors *g_Ae_8=new TGraphAsymmErrors(h_Events8, h_Init); g_Ae_8->SetTitle("nCJets 3 efficiency");
  TGraphAsymmErrors *g_Ae_10=new TGraphAsymmErrors(h_Events10, h_Init); g_Ae_10->SetTitle("Signal Acceptance #times Efficiency; m_X (GeV); Cumulative Efficiency"); g_Ae_10->SetLineColor(kGreen);
  TGraphAsymmErrors *g_Ae_14=new TGraphAsymmErrors(h_Events14, h_Init); g_Ae_14->SetTitle("btagging efficiency"); g_Ae_14->SetLineColor(kBlue);
  TGraphAsymmErrors *g_Ae_16=new TGraphAsymmErrors(h_Events16, h_Init); g_Ae_16->SetTitle("SR efficiency"); g_Ae_16->SetLineColor(kRed);
  
  TCanvas *c_Ae=new TCanvas("c_Ae", "c_Ae", 1000, 700);
  g_Ae_10->SetMaximum(0.5); g_Ae_10->SetMinimum(0);
  // g_Ae_2->Draw("AL*");
  // g_Ae_4->Draw("AL*");
  // g_Ae_6->Draw("L* same");
  // g_Ae_8->Draw("L* same");
  g_Ae_10->Draw("AL* same");
  g_Ae_14->Draw("L* same");
  g_Ae_16->Draw("L* same");
  // c_Ae->SetLogy();
  TLegend *leg=new TLegend(0.5, 0.75, 0.9, 0.9);
  leg->SetFillStyle(1); leg->SetFillColor(kWhite);
  // leg->AddEntry(g_Ae_4, "Trigger");
  // leg->AddEntry(g_Ae_6, "Vtype");
  // leg->AddEntry(g_Ae_8, "Jets > 3");
  leg->AddEntry(g_Ae_10, "HH candidate");
  leg->AddEntry(g_Ae_14, "b-tagging");
  leg->AddEntry(g_Ae_16, "Signal Region");
  leg->Draw();
  c_Ae->Update();
  c_Ae->SaveAs("c_Ae.png");
}
