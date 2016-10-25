#include <TROOT.h>
#include "TGraphAsymmErrors.h"


void CompareDifferentTTbar()
{

  TFile *f_ttbar_Madgraph=new TFile("InterpolDatacards_10GeV/UpperLimits_Interpol.root");
  TGraph *g_obs_Madgraph=(TGraph*)f_ttbar_Madgraph->Get("g_obs");
  TGraph *g_xsec=(TGraph*)f_ttbar_Madgraph->Get("g_xsec");
  TGraphAsymmErrors *g_xsec_1sigma=(TGraphAsymmErrors*)f_ttbar_Madgraph->Get("g_xsec_1sigma");
  TGraphAsymmErrors *g_xsec_2sigma=(TGraphAsymmErrors*)f_ttbar_Madgraph->Get("g_xsec_2sigma");
  
  TFile *f_ttbar_Herwig=new TFile("InterpolDatacards_ttbarHerwig/UpperLimits_Interpol.root");
  TGraph *g_obs_Herwig=(TGraph*)f_ttbar_Herwig->Get("g_obs");
  
  TCanvas *c_CompareDifferentTTbar=new TCanvas("c_CompareDifferentTTbar", "c_CompareDifferentTTbar", 700, 700);
  c_CompareDifferentTTbar->SetLogy();
  c_CompareDifferentTTbar->SetGridx(); c_CompareDifferentTTbar->SetGridy();
  g_xsec->SetMaximum(200); g_xsec->SetMinimum(7);
  g_xsec->Draw("AL*");
  g_xsec_2sigma->Draw("3");
  g_xsec_1sigma->Draw("3");
  g_xsec->Draw("L");
  g_obs_Madgraph->SetLineColor(kRed);
  g_obs_Madgraph->Draw("LP SAME");
  g_obs_Herwig->Draw("LP same");
  TLegend *leg=new TLegend(0.45, 0.6, 0.9, 0.85);
  leg->SetFillStyle(1); leg->SetFillColor(kWhite);
  leg->AddEntry(g_xsec, "Expected Upper Limit", "L");
  leg->AddEntry(g_xsec_1sigma, "Expected #pm 1 #sigma", "F");
  leg->AddEntry(g_xsec_2sigma, "Expected #pm 2 #sigma", "F");
  leg->AddEntry(g_obs_Madgraph, "Observed Upper Limit with MADGRAPH t#bar{t}", "LP");
  leg->AddEntry(g_obs_Herwig, "Observed Upper Limit with POWHEG t#bar{t}", "LP");
  leg->Draw();
  c_CompareDifferentTTbar->SaveAs("c_CompareDifferentTTbar.png");
  
}
