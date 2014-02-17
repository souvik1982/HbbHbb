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

std::string ftoa(double i) 
{
  char res[10];
  sprintf(res, "%f", i);
  std::string ret(res);
  return ret;
}

std::string itoa(int i) 
{
  char res[10];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

int DisplayDistributions_ttbar()
{

  double totalLuminosity=17928; // /pb
  
  double xsec_ttbar=225.197; // pb
  double br_ttbar_fulllept=(2.*0.108)*(2.*0.108);
  double br_ttbar_semilept=2*(2.*0.108)*0.67;
  double br_ttbar_hadronic=0.676*0.676;
  double xsec_ttbar_fulllept=xsec_ttbar*br_ttbar_fulllept;
  double xsec_ttbar_semilept=xsec_ttbar*br_ttbar_semilept;
  double xsec_ttbar_hadronic=xsec_ttbar*br_ttbar_hadronic;
  
  int rebin=1;
  
  xsec_ttbar_fulllept=24.56;
  xsec_ttbar_semilept=103.12;
  xsec_ttbar_hadronic=106.32;
  
  std::string tags="MMMM_nominal"; // MMMMbar, or MMMM_nominal
  std::string region="a_KinFit"; // a_KinFit, or b_KinFit
  
  TFile *ttbar_fulllept=new TFile((tags+"/"+region+"/Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root").c_str());
  TFile *ttbar_semilept=new TFile((tags+"/"+region+"/Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root").c_str());
  TFile *ttbar_hadronic=new TFile((tags+"/"+region+"/Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root").c_str());
  
  TFile *data=new TFile((tags+"/"+region+"/Histograms_BJetPlusX_Run2012BCD_Skim.root").c_str());
  
  double init_ttbar_fulllept=((TH1F*)ttbar_fulllept->Get("CountWithPU"))->GetBinContent(1);
  double init_ttbar_semilept=((TH1F*)ttbar_semilept->Get("CountWithPU"))->GetBinContent(1);
  double init_ttbar_hadronic=((TH1F*)ttbar_hadronic->Get("CountWithPU"))->GetBinContent(1);
  
  std::cout<<"init_ttbar_fulllept = "<<init_ttbar_fulllept<<std::endl;
  std::cout<<"init_ttbar_semilept = "<<init_ttbar_semilept<<std::endl;
  std::cout<<"init_ttbar_hadronic = "<<init_ttbar_hadronic<<std::endl;
  
  double scale_ttbar_fulllept=totalLuminosity*xsec_ttbar_fulllept/init_ttbar_fulllept;
  double scale_ttbar_semilept=totalLuminosity*xsec_ttbar_semilept/init_ttbar_semilept;
  double scale_ttbar_hadronic=totalLuminosity*xsec_ttbar_hadronic/init_ttbar_hadronic;
  
  std::cout<<"xsec_ttbar_fulllept = "<<xsec_ttbar_fulllept<<" pb"<<std::endl;
  std::cout<<"xsec_ttbar_semilept = "<<xsec_ttbar_semilept<<" pb"<<std::endl;
  std::cout<<"xsec_ttbar_hadronic = "<<xsec_ttbar_hadronic<<" pb"<<std::endl;
  
  std::cout<<"scale_ttbar_fulllept = "<<scale_ttbar_fulllept<<std::endl;
  std::cout<<"scale_ttbar_semilept = "<<scale_ttbar_semilept<<std::endl;
  std::cout<<"scale_ttbar_hadronic = "<<scale_ttbar_hadronic<<std::endl;
  
  TH1F *h_mX_SR_ttbar_fulllept=(TH1F*)ttbar_fulllept->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept=(TH1F*)ttbar_semilept->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic=(TH1F*)ttbar_hadronic->Get("h_mX_SR");
  
  TH1F *h_mX_CR2_ttbar_fulllept=(TH1F*)ttbar_fulllept->Get("h_mX_CR2");
  TH1F *h_mX_CR2_ttbar_semilept=(TH1F*)ttbar_semilept->Get("h_mX_CR2");
  TH1F *h_mX_CR2_ttbar_hadronic=(TH1F*)ttbar_hadronic->Get("h_mX_CR2");
  
  TH1F *h_mX_CR4_ttbar_fulllept=(TH1F*)ttbar_fulllept->Get("h_mX_CR4");
  TH1F *h_mX_CR4_ttbar_semilept=(TH1F*)ttbar_semilept->Get("h_mX_CR4");
  TH1F *h_mX_CR4_ttbar_hadronic=(TH1F*)ttbar_hadronic->Get("h_mX_CR4");
  
  TH1F *h_mX_SB_ttbar_fulllept=(TH1F*)h_mX_CR2_ttbar_fulllept->Clone("h_mX_SB"); h_mX_SB_ttbar_fulllept->Add(h_mX_CR4_ttbar_fulllept);
  TH1F *h_mX_SB_ttbar_semilept=(TH1F*)h_mX_CR2_ttbar_semilept->Clone("h_mX_SB"); h_mX_SB_ttbar_semilept->Add(h_mX_CR4_ttbar_semilept);
  TH1F *h_mX_SB_ttbar_hadronic=(TH1F*)h_mX_CR2_ttbar_hadronic->Clone("h_mX_SB"); h_mX_SB_ttbar_hadronic->Add(h_mX_CR4_ttbar_hadronic);
  
  // For data, if MMMM_nominal, take SB distribution and scale it to number of events in SR
  TH1F *h_mX_SR_data=(TH1F*)data->Get("h_mX_SR"); double n_SR=h_mX_SR_data->GetSumOfWeights();
  TH1F *h_mX_CR2_data=(TH1F*)data->Get("h_mX_CR2"); double n_CR2=h_mX_CR2_data->GetSumOfWeights();
  TH1F *h_mX_CR4_data=(TH1F*)data->Get("h_mX_CR4"); double n_CR4=h_mX_CR4_data->GetSumOfWeights();
  TH1F *h_mX_SB_data=(TH1F*)h_mX_CR2_data->Clone("h_mX_SB_data");
  h_mX_SB_data->Add(h_mX_CR4_data);
  if (tags=="MMMM_nominal")
  {
    h_mX_SR_data=(TH1F*)h_mX_SB_data->Clone("h_mX_SR_data");
    h_mX_SR_data->Scale(n_SR/(n_CR2+n_CR4));
  } 
  
  h_mX_SR_ttbar_fulllept->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic->Scale(scale_ttbar_hadronic);
  
  h_mX_SB_ttbar_fulllept->Scale(scale_ttbar_fulllept);
  h_mX_SB_ttbar_semilept->Scale(scale_ttbar_semilept);
  h_mX_SB_ttbar_hadronic->Scale(scale_ttbar_hadronic);
  
  h_mX_SR_ttbar_fulllept->Rebin(rebin);
  h_mX_SR_ttbar_semilept->Rebin(rebin);
  h_mX_SR_ttbar_hadronic->Rebin(rebin);
  h_mX_SB_ttbar_fulllept->Rebin(rebin);
  h_mX_SB_ttbar_semilept->Rebin(rebin);
  h_mX_SB_ttbar_hadronic->Rebin(rebin);
  h_mX_SB_data->Rebin(rebin);
  h_mX_SR_data->Rebin(rebin);
  
  h_mX_SR_ttbar_fulllept->SetFillColor(kRed);
  h_mX_SR_ttbar_semilept->SetFillColor(kBlue);
  h_mX_SR_ttbar_hadronic->SetFillColor(kGreen);
  
  h_mX_SB_ttbar_fulllept->SetFillColor(kRed);
  h_mX_SB_ttbar_semilept->SetFillColor(kBlue);
  h_mX_SB_ttbar_hadronic->SetFillColor(kGreen);
  
  h_mX_SB_data->GetXaxis()->SetRangeUser(200, 1200);
  h_mX_SR_data->GetXaxis()->SetRangeUser(200, 1200);
  if (region=="a_KinFit") 
  {
    h_mX_SB_data->SetTitle(("t#bar{t} Composition in SB; m_{X} GeV; Events / "+itoa(10*rebin)+" GeV").c_str());
    h_mX_SR_data->SetTitle(("t#bar{t} Composition of Data-driven Background Estimate in SR; m_{X} GeV; Events / "+itoa(10*rebin)+" GeV").c_str());
  }
  else if (region=="b_KinFit") 
  {
    h_mX_SB_data->SetTitle(("t#bar{t} Composition in VR-SB; m_{X} GeV; Events / "+itoa(10*rebin)+" GeV").c_str());
    h_mX_SR_data->SetTitle(("t#bar{t} Composition of VR; m_{X} GeV; Events / "+itoa(10*rebin)+" GeV").c_str());
  }
  h_mX_SB_data->SetLineColor(kBlack);
  h_mX_SR_data->SetLineColor(kBlack);
  h_mX_SB_data->SetLineWidth(2);
  h_mX_SR_data->SetLineWidth(2);
  
  THStack *s_mX_SR_ttbar=new THStack("s_mX_SR_ttbar", "s_mX_SR_ttbar");
  s_mX_SR_ttbar->Add(h_mX_SR_ttbar_hadronic, "hist");
  s_mX_SR_ttbar->Add(h_mX_SR_ttbar_semilept, "hist");
  s_mX_SR_ttbar->Add(h_mX_SR_ttbar_fulllept, "hist");
  
  THStack *s_mX_SB_ttbar=new THStack("s_mX_SB_ttbar", "s_mX_SB_ttbar");
  s_mX_SB_ttbar->Add(h_mX_SB_ttbar_hadronic, "hist");
  s_mX_SB_ttbar->Add(h_mX_SB_ttbar_semilept, "hist");
  s_mX_SB_ttbar->Add(h_mX_SB_ttbar_fulllept, "hist");
  
  // Subtracted data histograms
  TH1F *h_mX_SR_data_minusttbar=(TH1F*)h_mX_SR_data->Clone("h_mX_SR_data_minusttbar");
  h_mX_SR_data_minusttbar->Add(h_mX_SR_ttbar_hadronic, -1);
  h_mX_SR_data_minusttbar->Add(h_mX_SR_ttbar_semilept, -1);
  h_mX_SR_data_minusttbar->Add(h_mX_SR_ttbar_fulllept, -1);
  
  TH1F *h_mX_SB_data_minusttbar=(TH1F*)h_mX_SB_data->Clone("h_mX_SB_data_minusttbar");
  h_mX_SB_data_minusttbar->Add(h_mX_SB_ttbar_hadronic, -1);
  h_mX_SB_data_minusttbar->Add(h_mX_SB_ttbar_semilept, -1);
  h_mX_SB_data_minusttbar->Add(h_mX_SB_ttbar_fulllept, -1);
  
  h_mX_SR_data_minusttbar->SetLineColor(kBlue);
  h_mX_SB_data_minusttbar->SetLineColor(kBlue);
  
  TLegend *leg=new TLegend(0.5, 0.75, 0.9, 0.9);
  if (tags=="MMMM_nominal") leg->AddEntry(h_mX_SR_data, "Data from SB scaled to 17.928 /fb");
  else leg->AddEntry(h_mX_SR_data, "Data from SR in 17.928 /fb");
  leg->AddEntry(h_mX_SR_ttbar_fulllept, "Leptonic t#bar{t}");
  leg->AddEntry(h_mX_SR_ttbar_semilept, "Semi-leptonic t#bar{t}");
  leg->AddEntry(h_mX_SR_ttbar_hadronic, "Hadronic t#bar{t}");
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  
  TCanvas *c_mX_SR_ttbar=new TCanvas("c_mX_SR_ttbar", "c_mX_SR_ttbar", 700, 700);
  h_mX_SR_data->Draw("ep9");
  h_mX_SR_data_minusttbar->Draw("ep7 same");
  s_mX_SR_ttbar->Draw("same");
  double ks_SR=h_mX_SR_data->KolmogorovTest(h_mX_SR_data_minusttbar);
  std::cout<<"ks_SR = "<<ks_SR<<std::endl;
  TLegend *leg_SR=new TLegend(0.6, 0.55, 0.9, 0.6);
  leg_SR->AddEntry((TObject*)0, ("KS(data, data-ttbar) = "+ftoa(ks_SR)).c_str(), "");
  leg->Draw();
  leg_SR->Draw();
  s_mX_SR_ttbar->GetHistogram()->SetTitle("mX distribution in Central Region; m_{X} {GeV); Events / 40 GeV");
  if (tags=="MMMM_nominal") c_mX_SR_ttbar->SetLogy();
  c_mX_SR_ttbar->SaveAs(("c_mX_"+tags+"_"+region+"_SR_ttbar.png").c_str());
  
  TCanvas *c_mX_SB_ttbar=new TCanvas("c_mX_SB_ttbar", "c_mX_SB_ttbar", 700, 700);
  h_mX_SB_data->Draw("ep9");
  h_mX_SB_data_minusttbar->Draw("ep7 same");
  s_mX_SB_ttbar->Draw("same");
  double ks_SB=h_mX_SB_data->KolmogorovTest(h_mX_SB_data_minusttbar);
  std::cout<<"ks_SB = "<<ks_SB<<std::endl;
  TLegend *leg_SB=new TLegend(0.6, 0.55, 0.9, 0.6);
  leg_SB->AddEntry((TObject*)0, ("KS(data, data-ttbar) = "+ftoa(ks_SB)).c_str(), "");
  leg->Draw();
  leg_SB->Draw();
  s_mX_SB_ttbar->GetHistogram()->SetTitle("mX distribution in Sideband Region; m_{X} {GeV); Events / 40 GeV");
  if (tags=="MMMM_nominal") c_mX_SB_ttbar->SetLogy();
  c_mX_SB_ttbar->SaveAs(("c_mX_"+tags+"_"+region+"_SB_ttbar.png").c_str());
  
  double contrib_fulllept=h_mX_SR_ttbar_fulllept->GetSumOfWeights();
  double contrib_semilept=h_mX_SR_ttbar_semilept->GetSumOfWeights();
  double contrib_hadronic=h_mX_SR_ttbar_hadronic->GetSumOfWeights();
  double contrib_data=h_mX_SR_data->GetSumOfWeights();
  
  std::cout<<"Contribution of fully leptonic ttbar is "<<contrib_fulllept<<std::endl;
  std::cout<<"Contribution of semi leptonic ttbar is "<<contrib_semilept<<std::endl;
  std::cout<<"Contribution of hadronic ttbar is "<<contrib_hadronic<<std::endl;
  std::cout<<"Total contribution of ttbar is "<<contrib_fulllept+contrib_semilept+contrib_hadronic<<std::endl;
  std::cout<<"Number of events in data is "<<contrib_data<<std::endl; 
  
  return 1;
}
