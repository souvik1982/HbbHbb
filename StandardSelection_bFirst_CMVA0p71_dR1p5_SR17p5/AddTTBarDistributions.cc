#include <TFile.h>
#include <TH1F.h>

double totalLuminosity=17928; // /pb
double xsec_ttbar_fulllept=24.56;
double xsec_ttbar_semilept=103.12;
double xsec_ttbar_hadronic=106.32;

void AddTTBarDistributions()
{

  TFile *ttbar_fulllept=new TFile("Histograms_TTJets_FullLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_semilept=new TFile("Histograms_TTJets_SemiLeptMGDecays_8TeV-madgraph_Skim.root");
  TFile *ttbar_hadronic=new TFile("Histograms_TTJets_HadronicMGDecays_8TeV-madgraph_Skim.root");
  
  double init_ttbar_fulllept=((TH1F*)ttbar_fulllept->Get("CountWithPU"))->GetBinContent(1);
  double init_ttbar_semilept=((TH1F*)ttbar_semilept->Get("CountWithPU"))->GetBinContent(1);
  double init_ttbar_hadronic=((TH1F*)ttbar_hadronic->Get("CountWithPU"))->GetBinContent(1);
  
  double scale_ttbar_fulllept=totalLuminosity*xsec_ttbar_fulllept/init_ttbar_fulllept;
  double scale_ttbar_semilept=totalLuminosity*xsec_ttbar_semilept/init_ttbar_semilept;
  double scale_ttbar_hadronic=totalLuminosity*xsec_ttbar_hadronic/init_ttbar_hadronic;
  
  TH1F *h_mX_SR_ttbar_fulllept=(TH1F*)ttbar_fulllept->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_semilept=(TH1F*)ttbar_semilept->Get("h_mX_SR");
  TH1F *h_mX_SR_ttbar_hadronic=(TH1F*)ttbar_hadronic->Get("h_mX_SR");
  
  TH1F *h_mX_SB_ttbar_fulllept=(TH1F*)ttbar_fulllept->Get("h_mX_CR2"); h_mX_SB_ttbar_fulllept->Add((TH1F*)ttbar_fulllept->Get("h_mX_CR4"));
  TH1F *h_mX_SB_ttbar_semilept=(TH1F*)ttbar_semilept->Get("h_mX_CR2"); h_mX_SB_ttbar_semilept->Add((TH1F*)ttbar_semilept->Get("h_mX_CR4"));
  TH1F *h_mX_SB_ttbar_hadronic=(TH1F*)ttbar_hadronic->Get("h_mX_CR2"); h_mX_SB_ttbar_hadronic->Add((TH1F*)ttbar_hadronic->Get("h_mX_CR4"));
  
  h_mX_SR_ttbar_fulllept->Scale(scale_ttbar_fulllept);
  h_mX_SR_ttbar_semilept->Scale(scale_ttbar_semilept);
  h_mX_SR_ttbar_hadronic->Scale(scale_ttbar_hadronic);
  TH1F *h_mX_SR_ttbar=(TH1F*)h_mX_SR_ttbar_fulllept->Clone("h_mX_SR_ttbar");
  h_mX_SR_ttbar->Add(h_mX_SR_ttbar_semilept);
  h_mX_SR_ttbar->Add(h_mX_SR_ttbar_hadronic);
  h_mX_SR_ttbar->SetTitle("h_mX_SR_ttbar");
  
  h_mX_SB_ttbar_fulllept->Scale(scale_ttbar_fulllept);
  h_mX_SB_ttbar_semilept->Scale(scale_ttbar_semilept);
  h_mX_SB_ttbar_hadronic->Scale(scale_ttbar_hadronic);
  TH1F *h_mX_SB_ttbar=(TH1F*)h_mX_SB_ttbar_fulllept->Clone("h_mX_SB_ttbar");
  h_mX_SB_ttbar->Add(h_mX_SB_ttbar_semilept);
  h_mX_SB_ttbar->Add(h_mX_SB_ttbar_hadronic);
  h_mX_SB_ttbar->SetTitle("h_mX_SB_ttbar");
  
  TFile *outfile=new TFile("Histograms_TTJets_mX.root", "recreate");
  outfile->cd();
  h_mX_SR_ttbar->Write();
  h_mX_SB_ttbar->Write();
  outfile->Close();
}
