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

bool all=true;

bool signal_sets=true;
bool data_sets=true;
bool qcd_sets=false;

// QCD scale factors
double totalLuminosity=18600; // /pb
double nSignal_init=101421.;
double prodXsec_1=1.; // pb
/*
double xsec_QCDHT250To500=276000; // pb
double xsec_QCDHT500To1000=8426; // pb
double xsec_QCDHT1000ToInf=204; // pb
double init_QCDHT250To500=4669914;
double init_QCDHT500To1000=4779910;
double init_QCDHT1000ToInf=5192110;
double scale_QCDHT250To500=totalLuminosity*xsec_QCDHT250To500/init_QCDHT250To500;
double scale_QCDHT500To1000=totalLuminosity*xsec_QCDHT500To1000/init_QCDHT500To1000;
double scale_QCDHT1000ToInf=totalLuminosity*xsec_QCDHT1000ToInf/init_QCDHT1000ToInf;
*/
double xsec_QCDpT50To80=8148778.0;
double xsec_QCDpT80To120=1033680.0;
double xsec_QCDpT120To170=156293.3;
double xsec_QCDpT170To300=34138.15;
double xsec_QCDpT300To470=1759.549;
double xsec_QCDpT470To600=113.8791;
double xsec_QCDpT600To800=26.9921;
double xsec_QCDpT800To1000=3.550036;
double xsec_QCDpT1000To1400=0.737844;
double xsec_QCDpT1400To1800=0.03352235;
double xsec_QCDpT1800=0.001829005;
double init_QCDpT50To80=5788859;
double init_QCDpT80To120=5964863;
double init_QCDpT120To170=5745731;
double init_QCDpT170To300=5394397;
double init_QCDpT300To470=5558499;
double init_QCDpT470To600=2614847;
double init_QCDpT600To800=3966863;
double init_QCDpT800To1000=3968562;
double init_QCDpT1000To1400=1964087;
double init_QCDpT1400To1800=2000061;
double init_QCDpT1800=977585;
double scale_QCDpT50To80=totalLuminosity*xsec_QCDpT50To80/init_QCDpT50To80;
double scale_QCDpT80To120=totalLuminosity*xsec_QCDpT80To120/init_QCDpT80To120;
double scale_QCDpT120To170=totalLuminosity*xsec_QCDpT120To170/init_QCDpT120To170;
double scale_QCDpT170To300=totalLuminosity*xsec_QCDpT170To300/init_QCDpT170To300;
double scale_QCDpT300To470=totalLuminosity*xsec_QCDpT300To470/init_QCDpT300To470;
double scale_QCDpT470To600=totalLuminosity*xsec_QCDpT470To600/init_QCDpT470To600;
double scale_QCDpT600To800=totalLuminosity*xsec_QCDpT600To800/init_QCDpT600To800;
double scale_QCDpT800To1000=totalLuminosity*xsec_QCDpT800To1000/init_QCDpT800To1000;
double scale_QCDpT1000To1400=totalLuminosity*xsec_QCDpT1000To1400/init_QCDpT1000To1400;
double scale_QCDpT1400To1800=totalLuminosity*xsec_QCDpT1400To1800/init_QCDpT1400To1800;
double scale_QCDpT1800=totalLuminosity*xsec_QCDpT1800/init_QCDpT1800;

struct Samples
{
  std::vector<TH1F *> v_nJets;
  std::vector<TH1F *> v_nPV;
  std::vector<TH1F *> v_CSV1;
  std::vector<TH1F *> v_CSV2;
  std::vector<TH1F *> v_CSV3;
  std::vector<TH1F *> v_CSV4;
  std::vector<TH1F *> v_JetpT1;
  std::vector<TH1F *> v_JetpT2;
  std::vector<TH1F *> v_JetpT3;
  std::vector<TH1F *> v_JetpT4;
  std::vector<TH1F *> v_H1_mass;
  std::vector<TH1F *> v_H1_pT;
  std::vector<TH1F *> v_H2_mass;
  std::vector<TH1F *> v_H2_pT;
  std::vector<TH1F *> v_H1_CSV1;
  std::vector<TH1F *> v_H1_CSV2;
  std::vector<TH1F *> v_H2_CSV1;
  std::vector<TH1F *> v_H2_CSV2;
  std::vector<TH1F *> v_nTags;
  std::vector<TH1F *> v_X_mass;
  std::vector<TH1F *> v_X_pT;
  std::vector<TH1F *> v_HH_mass_diff;
  std::vector<TH1F *> v_HH_mass_mean;
  std::vector<TH1F *> v_HH_dPhi;
  std::vector<TH1F *> v_HH_deta;
  std::vector<TH1F *> v_Azimuth;
  std::vector<TH1F *> v_HH_mass_mean_tagged;
  std::vector<TH1F *> v_MET;
  std::vector<TH1F *> v_MET_sig;
  std::vector<TH1F *> v_HT;
  std::vector<TH1F *> v_genX_mass;
  std::vector<TH1F *> v_genH1_mass;
  std::vector<TH1F *> v_genH2_mass;
  std::vector<TH1F *> v_nTags_SB;
  std::vector<TH1F *> v_nTags_effSig;
  std::vector<TH1F *> v_Cuts;
  std::vector<TH1F *> v_H1_mass_bTagged;
  std::vector<TH1F *> v_H2_mass_bTagged;
  std::vector<TH1F *> v_H1_pT_bTagged;
  std::vector<TH1F *> v_H2_pT_bTagged;
  std::vector<TH1F *> v_X_mass_bTagged;
  std::vector<TH1F *> v_X_pT_bTagged;
  std::vector<TH1F *> v_mX_CR1;
  std::vector<TH1F *> v_mX_CR2;
  std::vector<TH1F *> v_mX_CR3;
  std::vector<TH1F *> v_mX_CR4;
  std::vector<TH1F *> v_mX_CR5;
  std::vector<TH1F *> v_mX_SR;
  std::vector<TH1F *> v_mX_CR24;
  std::vector<TH1F *> v_CountWithPU;
  std::vector<TH1F *> v_dR1;
  std::vector<TH1F *> v_dR1_KinFit;
} signals, QCDs, data;

std::string itoa(int i) 
{
  char res[4];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

double Normalize(TH1* h)
{
  double nEntries=h->GetSumOfWeights();
  h->Scale(1./nEntries);
  return nEntries;
}

void pushBackHistograms(Samples &sample, TFile *file, bool MC=true)
{
  sample.v_nJets.push_back((TH1F*)file->Get("h_nJets"));
  sample.v_nPV.push_back((TH1F*)file->Get("h_nPV"));
  sample.v_CSV1.push_back((TH1F*)file->Get("h_CSV1"));
  sample.v_CSV2.push_back((TH1F*)file->Get("h_CSV2"));
  sample.v_CSV3.push_back((TH1F*)file->Get("h_CSV3"));
  sample.v_CSV4.push_back((TH1F*)file->Get("h_CSV4"));
  sample.v_JetpT1.push_back((TH1F*)file->Get("h_JetpT1"));
  sample.v_JetpT2.push_back((TH1F*)file->Get("h_JetpT2"));
  sample.v_JetpT3.push_back((TH1F*)file->Get("h_JetpT3"));
  sample.v_JetpT4.push_back((TH1F*)file->Get("h_JetpT4"));
  sample.v_H1_mass.push_back((TH1F*)file->Get("h_H1_mass"));
  sample.v_H1_pT.push_back((TH1F*)file->Get("h_H1_pT"));
  sample.v_H2_mass.push_back((TH1F*)file->Get("h_H2_mass"));
  sample.v_H2_pT.push_back((TH1F*)file->Get("h_H2_pT"));
  sample.v_H1_CSV1.push_back((TH1F*)file->Get("h_H1_CSV1"));
  sample.v_H1_CSV2.push_back((TH1F*)file->Get("h_H1_CSV2"));
  sample.v_H2_CSV1.push_back((TH1F*)file->Get("h_H2_CSV1"));
  sample.v_H2_CSV2.push_back((TH1F*)file->Get("h_H2_CSV2"));
  sample.v_nTags.push_back((TH1F*)file->Get("h_nTags"));
  sample.v_X_mass.push_back((TH1F*)file->Get("h_X_mass"));
  sample.v_X_pT.push_back((TH1F*)file->Get("h_X_pT"));
  sample.v_HH_mass_diff.push_back((TH1F*)file->Get("h_HH_mass_diff_cand"));
  sample.v_HH_mass_mean.push_back((TH1F*)file->Get("h_HH_mass_mean"));
  sample.v_HH_dPhi.push_back((TH1F*)file->Get("h_HH_dPhi"));
  sample.v_HH_deta.push_back((TH1F*)file->Get("h_HH_deta"));
  sample.v_Azimuth.push_back((TH1F*)file->Get("h_Azimuth"));
  sample.v_HH_mass_mean_tagged.push_back((TH1F*)file->Get("h_HH_mass_mean_tagged"));
  sample.v_MET.push_back((TH1F*)file->Get("h_MET"));
  sample.v_MET_sig.push_back((TH1F*)file->Get("h_MET_sig"));
  sample.v_HT.push_back((TH1F*)file->Get("h_HT"));
  if (MC==true)
  {
    sample.v_genX_mass.push_back((TH1F*)file->Get("h_genX_mass"));
    sample.v_genH1_mass.push_back((TH1F*)file->Get("h_genH1_mass"));
    sample.v_genH2_mass.push_back((TH1F*)file->Get("h_genH2_mass"));
    sample.v_CountWithPU.push_back((TH1F*)file->Get("CountWithPU"));
  }
  sample.v_Cuts.push_back((TH1F*)file->Get("h_Cuts"));
  sample.v_H1_mass_bTagged.push_back((TH1F*)file->Get("h_H1_mass_bTagged"));
  sample.v_H2_mass_bTagged.push_back((TH1F*)file->Get("h_H2_mass_bTagged"));
  sample.v_H1_pT_bTagged.push_back((TH1F*)file->Get("h_H1_pT_bTagged"));
  sample.v_H2_pT_bTagged.push_back((TH1F*)file->Get("h_H2_pT_bTagged"));
  sample.v_X_mass_bTagged.push_back((TH1F*)file->Get("h_X_mass_bTagged"));
  sample.v_X_pT_bTagged.push_back((TH1F*)file->Get("h_X_pT_bTagged"));
  sample.v_mX_CR1.push_back((TH1F*)file->Get("h_mX_CR1"));
  sample.v_mX_CR2.push_back((TH1F*)file->Get("h_mX_CR2"));
  sample.v_mX_CR3.push_back((TH1F*)file->Get("h_mX_CR3"));
  sample.v_mX_CR4.push_back((TH1F*)file->Get("h_mX_CR4"));
  sample.v_mX_CR5.push_back((TH1F*)file->Get("h_mX_CR5"));
  sample.v_mX_SR.push_back((TH1F*)file->Get("h_mX_SR"));
  TH1F *h_mX_CR24=(TH1F*)file->Get("h_mX_CR2")->Clone("h_mX_CR24");
  h_mX_CR24->SetTitle("m_{X} in CR24");
  h_mX_CR24->Add((TH1F*)file->Get("h_mX_CR4"));
  sample.v_mX_CR24.push_back(h_mX_CR24);
  sample.v_dR1.push_back((TH1F*)file->Get("h_dR1"));
  sample.v_dR1_KinFit.push_back((TH1F*)file->Get("h_dR1_KinFit"));
  
  TH1F *h_nTags_SB=(TH1F*)((file->Get("h_nTags"))->Clone("h_nTags_SB")); h_nTags_SB->SetTitle("S/B function of nTags; nTags");
  TH1F *h_nTags_effSig=(TH1F*)((file->Get("h_nTags"))->Clone("h_nTags_effSig")); h_nTags_effSig->SetTitle("eff(S)/sqrt(eff(B)) function of nTags; nTags");
  sample.v_nTags_SB.push_back(h_nTags_SB);
  sample.v_nTags_effSig.push_back(h_nTags_effSig);
}

void EfficiencyTable(std::vector<TH1F *> v_Cuts, std::vector<TH1F *> v_Cuts_data)
{
  int massPts[15]={270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100};
  
  v_Cuts.at(3)->Draw();
  
  std::cout<<"EFFICIENCY RAW NUMBERS"<<std::endl;
  std::cout<<" === After Step 2 ==="<<std::endl;
  for (unsigned int i=0; i<v_Cuts.size(); ++i)
  {
    std::cout<<"--- mX = "<<massPts[i]<<" --- "<<std::endl;
    std::cout<<"#events after Step 2 = "<<v_Cuts.at(i)->GetBinContent(2)<<std::endl;
    std::cout<<"#events after trigger = "<<v_Cuts.at(i)->GetBinContent(4)<<std::endl;
    std::cout<<"#events after Vtype = "<<v_Cuts.at(i)->GetBinContent(6)<<std::endl;
    std::cout<<"#events after nCJets > 3 = "<<v_Cuts.at(i)->GetBinContent(8)<<std::endl;
    std::cout<<"#events after HH candidate = "<<v_Cuts.at(i)->GetBinContent(10)<<std::endl;
    std::cout<<"#events after bTagging = "<<v_Cuts.at(i)->GetBinContent(14)<<std::endl;
    std::cout<<"#events after SR = "<<v_Cuts.at(i)->GetBinContent(16)<<std::endl; 
  }
  std::cout<<"======================="<<std::endl;
  
  std::cout<<" === Step 2 Efficiencies wrt production ==="<<std::endl;
  for (unsigned int i=0; i<v_Cuts.size(); ++i)
  {
    double step2Eff=(v_Cuts.at(i)->GetBinContent(2))/(1e5);
    std::cout<<"mX = "<<massPts[i]<<" has Step 2 efficiency wrt production = "<<step2Eff*100.<<"%"<<std::endl;
  }
  std::cout<<"Data has events after Step 2 = "<<v_Cuts_data.at(0)->GetBinContent(2)<<std::endl;
  
  std::cout<<" === Trigger Efficiencies wrt Step 2 ==="<<std::endl;
  for (unsigned int i=0; i<v_Cuts.size(); ++i)
  {
    double triggerEff=(v_Cuts.at(i)->GetBinContent(4))/(v_Cuts.at(i)->GetBinContent(2));
    std::cout<<"mX = "<<massPts[i]<<" has trigger efficiency wrt step 2 = "<<triggerEff*100.<<"%"<<std::endl;
  }
  std::cout<<"Data has events after Trigger = "<<v_Cuts_data.at(0)->GetBinContent(4)<<std::endl;
  
  std::cout<<" === Vtype == Efficiency === "<<std::endl;
  for (unsigned int i=0; i<v_Cuts.size(); ++i)
  {
    double vtypeEff=(v_Cuts.at(i)->GetBinContent(6))/(v_Cuts.at(i)->GetBinContent(4));
    std::cout<<"mX = "<<massPts[i]<<" has vType efficiency = "<<vtypeEff*100.<<"%"<<std::endl;
  }
  std::cout<<"Data has events after vType = "<<v_Cuts_data.at(0)->GetBinContent(6)<<std::endl;
  
  std::cout<<" === nCJets > 3 Efficiency === "<<std::endl;
  for (unsigned int i=0; i<v_Cuts.size(); ++i)
  {
    double nCJetsEff=(v_Cuts.at(i)->GetBinContent(8))/(v_Cuts.at(i)->GetBinContent(6));
    std::cout<<"mX = "<<massPts[i]<<" has nCJets > 3 efficiency = "<<nCJetsEff*100.<<"%"<<std::endl;
  }
  std::cout<<"Data has events after nCJets = "<<v_Cuts_data.at(0)->GetBinContent(8)<<std::endl;
  
  std::cout<<" === HH Cand Efficiency === "<<std::endl;
  for (unsigned int i=0; i<v_Cuts.size(); ++i)
  {
    double hhCandEff=(v_Cuts.at(i)->GetBinContent(10))/(v_Cuts.at(i)->GetBinContent(8));
    std::cout<<"mX = "<<massPts[i]<<" has HH Cand efficiency = "<<hhCandEff*100.<<"%"<<std::endl;
  }
  std::cout<<"Data has events after HH Cand = "<<v_Cuts_data.at(0)->GetBinContent(10)<<std::endl;
  
  std::cout<<" === b-tagging Efficiency === "<<std::endl;
  for (unsigned int i=0; i<v_Cuts.size(); ++i)
  {
    double btagEff=(v_Cuts.at(i)->GetBinContent(14))/(v_Cuts.at(i)->GetBinContent(10));
    std::cout<<"mX = "<<massPts[i]<<" has b-tagging efficiency = "<<btagEff*100.<<"%"<<std::endl;
  }
  std::cout<<"Data has events after b-tagging = "<<v_Cuts_data.at(0)->GetBinContent(14)<<std::endl;
  
  std::cout<<" === SR Efficiency === "<<std::endl;
  for (unsigned int i=0; i<v_Cuts.size(); ++i)
  {
    double SREff=(v_Cuts.at(i)->GetBinContent(16))/(v_Cuts.at(i)->GetBinContent(14));
    std::cout<<"mX = "<<massPts[i]<<" has SR efficiency = "<<SREff*100.<<"%"<<std::endl;
  }
  std::cout<<"Data has events after SR = "<<v_Cuts_data.at(0)->GetBinContent(16)<<std::endl;
  
  // Cumulative eff wrt trigger
  std::cout<<" === Cumulative Efficiencies w.r.t. Trigger === "<<std::endl;
  for (unsigned int i=0; i<6; ++i)
  {
    double cumEff=(v_Cuts.at(i)->GetBinContent(16))/(v_Cuts.at(i)->GetBinContent(4));
    std::cout<<"mX = "<<massPts[i]<<" has cumulative efficiency wrt trigger = "<<cumEff*100.<<"%"<<std::endl;
  }
  
  // Cumulative eff wrt Production
  std::cout<<" === Cumulative Efficiencies w.r.t. Production === "<<std::endl;
  for (unsigned int i=0; i<v_Cuts.size(); ++i)
  {
    double cumEff=(v_Cuts.at(i)->GetBinContent(16))/(1e5);
    std::cout<<"mX = "<<massPts[i]<<" has cumulative efficiency wrt production = "<<cumEff*100.<<"%"<<std::endl;
  }
}

TCanvas* EfficiencyPlot(std::vector<TH1F*> v_Cuts, std::vector<TH1F*> v_CountWithPU)
{
  int massPts[15]={270, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100};
  
  TH1F *h_Init=new TH1F("h_Init", "h_Init", 88, 320, 1200);
  std::vector<TH1F *> h_Cut;
  
  TH1F *h_Events2=(TH1F*)h_Init->Clone("h_Events2");
  TH1F *h_Events4=(TH1F*)h_Init->Clone("h_Events4");
  TH1F *h_Events6=(TH1F*)h_Init->Clone("h_Events6");
  TH1F *h_Events8=(TH1F*)h_Init->Clone("h_Events8");
  TH1F *h_Events10=(TH1F*)h_Init->Clone("h_Events10");
  TH1F *h_Events14=(TH1F*)h_Init->Clone("h_Events14");
  TH1F *h_Events16=(TH1F*)h_Init->Clone("h_Events16");
  for (unsigned int i=0; i<v_Cuts.size(); ++i)
  {
    h_Init->SetBinContent(h_Init->FindBin(massPts[i]), v_CountWithPU.at(i)->GetBinContent(1));
    // std::cout<<"v_CountWithPU.at(i)->GetBinContent(1) = "<<v_CountWithPU.at(i)->GetBinContent(1)<<std::endl;
    // h_Events2->SetBinContent(i, v_Cuts.at(i)->GetBinContent(2));
    h_Events4->SetBinContent(h_Init->FindBin(massPts[i]), v_Cuts.at(i)->GetBinContent(4));
    h_Events6->SetBinContent(h_Init->FindBin(massPts[i]), v_Cuts.at(i)->GetBinContent(6));
    h_Events8->SetBinContent(h_Init->FindBin(massPts[i]), v_Cuts.at(i)->GetBinContent(8));
    h_Events10->SetBinContent(h_Init->FindBin(massPts[i]), v_Cuts.at(i)->GetBinContent(10));
    h_Events14->SetBinContent(h_Init->FindBin(massPts[i]), v_Cuts.at(i)->GetBinContent(14));
    h_Events16->SetBinContent(h_Init->FindBin(massPts[i]), v_Cuts.at(i)->GetBinContent(16));
    
    std::cout<<"m_X (GeV) = "<<massPts[i]<<", Efficiency = "<<v_Cuts.at(i)->GetBinContent(16)/v_CountWithPU.at(i)->GetBinContent(1)<<std::endl;
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
  
  return c_Ae;
} 

void CaterinaCompareTable(std::vector<TH1F *> v_mX_SR, std::vector<TH1F *> v_mX_CR2_data, std::vector<TH1F *> v_mX_CR4_data)
{
  int massPts[8]={300, 400, 500, 600, 700, 800, 0};
  
  // Cumulative eff wrt Production in +-30 GeV window around mass
  std::cout<<" === Cumulative Efficiencies in 30 GeV window around mass w.r.t. Production === "<<std::endl;
  for (unsigned int i=0; i<6; ++i)
  {
    // Draw and save each Signal SR plot
    TCanvas *c_mX_SR=new TCanvas("c_mX_SR", "", 700, 700);
    v_mX_SR.at(i)->Draw();
    c_mX_SR->SaveAs(("c_mX"+itoa(massPts[i])+"_SR_Signal.png").c_str());
    
    double massCenter;
    /*
    if (massPts[i]==300) massCenter=304;
    else if (massPts[i]==400) massCenter=406;
    else if (massPts[i]==500) massCenter=509;
    else if (massPts[i]==600) massCenter=614;
    else if (massPts[i]==700) massCenter=718;
    else if (massPts[i]==800) massCenter=820;
    */
    massCenter=v_mX_SR.at(i)->GetMean();
    
    double sigYield=v_mX_SR.at(i)->Integral(v_mX_SR.at(i)->FindBin(massCenter-30.), v_mX_SR.at(i)->FindBin(massCenter+30.));
    double bgYield=v_mX_CR2_data.at(0)->Integral(v_mX_CR2_data.at(0)->FindBin(massCenter-30.), v_mX_CR2_data.at(0)->FindBin(massCenter+30.))
                  +v_mX_CR4_data.at(0)->Integral(v_mX_CR4_data.at(0)->FindBin(massCenter-30.), v_mX_CR4_data.at(0)->FindBin(massCenter+30.));
    std::cout<<"mX = "<<massPts[i]<<", signal mean @ "<<v_mX_SR.at(i)->GetMean()
             <<" has 30 GeV-window #events in 18.6 GeV = "<<sigYield*totalLuminosity*prodXsec_1/nSignal_init
             <<" and background = "<<bgYield<<std::endl;
  }
}

void DrawSignal(std::vector<TH1F *> v, bool line=true, bool logPlot=false, double arrowX1=-999, double arrowX2=-999)
{
  std::string c_name=v.at(0)->GetName();
  c_name.replace(0, 1, "c");
  double max=0, min=1000.;
  if (v.size()>=6)
  {
    TCanvas *c_Signal=new TCanvas((c_name+"_Signal").c_str(), (c_name+"_Signal").c_str(), 800, 800);
    if (logPlot) c_Signal->SetLogy();
    if (line) {v.at(0)->Draw(""); v.at(0)->Draw("E SAME");} // C HIST
    else {v.at(0)->Draw();}
    max=v.at(0)->GetMaximum();
    min=v.at(0)->GetMinimum();
    for (unsigned int i=1; i<v.size(); ++i)
    {
      if (line) {v.at(i)->Draw("SAME"); v.at(i)->Draw("E SAME");} // C HIST
      else v.at(i)->Draw("SAME");
      if (v.at(i)->GetMaximum()>max) max=v.at(i)->GetMaximum();
      if (v.at(i)->GetMinimum()<min) min=v.at(i)->GetMinimum();
    }
    v.at(0)->SetMaximum(max*1.2);
    if (logPlot) v.at(0)->SetMinimum(min+10.);
    TLegend *leg_s=new TLegend(0.6, 0.6, 0.89, 0.89);
    leg_s->SetLineColor(0); leg_s->SetFillColor(0);
    leg_s->AddEntry((TObject*)0, "Signal Monte Carlo", "");
    if (all) leg_s->AddEntry(v.at(0), "m_{X}=450 GeV");
    leg_s->AddEntry(v.at(1), "m_{X}=500 GeV");
    if (all) leg_s->AddEntry(v.at(2), "m_{X}=550 GeV");
    leg_s->AddEntry(v.at(3), "m_{X}=600 GeV");
    if (all) leg_s->AddEntry(v.at(4), "m_{X}=650 GeV");
    leg_s->AddEntry(v.at(5), "m_{X}=700 GeV");
    if (all) leg_s->AddEntry(v.at(6), "m_{X}=800 GeV");
    if (all) leg_s->AddEntry(v.at(7), "m_{X}=900 GeV");
    leg_s->Draw();
    if (arrowX1!=-999)
    {
      TArrow *arrow=new TArrow(arrowX1, max*1.1, arrowX1, 0);
      arrow->Draw();
    }
    if (arrowX2!=-999)
    {
      TArrow *arrow=new TArrow(arrowX2, max*1.1, arrowX2, 0);
      arrow->Draw();
    }
    c_Signal->SaveAs((c_name+"_Signal.png").c_str());
  }
  else std::cout<<"ERROR: AT LEAST 6 SIGNAL SAMPLES NOT FOUND!"<<std::endl;
}

TCanvas* DrawQCD(std::vector<TH1F *> v, bool line=true, bool logPlot=false, double arrowX1=-999, double arrowX2=-999)
{
  std::string c_name=v.at(0)->GetName();
  c_name.replace(0, 1, "c");
  double max=0, min=1000.;
  
  TCanvas *c_QCD=new TCanvas((c_name+"_QCD").c_str(), (c_name+"_QCD").c_str(), 800, 800);
  THStack *s_QCD=new THStack((c_name+"_QCD").c_str(), (c_name+"_QCD").c_str());
  v.at(0)->Scale(scale_QCDpT50To80);
  v.at(1)->Scale(scale_QCDpT80To120);
  v.at(2)->Scale(scale_QCDpT120To170);
  v.at(3)->Scale(scale_QCDpT170To300);
  v.at(4)->Scale(scale_QCDpT300To470);
  v.at(5)->Scale(scale_QCDpT470To600);
  v.at(6)->Scale(scale_QCDpT600To800);
  v.at(7)->Scale(scale_QCDpT800To1000);
  v.at(8)->Scale(scale_QCDpT1000To1400);
  v.at(9)->Scale(scale_QCDpT1400To1800);
  v.at(10)->Scale(scale_QCDpT1800);
  for (unsigned int i=0; i<v.size(); ++i) s_QCD->Add(v.at(i));
  s_QCD->Draw();
  TLegend *leg=new TLegend(0.6, 0.9, 0.85, 0.7);
  leg->AddEntry(v.at(0), "QCD pT 50-80");
  leg->AddEntry(v.at(1), "QCD pT 80-120");
  leg->AddEntry(v.at(2), "QCD pT 120-170");
  leg->AddEntry(v.at(3), "QCD pT 170-300");
  leg->AddEntry(v.at(4), "QCD pT 300-470");
  leg->AddEntry(v.at(5), "QCD pT 470-600");
  leg->AddEntry(v.at(6), "QCD pT 600-800");
  leg->AddEntry(v.at(7), "QCD pT 800-1000");
  leg->AddEntry(v.at(8), "QCD pT 1000-1400");
  leg->AddEntry(v.at(9), "QCD pT 1400-1800");
  leg->AddEntry(v.at(10), "QCD pT 1800");
  leg->Draw();
  if (arrowX1!=-999)
  {
    TArrow *arrow=new TArrow(arrowX1, max*1.1, arrowX1, 0);
    arrow->Draw();
  }
  if (arrowX2!=-999)
  {
    TArrow *arrow=new TArrow(arrowX2, max*1.1, arrowX2, 0);
    arrow->Draw();
  }
  s_QCD->SetMinimum(100.);
  c_QCD->SetLogy();
  c_QCD->SaveAs((c_name+"_QCD.png").c_str());
  
  return c_QCD;
}

void fillSamples(Samples &s, string region)
{
  TFile *qcdpT50To80=new TFile((std::string("../")+region+std::string("/Histograms_QCD_Pt-50to80_TuneZ2star_8TeV_pythia6_Skim.root")).c_str());
  TFile *qcdpT80To120=new TFile((std::string("../")+region+std::string("Histograms_QCD_Pt-80to120_TuneZ2star_8TeV_pythia6_Skim.root")).c_str());
  TFile *qcdpT120To170=new TFile((std::string("../")+region+std::string("Histograms_QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Skim.root")).c_str());
  TFile *qcdpT170To300=new TFile((std::string("../")+region+std::string("Histograms_QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Skim.root")).c_str());
  TFile *qcdpT300To470=new TFile((std::string("../")+region+std::string("Histograms_QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Skim.root")).c_str());
  TFile *qcdpT470To600=new TFile((std::string("../")+region+std::string("Histograms_QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Skim.root")).c_str());
  TFile *qcdpT600To800=new TFile((std::string("../")+region+std::string("Histograms_QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Skim.root")).c_str());
  TFile *qcdpT800To1000=new TFile((std::string("../")+region+std::string("Histograms_QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Skim.root")).c_str());
  TFile *qcdpT1000To1400=new TFile((std::string("../")+region+std::string("Histograms_QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Skim.root")).c_str());
  TFile *qcdpT1400To1800=new TFile((std::string("../")+region+std::string("Histograms_QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Skim.root")).c_str());
  TFile *qcdpT1800=new TFile((std::string("../")+region+std::string("Histograms_QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Skim.root")).c_str());
  
  pushBackHistograms(s, qcdpT50To80);
  pushBackHistograms(s, qcdpT80To120);
  pushBackHistograms(s, qcdpT120To170);
  pushBackHistograms(s, qcdpT170To300);
  pushBackHistograms(s, qcdpT300To470);
  pushBackHistograms(s, qcdpT470To600);
  pushBackHistograms(s, qcdpT600To800);
  pushBackHistograms(s, qcdpT800To1000);
  pushBackHistograms(s, qcdpT1000To1400);
  pushBackHistograms(s, qcdpT1400To1800);
  pushBackHistograms(s, qcdpT1800);

}

void DrawQCDData(std::vector<TH1F *> v, std::vector<TH1F *> v_data, bool line=true, bool logPlot=false, double arrowX1=-999, double arrowX2=-999, double ratioMin=0., double ratioMax=0.)
{
  std::string c_name=v.at(0)->GetName();
  c_name.replace(0, 1, "c");
  double max=0, min=1000.;
  max=v_data.at(0)->GetMaximum();
  
  TCanvas *c_QCD=new TCanvas((c_name+"_QCD").c_str(), (c_name+"_QCD").c_str(), 800, 800);
  TPad *p_2=new TPad("p_2", "p_2", 0, 0, 1, 0.3);
  TPad *p_1=new TPad("p_1", "p_1", 0, 0.3, 1, 1);
  p_1->Draw();
  p_2->Draw();
  p_1->cd();
  if (logPlot) p_1->SetLogy();
  THStack *s_QCD=new THStack((c_name+"_QCD").c_str(), v.at(0)->GetTitle());
  v.at(0)->Scale(scale_QCDpT50To80);
  v.at(1)->Scale(scale_QCDpT80To120);
  v.at(2)->Scale(scale_QCDpT120To170);
  v.at(3)->Scale(scale_QCDpT170To300);
  v.at(4)->Scale(scale_QCDpT300To470);
  v.at(5)->Scale(scale_QCDpT470To600);
  v.at(6)->Scale(scale_QCDpT600To800);
  v.at(7)->Scale(scale_QCDpT800To1000);
  v.at(8)->Scale(scale_QCDpT1000To1400);
  v.at(9)->Scale(scale_QCDpT1400To1800);
  v.at(10)->Scale(scale_QCDpT1800);
  TH1F *h_sum=(TH1F*)v.at(0)->Clone("h_sum");
  h_sum->Sumw2();
  for (unsigned int i=0; i<v.size(); ++i) 
  {
    s_QCD->Add(v.at(i), "hist");
    if (i!=0) h_sum->Add(v.at(i));
  }
  s_QCD->Draw();
  s_QCD->SetMaximum(max*1.1);
  s_QCD->SetTitle((std::string("QCD MC & Data ")+v.at(0)->GetTitle()).c_str());
  s_QCD->GetXaxis()->SetTitle(v.at(0)->GetXaxis()->GetTitle());
  s_QCD->GetYaxis()->SetTitle(v.at(0)->GetYaxis()->GetTitle());
  v_data.at(0)->Draw("E SAME");
  TLegend *leg=new TLegend(0.6, 0.9, 0.85, 0.7);
  leg->AddEntry(v.at(0), "QCD pT 50-80");
  leg->AddEntry(v.at(1), "QCD pT 80-120");
  leg->AddEntry(v.at(2), "QCD pT 120-170");
  leg->AddEntry(v.at(3), "QCD pT 170-300");
  leg->AddEntry(v.at(4), "QCD pT 300-470");
  leg->AddEntry(v.at(5), "QCD pT 470-600");
  leg->AddEntry(v.at(6), "QCD pT 600-800");
  leg->AddEntry(v.at(7), "QCD pT 800-1000");
  leg->AddEntry(v.at(8), "QCD pT 1000-1400");
  leg->AddEntry(v.at(9), "QCD pT 1400-1800");
  leg->AddEntry(v.at(10), "QCD pT 1800");
  leg->AddEntry(v_data.at(0), "Data");
  leg->Draw();
  if (arrowX1!=-999)
  {
    TArrow *arrow=new TArrow(arrowX1, max*0.6, arrowX1, 0);
    arrow->Draw();
  }
  if (arrowX2!=-999)
  {
    TArrow *arrow=new TArrow(arrowX2, max*0.6, arrowX2, 0);
    arrow->Draw();
  }
  if (logPlot) s_QCD->SetMinimum(h_sum->GetMaximum()/1e5);
  p_2->cd();
  p_2->SetGridy();
  TH1F *h_ratio=(TH1F*)v_data.at(0)->Clone("h_ratio");
  h_ratio->SetTitle(("Data/MC Ratio "+std::string(h_ratio->GetTitle())+"; "+v.at(0)->GetXaxis()->GetTitle()).c_str());
  h_ratio->Divide(h_sum);
  if (ratioMin!=0. || ratioMax!=0.) 
  {
    h_ratio->SetMaximum(ratioMax); 
    h_ratio->SetMinimum(ratioMin);
  }
  h_ratio->Draw();
  p_1->cd();
  c_QCD->SaveAs((c_name+"_QCD.png").c_str());
}

void DrawData(std::vector<TH1F *> v, bool line=true, bool logPlot=false, double arrowX1=-999, double arrowX2=-999)
{
  std::string c_name=v.at(0)->GetName();
  c_name.replace(0, 1, "c");
  double max=0, min=1000.;
  
  TCanvas *c_Data=new TCanvas((c_name+"_Data").c_str(), (c_name+"_Data").c_str(), 800, 800);
  if (line) {v.at(0)->Draw(""); v.at(0)->Draw("E SAME");} // C HIST
  else v.at(0)->Draw();
  TLegend *leg_s=new TLegend(0.5, 0.6, 0.89, 0.89);
  leg_s->SetLineColor(0); leg_s->SetFillColor(0);
  leg_s->AddEntry((TObject*)0, "DATA", "");
  leg_s->AddEntry(v.at(0), "8 TeV pp Data 17.928 /pb");
  leg_s->Draw();
  max=v.at(0)->GetMaximum();
  if (arrowX1!=-999)
  {
    TArrow *arrow=new TArrow(arrowX1, max*0.9, arrowX1, 0);
    arrow->Draw();
  }
  if (arrowX2!=-999)
  {
    TArrow *arrow=new TArrow(arrowX2, max*0.9, arrowX2, 0);
    arrow->Draw();
  }
  if (logPlot) c_Data->SetLogy();
  if (logPlot) v.at(0)->SetMinimum(min+1.);
  c_Data->SaveAs((c_name+"_Data.png").c_str());
}

void lineColorSignal(std::vector<TH1F *> v)
{
  v.at(0)->SetLineColor(kGreen);
  v.at(1)->SetLineColor(kOrange);
  v.at(2)->SetLineColor(kMagenta);
  v.at(3)->SetLineColor(kBlue);
  v.at(4)->SetLineColor(kCyan);
  v.at(5)->SetLineColor(kRed);
  v.at(6)->SetLineColor(kBlack);
  v.at(6)->SetLineColor(kViolet);
  for (unsigned int i=0; i<v.size(); ++i)
  {
    v.at(i)->SetLineWidth(2);
  }
}

void lineColorQCD(std::vector<TH1F *> v)
{
  v.at(0)->SetFillColor(kGreen);
  v.at(1)->SetFillColor(kGreen+2);
  v.at(2)->SetFillColor(kOrange);
  v.at(3)->SetFillColor(kOrange+2);
  v.at(4)->SetFillColor(kMagenta);
  v.at(5)->SetFillColor(kMagenta+2);
  v.at(6)->SetFillColor(kBlue);
  v.at(7)->SetFillColor(kBlue+2);
  v.at(8)->SetFillColor(kCyan);
  v.at(9)->SetFillColor(kCyan+2);
  v.at(10)->SetFillColor(kRed);
  for (unsigned int i=0; i<v.size(); ++i) v.at(i)->SetLineWidth(2);
}

void lineColorData(std::vector<TH1F *> v)
{
  v.at(0)->SetLineColor(kBlack);
  v.at(0)->SetLineWidth(2);
  // v.at(0)->Sumw2();
}

void rebin(std::vector<TH1F *> v, int n)
{
  for (unsigned int i=0; i<v.size(); ++i)
    v.at(i)->Rebin(n);
}

void setTitle(std::vector<TH1F*> v, std::string title)
{
  for (unsigned int i=0; i<v.size(); ++i)
  {
    v.at(i)->SetTitle(title.c_str());
  }
}

int DisplayDistributions()
{
  if (signal_sets)
  {
    std::cout<<" === Opening Signal MC === "<<std::endl;
    
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
    
    // pushBackHistograms(signals, glugluToX270);
    // pushBackHistograms(signals, glugluToX300);
    // pushBackHistograms(signals, glugluToX350);
    // pushBackHistograms(signals, glugluToX400);
    pushBackHistograms(signals, glugluToX450);
    pushBackHistograms(signals, glugluToX500);
    pushBackHistograms(signals, glugluToX550);
    pushBackHistograms(signals, glugluToX600);
    pushBackHistograms(signals, glugluToX650);
    pushBackHistograms(signals, glugluToX700);
    pushBackHistograms(signals, glugluToX800);
    pushBackHistograms(signals, glugluToX900);
    // pushBackHistograms(signals, glugluToX1000);
    // pushBackHistograms(signals, glugluToX1100);
  }
  
  if (qcd_sets)
  {
    std::cout<<" === Opening QCD MC === "<<std::endl;
  
    TFile *qcdpT50To80=new TFile("/Users/souvik/HbbHbb/Analysis/FullData/MMMM_nominal/a_KinFit/Histograms_QCD_Pt-50to80_TuneZ2star_8TeV_pythia6_Skim.root");
    std::cout<<"Opened Histograms_QCD_Pt-50to80_TuneZ2star_8TeV_pythia6_Skim.root"<<std::endl;
    pushBackHistograms(QCDs, qcdpT50To80);
  
    TFile *qcdpT80To120=new TFile("/Users/souvik/HbbHbb/Analysis/FullData/MMMM_nominal/a_KinFit/Histograms_QCD_Pt-80to120_TuneZ2star_8TeV_pythia6_Skim.root");
    std::cout<<"Opened Histograms_QCD_Pt-80to120_TuneZ2star_8TeV_pythia6_Skim.root"<<std::endl;
    pushBackHistograms(QCDs, qcdpT80To120);
  
    TFile *qcdpT120To170=new TFile("/Users/souvik/HbbHbb/Analysis/FullData/MMMM_nominal/a_KinFit/Histograms_QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Skim.root");
    std::cout<<"Opened Histograms_QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Skim.root"<<std::endl;
    pushBackHistograms(QCDs, qcdpT120To170);
  
    TFile *qcdpT170To300=new TFile("/Users/souvik/HbbHbb/Analysis/FullData/MMMM_nominal/a_KinFit/Histograms_QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Skim.root");
    std::cout<<"Opened Histograms_QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Skim.root"<<std::endl;
    pushBackHistograms(QCDs, qcdpT170To300);
  
    TFile *qcdpT300To470=new TFile("/Users/souvik/HbbHbb/Analysis/FullData/MMMM_nominal/a_KinFit/Histograms_QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Skim.root");
    std::cout<<"Opened Histograms_QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Skim.root"<<std::endl;
    pushBackHistograms(QCDs, qcdpT300To470);
  
    TFile *qcdpT470To600=new TFile("/Users/souvik/HbbHbb/Analysis/FullData/MMMM_nominal/a_KinFit/Histograms_QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Skim.root");
    std::cout<<"Opened Histograms_QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Skim.root"<<std::endl;
    pushBackHistograms(QCDs, qcdpT470To600);
  
    TFile *qcdpT600To800=new TFile("/Users/souvik/HbbHbb/Analysis/FullData/MMMM_nominal/a_KinFit/Histograms_QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Skim.root");
    std::cout<<"Opened Histograms_QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Skim.root"<<std::endl;
    pushBackHistograms(QCDs, qcdpT600To800);
  
    TFile *qcdpT800To1000=new TFile("/Users/souvik/HbbHbb/Analysis/FullData/MMMM_nominal/a_KinFit/Histograms_QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Skim.root");
    std::cout<<"Opened Histograms_QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Skim.root"<<std::endl;
    pushBackHistograms(QCDs, qcdpT800To1000);
  
    TFile *qcdpT1000To1400=new TFile("/Users/souvik/HbbHbb/Analysis/FullData/MMMM_nominal/a_KinFit/Histograms_QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Skim.root");
    std::cout<<"Opened Histograms_QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6_Skim.root"<<std::endl;
    pushBackHistograms(QCDs, qcdpT1000To1400);
  
    TFile *qcdpT1400To1800=new TFile("/Users/souvik/HbbHbb/Analysis/FullData/MMMM_nominal/a_KinFit/Histograms_QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Skim.root");
    std::cout<<"Opened Histograms_QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6_Skim.root"<<std::endl;
    pushBackHistograms(QCDs, qcdpT1400To1800);
  
    TFile *qcdpT1800=new TFile("/Users/souvik/HbbHbb/Analysis/FullData/MMMM_nominal/a_KinFit/Histograms_QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Skim.root");
    std::cout<<"Opened Histograms_QCD_Pt-1800_TuneZ2star_8TeV_pythia6_Skim.root"<<std::endl;
    pushBackHistograms(QCDs, qcdpT1800);
    
    std::cout<<" === Scale Factors === "<<std::endl;
    std::cout<<" scale_QCDpT50To80 = "<<scale_QCDpT50To80<<std::endl;
    std::cout<<" scale_QCDpT80To120 = "<<scale_QCDpT80To120<<std::endl;
    std::cout<<" scale_QCDpT120To170 = "<<scale_QCDpT120To170<<std::endl;
    std::cout<<" scale_QCDpT170To300 = "<<scale_QCDpT170To300<<std::endl;
    std::cout<<" scale_QCDpT300To470 = "<<scale_QCDpT300To470<<std::endl;
    std::cout<<" scale_QCDpT470To600 = "<<scale_QCDpT470To600<<std::endl;
    std::cout<<" scale_QCDpT600To800 = "<<scale_QCDpT600To800<<std::endl;
    std::cout<<" scale_QCDpT800To1000 = "<<scale_QCDpT800To1000<<std::endl;
    std::cout<<" scale_QCDpT1000To1400 = "<<scale_QCDpT1000To1400<<std::endl;
    std::cout<<" scale_QCDpT1400To1800 = "<<scale_QCDpT1400To1800<<std::endl;
    std::cout<<" scale_QCDpT1800 = "<<scale_QCDpT1800<<std::endl;
  }
  
  if (data_sets)
  {
    std::cout<<" === Opening Data === "<<std::endl;
  
    TFile *data_8TeVData2012B_part1=new TFile("Histograms_BJetPlusX_Run2012BCD_Skim.root");
    std::cout<<"Opened Histograms_8TeVData2012BCD_Skim.root"<<std::endl;
    pushBackHistograms(data, data_8TeVData2012B_part1, false);
  }
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  
  if (signal_sets)
  {
    lineColorSignal(signals.v_nJets); setTitle(signals.v_nJets, "; n_{jets} with p_{T} > 40 GeV, |#eta| < 2.5, CMVA > 0.71"); 
    lineColorSignal(signals.v_nPV);
    lineColorSignal(signals.v_CSV1);
    lineColorSignal(signals.v_CSV2);
    lineColorSignal(signals.v_CSV3);
    lineColorSignal(signals.v_CSV4);
    lineColorSignal(signals.v_JetpT1);
    lineColorSignal(signals.v_JetpT2);
    lineColorSignal(signals.v_JetpT3);
    lineColorSignal(signals.v_JetpT4);
    lineColorSignal(signals.v_HH_mass_diff);
    lineColorSignal(signals.v_H1_mass); setTitle(signals.v_H1_mass, "Events/3 GeV; m_{H} (GeV)");
    lineColorSignal(signals.v_H2_mass);
    lineColorSignal(signals.v_H1_CSV1);
    lineColorSignal(signals.v_H1_CSV2);
    lineColorSignal(signals.v_H2_CSV1);
    lineColorSignal(signals.v_H2_CSV2);
    lineColorSignal(signals.v_nTags);
    lineColorSignal(signals.v_HH_dPhi);
    lineColorSignal(signals.v_HH_deta);
    lineColorSignal(signals.v_Azimuth);
    lineColorSignal(signals.v_H1_pT);
    lineColorSignal(signals.v_H2_pT);
    lineColorSignal(signals.v_X_mass);
    lineColorSignal(signals.v_X_pT);
    lineColorSignal(signals.v_MET);
    lineColorSignal(signals.v_MET_sig);
    lineColorSignal(signals.v_HT);
    lineColorSignal(signals.v_genX_mass); 
    lineColorSignal(signals.v_genH1_mass);
    lineColorSignal(signals.v_genH2_mass);
    lineColorSignal(signals.v_H1_mass_bTagged);
    lineColorSignal(signals.v_H2_mass_bTagged);
    lineColorSignal(signals.v_H1_pT_bTagged); 
    lineColorSignal(signals.v_H2_pT_bTagged); 
    lineColorSignal(signals.v_X_mass_bTagged);
    lineColorSignal(signals.v_X_pT_bTagged);
    lineColorSignal(signals.v_mX_CR1);
    lineColorSignal(signals.v_mX_CR2);
    lineColorSignal(signals.v_mX_CR3);
    lineColorSignal(signals.v_mX_CR4);
    lineColorSignal(signals.v_mX_SR);
    lineColorSignal(signals.v_mX_CR24);
    lineColorSignal(signals.v_dR1);
    lineColorSignal(signals.v_dR1_KinFit);
  }
  
  if (qcd_sets)
  {
    lineColorQCD(QCDs.v_nJets);           
    lineColorQCD(QCDs.v_nPV);             
    lineColorQCD(QCDs.v_CSV1);            
    lineColorQCD(QCDs.v_CSV2);            
    lineColorQCD(QCDs.v_CSV3);            
    lineColorQCD(QCDs.v_CSV4);            
    lineColorQCD(QCDs.v_JetpT1);          
    lineColorQCD(QCDs.v_JetpT2);          
    lineColorQCD(QCDs.v_JetpT3);          
    lineColorQCD(QCDs.v_JetpT4);          
    lineColorData(data.v_HH_mass_diff);
    lineColorQCD(QCDs.v_H1_mass);         
    lineColorQCD(QCDs.v_H2_mass);         
    lineColorQCD(QCDs.v_H1_CSV1);         
    lineColorQCD(QCDs.v_H1_CSV2);         
    lineColorQCD(QCDs.v_H2_CSV1);         
    lineColorQCD(QCDs.v_H2_CSV2);         
    lineColorQCD(QCDs.v_nTags);           
    lineColorQCD(QCDs.v_HH_dPhi);         
    lineColorQCD(QCDs.v_HH_deta);         
    lineColorQCD(QCDs.v_Azimuth);         
    lineColorQCD(QCDs.v_H1_pT);           
    lineColorQCD(QCDs.v_H2_pT);           
    lineColorQCD(QCDs.v_X_mass);          
    lineColorQCD(QCDs.v_X_pT);            
    lineColorQCD(QCDs.v_MET);             
    lineColorQCD(QCDs.v_MET_sig);         
    lineColorQCD(QCDs.v_HT);              
    lineColorQCD(QCDs.v_genX_mass);
    lineColorQCD(QCDs.v_genH1_mass);
    lineColorQCD(QCDs.v_genH2_mass);
    lineColorQCD(QCDs.v_H1_mass_bTagged); 
    lineColorQCD(QCDs.v_H2_mass_bTagged); 
    lineColorQCD(QCDs.v_H1_pT_bTagged);   
    lineColorQCD(QCDs.v_H2_pT_bTagged);   
    lineColorQCD(QCDs.v_X_mass_bTagged);  
    lineColorQCD(QCDs.v_X_pT_bTagged);    
    lineColorQCD(QCDs.v_mX_CR1);          
    lineColorQCD(QCDs.v_mX_CR2);          
    lineColorQCD(QCDs.v_mX_CR3);          
    lineColorQCD(QCDs.v_mX_CR4);          
    lineColorQCD(QCDs.v_mX_SR);           
    lineColorQCD(QCDs.v_mX_CR24);
  }
  
  if (data_sets)
  {
    lineColorData(data.v_nJets); setTitle(data.v_nJets, "; n_{jets} with p_{T} > 40 GeV, |#eta| < 2.5, CMVA > 0.71"); 
    lineColorData(data.v_nPV);
    lineColorData(data.v_CSV1);
    lineColorData(data.v_CSV2);
    lineColorData(data.v_CSV3);
    lineColorData(data.v_CSV4);
    lineColorData(data.v_JetpT1);
    lineColorData(data.v_JetpT2);
    lineColorData(data.v_JetpT3);
    lineColorData(data.v_JetpT4);
    lineColorData(data.v_H1_mass); setTitle(data.v_H1_mass, "Events/3 GeV; m_{H} (GeV)");
    lineColorData(data.v_H2_mass);
    lineColorData(data.v_H1_CSV1);
    lineColorData(data.v_H1_CSV2);
    lineColorData(data.v_H2_CSV1);
    lineColorData(data.v_H2_CSV2);
    lineColorData(data.v_nTags);
    lineColorData(data.v_HH_dPhi);
    lineColorData(data.v_HH_deta);
    lineColorData(data.v_Azimuth);
    lineColorData(data.v_H1_pT);
    lineColorData(data.v_H2_pT);
    lineColorData(data.v_X_mass);
    lineColorData(data.v_X_pT);
    lineColorData(data.v_MET);
    lineColorData(data.v_MET_sig);
    lineColorData(data.v_HT);
    lineColorData(data.v_H1_mass_bTagged);
    lineColorData(data.v_H2_mass_bTagged);
    lineColorData(data.v_H1_pT_bTagged);
    lineColorData(data.v_H2_pT_bTagged);
    lineColorData(data.v_X_mass_bTagged);
    lineColorData(data.v_X_pT_bTagged);
    lineColorData(data.v_mX_CR1);
    lineColorData(data.v_mX_CR2);
    lineColorData(data.v_mX_CR3);
    lineColorData(data.v_mX_CR4);
    lineColorData(data.v_mX_SR);
    lineColorData(data.v_mX_CR24);
    lineColorData(data.v_dR1);
    lineColorData(data.v_dR1_KinFit);
  }
  
  if (signal_sets)
  {
    DrawSignal(signals.v_genX_mass);
    DrawSignal(signals.v_genH1_mass);
    DrawSignal(signals.v_genH2_mass);
    DrawSignal(signals.v_nJets, false, true, 4);
    DrawSignal(signals.v_JetpT1, false, true); 
    DrawSignal(signals.v_JetpT2, false, true); 
    DrawSignal(signals.v_JetpT3, false, true); 
    DrawSignal(signals.v_JetpT4, false, true); 
    DrawSignal(signals.v_HH_mass_diff);
    DrawSignal(signals.v_H1_mass);
    DrawSignal(signals.v_H2_mass);
    DrawSignal(signals.v_H1_pT, false, true);
    DrawSignal(signals.v_H2_pT, false, true);
    DrawSignal(signals.v_H1_CSV1, true, false, 0.679); 
    DrawSignal(signals.v_H1_CSV2, true, false, 0.679); 
    DrawSignal(signals.v_H2_CSV1, true, false, 0.679); 
    DrawSignal(signals.v_H2_CSV2, true, false, 0.679); 
    DrawSignal(signals.v_HH_dPhi); 
    DrawSignal(signals.v_HH_deta); 
    DrawSignal(signals.v_Azimuth);
    
    // Draw H properties after b-tagging
    DrawSignal(signals.v_H1_mass_bTagged); 
    DrawSignal(signals.v_H2_mass_bTagged); 
    DrawSignal(signals.v_H1_pT_bTagged, false, true);
    DrawSignal(signals.v_H2_pT_bTagged, false, true);
    
    rebin(signals.v_X_mass, 2); DrawSignal(signals.v_X_mass); 
    rebin(signals.v_X_pT, 2); DrawSignal(signals.v_X_pT, false, true); 
    
    rebin(signals.v_X_pT_bTagged, 4); DrawSignal(signals.v_X_pT_bTagged, false, true); 
    rebin(signals.v_X_mass_bTagged, 1); DrawSignal(signals.v_X_mass_bTagged); 
    
    DrawSignal(signals.v_MET, false, true);
    DrawSignal(signals.v_MET_sig, false, true);  
    DrawSignal(signals.v_HT, false, true);
    
    rebin(signals.v_mX_CR24, 4); DrawSignal(signals.v_mX_CR24); 
    DrawSignal(signals.v_mX_SR, false, true);
    
    DrawSignal(signals.v_nTags, false, false, 3.);
    
    DrawSignal(signals.v_dR1); 
    DrawSignal(signals.v_dR1_KinFit);
  }
  
  if (data_sets)
  {
    DrawData(data.v_nJets, false, true, 4);
    DrawData(data.v_JetpT1, false, true);
    DrawData(data.v_JetpT2, false, true);
    DrawData(data.v_JetpT3, false, true);
    DrawData(data.v_JetpT4, false, true);
    DrawData(data.v_HH_mass_diff);
    DrawData(data.v_H1_mass); 
    DrawData(data.v_H2_mass); 
    DrawData(data.v_H1_pT);
    DrawData(data.v_H2_pT);
    DrawData(data.v_H1_CSV1, true, false, 0.679); 
    DrawData(data.v_H1_CSV2, true, false, 0.679); 
    DrawData(data.v_H2_CSV1, true, false, 0.679); 
    DrawData(data.v_H2_CSV2, true, false, 0.679); 
    DrawData(data.v_HH_dPhi);
    DrawData(data.v_HH_deta);
    DrawData(data.v_Azimuth);
    
    // Draw H properties after b-tagging
    rebin(data.v_H1_mass_bTagged, 2);
    DrawData(data.v_H1_mass_bTagged); 
    DrawData(data.v_H2_mass_bTagged); 
    rebin(data.v_H1_pT_bTagged, 2);
    DrawData(data.v_H1_pT_bTagged); 
    DrawData(data.v_H2_pT_bTagged); 
    
    rebin(data.v_X_mass, 2); DrawData(data.v_X_mass); 
    rebin(data.v_X_pT, 2); DrawData(data.v_X_pT, false, true);
    
    DrawData(data.v_MET);
    DrawData(data.v_MET_sig); 
    DrawData(data.v_HT);
    
    DrawData(data.v_nTags, false, false, 3.);
    
    DrawData(data.v_dR1);
    DrawData(data.v_dR1_KinFit);
  }
  
  if (qcd_sets)
  {
    DrawQCD(QCDs.v_nJets, false, false, 4);
    DrawQCDData(QCDs.v_JetpT1, data.v_JetpT1, false, true);
    DrawQCDData(QCDs.v_JetpT2, data.v_JetpT2, false, true);
    DrawQCDData(QCDs.v_JetpT3, data.v_JetpT3, false, true);
    DrawQCDData(QCDs.v_JetpT4, data.v_JetpT4, false, true);
    DrawQCDData(QCDs.v_H1_mass, data.v_H1_mass); 
    DrawQCDData(QCDs.v_H2_mass, data.v_H2_mass);
    DrawQCDData(QCDs.v_H1_pT, data.v_H1_pT, false, true);
    DrawQCDData(QCDs.v_H2_pT, data.v_H2_pT, false, true);
    DrawQCDData(QCDs.v_H1_CSV1, data.v_H1_CSV1, true, false, 0.679);
    DrawQCDData(QCDs.v_H1_CSV2, data.v_H1_CSV2, true, false, 0.679);
    DrawQCDData(QCDs.v_H2_CSV1, data.v_H2_CSV1, true, false, 0.679);
    DrawQCDData(QCDs.v_H2_CSV2, data.v_H2_CSV2, true, false, 0.679);
    DrawQCD(QCDs.v_HH_dPhi);
    DrawQCD(QCDs.v_HH_deta);
    DrawQCD(QCDs.v_Azimuth);
    
    // Draw H properties after b-tagging
    rebin(QCDs.v_H1_mass_bTagged, 2);
    DrawQCDData(QCDs.v_H1_mass_bTagged, data.v_H1_mass_bTagged);
    DrawQCDData(QCDs.v_H2_mass_bTagged, data.v_H2_mass_bTagged);
    rebin(QCDs.v_H1_pT_bTagged, 2);
    DrawQCDData(QCDs.v_H1_pT_bTagged, data.v_H1_pT_bTagged, false, true);
    DrawQCDData(QCDs.v_H2_pT_bTagged, data.v_H2_pT_bTagged, false, true);
    
    rebin(QCDs.v_X_mass, 2); rebin(data.v_X_mass, 2); DrawQCDData(QCDs.v_X_mass, data.v_X_mass); 
    rebin(QCDs.v_X_pT, 2); rebin(data.v_X_pT, 2); DrawQCDData(QCDs.v_X_pT, data.v_X_pT, false, true); 
    
    rebin(QCDs.v_X_pT_bTagged, 4); rebin(data.v_X_pT_bTagged, 4); DrawQCDData(QCDs.v_X_pT_bTagged, data.v_X_pT_bTagged, false, true, -999, -999, 0, 5);
    rebin(QCDs.v_X_mass_bTagged, 1); rebin(data.v_X_mass_bTagged, 1); DrawQCDData(QCDs.v_X_mass_bTagged, data.v_X_mass_bTagged, false, false, -999, -999, 0, 5);
    
    DrawQCDData(QCDs.v_MET, data.v_MET, false, true); 
    DrawQCDData(QCDs.v_MET_sig, data.v_MET_sig, false, true); 
    DrawQCDData(QCDs.v_HT, data.v_HT, false, true); 
    
    rebin(QCDs.v_mX_CR24, 4); rebin(data.v_mX_CR24, 4); DrawQCDData(QCDs.v_mX_CR24, data.v_mX_CR24);
    
    DrawQCD(QCDs.v_nTags, false, false, 3.);
  }
  
  
  // Efficiency outputs
  // EfficiencyTable(signals.v_Cuts, data.v_Cuts);
  
  // Efficiency plot output
  EfficiencyPlot(signals.v_Cuts, signals.v_CountWithPU);
  
  // Compare with Caterina
  // CaterinaCompareTable(signals.v_mX_SR, data.v_mX_CR2, data.v_mX_CR4);
  
  
  return 0;
}
  
  
  
  
  
  
  
  
  
