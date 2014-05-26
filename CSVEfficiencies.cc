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
#include <stdlib.h>

double bTagCSV_tightCut=0.898;
double bTagCSV_mediumCut=0.679;
double bTagCSV_looseCut=0.244;
double bTagCSV_noCut=0.;

double cut_CSV=bTagCSV_tightCut;
double jetpT_cut=40.;
double jeteta_cut=2.5;

// QCD scale factors
double totalLuminosity=13242; // /pb
double xsec_QCDHT250To500=276000; // pb
double xsec_QCDHT500To1000=8426; // pb
double xsec_QCDHT1000ToInf=204; // pb
double init_QCDHT250To500=4669914;
double init_QCDHT500To1000=4779910;
double init_QCDHT1000ToInf=5192110;
double scale_QCDHT250To500=totalLuminosity*xsec_QCDHT250To500/init_QCDHT250To500;
double scale_QCDHT500To1000=totalLuminosity*xsec_QCDHT500To1000/init_QCDHT500To1000;
double scale_QCDHT1000ToInf=totalLuminosity*xsec_QCDHT1000ToInf/init_QCDHT1000ToInf;

void CSVEfficiencies(std::string sample0,
                     std::string sample1,
                     std::string sample2,
                     std::string PUWeight0,
                     std::string PUWeight1,
                     std::string PUWeight2,
                     std::string csvTag="M")
{

  if (csvTag=="T") cut_CSV=bTagCSV_tightCut;
  if (csvTag=="M") cut_CSV=bTagCSV_mediumCut;
  if (csvTag=="L") cut_CSV=bTagCSV_looseCut;
  
  TH1F *h_udsg_JetpT_den=new TH1F("h_udsg_JetpT_den", "h_udsg_JetpT_den", 30, 0., 300.); h_udsg_JetpT_den->Sumw2();
  TH1F *h_udsg_Jeteta_den=new TH1F("h_udsg_Jeteta_den", "h_udsg_Jeteta_den", 30, 0., 2.5); h_udsg_Jeteta_den->Sumw2();
  TH2F *h_udsg_JetpTeta_den=new TH2F("h_udsg_JetpTeta_den", "h_udsg_JetpTeta_den", 30, 0., 300., 30, 0., 2.5); h_udsg_JetpTeta_den->Sumw2();
  TH1F *h_udsg_JetpT_num=new TH1F("h_udsg_JetpT_num", "h_udsg_JetpT_num", 30, 0., 300.); h_udsg_JetpT_num->Sumw2();
  TH1F *h_udsg_Jeteta_num=new TH1F("h_udsg_Jeteta_num", "h_udsg_Jeteta_num", 30, 0., 2.5); h_udsg_Jeteta_num->Sumw2();
  TH2F *h_udsg_JetpTeta_num=new TH2F("h_udsg_JetpTeta_num", "h_udsg_JetpTeta_num", 30, 0., 300., 30, 0., 2.5); h_udsg_JetpTeta_num->Sumw2();
  TH1F *h_c_JetpT_den=new TH1F("h_c_JetpT_den", "h_c_JetpT_den", 30, 0., 300.); h_c_JetpT_den->Sumw2();
  TH1F *h_c_Jeteta_den=new TH1F("h_c_Jeteta_den", "h_c_Jeteta_den", 30, 0., 2.5); h_c_Jeteta_den->Sumw2();
  TH2F *h_c_JetpTeta_den=new TH2F("h_c_JetpTeta_den", "h_c_JetpTeta_den", 30, 0., 300., 30, 0., 2.5); h_c_JetpTeta_den->Sumw2();
  TH1F *h_c_JetpT_num=new TH1F("h_c_JetpT_num", "h_c_JetpT_num", 30, 0., 300.); h_c_JetpT_num->Sumw2();
  TH1F *h_c_Jeteta_num=new TH1F("h_c_Jeteta_num", "h_c_Jeteta_num", 30, 0., 2.5); h_c_Jeteta_num->Sumw2();
  TH2F *h_c_JetpTeta_num=new TH2F("h_c_JetpTeta_num", "h_c_JetpTeta_num", 30, 0., 300., 30, 0., 2.5); h_c_JetpTeta_num->Sumw2();
  TH1F *h_b_JetpT_den=new TH1F("h_b_JetpT_den", "h_b_JetpT_den", 30, 0., 300.); h_b_JetpT_den->Sumw2();
  TH1F *h_b_Jeteta_den=new TH1F("h_b_Jeteta_den", "h_b_Jeteta_den", 30, 0., 2.5); h_b_Jeteta_den->Sumw2(); 
  TH2F *h_b_JetpTeta_den=new TH2F("h_b_JetpTeta_den", "h_b_JetpTeta_den", 30, 0., 300., 30, 0., 2.5); h_b_JetpTeta_den->Sumw2();
  TH1F *h_b_JetpT_num=new TH1F("h_b_JetpT_num", "h_b_JetpT_num", 30, 0., 300.); h_b_JetpT_num->Sumw2();
  TH1F *h_b_Jeteta_num=new TH1F("h_b_Jeteta_num", "h_b_Jeteta_num", 30, 0., 2.5); h_b_Jeteta_num->Sumw2();
  TH2F *h_b_JetpTeta_num=new TH2F("h_b_JetpTeta_num", "h_b_JetpTeta_num", 30, 0., 300., 30, 0., 2.5); h_b_JetpTeta_num->Sumw2();
  
  std::vector <std::string> samples;
  std::vector <std::string> PUWeights;
  
  samples.push_back(sample0); 
  samples.push_back(sample1);
  samples.push_back(sample2);
  PUWeights.push_back(PUWeight0);
  PUWeights.push_back(PUWeight1);
  PUWeights.push_back(PUWeight2);
  
  // Book variables
  int vType;
  bool triggerFlags[500];
  int nPV;
  int nhJets, naJets;
  int nJets;
  float hJetE[100], hJetpT[100], hJeteta[100], hJetphi[100], hJetCSV[100], hJetflavor[100], hJet_genPt[100];
  float aJetE[100], aJetpT[100], aJeteta[100], aJetphi[100], aJetCSV[100], aJetflavor[100], aJet_genPt[100];
  float jetE[100], jetpT[100], jeteta[100], jetphi[100], jetCSV[100], jetflavor[100], jet_genPt[100];
  
  for (unsigned int sample=0; sample<3; ++sample)
  {
    std::string inputfilename="/Users/souvik/HbbHbb/8TeV/OfficialStep2_"+samples.at(sample)+".root";
    TChain *tree=new TChain("tree");
    tree->Add(inputfilename.c_str());
    std::cout<<"Opened input file "<<inputfilename<<std::endl;
  
    std::string pufilename="/Users/souvik/HbbHbb/Analysis/"+PUWeights.at(sample);
    TFile *file_PUWeight=new TFile(pufilename.c_str());
    std::cout<<"Opened PU weight file = "<<pufilename<<std::endl;
    TH1F *h_PUWeight=(TH1F*)file_PUWeight->Get("h_PUWeight");
  
    // Retrieve variables
    tree->SetBranchAddress("Vtype", &(vType));
    tree->SetBranchAddress("triggerFlags", &(triggerFlags));
    tree->SetBranchAddress("nPVs", &(nPV));
    tree->SetBranchAddress("nhJets", &(nhJets));
    tree->SetBranchAddress("hJet_e", &(hJetE));
    tree->SetBranchAddress("hJet_pt", &(hJetpT));
    tree->SetBranchAddress("hJet_eta", &(hJeteta));
    tree->SetBranchAddress("hJet_phi", &(hJetphi));
    tree->SetBranchAddress("hJet_csv", &(hJetCSV));
    tree->SetBranchAddress("hJet_flavour", &(hJetflavor));
    tree->SetBranchAddress("hJet_genPt", &(hJet_genPt));
    tree->SetBranchAddress("naJets", &(naJets));
    tree->SetBranchAddress("aJet_e", &(aJetE));
    tree->SetBranchAddress("aJet_pt", &(aJetpT));
    tree->SetBranchAddress("aJet_eta", &(aJeteta));
    tree->SetBranchAddress("aJet_phi", &(aJetphi));
    tree->SetBranchAddress("aJet_csv", &(aJetCSV));
    tree->SetBranchAddress("aJet_flavour", &(aJetflavor));
    tree->SetBranchAddress("aJet_genPt", &(aJet_genPt));
    
    double scale_QCD;
    if (sample==0) scale_QCD=scale_QCDHT250To500;  
    if (sample==1) scale_QCD=scale_QCDHT500To1000;
    if (sample==2) scale_QCD=scale_QCDHT1000ToInf;
    std::cout<<"Scale QCD = "<<scale_QCD<<std::endl;
  
    // Loop over events
    int nEvents=tree->GetEntries();
    for (int i=0; i<nEvents; ++i)
    {
      tree->GetEvent(i);
      
      double weightPU=h_PUWeight->GetBinContent(h_PUWeight->FindBin(nPV));
      weightPU=weightPU*scale_QCD;
      
      // Collate hJets and aJets into Jets
      nJets=nhJets+naJets;
      for (int j=0; j<nhJets; ++j)
      {
        jetE[j]=hJetE[j];
        jetpT[j]=hJetpT[j];
        jeteta[j]=hJeteta[j];
        jetphi[j]=hJetphi[j];
        jetCSV[j]=hJetCSV[j];
        jetflavor[j]=hJetflavor[j];
        jet_genPt[j]=hJet_genPt[j];
      }
      for (int j=0; j<naJets; ++j)
      {
        jetE[j+nhJets]=aJetE[j];
        jetpT[j+nhJets]=aJetpT[j];
        jeteta[j+nhJets]=aJeteta[j];
        jetphi[j+nhJets]=aJetphi[j];
        jetCSV[j+nhJets]=aJetCSV[j];
        jetflavor[j+nhJets]=aJetflavor[j];
        jet_genPt[j+nhJets]=aJet_genPt[j];
      }
      
      int nCJets=0;
      for (int j=0; j<nJets; ++j)
      {
        if (fabs(jeteta[j])<jeteta_cut && jetpT[j]>jetpT_cut) ++nCJets;
      }
    
      if (vType==5 && nCJets>3)
      {
        for (int j=0; j<nJets; ++j)
        {
          if (fabs(jetflavor[j])==1 || fabs(jetflavor[j])==2 || fabs(jetflavor[j])==3 || fabs(jetflavor[j])==21) // d, u, s, g
          {
            if (fabs(jeteta[j])<jeteta_cut) h_udsg_JetpT_den->Fill(jet_genPt[j], weightPU); // h_udsg_JetpT_den->Fill(jetpT[j]);
            if (jetpT[j]>jetpT_cut) h_udsg_Jeteta_den->Fill(jeteta[j], weightPU);
            h_udsg_JetpTeta_den->Fill(jetpT[j], jeteta[j], weightPU);
            if (jetCSV[j]>cut_CSV)
            {
              if (fabs(jeteta[j])<jeteta_cut) h_udsg_JetpT_num->Fill(jet_genPt[j], weightPU); // h_udsg_JetpT_num->Fill(jetpT[j]);
              if (jetpT[j]>jetpT_cut) h_udsg_Jeteta_num->Fill(jeteta[j], weightPU);
              h_udsg_JetpTeta_num->Fill(jetpT[j], jeteta[j], weightPU);
            }
          }
          
          if (fabs(jetflavor[j])==4) // c
          {
            if (fabs(jeteta[j])<jeteta_cut) h_c_JetpT_den->Fill(jetpT[j], weightPU);
            if (jetpT[j]>jetpT_cut) h_c_Jeteta_den->Fill(jeteta[j], weightPU);
            h_c_JetpTeta_den->Fill(jetpT[j], jeteta[j], weightPU);
            if (jetCSV[j]>cut_CSV)
            {
              if (fabs(jeteta[j])<jeteta_cut) h_c_JetpT_num->Fill(jetpT[j], weightPU);
              if (jetpT[j]>jetpT_cut) h_c_Jeteta_num->Fill(jeteta[j], weightPU);
              h_c_JetpTeta_num->Fill(jetpT[j], jeteta[j], weightPU);
            }
          }
        
          if (fabs(jetflavor[j])==5) // b
          {
            if (fabs(jeteta[j])<jeteta_cut) h_b_JetpT_den->Fill(jetpT[j], weightPU);
            if (jetpT[j]>jetpT_cut) h_b_Jeteta_den->Fill(jeteta[j], weightPU);
            h_b_JetpTeta_den->Fill(jetpT[j], jeteta[j], weightPU);
            if (jetCSV[j]>cut_CSV)
            {
              if (fabs(jeteta[j])<jeteta_cut) h_b_JetpT_num->Fill(jetpT[j], weightPU);
              if (jetpT[j]>jetpT_cut) h_b_Jeteta_num->Fill(jeteta[j], weightPU);
              h_b_JetpTeta_num->Fill(jetpT[j], jeteta[j], weightPU);
            }
          }
        } // For loop over jets in the event
      } // If vType==5 && nCJets > 3
    } // For loop over events in this sample
  } // For loop over the 3 QCD samples
  
  std::string outputfilename="CSVFlavorEfficiencies_"+csvTag+".root";
  TFile *tFile=new TFile(outputfilename.c_str(), "RECREATE");
  h_udsg_JetpT_den->Write();
  h_udsg_Jeteta_den->Write();
  h_udsg_JetpTeta_den->Write();
  h_udsg_JetpT_num->Write();
  h_udsg_Jeteta_num->Write();
  h_udsg_JetpTeta_num->Write();
  h_c_JetpT_den->Write();
  h_c_Jeteta_den->Write();
  h_c_JetpTeta_den->Write();
  h_c_JetpT_num->Write();
  h_c_Jeteta_num->Write();
  h_c_JetpTeta_num->Write();
  h_b_JetpT_den->Write();
  h_b_Jeteta_den->Write();
  h_b_JetpTeta_den->Write();
  h_b_JetpT_num->Write();
  h_b_Jeteta_num->Write();
  h_b_JetpTeta_num->Write();
  tFile->Write();
  tFile->Close();
  
}
    
    
    
    
    
    
    
    
    
