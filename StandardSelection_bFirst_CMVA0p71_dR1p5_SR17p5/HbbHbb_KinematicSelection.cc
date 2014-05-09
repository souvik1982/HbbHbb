#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "/Users/souvik/HbbHbb/Analysis/bJetRegression/HelperFunctions.h"

double pi=3.14159265358979;

double jetpT_cut=40.;
double jeteta_cut=2.5;
double H_mass=160.0;
double mH_diff_cut=50.;
double mH_mean_cut=15.;

double bTagCSV_tightCut=0.898;
double bTagCSV_mediumCut=0.679;
double bTagCSV_looseCut=0.244;
double bTagCSV_noCut=0.;

typedef struct
{
  float et;
  float sumet;
  float sig;
  float phi;
} METInfo;

typedef struct
{
  float CSV;
  float E;
  float pT;
  float eta;
  float phi;
} JetInfo;

typedef struct 
{
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid;
} GenParticleInfo;

typedef struct
{
  int run;
  int lumi;
  int event;
  int json;
} EventInfo;

JetInfo jet_CSV1={-1,-1,-1,-1,-1};
JetInfo jet_CSV2={-1,-1,-1,-1,-1};
JetInfo jet_CSV3={-1,-1,-1,-1,-1};
JetInfo jet_CSV4={-1,-1,-1,-1,-1};
JetInfo jet_CSV5={-1,-1,-1,-1,-1};

typedef std::map<double, int> JetList;

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double E)
{
  TLorentzVector jet_p4;
  jet_p4.SetPtEtaPhiE(pT, eta, phi, E);
  return jet_p4;
}

int HbbHbb_KinematicSelection(std::string dir, std::string sample, 
                              std::string sigmaJECUnc_string="JEC", 
                              std::string sigmaJERUnc_string="JER")
{
  
  std::string inputfilename="/Users/souvik/HbbHbb/8TeV/"+dir+"/OfficialStep2_"+sample+".root";
  TChain *tree=new TChain("tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;
  
  double sigmaJECUnc; // -1, 0, +1
  double sigmaJERUnc; // -1, 0, +1
  if (sigmaJECUnc_string=="JEC") sigmaJECUnc=0;
  else if (sigmaJECUnc_string=="JECp1") sigmaJECUnc=1;
  else if (sigmaJECUnc_string=="JECm1") sigmaJECUnc=-1;
  
  if (sigmaJERUnc_string=="JER") sigmaJERUnc=0;
  else if (sigmaJERUnc_string=="JERp1") sigmaJERUnc=1;
  else if (sigmaJERUnc_string=="JERm1") sigmaJERUnc=-1;
  
  // Book variables
  EventInfo EVENT;
  int vType;
  bool triggerFlags[500];
  int nPV;
  int nhJets, naJets;
  int nJets, nCJets;
  float hJetE[100], hJetpT[100], hJeteta[100], hJetphi[100], hJetCSV[100], hJetflavor[100], hJetpTRaw[100], hJet_ptLeadTrack[100], hJet_vtx3dL[100], hJet_vtx3deL[100], hJet_vtxMass[100], hJet_vtxPt[100], hJet_cef[100], hJet_nconstituents[100], hJet_JECUnc[100], hJet_genpT[100];
  float aJetE[100], aJetpT[100], aJeteta[100], aJetphi[100], aJetCSV[100], aJetflavor[100], aJetpTRaw[100], aJet_ptLeadTrack[100], aJet_vtx3dL[100], aJet_vtx3deL[100], aJet_vtxMass[100], aJet_vtxPt[100], aJet_cef[100], aJet_nconstituents[100], aJet_JECUnc[100], aJet_genpT[100];
  float jetE[100], jetpT[100], jeteta[100], jetphi[100], jetCSV[100], jetflavor[100], jetpTRaw[100], jet_ptLeadTrack[100], jet_vtx3dL[100], jet_vtx3deL[100], jet_vtxMass[100], jet_vtxPt[100], jet_cef[100], jet_nconstituents[100], jet_JECUnc[100], jet_genpT[100];
  METInfo metObj;
  float met, metSig, ht;
  float weightPU;
  GenParticleInfo genX, genH1, genH2;
  
  // Retrieve variables
  tree->SetBranchAddress("EVENT", &EVENT);
  tree->SetBranchAddress("Vtype", &(vType));
  tree->SetBranchAddress("triggerFlags", &(triggerFlags));
  tree->SetBranchAddress("nPVs", &(nPV));
  tree->SetBranchAddress("nhJets", &(nhJets));
  tree->SetBranchAddress("hJet_e", &(hJetE)); 
  tree->SetBranchAddress("hJet_pt", &(hJetpT));
  tree->SetBranchAddress("hJet_eta", &(hJeteta));
  tree->SetBranchAddress("hJet_phi", &(hJetphi));
  // tree->SetBranchAddress("hJet_csv_nominal", &(hJetCSV));
  tree->SetBranchAddress("hJet_cmva", &(hJetCSV));
  tree->SetBranchAddress("hJet_JECUnc", &(hJet_JECUnc));
  tree->SetBranchAddress("hJet_genPt", &(hJet_genpT));
  tree->SetBranchAddress("naJets", &(naJets));
  tree->SetBranchAddress("aJet_e", &(aJetE));
  tree->SetBranchAddress("aJet_pt", &(aJetpT));
  tree->SetBranchAddress("aJet_eta", &(aJeteta));
  tree->SetBranchAddress("aJet_phi", &(aJetphi));
  // tree->SetBranchAddress("aJet_csv_nominal", &(aJetCSV));
  tree->SetBranchAddress("aJet_cmva", &(aJetCSV));
  tree->SetBranchAddress("aJet_JECUnc", &(aJet_JECUnc));
  tree->SetBranchAddress("aJet_genPt", &(aJet_genpT));
  tree->SetBranchAddress("MET", &(metObj));
  tree->SetBranchAddress("genX", &(genX));
  tree->SetBranchAddress("genH1", &(genH1));
  tree->SetBranchAddress("genH2", &(genH2));
  tree->SetBranchAddress("PUweight", &(weightPU));
  
  TH1F *h_nJets=new TH1F("h_nJets", "# Cleaned PAT Jets; n", 20, 0., 20.);
  TH1F *h_nPV=new TH1F("h_nPV", "# of Primary Vertices; nPV", 51, 0., 50.); h_nPV->Sumw2();
  TH1F *h_nPV_weighted=new TH1F("h_nPV_weighted", "# of Primary Vertices after Reweighting; nPV", 51, 0., 50.); h_nPV_weighted->Sumw2();
  
  TH1F *h_JetpT1=new TH1F("h_JetpT1", "JetpT1; Jet pT 1 (GeV); Events", 50, 0., 800.);
  TH1F *h_JetpT2=new TH1F("h_JetpT2", "JetpT2; Jet pT 2 (GeV); Events", 50, 0., 500.);
  TH1F *h_JetpT3=new TH1F("h_JetpT3", "JetpT3; Jet pT 3 (GeV); Events", 50, 0., 350.);
  TH1F *h_JetpT4=new TH1F("h_JetpT4", "JetpT4; Jet pT 4 (GeV); Events", 50, 0., 250.);
  
  TH1F *h_MET=new TH1F("h_MET", "MET; MET (GeV)", 100, 0, 200.);
  TH1F *h_MET_sig=new TH1F("h_MET_sig", "MET Significance; Sig", 20, 0., 20.);
  TH1F *h_HT=new TH1F("h_HT", "HT Distribution; HT (GeV); Events", 100, 0., 3000.);
  
  TH1F *h_H1_mass=new TH1F("h_H1_mass", "H1 mass; mass (GeV)", 50, 50., 200.);
  TH1F *h_H1_pT=new TH1F("h_H1_pT", "H1 p_{T}; p_{T} (GeV/c)", 50, 0., 800.);
  TH1F *h_H2_mass=new TH1F("h_H2_mass", "H2 mass; mass (GeV)", 50, 50., 200.);
  TH1F *h_H2_pT=new TH1F("h_H2_pT", "H2 p_{T}; p_{T} (GeV/c)", 50, 0., 800.);
  
  TH1F *h_HH_mass_diff_cand=new TH1F("h_HH_mass_diff_cand", "|#Deltam| between Higgs masses - all candidates", 50, 0., 200.);
  TH1F *h_HH_massNorm_diff=new TH1F("h_HH_massNorm_diff", "|#Deltam| between Higgs masses", 50, 0., 2.);
  
  TH1F *h_genX_mass=new TH1F("h_genX_mass", "Generator X mass; m_{X} GeV", 100, 0., 1000.);
  TH1F *h_genH1_mass=new TH1F("h_genH1_mass", "Generator H1 mass; m_{X} GeV", 100, 0., 1000.);
  TH1F *h_genH2_mass=new TH1F("h_genH2_mass", "Generator H2 mass; m_{X} GeV", 100, 0., 1000.);
  
  TH1F *h_trig=new TH1F("h_trig", "Trigger bits", 30, 0., 30.);
  
  TH1F *h_Cuts=new TH1F("h_Cuts", "Cut flow", 16, 0, 16);
  TAxis *a_Cuts=h_Cuts->GetXaxis();
  a_Cuts->SetBinLabel(2, "Step2");
  a_Cuts->SetBinLabel(4, "Trigger");
  a_Cuts->SetBinLabel(6, "Vtype");
  a_Cuts->SetBinLabel(8, "nCJets>3");
  a_Cuts->SetBinLabel(10, "HH Cand");
  a_Cuts->SetBinLabel(12, "massDiff<50");
  a_Cuts->SetBinLabel(14, "b-tagging");
  a_Cuts->SetBinLabel(16, "SR");
  
  std::string outfilename="OfficialStep2_KinematicallySelected_"+sample+".root";
  TFile *outfile=new TFile(outfilename.c_str(), "recreate");
  TTree *outtree=tree->CloneTree(0);
  int H1jet1_i, H1jet2_i;
  int H2jet1_i, H2jet2_i;
  int samepTHiggs=-1;
  outtree->Branch("H1jet1_i", &H1jet1_i, "H1jet1_i/I");
  outtree->Branch("H1jet2_i", &H1jet2_i, "H1jet2_i/I");
  outtree->Branch("H2jet1_i", &H2jet1_i, "H2jet1_i/I");
  outtree->Branch("H2jet2_i", &H2jet2_i, "H2jet2_i/I");
  outtree->Branch("samepTHiggs", &samepTHiggs, "samepTHiggs/I");
  
  bTagCSV_mediumCut=0.71;
  
  // Loop over events
  int nEvents=tree->GetEntries();
  double nCut0=0, nCut1=0, nCut2=0, nCut3=0, nCut4=0, nCut5=0, nCut6=0, nCut7=0;
  for (int i=0; i<nEvents; ++i)
  {
    if((207883<=EVENT.run && EVENT.run<=208307)) continue;
    ++nCut0;
    tree->GetEvent(i);
    
    h_nPV->Fill(nPV);
    h_nPV_weighted->Fill(nPV, weightPU);
    
    // Fill generator level info
    TLorentzVector h1, h2;
    h1.SetPtEtaPhiM(genH1.pt, genH1.eta, genH1.phi, genH1.mass);
    h2.SetPtEtaPhiM(genH2.pt, genH2.eta, genH2.phi, genH2.mass);
    h_genX_mass->Fill((h1+h2).M());
    // h_genX_mass->Fill(genX.mass);
    h_genH1_mass->Fill(genH1.mass);
    h_genH2_mass->Fill(genH2.mass);
    
    // Collate hJets and aJets into Jets
    nJets=nhJets+naJets;
    for (int j=0; j<nhJets; ++j)
    {
      jetpT[j]=smear_pt_resErr(hJetpT[j], hJet_genpT[j], hJeteta[j], sigmaJERUnc)+sigmaJECUnc*hJet_JECUnc[j];
      jetE[j]=hJetE[j]*jetpT[j]/hJetpT[j];
      jeteta[j]=hJeteta[j];
      jetphi[j]=hJetphi[j];
      jetCSV[j]=hJetCSV[j];
      jetflavor[j]=hJetflavor[j];
      jetpTRaw[j]=hJetpTRaw[j];
      jet_ptLeadTrack[j]=hJet_ptLeadTrack[j];
      jet_vtx3dL[j]=hJet_vtx3dL[j];
      jet_vtx3deL[j]=hJet_vtx3deL[j];
      jet_vtxMass[j]=hJet_vtxMass[j];
      jet_vtxPt[j]=hJet_vtxPt[j];
      jet_cef[j]=hJet_cef[j];
      jet_nconstituents[j]=hJet_nconstituents[j];
      jet_JECUnc[j]=hJet_JECUnc[j];
      jet_genpT[j]=hJet_genpT[j];
    }
    for (int j=0; j<naJets; ++j)
    {
      jetpT[j+nhJets]=smear_pt_resErr(aJetpT[j], aJet_genpT[j], aJeteta[j], sigmaJERUnc)+sigmaJECUnc*aJet_JECUnc[j];
      jetE[j+nhJets]=aJetE[j]*jetpT[j+nhJets]/aJetpT[j];
      jeteta[j+nhJets]=aJeteta[j];
      jetphi[j+nhJets]=aJetphi[j];
      jetCSV[j+nhJets]=aJetCSV[j];
      jetflavor[j+nhJets]=aJetflavor[j];
      jetpTRaw[j+nhJets]=aJetpTRaw[j];
      jet_ptLeadTrack[j+nhJets]=aJet_ptLeadTrack[j];
      jet_vtx3dL[j+nhJets]=aJet_vtx3dL[j];
      jet_vtx3deL[j+nhJets]=aJet_vtx3deL[j];
      jet_vtxMass[j+nhJets]=aJet_vtxMass[j];
      jet_vtxPt[j+nhJets]=aJet_vtxPt[j];
      jet_cef[j+nhJets]=aJet_cef[j];
      jet_nconstituents[j+nhJets]=aJet_nconstituents[j];
      jet_JECUnc[j+nhJets]=aJet_JECUnc[j];
      jet_genpT[j+nhJets]=aJet_genpT[j];
    }
    
    JetList jetList_pT;
    nCJets=0;
    ht=0;
    for (int j=0; j<nJets; ++j)
    {
      if (fabs(jeteta[j])<jeteta_cut) 
      {
        jetList_pT[jetpT[j]]=j;
        if (jetpT[j]>jetpT_cut && jetCSV[j]>bTagCSV_mediumCut)
	      {
          ++nCJets;
          ht+=jetpT[j];
        }
      }
    }
    
    // Trigger bits
    /*
    for (unsigned int trig=0; trig<30; ++trig)
    {
      if (triggerFlags[trig]) h_trig->Fill(trig);
    }
    */
    // Analysis begins here
    // if (triggerFlags[54] && !(207883<=EVENT.run && EVENT.run<=208307)) // For data
    if (triggerFlags[54])
    // if (triggerFlags[0]) // For signal MC
    {
      ++nCut1;
      if (vType==4 || vType==8 || vType==9 || vType==10)
      {
        ++nCut2;
        
        h_nJets->Fill(nCJets, weightPU);
        
        met=metObj.et;
        metSig=metObj.sig;
        // ht=metObj.sumet;
        h_MET->Fill(met);
        h_MET_sig->Fill(metSig);
        h_HT->Fill(ht);
        
        JetList::reverse_iterator iJet=jetList_pT.rbegin();
        if (iJet!=jetList_pT.rend())
        {
          h_JetpT1->Fill(iJet->first, weightPU);
          ++iJet; 
          if (iJet!=jetList_pT.rend())
          {
            h_JetpT2->Fill(iJet->first, weightPU);
            ++iJet; 
            if (iJet!=jetList_pT.rend()) 
            {
              h_JetpT3->Fill(iJet->first, weightPU);
              ++iJet; 
              if (iJet!=jetList_pT.rend()) 
              {
                h_JetpT4->Fill(iJet->first, weightPU);
              }
            }
          }
        }
      
        if (nCJets>3) // nCJets > 3
        {
          ++nCut3;
      
          bool foundHH=false;
          // int H1jet1_i, H1jet2_i;
          // int H2jet1_i, H2jet2_i;
          double m_diff_old=50.;
          for (int j=0; j<nJets; ++j)
          {
            TLorentzVector jet1_p4, jet2_p4, jet3_p4, jet4_p4;
            if (fabs(jeteta[j])<jeteta_cut && jetpT[j]>jetpT_cut && jetCSV[j]>bTagCSV_mediumCut)
            {
              jet1_p4=fillTLorentzVector(jetpT[j], jeteta[j], jetphi[j], jetE[j]);  
              for (int k=0; k<nJets; ++k)
              {
                if (k!=j && fabs(jeteta[k])<jeteta_cut && jetpT[k]>jetpT_cut && jetCSV[k]>bTagCSV_mediumCut)
                {
                  jet2_p4=fillTLorentzVector(jetpT[k], jeteta[k], jetphi[k], jetE[k]);
                  for (int l=0; l<nJets; ++l)
                  {
                    if (l!=k && l!=j && fabs(jeteta[l])<jeteta_cut && jetpT[l]>jetpT_cut && jetCSV[l]>bTagCSV_mediumCut)
                    {
                      jet3_p4=fillTLorentzVector(jetpT[l], jeteta[l], jetphi[l], jetE[l]);
                      for (int m=0; m<nJets; ++m)
                      {
                        if (m!=l && m!=k && m!=j && fabs(jeteta[m])<jeteta_cut && jetpT[m]>jetpT_cut && jetCSV[m]>bTagCSV_mediumCut)
                        {
                          jet4_p4=fillTLorentzVector(jetpT[m], jeteta[m], jetphi[m], jetE[m]);
                       
                          TLorentzVector diJet1_p4=jet1_p4+jet2_p4;
                          TLorentzVector diJet2_p4=jet3_p4+jet4_p4;
                          
                          double deltaR1=jet1_p4.DeltaR(jet2_p4);
                          double deltaR2=jet3_p4.DeltaR(jet4_p4);
                          // double etaSign=diJet1_p4.Eta()*diJet2_p4.Eta();
                        
                          double m_diff=fabs(diJet1_p4.M()-diJet2_p4.M());
                          if (m_diff<m_diff_old && deltaR1<1.5 && deltaR2<1.5) // && etaSign<0.)
                          {
                            H1jet1_i=j;
                            H1jet2_i=k;
                            H2jet1_i=l;
                            H2jet2_i=m;
                            m_diff_old=m_diff;
                            foundHH=true;
                            h_HH_mass_diff_cand->Fill(m_diff, weightPU);
                            h_HH_massNorm_diff->Fill(2.*m_diff/(diJet1_p4.M()+diJet2_p4.M()), weightPU);
                          }
                        } // Conditions on 4th jet
                      } // Loop over 4th jet
                    } // Conditions on 3rd jet
                  } // Loop over 3rd jet
                } // Conditions on 2nd jet
              } // Loop over 2nd jet
            } // Conditions on 1st jet
          } // Loop over 1st jet
        
          if (foundHH)
          {
            ++nCut4;
            
            TLorentzVector jet1_p4=fillTLorentzVector(jetpT[H1jet1_i], jeteta[H1jet1_i], jetphi[H1jet1_i], jetE[H1jet1_i]);
            TLorentzVector jet2_p4=fillTLorentzVector(jetpT[H1jet2_i], jeteta[H1jet2_i], jetphi[H1jet2_i], jetE[H1jet2_i]);    
            TLorentzVector jet3_p4=fillTLorentzVector(jetpT[H2jet1_i], jeteta[H2jet1_i], jetphi[H2jet1_i], jetE[H2jet1_i]);    
            TLorentzVector jet4_p4=fillTLorentzVector(jetpT[H2jet2_i], jeteta[H2jet2_i], jetphi[H2jet2_i], jetE[H2jet2_i]);
            
            // Randomization or ordering of which Higgs is which
            if (int((jet1_p4+jet2_p4).Pt()*100.) % 2 == 1) {swap(H1jet1_i, H2jet1_i); swap(H1jet2_i, H2jet2_i);} // swap if H pT is odd in second decimal place
            
            jet1_p4=fillTLorentzVector(jetpT[H1jet1_i], jeteta[H1jet1_i], jetphi[H1jet1_i], jetE[H1jet1_i]);
            jet2_p4=fillTLorentzVector(jetpT[H1jet2_i], jeteta[H1jet2_i], jetphi[H1jet2_i], jetE[H1jet2_i]); 
            jet3_p4=fillTLorentzVector(jetpT[H2jet1_i], jeteta[H2jet1_i], jetphi[H2jet1_i], jetE[H2jet1_i]); 
            jet4_p4=fillTLorentzVector(jetpT[H2jet2_i], jeteta[H2jet2_i], jetphi[H2jet2_i], jetE[H2jet2_i]);
            
            TLorentzVector H1_p4=jet1_p4+jet2_p4;
            TLorentzVector H2_p4=jet3_p4+jet4_p4;
            
            h_H1_mass->Fill(H1_p4.M(), weightPU);
            h_H1_pT->Fill(H1_p4.Pt(), weightPU);
            h_H2_mass->Fill(H2_p4.M(), weightPU);
            h_H2_pT->Fill(H2_p4.Pt(), weightPU);
            
            // Check if two highest tagged jets are in the same Higgs
            if ( ((jetpT[H1jet1_i]>jetpT[H2jet1_i] && jetpT[H1jet1_i]>jetpT[H2jet2_i]) &&
                  (jetpT[H1jet2_i]>jetpT[H2jet1_i] && jetpT[H1jet2_i]>jetpT[H2jet2_i])) ||
                 ((jetpT[H2jet1_i]>jetpT[H1jet1_i] && jetpT[H2jet1_i]>jetpT[H1jet2_i]) &&
                  (jetpT[H2jet2_i]>jetpT[H1jet1_i] && jetpT[H2jet2_i]>jetpT[H1jet2_i])) ) samepTHiggs=1;
            else samepTHiggs=0;
            
            // Write out tree
            outtree->Fill();
            
          } // found HH
        } // nCJets>3
      } // vType==5
    } // Trigger
  } // Event loop
  
  h_Cuts->Fill(1, nCut0);
  h_Cuts->Fill(3, nCut1);
  h_Cuts->Fill(5, nCut2);
  h_Cuts->Fill(7, nCut3);
  h_Cuts->Fill(9, nCut4);
  h_Cuts->Fill(11, nCut5);
  h_Cuts->Fill(13, nCut6);
  h_Cuts->Fill(15, nCut7);
  
  outtree->Write();
  outfile->Close();
  std::cout<<"Wrote output file "<<outfilename<<std::endl;
  
  TH1F *h_CountWithPU;
  // Pass the countWithPU forward
  // TFile *file=new TFile(inputfilename.c_str());
  // h_CountWithPU=(TH1F*)file->Get("CountWithPU");
                  
  std::string histfilename="Histograms_"+sample+".root";
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  // h_CountWithPU->Write();
  h_nJets->Write();
  h_nPV->Write();
  h_nPV_weighted->Write();
  h_JetpT1->Write();
  h_JetpT2->Write();
  h_JetpT3->Write();
  h_JetpT4->Write();
  h_MET->Write();
  h_MET_sig->Write();
  h_HT->Write();
  h_H1_mass->Write();
  h_H1_pT->Write();
  h_H2_mass->Write();
  h_H2_pT->Write();
  h_HH_mass_diff_cand->Write();
  h_HH_massNorm_diff->Write();
  h_genX_mass->Write();
  h_genH1_mass->Write();
  h_genH2_mass->Write();
  h_Cuts->Write();
  h_trig->Write();
  tFile->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  
  /*
  std::cout<<std::endl<<"=== Inference valid for QCD sample ==="<<std::endl;
  std::cout<<"n4Tags = "<<n4Tags<<std::endl;
  std::cout<<"n4Flavored = "<<n4Flavored<<std::endl;
  std::cout<<"n3Flavored = "<<n3Flavored<<std::endl;
  std::cout<<"n2Flavored = "<<n2Flavored<<std::endl;
  std::cout<<"n1Flavored = "<<n1Flavored<<std::endl;
  std::cout<<"n0Flavored = "<<n0Flavored<<std::endl;
  std::cout<<"n4Flavored/n4Tags = "<<n4Flavored/n4Tags<<std::endl;
  std::cout<<"n3Flavored/n4Tags = "<<n3Flavored/n4Tags<<std::endl;
  std::cout<<"n2Flavored/n4Tags = "<<n2Flavored/n4Tags<<std::endl;
  std::cout<<"n1Flavored/n4Tags = "<<n1Flavored/n4Tags<<std::endl;
  std::cout<<"n0Flavored/n4Tags = "<<n0Flavored/n4Tags<<std::endl;
  std::cout<<"=== x ==="<<std::endl;
  std::cout<<"nCheck = "<<nCheck<<std::endl;
  std::cout<<"Done."<<std::endl;
  */
  
  std::cout<<"=== Cut Efficiencies === "<<std::endl;
  std::cout<<"Number of events at the end of step 2 = "<<nCut0<<std::endl;
  std::cout<<"Number of events after trigger = "<<nCut1<<std::endl;
  std::cout<<"Number of events after Vtype==5 = "<<nCut2<<std::endl;
  std::cout<<"Number of events after nCJets>3 = "<<nCut3<<std::endl;
  std::cout<<"Number of events after finding HH kinematic candidate = "<<nCut4<<std::endl;
  std::cout<<"========================"<<std::endl;
  
  return 0;
}
            
      
  
  
  
  
