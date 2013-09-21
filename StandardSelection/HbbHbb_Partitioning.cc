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
#include <TCanvas.h>
#include <TProfile.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <TMath.h>

#include "/Users/souvik/HbbHbb/Analysis/bJetRegression/HelperFunctions.h"

#include "/Users/souvik/HbbHbb/Analysis/kinFit4b.h"

double sigmaJECUnc=0; // (-1, 0, 1)
double sigmaJERUnc=0; // (-1, 0, 1)

double pi=3.14159265358979;

double bTagCSV_tightCut=0.898;
double bTagCSV_mediumCut=0.679;
double bTagCSV_looseCut=0.244;
double bTagCSV_noCut=0.;

double jetpT_cut=40.;
double jeteta_cut=2.5;
double H_mass=125.0;
double mH_diff_cut=50.;
double mH_mean_cut=15.;
double HpT_cut=0.;
double jetCSV_cut=bTagCSV_mediumCut;

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

typedef std::map<double, int> JetList;

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double E)
{
  TLorentzVector jet_p4;
  jet_p4.SetPtEtaPhiE(pT, eta, phi, E);
  return jet_p4;
}

TLorentzVector fillTLorentzVector(GenParticleInfo gen)
{
  TLorentzVector gen_p4;
  gen_p4.SetPtEtaPhiM(gen.pt, gen.eta, gen.phi, gen.mass);
  return gen_p4;
}

/*
double deltaPhi(double phi1, double phi2)
{
  double dphi=phi1-phi2;
  if (dphi<-pi) dphi=2*pi+dphi;
  if (dphi>pi) dphi=2*pi-dphi;
  return dphi;
}
*/
double azimuthalAngle(TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, TLorentzVector v4)
{
  TLorentzVector v12=v1+v2;
  TVector3 v12_boost=v12.BoostVector();
  TLorentzVector v1_unboosted=v1; v1_unboosted.Boost(-v12_boost);
  TLorentzVector v2_unboosted=v2; v2_unboosted.Boost(-v12_boost);
  TVector3 v1v2Plane=v1_unboosted.Vect().Cross(v12_boost);
  
  TLorentzVector v34=v3+v4;
  TVector3 v34_boost=v34.BoostVector();
  TLorentzVector v3_unboosted=v3; v3_unboosted.Boost(-v34_boost);
  TLorentzVector v4_unboosted=v4; v4_unboosted.Boost(-v34_boost);
  TVector3 v3v4Plane=v3_unboosted.Vect().Cross(v34_boost);
  
  double angle=v1v2Plane.Angle(v3v4Plane);
  return angle;
}

TH2F *h_udsg_JetpTeta_den, *h_udsg_JetpTeta_num;
TH2F *h_c_JetpTeta_den, *h_c_JetpTeta_num;
TH2F *h_b_JetpTeta_den, *h_b_JetpTeta_num;

double CSVEfficiency(std::string j, double jetpT, double jeteta, double jetCSV)
{
  double weight, num, den;
  
  if (j=="j")
  {
    if (jetCSV>jetCSV_cut) weight=1.; else weight=0;
  }
  else if (j=="Q")
  {
    num=h_udsg_JetpTeta_num->GetBinContent(h_udsg_JetpTeta_num->FindBin(jetpT, jeteta));
    den=h_udsg_JetpTeta_den->GetBinContent(h_udsg_JetpTeta_den->FindBin(jetpT, jeteta));
    weight=num/den;
  }
  else if (j=="C")
  {
    num=h_c_JetpTeta_num->GetBinContent(h_c_JetpTeta_num->FindBin(jetpT, jeteta));
    den=h_c_JetpTeta_den->GetBinContent(h_c_JetpTeta_den->FindBin(jetpT, jeteta));
    weight=num/den;
  }
  else if (j=="B")
  {
    num=h_b_JetpTeta_num->GetBinContent(h_b_JetpTeta_num->FindBin(jetpT, jeteta));
    den=h_b_JetpTeta_den->GetBinContent(h_b_JetpTeta_den->FindBin(jetpT, jeteta));
    weight=num/den;
  }
  
  return weight;
}

bool returnRandom()
{
  double rand1=(double)rand()/double(RAND_MAX);
  bool returnValue;
  if (rand1<0.5) returnValue=false; else returnValue=true;
  return returnValue;
}

void fillStringsRandomly(std::string templ, std::string &j1, std::string &j2, std::string &j3, std::string &j4)
{
  if (templ=="QQjj") {j1="Q"; j2="Q"; j3="j"; j4="j";}
  if (templ=="QCjj") {j1="Q"; j2="C"; j3="j"; j4="j";}
  if (templ=="QBjj") {j1="Q"; j2="B"; j3="j"; j4="j";}
  if (templ=="CCjj") {j1="C"; j2="C"; j3="j"; j4="j";}
  if (templ=="CBjj") {j1="C"; j2="B"; j3="j"; j4="j";}
  if (templ=="BBjj") {j1="B"; j2="B"; j3="j"; j4="j";}
  if (templ=="QjQj") {j1="Q"; j2="j"; j3="Q"; j4="j";}
  if (templ=="QjCj") {j1="Q"; j2="j"; j3="C"; j4="j";}
  if (templ=="QjBj") {j1="Q"; j2="j"; j3="B"; j4="j";}
  if (templ=="CjCj") {j1="C"; j2="j"; j3="C"; j4="j";}
  if (templ=="CjBj") {j1="C"; j2="j"; j3="B"; j4="j";}
  if (templ=="BjBj") {j1="B"; j2="j"; j3="B"; j4="j";}
  
  // Now randomize between Higgs and between jets
  if (returnRandom())
  {
    swap(j1, j3);
    swap(j2, j4);
  }
  if (returnRandom())
  {
    swap(j1, j2);
  }
  if (returnRandom())
  {
    swap(j3, j4);
  }
}

std::string returnFlavorString(int flavor)
{
  std::string flav="j";
  if (fabs(flavor)==1 || fabs(flavor)==2 || fabs(flavor)==3 || fabs(flavor)==21) flav="Q";
  if (fabs(flavor)==4) flav="C";
  if (fabs(flavor)==5) flav="B";
  return flav;
}

bool truthMatch(TLorentzVector reco_p4, TLorentzVector gen_p4)
{
  bool returnValue=false;
  double match_cut=0.7;
  double dR=reco_p4.DeltaR(gen_p4); 
  if (dR<match_cut) returnValue=true;
  return returnValue;
}

std::string truthMatching(TLorentzVector jet1_p4, TLorentzVector jet2_p4,
                          TLorentzVector bH1_gen_p4, TLorentzVector bbH1_gen_p4, TLorentzVector bH2_gen_p4, TLorentzVector bbH2_gen_p4)
{
  std::string returnString="noMatch";
  
  if ((truthMatch(jet1_p4, bH1_gen_p4) && truthMatch(jet2_p4, bH2_gen_p4)) ||
	    (truthMatch(jet1_p4, bH1_gen_p4) && truthMatch(jet2_p4, bbH2_gen_p4)) ||
	    (truthMatch(jet1_p4, bbH1_gen_p4) && truthMatch(jet2_p4, bH2_gen_p4)) ||
	    (truthMatch(jet1_p4, bbH1_gen_p4) && truthMatch(jet2_p4, bbH2_gen_p4)) ||
	    (truthMatch(jet1_p4, bH2_gen_p4) && truthMatch(jet2_p4, bH1_gen_p4)) ||
	    (truthMatch(jet1_p4, bH2_gen_p4) && truthMatch(jet2_p4, bbH1_gen_p4)) ||
	    (truthMatch(jet1_p4, bbH2_gen_p4) && truthMatch(jet2_p4, bH1_gen_p4)) ||
	    (truthMatch(jet1_p4, bbH2_gen_p4) && truthMatch(jet2_p4, bbH1_gen_p4))) returnString="bH1bH2";
      
  if ((truthMatch(jet1_p4, bH1_gen_p4) && truthMatch(jet2_p4, bbH1_gen_p4)) || 
	    (truthMatch(jet1_p4, bbH1_gen_p4) && truthMatch(jet2_p4, bH1_gen_p4))) returnString="bH1bH1";
      
  if ((truthMatch(jet1_p4, bH2_gen_p4) && truthMatch(jet2_p4, bbH2_gen_p4)) ||
	    (truthMatch(jet1_p4, bbH2_gen_p4) && truthMatch(jet2_p4, bH2_gen_p4))) returnString="bH2bH2";
  
  return returnString;
}

int withinRegion(double mH1, double mH2, double r1=15., double r2=30., double mH1_c=H_mass, double mH2_c=H_mass)
{
  double r=pow(pow(mH1-mH1_c, 2)+pow(mH2-mH2_c, 2), 0.5);
  double angle=atan2(mH2-mH2_c, mH1-mH1_c);
  // std::cout<<"(mH1, mH2) = ("<<mH1<<", "<<mH2<<") lies in region ";
  int ret=-1;
  if (r<r1) ret=0;
  else if (r>r1 && r<r2)
  {
    if (angle>=0 && angle<pi/2.) ret=1;
    else if (angle>=pi/2. && angle<pi) ret=4;
    else if (angle<0 && angle>=-pi/2.) ret=2;
    else if (angle<pi/2.&& angle>=-pi) ret=3;
    else std::cout<<"This is within annulus but not within any CR!"<<std::endl;
  }
  else ret=5;
  // std::cout<<ret<<std::endl;
  return ret;
}
 

int HbbHbb_Partitioning(std::string sample, double h_mass, int kinConstraint=0)
{
  H_mass=h_mass;
  std::cout<<"H mass set at "<<H_mass<<std::endl;
	
	gSystem->Load("libPhysicsToolsKinFitter.so");
  
  std::string inputfilename="../OfficialStep2_KinematicallySelected_bTagged_"+sample+".root";
  // std::string inputfilename="../"+sample+"_selected_bTagged_.root";
  TChain *tree=new TChain("tree");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;
  
  // Book variables
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
  float eventWeight;
  GenParticleInfo genX, genH1, genH2, genH1B, genH1Bbar, genH2B, genH2Bbar;
  int H1jet1_i, H1jet2_i;
  int H2jet1_i, H2jet2_i;
  float regJet1E, regJet2E, regJet3E, regJet4E;
  float regJet1pT, regJet2pT, regJet3pT, regJet4pT;
  
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
  tree->SetBranchAddress("naJets", &(naJets));
  tree->SetBranchAddress("aJet_e", &(aJetE));
  tree->SetBranchAddress("aJet_pt", &(aJetpT));
  tree->SetBranchAddress("aJet_eta", &(aJeteta));
  tree->SetBranchAddress("aJet_phi", &(aJetphi));
  tree->SetBranchAddress("aJet_csv", &(aJetCSV));
  tree->SetBranchAddress("MET", &(metObj));
  tree->SetBranchAddress("genX", &(genX));
  tree->SetBranchAddress("genH1", &(genH1));
  tree->SetBranchAddress("genH2", &(genH2));
  tree->SetBranchAddress("genH1B", &(genH1B));
  tree->SetBranchAddress("genH1Bbar", &(genH1Bbar));
  tree->SetBranchAddress("genH2B", &(genH2B));
  tree->SetBranchAddress("genH2Bbar", &(genH2Bbar));
	tree->SetBranchAddress("H1jet1_i", &H1jet1_i);
  tree->SetBranchAddress("H1jet2_i", &H1jet2_i);
  tree->SetBranchAddress("H2jet1_i", &H2jet1_i);
  tree->SetBranchAddress("H2jet2_i", &H2jet2_i);
  tree->SetBranchAddress("regJet1E", &regJet1E);
  tree->SetBranchAddress("regJet2E", &regJet2E);
  tree->SetBranchAddress("regJet3E", &regJet3E);
  tree->SetBranchAddress("regJet4E", &regJet4E); 
  tree->SetBranchAddress("regJet1pT", &regJet1pT);
  tree->SetBranchAddress("regJet2pT", &regJet2pT);
  tree->SetBranchAddress("regJet3pT", &regJet3pT);
  tree->SetBranchAddress("regJet4pT", &regJet4pT); 
	tree->SetBranchAddress("eventWeight", &(eventWeight)); 
  
  TH1F *h_mX_CR1=new TH1F("h_mX_CR1", "h_mX_CR1", 200, 0., 2000.); h_mX_CR1->Sumw2();
  TH1F *h_mX_CR2=new TH1F("h_mX_CR2", "h_mX_CR2", 200, 0., 2000.); h_mX_CR2->Sumw2();
  TH1F *h_mX_CR3=new TH1F("h_mX_CR3", "h_mX_CR3", 200, 0., 2000.); h_mX_CR3->Sumw2();
  TH1F *h_mX_CR4=new TH1F("h_mX_CR4", "h_mX_CR4", 200, 0., 2000.); h_mX_CR4->Sumw2();
  TH1F *h_mX_CR5=new TH1F("h_mX_CR5", "h_mX_CR5", 200, 0., 2000.); h_mX_CR5->Sumw2();
  TH1F *h_mX_SR=new TH1F("h_mX_SR", "h_mX_SR", 200, 0., 2000.); h_mX_SR->Sumw2();
  TH1F *h_mX_SR_rightComb=new TH1F("h_mX_SR_rightComb", "h_mX_SR_rightComb", 200, 0., 2000.);
  TH1F *h_mX_SR_wrongComb=new TH1F("h_mX_SR_wrongComb", "h_mX_SR_wrongComb", 200, 0., 2000.);
  TH1F *h_mX_SR_noComb=new TH1F("h_mX_SR_noComb", "h_mX_SR_noComb", 200, 0., 2000.);
	TH1F *h_chi2=new TH1F("h_chi2", "Kinematic Square chi2", 100, 0., 100.);
	
	TH1F *h_jet1ET_noConst=new TH1F("h_jet1ET_noConst", "h_jet1ET_noConst", 50, 0., 500.);
	TH1F *h_jet2ET_noConst=new TH1F("h_jet2ET_noConst", "h_jet2ET_noConst", 50, 0., 250.);
	TH1F *h_jet1ET_kinConst=new TH1F("h_jet1ET_kinConst", "h_jet1ET_kinConst", 50, 0., 500.);
  TH1F *h_jet2ET_kinConst=new TH1F("h_jet2ET_kinConst", "h_jet2ET_kinConst", 50, 0., 250.);
	
	TH1F *h_jjTheta_noConst=new TH1F("h_jjTheta_noConst", "h_jjTheta_noConst", 50, 0., pi);
	TH1F *h_jjTheta_kinConst=new TH1F("h_jjTheta_kinConst", "h_jjTheta_kinConst", 50, 0., pi);
	
	TH1F *h_mH_noConst=new TH1F("h_mH_noConst", "Kinematically Unconstrained mH", 50, 50., 200.);
	TH1F *h_mH_kinConst=new TH1F("h_mH_kinConst", "Kinematically Constrained mH", 50, 50., 200.);
	
	TH1F *h_mX_SR_chi2=new TH1F("h_mX_SR_chi2", "h_mX_SR_chi2", 100, 0., 2000.); h_mX_SR_chi2->Sumw2();
  
  std::string histfilename="Histograms_"+sample+".root";
  gSystem->Exec(("cp ../"+histfilename+" "+histfilename).c_str());
  TFile *tFile1=new TFile(("../"+histfilename).c_str(), "READ");
  TH1F h_Cuts=*((TH1F*)((TH1F*)tFile1->Get("h_Cuts"))->Clone("h_Cuts"));
  tFile1->Close();
  
  // Loop over events
  int nEvents=tree->GetEntries();
  std::cout<<"Number of events in input file = "<<nEvents<<std::endl;
  double nCut7=0;
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);
    
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
    
    /*
    std::cout<<jetE[H1jet1_i]-regJet1E<<std::endl;
    std::cout<<jetE[H1jet2_i]-regJet2E<<std::endl;
    std::cout<<jetE[H2jet1_i]-regJet3E<<std::endl;
    std::cout<<jetE[H2jet2_i]-regJet4E<<std::endl;
    
    std::cout<<jetpT[H1jet1_i]-regJet1pT<<std::endl;
    std::cout<<jetpT[H1jet2_i]-regJet2pT<<std::endl;
    std::cout<<jetpT[H2jet1_i]-regJet3pT<<std::endl;
    std::cout<<jetpT[H2jet2_i]-regJet4pT<<std::endl;
    */
    
    jetE[H1jet1_i]=regJet1E;
    jetE[H1jet2_i]=regJet2E;
    jetE[H2jet1_i]=regJet3E;
    jetE[H2jet2_i]=regJet4E;
    jetpT[H1jet1_i]=regJet1pT;
    jetpT[H1jet2_i]=regJet2pT;
    jetpT[H2jet1_i]=regJet3pT;
    jetpT[H2jet2_i]=regJet4pT;
          
    TLorentzVector jet1_p4=fillTLorentzVector(jetpT[H1jet1_i], jeteta[H1jet1_i], jetphi[H1jet1_i], jetE[H1jet1_i]);      
    TLorentzVector jet2_p4=fillTLorentzVector(jetpT[H1jet2_i], jeteta[H1jet2_i], jetphi[H1jet2_i], jetE[H1jet2_i]);      
    TLorentzVector jet3_p4=fillTLorentzVector(jetpT[H2jet1_i], jeteta[H2jet1_i], jetphi[H2jet1_i], jetE[H2jet1_i]);      
    TLorentzVector jet4_p4=fillTLorentzVector(jetpT[H2jet2_i], jeteta[H2jet2_i], jetphi[H2jet2_i], jetE[H2jet2_i]);      
    
    TLorentzVector H1_p4=jet1_p4+jet2_p4;              
    TLorentzVector H2_p4=jet3_p4+jet4_p4;
    TLorentzVector X_p4=H1_p4+H2_p4;
          
    double H1_mass=H1_p4.M();
    double H2_mass=H2_p4.M();
    double X_mass=X_p4.M();
		
		// Put in kinematic constraining here
		h_jet1ET_noConst->Fill(jet1_p4.Et(), eventWeight);
		h_jet2ET_noConst->Fill(jet2_p4.Et(), eventWeight);
		h_mH_noConst->Fill((jet1_p4+jet2_p4).M(), eventWeight);
		h_jjTheta_noConst->Fill(jet1_p4.Angle(jet2_p4.Vect()), eventWeight);
		TLorentzVector jet1_new_p4, jet2_new_p4;
		if (kinConstraint==1)
		{
		  double chi2;
			TLorentzVector X_chi2_p4=H4b::Xchi2(jet1_p4, jet2_p4, jet3_p4, jet4_p4, H_mass, jet1_new_p4, jet2_new_p4);
			X_mass=X_chi2_p4.M();
			h_chi2->Fill(chi2, eventWeight);
			h_jet1ET_kinConst->Fill(jet1_new_p4.Et(), eventWeight);
			h_jet2ET_kinConst->Fill(jet2_new_p4.Et(), eventWeight);
			h_mH_kinConst->Fill((jet1_new_p4+jet2_new_p4).M(), eventWeight);
			h_jjTheta_kinConst->Fill(jet1_new_p4.Angle(jet2_new_p4.Vect()), eventWeight);
		}
    
    int region=withinRegion(H1_mass, H2_mass, 15., 35., H_mass, H_mass);
    if (region==0) // SR                                                  
    {                                                                     
      ++nCut7;                                                            
      h_mX_SR->Fill(X_mass, eventWeight);
      // Figure out combinatorics that went in
      TLorentzVector bH1_gen_p4=fillTLorentzVector(genH1B);
      TLorentzVector bbH1_gen_p4=fillTLorentzVector(genH1Bbar);
      TLorentzVector bH2_gen_p4=fillTLorentzVector(genH2B);
      TLorentzVector bbH2_gen_p4=fillTLorentzVector(genH2Bbar);
      std::string string_H1=truthMatching(jet1_p4, jet2_p4, bH1_gen_p4, bbH1_gen_p4, bH2_gen_p4, bbH2_gen_p4);
      std::string string_H2=truthMatching(jet3_p4, jet4_p4, bH1_gen_p4, bbH1_gen_p4, bH2_gen_p4, bbH2_gen_p4);
      if ((string_H1=="bH1bH1" && string_H2=="bH2bH2") || (string_H1=="bH2bH2" && string_H2=="bH1bH1")) h_mX_SR_rightComb->Fill(X_mass, eventWeight);
      else if (string_H1=="bH1bH2" || string_H2=="bH1bH2") h_mX_SR_wrongComb->Fill(X_mass, eventWeight);
      else if (string_H1=="noMatch" || string_H2=="noMatch") h_mX_SR_noComb->Fill(X_mass, eventWeight);
    }                                                                     
    else if (region==1) // CR1                                            
    {                                                                     
      h_mX_CR1->Fill(X_mass, eventWeight);                          
    }                                                                     
    else if (region==2) // CR2                                            
    {                                                                     
      h_mX_CR2->Fill(X_mass, eventWeight);                          
    }                                                                     
    else if (region==3) // CR3                                            
    {                                                                     
      h_mX_CR3->Fill(X_mass, eventWeight);                          
    }                                                                     
    else if (region==4) // CR4                                            
    {                                                                     
      h_mX_CR4->Fill(X_mass, eventWeight);                           
    }                                                                     
    else if (region==5) // CR5                                            
    {                                                                     
      h_mX_CR5->Fill(X_mass, eventWeight);                           
    }                                                                     
    else if (region==-1)                                                  
    {                                                                     
      std::cout<<"Didn't fall in any region!"<<std::endl;                 
    }
    
  } // Event loop
  
  std::cout<<"out of event loop"<<std::endl;
  std::cout<<"nCut7 = "<<nCut7<<std::endl;
  h_Cuts.Fill(15, nCut7);
  
  TFile *tFile2=new TFile(histfilename.c_str(), "UPDATE");
  tFile2->Delete("h_Cuts;1");
  h_mX_CR1->Write();
  h_mX_CR2->Write();
  h_mX_CR3->Write();
  h_mX_CR4->Write();
  h_mX_CR5->Write();
  h_mX_SR->Write();
  h_mX_SR_rightComb->Write();
  h_mX_SR_wrongComb->Write();
  h_mX_SR_noComb->Write();
	h_mX_SR_chi2->Write();
	h_chi2->Write();
	h_jet1ET_noConst->Write();
	h_jet2ET_noConst->Write();
	h_jet1ET_kinConst->Write();
  h_jet2ET_kinConst->Write();
	h_jjTheta_noConst->Write();
  h_jjTheta_kinConst->Write();
	h_mH_noConst->Write();
	h_mH_kinConst->Write();
  h_Cuts.Write();
  tFile2->Write();
  tFile2->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  
  std::cout<<"=== Cut Efficiencies === "<<std::endl;
  std::cout<<"Number of events in SR = "<<nCut7<<std::endl;
  std::cout<<"========================"<<std::endl;
  
  return 0;
}
            
      
  
  
  
  
