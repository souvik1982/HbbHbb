int CSVFlavorEfficiencies_Display()
{
  TFile *file=new TFile("CSVFlavorEfficiencies_M.root");
  TH1F *h_udsg_JetpT_den=(TH1F*)gDirectory->Get("h_udsg_JetpT_den");
  TH1F *h_udsg_Jeteta_den=(TH1F*)gDirectory->Get("h_udsg_Jeteta_den");
  TH2F *h_udsg_JetpTeta_den=(TH2F*)gDirectory->Get("h_udsg_JetpTeta_den");
  TH1F *h_udsg_JetpT_num=(TH1F*)gDirectory->Get("h_udsg_JetpT_num");
  TH1F *h_udsg_Jeteta_num=(TH1F*)gDirectory->Get("h_udsg_Jeteta_num");
  TH2F *h_udsg_JetpTeta_num=(TH2F*)gDirectory->Get("h_udsg_JetpTeta_num");
  TH1F *h_c_JetpT_den=(TH1F*)gDirectory->Get("h_c_JetpT_den");
  TH1F *h_c_Jeteta_den=(TH1F*)gDirectory->Get("h_c_Jeteta_den");
  TH2F *h_c_JetpTeta_den=(TH2F*)gDirectory->Get("h_c_JetpTeta_den");
  TH1F *h_c_JetpT_num=(TH1F*)gDirectory->Get("h_c_JetpT_num");
  TH1F *h_c_Jeteta_num=(TH1F*)gDirectory->Get("h_c_Jeteta_num");
  TH2F *h_c_JetpTeta_num=(TH2F*)gDirectory->Get("h_c_JetpTeta_num");
  TH1F *h_b_JetpT_den=(TH1F*)gDirectory->Get("h_b_JetpT_den");
  TH1F *h_b_Jeteta_den=(TH1F*)gDirectory->Get("h_b_Jeteta_den");
  TH2F *h_b_JetpTeta_den=(TH2F*)gDirectory->Get("h_b_JetpTeta_den");
  TH1F *h_b_JetpT_num=(TH1F*)gDirectory->Get("h_b_JetpT_num");
  TH1F *h_b_Jeteta_num=(TH1F*)gDirectory->Get("h_b_Jeteta_num");
  TH2F *h_b_JetpTeta_num=(TH2F*)gDirectory->Get("h_b_JetpTeta_num");
  
  h_udsg_JetpT_den->Rebin(2);
  h_udsg_JetpT_num->Rebin(2);
  h_udsg_JetpTeta_den->Rebin2D(2,2);
  h_udsg_Jeteta_den->Rebin(2);
  h_udsg_Jeteta_num->Rebin(2);
  h_udsg_JetpTeta_num->Rebin2D(2,2);
  
  h_c_JetpT_den->Rebin(3);
  h_c_JetpT_num->Rebin(3);
  h_c_JetpTeta_den->Rebin2D(3,3);
  h_c_Jeteta_den->Rebin(3);
  h_c_Jeteta_num->Rebin(3);
  h_c_JetpTeta_num->Rebin2D(3,3);
  
  h_b_JetpT_den->Rebin(2);
  h_b_JetpT_num->Rebin(2);
  h_b_JetpTeta_den->Rebin2D(2,3);
  h_b_Jeteta_den->Rebin(3);
  h_b_Jeteta_num->Rebin(3);
  h_b_JetpTeta_num->Rebin2D(2,3);
  
  
  TGraphAsymmErrors *g_udsg_JetpT=new TGraphAsymmErrors(h_udsg_JetpT_num, h_udsg_JetpT_den, "cp"); g_udsg_JetpT->SetTitle("b-tagging Efficiency for udsg Jets; p_{T}");
  TGraphAsymmErrors *g_udsg_Jeteta=new TGraphAsymmErrors(h_udsg_Jeteta_num, h_udsg_Jeteta_den, "cp"); g_udsg_Jeteta->SetTitle("b-tagging Efficiency for udsg Jets; #eta");
  
  TGraphAsymmErrors *g_c_JetpT=new TGraphAsymmErrors(h_c_JetpT_num, h_c_JetpT_den, "cp"); g_c_JetpT->SetTitle("b-tagging Efficiency for c Jets; p_{T}");
  TGraphAsymmErrors *g_c_Jeteta=new TGraphAsymmErrors(h_c_Jeteta_num, h_c_Jeteta_den, "cp"); g_c_Jeteta->SetTitle("b-tagging Efficiency for c Jets; #eta");
  
  TGraphAsymmErrors *g_b_JetpT=new TGraphAsymmErrors(h_b_JetpT_num, h_b_JetpT_den, "cp"); g_b_JetpT->SetTitle("b-tagging Efficiency for b Jets; p_{T}");
  TGraphAsymmErrors *g_b_Jeteta=new TGraphAsymmErrors(h_b_Jeteta_num, h_b_Jeteta_den, "cp"); g_b_Jeteta->SetTitle("b-tagging Efficiency for b Jets; #eta");
  
  TH2F *h_udsg_JetpTeta=h_udsg_JetpTeta_num->Clone("h_udsg_JetpTeta"); h_udsg_JetpTeta->SetTitle("b-tagging Efficiency for udsg Jets; p_{T} (GeV/c); #eta");
  h_udsg_JetpTeta->Divide(h_udsg_JetpTeta_den);
  
  TH2F *h_c_JetpTeta=h_c_JetpTeta_num->Clone("h_c_JetpTeta"); h_c_JetpTeta->SetTitle("b-tagging Efficiency for c Jets; p_{T} (GeV/c); #eta");
  h_c_JetpTeta->Divide(h_c_JetpTeta_den);
  
  TH2F *h_b_JetpTeta=h_b_JetpTeta_num->Clone("h_b_JetpTeta"); h_b_JetpTeta->SetTitle("b-tagging Efficiency for b Jets; p_{T} (GeV/c); #eta");
  h_b_JetpTeta->Divide(h_b_JetpTeta_den);
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(000000000);
  gStyle->SetPalette(1);
  
  TCanvas *c_udsg_Efficiency=new TCanvas("c_udsg_Efficiency", "", 1200, 400);
  c_udsg_Efficiency->Divide(3,1);
  c_udsg_Efficiency->cd(1);
  g_udsg_JetpT->SetMinimum(0.);
  g_udsg_JetpT->Draw("AP");
  c_udsg_Efficiency->cd(2);
  g_udsg_JetpT->SetMinimum(0.);
  g_udsg_Jeteta->GetYaxis()->SetRangeUser(0., 0.03);
  g_udsg_Jeteta->Draw("AP");
  c_udsg_Efficiency->cd(3);
  h_udsg_JetpTeta->SetMinimum(0.);
  h_udsg_JetpTeta->Draw("LEGO2");
  c_udsg_Efficiency->SaveAs("c_udsg_Efficiency.png");
  
  TCanvas *c_c_Efficiency=new TCanvas("c_c_Efficiency", "", 1200, 400);
  c_c_Efficiency->Divide(3,1);
  c_c_Efficiency->cd(1);
  g_c_JetpT->SetMinimum(0.);
  g_c_JetpT->Draw("AP");
  c_c_Efficiency->cd(2);
  g_c_Jeteta->SetMinimum(0.);
  g_c_Jeteta->GetYaxis()->SetRangeUser(0., 0.3);
  g_c_Jeteta->Draw("AP");
  c_c_Efficiency->cd(3);
  h_c_JetpTeta->SetMinimum(0.);
  h_c_JetpTeta->Draw("LEGO2");
  c_c_Efficiency->SaveAs("c_c_Efficiency.png");
  
  TCanvas *c_b_Efficiency=new TCanvas("c_b_Efficiency", "", 1200, 400);
  c_b_Efficiency->Divide(3,1);
  c_b_Efficiency->cd(1);
  g_b_JetpT->SetMinimum(0);
  g_b_JetpT->Draw("AP");
  c_b_Efficiency->cd(2);
  g_b_Jeteta->SetMinimum(0);
  g_b_Jeteta->Draw("AP");
  c_b_Efficiency->cd(3);
  h_b_JetpTeta->SetMinimum(0);
  h_b_JetpTeta->Draw("LEGO2");
  c_b_Efficiency->SaveAs("c_b_Efficiency.png");
}
  
