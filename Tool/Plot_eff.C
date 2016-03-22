#include "CMS_lumi.h"

void Plot_eff()
{
  TFile f_eff("v160302_ttbarSTWv6_Effs2dPlots.root");

  const std::string titre="CMS Supplementary";

  const bool do_c1=true;
  const bool do_c2=true;
  const bool do_c3=true;
  const bool do_c4=true;

  const int nxbin=7;
  Float_t xbins[nxbin+1]={10.0,20.0,30.0,40.0,50.0,70.0,100.0,120.0};
  //Float_t ybins[9]={0.0,5.0,10.0,20.0,40.0,60.0,80.0,100.0,120.0};
  const int nybin=5;
  Float_t ybins[nybin+1]={0.005,0.02,0.05,0.15,1.0,10.0};
 

  if (do_c1)
  {

  TCanvas *c1 = new TCanvas("c1", "c1",0,51,1920,1004);
  c1->SetFillColor(0);
  c1->cd();
  c1->SetLogy();

  TH2F *h_mu_iso_eff_ori = (TH2F*)f_eff.Get("mus_isoeffs");
  //h_mu_iso_eff_ori->Draw("colztext");

  TH2F *h_mu_iso_eff;
 h_mu_iso_eff = new TH2F("h_mu_iso_eff", "h_mu_iso_eff", nxbin, xbins, nybin, ybins);
  h_mu_iso_eff->SetTitle("");
  h_mu_iso_eff->SetXTitle("muon p_{T} [GeV]");
  h_mu_iso_eff->SetYTitle("activity");
  h_mu_iso_eff->SetStats(0);
  h_mu_iso_eff->SetMarkerSize(1.5);

  for (Int_t muonptc=1;muonptc<nxbin+1;++muonptc) {
    for (int activityc=1;activityc<nybin+1;++activityc)
    {
      h_mu_iso_eff->SetBinContent(muonptc,activityc,h_mu_iso_eff_ori->GetBinContent(muonptc,activityc+1));
      h_mu_iso_eff->SetBinError(muonptc,activityc,h_mu_iso_eff_ori->GetBinError(muonptc,activityc+1));
    }
  }

  gStyle->SetPaintTextFormat("1.2f");
  h_mu_iso_eff->Draw("colztexte");

//  TLatex *title = new TLatex(0.09770115,0.9194915,titre.c_str());
//  title->SetNDC();
//  title->SetTextSize(0.045);
//  title->Draw("same");
//
//   TLatex *   tex2 = new TLatex(0.8127615,0.9178645,"(13 TeV)");
//tex2->SetNDC();
//   tex2->SetTextSize(0.045);
//   tex2->SetLineWidth(2);
//   tex2->Draw();

   CMSStylePlot::CMS_lumi( c1, 0, 0 );

  c1->SaveAs( "_mu_2d_iso_eff.png" );
  c1->SaveAs( "2deffs_mus_iso.pdf" );
  c1->SaveAs( "_mu_2d_iso_eff.C" );
  }

  ///////////////////
  //     c2
  ///////////////////

  if (do_c2)
  {


  TCanvas *c2 = new TCanvas("c2", "c2",0,51,1920,1004);
  c2->SetFillColor(0);
  c2->cd();
  c2->SetLogy();

  TH2F *h_els_iso_eff_ori = (TH2F*)f_eff.Get("els_isoeffs");
  //h_els_iso_eff_ori->Draw("colztext");

  TH2F *h_els_iso_eff;
  //Float_t xbins[9]={5.0,10.0,20.0,30.0,40.0,50.0,70.0,100.0,120.0};
  //Float_t ybins[9]={0.0,5.0,10.0,20.0,40.0,60.0,80.0,100.0,120.0};
  h_els_iso_eff = new TH2F("h_els_iso_eff", "h_els_iso_eff", 7, xbins, nybin, ybins);
  h_els_iso_eff->SetTitle("");
  h_els_iso_eff->SetXTitle("electron p_{T} [GeV]");
  h_els_iso_eff->SetYTitle("activity");
  h_els_iso_eff->SetStats(0);
  h_els_iso_eff->SetMarkerSize(1.5);

  for (Int_t elsonptc=1;elsonptc<nxbin+1;++elsonptc) {
    for (int activityc=1;activityc<nybin+1;++activityc)
    {
      h_els_iso_eff->SetBinContent(elsonptc,activityc,h_els_iso_eff_ori->GetBinContent(elsonptc,activityc+1));
      h_els_iso_eff->SetBinError(elsonptc,activityc,h_els_iso_eff_ori->GetBinError(elsonptc,activityc+1));
    }
  }

  gStyle->SetPaintTextFormat("1.2f");
  h_els_iso_eff->Draw("colztexte");

//  title->SetNDC();
//  title->SetTextSize(0.045);
//  title->Draw("same");
//
//   TLatex *   tex2_c2 = new TLatex(0.8127615,0.9178645,"(13 TeV)");
//tex2_c2->SetNDC();
//   tex2_c2->SetTextSize(0.045);
//   tex2_c2->SetLineWidth(2);
//   tex2_c2->Draw();

   CMSStylePlot::CMS_lumi( c2, 0, 0 );

  c2->SaveAs( "_el_2d_iso_eff.png" );
  c2->SaveAs( "2deffs_els_iso.pdf" );
  c2->SaveAs( "_el_2d_iso_eff.C" );

  }

  ///////////////////
  //     c3
  ///////////////////
  if (do_c3)
  {

  TCanvas *c3 = new TCanvas("c3", "c3",0,51,1920,1004);
  c3->SetFillColor(0);
  c3->cd();
  c3->SetLogy();

  TH2F *h_mus_reco_eff_ori = (TH2F*)f_eff.Get("mus_recoeffs");
  //h_mus_reco_eff_ori->Draw("colztext");

  TH2F *h_mus_reco_eff;
  //Float_t xbins[9]={5.0,10.0,20.0,30.0,40.0,50.0,70.0,100.0,120.0};
  //Float_t ybins[9]={0.0,5.0,10.0,20.0,40.0,60.0,80.0,100.0,120.0};
  h_mus_reco_eff = new TH2F("h_mus_reco_eff", "h_mus_reco_eff", 7, xbins, nybin, ybins);
  h_mus_reco_eff->SetTitle("");
  h_mus_reco_eff->SetXTitle("muon p_{T} [GeV]");
  h_mus_reco_eff->SetYTitle("activity");
  h_mus_reco_eff->SetStats(0);
  h_mus_reco_eff->SetMarkerSize(1.5);

  for (Int_t musonptc=1;musonptc<nxbin+1;++musonptc) {
    for (int activityc=1;activityc<nybin+1;++activityc)
    {
      h_mus_reco_eff->SetBinContent(musonptc,activityc,h_mus_reco_eff_ori->GetBinContent(musonptc,activityc+1));
      h_mus_reco_eff->SetBinError(musonptc,activityc,h_mus_reco_eff_ori->GetBinError(musonptc,activityc+1));
    }
  }
  h_mus_reco_eff->SetMinimum(0.5);

  gStyle->SetPaintTextFormat("1.2f");
  h_mus_reco_eff->Draw("colztexte");

//  title->SetNDC();
//  title->SetTextSize(0.045);
//  title->Draw("same");
//
//   TLatex *   tex2_c3 = new TLatex(0.8127615,0.9178645,"(13 TeV)");
//tex2_c3->SetNDC();
//   tex2_c3->SetTextSize(0.045);
//   tex2_c3->SetLineWidth(2);
//   tex2_c3->Draw();

   CMSStylePlot::CMS_lumi( c3, 0, 0 );

  c3->SaveAs( "_mu_2d_reco_eff.png" );
  c3->SaveAs( "2deffs_mus_reco.pdf" );
  c3->SaveAs( "_mu_2d_reco_eff.C" );
  }
  ///////////////////
  //     c4
  ///////////////////
  if (do_c4)
  {

  TCanvas *c4 = new TCanvas("c4", "c4",0,51,1920,1004);
  c4->SetFillColor(0);
  c4->cd();
  c4->SetLogy();

  TH2F *h_els_reco_eff_ori = (TH2F*)f_eff.Get("els_recoeffs");
  //h_els_reco_eff_ori->Draw("colztext");

  TH2F *h_els_reco_eff;
  //Float_t xbins[9]={5.0,10.0,20.0,30.0,40.0,50.0,70.0,100.0,120.0};
  //Float_t ybins[9]={0.0,5.0,10.0,20.0,40.0,60.0,80.0,100.0,120.0};
  h_els_reco_eff = new TH2F("h_els_reco_eff", "h_els_reco_eff", 7, xbins, nybin, ybins);
  h_els_reco_eff->SetTitle("");
  h_els_reco_eff->SetXTitle("electron p_{T} [GeV]");
  h_els_reco_eff->SetYTitle("activity");
  h_els_reco_eff->SetStats(0);
  h_els_reco_eff->SetMarkerSize(1.5);

  for (Int_t elsonptc=1;elsonptc<nxbin+1;++elsonptc) {
    for (int activityc=1;activityc<nybin+1;++activityc)
    {
      h_els_reco_eff->SetBinContent(elsonptc,activityc,h_els_reco_eff_ori->GetBinContent(elsonptc,activityc+1));
      h_els_reco_eff->SetBinError(elsonptc,activityc,h_els_reco_eff_ori->GetBinError(elsonptc,activityc+1));
    }
  }

  gStyle->SetPaintTextFormat("1.2f");
  h_els_reco_eff->Draw("colztexte");

//  title->SetNDC();
//  title->SetTextSize(0.045);
//  title->Draw("same");
//
//   TLatex *   tex2_c4 = new TLatex(0.8127615,0.9178645,"(13 TeV)");
//tex2_c4->SetNDC();
//   tex2_c4->SetTextSize(0.045);
//   tex2_c4->SetLineWidth(2);
//   tex2_c4->Draw();

   CMSStylePlot::CMS_lumi( c4, 0, 0 );

  c4->SaveAs( "_el_2d_reco_eff.png" );
  c4->SaveAs( "2deffs_els_reco.pdf" );
  c4->SaveAs( "_el_2d_reco_eff.C" );
  }
}
