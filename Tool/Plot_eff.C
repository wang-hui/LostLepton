{
  TFile f_eff("Effs2dPlots.root");

  TCanvas *c1 = new TCanvas("c1", "c1",0,51,1920,1004);
  c1->SetFillColor(0);
  c1->cd();

  TH2F *h_mu_iso_eff_ori = (TH2F*)f_eff.Get("mus_isoeffs");
  //h_mu_iso_eff_ori->Draw("colztext");

  TH2F *h_mu_iso_eff;
  Float_t xbins[9]={5.0,10.0,20.0,30.0,40.0,50.0,70.0,100.0,120.0};
  Float_t ybins[9]={0.0,5.0,10.0,20.0,40.0,60.0,80.0,100.0,120.0};
  h_mu_iso_eff = new TH2F("h_mu_iso_eff", "h_mu_iso_eff", 8, xbins, 8, ybins);
  h_mu_iso_eff->SetTitle("");
  h_mu_iso_eff->SetXTitle("muon pT [GeV]");
  h_mu_iso_eff->SetYTitle("activity [GeV]");
  h_mu_iso_eff->SetStats(0);

  for (Int_t muonptc=1;muonptc<9;++muonptc) {
    for (int activityc=1;activityc<9;++activityc)
    {
      h_mu_iso_eff->SetBinContent(muonptc,activityc,h_mu_iso_eff_ori->GetBinContent(muonptc,activityc));
      h_mu_iso_eff->SetBinError(muonptc,activityc,h_mu_iso_eff_ori->GetBinError(muonptc,activityc));
    }
  }
  gStyle->SetPaintTextFormat("1.2f");
  h_mu_iso_eff->Draw("colztexte");


  ///////////////////
  //     c2
  ///////////////////

  TCanvas *c2 = new TCanvas("c2", "c2",0,51,1920,1004);
  c2->SetFillColor(0);
  c2->cd();

  TH2F *h_els_iso_eff_ori = (TH2F*)f_eff.Get("els_isoeffs");
  //h_els_iso_eff_ori->Draw("colztext");

  TH2F *h_els_iso_eff;
  //Float_t xbins[9]={5.0,10.0,20.0,30.0,40.0,50.0,70.0,100.0,120.0};
  //Float_t ybins[9]={0.0,5.0,10.0,20.0,40.0,60.0,80.0,100.0,120.0};
  h_els_iso_eff = new TH2F("h_els_iso_eff", "h_els_iso_eff", 8, xbins, 8, ybins);
  h_els_iso_eff->SetTitle("");
  h_els_iso_eff->SetXTitle("electron pT [GeV]");
  h_els_iso_eff->SetYTitle("activity [GeV]");
  h_els_iso_eff->SetStats(0);

  for (Int_t elsonptc=1;elsonptc<9;++elsonptc) {
    for (int activityc=1;activityc<9;++activityc)
    {
      h_els_iso_eff->SetBinContent(elsonptc,activityc,h_els_iso_eff_ori->GetBinContent(elsonptc,activityc));
      h_els_iso_eff->SetBinError(elsonptc,activityc,h_els_iso_eff_ori->GetBinError(elsonptc,activityc));
    }
  }
  gStyle->SetPaintTextFormat("1.2f");
  h_els_iso_eff->Draw("colztexte");

  ///////////////////
  //     c3
  ///////////////////

  TCanvas *c3 = new TCanvas("c3", "c3",0,51,1920,1004);
  c3->SetFillColor(0);
  c3->cd();

  TH2F *h_mus_reco_eff_ori = (TH2F*)f_eff.Get("mus_recoeffs");
  //h_mus_reco_eff_ori->Draw("colztext");

  TH2F *h_mus_reco_eff;
  //Float_t xbins[9]={5.0,10.0,20.0,30.0,40.0,50.0,70.0,100.0,120.0};
  //Float_t ybins[9]={0.0,5.0,10.0,20.0,40.0,60.0,80.0,100.0,120.0};
  h_mus_reco_eff = new TH2F("h_mus_reco_eff", "h_mus_reco_eff", 8, xbins, 8, ybins);
  h_mus_reco_eff->SetTitle("");
  h_mus_reco_eff->SetXTitle("muon pT [GeV]");
  h_mus_reco_eff->SetYTitle("activity [GeV]");
  h_mus_reco_eff->SetStats(0);

  for (Int_t musonptc=1;musonptc<9;++musonptc) {
    for (int activityc=1;activityc<9;++activityc)
    {
      h_mus_reco_eff->SetBinContent(musonptc,activityc,h_mus_reco_eff_ori->GetBinContent(musonptc,activityc));
      h_mus_reco_eff->SetBinError(musonptc,activityc,h_mus_reco_eff_ori->GetBinError(musonptc,activityc));
    }
  }
  gStyle->SetPaintTextFormat("1.2f");
  h_mus_reco_eff->Draw("colztexte");

  ///////////////////
  //     c4
  ///////////////////

  TCanvas *c4 = new TCanvas("c4", "c4",0,51,1920,1004);
  c4->SetFillColor(0);
  c4->cd();

  TH2F *h_els_reco_eff_ori = (TH2F*)f_eff.Get("els_recoeffs");
  //h_els_reco_eff_ori->Draw("colztext");

  TH2F *h_els_reco_eff;
  //Float_t xbins[9]={5.0,10.0,20.0,30.0,40.0,50.0,70.0,100.0,120.0};
  //Float_t ybins[9]={0.0,5.0,10.0,20.0,40.0,60.0,80.0,100.0,120.0};
  h_els_reco_eff = new TH2F("h_els_reco_eff", "h_els_reco_eff", 8, xbins, 8, ybins);
  h_els_reco_eff->SetTitle("");
  h_els_reco_eff->SetXTitle("electron pT [GeV]");
  h_els_reco_eff->SetYTitle("activity [GeV]");
  h_els_reco_eff->SetStats(0);

  for (Int_t elsonptc=1;elsonptc<9;++elsonptc) {
    for (int activityc=1;activityc<9;++activityc)
    {
      h_els_reco_eff->SetBinContent(elsonptc,activityc,h_els_reco_eff_ori->GetBinContent(elsonptc,activityc));
      h_els_reco_eff->SetBinError(elsonptc,activityc,h_els_reco_eff_ori->GetBinError(elsonptc,activityc));
    }
  }
  gStyle->SetPaintTextFormat("1.2f");
  h_els_reco_eff->Draw("colztexte");


}
