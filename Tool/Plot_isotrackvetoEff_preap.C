{
  double isoTrackEff_SB[37] = {
0.614593, 0.617546, 0.634198, 0.587526, 0.619124, 0.638449, 0.595616, 0.56704, 0.612878, 0.64588, 0.729552, 0.609632, 0.628899, 0.605931, 0.626885, 0.605384, 0.593897, 0.628592, 0.617852, 0.600823, 0.789169, 0.565181, 0.571983, 0.567818, 0.613757, 0.603237, 0.54181, 0.668762, 0.596302, 0.581645, 0.5509, 0.653592, 0.608133, 0.603655, 0.599706, 0.514887, 0.69194, };

double isoTrackErr[37];
 isoTrackErr[0] = 0.0064399450138;
 isoTrackErr[1] = 0.011496246345;
 isoTrackErr[2] = 0.019285668138;
 isoTrackErr[3] = 0.037330456215;
 isoTrackErr[4] = 0.012420347007;
 isoTrackErr[5] = 0.015367213817;
 isoTrackErr[6] = 0.026068904857;
 isoTrackErr[7] = 0.043041328322;
 isoTrackErr[8] = 0.036010157058;
 isoTrackErr[9] = 0.025627383556;
 isoTrackErr[10] = 0.031484245268;
 isoTrackErr[11] = 0.005816864258;
 isoTrackErr[12] = 0.010140718489;
 isoTrackErr[13] = 0.017111155375;
 isoTrackErr[14] = 0.029890178181;
 isoTrackErr[15] = 0.01576539081;
 isoTrackErr[16] = 0.018700075717;
 isoTrackErr[17] = 0.023491214954;
 isoTrackErr[18] = 0.04294318543;
 isoTrackErr[19] = 0.033735096675;
 isoTrackErr[20] = 0.042189978275;
 isoTrackErr[21] = 0.010359565724;
 isoTrackErr[22] = 0.020856206271;
 isoTrackErr[23] = 0.034100135733;
 isoTrackErr[24] = 0.014479372271;
 isoTrackErr[25] = 0.022201963212;
 isoTrackErr[26] = 0.030778790437;
 isoTrackErr[27] = 0.03788246376;
 isoTrackErr[28] = 0.033413122588;
 isoTrackErr[29] = 0.010833724279;
 isoTrackErr[30] = 0.022349160002;
 isoTrackErr[31] = 0.036803132395;
 isoTrackErr[32] = 0.015652160221;
 isoTrackErr[33] = 0.023558335526;
 isoTrackErr[34] = 0.034077324449;
 isoTrackErr[35] = 0.045691322437;
 isoTrackErr[36] = 0.032878319254;

  const std::string titre="CMS Simulation 2016, #sqrt{s} = 13 TeV";

  TCanvas *c1 = new TCanvas("c1", "c1",0,51,1920,1004);
  c1->SetFillColor(0);
  c1->cd();

  TH1D *h_mu_iso_eff = new TH1D("h_mu_iso_eff", "h_mu_iso_eff", 37 , 0 , 37);
  h_mu_iso_eff->SetTitle("");
  h_mu_iso_eff->SetXTitle("Search Bins");
  h_mu_iso_eff->SetYTitle("iso track veto efficiency");
  h_mu_iso_eff->SetStats(0);
  h_mu_iso_eff->SetLineWidth(3);

  for (Int_t sbc=1;sbc<=37;++sbc) {
    h_mu_iso_eff->SetBinContent(sbc,isoTrackEff_SB[sbc-1]);
    h_mu_iso_eff->SetBinError(sbc,isoTrackErr[sbc-1]);
  }

  gStyle->SetPaintTextFormat("1.2f");
  h_mu_iso_eff->Draw();

  TLatex *title = new TLatex(0.09770115,0.9194915,titre.c_str());
  title->SetNDC();
  title->SetTextSize(0.045);
  title->Draw("same");

  c1->SaveAs( "h_isotrackvetoEff.png" );

}
