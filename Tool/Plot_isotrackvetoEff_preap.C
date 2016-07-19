#include "CMS_lumi.h"

void Plot_isotrackvetoEff_preap()
{
  double isoTrackEff_SB[37] = {0.61438, 0.596249, 0.621453, 0.62781, 0.608466, 0.610716, 0.662686, 0.661158, 0.628881, 0.667907, 0.675414, 0.608604, 0.636072, 0.620656, 0.664064, 0.628216, 0.648727, 0.607444, 0.597922, 0.689885, 0.798671, 0.548954, 0.539646, 0.575433, 0.599655, 0.60492, 0.577339, 0.641311, 0.523065, 0.587981, 0.583319, 0.679883, 0.582554, 0.605587, 0.606145, 0.532442, 0.681606};

double isoTrackErr[37];
isoTrackErr[0] = 0.0065161768148;
isoTrackErr[1] = 0.011970003298;
isoTrackErr[2] = 0.021341294997;
isoTrackErr[3] = 0.043182665314;
isoTrackErr[4] = 0.012646661427;
isoTrackErr[5] = 0.015483745016;
isoTrackErr[6] = 0.026132731854;
isoTrackErr[7] = 0.050315177651;
isoTrackErr[8] = 0.039410194819;
isoTrackErr[9] = 0.027308312153;
isoTrackErr[10] = 0.03333558055;
isoTrackErr[11] = 0.0060849980935;
isoTrackErr[12] = 0.010894531834;
isoTrackErr[13] = 0.018780241888;
isoTrackErr[14] = 0.035394571365;
isoTrackErr[15] = 0.016673859137;
isoTrackErr[16] = 0.019516455798;
isoTrackErr[17] = 0.026650595421;
isoTrackErr[18] = 0.047586995334;
isoTrackErr[19] = 0.03809529272;
isoTrackErr[20] = 0.048463259589;
isoTrackErr[21] = 0.010726196485;
isoTrackErr[22] = 0.022787744916;
isoTrackErr[23] = 0.040441213778;
isoTrackErr[24] = 0.015178894685;
isoTrackErr[25] = 0.023896511871;
isoTrackErr[26] = 0.036463666443;
isoTrackErr[27] = 0.044898492361;
isoTrackErr[28] = 0.03787541448;
isoTrackErr[29] = 0.011482565913;
isoTrackErr[30] = 0.024669260564;
isoTrackErr[31] = 0.041694741508;
isoTrackErr[32] = 0.01713434815;
isoTrackErr[33] = 0.026187526589;
isoTrackErr[34] = 0.040168741147;
isoTrackErr[35] = 0.051369476129;
isoTrackErr[36] = 0.040775186172;

  const std::string titre="CMS Supplementary";

  TCanvas *c1 = new TCanvas("c1", "c1",0,51,1920,1004);
  c1->SetFillColor(0);
  c1->cd();

  TH1D *h_mu_iso_eff = new TH1D("h_mu_iso_eff", "h_mu_iso_eff", 37 , 0 , 37);
  h_mu_iso_eff->SetTitle("");
  h_mu_iso_eff->SetXTitle("Search region bin number");
  h_mu_iso_eff->SetYTitle("iso track veto efficiency");
  h_mu_iso_eff->SetStats(0);
  h_mu_iso_eff->SetLineWidth(3);
  h_mu_iso_eff->SetMinimum(0.0);

  for (Int_t sbc=1;sbc<=37;++sbc)
  {
    h_mu_iso_eff->SetBinContent(sbc,isoTrackEff_SB[sbc-1]);
    h_mu_iso_eff->SetBinError(sbc,isoTrackErr[sbc-1]);
  }

  gStyle->SetPaintTextFormat("1.2f");
  h_mu_iso_eff->Draw();

//  TLatex *title = new TLatex(0.09770115,0.9194915,titre.c_str());
//  title->SetNDC();
//  title->SetTextSize(0.045);
//  title->Draw("same");
//  TLatex *   tex2 = new TLatex(0.8127615,0.9178645,"(13 TeV)");
//  tex2->SetNDC();
//  tex2->SetTextSize(0.045);
//  tex2->SetLineWidth(2);
//  tex2->Draw();

  CMSStylePlot::CMS_lumi( c1, 0, 0 );
  c1->SaveAs( "h_isotrackvetoEff.png" );
}
