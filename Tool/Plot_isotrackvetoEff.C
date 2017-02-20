#include "CMS_lumi.h"

void Plot_isotrackvetoEff()
{
double isoTrackEff_SB[84] = {0.580997, 0.58417, 0.620865, 0.603369, 0.442557, 0.586969, 0.587196, 0.610037, 0.716258, 0.572515, 0.639494, 0.500386, 0.563177, 0, 0.408404, 0.672944, 0.62005, 1, 0, 0.674076, 1, 0.570958, 0.589226, 0.666014, 0.577494, 0.5728, 0.563938, 0.630374, 0.642516, 1, 0.675154, 0.524279, 0.499238, 0.599968, 0.5, 1, 0, 0.585174, 0.577849, 0.582268, 0.557018, 0.60262, 0.559008, 0.667016, 0.635809, 0.629136, 0.5, 0.672397, 0.529846, 0.535705, 0.685794, 0.504156, 0.602791, 0.628167, 0.568111, 0.625537, 0.606724, 0.75, 0.518394, 0.573551, 0.577222, 0.802009, 0.608125, 0.584275, 0.751041, 0.684588, 0.546287, 0.661034, 0.837938, 0.768531, 0.519388, 0.644262, 0, 0.684513, 0.640202, 0.784768, 0.33588, 1, 0.498607, 0, 0.548405, 1, 0.662962, 0};

double isoTrackErr[84];

isoTrackErr[0] = 0.0050850785268;
isoTrackErr[1] = 0.025391769075;
isoTrackErr[2] = 0.067722243068;
isoTrackErr[3] = 0.10377832429;
isoTrackErr[4] = 0.21848399055;
isoTrackErr[5] = 0.006831792995;
isoTrackErr[6] = 0.026765783973;
isoTrackErr[7] = 0.070476377954;
isoTrackErr[8] = 0.10214497708;
isoTrackErr[9] = 0.021875646876;
isoTrackErr[10] = 0.024561638021;
isoTrackErr[11] = 0.058492207328;
isoTrackErr[12] = 0.12223033488;
isoTrackErr[13] = 1.356854731;
isoTrackErr[14] = 0.14582713615;
isoTrackErr[15] = 0.078578886874;
isoTrackErr[16] = 0.10512827005;
isoTrackErr[17] = 0.33921368274;
isoTrackErr[18] = 1.356854731;
isoTrackErr[19] = 0.46925458285;
isoTrackErr[20] = 0.67842736548;
isoTrackErr[21] = 0.0070625703729;
isoTrackErr[22] = 0.029065366792;
isoTrackErr[23] = 0.074436250981;
isoTrackErr[24] = 0.11346989569;
isoTrackErr[25] = 0.25746433808;
isoTrackErr[26] = 0.022968600857;
isoTrackErr[27] = 0.057503837057;
isoTrackErr[28] = 0.16188373019;
isoTrackErr[29] = 0.22614245516;
isoTrackErr[30] = 0.085061106842;
isoTrackErr[31] = 0.070520823175;
isoTrackErr[32] = 0.12619924503;
isoTrackErr[33] = 0.32351608395;
isoTrackErr[34] = 1.8;
isoTrackErr[35] = 0.67842736548;
isoTrackErr[36] = 0.67842736548;
isoTrackErr[37] = 0.018813319705;
isoTrackErr[38] = 0.043877590506;
isoTrackErr[39] = 0.10235501523;
isoTrackErr[40] = 0.21848399055;
isoTrackErr[41] = 0.047741152561;
isoTrackErr[42] = 0.085329724693;
isoTrackErr[43] = 0.15363447457;
isoTrackErr[44] = 0.1890564512;
isoTrackErr[45] = 0.093992342889;
isoTrackErr[46] = 0.18241924074;
isoTrackErr[47] = 0.46925458285;
isoTrackErr[48] = 0.021325293663;
isoTrackErr[49] = 0.064804110006;
isoTrackErr[50] = 0.14614296059;
isoTrackErr[51] = 0.38070289019;
isoTrackErr[52] = 0.025645377517;
isoTrackErr[53] = 0.054166118332;
isoTrackErr[54] = 0.12223033488;
isoTrackErr[55] = 0.099123160108;
isoTrackErr[56] = 0.1209665176;
isoTrackErr[57] = 0.37298610151;
isoTrackErr[58] = 0.023279923811;
isoTrackErr[59] = 0.06489050225;
isoTrackErr[60] = 0.13660848589;
isoTrackErr[61] = 0.30972485952;
isoTrackErr[62] = 0.030768067593;
isoTrackErr[63] = 0.075426859052;
isoTrackErr[64] = 0.16973026465;
isoTrackErr[65] = 0.056531331451;
isoTrackErr[66] = 0.19229119818;
isoTrackErr[67] = 0.46925458285;
isoTrackErr[68] = 0.26479636657;
isoTrackErr[69] = 0.15879914585;
isoTrackErr[70] = 0.062080635024;
isoTrackErr[71] = 0.10637044735;
isoTrackErr[72] = 0.45228491032;
isoTrackErr[73] = 0.087120602182;
isoTrackErr[74] = 0.1890564512;
isoTrackErr[75] = 0.14917313501;
isoTrackErr[76] = 0.21358743205;
isoTrackErr[77] = 1.356854731;
isoTrackErr[78] = 0.18241924074;
isoTrackErr[79] = 0.67842736548;
isoTrackErr[80] = 0.13311390862;
isoTrackErr[81] = 0.45228491032;
isoTrackErr[82] = 0.46925458285;
isoTrackErr[83] = 1.356854731;

  const std::string titre="CMS Supplementary";

  TCanvas *c1 = new TCanvas("c1", "c1",0,51,1920,1004);
  c1->SetFillColor(0);
  c1->cd();

  TH1D *h_mu_iso_eff = new TH1D("h_mu_iso_eff", "h_mu_iso_eff", 84 , 0 , 84);
  h_mu_iso_eff->SetTitle("");
  h_mu_iso_eff->SetXTitle("Search region bin number");
  h_mu_iso_eff->SetYTitle("iso track veto efficiency");
  h_mu_iso_eff->SetStats(0);
  h_mu_iso_eff->SetLineWidth(3);
  h_mu_iso_eff->SetMinimum(0.0);

  for (Int_t sbc=1;sbc<=84;++sbc) {
    h_mu_iso_eff->SetBinContent(sbc,isoTrackEff_SB[sbc-1]);
    h_mu_iso_eff->SetBinError(sbc,isoTrackErr[sbc-1]);
  }

  gStyle->SetPaintTextFormat("1.2f");
  h_mu_iso_eff->Draw();

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


  c1->SaveAs( "v5_isotrackvetoEff_tt.pdf" );

}
