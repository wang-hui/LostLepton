#include "CMS_lumi.h"

void Plot_isotrackvetoEff()
{

  double isoTrackEff_SB[84] = {0.590493, 0.561236, 0.695239, 0.644243, 0.456505, 0.589549, 0.590079, 0.605812, 0.640944, 0.567335, 0.620706, 0.553416, 0.590924, 0.52073, 0.457827, 0.639981, 0.601433, 0.694462, 0.623048, 0.544614, 0.597239, 0.575559, 0.575305, 0.637019, 0.639435, 0.656953, 0.569597, 0.624, 0.522413, 0.801731, 0.635865, 0.53689, 0.456575, 0.478029, 0.718798, 0.577477, 0.360426, 0.583126, 0.587228, 0.619752, 0.640475, 0.593927, 0.621907, 0.661975, 0.542692, 0.665093, 0.619586, 0.681133, 0.540293, 0.524838, 0.721528, 0.547977, 0.616599, 0.612173, 0.63301, 0.687242, 0.630832, 0.531038, 0.499387, 0.534385, 0.63721, 0.68155, 0.593246, 0.602999, 0.703103, 0.630304, 0.506049, 0.408792, 0.698086, 0.772164, 0.518214, 0.539126, 0.153339, 0.676463, 0.5882, 0.588093, 0.575513, 0.625436, 0.556933, 0.32517, 0.554402, 0.602495, 0.346128, 0.51129};

double isoTrackErr[84];

isoTrackErr[0] = 0.0030355236207;
isoTrackErr[1] = 0.0095990611522;
isoTrackErr[2] = 0.017589254484;
isoTrackErr[3] = 0.029451681831;
isoTrackErr[4] = 0.025842069373;
isoTrackErr[5] = 0.0042220123835;
isoTrackErr[6] = 0.010517122589;
isoTrackErr[7] = 0.017884725852;
isoTrackErr[8] = 0.025403657432;
isoTrackErr[9] = 0.0093519057732;
isoTrackErr[10] = 0.011671161591;
isoTrackErr[11] = 0.019960807908;
isoTrackErr[12] = 0.02972734912;
isoTrackErr[13] = 0.047612659855;
isoTrackErr[14] = 0.033500455575;
isoTrackErr[15] = 0.029057145171;
isoTrackErr[16] = 0.030875278685;
isoTrackErr[17] = 0.046247897914;
isoTrackErr[18] = 0.085557048661;
isoTrackErr[19] = 0.079908568382;
isoTrackErr[20] = 0.044042328808;
isoTrackErr[21] = 0.0040640009295;
isoTrackErr[22] = 0.010740085163;
isoTrackErr[23] = 0.018640586045;
isoTrackErr[24] = 0.025779779792;
isoTrackErr[25] = 0.03311381757;
isoTrackErr[26] = 0.012014548689;
isoTrackErr[27] = 0.023080446328;
isoTrackErr[28] = 0.041425540606;
isoTrackErr[29] = 0.048386140509;
isoTrackErr[30] = 0.027516624129;
isoTrackErr[31] = 0.029788035228;
isoTrackErr[32] = 0.041166162798;
isoTrackErr[33] = 0.054450646687;
isoTrackErr[34] = 0.094041677246;
isoTrackErr[35] = 0.078992098052;
isoTrackErr[36] = 0.068218350187;
isoTrackErr[37] = 0.015165453189;
isoTrackErr[38] = 0.031584731714;
isoTrackErr[39] = 0.059184949614;
isoTrackErr[40] = 0.12481031216;
isoTrackErr[41] = 0.013394856838;
isoTrackErr[42] = 0.021677648391;
isoTrackErr[43] = 0.036406460646;
isoTrackErr[44] = 0.048286961087;
isoTrackErr[45] = 0.014820533231;
isoTrackErr[46] = 0.019200409301;
isoTrackErr[47] = 0.032211443048;
isoTrackErr[48] = 0.0098838534385;
isoTrackErr[49] = 0.018915862407;
isoTrackErr[50] = 0.031395286661;
isoTrackErr[51] = 0.040201484052;
isoTrackErr[52] = 0.013557829641;
isoTrackErr[53] = 0.02264762426;
isoTrackErr[54] = 0.031710245788;
isoTrackErr[55] = 0.035444531952;
isoTrackErr[56] = 0.044143776548;
isoTrackErr[57] = 0.061911426366;
isoTrackErr[58] = 0.01067372206;
isoTrackErr[59] = 0.020559239771;
isoTrackErr[60] = 0.035254632145;
isoTrackErr[61] = 0.044620576516;
isoTrackErr[62] = 0.017204833394;
isoTrackErr[63] = 0.027516624129;
isoTrackErr[64] = 0.040169343928;
isoTrackErr[65] = 0.028734913507;
isoTrackErr[66] = 0.062044314554;
isoTrackErr[67] = 0.061100223984;
isoTrackErr[68] = 0.078024918783;
isoTrackErr[69] = 0.067224591949;
isoTrackErr[70] = 0.054197430254;
isoTrackErr[71] = 0.090616460245;
isoTrackErr[72] = 0.12489542248;
isoTrackErr[73] = 0.031963203762;
isoTrackErr[74] = 0.051025896027;
isoTrackErr[75] = 0.022057444731;
isoTrackErr[76] = 0.03309269208;
isoTrackErr[77] = 0.051779457135;
isoTrackErr[78] = 0.046062197559;
isoTrackErr[79] = 0.058447827815;
isoTrackErr[80] = 0.042892705907;
isoTrackErr[81] = 0.069484686791;
isoTrackErr[82] = 0.071710279109;
isoTrackErr[83] = 0.10088421809;

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


  c1->SaveAs( "v4_isotrackvetoEff.pdf" );

}
