#include "CMS_lumi.h"

void Plot_isotrackvetoEff()
{
double isoTrackEff_SB[84] = {0.594375, 0.546478, 0.7282, 0.688987, 0.595184, 0.591746, 0.630295, 0.563533, 0.614842, 0.633231, 0.599336, 0.608467, 0.574718, 0.49790383114, 0.519775, 0.697071, 0.524727, 0.444574, 0.79088364396, 0.79088364396, 0.79088364396, 0.592583, 0.577918, 0.603168, 0.440157, 0.577365, 0.591467, 0.482023, 0.643849, 0.63944231766, 0.6993, 0.570809, 0.45026, 0.460367, 0.71155236655 , 0.71155236655, 0.71155236655, 0.568682, 0.627962, 0.71398, 0.6915667528, 0.541915, 0.525106, 0.730244, 0.710579, 0.709875, 0.653734, 0.864201, 0.546497, 0.545002, 0.836265, 0.80460967188, 0.636125, 0.567026, 0.441162, 0.529601, 0.693581, 0.70000631589, 0.526511, 0.548097, 0.50195, 0.52638378831, 0.561576, 0.585059, 0.609084, 0.561132, 0.593056, 0.56913185298, 0.906777, 0.475208, 0.490372, 0.660005, 0.5, 0.501237, 0.64193, 0.700626, 0.606326, 0.55556847701, 0.578813, 0.61165800818, 0.621471, 0.59577751302, 0.604189, 0.55815360214};

double isoTrackErr[84];
isoTrackErr[0] = 0.0047617046468;
isoTrackErr[1] = 0.020504840014;
isoTrackErr[2] = 0.045746934432;
isoTrackErr[3] = 0.060171579346;
isoTrackErr[4] = 0.059189652946;
isoTrackErr[5] = 0.0059259518815;
isoTrackErr[6] = 0.020484422649;
isoTrackErr[7] = 0.040790801523;
isoTrackErr[8] = 0.062236648879;
isoTrackErr[9] = 0.014992538047;
isoTrackErr[10] = 0.017500369138;
isoTrackErr[11] = 0.028320976673;
isoTrackErr[12] = 0.051105242061;
isoTrackErr[13] = 0.058045087455;
isoTrackErr[14] = 0.067566914065;
isoTrackErr[15] = 0.041899036557;
isoTrackErr[16] = 0.042570906097;
isoTrackErr[17] = 0.074664344795;
isoTrackErr[18] = 0.055465807665;
isoTrackErr[19] = 0.055465807665;
isoTrackErr[20] = 0.055465807665;
isoTrackErr[21] = 0.0067496172598;
isoTrackErr[22] = 0.025245535039;
isoTrackErr[23] = 0.054421027237;
isoTrackErr[24] = 0.089349441538;
isoTrackErr[25] = 0.098483843292;
isoTrackErr[26] = 0.020084140666;
isoTrackErr[27] = 0.042780782671;
isoTrackErr[28] = 0.10012215508;
isoTrackErr[29] = 0.083585738567;
isoTrackErr[30] = 0.063128451837;
isoTrackErr[31] = 0.055927489941;
isoTrackErr[32] = 0.070722984586;
isoTrackErr[33] = 0.10088421809;
isoTrackErr[34] = 0.084607176399;
isoTrackErr[35] = 0.084607176399;
isoTrackErr[36] = 0.084607176399;
isoTrackErr[37] = 0.016734810319;
isoTrackErr[38] = 0.035639317182;
isoTrackErr[39] = 0.091730235211;
isoTrackErr[40] = 0.083455452178;
isoTrackErr[41] = 0.04211230423;
isoTrackErr[42] = 0.068885277539;
isoTrackErr[43] = 0.10857592821;
isoTrackErr[44] = 0.12049747516;
isoTrackErr[45] = 0.060197559742;
isoTrackErr[46] = 0.10228384065;
isoTrackErr[47] = 0.14767475903;
isoTrackErr[48] = 0.021827644143;
isoTrackErr[49] = 0.06227544862;
isoTrackErr[50] = 0.10833011818;
isoTrackErr[51] = 0.088934994237;
isoTrackErr[52] = 0.024887298837;
isoTrackErr[53] = 0.046415145194;
isoTrackErr[54] = 0.1009918555;
isoTrackErr[55] = 0.09274668227;
isoTrackErr[56] = 0.093361883806;
isoTrackErr[57] = 0.075234675626;
isoTrackErr[58] = 0.022757674081;
isoTrackErr[59] = 0.062402245208;
isoTrackErr[60] = 0.13354323019;
isoTrackErr[61] = 0.1296285281;
isoTrackErr[62] = 0.030009213409;
isoTrackErr[63] = 0.070476377954;
isoTrackErr[64] = 0.13660848589;
isoTrackErr[65] = 0.055728041731;
isoTrackErr[66] = 0.18108939694;
isoTrackErr[67] = 0.13311390862;
isoTrackErr[68] = 0.14767475903;
isoTrackErr[69] = 0.17325608641;
isoTrackErr[70] = 0.04887242706;
isoTrackErr[71] = 0.12304085569;
isoTrackErr[72] = 0.2867016056;
isoTrackErr[73] = 0.093361734968;
isoTrackErr[74] = 0.16188373019;
isoTrackErr[75] = 0.11710383051;
isoTrackErr[76] = 0.18108939694;
isoTrackErr[77] = 0.17325608641;
isoTrackErr[78] = 0.17099598839;
isoTrackErr[79] = 0.14582713615;
isoTrackErr[80] = 0.12959976516;
isoTrackErr[81] = 0.11385852291;
isoTrackErr[82] = 0.28171901436;
isoTrackErr[83] = 0.20278085194;

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


  c1->SaveAs( "v5_isotrackvetoEff.pdf" );

}
