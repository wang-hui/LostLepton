#include "CMS_lumi.h"

void Plot_isotrackvetoEff()
{
  double isoTrackEff_SB[84] = {0.595698, 0.559767, 0.690506, 0.598029, 0.483934, 0.588507, 0.591973, 0.644862, 0.646893, 0.55474, 0.622107, 0.526077, 0.556589, 0.51727344274, 0.59668, 0.618632, 0.599034, 0.63984638821, 0.72090258607, 0.74720555884, 0.68898843912, 0.578413, 0.576134, 0.618529, 0.656065, 0.653417, 0.581615, 0.673862, 0.539743, 0.61501362449, 0.649015, 0.58636, 0.483397, 0.47836, 0.82606110674, 0.718682, 0.84149043743, 0.578493, 0.572225, 0.591064, 0.597076, 0.616388, 0.622267, 0.645792, 0.546865, 0.647166, 0.653966, 0.71767282287, 0.531598, 0.494194, 0.627879, 0.6694277283, 0.627342, 0.647807, 0.695249, 0.646222, 0.735787, 0.58018408413, 0.492922, 0.517628, 0.589567, 0.68693268678, 0.540302, 0.641796, 0.714931, 0.586757, 0.436503, 0.65469310577, 0.60647246182, 0.723045, 0.508622, 0.455777, 0.74709949247, 0.704111, 0.614436, 0.542777, 0.443117, 0.56953212606, 0.647991, 0.72521427671, 0.406943, 0.51864099156, 0.57573498944, 0.54163554253};

double isoTrackErr[84];

isoTrackErr[0] = 0.0032675293939;
isoTrackErr[1] = 0.010235943347;
isoTrackErr[2] = 0.019306485726;
isoTrackErr[3] = 0.033511076395;
isoTrackErr[4] = 0.027691346059;
isoTrackErr[5] = 0.0044646359285;
isoTrackErr[6] = 0.010917287899;
isoTrackErr[7] = 0.018619774848;
isoTrackErr[8] = 0.026802257147;
isoTrackErr[9] = 0.0096153962283;
isoTrackErr[10] = 0.012055376461;
isoTrackErr[11] = 0.020678292933;
isoTrackErr[12] = 0.030615293193;
isoTrackErr[13] = 0.050371388614;
isoTrackErr[14] = 0.034867802637;
isoTrackErr[15] = 0.029079152987;
isoTrackErr[16] = 0.030844864188;
isoTrackErr[17] = 0.049860586809;
isoTrackErr[18] = 0.088934994237;
isoTrackErr[19] = 0.080501696832;
isoTrackErr[20] = 0.044967359616;
isoTrackErr[21] = 0.0044487495078;
isoTrackErr[22] = 0.011385236097;
isoTrackErr[23] = 0.020070425108;
isoTrackErr[24] = 0.027328543144;
isoTrackErr[25] = 0.034861539127;
isoTrackErr[26] = 0.013250062207;
isoTrackErr[27] = 0.02453850779;
isoTrackErr[28] = 0.041792229926;
isoTrackErr[29] = 0.052834755284;
isoTrackErr[30] = 0.029737284165;
isoTrackErr[31] = 0.031798689812;
isoTrackErr[32] = 0.047240986292;
isoTrackErr[33] = 0.05430731066;
isoTrackErr[34] = 0.089285443921;
isoTrackErr[35] = 0.078894617102;
isoTrackErr[36] = 0.073103312407;
isoTrackErr[37] = 0.016671863221;
isoTrackErr[38] = 0.033396373533;
isoTrackErr[39] = 0.062130318796;
isoTrackErr[40] = 0.13720555325;
isoTrackErr[41] = 0.014094370208;
isoTrackErr[42] = 0.022817320107;
isoTrackErr[43] = 0.038002427201;
isoTrackErr[44] = 0.051409969483;
isoTrackErr[45] = 0.015429628607;
isoTrackErr[46] = 0.019960592951;
isoTrackErr[47] = 0.032639622881;
isoTrackErr[48] = 0.011385660653;
isoTrackErr[49] = 0.021948759179;
isoTrackErr[50] = 0.036571498164;
isoTrackErr[51] = 0.047820324062;
isoTrackErr[52] = 0.015966939527;
isoTrackErr[53] = 0.026359001942;
isoTrackErr[54] = 0.036006172134;
isoTrackErr[55] = 0.038608882224;
isoTrackErr[56] = 0.051094485831;
isoTrackErr[57] = 0.069647212813;
isoTrackErr[58] = 0.012036688472;
isoTrackErr[59] = 0.023122330496;
isoTrackErr[60] = 0.041168591417;
isoTrackErr[61] = 0.050447876804;
isoTrackErr[62] = 0.019648102934;
isoTrackErr[63] = 0.030144776468;
isoTrackErr[64] = 0.044600864237;
isoTrackErr[65] = 0.034024555957;
isoTrackErr[66] = 0.069275869973;
isoTrackErr[67] = 0.072276700225;
isoTrackErr[68] = 0.085703470511;
isoTrackErr[69] = 0.07060906489;
isoTrackErr[70] = 0.062977158268;
isoTrackErr[71] = 0.11417714726;
isoTrackErr[72] = 0.12639978429;
isoTrackErr[73] = 0.034324995512;
isoTrackErr[74] = 0.057390856228;
isoTrackErr[75] = 0.024646537796;
isoTrackErr[76] = 0.035760833066;
isoTrackErr[77] = 0.059221019011;
isoTrackErr[78] = 0.055427009687;
isoTrackErr[79] = 0.0704316337;
isoTrackErr[80] = 0.051911144961;
isoTrackErr[81] = 0.086239170422;
isoTrackErr[82] = 0.091289794614;
isoTrackErr[83] = 0.10925437915;

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


  c1->SaveAs( "isotrackvetoEff.pdf" );

}
