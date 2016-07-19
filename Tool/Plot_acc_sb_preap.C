#include "CMS_lumi.h"

void Plot_acc_sb_preap()
{
  double err_acc[37];
  double err_acc_el[37];
  const double ttbar_mus_acc[37] = {0.833188,0.813335,0.830735,0.839931,0.846007,0.808412,0.791037,0.863016,0.856822,0.745926,0.741275,0.837696,0.819089,0.811669,0.83177,0.8392,0.802207,0.759466,0.838426,0.839024,0.78539,0.907118,0.898193,0.945786,0.849418,0.840731,0.897796,0.809444,0.850811,0.910869,0.910453,0.948242,0.857614,0.842478,0.902275,0.830653,0.836152};

  const double ttbar_els_acc[37] = {0.847301,0.81815,0.823032,0.849582,0.847771,0.797805,0.818129,0.855029,0.870859,0.814914,0.759492,0.831736,0.821877,0.799215,0.751261,0.844663,0.811157,0.792653,0.765419,0.820682,0.701654,0.91088,0.905967,0.931121,0.843057,0.848582,0.895733,0.882179,0.828121,0.903091,0.921523,0.930069,0.844546,0.839703,0.863359,0.875123,0.788378};

err_acc[0] = 0.0033252876063;
err_acc[1] = 0.0067599878464;
err_acc[2] = 0.012089845381;
err_acc[3] = 0.02333468668;
err_acc[4] = 0.006744380329;
err_acc[5] = 0.0089445130746;
err_acc[6] = 0.015945966545;
err_acc[7] = 0.02914241673;
err_acc[8] = 0.017842694331;
err_acc[9] = 0.018407631864;
err_acc[10] = 0.02316708597;
err_acc[11] = 0.0032491967865;
err_acc[12] = 0.0063191529388;
err_acc[13] = 0.011356334581;
err_acc[14] = 0.02293372401;
err_acc[15] = 0.008933290665;
err_acc[16] = 0.01177060416;
err_acc[17] = 0.019105386469;
err_acc[18] = 0.03000528289;
err_acc[19] = 0.020466186058;
err_acc[20] = 0.035929766823;
err_acc[21] = 0.0036228774134;
err_acc[22] = 0.0077729903308;
err_acc[23] = 0.010481042898;
err_acc[24] = 0.0074724697595;
err_acc[25] = 0.011603837296;
err_acc[26] = 0.015245992461;
err_acc[27] = 0.025447113373;
err_acc[28] = 0.020870663489;
err_acc[29] = 0.0038995968456;
err_acc[30] = 0.0081437906046;
err_acc[31] = 0.011393920452;
err_acc[32] = 0.0083704581724;
err_acc[33] = 0.01407686086;
err_acc[34] = 0.016553661421;
err_acc[35] = 0.02814771485;
err_acc[36] = 0.02290746571;
err_acc_el[0] = 0.0033327246607;
err_acc_el[1] = 0.006690741288;
err_acc_el[2] = 0.012382618751;
err_acc_el[3] = 0.022435646439;
err_acc_el[4] = 0.0062952683646;
err_acc_el[5] = 0.0089453313795;
err_acc_el[6] = 0.014491817514;
err_acc_el[7] = 0.030130021412;
err_acc_el[8] = 0.018305872229;
err_acc_el[9] = 0.016517732116;
err_acc_el[10] = 0.021666726668;
err_acc_el[11] = 0.0032891629494;
err_acc_el[12] = 0.0064493415074;
err_acc_el[13] = 0.011773424644;
err_acc_el[14] = 0.023451554913;
err_acc_el[15] = 0.0088607604077;
err_acc_el[16] = 0.011803412898;
err_acc_el[17] = 0.017906571206;
err_acc_el[18] = 0.03181950751;
err_acc_el[19] = 0.021029153012;
err_acc_el[20] = 0.035531027944;
err_acc_el[21] = 0.0036594484786;
err_acc_el[22] = 0.0076775132589;
err_acc_el[23] = 0.012234102026;
err_acc_el[24] = 0.0076411030843;
err_acc_el[25] = 0.011445341113;
err_acc_el[26] = 0.014889646039;
err_acc_el[27] = 0.024191449634;
err_acc_el[28] = 0.021609389986;
err_acc_el[29] = 0.0041457295382;
err_acc_el[30] = 0.0079653292382;
err_acc_el[31] = 0.012489920279;
err_acc_el[32] = 0.0086735649021;
err_acc_el[33] = 0.01332465423;
err_acc_el[34] = 0.018418906272;
err_acc_el[35] = 0.026835303461;
err_acc_el[36] = 0.025125889231;

  const std::string titre="CMS Supplementary";

  TCanvas *c1 = new TCanvas("c1", "c1",0,51,1920,1004);
  c1->SetFillColor(0);
  c1->cd();

  TH1D *h_mu_iso_eff = new TH1D("h_mu_iso_eff", "h_mu_iso_eff", 37 , 0 , 37);
  h_mu_iso_eff->SetTitle("");
  h_mu_iso_eff->SetXTitle("Search region bin number");
  h_mu_iso_eff->SetYTitle("acceptance");
  h_mu_iso_eff->SetStats(0);
  h_mu_iso_eff->SetLineWidth(3);
  h_mu_iso_eff->SetMinimum(0.0);

  for (Int_t sbc=1;sbc<=37;++sbc)
  {
    h_mu_iso_eff->SetBinContent(sbc,ttbar_mus_acc[sbc-1]);
    h_mu_iso_eff->SetBinError(sbc,err_acc[sbc-1]);
    //h_mu_iso_eff->SetBinContent(sbc,sbc);
  }

  gStyle->SetPaintTextFormat("1.2f");
  h_mu_iso_eff->Draw();

  TH1D *h_e_iso_eff = new TH1D("h_e_iso_eff", "h_e_iso_eff", 37 , 0 , 37);
  h_e_iso_eff->SetLineColor(2);
  h_e_iso_eff->SetLineWidth(3);

  for (Int_t sbc=1;sbc<=37;++sbc) 
  {
    h_e_iso_eff->SetBinContent(sbc,ttbar_els_acc[sbc-1]);
    h_e_iso_eff->SetBinError(sbc,err_acc_el[sbc-1]);
    //h_mu_iso_eff->SetBinContent(sbc,sbc);
  }
  h_e_iso_eff->Draw("same");

//  TLatex *title = new TLatex(0.09770115,0.9194915,titre.c_str());
//  title->SetNDC();
//  title->SetTextSize(0.045);
//  title->Draw("same");
//  TLatex *tex2 = new TLatex(0.8127615,0.9178645,"(13 TeV)");
//  tex2->SetNDC();
//  tex2->SetTextSize(0.045);
//  tex2->SetLineWidth(2);
//  tex2->Draw();
  CMSStylePlot::CMS_lumi( c1, 0, 0 );


  TLegend* leg = new TLegend(0.6,0.15,0.9,0.27);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillColor(0);
  leg->AddEntry(h_mu_iso_eff,"Muons","l");
  leg->AddEntry(h_e_iso_eff,"Electrons","l");
  leg->Draw("same");

  c1->SaveAs( "h_isotrackvetoEff.png" );
}
