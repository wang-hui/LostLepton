
void Plot_pred_preap()
{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Fri Jan 29 09:59:25 2016) by ROOT version6.02/05
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",260,124,538,327);
   Canvas_1->Range(-5.75,-2.241157,51.75,20.17042);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   gStyle->SetOptStat(0);

   TH1D *h_pred_lept_isotrk1 = new TH1D("h_pred_lept_isotrk1","",38,0,38);
   h_pred_lept_isotrk1->SetBinContent(1,15.41688);
   h_pred_lept_isotrk1->SetBinContent(2,2.689249);
   h_pred_lept_isotrk1->SetBinContent(3,1.770631);
   h_pred_lept_isotrk1->SetBinContent(5,3.025187);
   h_pred_lept_isotrk1->SetBinContent(6,5.676313);
   h_pred_lept_isotrk1->SetBinContent(7,1.22127);
   h_pred_lept_isotrk1->SetBinContent(9,0.8200502);
   h_pred_lept_isotrk1->SetBinContent(10,1.764622);
   h_pred_lept_isotrk1->SetBinContent(11,1.859962);
   h_pred_lept_isotrk1->SetBinContent(12,16.48773);
   h_pred_lept_isotrk1->SetBinContent(13,5.359042);
   h_pred_lept_isotrk1->SetBinContent(15,0.416388);
   h_pred_lept_isotrk1->SetBinContent(16,2.17177);
   h_pred_lept_isotrk1->SetBinContent(17,1.205399);
   h_pred_lept_isotrk1->SetBinContent(18,0.5830355);
   h_pred_lept_isotrk1->SetBinContent(20,0.4255685);
   h_pred_lept_isotrk1->SetBinContent(22,5.542776);
   h_pred_lept_isotrk1->SetBinContent(23,1.976001);
   h_pred_lept_isotrk1->SetBinContent(24,0.3044756);
   h_pred_lept_isotrk1->SetBinContent(25,3.625229);
   h_pred_lept_isotrk1->SetBinContent(26,2.021866);
   h_pred_lept_isotrk1->SetBinContent(27,1.036188);
   h_pred_lept_isotrk1->SetBinContent(29,0.5567768);
   h_pred_lept_isotrk1->SetBinContent(30,7.765275);
   h_pred_lept_isotrk1->SetBinContent(31,1.531608);
   h_pred_lept_isotrk1->SetBinContent(32,0.1719414);
   h_pred_lept_isotrk1->SetBinContent(33,3.271624);
   h_pred_lept_isotrk1->SetBinContent(34,0.5379648);
   h_pred_lept_isotrk1->SetBinError(1,2.613118);
   h_pred_lept_isotrk1->SetBinError(2,1.147244);
   h_pred_lept_isotrk1->SetBinError(3,1.076659);
   h_pred_lept_isotrk1->SetBinError(4,0.79866606148);
   h_pred_lept_isotrk1->SetBinError(5,1.175504);
   h_pred_lept_isotrk1->SetBinError(6,1.783231);
   h_pred_lept_isotrk1->SetBinError(7,0.7222653);
   h_pred_lept_isotrk1->SetBinError(8,0.94950725318);
   h_pred_lept_isotrk1->SetBinError(9,0.4752271);
   h_pred_lept_isotrk1->SetBinError(10,1.029543);
   h_pred_lept_isotrk1->SetBinError(11,1.091137);
   h_pred_lept_isotrk1->SetBinError(12,3.024688);
   h_pred_lept_isotrk1->SetBinError(13,1.681211);
   h_pred_lept_isotrk1->SetBinError(14,1.1494708018);
   h_pred_lept_isotrk1->SetBinError(15,0.416388);
   h_pred_lept_isotrk1->SetBinError(16,0.8928891);
   h_pred_lept_isotrk1->SetBinError(17,0.6985262);
   h_pred_lept_isotrk1->SetBinError(18,0.5830355);
   h_pred_lept_isotrk1->SetBinError(19,1.1476235882);
   h_pred_lept_isotrk1->SetBinError(20,0.4255685);
   h_pred_lept_isotrk1->SetBinError(21,1.2709769138);
   h_pred_lept_isotrk1->SetBinError(22,1.627558);
   h_pred_lept_isotrk1->SetBinError(23,0.835207);
   h_pred_lept_isotrk1->SetBinError(24,0.3044756);
   h_pred_lept_isotrk1->SetBinError(25,1.475125);
   h_pred_lept_isotrk1->SetBinError(26,0.9786622);
   h_pred_lept_isotrk1->SetBinError(27,0.6186725);
   h_pred_lept_isotrk1->SetBinError(28,0.96622376884);
   h_pred_lept_isotrk1->SetBinError(29,0.5567768);
   h_pred_lept_isotrk1->SetBinError(30,2.485114);
   h_pred_lept_isotrk1->SetBinError(31,0.8225499);
   h_pred_lept_isotrk1->SetBinError(32,0.1719414);
   h_pred_lept_isotrk1->SetBinError(33,1.113197);
   h_pred_lept_isotrk1->SetBinError(34,0.5379648);
   h_pred_lept_isotrk1->SetBinError(35,0.78565430477);
   h_pred_lept_isotrk1->SetBinError(36,0.82004298885);
   h_pred_lept_isotrk1->SetBinError(37,1.1692082201);
   h_pred_lept_isotrk1->SetEntries(37);

   h_pred_lept_isotrk1->GetXaxis()->SetTitle("Search Bins");
   h_pred_lept_isotrk1->SetMaximum(20.0);
   //h_pred_lept_isotrk1->SetMinimum(0.0);
   h_pred_lept_isotrk1->SetMinimum(-0.0001);
   h_pred_lept_isotrk1->SetLineWidth(3);

//   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
//   ptstats->SetName("stats");
//   ptstats->SetBorderSize(1);
//   ptstats->SetFillColor(0);
//   ptstats->SetTextAlign(12);
//   ptstats->SetTextFont(42);
//   TText *AText = ptstats->AddText("h_pred_lept_isotrk");
//   AText->SetTextSize(0.0368);
//   AText = ptstats->AddText("Entries = 45     ");
//   AText = ptstats->AddText("Mean  =  15.47");
//   AText = ptstats->AddText("RMS   =  12.31");
//   ptstats->SetOptStat(1111);
//   ptstats->SetOptFit(0);
//   ptstats->Draw();
//   h_pred_lept_isotrk1->GetListOfFunctions()->Add(ptstats);
//   ptstats->SetParent(h_pred_lept_isotrk1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   h_pred_lept_isotrk1->SetLineColor(ci);
   h_pred_lept_isotrk1->GetXaxis()->SetLabelFont(42);
   h_pred_lept_isotrk1->GetXaxis()->SetLabelSize(0.035);
   h_pred_lept_isotrk1->GetXaxis()->SetTitleSize(0.035);
   h_pred_lept_isotrk1->GetXaxis()->SetTitleFont(42);
   h_pred_lept_isotrk1->GetYaxis()->SetLabelFont(42);
   h_pred_lept_isotrk1->GetYaxis()->SetLabelSize(0.035);
   h_pred_lept_isotrk1->GetYaxis()->SetTitleSize(0.035);
   h_pred_lept_isotrk1->GetYaxis()->SetTitleFont(42);
   h_pred_lept_isotrk1->GetZaxis()->SetLabelFont(42);
   h_pred_lept_isotrk1->GetZaxis()->SetLabelSize(0.035);
   h_pred_lept_isotrk1->GetZaxis()->SetTitleSize(0.035);
   h_pred_lept_isotrk1->GetZaxis()->SetTitleFont(42);
   h_pred_lept_isotrk1->Draw("");

   //drawSBregionDef(0.0, 19.0, false);

   TLatex *   tex = new TLatex(0.09770115,0.9194915,"CMS Preliminary 2016, 2.3 fb^{-1}, #sqrt{s} = 13 TeV");
   tex->SetNDC();
   tex->SetTextSize(0.045);
   tex->SetLineWidth(2);
   tex->Draw();

   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
