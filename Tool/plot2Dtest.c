int plot2Dtest()
{
TCanvas* mycanvas = new TCanvas();
gStyle->SetOptStat(kFALSE);
gStyle->SetPaintTextFormat("4.4f");
TFile f1("v1_PredLL_gen_el.root");
TH2F *h1 = (TH2F*)f1.Get("h_pred_el_all_2d_met_njets");
h1->Rebin2D(20,4);
h1->GetYaxis()->SetRangeUser(4,12);
h1->GetXaxis()->SetRangeUser(200,1000);
double binC = h1->GetBinContent(5,3);
cout << binC << endl;
double binE = h1->GetBinError(5,3);
cout << binE << endl;

TFile f2("v1_PredLL_el.root");
TH2F *h2 = (TH2F*)f2.Get("h_pred_el_all_2d_met_njets");
h2->Rebin2D(20,4);
h2->GetYaxis()->SetRangeUser(4,12);
h2->GetXaxis()->SetRangeUser(200,1000);
binC = h2->GetBinContent(5,3);
cout << binC << endl;
binE = h2->GetBinError(5,3);
cout << binE << endl;

h1->Sumw2();
h2->Sumw2();
TH2F *rat = (TH2F*)h1->Clone();
rat->SetName("Ratio");
rat->Divide(h1,h2,1.,1.,"B");
//rat->Divide(h1,h2,1.,1.);
rat->SetMarkerSize(2);
rat->GetYaxis()->SetTitle("Njets");
rat->GetXaxis()->SetTitle("MET");
rat->DrawClone("Colztexte");
binC = rat->GetBinContent(5,3);
cout << binC << endl;
binE = rat->GetBinError(5,3);
cout << binE << endl;


/* double y_range = h1->GetMaximum();
  if (h1->GetMaximum() < h2->GetMaximum())
  {y_range = h2->GetMaximum();}
  y_range = y_range * 1.1;
  h1->GetYaxis()->SetRangeUser(0,y_range);
  h2->GetYaxis()->SetRangeUser(0,y_range);
*/

mycanvas->Print("gen_el_all_2D_ratio.png");

return 0;
}
