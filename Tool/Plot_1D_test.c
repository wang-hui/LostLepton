int Plot_1D_test()
{
TCanvas* mycanvas = new TCanvas();
gStyle->SetOptStat(kFALSE);
TFile f1("CalLL.root");

if (true)
{
TH1D *h1 = (TH1D*)f1.Get("h_iso_track_veto_nb");
//h1->Rebin2D(20,4);
h1->GetXaxis()->SetRangeUser(0,5);

TH1D *h2 = (TH1D*)f1.Get("h_iso_track_veto_nb_all");
//h2->Rebin2D(20,4);
//h2->GetYaxis()->SetRangeUser(4,12);
//h2->GetXaxis()->SetRangeUser(200,1000);

h1->Sumw2();
h2->Sumw2();
TH1D *rat = (TH1D*)h1->Clone();
rat->Divide(h1,h2,1.,1.,"B");
//rat->Divide(h1,h2,1.,1.);
rat->SetMarkerSize(2);
rat->GetYaxis()->SetTitle("iso reack veto eff");
rat->GetXaxis()->SetTitle("number of bottom");
rat->GetYaxis()->SetRangeUser(0,1);
rat->Draw();

mycanvas->Print("iso_track_veto_nb.png");
}
if (true)
{
TH1D *h1 = (TH1D*)f1.Get("h_iso_track_veto_nt");
//h1->Rebin2D(20,4);
h1->GetXaxis()->SetRangeUser(0,5);

TH1D *h2 = (TH1D*)f1.Get("h_iso_track_veto_nt_all");
//h2->Rebin2D(20,4);
//h2->GetYaxis()->SetRangeUser(4,12);
//h2->GetXaxis()->SetRangeUser(200,1000);

h1->Sumw2();
h2->Sumw2();
TH1D *rat = (TH1D*)h1->Clone();
//rat->SetName("Ratio");
rat->Divide(h1,h2,1.,1.,"B");
//rat->Divide(h1,h2,1.,1.);
rat->SetMarkerSize(2);
rat->GetYaxis()->SetTitle("iso reack veto eff");
rat->GetXaxis()->SetTitle("number of top");
rat->GetYaxis()->SetRangeUser(0,1);
rat->Draw();

mycanvas->Print("iso_track_veto_nt.png");
}
if (true)
{
TH1D *h1 = (TH1D*)f1.Get("h_iso_track_veto_ht");
//h1->Rebin2D(20,4);
//h1->GetXaxis()->SetRangeUser(200,1000);

TH1D *h2 = (TH1D*)f1.Get("h_iso_track_veto_ht_all");
//h2->Rebin2D(20,4);
//h2->GetYaxis()->SetRangeUser(4,12);
//h2->GetXaxis()->SetRangeUser(200,1000);

h1->Sumw2();
h2->Sumw2();
TH1D *rat = (TH1D*)h1->Clone();
//rat->SetName("Ratio");
rat->Divide(h1,h2,1.,1.,"B");
//rat->Divide(h1,h2,1.,1.);
rat->SetMarkerSize(2);
rat->GetYaxis()->SetTitle("iso reack veto eff");
rat->GetXaxis()->SetTitle("ht");
rat->GetYaxis()->SetRangeUser(0,1);
rat->Draw();

mycanvas->Print("iso_track_veto_ht.png");
}
if (true)
{
TH1D *h1 = (TH1D*)f1.Get("h_iso_track_veto_mt2");
//h1->Rebin2D(20,4);
//h1->GetXaxis()->SetRangeUser(200,1000);

TH1D *h2 = (TH1D*)f1.Get("h_iso_track_veto_mt2_all");
//h2->Rebin2D(20,4);
//h2->GetYaxis()->SetRangeUser(4,12);
//h2->GetXaxis()->SetRangeUser(200,1000);

h1->Sumw2();
h2->Sumw2();
TH1D *rat = (TH1D*)h1->Clone();
//rat->SetName("Ratio");
rat->Divide(h1,h2,1.,1.,"B");
//rat->Divide(h1,h2,1.,1.);
rat->SetMarkerSize(2);
rat->GetYaxis()->SetTitle("iso reack veto eff");
rat->GetXaxis()->SetTitle("mt2");
rat->GetYaxis()->SetRangeUser(0,1);
rat->Draw();

mycanvas->Print("iso_track_veto_mt2.png");
}
if (true)
{
TH1D *h1 = (TH1D*)f1.Get("h_iso_track_veto_met");
//h1->Rebin2D(20,4);
//h1->GetXaxis()->SetRangeUser(200,1000);

TH1D *h2 = (TH1D*)f1.Get("h_iso_track_veto_met_all");
//h2->Rebin2D(20,4);
//h2->GetYaxis()->SetRangeUser(4,12);
//h2->GetXaxis()->SetRangeUser(200,1000);

h1->Sumw2();
h2->Sumw2();
TH1D *rat = (TH1D*)h1->Clone();
//rat->SetName("Ratio");
rat->Divide(h1,h2,1.,1.,"B");
//rat->Divide(h1,h2,1.,1.);
rat->SetMarkerSize(2);
rat->GetYaxis()->SetTitle("iso reack veto eff");
rat->GetXaxis()->SetTitle("met");
rat->GetYaxis()->SetRangeUser(0,1);
rat->Draw();

mycanvas->Print("iso_track_veto_met.png");
}

return 0;
}
