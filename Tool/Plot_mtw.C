{ 
  TFile *f = new TFile("RootForPlotting/123.root"); 

  //f.ls();
  TCanvas *c = new TCanvas("c","A Simple Graph Example",200,10,700,500); 
  c->SetLogy();
  gStyle->SetOptStat(0);

  TH1D *h1 = (TH1D*)f.Get("h_b_mtw_ttbar"); 
  h1->GetXaxis()->SetRangeUser(0,500);
  h1->GetXaxis()->SetTitle("MtW [GeV]");
  h1->SetLineColor(1);

  TH1D *h2 = (TH1D*)f.Get("h_b_mtw_t2tt_1");                
  h2->GetXaxis()->SetRangeUser(0,500);
  h2->GetXaxis()->SetTitle("MtW [GeV]");
  h2->SetLineColor(2);

  TH1D *h3 = (TH1D*)f.Get("h_b_mtw_t2tt_2");                
  h3->GetXaxis()->SetRangeUser(0,500);
  h3->GetXaxis()->SetTitle("MtW [GeV]");
  h3->SetLineColor(3);

  TH1D *h4 = (TH1D*)f.Get("h_b_mtw_t2tt_3");
  h4->GetXaxis()->SetRangeUser(0,500);
  h4->GetXaxis()->SetTitle("MtW [GeV]");
  h4->SetLineColor(4);

  TH1D *h5 = (TH1D*)f.Get("h_b_mtw_t2tt_4");
  h5->GetXaxis()->SetRangeUser(0,500);
  h5->GetXaxis()->SetTitle("MtW [GeV]");
  h5->SetLineColor(5);

  //TLine *line1 = new TLine(0.02,0,0.02,250);
  //line1->SetLineColor(2);
  //TLine *line2 = new TLine(0.2,0,0.2,100);
  //line2->SetLineColor(2);
  
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  h4->Draw("same");
  h5->Draw("same");

  //line1->Draw("same");
  //line2->Draw("same");

  const std::string titre="CMS Preliminary 2015";
  TLatex *title = new TLatex(0.09770115,0.9194915,titre.c_str());
  title->SetNDC();
  title->SetTextSize(0.045);
  title->Draw("same");

  TLegend* leg = new TLegend(0.6,0.75,0.9,0.87);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillColor(0);
  leg->AddEntry(h1,"TTBar","l");
  leg->AddEntry(h2,"T2tt_425_325","l");
  leg->AddEntry(h3,"T2tt_500_325","l");
  leg->AddEntry(h4,"T2tt_650_325","l");
  leg->AddEntry(h5,"T2tt_850_100","l");
  leg->Draw("same");
}
