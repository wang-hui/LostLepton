#include <vector>
#include <iostream>
#include <string>
#include <cstring>

#include "TFile.h"
#include "TList.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TExec.h"
#include "TLine.h"

#include "SusyAnaTools/Tools/searchBins.h"

class ClosurePlots
{
 public:
  TFile * fin;
  TFile * fin2;
  TList * list;
  TList * list2;

  double scale;

  void Initialization(); 
  void PrintPlotsName();
  void SetScale(
               double Nevents,
               double XSec,
               double Lumi
               );
  void ClosureTemplate(
                       TString hist_tag,
                       TString XTitle,
                       double min,
                       double max
                       );
  void DiLeptonPlots(
                     TString SLhist,
                     TString DLhist,
                     TString XTitle,
                     double min,
                     double max
                    );

};

void ClosurePlots::Initialization()
{
  //fin = TFile::Open("RootForPlotting/test.root");
  //fin = TFile::Open("RootForPlotting/v151201_0p5fb_ExpLL.root");
  //fin = TFile::Open("RootForPlotting/v151201_2p1fb_ExpLL.root");
  //fin = TFile::Open("RootForPlotting/v151216_ExpLL.root");
  //fin = TFile::Open("RootForPlotting/v160106_ttbar_ExpLL.root");
  //fin = TFile::Open("RootForPlotting/v160112_notrigeff_ExpLL.root");
  //fin = TFile::Open("RootForPlotting/v160113_notrigeff_ExpLL.root");
  //fin = TFile::Open("RootForPlotting/v160217_ttbar_ExpLL.root");
  //fin = TFile::Open("RootForPlotting/v160217_ttbar_v3_ExpLL.root");
  //fin = TFile::Open("RootForPlotting/v160217_ttbar_v3_invertedDPhi_ExpLL.root");
  //fin = TFile::Open("RootForPlotting/v160120_invertedDPhi_ttbar_ExpLL.root");
  //fin = TFile::Open("RootForPlotting/v160113_nodilepton_ExpLL.root");
  //fin = TFile::Open("RootForPlotting/v160309_ttbarSingletopW_45bins_ExpLL.root");
  fin = TFile::Open("RootForPlotting/v20160719_HuaTest_ExpLL.root");
  //fin2 = TFile::Open("RootForPlotting/v151204_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v151209_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v151216_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/ttbar_v151217_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160106_ttbar_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160106_data_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160106_data_trigSel_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160112_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160112_nomtw_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160112_nodilepton_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160114_accSB_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160114_ttbar_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160116_data_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160217_ttbar_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160217_ttbar_v3_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160309_ttbarSingletopW_45bins_PredLL.root");
  fin2 = TFile::Open("RootForPlotting/v20160719_HuaTest_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160217_ttbar_v3_invertedDPhi_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160120_invertedDPhi_ttbar_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160122_accfromInvertedDPhi_PredLL.root");
  //fin2 = TFile::Open("RootForPlotting/v160122_efffromInvertedDPhi_PredLL.root");
  list = fin->GetListOfKeys();
  list2 = fin2->GetListOfKeys();
  scale=1.0;
}

void ClosurePlots::PrintPlotsName()
{
  for(int i  = 0 ; i < list->GetSize() ; i++)
  {
    std::cout<<"Name: "<< list->At(i)->GetName() << "("<< i <<")"<<std::endl;
  }
  
  return ;
}


void ClosurePlots::SetScale(
                            double Nevents,
                            double XSec,
                            double Lumi
                           )
{
  scale = XSec*Lumi/Nevents;
}



void ClosurePlots::DiLeptonPlots(
                                 TString SLhist,
                                 TString DLhist,
                                 TString XTitle,
                                 double min,
                                 double max
                                )
{
  TH1D * h_exp_sl = (TH1D*)fin->Get(SLhist)->Clone();
  TH1D * h_exp_dl = (TH1D*)fin->Get(DLhist)->Clone();

  TCanvas *c = new TCanvas("c","A Simple Graph Example",200,10,700,500);
  gStyle->SetOptStat(0);

  h_exp_dl->GetXaxis()->SetRangeUser(min,max);
  h_exp_dl->GetXaxis()->SetTitle(XTitle);
  h_exp_dl->SetLineColor(1);
  h_exp_dl->SetLineWidth(3);
  h_exp_dl->Sumw2();
  h_exp_dl->Scale(scale);

  h_exp_sl->GetXaxis()->SetRangeUser(min,max);
  h_exp_sl->SetLineColor(2);
  h_exp_sl->SetLineWidth(3);
  h_exp_sl->Sumw2();
  h_exp_sl->Scale(scale);

  h_exp_dl->Draw();
  h_exp_sl->Draw("same");

  const std::string titre="CMS Preliminary 2016, 2.3 fb^{-1}, #sqrt{s} = 13 TeV";
  TLatex *title = new TLatex(0.09770115,0.9194915,titre.c_str());
  title->SetNDC();
  title->SetTextSize(0.045);
  title->Draw("same");

  TLegend* leg = new TLegend(0.6,0.75,0.85,0.85);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillColor(0);
  leg->AddEntry(h_exp_dl,"Expectation 1 or 2 lepton","l");
  leg->AddEntry(h_exp_sl,"Expectation 1 lepton","l");
  leg->Draw("same");

  c->SaveAs( TString("Plotsc/") + DLhist + TString("_compare.png") );
  c->SaveAs( TString("Plotsc/") + DLhist + TString("_compare.C") );
}

void ClosurePlots::ClosureTemplate(
                                   TString hist_tag,
                                   TString XTitle,
                                   double min,
                                   double max
                                  )
{ 
  TH1D * h_pred;
  TH1D * h_exp;

  for(int i  = 0 ; i < list->GetSize() ; i++)
  {
    if( TString(list->At(i)->GetName()).Contains( hist_tag ) )
    {
      if( TString(list->At(i)->GetName()).Contains( "_exp_" ) )
      {
        h_exp = (TH1D*)fin->Get(list->At(i)->GetName())->Clone();
      }
    }
  }
  for(int i  = 0 ; i < list2->GetSize() ; i++)
  {
    if( TString(list2->At(i)->GetName()).Contains( hist_tag ) )
    {
      if( TString(list2->At(i)->GetName()).Contains( "_pred_" ) )
      {
        h_pred = (TH1D*)fin2->Get(list2->At(i)->GetName())->Clone();
      }
    }
  }

  //Set Style for h_exp and h_pred
  h_pred->SetMarkerStyle(20);
  h_pred->SetMarkerColor(kBlue);
  h_pred->SetMarkerSize(0);
  h_pred->SetLineColor(h_pred->GetMarkerColor());
  h_pred->SetLineWidth(3);
  h_pred->GetXaxis()->SetRangeUser(min,max);
  h_pred->GetXaxis()->SetTitle(XTitle);
  h_pred->GetYaxis()->SetTitleOffset(0.6);
  h_pred->GetYaxis()->SetTitleFont(42);
  h_pred->GetYaxis()->SetTitleSize(0.065);
  h_pred->GetYaxis()->SetLabelSize(0.04);
  h_pred->GetYaxis()->SetLabelFont(42);
  h_pred->GetYaxis()->SetTitle("Events");
  //h_pred->Sumw2();
  //h_pred->Scale(scale);

  h_exp->SetMarkerStyle(20);
  h_exp->SetMarkerColor(kRed);
  h_exp->SetMarkerSize(0.9);
  h_exp->SetLineColor(h_exp->GetMarkerColor());
  h_exp->SetLineWidth(3);
  h_exp->GetXaxis()->SetRangeUser(min,max);
  h_exp->GetXaxis()->SetTitle(XTitle);
  h_exp->GetYaxis()->SetTitleOffset(0.6);
  h_exp->GetYaxis()->SetTitleFont(42);
  h_exp->GetYaxis()->SetTitleSize(0.065);
  h_exp->GetYaxis()->SetLabelSize(0.04);
  h_exp->GetYaxis()->SetLabelFont(42);
  h_exp->GetYaxis()->SetTitle("Events");
  //h_exp->Sumw2();
  //h_exp->Scale(scale);

  // Ratio plots
  TH1* h_ratioFrame;
  TH1* h_ratio;
    
  h_ratio = static_cast<TH1*>(h_exp->Clone("Ratio"));
  h_ratio->Divide(h_pred);
  h_ratio->SetMarkerSize(1);
  h_ratio->GetYaxis()->SetTitle("#frac{Direct}{Prediction}");
  //h_ratio->GetYaxis()->SetRangeUser(0.0,5.1);
  h_ratio->SetTitle("");
  h_ratio->SetStats(0);
  h_ratio->SetLineWidth(1);
  h_ratio->GetYaxis()->SetTitleSize(0.15);
  h_ratio->GetYaxis()->SetTitleOffset(0.3);
  h_ratio->GetYaxis()->SetTitleFont(42);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetLabelFont(42);
  h_ratio->GetXaxis()->SetLabelOffset(0.01);
  h_ratio->GetXaxis()->SetLabelFont(42);
  h_ratio->GetXaxis()->SetLabelSize(0.08);
  h_ratio->GetXaxis()->SetTitleSize(0.16);
  h_ratio->GetXaxis()->SetTitleFont(42);
  h_ratio->GetXaxis()->SetTitleOffset(0.6);
  if ( hist_tag.Contains("_sb") )
  {
    h_ratio->GetXaxis()->SetTitle("Search region bin number");
  }

  //Create LUMI stamp
  //const std::string titre="CMS Preliminary 2015, "+ lumi_str + " fb^{-1}, #sqrt{s} = 13 TeV";
  //const std::string titre="CMS Preliminary 2016, 2.3 fb^{-1}, #sqrt{s} = 13 TeV";
  const std::string titre="CMS Supplementary                                                             2.3 fb^{-1}(13 TeV)";

  TLatex *title = new TLatex(0.09770115,0.9194915,titre.c_str());
  title->SetNDC();
  title->SetTextSize(0.045);

  //Create Legend
  TLegend* leg = new TLegend(0.55,0.75,0.90,0.90);
  leg->SetBorderSize(1);
  leg->SetLineColor(1);
  leg->SetLineWidth(2);
  leg->SetFillColor(0);
  //leg->SetFillStyle();
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader("LostLepton background");
  leg->AddEntry(h_exp,"Direct from simulation","P");
  //leg->AddEntry(hPred[0],"Treat simulation like data","L");
  leg->AddEntry(h_pred,"Treat simulation like data");

  //Draw plots on Canvas
  TCanvas *c = new TCanvas("c","",50,50,800,600); 
  //HistStyle::init();
  gStyle->SetOptStat(0);

  TPad *pad = (TPad*) c->GetPad(0); 
  pad->Clear();
  pad->Divide(2, 1);

  double divRatio = 0.20;
  double labelRatio = (1-divRatio)/divRatio;
  double small = 0;

  pad->cd(1); 
  TPad *pad1 = (TPad*) pad->GetPad(1);
  pad1->SetPad("", "", 0, divRatio, 1.0, 1.0, kWhite);
  pad1->SetBottomMargin(0.005);
  pad1->SetBorderMode(0);
  
  TExec *setex = new TExec("setex", "gStyle->SetErrorX(0.0)");

  if( hist_tag.Contains("_sb") )
  { 
    pad1->SetLogy(); 
    h_exp->GetXaxis()->SetRangeUser(0.,45);
    //h_exp->GetYaxis()->SetRangeUser(0.,100);

    Double_t pred,exp,pred_err,exp_err;
    double non_closure_unc[45] ={-10};
    for (Int_t i = 1; i < h_pred->GetNbinsX(); i++)
    {
      pred = h_pred->GetBinContent(i);
      exp = h_exp->GetBinContent(i);
      pred_err = h_pred->GetBinError(i);
      exp_err = h_exp->GetBinError(i);
      //std::cout << "i: " << i << " pred_err: " << pred_err << " exp_err: " << exp_err << std::endl;
      double r = 100;
      double e = 5;
      if ( (pred > 0) && (exp > 0) ) 
      { 
        r = exp/pred;
        e = std::sqrt( exp_err*exp_err + pred_err*pred_err*r*r ) / pred;
        non_closure_unc[i-1] = e;
        //double percent = (exp-pred)/pred * 100;
        std::cout << "i: " << i << " Pred: "<< pred << " Pred_Err: "<< pred_err << std::endl;
        //std::cout << "i: " << i << " Pred: "<< pred << " Exp: "<< exp << " Ratio: " << r-1 << " Error: " << e << std::endl;
      }
    }
  }
  
  setex->Draw();
  h_exp->Draw("PE1");
  h_pred->DrawCopy("hist same");
  h_pred->SetFillColor(kBlue-4);
  h_pred->SetFillStyle(3001);
  h_pred->Draw("E2 same");

  SearchBins theSearchBins("SB_69_2016");
  if( hist_tag.Contains("_sb") ){ theSearchBins.drawSBregionDef(0.0, 75.0, true); }
  title->Draw("same");
  leg->Draw("same");

  c->Update(); 
  
  pad->cd(2);
  TPad *pad2 = (TPad*) pad->GetPad(2);
  pad2->SetPad("ratio", "", 0, 0, 1.0, divRatio, kWhite);
  pad2->SetBottomMargin(0.3);
  pad2->SetTopMargin(small);
  pad2->SetBorderMode(0);

  h_ratio->SetMaximum(1.8);
  h_ratio->SetMinimum(0.4);
  if( hist_tag.Contains("_sb") )
  {
    h_ratio->GetXaxis()->SetLabelSize(0.1);
    h_ratio->GetXaxis()->SetTitleOffset(0.8);
    h_ratio->GetXaxis()->SetRangeUser(0.,45);
  }

  TLine *tl_one = new TLine();
  tl_one->SetLineStyle(2);
  tl_one->SetLineColor(1);
  tl_one->SetLineWidth(2);
  
  h_ratio->Draw("PE1");
  tl_one->DrawLine(min,1.,max,1.);

  c->SaveAs( TString("Plotsc/") + hist_tag + TString(".png") );
  c->SaveAs( TString("Plotsc/") + hist_tag + TString(".pdf") );
  c->SaveAs( TString("Plotsc/") + hist_tag + TString(".C") );

  /*
  TCanvas *c = new TCanvas("c","A Simple Graph Example",200,10,700,500); 
  gStyle->SetOptStat(0);

  TPad *pad = (TPad*) c->GetPad(0); 
  pad->Clear();
  pad->Divide(2, 1);

  double divRatio = 0.20;
  double labelRatio = (1-divRatio)/divRatio;
  double small = 0;

  pad->cd(1); 
  TPad *pad1 = (TPad*) pad->GetPad(1);
  pad1->SetPad("", "", 0, divRatio, 1.0, 1.0, kWhite);
  pad1->SetBottomMargin(0.005);
  pad1->SetBorderMode(0);

  h_pred->GetXaxis()->SetRangeUser(min,max);
  h_pred->GetXaxis()->SetTitle(XTitle);
  h_pred->SetLineColor(1);
  h_pred->SetLineWidth(3);
  //h_pred->Sumw2();
  h_pred->Scale(scale);

  h_exp->GetXaxis()->SetRangeUser(min,max);
  h_exp->GetXaxis()->SetTitle(XTitle);
  h_exp->SetLineColor(2);
  h_exp->SetLineWidth(3);
  //h_exp->Sumw2();
  h_exp->Scale(scale);

  h_pred->Draw(); 
  h_exp->Draw("same");

  const std::string titre="CMS Preliminary 2016, 2.3 fb^{-1}, #sqrt{s} = 13 TeV";
  TLatex *title = new TLatex(0.09770115,0.9194915,titre.c_str());
  title->SetNDC();
  title->SetTextSize(0.045);
  title->Draw("same");

  TLegend* leg = new TLegend(0.6,0.75,0.85,0.85);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetFillColor(0);
  leg->AddEntry(h_pred,"Prediction","l");
  leg->AddEntry(h_exp,"Expectation","l");
  leg->Draw("same");

  c->Update(); 
 
  pad->cd(2); 
  TPad *pad2 = (TPad*) pad->GetPad(2);
  pad2->SetPad("ratio", "", 0, 0, 1.0, divRatio, kWhite);
  pad2->SetBottomMargin(0.3);
  pad2->SetTopMargin(small);
  pad2->SetBorderMode(0);

  TH1D *ratio = (TH1D*) h_exp->Clone();
  TH1D *allmc = (TH1D*) h_pred->Clone();

  ratio->Add(allmc, -1);
  ratio->Divide(allmc);
  ratio->GetYaxis()->SetTitle( "(exp - pred)/pred" );

  TAxis* xHT = ratio->GetXaxis();

  xHT->SetTickLength(xHT->GetTickLength()*labelRatio);
  xHT->SetLabelSize(xHT->GetLabelSize()*labelRatio);
  xHT->SetLabelOffset(xHT->GetLabelOffset()*labelRatio);
  ratio->SetMinimum(-1.0);
  ratio->SetMaximum(1.0);

  TAxis* yHT = ratio->GetYaxis();
  yHT->SetNdivisions(010);
  yHT->SetLabelSize(yHT->GetLabelSize()*2.0);
  yHT->SetTitleOffset(0.3);
  yHT->SetTitleSize(0.08);

  ratio->SetTitleSize(0.15);
  ratio->SetStats(kFALSE);
  ratio->SetMarkerStyle(kFullDotMedium);
  //ratio->Sumw2();
  ratio->DrawCopy();

  TH1D *zero = (TH1D*)ratio->Clone(); 
  zero->Reset();
  for(int ib=0; ib<ratio->GetNbinsX(); ib++){ zero->SetBinContent(ib+1, 0); }
  zero->SetLineColor(kRed); zero->SetLineWidth(1);
  zero->DrawCopy("same");
  
  c->SaveAs( TString("Plotsc/") + hist_tag + TString(".png") );
  c->SaveAs( TString("Plotsc/") + hist_tag + TString(".C") );
  */
}


struct Plotting_Parameter
{
  TString hist_tag;
  TString XTitle;
  double min;
  double max;
};
