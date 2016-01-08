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
  //fin = TFile::Open("test.root");
  //fin = TFile::Open("v151201_0p5fb_ExpLL.root");
  //fin = TFile::Open("v151201_2p1fb_ExpLL.root");
  //fin = TFile::Open("v151216_ExpLL.root");
  fin = TFile::Open("v160106_ttbar_ExpLL.root");
  //fin2 = TFile::Open("v151204_PredLL.root");
  //fin2 = TFile::Open("v151209_PredLL.root");
  //fin2 = TFile::Open("v151216_PredLL.root");
  //fin2 = TFile::Open("ttbar_v151217_PredLL.root");
  //fin2 = TFile::Open("v160106_ttbar_PredLL.root");
  //fin2 = TFile::Open("v160106_data_PredLL.root");
  fin2 = TFile::Open("v160106_data_trigSel_PredLL.root");
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

  const std::string titre="CMS Preliminary 2016, 2.2 fb^{-1}, #sqrt{s} = 13 TeV";
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
  h_pred->Sumw2();
  h_pred->Scale(scale);

  h_exp->GetXaxis()->SetRangeUser(min,max);
  h_exp->GetXaxis()->SetTitle(XTitle);
  h_exp->SetLineColor(2);
  h_exp->SetLineWidth(3);
  h_exp->Sumw2();
  h_exp->Scale(scale);

  h_pred->Draw(); 
  h_exp->Draw("same");

  const std::string titre="CMS Preliminary 2016, 2.2 fb^{-1}, #sqrt{s} = 13 TeV";
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
  ratio->Sumw2();
  ratio->DrawCopy();

  TH1D *zero = (TH1D*)ratio->Clone(); 
  zero->Reset();
  for(int ib=0; ib<ratio->GetNbinsX(); ib++){ zero->SetBinContent(ib+1, 0); }
  zero->SetLineColor(kRed); zero->SetLineWidth(1);
  zero->DrawCopy("same");

  c->SaveAs( TString("Plotsc/") + hist_tag + TString(".png") );
  c->SaveAs( TString("Plotsc/") + hist_tag + TString(".C") );
}


struct Plotting_Parameter
{
  TString hist_tag;
  TString XTitle;
  double min;
  double max;
};
