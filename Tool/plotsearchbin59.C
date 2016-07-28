#include "HistStyle.h"
#include "HistReader.h"

void plotsearchbin59(const TString &fileName = "app_ExpLL.root",
		     const TString &graphicsFormat = "pdf") {
  gROOT->ProcessLine(".L HistStyle.h+");
  gROOT->ProcessLine(".L HistReader.h+");
  HistStyle::init();
  int canSizeX = 600;
  int canSizeY = 600;

  double ymax_Yields = 200000000.;
  double ymin_Yields = 0.0002;
  double ymax_Ratio = 3;
  double ymin_Ratio = 0;
  
  // Get histograms from file
  const unsigned int kNDists = 6;
  TH1* hPred[kNDists];
  TH1* hTrue[kNDists];
  for(unsigned int i = 0; i < kNDists; ++i) {
    TString name = "";
    if(      i == 0 ) name = "Yields";
    else if( i == 1 ) name = "met";
    else if( i == 2 ) name = "NJets";
    else if( i == 3) name = "NbJets";
    else if( i == 4) name = "NTops";
    else if( i == 5) name = "MT2";
    else if (i == 6) name = "Ht"; 
    else if (i == 7) name = "dPhi0"; 
    else if (i == 8) name = "dPhi1"; 
    else if (i == 9) name = "dPhi2"; 
    TString fileName2 = "app_PredLL.root";
   hPred[i] = HistReader::get(fileName2,"h_pred_lept_isotrk");
   hTrue[i] = HistReader::get(fileName,"h_exp_lept_isotrk");
  }
  
  // Set style
  for(unsigned int i = 0; i < kNDists; ++i) {
   hTrue[i]->SetTitle("");
   hTrue[i]->SetMarkerStyle(20);
   hTrue[i]->SetMarkerColor(kRed);
   hTrue[i]->SetLineColor(hTrue[i]->GetMarkerColor());
   hTrue[i]->SetMarkerSize(0.9);
   hTrue[i]->GetYaxis()->SetTitleOffset(0.6);
   hTrue[i]->GetYaxis()->SetTitleFont(42);
   hTrue[i]->GetYaxis()->SetLabelSize(0.065);
   hTrue[i]->GetYaxis()->SetLabelSize(0.04);
   hTrue[i]->GetYaxis()->SetLabelFont(42);
   hPred[i]->SetLineColor(kBlue);
   hPred[i]->SetMarkerStyle(20);
   hPred[i]->SetMarkerSize(0);
    if (i == 0){
     hTrue[i]->SetAxisRange(0,59,"X");
     hPred[i]->SetAxisRange(0,59,"X");
     hTrue[i]->SetAxisRange(ymin_Yields,ymax_Yields,"Y");      
     hPred[i]->SetAxisRange(ymin_Yields,ymax_Yields,"Y");      
    }
   }

  // Ratio plots
  TH1* hRatioFrame[kNDists];
  TH1* hRatio[kNDists];
  for(unsigned int i = 0; i < kNDists; ++i) {
    hRatio[i] = static_cast<TH1*>(hTrue[i]->Clone("Ratio"));
    hRatio[i]->Divide(hPred[i]);
    hRatio[i]->GetYaxis()->SetTitle("#frac{Direct}{Prediction}");
    hRatio[i]->GetYaxis()->SetRangeUser(0.2,1.9);
    hRatio[i]->SetTitle("");
    hRatio[i]->SetStats(0);
    hRatio[i]->SetLineWidth(1);
    hRatio[i]->GetYaxis()->SetTitleSize(0.15);
    hRatio[i]->GetYaxis()->SetTitleOffset(0.3);
    hRatio[i]->GetYaxis()->SetTitleFont(42);
    hRatio[i]->GetYaxis()->SetLabelSize(0.1);
    hRatio[i]->GetYaxis()->SetLabelFont(42);
    hRatio[i]->GetXaxis()->SetLabelOffset(0.01);
    hRatio[i]->GetXaxis()->SetLabelFont(42);    
    hRatio[i]->GetXaxis()->SetLabelSize(0.08);
    hRatio[i]->GetXaxis()->SetTitleSize(0.16);
    hRatio[i]->GetXaxis()->SetTitleFont(42);
    hRatio[i]->GetXaxis()->SetTitleOffset(0.6);
    if (i == 0){
      hRatio[i]->GetXaxis()->SetTitle("Search region bin number");
      hRatio[i]->SetAxisRange(0,59,"X");
      hRatio[i]->SetAxisRange(0,3,"Y");
    }
    
    hRatioFrame[i] = static_cast<TH1*>(hRatio[i]->Clone("RatioFrame"));
    for(int bin = 1; bin <= hRatioFrame[i]->GetNbinsX(); ++bin) {
      hRatioFrame[i]->SetBinContent(bin,1.);
    }
    hRatioFrame[i]->SetLineStyle(2);
    hRatioFrame[i]->SetLineWidth(2);
    hRatioFrame[i]->SetLineColor(kBlack);
  }

  // Create legend
  TLegend* leg = new TLegend(0.56,0.59,0.94,0.85);
  leg->SetBorderSize(1);
  leg->SetLineColor(1);
  leg->SetLineWidth(2);
  leg->SetFillColor(0);
  //  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader("Lost-lepton background");
  leg->AddEntry(hTrue[0],"Direct from simulation","P");
  leg->AddEntry(hPred[0],"Treat simulation like data");


  // Draw
  for(unsigned int i = 0; i < kNDists; ++i) {
    TString name = "";
    if(      i == 0 ) {name = "Yields";canSizeX=900.;}
    else if( i == 1 ) name = "MET";
    else if( i == 2 ) name = "NJets";
    else if( i == 3) name = "NbJets";
    else if( i == 4) name = "NTops";
    else if( i == 5) name = "MT2";
    else if( i == 6) name = "HT";
    else if( i == 7) name = "dPhi0";
    else if( i == 8) name = "dPhi1";
    else if( i == 9) name = "dPhi2";
    TCanvas* can = new TCanvas(name,name,canSizeX,canSizeY);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetTopMargin(0.12);   
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetLeftMargin(0.10);
    pad1->SetRightMargin(0.05);
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
        
    TExec * setex2 = new TExec("setex2", "gStyle->SetErrorX(0.0)");
    TExec *setex = new TExec("setex", "gStyle->SetErrorX(0.5)");
    hTrue[i]->GetXaxis()->SetRangeUser(0., 59.);   
    setex2->Draw();
    hTrue[i]->Draw("PE1");
    setex->Draw();
    hPred[i]->DrawCopy("histsame");
    hPred[i]->SetFillColor(kBlue-4);
    hPred[i]->SetFillStyle(3001);
    hPred[i]->Draw("E2Lsame");
    setex2->Draw();
    hTrue[i]->Draw("PE1same");    
       if (i==0){
      //-----------------------------------------------------------
      // Putting lines and labels explaining search region definitions
      //-----------------------------------------------------------
      // Ntop separation lines
      TLine *tl_ntop = new TLine();
      tl_ntop->SetLineStyle(2);
      tl_ntop->DrawLine(28,ymin_Yields,28,ymax_Yields);
      tl_ntop->DrawLine(53,ymin_Yields,53,ymax_Yields);

      // Ntop labels
      TLatex * ttext_ntop = new TLatex();
      ttext_ntop->SetTextFont(42);
      ttext_ntop->SetTextSize(0.07);
      ttext_ntop->SetTextAlign(22);
      ttext_ntop->DrawLatex(11.5, ymax_Yields/3 ,"N_{t} = 1");
      ttext_ntop->DrawLatex(39.5, ymax_Yields/40000. ,"N_{t} = 2");
      ttext_ntop->SetTextAngle(90.);
      ttext_ntop->DrawLatex(55.5 , ymax_Yields/200000. ,"N_{t} #geq 3");
      
      // Nb separation lines
      TLine *tl_nb = new TLine();
      tl_nb->SetLineStyle(3);
      //      tl_nb->SetLineColor(32);
      tl_nb->DrawLine(12,ymin_Yields,12,ymax_Yields/20.); 
      tl_nb->DrawLine(23,ymin_Yields,23,ymax_Yields/20.); 
      tl_nb->DrawLine(40,ymin_Yields,40,ymax_Yields/200000.); 
      tl_nb->DrawLine(49,ymin_Yields,49,ymax_Yields/200000.); 
      tl_nb->DrawLine(55,ymin_Yields,55,ymax_Yields/4000000.); 
      tl_nb->DrawLine(57,ymin_Yields,57,ymax_Yields/4000000.); 

      // Nb labels
      TLatex * ttext2 = new TLatex();
      ttext2->SetTextFont(42);
      // ttext2->SetTextColor(32);
      ttext2->SetTextSize(0.07);
      ttext2->SetTextAlign(22);
      ttext2->Draw();

      ttext2->DrawLatex( 4.5, ymax_Yields/20., "N_{b} = 1");
      ttext2->DrawLatex(17.5, ymax_Yields/20., "N_{b} = 2");
      // ttext2->DrawLatex( 28, ymax_Yields/200., "N_{b} = 1");
      // ttext2->DrawLatex(36, ymax_Yields/200., "N_{b} = 2");
      ttext2->SetTextAngle(90.);
      ttext2->DrawLatex(26, ymax_Yields/200. , "N_{b} #geq 3");

      // MT2 separation lines
      TLine *tl_mt2 = new TLine();
      tl_mt2->SetLineStyle(4);
      // tl_mt2->SetLineColor(49);
      tl_mt2->DrawLine(4,ymin_Yields,4,ymax_Yields/200.); 
      tl_mt2->DrawLine(8,ymin_Yields,8,ymax_Yields/200.); 
      tl_mt2->DrawLine(16,ymin_Yields,16,ymax_Yields/200.); 
      tl_mt2->DrawLine(20,ymin_Yields,20,ymax_Yields/200.); 

      tl_mt2->DrawLine(26,ymin_Yields,26,ymax_Yields/2000.); 
      tl_mt2->DrawLine(32,ymin_Yields,32,ymax_Yields/400000.); 
      tl_mt2->DrawLine(36,ymin_Yields,36,ymax_Yields/400000.); 
      tl_mt2->DrawLine(43,ymin_Yields,43,ymax_Yields/400000.); 
      tl_mt2->DrawLine(46,ymin_Yields,46,ymax_Yields/400000.); 
      tl_mt2->DrawLine(51,ymin_Yields,51,ymax_Yields/2000000.); 
     
      // MT2 labels
      TLatex * ttextmt2 = new TLatex();
      ttextmt2->SetTextFont(42);
      // ttextmt2->SetTextColor(49);
      ttextmt2->SetTextSize(0.045);
      ttextmt2->SetTextAlign(12);
      ttextmt2->SetTextAngle(90);
      ttextmt2->DrawLatex( 2.5, ymax_Yields/2000000. , "M_{T2}=[200,300]");
      ttextmt2->DrawLatex( 6.0, ymax_Yields/2000000. , "M_{T2}=[300,400]");
      ttextmt2->DrawLatex(9.5, ymax_Yields/2000000. , "M_{T2}#geq400 GeV");
      ttextmt2->DrawLatex( 27, ymax_Yields/10000000. , "M_{T2}#geq350 GeV");
      //-----------------------------------------------------------
 
    }
   
    leg->Draw("same");
    name = "hClosure"+name;
    gPad->SetLogy();
    TLatex mark;
    mark.SetNDC(true);
    double fontScale = 1.6;
    char lumistamp[128];
    //sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", 8000/ 1000.0);
    sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", 12900/ 1000.0);
    //Draw CMS mark                                                                                                                    
    mark.SetTextAlign(11);
    mark.SetTextSize(0.042 * fontScale * 1.25);
    mark.SetTextFont(61);
    mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS"); // #scale[0.8]{#it{Preliminary}}");        
    mark.SetTextSize(0.042 * fontScale);
    mark.SetTextFont(52);
    mark.DrawLatex(gPad->GetLeftMargin() + 0.09, 1 - (gPad->GetTopMargin() - 0.017), "Simulation");
    //Draw lumistamp                                                                                                                         
    mark.SetTextFont(42);
    mark.SetTextAlign(31);
    mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);
    
    can->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);
    pad2->SetLeftMargin(0.1);
    pad2->SetRightMargin(0.05);
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    hRatio[i]->SetMarkerSize(1);
    hRatio[i]->SetMaximum(2.9);
    if(i==0){
      hRatio[i]->GetXaxis()->SetLabelSize(0.15);
      hRatio[i]->GetXaxis()->SetTitleOffset(0.8);
    }
    hRatio[i]->GetXaxis()->SetRangeUser(0., 59.);
    //hRatio[i]->Draw("PE1");
    hRatio[i]->Draw("PE0");
    hRatioFrame[i]->Draw("HISTsame");
    hRatio[i]->Draw("PE1same");   
   
    if (i==0){

      //-----------------------------------------------------------
      // Putting lines and labels explaining search region definitions
      //-----------------------------------------------------------
      //Ntop separation lines
      TLine *tl_ntop = new TLine();
      tl_ntop->SetLineStyle(2);
      tl_ntop->DrawLine(28,ymin_Ratio,28,ymax_Ratio); 
      tl_ntop->DrawLine(53,ymin_Ratio,53,ymax_Ratio); 

      // Nb separation lines
      TLine *tl_nb = new TLine();
      tl_nb->SetLineStyle(3);
      tl_nb->SetLineColor(1);
      tl_nb->DrawLine(12,ymin_Ratio,12,ymax_Ratio); 
      tl_nb->DrawLine(23,ymin_Ratio,23,ymax_Ratio); 
      tl_nb->DrawLine(40,ymin_Ratio,40,ymax_Ratio); 
      tl_nb->DrawLine(49,ymin_Ratio,49,ymax_Ratio); 
      tl_nb->DrawLine(55,ymin_Ratio,55,ymax_Ratio); 
      tl_nb->DrawLine(57,ymin_Ratio,57,ymax_Ratio); 
      
      // MT2 separation lines
      TLine *tl_mt2 = new TLine();
      tl_mt2->SetLineStyle(4);
      tl_mt2->SetLineColor(1);
      tl_mt2->DrawLine(4,ymin_Ratio,4,ymax_Ratio); 
      tl_mt2->DrawLine(8,ymin_Ratio,8,ymax_Ratio); 
      tl_mt2->DrawLine(16,ymin_Ratio,16,ymax_Ratio); 
      tl_mt2->DrawLine(20,ymin_Ratio,20,ymax_Ratio); 

      tl_mt2->DrawLine(26,ymin_Ratio,26,ymax_Ratio); 
      tl_mt2->DrawLine(32,ymin_Ratio,32,ymax_Ratio); 
      tl_mt2->DrawLine(36,ymin_Ratio,36,ymax_Ratio); 
      tl_mt2->DrawLine(43,ymin_Ratio,43,ymax_Ratio); 
      tl_mt2->DrawLine(46,ymin_Ratio,46,ymax_Ratio); 
      tl_mt2->DrawLine(51,ymin_Ratio,51,ymax_Ratio); 

      //-----------------------------------------------------------
    }
   
    //    can->SaveAs(name+"_hadtau."+graphicsFormat);
   
  }

  //Some checking
  /*  std::cout << std::setprecision(3) << std::fixed;
  cout<<"Exp: "<<hTrue[0]->Integral()<<" Pred: "<<hPred[0]->Integral()<<endl;
  double ClosureUnc[59]={0};
  for(unsigned ib =0; ib<hTrue[0]->GetNbinsX(); ib++){
    //    cout<<ib<<"\t"<<hTrue[0]->GetBinContent(ib+1)<<" ± "<<hTrue[0]->GetBinError(ib+1)<<"\t"<<hPred[0]->GetBinContent(ib+1)<<" ± "<<hPred[0]->GetBinError(ib+1)<<"\t"<<hRatio[0]->GetBinContent(ib+1)<<" ± "<<hRatio[0]->GetBinError(ib+1)<<"\t"<<hTrue[0]->GetBinContent(ib+1)/hPred[0]->GetBinContent(ib+1)<<"\t"<<fabs(hTrue[0]->GetBinContent(ib+1)-hPred[0]->GetBinContent(ib+1))/hPred[0]->GetBinContent(ib+1)<<endl;
    double RatioError = hRatio[0]->GetBinError(ib+1);
    double RelnonClosure = std::fabs((hPred[0]->GetBinContent(ib+1)-hPred[0]->GetBinContent(ib+1))/hPred[0]->GetBinContent(ib+1));
    ClosureUnc[ib] = RelnonClosure > RatioError? RelnonClosure : RatioError;
    cout<<ClosureUnc[ib]<<", ";
    }*/
}
