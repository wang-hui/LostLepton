#include "CMS_lumi.h"

void drawSBregionDef(const double ymin_Yields = 0.05, const double ymax_Yields = 500., const bool logscale=true){
    int NSB = 37;

    const double adjHalfBin = 0.5;
    //-----------------------------------------------------------
    // Putting lines and labels explaining search region definitions
    //-----------------------------------------------------------
    // Ntop separation lines
    TLine *tl_ntop = new TLine();
    tl_ntop->SetLineStyle(2);
    if(NSB == 45)
    {
	tl_ntop->DrawLine(23.5 + adjHalfBin,ymin_Yields,23.5 + adjHalfBin,ymax_Yields);
	tl_ntop->DrawLine(41.5 + adjHalfBin,ymin_Yields,41.5 + adjHalfBin,ymax_Yields);
    } else if(NSB == 37)
    {
	tl_ntop->DrawLine(20.5 + adjHalfBin,ymin_Yields,20.5 + adjHalfBin,ymax_Yields);
    }
   // Ntop labels
   TLatex * ttext_ntop = new TLatex();
   ttext_ntop->SetTextFont(42);
   ttext_ntop->SetTextSize(0.06);
   ttext_ntop->SetTextAlign(22);
   if(logscale)
   {
       if(NSB == 45)
       {
	   ttext_ntop->DrawLatex(11.5 + adjHalfBin, ymax_Yields/1.8 ,"N_{top} = 1");
	   ttext_ntop->DrawLatex(35.5 + adjHalfBin, ymax_Yields/35. ,"N_{top} = 2");
	   ttext_ntop->SetTextAngle(90.);
	   ttext_ntop->DrawLatex(43. + adjHalfBin, ymax_Yields/15. ,"N_{top} #geq 3"); 
       } else if (NSB == 37)
       {
	   ttext_ntop->DrawLatex(10 + adjHalfBin, ymax_Yields/1.8 ,"N_{top} = 1");
	   ttext_ntop->DrawLatex(29.5 + adjHalfBin, ymax_Yields/35. ,"N_{top} = 2");
       }
   } else
   {
       if(NSB == 45)
       {
	   ttext_ntop->DrawLatex(11.5 + adjHalfBin, ymax_Yields*0.92 ,"N_{top} = 1");
	   ttext_ntop->DrawLatex(35.5 + adjHalfBin, ymax_Yields*0.92 ,"N_{top} = 2");
	   ttext_ntop->SetTextAngle(90.);
	   ttext_ntop->DrawLatex(43. + adjHalfBin, ymax_Yields*0.85 ,"N_{top} #geq 3");
       } else if(NSB == 37)
       {
	   ttext_ntop->DrawLatex(10 + adjHalfBin, ymax_Yields*0.92 ,"N_{top} = 1");
	   ttext_ntop->DrawLatex(29.5 + adjHalfBin, ymax_Yields*0.92 ,"N_{top} = 2");
       }
   }

   // Nb separation lines
   TLine *tl_nb = new TLine();
   tl_nb->SetLineStyle(3);
   tl_nb->SetLineColor(32);
   if(logscale)
   {
       if(NSB == 45)
       {
	   tl_nb->DrawLine(10.5 + adjHalfBin,ymin_Yields,10.5 + adjHalfBin,ymax_Yields/4.);
	   tl_nb->DrawLine(20.5 + adjHalfBin,ymin_Yields,20.5 + adjHalfBin,ymax_Yields/4.);
	   tl_nb->DrawLine(31.5 + adjHalfBin,ymin_Yields,31.5 + adjHalfBin,ymax_Yields/40.);
	   tl_nb->DrawLine(39.5 + adjHalfBin,ymin_Yields,39.5 + adjHalfBin,ymax_Yields/40.);
	   tl_nb->DrawLine(42.5 + adjHalfBin,ymin_Yields,42.5 + adjHalfBin,ymax_Yields/40.);
	   tl_nb->DrawLine(43.5 + adjHalfBin,ymin_Yields,43.5 + adjHalfBin,ymax_Yields/40.);
       } else if (NSB == 37)
       {
	   tl_nb->DrawLine(10.5 + adjHalfBin,ymin_Yields,10.5 + adjHalfBin,ymax_Yields/4.);
	   tl_nb->DrawLine(28.5 + adjHalfBin,ymin_Yields,28.5 + adjHalfBin,ymax_Yields/40.);
       }
   } else
   {
       if(NSB == 45)
       {
	   tl_nb->DrawLine(10.5 + adjHalfBin,ymin_Yields,10.5 + adjHalfBin,ymax_Yields*0.8);
	   tl_nb->DrawLine(20.5 + adjHalfBin,ymin_Yields,20.5 + adjHalfBin,ymax_Yields*0.8);
	   tl_nb->DrawLine(31.5 + adjHalfBin,ymin_Yields,31.5 + adjHalfBin,ymax_Yields*0.8);
	   tl_nb->DrawLine(39.5 + adjHalfBin,ymin_Yields,39.5 + adjHalfBin,ymax_Yields*0.8);
	   tl_nb->DrawLine(42.5 + adjHalfBin,ymin_Yields,42.5 + adjHalfBin,ymax_Yields*0.8);
	   tl_nb->DrawLine(43.5 + adjHalfBin,ymin_Yields,43.5 + adjHalfBin,ymax_Yields*0.8);
       } else if (NSB == 37)
       {
	   tl_nb->DrawLine(10.5 + adjHalfBin,ymin_Yields,10.5 + adjHalfBin,ymax_Yields*0.8);
	   tl_nb->DrawLine(28.5 + adjHalfBin,ymin_Yields,28.5 + adjHalfBin,ymax_Yields*0.8);
       }
   }
   // Nb labels
   TLatex * ttext2 = new TLatex();
   ttext2->SetTextFont(42);
   ttext2->SetTextColor(32);
   ttext2->SetTextSize(0.05);
   ttext2->SetTextAlign(22);
   ttext2->Draw();

   if(logscale)
   {
       ttext2->DrawLatex( 4.5 + adjHalfBin, ymax_Yields/3., "N_{b} = 1");
       ttext2->DrawLatex(14.5 + adjHalfBin, ymax_Yields/3., "N_{b} = 2");
       if(NSB == 45)
       {
	   ttext2->SetTextAngle(90.);
	   ttext2->DrawLatex(21.5 + adjHalfBin, ymax_Yields/3. , "N_{b} #geq 3");
       }
   } else
   {
       ttext2->DrawLatex( 4.5 + adjHalfBin, ymax_Yields*0.8, "N_{b} = 1");
       ttext2->DrawLatex(14.5 + adjHalfBin, ymax_Yields*0.8, "N_{b} = 2");
       if(NSB == 45)
       {
	   ttext2->SetTextAngle(90.);
	   ttext2->DrawLatex(21.5 + adjHalfBin + 0.5, ymax_Yields*0.8, "N_{b} #geq 3");
       }
   }
   // MT2 separation lines
   TLine *tl_mt2 = new TLine();
   tl_mt2->SetLineStyle(4);
   tl_mt2->SetLineColor(49);
   if(logscale)
   {
       tl_mt2->DrawLine(3.5 + adjHalfBin,ymin_Yields,3.5 + adjHalfBin,ymax_Yields/20.);
       tl_mt2->DrawLine(7.5 + adjHalfBin,ymin_Yields,7.5 + adjHalfBin,ymax_Yields/20.);
       tl_mt2->DrawLine(14.5 + adjHalfBin,ymin_Yields,14.5 + adjHalfBin,ymax_Yields/20.);
       tl_mt2->DrawLine(18.5 + adjHalfBin,ymin_Yields,18.5 + adjHalfBin,ymax_Yields/20.);
       if(NSB == 45)
       {
	   tl_mt2->DrawLine(26.5 + adjHalfBin,ymin_Yields,26.5 + adjHalfBin,ymax_Yields/80.);
	   tl_mt2->DrawLine(29.5 + adjHalfBin,ymin_Yields,29.5 + adjHalfBin,ymax_Yields/80.);
	   tl_mt2->DrawLine(34.5 + adjHalfBin,ymin_Yields,34.5 + adjHalfBin,ymax_Yields/80.);
	   tl_mt2->DrawLine(37.5 + adjHalfBin,ymin_Yields,37.5 + adjHalfBin,ymax_Yields/80.);
       } else if(NSB == 37)
       {
	   tl_mt2->DrawLine(23.5 + adjHalfBin,ymin_Yields,23.5 + adjHalfBin,ymax_Yields/80.);
	   tl_mt2->DrawLine(26.5 + adjHalfBin,ymin_Yields,26.5 + adjHalfBin,ymax_Yields/80.);
	   tl_mt2->DrawLine(31.5 + adjHalfBin,ymin_Yields,31.5 + adjHalfBin,ymax_Yields/80.);
	   tl_mt2->DrawLine(34.5 + adjHalfBin,ymin_Yields,34.5 + adjHalfBin,ymax_Yields/80.);
       }
   } else
   {
       tl_mt2->DrawLine(3.5 + adjHalfBin,ymin_Yields,3.5 + adjHalfBin,ymax_Yields*0.6);
       tl_mt2->DrawLine(7.5 + adjHalfBin,ymin_Yields,7.5 + adjHalfBin,ymax_Yields*0.6);
       tl_mt2->DrawLine(14.5 + adjHalfBin,ymin_Yields,14.5 + adjHalfBin,ymax_Yields*0.6);
       tl_mt2->DrawLine(18.5 + adjHalfBin,ymin_Yields,18.5 + adjHalfBin,ymax_Yields*0.6);
       if(NSB == 45)
       {       
	   tl_mt2->DrawLine(26.5 + adjHalfBin,ymin_Yields,26.5 + adjHalfBin,ymax_Yields*0.6);
	   tl_mt2->DrawLine(29.5 + adjHalfBin,ymin_Yields,29.5 + adjHalfBin,ymax_Yields*0.6);
	   tl_mt2->DrawLine(34.5 + adjHalfBin,ymin_Yields,34.5 + adjHalfBin,ymax_Yields*0.6);
	   tl_mt2->DrawLine(37.5 + adjHalfBin,ymin_Yields,37.5 + adjHalfBin,ymax_Yields*0.6);
       } else if(NSB == 37)
       {
	   tl_mt2->DrawLine(23.5 + adjHalfBin,ymin_Yields,23.5 + adjHalfBin,ymax_Yields*0.6);
	   tl_mt2->DrawLine(26.5 + adjHalfBin,ymin_Yields,26.5 + adjHalfBin,ymax_Yields*0.6);
	   tl_mt2->DrawLine(31.5 + adjHalfBin,ymin_Yields,31.5 + adjHalfBin,ymax_Yields*0.6);
	   tl_mt2->DrawLine(34.5 + adjHalfBin,ymin_Yields,34.5 + adjHalfBin,ymax_Yields*0.6);
       }
   }
   // MT2 labels
   TLatex * ttextmt2 = new TLatex();
   ttextmt2->SetTextFont(42);
   ttextmt2->SetTextColor(49);
   ttextmt2->SetTextSize(0.02);
   ttextmt2->SetTextAlign(22);
   ttextmt2->SetTextAngle(90);
   if(logscale)
   {
       ttextmt2->DrawLatex( 1.5 + adjHalfBin, ymax_Yields/10. , "M_{T2}=[200,300]");
       ttextmt2->DrawLatex( 5.0 + adjHalfBin, ymax_Yields/10. , "M_{T2}=[300,400]");
       ttextmt2->DrawLatex(8.5 + adjHalfBin, ymax_Yields/10. , "M_{T2}#geq400");
   } else
   {
       ttextmt2->DrawLatex( 1.5 + adjHalfBin, ymax_Yields*0.55 , "M_{T2}=[200,300]");
       ttextmt2->DrawLatex( 5.0 + adjHalfBin, ymax_Yields*0.55 , "M_{T2}=[300,400]");
       ttextmt2->DrawLatex(8.5 + adjHalfBin, ymax_Yields*0.55 , "M_{T2}#geq400");
   }
   //-----------------------------------------------------------
}


void Plot_pred_v6()
{

  double Pred[37];

Pred[0] = 14.966986206;
Pred[1] = 2.8397569373;
Pred[2] = 1.7515799273;
Pred[3] = 0;
Pred[4] = 3.0238918586;
Pred[5] = 5.7521929994;
Pred[6] = 1.6478511195;
Pred[7] = 0;
Pred[8] = 0.90563180631;
Pred[9] = 2.1528947325;
Pred[10] = 2.331821985;
Pred[11] = 16.254900066;
Pred[12] = 5.7242718946;
Pred[13] = 0.39153569034;
Pred[14] = 0.47626785819;
Pred[15] = 2.062977317;
Pred[16] = 1.3628468412;
Pred[17] = 0.5377642659;
Pred[18] = 0;
Pred[19] = 0.48624384132;
Pred[20] = 0;
Pred[21] = 5.4188607609;
Pred[22] = 1.6879924485;
Pred[23] = 0.32105539163;
Pred[24] = 2.8844669648;
Pred[25] = 2.1005053514;
Pred[26] = 1.0788833013;
Pred[27] = 0;
Pred[28] = 0.45345327562;
Pred[29] = 7.6172217543;
Pred[30] = 1.5792066031;
Pred[31] = 0.18283642101;
Pred[32] = 3.5135658247;
Pred[33] = 0.59907219402;
Pred[34] = 0;
Pred[35] = 0;
Pred[36] = 0;

  double uncUpPred[37];
  double uncDownPred[37];

uncUpPred[0] = 2.7967550481;
uncDownPred[0] = 2.6350282405;
uncUpPred[1] = 1.5679315516;
uncDownPred[1] = 1.2056304163;
uncUpPred[2] = 1.4786544011;
uncDownPred[2] = 1.0678053218;
uncUpPred[3] = 0.87514098028;
uncDownPred[3] = 0;
uncUpPred[4] = 1.4666661274;
uncDownPred[4] = 1.1803413142;
uncUpPred[5] = 2.0723119241;
uncDownPred[5] = 1.8023373233;
uncUpPred[6] = 1.5126382404;
uncDownPred[6] = 0.96874970321;
uncUpPred[7] = 0.91438010697;
uncDownPred[7] = 0;
uncUpPred[8] = 0.99002726007;
uncDownPred[8] = 0.5239046075;
uncUpPred[9] = 1.853340723;
uncDownPred[9] = 1.2593377383;
uncUpPred[10] = 1.9765446398;
uncDownPred[10] = 1.3696220033;
uncUpPred[11] = 3.1706899669;
uncDownPred[11] = 3.0237047761;
uncUpPred[12] = 2.0951654142;
uncDownPred[12] = 1.7955762707;
uncUpPred[13] = 1.1974652474;
uncDownPred[13] = 0.39153569034;
uncUpPred[14] = 1.2980065582;
uncDownPred[14] = 0.47626785819;
uncUpPred[15] = 1.2480925093;
uncDownPred[15] = 0.84729315349;
uncUpPred[16] = 1.3932840225;
uncDownPred[16] = 0.78907764493;
uncUpPred[17] = 1.365775925;
uncDownPred[17] = 0.5377642659;
uncUpPred[18] = 1.1345737166;
uncDownPred[18] = 0;
uncUpPred[19] = 1.2209443722;
uncDownPred[19] = 0.48624384132;
uncUpPred[20] = 1.6786440997;
uncDownPred[20] = 0;
uncUpPred[21] = 1.6800115552;
uncDownPred[21] = 1.5864322181;
uncUpPred[22] = 0.94129132422;
uncDownPred[22] = 0.74460084481;
uncUpPred[23] = 0.5900581438;
uncDownPred[23] = 0.32105539163;
uncUpPred[24] = 1.5121431929;
uncDownPred[24] = 1.2571422268;
uncUpPred[25] = 1.3158400892;
uncDownPred[25] = 0.99148872881;
uncUpPred[26] = 0.91890012085;
uncDownPred[26] = 0.65224265068;
uncUpPred[27] = 0.92734282761;
uncDownPred[27] = 0;
uncUpPred[28] = 0.92883834591;
uncDownPred[28] = 0.45345327562;
uncUpPred[29] = 2.5023926895;
uncDownPred[29] = 2.43051492;
uncUpPred[30] = 1.0290955668;
uncDownPred[30] = 0.84968724761;
uncUpPred[31] = 0.58181653886;
uncDownPred[31] = 0.18283642101;
uncUpPred[32] = 1.4094456072;
uncDownPred[32] = 1.1367411199;
uncUpPred[33] = 1.1071648503;
uncDownPred[33] = 0.59907219402;
uncUpPred[34] = 0.76583341382;
uncDownPred[34] = 0;
uncUpPred[35] = 0.81898299273;
uncDownPred[35] = 0;
uncUpPred[36] = 1.1367889905;
uncDownPred[36] = 0;


  const std::string titre="CMS Supplementary";

  TCanvas *c1 = new TCanvas("c1", "c1",0,51,1920,1004);
  c1->SetFillColor(0);
  c1->cd();

//  TH1D *h_mu_iso_eff = new TH1D("h_mu_iso_eff", "h_mu_iso_eff", 37 , 0 , 37);
//  h_mu_iso_eff->SetTitle("");
//  h_mu_iso_eff->SetXTitle("Search Bins");
//  //h_mu_iso_eff->SetYTitle("iso track veto efficiency");
//  h_mu_iso_eff->SetStats(0);
//  h_mu_iso_eff->SetLineWidth(3);
//
//  for (Int_t sbc=1;sbc<=37;++sbc) {
//    h_mu_iso_eff->SetBinContent(sbc,Pred[sbc-1]);
//    h_mu_iso_eff->SetBinError(sbc,isoTrackErr[sbc-1]);
//  }
//
//  gStyle->SetPaintTextFormat("1.2f");
//  h_mu_iso_eff->Draw();

  Double_t x[37];
  Double_t exl[37];
  Double_t exh[37];
  for (int xc=0;xc<37;++xc)
  {
    x[xc]=xc+0.5;
    exl[xc]=0.5;
    exh[xc]=0.5;
  }

   gr = new TGraphAsymmErrors(37,x,Pred,exl,exh,uncDownPred,uncUpPred);
   gr->SetTitle("");
   gr->SetLineWidth(3);
   gr->SetLineColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");
   gr->SetMaximum(22.0);
   gr->GetXaxis()->SetTitle("Search region bin number");

  drawSBregionDef(0.0, 22.0, false);

//  TLatex *title = new TLatex(0.09770115,0.9194915,titre.c_str());
//  title->SetNDC();
//  title->SetTextSize(0.045);
//  title->Draw("same");
//
//   TLatex *   tex2 = new TLatex(0.7405858,0.9178645,"2.3 fb^{-1} (13 TeV)");
//tex2->SetNDC();
//   tex2->SetTextSize(0.045);
//   tex2->SetLineWidth(2);
//   tex2->Draw();

  CMSStylePlot::CMS_lumi( c1, 4, 0 );


  c1->SaveAs( "h_isotrackvetoEff.png" );

}
