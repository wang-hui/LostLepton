#include "CMS_lumi.h"

// Function to draw the signal bin definition
//void drawSBregionDef(const double ymin_Yields, const double ymax_Yields, const bool logscale){
    int NSB = 84;

    const double adjHalfBin = 0.5;
    const double deltaY = ymax_Yields - ymin_Yields;
    //-----------------------------------------------------------
    // Putting lines and labels explaining search region definitions
    //-----------------------------------------------------------
    // Ntop separation lines
    TLine *tl_ntop = new TLine();
    tl_ntop->SetLineStyle(2);
    if(NSB == 59)
      {
	tl_ntop->DrawLine(27.5 + adjHalfBin,ymin_Yields,27.5 + adjHalfBin,ymax_Yields);
        tl_ntop->DrawLine(52.5 + adjHalfBin,ymin_Yields,52.5 + adjHalfBin,ymax_Yields);

      }
    else if(NSB == 45)
      {
        tl_ntop->DrawLine(23.5 + adjHalfBin,ymin_Yields,23.5 + adjHalfBin,ymax_Yields);
        tl_ntop->DrawLine(41.5 + adjHalfBin,ymin_Yields,41.5 + adjHalfBin,ymax_Yields);
      } 
    else if(NSB == 37)
      {
        tl_ntop->DrawLine(20.5 + adjHalfBin,ymin_Yields,20.5 + adjHalfBin,ymax_Yields);
      }
    // Ntop labels
    TLatex * ttext_ntop = new TLatex();
    ttext_ntop->SetTextFont(42);
    ttext_ntop->SetTextSize(0.045);
    ttext_ntop->SetTextAlign(22);
    if(logscale)
      {
      
	if(NSB == 59)
	{
	  ttext_ntop->DrawLatex(11.5 + adjHalfBin, ymax_Yields/2 ,"N_{t} = 1");
	  ttext_ntop->DrawLatex(39.5 + adjHalfBin, ymax_Yields/200. ,"N_{t} = 2");
	  ttext_ntop->SetTextAngle(90.);
	  ttext_ntop->DrawLatex(54.5 + adjHalfBin, ymax_Yields/20. ,"N_{t} #geq 3");

	}

        else if(NSB == 45)
        {
            ttext_ntop->DrawLatex(11.5 + adjHalfBin, ymax_Yields/3 ,"N_{t} = 1");
            ttext_ntop->DrawLatex(32.5 + adjHalfBin, ymax_Yields/300. ,"N_{t} = 2");
            ttext_ntop->SetTextAngle(90.);
            ttext_ntop->DrawLatex(42.5 + adjHalfBin, ymax_Yields/30. ,"N_{t} #geq 3"); 
        } 
        else if (NSB == 37)
        {
            ttext_ntop->DrawLatex(10 + adjHalfBin, ymax_Yields/1.8 ,"N_{t} = 1");
            ttext_ntop->DrawLatex(29.5 + adjHalfBin, ymax_Yields/70. ,"N_{t} #geq 2");
        }
    }
    else
      {
	if(NSB == 59)
	  {
	    ttext_ntop->DrawLatex(11.5 + adjHalfBin, ymax_Yields*0.93 ,"N_{t} = 1");
            ttext_ntop->DrawLatex(38.5 + adjHalfBin, ymax_Yields*0.93 ,"N_{t} = 2");
            ttext_ntop->SetTextAngle(90.);
	    ttext_ntop->DrawLatex(55.5 + adjHalfBin, ymax_Yields*0.85 ,"N_{t} #geq 3");
	  }
        else if(NSB == 45)
	  {
            ttext_ntop->DrawLatex(11.5 + adjHalfBin, ymax_Yields*0.93 ,"N_{t} = 1");
            ttext_ntop->DrawLatex(35.5 + adjHalfBin, ymax_Yields*0.93 ,"N_{t} = 2");
            ttext_ntop->SetTextAngle(90.);
            ttext_ntop->DrawLatex(43. + adjHalfBin, ymax_Yields*0.85 ,"N_{t} #geq 3");
	  }
        else if(NSB == 37)
	  {
            ttext_ntop->DrawLatex(10 + adjHalfBin, ymax_Yields*0.93 ,"N_{t} = 1");
            ttext_ntop->DrawLatex(29.5 + adjHalfBin, ymax_Yields*0.93 ,"N_{t} #geq 2");
	  }
    }

    // Nb separation lines
    TLine *tl_nb = new TLine();
    tl_nb->SetLineStyle(3);
    tl_nb->SetLineColor(32);
    if(logscale)
      {
	if(NSB == 59)
	  {
	    tl_nb->DrawLine(11.5 + adjHalfBin,ymin_Yields,11.5 + adjHalfBin,ymax_Yields/6);
            tl_nb->DrawLine(22.5 + adjHalfBin,ymin_Yields,22.5 + adjHalfBin,ymax_Yields/6.);
            tl_nb->DrawLine(39.5 + adjHalfBin,ymin_Yields,39.5 + adjHalfBin,ymax_Yields/160.);
            tl_nb->DrawLine(48.5 + adjHalfBin,ymin_Yields,48.5 + adjHalfBin,ymax_Yields/160.);
            tl_nb->DrawLine(54.5 + adjHalfBin,ymin_Yields,54.5 + adjHalfBin,ymax_Yields/160.);
            tl_nb->DrawLine(56.5 + adjHalfBin,ymin_Yields,56.5 + adjHalfBin,ymax_Yields/160.);  
	    
	  }
      
        else if(NSB == 45)
	  {
            tl_nb->DrawLine(10.5 + adjHalfBin,ymin_Yields,10.5 + adjHalfBin,ymax_Yields/40.);
            tl_nb->DrawLine(20.5 + adjHalfBin,ymin_Yields,20.5 + adjHalfBin,ymax_Yields/40.);
            tl_nb->DrawLine(31.5 + adjHalfBin,ymin_Yields,31.5 + adjHalfBin,ymax_Yields/400.);
            tl_nb->DrawLine(39.5 + adjHalfBin,ymin_Yields,39.5 + adjHalfBin,ymax_Yields/400.);
            tl_nb->DrawLine(42.5 + adjHalfBin,ymin_Yields,42.5 + adjHalfBin,ymax_Yields/400.);
            tl_nb->DrawLine(43.5 + adjHalfBin,ymin_Yields,43.5 + adjHalfBin,ymax_Yields/400.);
	  }
        else if (NSB == 37)
	  {
            tl_nb->DrawLine(10.5 + adjHalfBin,ymin_Yields,10.5 + adjHalfBin,ymax_Yields/4.);
            tl_nb->DrawLine(28.5 + adjHalfBin,ymin_Yields,28.5 + adjHalfBin,ymax_Yields/160.);
	  }
    }
    else
      {
	if(NSB == 59)
	  {
	    tl_nb->DrawLine(11.5 + adjHalfBin,ymin_Yields,11.5 + adjHalfBin,ymax_Yields*0.8);
            tl_nb->DrawLine(22.5 + adjHalfBin,ymin_Yields,22.5 + adjHalfBin,ymax_Yields*0.8);
            tl_nb->DrawLine(39.5 + adjHalfBin,ymin_Yields,39.5 + adjHalfBin,ymax_Yields*0.8);
            tl_nb->DrawLine(48.5 + adjHalfBin,ymin_Yields,48.5 + adjHalfBin,ymax_Yields*0.8);
            tl_nb->DrawLine(54.5 + adjHalfBin,ymin_Yields,54.5 + adjHalfBin,ymax_Yields*0.8);
            tl_nb->DrawLine(56.5 + adjHalfBin,ymin_Yields,56.5 + adjHalfBin,ymax_Yields*0.8);
	  }

        else if(NSB == 45)
	  {
            tl_nb->DrawLine(10.5 + adjHalfBin,ymin_Yields,10.5 + adjHalfBin,ymax_Yields*0.8);
            tl_nb->DrawLine(20.5 + adjHalfBin,ymin_Yields,20.5 + adjHalfBin,ymax_Yields*0.8);
            tl_nb->DrawLine(31.5 + adjHalfBin,ymin_Yields,31.5 + adjHalfBin,ymax_Yields*0.8);
            tl_nb->DrawLine(39.5 + adjHalfBin,ymin_Yields,39.5 + adjHalfBin,ymax_Yields*0.8);
            tl_nb->DrawLine(42.5 + adjHalfBin,ymin_Yields,42.5 + adjHalfBin,ymax_Yields*0.8);
            tl_nb->DrawLine(43.5 + adjHalfBin,ymin_Yields,43.5 + adjHalfBin,ymax_Yields*0.8);
	  } 
        else if (NSB == 37)
	  {
	    tl_nb->DrawLine(10.5 + adjHalfBin,ymin_Yields,10.5 + adjHalfBin,ymax_Yields*0.8);
            tl_nb->DrawLine(28.5 + adjHalfBin,ymin_Yields,28.5 + adjHalfBin,ymax_Yields*0.8);
	  }
    }
    // Nb labels
    TLatex * ttext2 = new TLatex();
    ttext2->SetTextFont(42);
    //ttext2->SetTextColor(32);
    ttext2->SetTextColor(kBlack);
    ttext2->SetTextSize(0.045);
    ttext2->SetTextAlign(22);
    ttext2->Draw();

    if(logscale)
      {
        ttext2->DrawLatex( 4.5 + adjHalfBin, ymax_Yields/4., "N_{b} = 1");
        if(NSB == 37) ttext2->DrawLatex(14.5 + adjHalfBin, ymax_Yields/3.5, "N_{b} #geq 2");
        else if(NSB == 45)
	  {
            ttext2->DrawLatex(14.5 + adjHalfBin, ymax_Yields/7, "N_{b} = 2");
            ttext2->SetTextAngle(90.);
            ttext2->DrawLatex(22.5  + adjHalfBin, ymax_Yields/20 , "N_{b} #geq 3");
	  }
	else if(NSB == 59)
	  {
	    ttext2->DrawLatex(17.5 + adjHalfBin, ymax_Yields/4., "N_{b} = 2");
            ttext2->SetTextAngle(90.);
            ttext2->DrawLatex(23.5  + adjHalfBin, ymax_Yields/30 , "N_{b} #geq 3");
	  }
      } 
    else
      {
        ttext2->DrawLatex( 4.5 + adjHalfBin, ymax_Yields*0.87, "N_{b} = 1");
        if(NSB == 37) ttext2->DrawLatex(14.5 + adjHalfBin, ymax_Yields*0.87, "N_{b} #geq 2");
        else if(NSB == 45)
	  {
            ttext2->DrawLatex(14.5 + adjHalfBin, ymax_Yields*0.87, "N_{b} = 2");
            ttext2->SetTextAngle(90.);
            ttext2->DrawLatex(21.5 + adjHalfBin + 0.5, ymax_Yields*0.87, "N_{b} #geq 3");
	  }
	else if(NSB == 59)
          {
            ttext2->DrawLatex(17.5 + adjHalfBin, ymax_Yields*0.87, "N_{b} = 2");
            ttext2->SetTextAngle(90.);
            ttext2->DrawLatex(25.5 + adjHalfBin + 0.5, ymax_Yields*0.87, "N_{b} #geq 3");
          }

      }
    // MT2 separation lines
    TLine *tl_mt2 = new TLine();
    tl_mt2->SetLineStyle(4);
    tl_mt2->SetLineColor(49);
    if(logscale)
    {
        tl_mt2->DrawLine(3.5 + adjHalfBin,ymin_Yields,3.5 + adjHalfBin,ymax_Yields/40.);
        tl_mt2->DrawLine(7.5 + adjHalfBin,ymin_Yields,7.5 + adjHalfBin,ymax_Yields/40.);
        

        if(NSB == 45)
	  {
	    tl_mt2->DrawLine(14.5 + adjHalfBin,ymin_Yields,14.5 + adjHalfBin,ymax_Yields/300.);
	    tl_mt2->DrawLine(18.5 + adjHalfBin,ymin_Yields,18.5 + adjHalfBin,ymax_Yields/300.);
            tl_mt2->DrawLine(26.5 + adjHalfBin,ymin_Yields,26.5 + adjHalfBin,ymax_Yields/300.);
            tl_mt2->DrawLine(29.5 + adjHalfBin,ymin_Yields,29.5 + adjHalfBin,ymax_Yields/300.);
            tl_mt2->DrawLine(34.5 + adjHalfBin,ymin_Yields,34.5 + adjHalfBin,ymax_Yields/300.);
            tl_mt2->DrawLine(37.5 + adjHalfBin,ymin_Yields,37.5 + adjHalfBin,ymax_Yields/300.);
        } 
        else if(NSB == 37)
	  {
	    tl_mt2->DrawLine(14.5 + adjHalfBin,ymin_Yields,14.5 + adjHalfBin,ymax_Yields/500.);
	    tl_mt2->DrawLine(18.5 + adjHalfBin,ymin_Yields,18.5 + adjHalfBin,ymax_Yields/500.);
	    tl_mt2->DrawLine(23.5 + adjHalfBin,ymin_Yields,23.5 + adjHalfBin,ymax_Yields/320.);
            tl_mt2->DrawLine(26.5 + adjHalfBin,ymin_Yields,26.5 + adjHalfBin,ymax_Yields/320.);
            tl_mt2->DrawLine(31.5 + adjHalfBin,ymin_Yields,31.5 + adjHalfBin,ymax_Yields/320.);
            tl_mt2->DrawLine(34.5 + adjHalfBin,ymin_Yields,34.5 + adjHalfBin,ymax_Yields/320.);
        }

	else if(NSB == 59)
	  {
	    tl_mt2->DrawLine(15.5 + adjHalfBin,ymin_Yields,15.5 + adjHalfBin,ymax_Yields/40.);
	    tl_mt2->DrawLine(19.5 + adjHalfBin,ymin_Yields,19.5 + adjHalfBin,ymax_Yields/40.);
	    tl_mt2->DrawLine(25.5 + adjHalfBin,ymin_Yields,25.5 + adjHalfBin,ymax_Yields/100.);
	    tl_mt2->DrawLine(31.5 + adjHalfBin,ymin_Yields,31.5 + adjHalfBin,ymax_Yields/100.);
            tl_mt2->DrawLine(35.5 + adjHalfBin,ymin_Yields,35.5 + adjHalfBin,ymax_Yields/320.);
            tl_mt2->DrawLine(42.5 + adjHalfBin,ymin_Yields,42.5 + adjHalfBin,ymax_Yields/320.);
            tl_mt2->DrawLine(45.5 + adjHalfBin,ymin_Yields,45.5 + adjHalfBin,ymax_Yields/320.);
            tl_mt2->DrawLine(50.5 + adjHalfBin,ymin_Yields,50.5 + adjHalfBin,ymax_Yields/320.);
	  }

    } 
    else
    {
        tl_mt2->DrawLine( 3.5 + adjHalfBin, ymin_Yields,  3.5 + adjHalfBin, ymax_Yields*0.6);
        tl_mt2->DrawLine( 7.5 + adjHalfBin, ymin_Yields,  7.5 + adjHalfBin, ymax_Yields*0.6);

        if(NSB == 45)
	  {
	    tl_mt2->DrawLine(14.5 + adjHalfBin, ymin_Yields, 14.5 + adjHalfBin, ymax_Yields*0.6);
	    tl_mt2->DrawLine(18.5 + adjHalfBin, ymin_Yields, 18.5 + adjHalfBin, ymax_Yields*0.6);  
            tl_mt2->DrawLine(26.5 + adjHalfBin, ymin_Yields, 26.5 + adjHalfBin, ymax_Yields*0.6);
            tl_mt2->DrawLine(29.5 + adjHalfBin, ymin_Yields, 29.5 + adjHalfBin, ymax_Yields*0.6);
            tl_mt2->DrawLine(34.5 + adjHalfBin, ymin_Yields, 34.5 + adjHalfBin, ymax_Yields*0.6);
            tl_mt2->DrawLine(37.5 + adjHalfBin, ymin_Yields, 37.5 + adjHalfBin, ymax_Yields*0.6);
	  } 
        else if(NSB == 37)
	  {
	    tl_mt2->DrawLine(14.5 + adjHalfBin, ymin_Yields, 14.5 + adjHalfBin, ymax_Yields*0.6);
	    tl_mt2->DrawLine(18.5 + adjHalfBin, ymin_Yields, 18.5 + adjHalfBin, ymax_Yields*0.6);
            tl_mt2->DrawLine(23.5 + adjHalfBin, ymin_Yields, 23.5 + adjHalfBin, ymax_Yields*0.6);
            tl_mt2->DrawLine(26.5 + adjHalfBin, ymin_Yields, 26.5 + adjHalfBin, ymax_Yields*0.6);
            tl_mt2->DrawLine(31.5 + adjHalfBin, ymin_Yields, 31.5 + adjHalfBin, ymax_Yields*0.6);
            tl_mt2->DrawLine(34.5 + adjHalfBin, ymin_Yields, 34.5 + adjHalfBin, ymax_Yields*0.6);
	  }
	else if(NSB == 59)
          {
            tl_mt2->DrawLine(15.5 + adjHalfBin, ymin_Yields, 15.5 + adjHalfBin, ymax_Yields*0.6);
	    tl_mt2->DrawLine(19.5 + adjHalfBin, ymin_Yields, 19.5 + adjHalfBin, ymax_Yields*0.6);
	    tl_mt2->DrawLine(25.5 + adjHalfBin, ymin_Yields, 25.5 + adjHalfBin, ymax_Yields*0.6);
	    tl_mt2->DrawLine(31.5 + adjHalfBin, ymin_Yields, 31.5 + adjHalfBin, ymax_Yields*0.6);
            tl_mt2->DrawLine(35.5 + adjHalfBin, ymin_Yields, 35.5 + adjHalfBin, ymax_Yields*0.6);
            tl_mt2->DrawLine(42.5 + adjHalfBin, ymin_Yields, 42.5 + adjHalfBin, ymax_Yields*0.6);
            tl_mt2->DrawLine(45.5 + adjHalfBin, ymin_Yields, 45.5 + adjHalfBin, ymax_Yields*0.6);
            tl_mt2->DrawLine(50.5 + adjHalfBin, ymin_Yields, 50.5 + adjHalfBin, ymax_Yields*0.6);
          }
    }
    // MT2 labels
    TLatex * ttextmt2 = new TLatex();
    ttextmt2->SetTextFont(42);
    //ttextmt2->SetTextColor(49);
    ttextmt2->SetTextColor(kBlack);
    ttextmt2->SetTextSize(0.040);
    ttextmt2->SetTextAlign(12);
    ttextmt2->SetTextAngle(90);
    if(logscale)
      {
	//Only done upto here
	if(NSB == 59)
	  {
	    ttextmt2->DrawLatex( 2.0, ymax_Yields/40. , "M_{T2}=[200,350]");
	    ttextmt2->DrawLatex( 6.0, ymax_Yields/100. , "M_{T2}=[350,450]");
	    ttextmt2->DrawLatex( 10.5, ymax_Yields/100. , "M_{T2}#geq450 GeV");
	    ttextmt2->DrawLatex( 27.5, ymax_Yields/100. , "M_{T2}#geq350 GeV");
	  }
	else if(NSB == 37 || NSB == 45)
	  {
	    ttextmt2->DrawLatex( 2.0, ymax_Yields/2000. , "M_{T2}=[200,300]");
	    ttextmt2->DrawLatex( 6.0, ymax_Yields/2000. , "M_{T2}=[300,400]");
	    ttextmt2->DrawLatex( 9.5, ymax_Yields/2000. , "M_{T2}#geq400 GeV");
	  } 
      }
    else
      {
	if(NSB == 59)
	  {
	    ttextmt2->DrawLatex( 2.0, deltaY*0.52 + ymin_Yields , "M_{T2}=[200,350]");
	    ttextmt2->DrawLatex( 6.0, deltaY*0.52 + ymin_Yields , "M_{T2}=[350,450]");
	    ttextmt2->DrawLatex( 10.5, deltaY*0.52 + ymin_Yields , "M_{T2}#geq450 GeV");
	    ttextmt2->DrawLatex( 27.5, deltaY*0.52 + ymin_Yields , "M_{T2}#geq350 GeV");
	  }
	
	else if(NSB == 37 || NSB == 45)
	  {
	    ttextmt2->DrawLatex( 2.0, deltaY*0.52 + ymin_Yields , "M_{T2}=[200,300]");
	    ttextmt2->DrawLatex( 6.0, deltaY*0.52 + ymin_Yields , "M_{T2}=[300,400]");
	    ttextmt2->DrawLatex( 9.5, deltaY*0.52 + ymin_Yields , "M_{T2}#geq400 GeV");
	  }
      }


}

void Plot_pred_ICHEP()
{

  const unsigned nsb=84;
  double Pred[nsb];
Pred[0] = 205.62219;
Pred[1] = 21.18378;
Pred[2] = 3.69890;
Pred[3] = 1.27805;
Pred[4] = 22.42201;
Pred[5] = 5.57040;
Pred[6] = 1.34808;
Pred[7] = 0.80865;
Pred[8] = 0.81858;
Pred[9] = 5.06885;
Pred[10] = 0.37940;
Pred[11] = 1.16050;
Pred[12] = 155.01298;
Pred[13] = 14.40568;
Pred[14] = 0.00000;
Pred[15] = 0.23998;
Pred[16] = 13.26749;
Pred[17] = 4.04186;
Pred[18] = 0.00000;
Pred[19] = 0.75520;
Pred[20] = 0.53093;
Pred[21] = 1.21581;
Pred[22] = 0.00000;
Pred[23] = 25.15823;
Pred[24] = 2.47372;
Pred[25] = 0.00000;
Pred[26] = 4.39193;
Pred[27] = 0.67475;
Pred[28] = 54.41609;
Pred[29] = 5.27224;
Pred[30] = 0.24601;
Pred[31] = 0.00000;
Pred[32] = 6.12253;
Pred[33] = 0.78582;
Pred[34] = 0.52118;
Pred[35] = 0.00000;
Pred[36] = 0.00000;
Pred[37] = 0.83659;
Pred[38] = 0.75205;
Pred[39] = 1.31496;
Pred[40] = 37.68324;
Pred[41] = 5.21535;
Pred[42] = 0.56579;
Pred[43] = 2.65892;
Pred[44] = 1.23919;
Pred[45] = 0.68652;
Pred[46] = 0.97852;
Pred[47] = 0.95411;
Pred[48] = 0.00000;
Pred[49] = 7.27214;
Pred[50] = 0.35616;
Pred[51] = 0.00000;
Pred[52] = 0.00000;
Pred[53] = 1.63515;
Pred[54] = 0.00000;
Pred[55] = 1.19445;
Pred[56] = 0.00000;
Pred[57] = 0.00000;
Pred[58] = 0.00000;



  double uncUpPred[nsb];
  double uncDownPred[nsb];

uncUpPred[0] = 11.97116;
uncDownPred[0] = 11.93159;
uncUpPred[1] = 4.23288;
uncDownPred[1] = 4.04097;
uncUpPred[2] = 2.09447;
uncDownPred[2] = 1.71820;
uncUpPred[3] = 1.60135;
uncDownPred[3] = 0.95535;
uncUpPred[4] = 4.29821;
uncDownPred[4] = 4.18929;
uncUpPred[5] = 2.15959;
uncDownPred[5] = 1.79942;
uncUpPred[6] = 1.28087;
uncDownPred[6] = 0.78208;
uncUpPred[7] = 1.22960;
uncDownPred[7] = 0.80865;
uncUpPred[8] = 0.77766;
uncDownPred[8] = 0.60703;
uncUpPred[9] = 2.66190;
uncDownPred[9] = 2.42326;
uncUpPred[10] = 1.13830;
uncDownPred[10] = 0.37940;
uncUpPred[11] = 1.42336;
uncDownPred[11] = 0.67095;
uncUpPred[12] = 10.76360;
uncDownPred[12] = 10.71805;
uncUpPred[13] = 3.32287;
uncDownPred[13] = 3.08320;
uncUpPred[14] = 1.27473;
uncDownPred[14] = 0.00000;
uncUpPred[15] = 0.54177;
uncDownPred[15] = 0.23998;
uncUpPred[16] = 4.17646;
uncDownPred[16] = 4.06719;
uncUpPred[17] = 1.77479;
uncDownPred[17] = 1.51024;
uncUpPred[18] = 1.80059;
uncDownPred[18] = 0.00000;
uncUpPred[19] = 2.06521;
uncDownPred[19] = 0.75520;
uncUpPred[20] = 1.34473;
uncDownPred[20] = 0.53093;
uncUpPred[21] = 1.69991;
uncDownPred[21] = 0.86333;
uncUpPred[22] = 0.89120;
uncDownPred[22] = 0.00000;
uncUpPred[23] = 4.58035;
uncDownPred[23] = 4.47277;
uncUpPred[24] = 1.76191;
uncDownPred[24] = 1.31721;
uncUpPred[25] = 1.29548;
uncDownPred[25] = 0.00000;
uncUpPred[26] = 3.17373;
uncDownPred[26] = 2.86906;
uncUpPred[27] = 1.65726;
uncDownPred[27] = 0.67475;
uncUpPred[28] = 4.68248;
uncDownPred[28] = 4.63731;
uncUpPred[29] = 1.54289;
uncDownPred[29] = 1.40420;
uncUpPred[30] = 0.80041;
uncDownPred[30] = 0.24601;
uncUpPred[31] = 0.41258;
uncDownPred[31] = 0.00000;
uncUpPred[32] = 2.18757;
uncDownPred[32] = 1.87757;
uncUpPred[33] = 1.17711;
uncDownPred[33] = 0.55566;
uncUpPred[34] = 1.41961;
uncDownPred[34] = 0.36942;
uncUpPred[35] = 0.95377;
uncDownPred[35] = 0.00000;
uncUpPred[36] = 0.86950;
uncDownPred[36] = 0.00000;
uncUpPred[37] = 1.04210;
uncDownPred[37] = 0.61258;
uncUpPred[38] = 1.41619;
uncDownPred[38] = 0.75205;
uncUpPred[39] = 1.49657;
uncDownPred[39] = 0.98803;
uncUpPred[40] = 4.44353;
uncDownPred[40] = 4.40273;
uncUpPred[41] = 2.33520;
uncDownPred[41] = 2.26295;
uncUpPred[42] = 0.89776;
uncDownPred[42] = 0.40033;
uncUpPred[43] = 1.49570;
uncDownPred[43] = 1.10226;
uncUpPred[44] = 1.20917;
uncDownPred[44] = 0.71943;
uncUpPred[45] = 1.24124;
uncDownPred[45] = 0.68652;
uncUpPred[46] = 1.01914;
uncDownPred[46] = 0.57596;
uncUpPred[47] = 1.60704;
uncDownPred[47] = 0.95411;
uncUpPred[48] = 1.35798;
uncDownPred[48] = 0.00000;
uncUpPred[49] = 2.06321;
uncDownPred[49] = 1.94062;
uncUpPred[50] = 0.50773;
uncDownPred[50] = 0.25424;
uncUpPred[51] = 1.07690;
uncDownPred[51] = 0.00000;
uncUpPred[52] = 0.91775;
uncDownPred[52] = 0.00000;
uncUpPred[53] = 1.11923;
uncDownPred[53] = 0.89839;
uncUpPred[54] = 0.65075;
uncDownPred[54] = 0.00000;
uncUpPred[55] = 1.02230;
uncDownPred[55] = 0.75590;
uncUpPred[56] = 0.47961;
uncDownPred[56] = 0.00000;
uncUpPred[57] = 0.81399;
uncDownPred[57] = 0.00000;
uncUpPred[58] = 0.32653;
uncDownPred[58] = 0.00000;


  const std::string titre="CMS Supplementary";

  TCanvas *c1 = new TCanvas("c1", "c1",0,51,1920,1004);
  c1->SetFillColor(0);
  c1->cd();
  //c1 -> SetLogy();

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

  Double_t x[nsb];
  Double_t exl[nsb];
  Double_t exh[nsb];
  for (int xc=0;xc<nsb;++xc)
  {
    x[xc]=xc+0.5;
    exl[xc]=0.5;
    exh[xc]=0.5;
  }

   TGraphAsymmErrors* gr = new TGraphAsymmErrors(nsb,x,Pred,exl,exh,uncDownPred,uncUpPred);
   gr->SetTitle("");
   gr->SetLineWidth(3);
   gr->SetLineColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("AP");
   gr->SetMaximum(220.0);
   gr->SetMinimum(0.0);
   //gr->SetMaximum(1000.0);
   //gr->SetMinimum(0.1);
   gr->GetXaxis()->SetTitle("Search region bin number");
   gr->GetXaxis()->SetTitleSize(0.045);
   gr->GetYaxis()->SetTitle("Events");
   gr->GetYaxis()->SetTitleSize(0.045);
   gr->GetXaxis()->SetLimits(0.0,59.0);

   drawSBregionDef(0.0, 220.0, false);
   //drawSBregionDef(0.1, 1000.0, true);

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
