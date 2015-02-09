#include <vector>
#include <iostream>
#include <string>
#include <cstring>
#include <set>

#include "TFile.h"
#include "TList.h"
#include "TString.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include "ClosurePlots.h"

using namespace std;

int main()
{
  ClosurePlots myClosurePlots;
  //initialize the closure plots we want to investigate
  vector<Plotting_Parameter> myPlotting_Paramete = 
  { 
    {"mu_all_met"     , 200 , 600 },
    {"mu_all_njets"   , 2   , 17  }, 
    {"mu_all_mt2"     , 100 , 600 },                                               
    {"mu_all_topmass" , 50  , 300 }                                                
  };

  
  vector<Plotting_Parameter>::iterator iter_p;

  for( iter_p = myPlotting_Paramete.begin() ; iter_p != myPlotting_Paramete.end() ; iter_p ++)
  {

    myClosurePlots.ClosureTemplate(
                                 (*iter_p).hist_tag,
                                 (*iter_p).min,
                                 (*iter_p).max 
                                );
  }

  return 0;
}

