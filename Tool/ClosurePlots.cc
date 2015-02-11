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
  myClosurePlots.Initialization();
  //initialize the closure plots parameter we want to investigate
  vector<Plotting_Parameter> myPlotting_Paramete = 
  { 
    //muon closure plots, all
    {"_mu_all_met"     , "MET[GeV]"     ,200 , 600 },
    {"_mu_all_njets"   , "NJets30"      ,2   , 17  }, 
    {"_mu_all_mt2"     , "MT2[GeV]"     ,100 , 600 },                                               
    {"_mu_all_topmass" , "TopMass[GeV]" ,50  , 300 },
    //muon closure plots, acc
    {"_mu_acc_met"     , "MET[GeV]"     ,200 , 600 },
    {"_mu_acc_njets"   , "NJets30"      ,2   , 17  },
    {"_mu_acc_mt2"     , "MT2[GeV]"     ,100 , 600 },
    {"_mu_acc_topmass" , "TopMass[GeV]" ,50  , 300 },
    //muon closure plots, id
    {"_mu_id_met"     , "MET[GeV]"     ,200 , 600 },
    {"_mu_id_njets"   , "NJets30"      ,2   , 17  },
    {"_mu_id_mt2"     , "MT2[GeV]"     ,100 , 600 },
    {"_mu_id_topmass" , "TopMass[GeV]" ,50  , 300 },
    //muon closure plots, iso
    {"_mu_iso_met"     , "MET[GeV]"     ,200 , 600 },
    {"_mu_iso_njets"   , "NJets30"      ,2   , 17  },
    {"_mu_iso_mt2"     , "MT2[GeV]"     ,100 , 600 },
    {"_mu_iso_topmass" , "TopMass[GeV]" ,50  , 300 },
    //electron closure plots, all
    {"_el_all_met"     , "MET[GeV]"     ,200 , 600 },
    {"_el_all_njets"   , "NJets30"      ,2   , 17  },
    {"_el_all_mt2"     , "MT2[GeV]"     ,100 , 600 },
    {"_el_all_topmass" , "TopMass[GeV]" ,50  , 300 },
    //electron closure plots, acc
    {"_el_acc_met"     , "MET[GeV]"     ,200 , 600 },
    {"_el_acc_njets"   , "NJets30"      ,2   , 17  },
    {"_el_acc_mt2"     , "MT2[GeV]"     ,100 , 600 },
    {"_el_acc_topmass" , "TopMass[GeV]" ,50  , 300 },
    //electron closure plots, id
    {"_el_id_met"     , "MET[GeV]"     ,200 , 600 },
    {"_el_id_njets"   , "NJets30"      ,2   , 17  },
    {"_el_id_mt2"     , "MT2[GeV]"     ,100 , 600 },
    {"_el_id_topmass" , "TopMass[GeV]" ,50  , 300 },
    //electron closure plots, iso
    {"_el_iso_met"     , "MET[GeV]"     ,200 , 600 },
    {"_el_iso_njets"   , "NJets30"      ,2   , 17  },
    {"_el_iso_mt2"     , "MT2[GeV]"     ,100 , 600 },
    {"_el_iso_topmass" , "TopMass[GeV]" ,50  , 300 },
  };

  
  vector<Plotting_Parameter>::iterator iter_p;

  for( iter_p = myPlotting_Paramete.begin() ; iter_p != myPlotting_Paramete.end() ; iter_p ++)
  {

    myClosurePlots.ClosureTemplate(
                                   (*iter_p).hist_tag,
                                   (*iter_p).XTitle,
                                   (*iter_p).min,
                                   (*iter_p).max 
                                  );
  }

  return 0;
}

