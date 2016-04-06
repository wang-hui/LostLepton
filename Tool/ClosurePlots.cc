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
  //myClosurePlots.SetScale(11339232, 831.76, 2154.5);
  //initialize the closure plots parameter we want to investigate
  vector<Plotting_Parameter> myPlotting_Paramete = 
  { 
    //muon closure plots, all
    {"_mu_all_met"     , "MET [GeV]"             ,200 , 600  },
    {"_mu_all_njets"   , "Number of 30 GeV jets" ,2   , 17   }, 
    {"_mu_all_mt2"     , "MT2 [GeV]"             ,100 , 600  },                                               
    {"_mu_all_ht"      , "HT [GeV]"              ,400 , 2000 },
    {"_mu_all_mht"     , "MHT [GeV]"             ,100 , 500  },
    {"_mu_all_ntopjets", "Number of Top jets"    ,0   , 5    },
    //muon closure plots, acc
    {"_mu_acc_met"     , "MET [GeV]"             ,200 , 600  },
    {"_mu_acc_njets"   , "Number of 30 GeV jets" ,2   , 17   },
    {"_mu_acc_mt2"     , "MT2 [GeV]"             ,100 , 600  },
    //{"_mu_acc_ht"      , "HT [GeV]"              ,100 , 1100 },
    {"_mu_acc_ht"      , "HT [GeV]"              ,400 , 2000 },
    {"_mu_acc_mht"     , "MHT [GeV]"             ,100 , 500  },
    {"_mu_acc_ntopjets", "Number of Top jets"    ,0   , 5    },
    //muon closure plots, id
    {"_mu_id_met"     , "MET [GeV]"             ,200 , 600  },
    {"_mu_id_njets"   , "Number of 30 GeV jets" ,2   , 17   },
    {"_mu_id_mt2"     , "MT2 [GeV]"             ,100 , 600  },
    {"_mu_id_ht"      , "HT [GeV]"              ,400 , 2000 },
    {"_mu_id_mht"     , "MHT [GeV]"             ,100 , 500  },
    {"_mu_id_ntopjets", "Number of Top jets"    ,0   , 5    },
    //muon closure plots, iso
    {"_mu_iso_met"     , "MET [GeV]"             ,200 , 600  },
    {"_mu_iso_njets"   , "Number of 30 GeV jets" ,2   , 17   },
    {"_mu_iso_mt2"     , "MT2 [GeV]"             ,100 , 600  },
    {"_mu_iso_ht"      , "HT [GeV]"              ,400 , 2000 },
    {"_mu_iso_mht"     , "MHT [GeV]"             ,100 , 500  },
    {"_mu_iso_ntopjets", "Number of Top jets"    ,0   , 5    },
    //electron closure plots, all
    {"_el_all_met"     , "MET [GeV]"             ,200 , 600  },
    {"_el_all_njets"   , "Number of 30 GeV jets" ,2   , 17   },
    {"_el_all_mt2"     , "MT2 [GeV]"             ,100 , 600  },
    {"_el_all_ht"      , "HT [GeV]"              ,400 , 2000 },
    {"_el_all_mht"     , "MHT [GeV]"             ,100 , 500  },
    {"_el_all_ntopjets", "Number of Top jets"    ,0   , 5    },
    //electron closure plots, acc
    {"_el_acc_met"     , "MET [GeV]"             ,200 , 600  },
    {"_el_acc_njets"   , "Number of 30 GeV jets" ,2   , 17   },
    {"_el_acc_mt2"     , "MT2 [GeV]"             ,100 , 600  },
    {"_el_acc_ht"      , "HT [GeV]"              ,400 , 2000 },
    {"_el_acc_mht"     , "MHT [GeV]"             ,100 , 500  },
    {"_el_acc_ntopjets", "Number of Top jets"    ,0   , 5    },
    //electron closure plots, id
    {"_el_id_met"     , "MET [GeV]"             ,200 , 600  },
    {"_el_id_njets"   , "Number of 30 GeV jets" ,2   , 17   },
    {"_el_id_mt2"     , "MT2 [GeV]"             ,100 , 600  },
    {"_el_id_ht"      , "HT [GeV]"              ,400 , 2000 },
    {"_el_id_mht"     , "MHT [GeV]"             ,100 , 500  },
    {"_el_id_ntopjets", "Number of Top jets"    ,0   , 5    },
    //electron closure plots, iso
    {"_el_iso_met"     , "MET [GeV]"             ,200 , 600  },
    {"_el_iso_njets"   , "Number of 30 GeV jets" ,2   , 17   },
    {"_el_iso_mt2"     , "MT2 [GeV]"             ,100 , 600  },
    {"_el_iso_ht"      , "HT [GeV]"              ,400 , 2000 },
    {"_el_iso_mht"     , "MHT [GeV]"             ,100 , 500  },
    {"_el_iso_ntopjets", "Number of Top jets"    ,0   , 5    },
    //search bin closure plots
    {"_mu_sb"          , "Search Bins"           ,0   , 45  },
    {"_mu_iso_sb"          , "Search Bins"           ,0   , 45  },
    {"_mu_reco_sb"          , "Search Bins"           ,0   , 45  },
    {"_mu_acc_sb"          , "Search Bins"           ,0   , 45  },
    {"_el_sb"          , "Search Bins"           ,0   , 45  },
    {"_lept_sb"        , "Search Bins"           ,0   , 45  },
    {"_lept_isotrk" , "Search Bins"           ,0   , 45  },
    // mu+ele all closure plots
    {"_lept_all_met"     , "MET [GeV]"             ,200 , 600  },
    {"_lept_all_njets"   , "Number of 30 GeV jets" ,2   , 17   }, 
    {"_lept_all_mt2"     , "MT2 [GeV]"             ,100 , 600  },                                               
    {"_lept_all_ht"      , "HT [GeV]"              ,400 , 2000 },
    {"_lept_all_mht"     , "MHT [GeV]"             ,100 , 500  },
    {"_lept_all_ntopjets", "Number of Top jets"    ,0   , 5    },

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
  
  myClosurePlots.DiLeptonPlots(
                              "h_exp_musingle_all_met",
                              "h_exp_mu_all_met",
                              "MET [GeV]",
                              200,
                              600
                              );
  myClosurePlots.DiLeptonPlots(
                              "h_exp_musingle_all_njets",
                              "h_exp_mu_all_njets",
                              "Number of 30 GeV jets",
                              2,
                              17
                              );
  myClosurePlots.DiLeptonPlots(
                              "h_exp_musingle_all_mt2",
                              "h_exp_mu_all_mt2",
                              "MT2 [GeV]",
                              100,
                              600
                              );
  
  myClosurePlots.DiLeptonPlots(
                              "h_exp_musingle_all_ht",
                              "h_exp_mu_all_ht",
                              "HT [GeV]",
                              400,
                              2000
                              );
  myClosurePlots.DiLeptonPlots(
                              "h_exp_musingle_all_mht",
                              "h_exp_mu_all_mht",
                              "MHT [GeV]",
                              100,
                              500
                              );

  myClosurePlots.DiLeptonPlots(
                              "h_exp_elsingle_all_met",
                              "h_exp_el_all_met",
                              "MET [GeV]",
                              200,
                              600
                              );
  myClosurePlots.DiLeptonPlots(
                              "h_exp_elsingle_all_njets",
                              "h_exp_el_all_njets",
                              "Number of 30 GeV jets",
                              2,
                              17
                              );
  myClosurePlots.DiLeptonPlots(
                              "h_exp_elsingle_all_mt2",
                              "h_exp_el_all_mt2",
                              "MT2 [GeV]",
                              100,
                              600
                              );

  myClosurePlots.DiLeptonPlots(
                              "h_exp_elsingle_all_ht",
                              "h_exp_el_all_ht",
                              "HT [GeV]",
                              400,
                              2000
                              );
  myClosurePlots.DiLeptonPlots(
                              "h_exp_elsingle_all_mht",
                              "h_exp_el_all_mht",
                              "MHT [GeV]",
                              100,
                              500
                              );

  return 0;
}

