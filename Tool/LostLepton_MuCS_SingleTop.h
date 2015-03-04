#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"

#include "SusyAnaTools/Tools/NTupleReader.h"

#define PT_BINS 8
//#define AC_BINS 1

//############################begin to defin class AccRecoIsoEffs###################//


class AccRecoIsoEffs
{
 public:
  void GetAccRecoIsoEffs();
  void GetDiLeptonFactor();
  void printOverview();
  void printAccRecoIsoEffs();
  void printNormalizeFlowNumber();

  //define the number of Pt bins we need in the calculation
  //static int PT_BINS = 8;

  //here we define the overall information for ttbar sample
  int nevents_tot = 0;
  int nevents_sel_base = 0;
  int nevents_sel_mus = 0;
  int nevents_sel_els = 0;

  //define acceptance, reco eff and iso eff to be passed from ttbar result
  double mus_acc = 0, els_acc = 0;
  double mus_acc_err = 0, els_acc_err = 0;
  double mus_recoeff[PT_BINS] = {0}, els_recoeff[PT_BINS] = {0};
  double mus_isoeff[PT_BINS] = {0}, els_isoeff[PT_BINS] = {0};
  double mus_recoeff_err[PT_BINS] = {0}, els_recoeff_err[PT_BINS] = {0};
  double mus_isoeff_err[PT_BINS] = {0}, els_isoeff_err[PT_BINS] = {0};

  //here we define the muon/electron number we need to count in the loop
  double nmus = 0, nmus_acc = 0, nels = 0, nels_acc =0;

  //here we define the event weight we are going to use in the second loop ( muon/electron CS and prediction plots)
  double mus_EventWeight_iso[PT_BINS] = {0}, mus_EventWeight_reco[PT_BINS] = {0}, mus_EventWeight_acc[PT_BINS] = {0};
  double els_EventWeight_iso[PT_BINS] = {0}, els_EventWeight_reco[PT_BINS] = {0}, els_EventWeight_acc[PT_BINS] = {0};

  //di-lepton correction
  double nevents_single_mus = 0, nevents_di_mus = 0;
  double nevents_single_els = 0, nevents_di_els = 0;
  
  //flow number table
  double nevents_exp_all_mus = 0, nevents_exp_all_els = 0;
  double nevents_exp_acc_mus = 0, nevents_exp_acc_els = 0;
  double nevents_exp_id_mus = 0, nevents_exp_id_els = 0;
  double nevents_exp_iso_mus = 0, nevents_exp_iso_els = 0;
  double nevents_exp_all_mus_err = 0, nevents_exp_all_els_err = 0;
  double nevents_exp_acc_mus_err = 0, nevents_exp_acc_els_err = 0;
  double nevents_exp_id_mus_err = 0, nevents_exp_id_els_err = 0;
  double nevents_exp_iso_mus_err = 0, nevents_exp_iso_els_err = 0;

  double nevents_pred_all_mus = 0, nevents_pred_all_els = 0;
  double nevents_pred_acc_mus = 0, nevents_pred_acc_els = 0;
  double nevents_pred_id_mus = 0, nevents_pred_id_els = 0;
  double nevents_pred_iso_mus = 0, nevents_pred_iso_els = 0;
  double nevents_pred_all_mus_err = 0, nevents_pred_all_els_err = 0;
  double nevents_pred_acc_mus_err = 0, nevents_pred_acc_els_err = 0;
  double nevents_pred_id_mus_err = 0, nevents_pred_id_els_err = 0;
  double nevents_pred_iso_mus_err = 0, nevents_pred_iso_els_err = 0;

  //di lepton correction factor
  double corrfactor_di_mus = 0;
  double corrfactor_di_els = 0;
  double corrfactor_di_mus_err = 0;
  double corrfactor_di_els_err = 0;

  void EffstoWeights();
  void NormalizeFlowNumber();

 private:
  double get_stat_Error(double a,
                        double an
                       );
  
  double get_sys_Error(double r,
                       double p
                      );
};


double AccRecoIsoEffs::get_stat_Error(double a,
                                      double an
                                     )
{
  double n;
  n = an - a;

  double err;
  err = 1000;

  double alpha;
  alpha = 1-0.6827;

  if( a>=0 && n>=0 )
  {
    err = std::sqrt(n/(a+n)/(a+n)*n/(a+n)/(a+n)*ROOT::Math::gamma_quantile_c(alpha/2,a+1,1)+a/(a+n)/(a+n)*a/(a+n)/(a+n)*ROOT::Math::gamma_quantile_c(alpha/2,n+1,1));
    return err;
  }
  else
  {
    return -1;
  }
}

double AccRecoIsoEffs::get_sys_Error(double r,
                                     double p
                                    )
{
  double err;
  err = 1000;

  err = r - (1-(1-r)*(1+p));

  if(err < 0)
  {
    return -1;
  }
  return err;
}

//############finish the definition of class AccRecoEffs######################

class BaseHistgram
{
 public:
  void BookHistgram(const char *);

  TFile *oFile;
  TH1D *h_b_all_MET;
  TH1D *h_b_baseline_nMuons, *h_b_baseline_njets, *h_b_baseline_nbjetsCSVM, *h_b_baseline_bestTopMass, *h_b_baseline_MET, *h_b_baseline_jetpt2, *h_b_baseline_jetpt4, *h_b_baseline_jet1_met_phi_diff, *h_b_baseline_jet2_met_phi_diff, *h_b_baseline_jet3_met_phi_diff;
  TH1D *h_b_acc_njets, *h_b_acc_nbjetsCSVM, *h_b_acc_bestTopMass, *h_b_acc_MET, *h_b_acc_jetpt2, *h_b_acc_jetpt4, *h_b_acc_jet1_met_phi_diff, *h_b_acc_jet2_met_phi_diff, *h_b_acc_jet3_met_phi_diff;
  TH1D *h_b_reco_nMuons, *h_b_reco_njets, *h_b_reco_nbjetsCSVM, *h_b_reco_bestTopMass, *h_b_reco_MET, *h_b_reco_jetpt2, *h_b_reco_jetpt4, *h_b_reco_jet1_met_phi_diff, *h_b_reco_jet2_met_phi_diff, *h_b_reco_jet3_met_phi_diff;
  TH1D *h_b_deltaR_mus, *h_b_deltaR_els;

  //closure plots definition
  TH1D *h_pred_mu_iso_met, *h_pred_mu_iso_njets, *h_pred_mu_iso_mt2, *h_pred_mu_iso_topmass, *h_pred_mu_iso_ht, *h_pred_mu_iso_mht, *h_pred_mu_iso_ntopjets;
  TH1D *h_pred_mu_id_met, *h_pred_mu_id_njets, *h_pred_mu_id_mt2, *h_pred_mu_id_topmass, *h_pred_mu_id_ht, *h_pred_mu_id_mht, *h_pred_mu_id_ntopjets;
  TH1D *h_pred_mu_acc_met, *h_pred_mu_acc_njets, *h_pred_mu_acc_mt2, *h_pred_mu_acc_topmass, *h_pred_mu_acc_ht, *h_pred_mu_acc_mht, *h_pred_mu_acc_ntopjets;
  TH1D *h_pred_mu_all_met, *h_pred_mu_all_njets, *h_pred_mu_all_mt2, *h_pred_mu_all_topmass, *h_pred_mu_all_ht, *h_pred_mu_all_mht, *h_pred_mu_all_ntopjets;

  TH1D *h_exp_mu_iso_met, *h_exp_mu_iso_njets, *h_exp_mu_iso_mt2, *h_exp_mu_iso_topmass, *h_exp_mu_iso_ht, *h_exp_mu_iso_mht, *h_exp_mu_iso_ntopjets;
  TH1D *h_exp_mu_id_met, *h_exp_mu_id_njets, *h_exp_mu_id_mt2, *h_exp_mu_id_topmass, *h_exp_mu_id_ht, *h_exp_mu_id_mht, *h_exp_mu_id_ntopjets;
  TH1D *h_exp_mu_acc_met, *h_exp_mu_acc_njets, *h_exp_mu_acc_mt2, *h_exp_mu_acc_topmass, *h_exp_mu_acc_ht, *h_exp_mu_acc_mht, *h_exp_mu_acc_ntopjets;
  TH1D *h_exp_mu_all_met, *h_exp_mu_all_njets, *h_exp_mu_all_mt2, *h_exp_mu_all_topmass, *h_exp_mu_all_ht, *h_exp_mu_all_mht, *h_exp_mu_all_ntopjets;

  TH1D *h_exp_musingle_all_met, *h_exp_musingle_all_njets, *h_exp_musingle_all_mt2, *h_exp_musingle_all_topmass, *h_exp_musingle_all_ht, *h_exp_musingle_all_mht, *h_exp_musingle_all_ntopjets;

  TH1D *h_pred_el_iso_met, *h_pred_el_iso_njets, *h_pred_el_iso_mt2, *h_pred_el_iso_topmass, *h_pred_el_iso_ht, *h_pred_el_iso_mht, *h_pred_el_iso_ntopjets;
  TH1D *h_pred_el_id_met, *h_pred_el_id_njets, *h_pred_el_id_mt2, *h_pred_el_id_topmass, *h_pred_el_id_ht, *h_pred_el_id_mht, *h_pred_el_id_ntopjets;
  TH1D *h_pred_el_acc_met, *h_pred_el_acc_njets, *h_pred_el_acc_mt2, *h_pred_el_acc_topmass, *h_pred_el_acc_ht, *h_pred_el_acc_mht, *h_pred_el_acc_ntopjets;
  TH1D *h_pred_el_all_met, *h_pred_el_all_njets, *h_pred_el_all_mt2, *h_pred_el_all_topmass, *h_pred_el_all_ht, *h_pred_el_all_mht, *h_pred_el_all_ntopjets;

  TH1D *h_exp_el_iso_met, *h_exp_el_iso_njets, *h_exp_el_iso_mt2, *h_exp_el_iso_topmass, *h_exp_el_iso_ht, *h_exp_el_iso_mht, *h_exp_el_iso_ntopjets;
  TH1D *h_exp_el_id_met, *h_exp_el_id_njets, *h_exp_el_id_mt2, *h_exp_el_id_topmass, *h_exp_el_id_ht, *h_exp_el_id_mht, *h_exp_el_id_ntopjets;
  TH1D *h_exp_el_acc_met, *h_exp_el_acc_njets, *h_exp_el_acc_mt2, *h_exp_el_acc_topmass, *h_exp_el_acc_ht, *h_exp_el_acc_mht, *h_exp_el_acc_ntopjets;
  TH1D *h_exp_el_all_met, *h_exp_el_all_njets, *h_exp_el_all_mt2, *h_exp_el_all_topmass, *h_exp_el_all_ht, *h_exp_el_all_mht, *h_exp_el_all_ntopjets;

  TH1D *h_exp_elsingle_all_met, *h_exp_elsingle_all_njets, *h_exp_elsingle_all_mt2, *h_exp_elsingle_all_topmass, *h_exp_elsingle_all_ht, *h_exp_elsingle_all_mht, *h_exp_elsingle_all_ntopjets;
};

void BaseHistgram::BookHistgram(const char *outFileName)
{

  oFile = new TFile(outFileName, "recreate");

  h_b_all_MET = new TH1D("h_b_all_MET","",1000,0,1000);

  h_b_baseline_nMuons = new TH1D("h_b_baseline_nMuons","",10,0,10);
  h_b_baseline_njets = new TH1D("h_b_baseline_njets","",10,0,10);
  h_b_baseline_nbjetsCSVM = new TH1D("h_b_baseline_nbjetsCSVM","",10,0,10);
  h_b_baseline_bestTopMass = new TH1D("h_b_baseline_bestTopMass","",1000,0,500);
  h_b_baseline_MET = new TH1D("h_b_baseline_MET","",1000,0,1000);
  h_b_baseline_jetpt4 = new TH1D("h_b_baseline_jetpt4","",1000,0,1000);
  h_b_baseline_jetpt2 = new TH1D("h_b_baseline_jetpt2","",1000,0,1000);
  h_b_baseline_jet1_met_phi_diff = new TH1D("h_b_baseline_jet1_met_phi_diff","",1000,-5,5);
  h_b_baseline_jet2_met_phi_diff = new TH1D("h_b_baseline_jet2_met_phi_diff","",1000,-5,5);
  h_b_baseline_jet3_met_phi_diff = new TH1D("h_b_baseline_jet3_met_phi_diff","",1000,-5,5);

  h_b_acc_njets = new TH1D("h_b_acc_njets","",10,0,10);
  h_b_acc_nbjetsCSVM = new TH1D("h_b_acc_nbjetsCSVM","",10,0,10);
  h_b_acc_bestTopMass = new TH1D("h_b_acc_bestTopMass","",1000,0,500);
  h_b_acc_MET = new TH1D("h_b_acc_MET","",1000,0,1000);
  h_b_acc_jetpt4 = new TH1D("h_b_acc_jetpt4","",1000,0,1000);
  h_b_acc_jetpt2 = new TH1D("h_b_acc_jetpt2","",1000,0,1000);
  h_b_acc_jet1_met_phi_diff = new TH1D("h_b_acc_jet1_met_phi_diff","",1000,-5,5);
  h_b_acc_jet2_met_phi_diff = new TH1D("h_b_acc_jet2_met_phi_diff","",1000,-5,5);
  h_b_acc_jet3_met_phi_diff = new TH1D("h_b_acc_jet3_met_phi_diff","",1000,-5,5);

  h_b_reco_nMuons = new TH1D("h_b_reco_nMuons","",10,0,10);
  h_b_reco_njets = new TH1D("h_b_reco_njets","",10,0,10);
  h_b_reco_nbjetsCSVM = new TH1D("h_b_reco_nbjetsCSVM","",10,0,10);
  h_b_reco_bestTopMass = new TH1D("h_b_reco_bestTopMass","",1000,0,500);
  h_b_reco_MET = new TH1D("h_b_reco_MET","",1000,0,1000);
  h_b_reco_jetpt4 = new TH1D("h_b_reco_jetpt4","",1000,0,1000);
  h_b_reco_jetpt2 = new TH1D("h_b_reco_jetpt2","",1000,0,1000);
  h_b_reco_jet1_met_phi_diff = new TH1D("h_b_reco_jet1_met_phi_diff","",1000,-5,5);
  h_b_reco_jet2_met_phi_diff = new TH1D("h_b_reco_jet2_met_phi_diff","",1000,-5,5);
  h_b_reco_jet3_met_phi_diff = new TH1D("h_b_reco_jet3_met_phi_diff","",1000,-5,5);

  h_b_deltaR_mus = new TH1D("h_b_deltaR_mus","",1000,0,0.5);
  h_b_deltaR_els = new TH1D("h_b_deltaR_els","",1000,0,0.5);

  //start closure plots
  h_pred_mu_iso_met = new TH1D("h_pred_mu_iso_met","",100,0,1000);
  h_pred_mu_iso_njets = new TH1D("h_pred_mu_iso_njets","",20,0,20);
  h_pred_mu_iso_mt2 = new TH1D("h_pred_mu_iso_mt2","",100,0,1000);
  h_pred_mu_iso_topmass = new TH1D("h_pred_mu_iso_topmass","",100,50,300);
  h_pred_mu_iso_ht = new TH1D("h_pred_mu_iso_ht","",300,0,3000);
  h_pred_mu_iso_mht = new TH1D("h_pred_mu_iso_mht","",100,0,1000);
  h_pred_mu_iso_ntopjets = new TH1D("h_pred_mu_iso_ntopjets","",20,0,20);

  h_pred_mu_id_met = new TH1D("h_pred_mu_id_met","",100,0,1000);
  h_pred_mu_id_njets = new TH1D("h_pred_mu_id_njets","",20,0,20);
  h_pred_mu_id_mt2 = new TH1D("h_pred_mu_id_mt2","",100,0,1000);
  h_pred_mu_id_topmass = new TH1D("h_pred_mu_id_topmass","",100,50,300);
  h_pred_mu_id_ht = new TH1D("h_pred_mu_id_ht","",300,0,3000);
  h_pred_mu_id_mht = new TH1D("h_pred_mu_id_mht","",100,0,1000);
  h_pred_mu_id_ntopjets = new TH1D("h_pred_mu_id_ntopjets","",20,0,20);

  h_pred_mu_acc_met = new TH1D("h_pred_mu_acc_met","",100,0,1000);
  h_pred_mu_acc_njets = new TH1D("h_pred_mu_acc_njets","",20,0,20);
  h_pred_mu_acc_mt2 = new TH1D("h_pred_mu_acc_mt2","",100,0,1000);
  h_pred_mu_acc_topmass = new TH1D("h_pred_mu_acc_topmass","",100,50,300);
  h_pred_mu_acc_ht = new TH1D("h_pred_mu_acc_ht","",300,0,3000);
  h_pred_mu_acc_mht = new TH1D("h_pred_mu_acc_mht","",100,0,1000);
  h_pred_mu_acc_ntopjets = new TH1D("h_pred_mu_acc_ntopjets","",20,0,20);

  h_pred_mu_all_met = new TH1D("h_pred_mu_all_met","",100,0,1000);
  h_pred_mu_all_njets = new TH1D("h_pred_mu_all_njets","",20,0,20);
  h_pred_mu_all_mt2 = new TH1D("h_pred_mu_all_mt2","",100,0,1000);
  h_pred_mu_all_topmass = new TH1D("h_pred_mu_all_topmass","",100,50,300);
  h_pred_mu_all_ht = new TH1D("h_pred_mu_all_ht","",300,0,3000);
  h_pred_mu_all_mht = new TH1D("h_pred_mu_all_mht","",100,0,1000);
  h_pred_mu_all_ntopjets = new TH1D("h_pred_mu_all_ntopjets","",20,0,20);

  h_exp_mu_iso_met = new TH1D("h_exp_mu_iso_met","",100,0,1000);
  h_exp_mu_iso_njets = new TH1D("h_exp_mu_iso_njets","",20,0,20);
  h_exp_mu_iso_mt2 = new TH1D("h_exp_mu_iso_mt2","",100,0,1000);
  h_exp_mu_iso_topmass = new TH1D("h_exp_mu_iso_topmass","",100,50,300);
  h_exp_mu_iso_ht = new TH1D("h_exp_mu_iso_ht","",300,0,3000);
  h_exp_mu_iso_mht = new TH1D("h_exp_mu_iso_mht","",100,0,1000);
  h_exp_mu_iso_ntopjets = new TH1D("h_exp_mu_iso_ntopjets","",20,0,20);

  h_exp_mu_id_met = new TH1D("h_exp_mu_id_met","",100,0,1000);
  h_exp_mu_id_njets = new TH1D("h_exp_mu_id_njets","",20,0,20);
  h_exp_mu_id_mt2 = new TH1D("h_exp_mu_id_mt2","",100,0,1000);
  h_exp_mu_id_topmass = new TH1D("h_exp_mu_id_topmass","",100,50,300);
  h_exp_mu_id_ht = new TH1D("h_exp_mu_id_ht","",300,0,3000);
  h_exp_mu_id_mht = new TH1D("h_exp_mu_id_mht","",100,0,1000);
  h_exp_mu_id_ntopjets = new TH1D("h_exp_mu_id_ntopjets","",20,0,20);

  h_exp_mu_acc_met = new TH1D("h_exp_mu_acc_met","",100,0,1000);
  h_exp_mu_acc_njets = new TH1D("h_exp_mu_acc_njets","",20,0,20);
  h_exp_mu_acc_mt2 = new TH1D("h_exp_mu_acc_mt2","",100,0,1000);
  h_exp_mu_acc_topmass = new TH1D("h_exp_mu_acc_topmass","",100,50,300);
  h_exp_mu_acc_ht = new TH1D("h_exp_mu_acc_ht","",300,0,3000);
  h_exp_mu_acc_mht = new TH1D("h_exp_mu_acc_mht","",100,0,1000);
  h_exp_mu_acc_ntopjets = new TH1D("h_exp_mu_acc_ntopjets","",20,0,20);

  h_exp_mu_all_met = new TH1D("h_exp_mu_all_met","",100,0,1000);
  h_exp_mu_all_njets = new TH1D("h_exp_mu_all_njets","",20,0,20);
  h_exp_mu_all_mt2 = new TH1D("h_exp_mu_all_mt2","",100,0,1000);
  h_exp_mu_all_topmass = new TH1D("h_exp_mu_all_topmass","",100,50,300);
  h_exp_mu_all_ht = new TH1D("h_exp_mu_all_ht","",300,0,3000);
  h_exp_mu_all_mht = new TH1D("h_exp_mu_all_mht","",100,0,1000);
  h_exp_mu_all_ntopjets = new TH1D("h_exp_mu_all_ntopjets","",20,0,20);

  h_exp_musingle_all_met = new TH1D("h_exp_musingle_all_met","",100,0,1000);
  h_exp_musingle_all_njets = new TH1D("h_exp_musingle_all_njets","",20,0,20);
  h_exp_musingle_all_mt2 = new TH1D("h_exp_musingle_all_mt2","",100,0,1000);
  h_exp_musingle_all_topmass = new TH1D("h_exp_musingle_all_topmass","",100,50,300);
  h_exp_musingle_all_ht = new TH1D("h_exp_musingle_all_ht","",300,0,3000);
  h_exp_musingle_all_mht = new TH1D("h_exp_musingle_all_mht","",100,0,1000);
  h_exp_musingle_all_ntopjets = new TH1D("h_exp_musingle_all_ntopjets","",20,0,20);

  h_pred_el_iso_met = new TH1D("h_pred_el_iso_met","",100,0,1000);
  h_pred_el_iso_njets = new TH1D("h_pred_el_iso_njets","",20,0,20);
  h_pred_el_iso_mt2 = new TH1D("h_pred_el_iso_mt2","",100,0,1000);
  h_pred_el_iso_topmass = new TH1D("h_pred_el_iso_topmass","",100,50,300);
  h_pred_el_iso_ht = new TH1D("h_pred_el_iso_ht","",300,0,3000);
  h_pred_el_iso_mht = new TH1D("h_pred_el_iso_mht","",100,0,1000);
  h_pred_el_iso_ntopjets = new TH1D("h_pred_el_iso_ntopjets","",20,0,20);

  h_pred_el_id_met = new TH1D("h_pred_el_id_met","",100,0,1000);
  h_pred_el_id_njets = new TH1D("h_pred_el_id_njets","",20,0,20);
  h_pred_el_id_mt2 = new TH1D("h_pred_el_id_mt2","",100,0,1000);
  h_pred_el_id_topmass = new TH1D("h_pred_el_id_topmass","",100,50,300);
  h_pred_el_id_ht = new TH1D("h_pred_el_id_ht","",300,0,3000);
  h_pred_el_id_mht = new TH1D("h_pred_el_id_mht","",100,0,1000);
  h_pred_el_id_ntopjets = new TH1D("h_pred_el_id_ntopjets","",20,0,20);

  h_pred_el_acc_met = new TH1D("h_pred_el_acc_met","",100,0,1000);
  h_pred_el_acc_njets = new TH1D("h_pred_el_acc_njets","",20,0,20);
  h_pred_el_acc_mt2 = new TH1D("h_pred_el_acc_mt2","",100,0,1000);
  h_pred_el_acc_topmass = new TH1D("h_pred_el_acc_topmass","",100,50,300);
  h_pred_el_acc_ht = new TH1D("h_pred_el_acc_ht","",300,0,3000);
  h_pred_el_acc_mht = new TH1D("h_pred_el_acc_mht","",100,0,1000);
  h_pred_el_acc_ntopjets = new TH1D("h_pred_el_acc_ntopjets","",20,0,20);

  h_pred_el_all_met = new TH1D("h_pred_el_all_met","",100,0,1000);
  h_pred_el_all_njets = new TH1D("h_pred_el_all_njets","",20,0,20);
  h_pred_el_all_mt2 = new TH1D("h_pred_el_all_mt2","",100,0,1000);
  h_pred_el_all_topmass = new TH1D("h_pred_el_all_topmass","",100,50,300);
  h_pred_el_all_ht = new TH1D("h_pred_el_all_ht","",300,0,3000);
  h_pred_el_all_mht = new TH1D("h_pred_el_all_mht","",100,0,1000);
  h_pred_el_all_ntopjets = new TH1D("h_pred_el_all_ntopjets","",20,0,20);

  h_exp_el_iso_met = new TH1D("h_exp_el_iso_met","",100,0,1000);
  h_exp_el_iso_njets = new TH1D("h_exp_el_iso_njets","",20,0,20);
  h_exp_el_iso_mt2 = new TH1D("h_exp_el_iso_mt2","",100,0,1000);
  h_exp_el_iso_topmass = new TH1D("h_exp_el_iso_topmass","",100,50,300);
  h_exp_el_iso_ht = new TH1D("h_exp_el_iso_ht","",300,0,3000);
  h_exp_el_iso_mht = new TH1D("h_exp_el_iso_mht","",100,0,1000);
  h_exp_el_iso_ntopjets = new TH1D("h_exp_el_iso_ntopjets","",20,0,20);

  h_exp_el_id_met = new TH1D("h_exp_el_id_met","",100,0,1000);
  h_exp_el_id_njets = new TH1D("h_exp_el_id_njets","",20,0,20);
  h_exp_el_id_mt2 = new TH1D("h_exp_el_id_mt2","",100,0,1000);
  h_exp_el_id_topmass = new TH1D("h_exp_el_id_topmass","",100,50,300);
  h_exp_el_id_ht = new TH1D("h_exp_el_id_ht","",300,0,3000);
  h_exp_el_id_mht = new TH1D("h_exp_el_id_mht","",100,0,1000);
  h_exp_el_id_ntopjets = new TH1D("h_exp_el_id_ntopjets","",20,0,20);

  h_exp_el_acc_met = new TH1D("h_exp_el_acc_met","",100,0,1000);
  h_exp_el_acc_njets = new TH1D("h_exp_el_acc_njets","",20,0,20);
  h_exp_el_acc_mt2 = new TH1D("h_exp_el_acc_mt2","",100,0,1000);
  h_exp_el_acc_topmass = new TH1D("h_exp_el_acc_topmass","",100,50,300);
  h_exp_el_acc_ht = new TH1D("h_exp_el_acc_ht","",300,0,3000);
  h_exp_el_acc_mht = new TH1D("h_exp_el_acc_mht","",100,0,1000);
  h_exp_el_acc_ntopjets = new TH1D("h_exp_el_acc_ntopjets","",20,0,20);

  h_exp_el_all_met = new TH1D("h_exp_el_all_met","",100,0,1000);
  h_exp_el_all_njets = new TH1D("h_exp_el_all_njets","",20,0,20);
  h_exp_el_all_mt2 = new TH1D("h_exp_el_all_mt2","",100,0,1000);
  h_exp_el_all_topmass = new TH1D("h_exp_el_all_topmass","",100,50,300);
  h_exp_el_all_ht = new TH1D("h_exp_el_all_ht","",300,0,3000);
  h_exp_el_all_mht = new TH1D("h_exp_el_all_mht","",100,0,1000);
  h_exp_el_all_ntopjets = new TH1D("h_exp_el_all_ntopjets","",20,0,20);

  h_exp_elsingle_all_met = new TH1D("h_exp_elsingle_all_met","",100,0,1000);
  h_exp_elsingle_all_njets = new TH1D("h_exp_elsingle_all_njets","",20,0,20);
  h_exp_elsingle_all_mt2 = new TH1D("h_exp_elsingle_all_mt2","",100,0,1000);
  h_exp_elsingle_all_topmass = new TH1D("h_exp_elsingle_all_topmass","",100,50,300);
  h_exp_elsingle_all_ht = new TH1D("h_exp_elsingle_all_ht","",300,0,3000);
  h_exp_elsingle_all_mht = new TH1D("h_exp_elsingle_all_mht","",100,0,1000);
  h_exp_elsingle_all_ntopjets = new TH1D("h_exp_elsingle_all_ntopjets","",20,0,20);
}

//Fill chain from txt file
bool FillChain(TChain *chain, const TString &inputFileList)
{
  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open())
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return false;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1)
  {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return true;
}


//############determine the pt bin number from generation level pt############
int Set_ptbin_number(
                     double gen_pt
                    )
{
  int ptbin_num;

  if(gen_pt < 10)
  {
    ptbin_num = 0;
  }
  else if(gen_pt >= 10 && gen_pt < 20)
  {
    ptbin_num = 1;
  }
  else if(gen_pt >= 20 && gen_pt < 30)
  {
    ptbin_num = 2;
  }
  else if(gen_pt >= 30 && gen_pt < 40)
  {
    ptbin_num = 3;
  }
  else if(gen_pt >= 40 && gen_pt < 50)
  {
    ptbin_num = 4;
  }
  else if(gen_pt >= 50 && gen_pt < 70)
  {
    ptbin_num = 5;
  }
  else if(gen_pt >= 70 && gen_pt < 100)
  {
    ptbin_num = 6;
  }
  else if(gen_pt >= 100 )
  {
    ptbin_num = 7;
  }
  
  return ptbin_num;
}


//##########functions to calculate Delta_R and Delta Phi###############
double DeltaPhi(double phi1, double phi2) 
{
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}

double DeltaR(double eta1, double phi1, double eta2, double phi2) 
{
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}

/*
bool isGoodMuonID(double id_GlobalMuonPromptTight,
                  //double isGlobalMuon,
                  double isPFMuon,
                  //double tk_chi2,
                  //double tk_numvalhits,
                  double numberOfMatchedStations,
                  double id_d0,
                  double id_tk_dz,
                  double tk_numvalPixelhits,
                  double tk_LayersWithMeasurement
                  )
{

  if((int)id_GlobalMuonPromptTight==0)
    return false;
  //if((int)isGlobalMuon==0)
    //return false;
  if((int)isPFMuon==0)
    return false;
  
  //if(tk_chi2>=10)
    //return false;
  //if(tk_numvalhits<=0)
    //return false;
  if(numberOfMatchedStations<=1)
    return false;
  if(id_d0>=0.2 || id_d0<=-0.2)
    return false;
  if(id_tk_dz>=0.5 || id_tk_dz<=-0.5)
    return false;
  if(tk_numvalPixelhits<=0)
    return false;
  if(tk_LayersWithMeasurement<=5)
    return false;
  
  return true;
}
*/
/*
bool isGoodElectronID(double isEB,
                      double isEE,
                      double dEtaIn,
                      double dPhiIn,
                      double full5x5_sigmaIetaIeta,
                      double hadOverEm
                      )
{ 
  if(( (int)isEB==1 && (int)isEE==1 )||( (int)isEB==0 && (int)isEE==0 ) )
    return false;

  if((int)isEB==1)
  {
    if(dEtaIn > 0.021)
      return false;
    if(dPhiIn > 0.25)
      return false;
    if(full5x5_sigmaIetaIeta > 0.012)
      return false;
    if(hadOverEm > 0.24)
      return false;
  }
  
  if((int)isEE==1)
  {
    if(dEtaIn > 0.028)
      return false;
    if(dPhiIn > 0.23)
      return false;
    if(full5x5_sigmaIetaIeta > 0.035)
      return false;
    if(hadOverEm > 0.19)
      return false;
  }

  return true;
}
*/
/*
bool isMatchedDeltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaR;
  deltaR = DeltaR(eta1,phi1,eta2,phi2);

  if(deltaR<0.01)
    return true;
  else
    return false;
}
*/

/*
bool isPassMuonIso(double mus_pfIsoR04_sumNeutralHadronEt,
                   double mus_pfIsoR04_sumPhotonEt,
                   double mus_pfIsoR04_sumPUPt,
                   double mus_pfIsoR04_sumChargedHadronPt,
                   double reco_mus_pt
                   )
{
  double isoNeutral = mus_pfIsoR04_sumNeutralHadronEt + mus_pfIsoR04_sumPhotonEt - 0.5*mus_pfIsoR04_sumPUPt;
  isoNeutral = ( isoNeutral > 0) ? isoNeutral : 0;

  double  muRelIso = (mus_pfIsoR04_sumChargedHadronPt + isoNeutral) / reco_mus_pt;
  if( muRelIso < 0.2 )
  {
    return true;
  }
  else
    return false;
}
*/


