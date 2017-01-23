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
#include "SusyAnaTools/Tools/baselineDef.h"

#include "LLBinFunction.h"

#include "SusyAnaTools/Tools/PDFUncertainty.h"

//using namespace std;

//static PDFUncertainty pdf;
//void mypdf(NTupleReader& tr)
//{
//  pdf(tr);
//}


//############################begin to defin class AccRecoIsoEffs###################
static BaselineVessel *myBaselineVessel;
void mypassBaselineFunc(NTupleReader& tr)
{
  (*myBaselineVessel)(tr);
}

class AccRecoIsoEffs
{
 public:
  void printOverview();
  void printAccRecoIsoEffs();
  void printNormalizeFlowNumber();
  void printEffsHeader();

  //here we define the overall information for ttbar sample
  int nevents_tot = 0;
  int nevents_sel_base = 0;
  int nevents_sel_mus = 0;
  int nevents_sel_els = 0;
  int nevents_cs_mus = 0;
  int nevents_cs_mus_sb[NSEARCH_BINS] = {0};

  //define acceptance, reco eff and iso eff to be calculated
  double mus_acc[NJETS_BINS][NHT_BINS] = {{0}}, els_acc[NSEARCH_BINS] = {0};
  double mus_acc_err[NJETS_BINS][NHT_BINS] = {{0}}, els_acc_err[NSEARCH_BINS] = {0};

  double mus_recoeff[PT_BINS][AC_BINS] = {{0}}, els_recoeff[PT_BINS][AC_BINS] = {{0}};
  double mus_isoeff[PT_BINS][AC_BINS] = {{0}}, els_isoeff[PT_BINS][AC_BINS] = {{0}};
  double mus_isoeff_allreco[PT_BINS][AC_BINS] = {{0}}, els_isoeff_allreco[PT_BINS][AC_BINS] = {{0}};
 
  double mus_recoeff_err[PT_BINS][AC_BINS] = {{0}}, els_recoeff_err[PT_BINS][AC_BINS] = {{0}};
  double mus_isoeff_err[PT_BINS][AC_BINS] = {{0}}, els_isoeff_err[PT_BINS][AC_BINS] = {{0}};
  double mus_isoeff_err_allreco[PT_BINS][AC_BINS] = {{0}}, els_isoeff_err_allreco[PT_BINS][AC_BINS] = {{0}};

  double mtwcorrfactor[PT_BINS] = {0}, mtwcorrfactor_err[PT_BINS] = {0};
 
  //here we define the muon/electron number we need to count in the loop
  double nmus_MC[NJETS_BINS][NHT_BINS] = {{0}}, nmus_acc_MC[NJETS_BINS][NHT_BINS] = {{0}}, nels_MC[NSEARCH_BINS] = {0}, nels_acc_MC[NSEARCH_BINS] = {0};
  double nmus_acc_bin_MC[PT_BINS][AC_BINS] = {{0}}, nels_acc_bin_MC[PT_BINS][AC_BINS] = {{0}};
  double nmus_reco_MC[PT_BINS][AC_BINS] = {{0}}, nels_reco_MC[PT_BINS][AC_BINS] = {{0}};
  double nmus_iso_MC[PT_BINS][AC_BINS] = {{0}}, nels_iso_MC[PT_BINS][AC_BINS] = {{0}};

  double nmus_MC_sb[NSEARCH_BINS] = {0};
  double nmus_sb[NSEARCH_BINS] = {0};
  double nmus_acc_MC_sb[NSEARCH_BINS] = {0};
  double nmus_acc_sb[NSEARCH_BINS] = {0};

  double nmus_reco_MC_allreco[PT_BINS][AC_BINS] = {{0}}, nels_reco_MC_allreco[PT_BINS][AC_BINS] = {{0}};
  double nmus_iso_MC_allreco[PT_BINS][AC_BINS] = {{0}}, nels_iso_MC_allreco[PT_BINS][AC_BINS] = {{0}};

  double mtwall_MC[PT_BINS] = {0}, mtw100_MC[PT_BINS] = {0};

  double nmus[NJETS_BINS][NHT_BINS] = {{0}}, nmus_acc[NJETS_BINS][NHT_BINS] = {{0}}, nels[NSEARCH_BINS] = {0}, nels_acc[NSEARCH_BINS] = {0};
  double nmus_acc_bin[PT_BINS][AC_BINS] = {{0}}, nels_acc_bin[PT_BINS][AC_BINS] = {{0}};
  double nmus_reco[PT_BINS][AC_BINS] = {{0}}, nels_reco[PT_BINS][AC_BINS] = {{0}};
  double nmus_iso[PT_BINS][AC_BINS] = {{0}}, nels_iso[PT_BINS][AC_BINS] = {{0}};

  double nmus_reco_allreco[PT_BINS][AC_BINS] = {{0}}, nels_reco_allreco[PT_BINS][AC_BINS] = {{0}};
  double nmus_iso_allreco[PT_BINS][AC_BINS] = {{0}}, nels_iso_allreco[PT_BINS][AC_BINS] = {{0}};

  double mtwall[PT_BINS] = {0}, mtw100[PT_BINS] = {0};

  //here we define the event weight we are going to use in the second loop ( muon/electron CS and prediction plots)
  double mus_EventWeight_iso[NJETS_BINS][PT_BINS][AC_BINS][NHT_BINS][NSEARCH_BINS] = {{{{{0}}}}}, mus_EventWeight_reco[NJETS_BINS][PT_BINS][AC_BINS][NHT_BINS][NSEARCH_BINS] = {{{{{0}}}}}, mus_EventWeight_acc[NJETS_BINS][PT_BINS][AC_BINS][NHT_BINS][NSEARCH_BINS] = {{{{{0}}}}};
  double els_EventWeight_iso[NJETS_BINS][PT_BINS][AC_BINS][NHT_BINS][NSEARCH_BINS] = {{{{{0}}}}}, els_EventWeight_reco[NJETS_BINS][PT_BINS][AC_BINS][NHT_BINS][NSEARCH_BINS] = {{{{{0}}}}}, els_EventWeight_acc[NJETS_BINS][PT_BINS][AC_BINS][NHT_BINS][NSEARCH_BINS] = {{{{{0}}}}};

  //here we define the search bin variables
  double nevents_mus_CS_SB_MC[NSEARCH_BINS] = {0}, nevents_mus_CS_SB_Normalized[NSEARCH_BINS] = {0};
  double nevents_mus_pred_SB_MC[NSEARCH_BINS] = {0}, nevents_mus_pred_SB_Normalized[NSEARCH_BINS] = {0};
  double nevents_mus_pred_SB_MC_err[NSEARCH_BINS] = {0};
  double nevents_mus_exp_SB_MC[NSEARCH_BINS] = {0}, nevents_mus_exp_SB_Normalized[NSEARCH_BINS] = {0};
  double nevents_els_pred_SB_MC[NSEARCH_BINS] = {0}, nevents_els_pred_SB_Normalized[NSEARCH_BINS] = {0};
  double nevents_els_exp_SB_MC[NSEARCH_BINS] = {0}, nevents_els_exp_SB_Normalized[NSEARCH_BINS] = {0};
  double nevents_lept_pred_SB_MC[NSEARCH_BINS] = {0}, nevents_lept_pred_SB_Normalized[NSEARCH_BINS] = {0};
  double nevents_lept_exp_SB_MC[NSEARCH_BINS] = {0}, nevents_lept_exp_SB_Normalized[NSEARCH_BINS] = {0};
  double nevents_lept_exp_SB_MC_isotrk[NSEARCH_BINS] = {0}, nevents_lept_exp_SB_Normalized_isotrk[NSEARCH_BINS] = {0};
  double nevents_mus_pred_iso_SB_Normalized[NSEARCH_BINS] = {0}, nevents_mus_pred_reco_SB_Normalized[NSEARCH_BINS] = {0}, nevents_mus_pred_acc_SB_Normalized[NSEARCH_BINS] = {0};
  double nevents_mus_exp_iso_SB_Normalized[NSEARCH_BINS] = {0}, nevents_mus_exp_reco_SB_Normalized[NSEARCH_BINS] = {0}, nevents_mus_exp_acc_SB_Normalized[NSEARCH_BINS] = {0};
  double nevents_mus_exp_iso_SB_Normalized_isotrk[NSEARCH_BINS] = {0}, nevents_mus_exp_reco_SB_Normalized_isotrk[NSEARCH_BINS] = {0}, nevents_mus_exp_acc_SB_Normalized_isotrk[NSEARCH_BINS] = {0};
  double w2_mus_exp_iso_SB_Normalized_isotrk[NSEARCH_BINS] = {0}, w2_mus_exp_reco_SB_Normalized_isotrk[NSEARCH_BINS] = {0}, w2_mus_exp_acc_SB_Normalized_isotrk[NSEARCH_BINS] = {0};
  double nevents_els_exp_iso_SB_Normalized_isotrk[NSEARCH_BINS] = {0}, nevents_els_exp_reco_SB_Normalized_isotrk[NSEARCH_BINS] = {0}, nevents_els_exp_acc_SB_Normalized_isotrk[NSEARCH_BINS] = {0};
  double w2_els_exp_iso_SB_Normalized_isotrk[NSEARCH_BINS] = {0}, w2_els_exp_reco_SB_Normalized_isotrk[NSEARCH_BINS] = {0}, w2_els_exp_acc_SB_Normalized_isotrk[NSEARCH_BINS] = {0};

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

/*  TFile *Effs2dPlots = new TFile("v3_Effs2dPlots.root", "recreate");
  double ptbins[PT_BINS+1]={10.0,20.0,30.0,40.0,50.0,70.0,100.0,120.0};
  double acbins[AC_BINS+2]={0.0,0.005,0.02,0.05,0.15,1.0,10.0};
  double njetbins[NJETS_BINS+1]={3.5,4.5,5.5,6.5,7.5,8.5,9.5};
  double htbins[NHT_BINS+1]={500.0,650.0,900.0,1400.0};
  TH2D *mus_recoeffs2d  = new TH2D("mus_recoeffs","Muon RecoEffs",PT_BINS,ptbins,AC_BINS+1,acbins);
  TH2D *mus_isoeffs2d  = new TH2D("mus_isoeffs","Muon IsoEffs",PT_BINS,ptbins,AC_BINS+1,acbins);
  TH2D *els_recoeffs2d  = new TH2D("els_recoeffs","Electron RecoEffs",PT_BINS,ptbins,AC_BINS+1,acbins);
  TH2D *els_isoeffs2d  = new TH2D("els_isoeffs","Electron IsoEffs",PT_BINS,ptbins,AC_BINS+1,acbins);
  TH2D *mus_acc2d  = new TH2D("mus_acc","mus_acc",NJETS_BINS,njetbins,NHT_BINS,htbins);
*/
  //TH2D *mus_recoeffs2d  = new TH2D("mus_recoeffs","Muon RecoEffs",PT_BINS,0,PT_BINS,AC_BINS,0,AC_BINS);
  //TH2D *mus_isoeffs2d  = new TH2D("mus_isoeffs","Muon IsoEffs",PT_BINS,0,PT_BINS,AC_BINS,0,AC_BINS);
  //TH2D *els_recoeffs2d  = new TH2D("els_recoeffs","Electron RecoEffs",PT_BINS,0,PT_BINS,AC_BINS,0,AC_BINS);
  //TH2D *els_isoeffs2d  = new TH2D("els_isoeffs","Electron IsoEffs",PT_BINS,0,PT_BINS,AC_BINS,0,AC_BINS);

  void NumberstoEffs();
  void EffsPlotsGen();
  void EffstoWeights();
  void EffstoWeights_fromH();
  void GetDiLeptonFactor();
  void NormalizeFlowNumber();
  //void printSearchBin(ClosureHistgram& myClosureHistgram);
  double get_stat_Error(double a, double an);

 private:
  double scale = 1;  
  double get_stat_Error_APNOA(double a, double an,double ua, double un);
  double get_sys_Error(double r, double p);
};


double AccRecoIsoEffs::get_stat_Error(double a, double an)
{
  // give the uncertainty for R=A/an=A/(A+N)
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

double AccRecoIsoEffs::get_stat_Error_APNOA(double a, double an,double ua, double un)
{
  // give the uncertainty for R=an/A=(A+N)/A
  double n;
  n = an - a;

  double err;
  err = 1000;

  double alpha;
  alpha = 1-0.6827;

  if( a>=0 && n>=0 )
  {
    err = std::sqrt(n/(a*a)*n/(a*a)*ua+1.0/(a*a)*un);
    return err;
  }
  else
  {
    return -1;
  }
}



double AccRecoIsoEffs::get_sys_Error(
                                     double r,
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
  TH1D *h_b_baseline_nMuons, *h_b_baseline_njets, *h_b_baseline_nbjetsCSVM, *h_b_baseline_MET, *h_b_baseline_jetpt2, *h_b_baseline_jetpt4, *h_b_baseline_jet1_met_phi_diff, *h_b_baseline_jet2_met_phi_diff, *h_b_baseline_jet3_met_phi_diff;
  TH1D *h_b_acc_njets, *h_b_acc_nbjetsCSVM, *h_b_acc_MET, *h_b_acc_jetpt2, *h_b_acc_jetpt4, *h_b_acc_jet1_met_phi_diff, *h_b_acc_jet2_met_phi_diff, *h_b_acc_jet3_met_phi_diff;
  TH1D *h_b_reco_nMuons, *h_b_reco_njets, *h_b_reco_nbjetsCSVM, *h_b_reco_MET, *h_b_reco_jetpt2, *h_b_reco_jetpt4, *h_b_reco_jet1_met_phi_diff, *h_b_reco_jet2_met_phi_diff, *h_b_reco_jet3_met_phi_diff;
  TH1D *h_b_deltaR_mus, *h_b_deltaR_els;
  TH1D *h_id_genactivity_mus, *h_id_genactivity_els, *h_id_recoactivity_mus, *h_id_recoactivity_els;
  TH1D *h_b_njets_mus, *h_b_njets_els;
  TH1D *h_b_jet_pt;

  TH1D *h_b_deltaR_genup_mus, *h_b_deltaR_genup_els;
  TH2D *h_b_deltaR_pt_mus, *h_b_deltaR_pt_els;
  TH2D *h_b_njets30_pt_mus, *h_b_njets30_pt_els, *h_b_njets30_eta_mus, *h_b_njets30_eta_els;
  TH1D *h_b_njets30_4_pt_mus, *h_b_njets30_5_pt_mus, *h_b_njets30_6_pt_mus, *h_b_njets30_7_pt_mus, *h_b_njets30_8_pt_mus, *h_b_njets30_9_pt_mus;
  TH1D *h_b_njets30_4_eta_mus, *h_b_njets30_5_eta_mus, *h_b_njets30_6_eta_mus, *h_b_njets30_7_eta_mus, *h_b_njets30_8_eta_mus, *h_b_njets30_9_eta_mus;
  TH1D *h_b_njets30_4_pt_els, *h_b_njets30_5_pt_els, *h_b_njets30_6_pt_els, *h_b_njets30_7_pt_els, *h_b_njets30_8_pt_els, *h_b_njets30_9_pt_els;
  TH1D *h_b_njets30_4_eta_els, *h_b_njets30_5_eta_els, *h_b_njets30_6_eta_els, *h_b_njets30_7_eta_els, *h_b_njets30_8_eta_els, *h_b_njets30_9_eta_els;
  TH1D *h_b_njets30_4_ht_mus, *h_b_njets30_5_ht_mus, *h_b_njets30_6_ht_mus, *h_b_njets30_7_ht_mus, *h_b_njets30_8_ht_mus, *h_b_njets30_9_ht_mus;

  TH1D *h_mtw_mus;
};

void BaseHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");

  h_b_all_MET = new TH1D("h_b_all_MET","",1000,0,1000);

  h_b_baseline_nMuons = new TH1D("h_b_baseline_nMuons","",10,0,10);
  h_b_baseline_njets = new TH1D("h_b_baseline_njets","",10,0,10);
  h_b_baseline_nbjetsCSVM = new TH1D("h_b_baseline_nbjetsCSVM","",10,0,10);
  h_b_baseline_MET = new TH1D("h_b_baseline_MET","",1000,0,1000);
  h_b_baseline_jetpt4 = new TH1D("h_b_baseline_jetpt4","",1000,0,1000);
  h_b_baseline_jetpt2 = new TH1D("h_b_baseline_jetpt2","",1000,0,1000);
  h_b_baseline_jet1_met_phi_diff = new TH1D("h_b_baseline_jet1_met_phi_diff","",1000,-5,5);
  h_b_baseline_jet2_met_phi_diff = new TH1D("h_b_baseline_jet2_met_phi_diff","",1000,-5,5);
  h_b_baseline_jet3_met_phi_diff = new TH1D("h_b_baseline_jet3_met_phi_diff","",1000,-5,5);

  h_b_acc_njets = new TH1D("h_b_acc_njets","",10,0,10);
  h_b_acc_nbjetsCSVM = new TH1D("h_b_acc_nbjetsCSVM","",10,0,10);
  h_b_acc_MET = new TH1D("h_b_acc_MET","",1000,0,1000);
  h_b_acc_jetpt4 = new TH1D("h_b_acc_jetpt4","",1000,0,1000);
  h_b_acc_jetpt2 = new TH1D("h_b_acc_jetpt2","",1000,0,1000);
  h_b_acc_jet1_met_phi_diff = new TH1D("h_b_acc_jet1_met_phi_diff","",1000,-5,5);
  h_b_acc_jet2_met_phi_diff = new TH1D("h_b_acc_jet2_met_phi_diff","",1000,-5,5);
  h_b_acc_jet3_met_phi_diff = new TH1D("h_b_acc_jet3_met_phi_diff","",1000,-5,5);

  h_b_reco_nMuons = new TH1D("h_b_reco_nMuons","",10,0,10);
  h_b_reco_njets = new TH1D("h_b_reco_njets","",10,0,10);
  h_b_reco_nbjetsCSVM = new TH1D("h_b_reco_nbjetsCSVM","",10,0,10);
  h_b_reco_MET = new TH1D("h_b_reco_MET","",1000,0,1000);
  h_b_reco_jetpt4 = new TH1D("h_b_reco_jetpt4","",1000,0,1000);
  h_b_reco_jetpt2 = new TH1D("h_b_reco_jetpt2","",1000,0,1000);
  h_b_reco_jet1_met_phi_diff = new TH1D("h_b_reco_jet1_met_phi_diff","",1000,-5,5);
  h_b_reco_jet2_met_phi_diff = new TH1D("h_b_reco_jet2_met_phi_diff","",1000,-5,5);
  h_b_reco_jet3_met_phi_diff = new TH1D("h_b_reco_jet3_met_phi_diff","",1000,-5,5);

  h_b_deltaR_mus = new TH1D("h_b_deltaR_mus","",1000,0,0.5);
  h_b_deltaR_els = new TH1D("h_b_deltaR_els","",1000,0,0.5);
  
  h_b_deltaR_genup_mus = new TH1D("h_b_deltaR_genup_mus","",50,0,5);
  h_b_deltaR_genup_els = new TH1D("h_b_deltaR_genup_els","",50,0,5);

  h_b_deltaR_pt_mus = new TH2D("h_b_deltaR_pt_mus","",100,0,1,100,0,500);
  h_b_deltaR_pt_els = new TH2D("h_b_deltaR_pt_els","",100,0,1,100,0,500);

  h_b_njets30_pt_mus = new TH2D("h_b_njets30_pt_mus","",15,0,15,20,0,200);
  h_b_njets30_pt_els = new TH2D("h_b_njets30_pt_els","",15,0,15,20,0,200);
  h_b_njets30_eta_mus = new TH2D("h_b_njets30_eta_mus","",15,0,15,12,-3,3);
  h_b_njets30_eta_els = new TH2D("h_b_njets30_eta_els","",15,0,15,12,-3,3);

  h_b_njets30_4_pt_mus = new TH1D("h_b_njets30_4_pt_mus","",200,0,200);
  h_b_njets30_5_pt_mus = new TH1D("h_b_njets30_5_pt_mus","",200,0,200);
  h_b_njets30_6_pt_mus = new TH1D("h_b_njets30_6_pt_mus","",200,0,200);
  h_b_njets30_7_pt_mus = new TH1D("h_b_njets30_7_pt_mus","",200,0,200);
  h_b_njets30_8_pt_mus = new TH1D("h_b_njets30_8_pt_mus","",200,0,200);
  h_b_njets30_9_pt_mus = new TH1D("h_b_njets30_9_pt_mus","",200,0,200);

  h_b_njets30_4_ht_mus = new TH1D("h_b_njets30_4_ht_mus","",200,0,2000);
  h_b_njets30_5_ht_mus = new TH1D("h_b_njets30_5_ht_mus","",200,0,2000);
  h_b_njets30_6_ht_mus = new TH1D("h_b_njets30_6_ht_mus","",200,0,2000);
  h_b_njets30_7_ht_mus = new TH1D("h_b_njets30_7_ht_mus","",200,0,2000);
  h_b_njets30_8_ht_mus = new TH1D("h_b_njets30_8_ht_mus","",200,0,2000);
  h_b_njets30_9_ht_mus = new TH1D("h_b_njets30_9_ht_mus","",200,0,2000);

  h_b_njets30_4_eta_mus = new TH1D("h_b_njets30_4_eta_mus","",60,-3,3);
  h_b_njets30_5_eta_mus = new TH1D("h_b_njets30_5_eta_mus","",60,-3,3);
  h_b_njets30_6_eta_mus = new TH1D("h_b_njets30_6_eta_mus","",60,-3,3);
  h_b_njets30_7_eta_mus = new TH1D("h_b_njets30_7_eta_mus","",60,-3,3);
  h_b_njets30_8_eta_mus = new TH1D("h_b_njets30_8_eta_mus","",60,-3,3);
  h_b_njets30_9_eta_mus = new TH1D("h_b_njets30_9_eta_mus","",60,-3,3);

  h_b_njets30_4_pt_els = new TH1D("h_b_njets30_4_pt_els","",200,0,200);
  h_b_njets30_5_pt_els = new TH1D("h_b_njets30_5_pt_els","",200,0,200);
  h_b_njets30_6_pt_els = new TH1D("h_b_njets30_6_pt_els","",200,0,200);
  h_b_njets30_7_pt_els = new TH1D("h_b_njets30_7_pt_els","",200,0,200);
  h_b_njets30_8_pt_els = new TH1D("h_b_njets30_8_pt_els","",200,0,200);
  h_b_njets30_9_pt_els = new TH1D("h_b_njets30_9_pt_els","",200,0,200);

  h_b_njets30_4_eta_els = new TH1D("h_b_njets30_4_eta_els","",60,-3,3);
  h_b_njets30_5_eta_els = new TH1D("h_b_njets30_5_eta_els","",60,-3,3);
  h_b_njets30_6_eta_els = new TH1D("h_b_njets30_6_eta_els","",60,-3,3);
  h_b_njets30_7_eta_els = new TH1D("h_b_njets30_7_eta_els","",60,-3,3);
  h_b_njets30_8_eta_els = new TH1D("h_b_njets30_8_eta_els","",60,-3,3);
  h_b_njets30_9_eta_els = new TH1D("h_b_njets30_9_eta_els","",60,-3,3);

  h_mtw_mus = new TH1D("h_mtw_mus","",200,0,200);

  h_id_genactivity_mus = new TH1D("h_id_genactivity_mus","",200,0,1);
  h_id_genactivity_els = new TH1D("h_id_genactivity_els","",200,0,1);
  h_id_recoactivity_mus = new TH1D("h_id_recoactivity_mus","",200,0,1);
  h_id_recoactivity_els = new TH1D("h_id_recoactivity_els","",200,0,1);

  h_b_njets_mus = new TH1D("h_b_njets_mus","",40,0,40);
  h_b_njets_els = new TH1D("h_b_njets_els","",40,0,40);

  h_b_jet_pt = new TH1D("h_b_jet_pt","",1000,0,200);
}


class ClosureHistgram
{
 public:
  void BookHistgram(const char *);

  TFile *oFile;
  //closure plots definition
  TH1D *h_pred_mu_iso_met, *h_pred_mu_iso_njets, *h_pred_mu_iso_mt2, *h_pred_mu_iso_ht, *h_pred_mu_iso_mht, *h_pred_mu_iso_ntopjets;
  TH1D *h_pred_mu_id_met, *h_pred_mu_id_njets, *h_pred_mu_id_mt2, *h_pred_mu_id_ht, *h_pred_mu_id_mht, *h_pred_mu_id_ntopjets;
  TH1D *h_pred_mu_acc_met, *h_pred_mu_acc_njets, *h_pred_mu_acc_mt2, *h_pred_mu_acc_ht, *h_pred_mu_acc_mht, *h_pred_mu_acc_ntopjets;
  TH1D *h_pred_mu_all_met, *h_pred_mu_all_njets, *h_pred_mu_all_mt2, *h_pred_mu_all_ht, *h_pred_mu_all_mht, *h_pred_mu_all_ntopjets;

  TH1D *h_exp_mu_iso_met, *h_exp_mu_iso_njets, *h_exp_mu_iso_mt2, *h_exp_mu_iso_ht, *h_exp_mu_iso_mht, *h_exp_mu_iso_ntopjets;
  TH1D *h_exp_mu_id_met, *h_exp_mu_id_njets, *h_exp_mu_id_mt2, *h_exp_mu_id_ht, *h_exp_mu_id_mht, *h_exp_mu_id_ntopjets;
  TH1D *h_exp_mu_acc_met, *h_exp_mu_acc_njets, *h_exp_mu_acc_mt2, *h_exp_mu_acc_ht, *h_exp_mu_acc_mht, *h_exp_mu_acc_ntopjets;
  TH1D *h_exp_mu_all_met, *h_exp_mu_all_njets, *h_exp_mu_all_mt2, *h_exp_mu_all_ht, *h_exp_mu_all_mht, *h_exp_mu_all_ntopjets;

  TH1D *h_exp_musingle_all_met, *h_exp_musingle_all_njets, *h_exp_musingle_all_mt2, *h_exp_musingle_all_ht, *h_exp_musingle_all_mht, *h_exp_musingle_all_ntopjets;

  TH1D *h_pred_el_iso_met, *h_pred_el_iso_njets, *h_pred_el_iso_mt2, *h_pred_el_iso_ht, *h_pred_el_iso_mht, *h_pred_el_iso_ntopjets;
  TH1D *h_pred_el_id_met, *h_pred_el_id_njets, *h_pred_el_id_mt2, *h_pred_el_id_ht, *h_pred_el_id_mht, *h_pred_el_id_ntopjets;
  TH1D *h_pred_el_acc_met, *h_pred_el_acc_njets, *h_pred_el_acc_mt2, *h_pred_el_acc_ht, *h_pred_el_acc_mht, *h_pred_el_acc_ntopjets;
  TH1D *h_pred_el_all_met, *h_pred_el_all_njets, *h_pred_el_all_mt2, *h_pred_el_all_ht, *h_pred_el_all_mht, *h_pred_el_all_ntopjets;

  TH1D *h_exp_el_iso_met, *h_exp_el_iso_njets, *h_exp_el_iso_mt2, *h_exp_el_iso_ht, *h_exp_el_iso_mht, *h_exp_el_iso_ntopjets;
  TH1D *h_exp_el_id_met, *h_exp_el_id_njets, *h_exp_el_id_mt2, *h_exp_el_id_ht, *h_exp_el_id_mht, *h_exp_el_id_ntopjets;
  TH1D *h_exp_el_acc_met, *h_exp_el_acc_njets, *h_exp_el_acc_mt2, *h_exp_el_acc_ht, *h_exp_el_acc_mht, *h_exp_el_acc_ntopjets;
  TH1D *h_exp_el_all_met, *h_exp_el_all_njets, *h_exp_el_all_mt2, *h_exp_el_all_ht, *h_exp_el_all_mht, *h_exp_el_all_ntopjets;

  TH1D *h_pred_lept_all_met, *h_pred_lept_all_njets, *h_pred_lept_all_mt2, *h_pred_lept_all_ht, *h_pred_lept_all_mht, *h_pred_lept_all_ntopjets;
  TH1D *h_exp_lept_all_met, *h_exp_lept_all_njets, *h_exp_lept_all_mt2, *h_exp_lept_all_ht, *h_exp_lept_all_mht, *h_exp_lept_all_ntopjets;
  TH1D *h_exp_lept_all_Z_sb, *h_exp_lept_all_Z_njets30, *h_exp_lept_all_Z_njets50, *h_exp_lept_all_Z_ntops, *h_exp_lept_all_Z_nbjets, *h_exp_lept_all_Z_MET, *h_exp_lept_all_Z_MT2, *h_exp_lept_all_Z_HT;

  TH2D *h_pred_lept_all_2d_met_mupt;
  TH2D *h_pred_el_all_2d_met_njets;

  TH1D *h_exp_elsingle_all_met, *h_exp_elsingle_all_njets, *h_exp_elsingle_all_mt2, *h_exp_elsingle_all_ht, *h_exp_elsingle_all_mht, *h_exp_elsingle_all_ntopjets;

  //closure for search bin
  TH1D *h_exp_mu_sb, *h_pred_mu_sb, *h_exp_el_sb, *h_pred_el_sb, *h_exp_lept_sb, *h_pred_lept_sb, *h_exp_lept_sb_isotrk, *h_pred_lept_sb_isotrk;
  TH1D * h_pred_mu_iso_sb, *h_pred_mu_reco_sb, *h_pred_mu_acc_sb, *h_exp_iso_mu_sb, *h_exp_reco_mu_sb, *h_exp_acc_mu_sb;


};

void ClosureHistgram::BookHistgram(const char *outFileName)
{

  oFile = new TFile(outFileName, "recreate");
  //start closure plots
  h_pred_mu_iso_met = new TH1D("h_pred_mu_iso_met","",100,0,1000);
  h_pred_mu_iso_njets = new TH1D("h_pred_mu_iso_njets","",20,0,20);
  h_pred_mu_iso_mt2 = new TH1D("h_pred_mu_iso_mt2","",100,0,1000);
  h_pred_mu_iso_ht = new TH1D("h_pred_mu_iso_ht","",300,0,3000);
  h_pred_mu_iso_mht = new TH1D("h_pred_mu_iso_mht","",100,0,1000);
  h_pred_mu_iso_ntopjets = new TH1D("h_pred_mu_iso_ntopjets","",20,0,20);

  h_pred_mu_id_met = new TH1D("h_pred_mu_id_met","",100,0,1000);
  h_pred_mu_id_njets = new TH1D("h_pred_mu_id_njets","",20,0,20);
  h_pred_mu_id_mt2 = new TH1D("h_pred_mu_id_mt2","",100,0,1000);
  h_pred_mu_id_ht = new TH1D("h_pred_mu_id_ht","",300,0,3000);
  h_pred_mu_id_mht = new TH1D("h_pred_mu_id_mht","",100,0,1000);
  h_pred_mu_id_ntopjets = new TH1D("h_pred_mu_id_ntopjets","",20,0,20);

  h_pred_mu_acc_met = new TH1D("h_pred_mu_acc_met","",100,0,1000);
  h_pred_mu_acc_njets = new TH1D("h_pred_mu_acc_njets","",20,0,20);
  h_pred_mu_acc_mt2 = new TH1D("h_pred_mu_acc_mt2","",100,0,1000);
  h_pred_mu_acc_ht = new TH1D("h_pred_mu_acc_ht","",300,0,3000);
  h_pred_mu_acc_mht = new TH1D("h_pred_mu_acc_mht","",100,0,1000);
  h_pred_mu_acc_ntopjets = new TH1D("h_pred_mu_acc_ntopjets","",20,0,20);

  h_pred_mu_all_met = new TH1D("h_pred_mu_all_met","",100,0,1000);
  h_pred_mu_all_njets = new TH1D("h_pred_mu_all_njets","",20,0,20);
  h_pred_mu_all_mt2 = new TH1D("h_pred_mu_all_mt2","",100,0,1000);
  h_pred_mu_all_ht = new TH1D("h_pred_mu_all_ht","",300,0,3000);
  h_pred_mu_all_mht = new TH1D("h_pred_mu_all_mht","",100,0,1000);
  h_pred_mu_all_ntopjets = new TH1D("h_pred_mu_all_ntopjets","",20,0,20);

  h_exp_mu_iso_met = new TH1D("h_exp_mu_iso_met","",100,0,1000);
  h_exp_mu_iso_njets = new TH1D("h_exp_mu_iso_njets","",20,0,20);
  h_exp_mu_iso_mt2 = new TH1D("h_exp_mu_iso_mt2","",100,0,1000);
  h_exp_mu_iso_ht = new TH1D("h_exp_mu_iso_ht","",300,0,3000);
  h_exp_mu_iso_mht = new TH1D("h_exp_mu_iso_mht","",100,0,1000);
  h_exp_mu_iso_ntopjets = new TH1D("h_exp_mu_iso_ntopjets","",20,0,20);

  h_exp_mu_id_met = new TH1D("h_exp_mu_id_met","",100,0,1000);
  h_exp_mu_id_njets = new TH1D("h_exp_mu_id_njets","",20,0,20);
  h_exp_mu_id_mt2 = new TH1D("h_exp_mu_id_mt2","",100,0,1000);
  h_exp_mu_id_ht = new TH1D("h_exp_mu_id_ht","",300,0,3000);
  h_exp_mu_id_mht = new TH1D("h_exp_mu_id_mht","",100,0,1000);
  h_exp_mu_id_ntopjets = new TH1D("h_exp_mu_id_ntopjets","",20,0,20);

  h_exp_mu_acc_met = new TH1D("h_exp_mu_acc_met","",100,0,1000);
  h_exp_mu_acc_njets = new TH1D("h_exp_mu_acc_njets","",20,0,20);
  h_exp_mu_acc_mt2 = new TH1D("h_exp_mu_acc_mt2","",100,0,1000);
  h_exp_mu_acc_ht = new TH1D("h_exp_mu_acc_ht","",300,0,3000);
  h_exp_mu_acc_mht = new TH1D("h_exp_mu_acc_mht","",100,0,1000);
  h_exp_mu_acc_ntopjets = new TH1D("h_exp_mu_acc_ntopjets","",20,0,20);

  h_exp_mu_all_met = new TH1D("h_exp_mu_all_met","",100,0,1000);
  h_exp_mu_all_njets = new TH1D("h_exp_mu_all_njets","",20,0,20);
  h_exp_mu_all_mt2 = new TH1D("h_exp_mu_all_mt2","",100,0,1000);
  h_exp_mu_all_ht = new TH1D("h_exp_mu_all_ht","",300,0,3000);
  h_exp_mu_all_mht = new TH1D("h_exp_mu_all_mht","",100,0,1000);
  h_exp_mu_all_ntopjets = new TH1D("h_exp_mu_all_ntopjets","",20,0,20);

  h_exp_lept_all_met = new TH1D("h_exp_lept_all_met","",100,0,1000);
  h_exp_lept_all_njets = new TH1D("h_exp_lept_all_njets","",20,0,20);
  h_exp_lept_all_mt2 = new TH1D("h_exp_lept_all_mt2","",100,0,1000);
  h_exp_lept_all_ht = new TH1D("h_exp_lept_all_ht","",300,0,3000);
  h_exp_lept_all_mht = new TH1D("h_exp_lept_all_mht","",100,0,1000);
  h_exp_lept_all_ntopjets = new TH1D("h_exp_lept_all_ntopjets","",20,0,20);

  h_exp_lept_all_Z_sb = new TH1D("hSearchBins"      , "Search Bins;Search Bin;Events"        , 45 , 0 , 45);
  h_exp_lept_all_Z_njets30 = new TH1D("hNJets30"    , "NJets30;N_{jets} (p_{T} > 30);Events" , 10 , 0    , 10);   // "cntNJetsPt30Eta24"
  h_exp_lept_all_Z_njets50 = new TH1D("hNJets50"    , "NJets50;N_{jets} (p_{T} > 50);Events" , 10 , 0    , 10);   // "cntNJetsPt50Eta24"
  h_exp_lept_all_Z_ntops = new TH1D("hNTops"        , "NTops;N_{tops};Events"                , 5  , 0    , 5);    // "nTopCandSortedCnt"
  h_exp_lept_all_Z_nbjets = new TH1D("hNbJets"      , "NbJets;N_{bjets};Events"              , 5  , 0    , 5);    // "cntCSVS"
  h_exp_lept_all_Z_MET = new TH1D("hMET"            , "MET;#slash{E}_{T} [GeV];Events"       , 24 , 200  , 800);  // "met"
  h_exp_lept_all_Z_MT2 = new TH1D("hMT2"            , "MT2;M_{T2} [GeV];Events"              , 24 , 200  , 800);  // "best_had_brJet_MT2"
  h_exp_lept_all_Z_HT = new TH1D("hHT"              , "HT;H_{T} [GeV];Events"                , 20 , 500  , 1000); // "HT"

  h_exp_musingle_all_met = new TH1D("h_exp_musingle_all_met","",100,0,1000);
  h_exp_musingle_all_njets = new TH1D("h_exp_musingle_all_njets","",20,0,20);
  h_exp_musingle_all_mt2 = new TH1D("h_exp_musingle_all_mt2","",100,0,1000);
  h_exp_musingle_all_ht = new TH1D("h_exp_musingle_all_ht","",300,0,3000);
  h_exp_musingle_all_mht = new TH1D("h_exp_musingle_all_mht","",100,0,1000);
  h_exp_musingle_all_ntopjets = new TH1D("h_exp_musingle_all_ntopjets","",20,0,20);

  h_pred_el_iso_met = new TH1D("h_pred_el_iso_met","",100,0,1000);
  h_pred_el_iso_njets = new TH1D("h_pred_el_iso_njets","",20,0,20);
  h_pred_el_iso_mt2 = new TH1D("h_pred_el_iso_mt2","",100,0,1000);
  h_pred_el_iso_ht = new TH1D("h_pred_el_iso_ht","",300,0,3000);
  h_pred_el_iso_mht = new TH1D("h_pred_el_iso_mht","",100,0,1000);
  h_pred_el_iso_ntopjets = new TH1D("h_pred_el_iso_ntopjets","",20,0,20);

  h_pred_el_id_met = new TH1D("h_pred_el_id_met","",100,0,1000);
  h_pred_el_id_njets = new TH1D("h_pred_el_id_njets","",20,0,20);
  h_pred_el_id_mt2 = new TH1D("h_pred_el_id_mt2","",100,0,1000);
  h_pred_el_id_ht = new TH1D("h_pred_el_id_ht","",300,0,3000);
  h_pred_el_id_mht = new TH1D("h_pred_el_id_mht","",100,0,1000);
  h_pred_el_id_ntopjets = new TH1D("h_pred_el_id_ntopjets","",20,0,20);

  h_pred_el_acc_met = new TH1D("h_pred_el_acc_met","",100,0,1000);
  h_pred_el_acc_njets = new TH1D("h_pred_el_acc_njets","",20,0,20);
  h_pred_el_acc_mt2 = new TH1D("h_pred_el_acc_mt2","",100,0,1000);
  h_pred_el_acc_ht = new TH1D("h_pred_el_acc_ht","",300,0,3000);
  h_pred_el_acc_mht = new TH1D("h_pred_el_acc_mht","",100,0,1000);
  h_pred_el_acc_ntopjets = new TH1D("h_pred_el_acc_ntopjets","",20,0,20);

  h_pred_el_all_met = new TH1D("h_pred_el_all_met","",100,0,1000);
  h_pred_el_all_njets = new TH1D("h_pred_el_all_njets","",20,0,20);
  h_pred_el_all_mt2 = new TH1D("h_pred_el_all_mt2","",100,0,1000);
  h_pred_el_all_ht = new TH1D("h_pred_el_all_ht","",300,0,3000);
  h_pred_el_all_mht = new TH1D("h_pred_el_all_mht","",100,0,1000);
  h_pred_el_all_ntopjets = new TH1D("h_pred_el_all_ntopjets","",20,0,20);

  h_pred_lept_all_met = new TH1D("h_pred_lept_all_met","",100,0,1000);
  h_pred_lept_all_2d_met_mupt = new TH2D("h_pred_lept_all_2d_met_mupt","",100,-10,10,100,0,100);
  h_pred_el_all_2d_met_njets = new TH2D("h_pred_el_all_2d_met_njets","",100,0,1000,20,0,20);
  h_pred_lept_all_njets = new TH1D("h_pred_lept_all_njets","",20,0,20);
  h_pred_lept_all_mt2 = new TH1D("h_pred_lept_all_mt2","",100,0,1000);
  h_pred_lept_all_ht = new TH1D("h_pred_lept_all_ht","",300,0,3000);
  h_pred_lept_all_mht = new TH1D("h_pred_lept_all_mht","",100,0,1000);
  h_pred_lept_all_ntopjets = new TH1D("h_pred_lept_all_ntopjets","",20,0,20);

  h_exp_el_iso_met = new TH1D("h_exp_el_iso_met","",100,0,1000);
  h_exp_el_iso_njets = new TH1D("h_exp_el_iso_njets","",20,0,20);
  h_exp_el_iso_mt2 = new TH1D("h_exp_el_iso_mt2","",100,0,1000);
  h_exp_el_iso_ht = new TH1D("h_exp_el_iso_ht","",300,0,3000);
  h_exp_el_iso_mht = new TH1D("h_exp_el_iso_mht","",100,0,1000);
  h_exp_el_iso_ntopjets = new TH1D("h_exp_el_iso_ntopjets","",20,0,20);

  h_exp_el_id_met = new TH1D("h_exp_el_id_met","",100,0,1000);
  h_exp_el_id_njets = new TH1D("h_exp_el_id_njets","",20,0,20);
  h_exp_el_id_mt2 = new TH1D("h_exp_el_id_mt2","",100,0,1000);
  h_exp_el_id_ht = new TH1D("h_exp_el_id_ht","",300,0,3000);
  h_exp_el_id_mht = new TH1D("h_exp_el_id_mht","",100,0,1000);
  h_exp_el_id_ntopjets = new TH1D("h_exp_el_id_ntopjets","",20,0,20);

  h_exp_el_acc_met = new TH1D("h_exp_el_acc_met","",100,0,1000);
  h_exp_el_acc_njets = new TH1D("h_exp_el_acc_njets","",20,0,20);
  h_exp_el_acc_mt2 = new TH1D("h_exp_el_acc_mt2","",100,0,1000);
  h_exp_el_acc_ht = new TH1D("h_exp_el_acc_ht","",300,0,3000);
  h_exp_el_acc_mht = new TH1D("h_exp_el_acc_mht","",100,0,1000);
  h_exp_el_acc_ntopjets = new TH1D("h_exp_el_acc_ntopjets","",20,0,20);

  h_exp_el_all_met = new TH1D("h_exp_el_all_met","",100,0,1000);
  h_exp_el_all_njets = new TH1D("h_exp_el_all_njets","",20,0,20);
  h_exp_el_all_mt2 = new TH1D("h_exp_el_all_mt2","",100,0,1000);
  h_exp_el_all_ht = new TH1D("h_exp_el_all_ht","",300,0,3000);
  h_exp_el_all_mht = new TH1D("h_exp_el_all_mht","",100,0,1000);
  h_exp_el_all_ntopjets = new TH1D("h_exp_el_all_ntopjets","",20,0,20);

  h_exp_elsingle_all_met = new TH1D("h_exp_elsingle_all_met","",100,0,1000);
  h_exp_elsingle_all_njets = new TH1D("h_exp_elsingle_all_njets","",20,0,20);
  h_exp_elsingle_all_mt2 = new TH1D("h_exp_elsingle_all_mt2","",100,0,1000);
  h_exp_elsingle_all_ht = new TH1D("h_exp_elsingle_all_ht","",300,0,3000);
  h_exp_elsingle_all_mht = new TH1D("h_exp_elsingle_all_mht","",100,0,1000);
  h_exp_elsingle_all_ntopjets = new TH1D("h_exp_elsingle_all_ntopjets","",20,0,20);

  //closure for search bin
  h_exp_mu_sb = new TH1D("h_exp_mu_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
  h_pred_mu_sb = new TH1D("h_pred_mu_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
  h_exp_el_sb = new TH1D("h_exp_el_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
  h_pred_el_sb = new TH1D("h_pred_el_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
  h_exp_lept_sb = new TH1D("h_exp_lept_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
  h_exp_lept_sb_isotrk = new TH1D("h_exp_lept_isotrk","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
  h_pred_lept_sb = new TH1D("h_pred_lept_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
  h_pred_lept_sb_isotrk = new TH1D("h_pred_lept_isotrk","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);

  h_pred_mu_iso_sb = new TH1D("h_pred_mu_iso_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
  h_pred_mu_reco_sb = new TH1D("h_pred_mu_reco_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
  h_pred_mu_acc_sb = new TH1D("h_pred_mu_acc_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);

  h_exp_iso_mu_sb = new TH1D("h_exp_mu_iso_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
  h_exp_reco_mu_sb = new TH1D("h_exp_mu_reco_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
  h_exp_acc_mu_sb = new TH1D("h_exp_mu_acc_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);

}

//Common function that will be used in both calculation and expectation loop
class LostLeptonObj
{
 public:
  //initialize those variables with some crazy number
  bool isMu = false, isEl = false;
  double gen_eta = -10, gen_phi = -10, gen_pt = -10, gen_activity = -10;
  double reco_eta = -10, reco_phi = -10, reco_pt = -10, reco_activity = -10;
  bool passAcc = false, passId = false, passIso = false;
  bool doneAcc = false, doneId = false, doneIso = false;
  int reco_index = -10;  

  void SetMyLL(
                //LL variables that will be always useful 
                int pid,
                TLorentzVector onegenlept,
                double genpfActivity,
                std::vector<int> IdFlag,
                std::vector<TLorentzVector> recoleptLVec,
                std::vector<double> recoleptactivityVec,
                std::vector<double> MiniIso
              );

 private:
  void SetFlavor(int pid);
  void genLeptonSetup(TLorentzVector onegenlept, double activity);
  void gogoAcc();
  void gogoId(std::vector<int> IdFlag, std::vector<TLorentzVector> leptLVec);
  void recoLeptonSetup(TLorentzVector onerecolept, double activity);
  void gogoIso(std::vector<double> iso);
};


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
