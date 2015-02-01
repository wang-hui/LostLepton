#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"

#include "NTupleReader.h"

#define PT_BINS 8

class AccRecoIsoEffs
{
 public:
  void printOverview();
  void printAccRecoIsoEffs();

  //define the number of Pt bins we need in the calculation
  //static int PT_BINS = 8;

  //here we define the overall information for ttbar sample
  int nevents_tot = 0;
  int nevents_sel_base = 0;
  int nevents_sel_mus = 0;
  int nevents_sel_els = 0;

  //define acceptance, reco eff and iso eff to be calculated
  double mus_acc = 0, els_acc = 0;
  double mus_acc_err = 0, els_acc_err = 0;
  double mus_recoeff[PT_BINS] = {0}, els_recoeff[PT_BINS] = {0};
  double mus_isoeff[PT_BINS] = {0}, els_isoeff[PT_BINS] = {0};
  double mus_recoeff_err[PT_BINS] = {0}, els_recoeff_err[PT_BINS] = {0};
  double mus_isoeff_err[PT_BINS] = {0}, els_isoeff_err[PT_BINS] = {0};

  //here we defin the muon/electron number we need to count in the loop
  double nmus = 0, nmus_acc = 0, nels = 0, nels_acc =0;
  double nmus_acc_bin[PT_BINS] = {0}, nels_acc_bin[PT_BINS] = {0};
  double nmus_reco[PT_BINS] = {0}, nels_reco[PT_BINS] = {0};
  double nmus_iso[PT_BINS] = {0}, nels_iso[PT_BINS] = {0};

  void NumberstoEffs();

 private:
  double get_stat_Error(double a,
                        double an
                       );
  
  double get_sys_Error(double r,
                       double p
                      );
};



int Set_ptbin_number(double gen_pt
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

